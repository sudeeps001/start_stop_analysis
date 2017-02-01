#!/usr/bin/env python
import argparse
import cPickle as pickle
import gzip
import jinja2
import json
import logging
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import re
import random
import sys
import time
from functools import wraps

'''
Created on Jan 11, 2017
Python helper file for start and stop codon analysis from ribosome profiling data
This module is designed to parse the results from the run
bamToBed -i RPF.bam -bed12 -split | \
    windowBed -w 100 -sm -b stdin -a start_stop.bed | \
    cut -f 7- |sort -k1,1 -k2,2g | \
    closestBed -s -t "last" -a stdin -b start_stop.bed > codon_analysis.csv

Copyright (C) 2016  Sudeep Sahadevan
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>
'''
logger = logging.getLogger()

def logfn(fn):
    @wraps(fn)
    def _wrapper(*args,**kwargs):
        if len(kwargs) > 0:
            logging.debug('function:  {}{}{}'.format(fn.__name__,args,kwargs))
        else:
            logging.debug('function:  {}{}'.format(fn.__name__,args))
        _fn = fn(*args,**kwargs)
        return _fn
    return _wrapper

def timefn(fn):
    @wraps(fn)
    def _wrapper(*args,**kwargs):
        logging.debug('Function: {} starting run @ {}'.format(fn.__name__,time.strftime('%Y-%m-%d %H:%M:%S')))
        res =  fn(*args,**kwargs)
        logging.debug('Function: {} finished run @ {}'.format(fn.__name__,time.strftime('%Y-%m-%d %H:%M:%S')))
        return res
    return _wrapper

class StartStopData():
    
    def __init__(self):
        self.startcodons = {}
        self.stopcodons = {}
        self.__start_len_dict = {} # read length to distance dictionary for start codons
        self.__start_dist = set() # distance set for start codons
        self.__stop_len_dict = {} # read length to distance dictionary for stop codons
        self.__stop_dist = set() # distance set for stop codons
    
    def add_start_codon(self,chr,startpos,strand,readstart,readlen):
        '''
        Given chr id/name, start positon of start codon (5'/3') w.r.t strand, strand information,
        read start position (w.r.t strand) and read length add this data to startcodon dictionary
        '''
        ids = '|'.join([chr,str(startpos),strand])
        dist = startpos-readstart if strand =='-' else (readstart-startpos)
        try:
            self.startcodons[ids] = np.concatenate((self.startcodons[ids],[[startpos,readstart,dist,readlen]]))
        except KeyError:
            self.startcodons[ids] = np.array([[startpos,readstart,dist,readlen]])
#        read length to start distance
        self.__start_dist.add(dist)
            
    def add_stop_codon(self,chr,stoppos,strand,readstop,readlen):
        '''
        Given chr id/name, start positon of start codon (5'/3') w.r.t strand, strand information,
        read stop position (w.r.t strand) and read length add this data to stopcodon dictionary
        '''
        ids = '|'.join([chr,str(stoppos),strand])
        dist = stoppos-readstop if strand =='-' else (readstop-stoppos)
        try:
            self.stopcodons[ids] = np.concatenate((self.stopcodons[ids],[[stoppos,readstop,dist,readlen]]))
        except KeyError:
            self.stopcodons[ids] = np.array([[stoppos,readstop,dist,readlen]])
#        read length to stop distance
        self.__stop_dist.add(dist)
    
    def get_start_codon_dist_count(self,mincount=20):
        '''
        return a pandas data.frame where columns are read 5' positions to start codon 5' distance ,
        indices are read lengths, and values are the coverage per read 5' - start codon 5' distance
        '''
        start_len_dict = {}
        for atg in self.startcodons.keys():
            dat = self.startcodons[atg]
            if dat.shape[0]>=mincount:
                for e in dat:
                    try:
                        start_len_dict[int(e[3])].append(int(e[2]))
                    except KeyError:
                        start_len_dict[int(e[3])] = [int(e[2])]
                
        cols = sorted(self.__start_dist)
        sdf = pd.DataFrame(columns=cols)
        readlens = sorted(start_len_dict)
        for rl in readlens:
            distuniq,distcount = np.unique(start_len_dict[rl], return_counts=True)
            tmpdf = pd.DataFrame(columns=distuniq)
            tmpdf.loc[0] = distcount
            tmpdf.index=[rl]
            sdf = sdf.append(tmpdf)
        sdf = sdf.fillna(0)
        return sdf
    
    def get_stop_codon_dist_count(self,mincount=20):
        '''
        return a pandas data.frame where columns are read 3' positions to stop codon 3' distance ,
        indices are read lengths, and values are the coverage per read 5' - stop codon 3' distance
        '''
        stop_len_dict = {}
        for tga in self.stopcodons.keys():
            dat = self.stopcodons[tga]
            if dat.shape[0]>=mincount:
                for e in dat:
                    try:
                        stop_len_dict[int(e[3])].append(int(e[2]))
                    except KeyError:
                        stop_len_dict[int(e[3])] = [int(e[2])]
                
        cols = sorted(self.__stop_dist)
        sdf = pd.DataFrame(columns=cols)
        readlens = sorted(stop_len_dict)
        for rl in readlens:
            distuniq,distcount = np.unique(stop_len_dict[rl], return_counts=True)
            tmpdf = pd.DataFrame(columns=distuniq)
            tmpdf.loc[0] = distcount
            tmpdf.index=[rl]
            sdf = sdf.append(tmpdf)
        sdf = sdf.fillna(0)
        return sdf

@logfn
@timefn
def get_cmap(N):
    '''
    Returns a function that maps each index in 0, 1, ... N-1 to a distinct 
    RGB color. colmap list generated from http://tools.medialab.sciences-po.fr/iwanthue/
    '''
    colmap = ["#60002e","#56ee8f","#7d2194","#a5dc57","#6938aa","#58b845","#a53cac","#31b44e","#ab63d9","#6a9500","#0054c5","#fad048",
            "#4886ff","#9b8f00","#00348d","#d7da6d","#7f0076","#00dca6","#c00f79","#71eba8","#650061","#45edc6","#d11b58","#01d3b5",
            "#cf2445","#008429","#e18fff","#567500","#ff8df1","#004d04","#ea4da7","#00965a","#ff74c7","#7c6800","#ba9fff","#ffa741",
            "#025db1","#ff9b4f","#003978","#d15c1d","#0067a6","#ff8852","#450353","#e8d487","#3b124f","#ffaf7f","#9caeff","#963c00",
            "#e3b3ff","#684b00","#ff9dd0","#580300","#ff86aa","#782300","#612b63","#ff7b5b","#88005a","#bd7451","#8d0039","#ff999a",
            "#91000f","#ff7081","#8d0020","#ff726a"]
    if N<=len(colmap):
        return [colmap[i] for i in random.sample(range(0,len(colmap)),N)]
    else:
        max_value = 16581375 #255**3
        interval = int(max_value / N)
        colors = ['#'+hex(I)[2:].zfill(6) for I in range(0, max_value, interval)]
        return colors

@logfn
@timefn
def plot_pdf(df,filename,title='Figure title',cols = None):
    '''
    Given a pandas dataFrame object plot each index values as line plot
    '''
    if cols is None:
        cols = get_cmap(len(df))
    fg,subp = plt.subplots()
    subp.tick_params(axis='both', which='major', labelsize=5)
    subp.set_ylim([0,int(round(df.values.max()+100,-1))])
    subp.set_xticks(df.columns.values)
    subax = subp.axes
    subax.xaxis.set_ticks_position('none')
    subax.xaxis.grid(which='major',alpha=0.2,linestyle='-')
    for i in range(0,len(df)):
        subp.plot(df.columns.values,df.iloc[i],linewidth=1,label=df.index[i],color=cols[i])
    subp.set_xlabel('Distance to codon')
    subp.set_ylabel('Coverage')
    subp.set_title(title)
    subp.legend(loc=0,ncol=3,prop={'size':6},borderaxespad=0.)
    plt.savefig(filename, format='pdf', bbox_inches='tight')

@logfn
@timefn
def plot_d3(coveragepd,colmap,basename,plottype):
    '''
    return D3 library based js/html for interactive visualization.
    D3: https://d3js.org/
    '''
    if colmap is None:
        colmap = get_cmap(len(coveragepd))
    readLenbins = json.dumps(list(coveragepd.index.astype(np.str)))
    lenColors = json.dumps(colmap)
    distBins = json.dumps(list(coveragepd.columns))
    covData = json.dumps(np.array(coveragepd).tolist()) # readLenbins,lenColors,distBins,covData
    title = basename+' '+plottype
    header = basename+' '+plottype+' distance plot'
    xlab = plottype+' distance'
    htmltemplate = jinja2.Template('''
    <html>
    <head>
        <script type="text/javascript" src="http://d3js.org/d3.v3.js" charset="utf-8"></script>
        <link href='https://fonts.googleapis.com/css?family=Open+Sans' rel='stylesheet' type='text/css'>
        <!--
        Copyright (C) 2016  Sudeep Sahadevan
        This program is free software: you can redistribute it and/or modify
        it under the terms of the GNU General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.

        This program is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        GNU General Public License for more details.
        
        You should have received a copy of the GNU General Public License
        along with this program.  If not, see <http://www.gnu.org/licenses/>
        -->
        <meta charset="utf-8"/>
        <title> {{title}} </title>
    </head>
    <body>
        <div id="cbselector">
            <div id="formSelector">
                <div id="dataSelector"></div>
                <div id="plotOptions">
                    <form id="inputf1" autocomplete="off"> Select plot type:<br>
                        <input value="line" name="plotChbx" id="pC1" type="checkbox" checked="checked">line
                        <input value="dot" name="plotChbx" id="pC2" type="checkbox" checked="checked">dot
                        <input value="area" name="plotChbx" id="pC3" type="checkbox">area<br>
                    </form>
                    <label for="inputf2"> Select Y axis scale:<br>
                    <select id="inputf2" autocomplete="off">
                        <option value="linear"  selected="selected">linear</option>
                        <option value="log">log</option>
                    </select>
                </div>
            </div>
            <div id="plotter"></div>
            <style type="text/css">
                head, body {font-family: 'Open Sans', sans-serif;}
                .axis{font-size: 14px;}
                .axis path, .axis line {fill: none; stroke: grey; stroke-width: 2; shape-rendering: crispEdges;}
                .area {opacity: 0.30;}                
                .line {fill: none; stroke-width: 2; opacity: 0.60}
            </style>
            <script type="text/javascript" >
                var ids1 = {{readLenbins}};
                var idSelector = [];
                var idColors = {{lenColors}};
                var dataIndices = {{distBins}};
                var dataMat = {{covData}};
                function plotRender (dataIndices,dataMat,idSel,plotL,plotD,plotA,ySel) {
                    var plotData = [];
                    var plotCol = [];
                    var plotNames = [];
                    var minMax = [];
                    var firstDat = true;
                    for (i=0;i<idSel.length;i++) {
                            plotData[i] = dataMat[idSel[i]];
                            if (ySel==="log"){ // add 1 to all values in case there are zeros
                                plotData[i] = plotData[i].map(function (x) {return x+1});
                            }
                        plotCol[i] = idColors[idSel[i]];
                        plotNames.push(ids1[idSel[i]]);
                        if (firstDat) {
                            minMax = d3.extent(dataMat[idSel[i]])
                            firstDat = false;
                        }else {
                            tmpDat = d3.extent(dataMat[idSel[i]])
                            if (tmpDat[0] < minMax[0]) {
                                minMax[0] = tmpDat[0]
                            }
                            if (tmpDat[1] > minMax[1]) {
                                minMax[1] = tmpDat[1]
                            }
                        }
                    }
                    var yScale;
                    var yaxis;
                    var yLabelString = "Coverage";
                    if (ySel==="linear") {
                        yScale = d3.scale.linear().domain(minMax).range([height,margin.top]).nice();
                        yaxis = d3.svg.axis().scale(yScale).orient("left");
                    }else {
                        yScale = d3.scale.log().range([height,margin.top]).domain([minMax[0]+1,minMax[1]+1]).nice();
                        yaxis = d3.svg.axis().scale(yScale).orient("left");
                        yaxis.tickFormat(function (d) {return yScale.tickFormat(4,d3.format(",d"))(d) }).tickSize(5,0);
                        yLabelString += " (log scale)"
                    }
                    var area = d3.svg.area()
                                    .x(function(d,i) {return xScale(dataIndices[i]);})
                                    .y0(height)
                                    .y1(function(d) {return yScale(d);});
                    minMaxPos = d3.extent(dataIndices)
                    var xScale = d3.scale.linear().domain(minMaxPos).range([10, width]);
                    var xaxis = d3.svg.axis().scale(xScale).ticks(5);                        
                    var line = d3.svg.line()
                            .x(function(d,i) {return xScale(dataIndices[i]);})
                            .y(function(d) {return yScale(d)});                                   
                    var lines = svg.selectAll( "g" ) .data( plotData );
                    var dots = svg.selectAll( "g" ) .data( plotData );
                    var areas = svg.selectAll( "g" ) .data( plotData );
                    svg.append("g")
                        .attr("class", "x axis")
                        .attr("transform", "translate(0," + height + ")")
                        .call(xaxis);
                    svg.append("g")
                        .attr("class", "y axis")
                        .call(yaxis);
                    if (plotL) {
                        var aLineContainer = lines.enter().append("g");
                        aLineContainer.append("path")
                            .attr("class", "line")
                            .attr("d", line)
                            .style("stroke", function(d,i) { return plotCol[i]; });
                    }
                    if (plotD) {// last so elements would be on top
                        var aDotContainer = dots.enter().append("g");
                        aDotContainer.selectAll("circle")
                        .data(function (d,i) {return d;})
                        .enter()
                        .append("circle")
                        .attr("class","dots")
                        .attr("r",4)
                        .attr("cx",function (d,i) {return xScale(dataIndices[i]);})
                        .attr("cy",function (d,i,j) {return yScale(d);})
                        .style("fill",function (d) {
                            var index = aDotContainer.data().indexOf(d3.select(this.parentNode).datum());
                            return plotCol[index];
                        })//console.log(aDotContainer.data().indexOf(d3.select(this.parentNode).datum())); access parent node index
                        .on("click", handleClick)
                        .on("mouseout", handleMouseOut);
//                        handle clicking events on data points
                        function handleClick(d, i) {  // Add interactivity
//                         Use D3 to select element, change color and size
                        d3.select(this).attr({r: 6});
                        var textAreaUpdate = document.getElementById("dotInfo1");
                        textAreaUpdate.style.fontSize = "13px"
                        var dotInfo = "<table id=infoTable><tr><td style=\\"color:#47484a;\\">Distance to codon:</td><td style=\\"color:#47484a;\\">"+dataIndices[i]+"</td></tr>";
                        dotInfo += "<table id=infoTable><tr><th style=\\"color:#47484a;\\">Read length</td><th style=\\"color:#47484a;\\">Count</td></tr>";
                        for (d=0;d<plotData.length;d++) {
                            var tmpI = plotData[d];
                            var tmpVal = tmpI[i];
                            if (ySel==='log'){ tmpVal = tmpVal-1 }
                            var tmpCount = d3.format(",d")(tmpVal);
                            dotInfo += "<tr><td style=\\"color:"+plotCol[d]+";\\" align=\\"center\\" >"+plotNames[d]+"</td><td style=\\"color:#47484a;\\">"+tmpCount+"</td></tr>";
                        }
                        dotInfo +="</table>";
                        textAreaUpdate.innerHTML = dotInfo;
                    }
//                    resize data point on mouse out
                    function handleMouseOut(d, i) {  // Add interactivity
//                         Use D3 to select element, change color and size
                        d3.select(this).attr({r: 4});
                    }
                    }
                    if (plotA) {
                        var aAreasContainer = areas.enter().append("g");
                        aAreasContainer.append("path")
                            .attr("class", "area")
                            .attr("d", area)
                            .style("fill", function(d,i) { return plotCol[i]; });
                    }
//                    add title
                    svg.append("text")
                        .attr("class","title")
                        .attr("x", (width/2))
                        .attr("y", margin.top/2)
                        .attr("text-anchor", "middle")
                        .style("font-size", "18px")
                        .style("font-weight", "bold")
                        .text(\"{{header}}\");
//                    add x axis label 
                    svg.append("text")
                        .attr("class","x label")
                        .attr("x", (width / 2))
                        .attr("y", height+(margin.bottom*0.8))
                        .attr("text-anchor", "middle")
                        .style("font-size", "16px")
                        .text("{{xlab}}");
//                    add y axis label 
                    svg.append("text")
                        .attr("class","y label")
                        .attr("x", 0-height/2)
                        .attr("y", 0-margin.left*0.8)
                        .attr("text-anchor", "middle")
                        .attr("transform","rotate(-90)")
                        .style("font-size", "16px")
                        .text(yLabelString);        
                }
//                dynamically add dataset names to form selector
                var container = document.getElementById("dataSelector");
                container.class="container";
                contLabel = document.createElement("label")
                contLabel.htmlfor="dataSelector";
                contLabel.appendChild(document.createTextNode("Select read lengths to plot:"));
                container.appendChild(contLabel);
                container.appendChild(document.createElement("br"));
//                .container {width:300px; height: 50px; overflow-y: scroll; }
                var finput = document.createElement("form");
                    finput.id = "cboxInput1";
                    finput.autocomplete="off";
                    finput.class="container";
                    finput.style.paddingTop="5px";
                    finput.style.width="360px";
                    finput.style.height="75px";
                    finput.style.overflowY="scroll";
                    finput.style.overflowX="scroll";
                    finput.style.paddingBottom="5px";
                container.appendChild(finput);
                for(i=0;i<ids1.length;i++ ){
//                    idColors[i] = color1.get(true);
                    var cbox = document.createElement("input");
                    cbox.type = "checkbox";
                    cbox.name = "datCB";
                    cbox.value = i;
                    cbox.style.display = "inline-block";
                    cbox.id = ids1[i].replace(/\s{1,}/g,"");
                    var label = document.createElement('label');
                    label.htmlfor = "cboxid"+i; 
                    label.style.color = idColors[i];
                    label.style.display = "inline-block";
                    label.appendChild(document.createTextNode(ids1[i]));
                    finput.appendChild(cbox);// add checkbox before label to avoid dealing with string length
                    finput.appendChild(label); 
//                    finput.appendChild(document.createElement("br"));
                    if (i===0) {
                        document.getElementById(ids1[0].replace(/\s{1,}/g,"")).checked=true;
                        idSelector.push(0);
                    }
                    if ((i+1)%5===0){finput.appendChild(document.createElement("br"));}
                    else { label.appendChild(document.createTextNode(' '));}
                }
//                add update button
                var but1 = document.createElement("BUTTON");
                but1.id = "updater1";
                but1.appendChild(document.createTextNode("Update"));
                but1.onclick = callUpdater;
//                add information area
                var textArea = document.createElement("div");
                textArea.id="dotInfo1";
                textArea.contenteditable="true";
                textArea.style.fontSize = "13px";
                textArea.style.color = "#47484a";
                textArea.style.paddingTop="10px";
                textArea.innerHTML="Click on a data dot (from dot plot) for details.";
                var optDiv = document.getElementById("formSelector")
                optDiv.appendChild(but1);
                optDiv.appendChild(textArea);
//                plotting stuff
                W1 = 1050;
                H1 = 600;                
                var margin = {top: 50, right: 20, bottom: 50, left: 100},
                                    width = W1 - margin.left - margin.right,
                                    height = H1 - margin.top - margin.bottom;
                var svg = d3.select("#plotter")
                                .append("svg")
                                .attr("height", H1)
                                .attr("width", W1)
                                .append("g")
                                .attr("transform", "translate(" + margin.left + "," + margin.top + ")");
                var formPos = document.getElementById("formSelector");
                formPos.style.position="absolute";
                formPos.style.width="400px";
                formPos.style.left=W1+5;
                formPos.style.top=margin.top+30;
//                predefine y scale
                var ySelectorVal = "linear"
//                predefine what to plot
                var plotLine = true;
                var plotDot = true;
                var plotArea = false;
                plotRender (dataIndices,dataMat,idSelector,plotLine,plotDot,plotArea,ySelectorVal);
//                function to redraw the plots on user input
                function callUpdater () {
                    idSelector = [];
                    var cbVals = document.getElementsByName("datCB");
                    for (var j=0; j<cbVals.length; j++) {
                        if (cbVals[j].checked) {
                            idSelector.push(parseInt(cbVals[j].value));
                        }
                    }
//                     if all the datasets are unselected, select the first one by default
                    if (idSelector.length==0) {
                        document.getElementById(ids1[0].replace(/\s{1,}/g,"")).checked=true;
                        idSelector.push(0);
                    }
                    var plotChecker = document.getElementsByName("plotChbx");
                    for (var p=0;p<plotChecker.length;p++) {
                        if (plotChecker[p].value==="dot") {
                            plotDot = plotChecker[p].checked;
                        }else if (plotChecker[p].value==="area") {
                            plotArea = plotChecker[p].checked;
                        }else {
                            plotLine = plotChecker[p].checked;
                        }
                    }
//                    if both the plot options are unselected, select line plot by default 
                    if (plotLine===false && plotDot===false && plotArea===false){ 
                        plotLine = true;
                        plotDot = true;
                        document.getElementById("pC1").checked=plotLine;
                        document.getElementById("pC2").checked=plotDot;
                    }
                    yVals = document.getElementById("inputf2");
                    ySelectorVal = yVals.options[yVals.selectedIndex].value;
                    svg.selectAll("g > *").remove();
                    plotRender(dataIndices,dataMat,idSelector,plotLine,plotDot,plotArea,ySelectorVal);
                }
            </script>
        </div>
    </body>
</html>
    ''')
    return htmltemplate.render(readLenbins=readLenbins,lenColors=lenColors,distBins=distBins,covData=covData,title=title,header=header,xlab=xlab)

@logfn
@timefn
def parse_start_stop_codon(fin,startcodon='start_codon',stopcodon='stop_codon',minQ=10,minL=20):
    '''
    Given an output from closestBed in bedtools suite, parse the file.
    Expected columns in order: 
    '''
    readids = set() # do not count reads more than once
    startstop = StartStopData()
    featPat = re.compile('^.*\t('+startcodon+'|'+stopcodon+')\t.*$',re.IGNORECASE)
    fh = gzip.open(fin, 'rb') if re.search('^.*.gz',fin) else open(fin,'r')
    for f in fh:
        f = f.strip('\n')
        if re.search(featPat, f):
            fdat = f.split('\t')
            qual = float(fdat[4])
            if (fdat[0]==fdat[12]) and (fdat[5]==fdat[17]) and qual >=minQ: # just make sure
                stype = fdat[15]
                if (stype==startcodon and fdat[5]=='+') or (stype==stopcodon and fdat[5]=='-'):
                    fdist = int(fdat[13]) - int(fdat[6]) # +ve for start pos - read start pos in +ve strand, 
                elif (stype==startcodon and fdat[5]=='-') or (stype==stopcodon and fdat[5]=='+'):
                    fdist = int(fdat[14]) - int(fdat[7]) # -ve for stop pos - read stop pos in -ve strand
                readlen = sum([int(s) for s in fdat[10].split(',')])
                if not readids.__contains__(fdat[3]):
                    if abs(fdist)<=readlen and fdist>=0 and stype==startcodon and fdat[5]=='+' and readlen>=minL:
                        startstop.add_start_codon(fdat[0], int(fdat[13]), fdat[5], int(fdat[6]), readlen)
                    elif abs(fdist)<=readlen and fdist<=0 and stype==stopcodon and fdat[5]=='+' and readlen>=minL:
                        startstop.add_stop_codon(fdat[0], int(fdat[14]), fdat[5], int(fdat[7]), readlen)
                    elif abs(fdist)<=readlen and fdist<=0 and stype==startcodon and fdat[5]=='-' and readlen>=minL:
                        startstop.add_start_codon(fdat[0], int(fdat[14]), fdat[5], int(fdat[7]), readlen)
                    elif abs(fdist)<=readlen and fdist>=0 and stype==stopcodon and fdat[5]=='-' and readlen>=minL:
                        startstop.add_stop_codon(fdat[0], int(fdat[13]), fdat[5], int(fdat[6]), readlen)
                    readids.add(fdat[3])
    fh.close()
    return startstop

@logfn
@timefn
def check_file(fname):
    '''
    If file exists, log a file rewrite warning
    '''
    if os.path.exists(fname):
        logging.warn('Rewriting existing file: %s'%fname)
    return fname

def main(argv):
    prog = re.sub('^.*\/','',argv[0])
    version = '0.1.1a'
    description = '''Script for start/stop codon coverage analysis
    The input file for this script is the output file from the run:
    bamToBed -i RPF.bam -bed12 -split | windowBed -w 100 -sm -b stdin -a start_stop.bed | cut -f 7- |sort -k1,1 -k2,2g | closestBed -s -t "last" -a stdin -b start_stop.bed > codon_analysis.csv
    The output file is expected to have columns in the following order: 
    chr, start, end, read.id, map.quality, strand, .1, .2, .3, span.exons, exon.length, intron.length, chr.start_stop, begin.start_stop, end.start_stop, type.start_stop, gene_id.start_stop, strand.start_stop
    
    Plot options:
        pdf: Generate static pdf plot
        html: Generate interactive D3 plot
    
    Output file(s) naming schema:
        start codon analysis files: Base_name.start_codon.csv (.pdf|.html)
        stop codon analysis files: Base_name.stop_codon.csv (.pdf|.html)
    '''
    plottype=['pdf','html']
    loglevel = ['debug','info','warn','quiet']
    rd = argparse.ArgumentParser(prog, description=description,version=version,formatter_class = argparse.RawTextHelpFormatter)
    rd.add_argument('--f',metavar='Codon coverage',dest='in',help='Result file from: bamToBed -i RPF.bam -bed12 -split | windowBed -w 100 -sm -b stdin -a start_stop.bed | cut -f 7- |sort -k1,1 -k2,2g | closestBed -s -t "last" -a stdin -b start_stop.bed > codon_analysis.csv (supports .gz files)',required=True)
    rd.add_argument('--start',metavar='Start codon name tag',dest='start',help='Start codon name tag used in BED file (default: start_codon)',default='start_codon')
    rd.add_argument('--stop',metavar='Stop codon name tag',dest='stop',help='Stop codon name tag used in BED file (default: stop_codon)',default='stop_codon')
    rd.add_argument('--l',metavar='Min. read length',dest='minl',help='Prune alignments with readlength < min. read length (default: 20)',default=20,type=int)
    rd.add_argument('--q',metavar='Min. quality',dest='minq',help='Prune alignments with alignment quality < min. quality  (default: 10)',default=10,type=int)
    rd.add_argument('--c',metavar='Min. count',dest='minc',help='Prune start/stop codons with < min. count reads mapped (default: 0)',default=0,type=int)
    rd.add_argument('--plot',dest='plot',metavar='Choose plot type. Allowed options: '+", ".join(plottype)+' (default:pdf).',default=plottype[0])
    rd.add_argument('--name',metavar='Base name',dest='name',help='Base name to save output files',required=True)
    rd.add_argument('--verbose',metavar='Log level',dest='verbose',help='Allowed choices: '+", ".join(loglevel)+" (default: info)",default=loglevel[1])
    rvars = vars(rd.parse_args())
    fin = rvars['in']
    start = rvars['start']
    stop = rvars['stop']
    minL = rvars['minl']
    minQ = rvars['minq']
    minC = rvars['minc']
    plotType = rvars['plot']
    basename = rvars['name']
    verbose = rvars['verbose']
    if verbose=='quiet':
        ch = logging.StreamHandler(logging.NullHandler)
        logger.addHandler(ch)
    else:
        logger.setLevel(logging.getLevelName(verbose.upper()))
        ch = logging.StreamHandler(sys.stdout)
        ch.setLevel(logging.getLevelName(verbose.upper()))
        formatter = logging.Formatter(' [%(levelname)s]  %(message)s')
        ch.setFormatter(formatter)
        logger.addHandler(ch)
    basename = re.sub('\s{1,}', '_', basename)
    logging.info(' Parsing file: %s'%os.path.abspath(fin))
    logging.info(' Start codon tag: %s'%start)
    logging.info(' Stop codon tag: %s'%stop)
    logging.info(' Minimum read length: %d'%minL)
    logging.info(' Minimum quality: %f'%minQ)
    logging.info(' Minimum read count: %d'%minC)
    logging.info(' Basename: %s'%basename)
    logging.info(' Plot type: %s'%plotType)
    codon_analysis = parse_start_stop_codon(fin, start, stop, minQ, minL)
    start_pd = codon_analysis.get_start_codon_dist_count(minC)
    stop_pd = codon_analysis.get_stop_codon_dist_count(minC)
    start_csv = check_file(basename+'_start_codon.csv')
    stop_csv = check_file(basename+'_stop_codon.csv')
    start_pd.to_csv(start_csv, "\t", encoding="utf-8")
    stop_pd.to_csv(stop_csv, "\t", encoding="utf-8")
    cols = get_cmap(len(start_pd)) if len(start_pd) == len(stop_pd)  else None # so that both plots have same color
    logging.info(' *** Output files ***')
    logging.info(' Start codon analysis csv file: %s '%start_csv)
    logging.info(' Stop codon analysis csv file: %s '%stop_csv)
    if plotType=='pdf':
        start_pdf =  check_file(basename+'_start_codon.pdf')
        plot_pdf(start_pd, start_pdf, basename+' start codon analysis',cols)
        stop_pdf = check_file(basename+'_stop_codon.pdf')
        plot_pdf(stop_pd, stop_pdf, basename+' stop codon analysis',cols)
        logging.info(' Start codon pdf plot: %s '%start_pdf)
        logging.info(' Stop codon pdf plot: %s '%stop_pdf)
    elif plotType=='html':
        startcodonplot = basename+'_start_codon.html'
        fightml = open(startcodonplot,'w')
        fightml.write(plot_d3(start_pd, cols, basename, 'start codon'))
        fightml.close()
        stopcodonplot = basename+'_stop_codon.html'
        fightml = open(stopcodonplot,'w')
        fightml.write(plot_d3(stop_pd, cols, basename, 'stop codon'))
        fightml.close()
        logging.info(' Start codon html plot: %s '%startcodonplot)
        logging.info(' Stop codon html plot: %s '%stopcodonplot)
    
if __name__ == '__main__':
    main(sys.argv)