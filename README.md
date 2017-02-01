
Start - Stop codon analysis
---------------------------

Workflow for start, stop codon analysis for Ribosome Protected Fragments (RPFs). The pipeline used here follows part of the method adopted by 
[RiboTaper workflow](https://ohlerlab.mdc-berlin.de/software/RiboTaper_126/) for start and stop codon analysis.

**Description**

For an input .bam file, count the number of mapped reads and sub-sample if the number of mapped reads is greater than the given threshold.
Generate a list of reads that maps to a window of 100nt around start and stop codons from BED file. Bin reads according to read lengths and calculate the distance between 5'/3' 
ends of the reads to start or stop codon positions. Compute coverage data per read length bin for read 5'/3' ends to start - stop positions.

**Dependencies**
+ [SAMtools](http://samtools.sourceforge.net/)  
+ [bedtools](http://bedtools.readthedocs.io/en/latest/)
+ Python (>=2.7)
  + [jinja2](https://pypi.python.org/pypi/Jinja2) (tested version: 2.7.2)
  + [matplotlib](https://pypi.python.org/pypi/matplotlib) (tested version: 1.4.3)
  + [numpy](https://pypi.python.org/pypi/numpy) (tested vesion: 1.9.1)
  + [pandas](https://pypi.python.org/pypi/pandas) (tested verion: 0.18.1)
 
**Usage**  
Complete workflow:  
`start_stop_analysis.sh <RPF.bam>  <start_stop.bed> <sample_name>  <downsample_size>`  
Python helper script options:  
`python start_stop_analysis.py -h`

