#### About the program

#### Input files to be provided: 
A **SAM file** where the following columns contain the following information: <br>
•	(3) RNAME: the chromosome to which the read aligned <br>
•	(4) POS: the position on the chromosome where the alignment starts <br>
•	(6) CIGAR: a string describing the alignment in cigar format, see below <br>
•	Last column is NH:i:x: an integer describing how many times this read aligned to the reference<br>
<br>
A **tab-separated file**, containing a header and 3 columns <br>
(1)	Gene ID <br>
(2)	Transcript ID <br>
(3)	Location of the gene <br>
in the following format: GeneID_chrNum:startposition..endposition(strand) <br>
<br>
For instance TGME49_chrVIII:6,793,066..6,795,596(-) where:<br>
TGME49_chrVIII is name of the chromosome where the gene is encoded, <br>
6793066 is the start position, <br>
6795596 is the end position and <br>
(-) is the strand the gene is located <br>

#### Input format: <br>
The input files are to be given at the command line in the following order: <br>
````
$ python ./<Script.py> <samFile> <locations>
````


#### Output 

A **tab separated file** showing the locations of all splice junctions found withing the boundaries of each gene. An empty row is inserted at the end of every gene. 
The columns present in the output file are: <br>
(1)	Gene id <br>
(2)	Junction Start <br>
(3)	Junction End <br>
(4)	Number of reads supporting the junction <br>
<br>



