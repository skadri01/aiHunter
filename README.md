# Amplicon Indel Hunter (AIH/AID) - Indel detection in amplicon-based, paired-end, NGS data

## aiHunter v1.1.0

### Contact
Sabah Kadri: skadri@bsd.uchicago.edu

Jeremy P Segal: jsegal@bsd.uchicago.edu

### **Publication (please cite if using software):**
Sabah Kadri, Chao J. Zhen, Michelle N. Wurst, Bradley C. Long, Zi-Feng Jiang, Y. Lynn Wang, Larissa V. Furtado, Jeremy P. Segal, **Amplicon Indel Hunter Is a Novel Bioinformatics Tool to Detect Large Somatic Insertion/Deletion Mutations in Amplicon-Based Next-Generation Sequencing Data**, The Journal of Molecular Diagnostics, Volume 17, Issue 6, November 2015, Pages 635-643, ISSN 1525-1578, 
http://dx.doi.org/10.1016/j.jmoldx.2015.06.005.\newline
(http://www.sciencedirect.com/science/article/pii/S1525157815001506)


## Introduction 
Accurate detection of large insertions and deletions (indels) via amplicon-based targeted NGS assays remains a challenge when depending on alignment-based methods. Sequencing reads that cover these indels are, by definition, different from the reference sequence, and lead to variable performance of alignment algorithms. Amplicon Indel Hunter (AIH) is a large (>5-bp) indel detection method that is reference genome independent and highly sensitive for the identification of somatic indels in amplicon-based, paired-end, NGS data. The software (aiHunter) takes as input paired end fastq files and information about the amplicons in the assay, and detects amplicons with potential large indel. The output from AIH is given as input to a helper tool, Amplicon Indel Diagnoser (AID) which tries to annotate the exact indel sequence and returns output in VCF format.


## Release Notes
v1.1.0 has the following changes over v1.0.0:

1. **Bug fix in VCF file generation:**  
If AID fails to decode the exact indel, due to bad quality or the size of the indel, it returns a  
``Failed indel identification - manual review required for <amplicon name>``  
line in the \<read1filename\>.R2fastq.finalindelstats file.  
In v1.0.0, the VCF file generation was failing if there was a decoding failure by AID. This is now fixed. The VCF file will only contain entries for the indels which are successfully decoded by AID. However, we encourage all users to also investigate the \<read1 filename\>.R2fastq.finalindelstats file to manually check entries which failed with AID.    
2. **Example data:**  
Example data, and its associated documentation is now included in the aiHunter folder. Please see the folder aiHunter_example_data.
3. **New FAQ section**
We have added a new FAQ section to this documentation, based on common issues faced by some of our users. Please refer to the FAQ section for better understanding of any errors encounters. Please feel free to reach out Sabah Kadri (skadri@bsd.uchicago.edu) for any additional trouble- shooting.


## Quickstart
A quick way to start using aiHunter is to have paired end Fastq files, an infosheet with details of the assay, and a fasta file with the insert sequences. Please see "Input Data" below for more details.  

``python aiHunter.py --read1 <Read1 Fastq> --read2 <Read2 Fastq> --amp <Amplicon information sheet> --inserts <Inserts Fasta file> [options]``


## Installation
aiHunter has been tested with python v2.7.6. It requires the following packages in the python PATH to run it: (1) argparse (2) os (3) sys (4) time (5) Bio (biopython) (6) math  
An executable binary is also provided with its own documentation in the binary_distribution folder. This can only be run on the Linux OS.


## License 
See License.txt

## Input Data
aiHunter requires the following input files in order to run:  

|     |     |
|-----|-----|
|``--read1``  |Read1 Fastq file.   |
|``--read2``  |Read2 Fastq file.  |
|``--amp``  |Tab-delimited amplicon info file. See below.   |
|``--inserts``  |Fasta file with insert sequences. See below.   |
  
     
   
       
   
