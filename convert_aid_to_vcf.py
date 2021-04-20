#!/usr/bin/env python

import argparse
import os
import sys
import time
from time import gmtime, strftime
from math import ceil
###########################
# Convert AID output to VCF
###########################
# @ authors: skadri
# Copyright: The University of Chicago
# This script is part of the variant annotation pipeline. We assume the data in the following format:

#Chr     Pos     OrigDepth       Depth>30        RefAllele       RefAllele%      HighVarAllele   HighVarAllele%  Comments
#chr1    115258747       5557    5082    C       0.941755214     G       0.056867375
#chr5    170837543       3620    4286    C       0.755482968     +4TCTG  0.24265049      TCTG inserted after C
#chr12   121432114       1126    1095    C       0.946118721     -1G     0.052968037     G deleted after C
#chr11   108123558       604     590     T       0.913559322     *       0.083050847     skip; print those cases in a separate output file
#chr13   28608263        4191    4191    T       0.620138392     +24TCATATTCTCTGAAATCAACGTAG     0.379861608     long indel

def main(inputfile,outputf):

	# Error checks for files existence
	if (not os.path.exists(inputfile)):
		print inputfile + " does not exist!"
		sys.exit(1)

	varsfile=open(inputfile,"r")
	# THERE IS NO HEADER
	# line=varsfile.readline()

	outputfile=open(outputf,"w")
	# Write the header
	outputfile.write("##fileformat=VCFv4.1\n")
	outputfile.write("##fileDate="+strftime("%Y-%m-%d %H:%M:%S", gmtime())+"\n")
	outputfile.write("##source="+inputfile+"\n")
	outputfile.write("##reference=file:///"+"\n")
	outputfile.write("##contig=\n")
	outputfile.write("##phasing=partial\n")
	outputfile.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n")
	outputfile.write("##INFO=<ID=DP30,Number=1,Type=Integer,Description=\"Total Depth higher than 30x\">\n")
	outputfile.write("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n")
	outputfile.write("##FILTER=<ID=q10,Description=\"Quality below 10\">\n")
	outputfile.write("##FILTER=<ID=s50,Description=\"Less than 50% of samples have data\">\n")
	outputfile.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
	outputfile.write("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n")
	outputfile.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n")
	outputfile.write("##FORMAT=<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">\n")
	outputfile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmy_sample\n")

	line=varsfile.readline()
	while(line):
		if("Failed indel identification" not in line):
			line=line.strip()
			entries=line.split("\t")
			# chr  Pos  OrigDepth  Depth_30  RefAllele  RefAllele_p  HighVarAllele  HighVarAllele_p
			# 0   1	 2	  3	 4		 5	   6		  7

			chrom=entries[0]
			pos=entries[1]
			DP=entries[2]
			DP30=entries[3]
			maf=entries[7]

			# ROUND THE MAF to 3 decimal places
			maf=str(ceil(float(maf) * 1000) / 1000)

			ref_allele=entries[4]
			var_allele=None
			high_varallele=entries[6]
			# if insertion, add the REF to the ALT and remove the +
			if("+" in entries[6]):
				# REF_ALLELE IS AS IS
				# VAR_ALLELE is REF_ALLELE+HIGH_VARALLELE

				# ignore the "+"
				high_varallele=high_varallele[1:]
				# ignore the numbers
				while(high_varallele[0].isdigit()):
					high_varallele=high_varallele[1:]
				var_allele=ref_allele+high_varallele
			else:
				# if deletion
				if("-" in entries[6]):
					# REF_ALLELE IS REF_ALLELE+HIGH_VARALLELE
					# VAR_ALLELE IS REF_ALLELE
					var_allele=ref_allele

					# ignore the "-"
					high_varallele=high_varallele[1:]
					# ignore the numbers
					while(high_varallele[0].isdigit()):
						high_varallele=high_varallele[1:]
					ref_allele=ref_allele+high_varallele
				# NO INSERTION OR DELETION
				else:
					var_allele=high_varallele

			# OUTPUT:
			# chr  Pos   .  RefAllele  VarAllele  .  .  DP=OrigDepth;DP30=Depth_30;AF=HighVarAllele_p;REF=RefAllele;ALT=VarAllele;  GT:DP   :
			outputfile.write(chrom+"\t"+pos+"\t"+"."+"\t"+ref_allele+"\t"+var_allele+"\t"+"."+"\t"+"."+"\t")
			outputfile.write("DP="+DP+";DP30="+DP30+";AF="+maf+";REF="+ref_allele+";ALT="+var_allele+";"+"\t"+"GT:DP"+"\t"+":"+"\n")

		line=varsfile.readline()


	outputfile.close()
	varsfile.close()


# Returns the directory of the script
#def getScriptPath():
#	return os.path.dirname(os.path.realpath(sys.argv[0]))

# Adds two classes to the formatter_class
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
	pass

if __name__ == "__main__":
	version_string="v1.0.1"

	description_string="##########################################################\n| Convert variant stats output to VCF |\n##########################################################"

	usage_short="convert_aid_to_vcf.py -i <Input variantstats file. Required> -o <Output VCF name.Required>"

	# RawDescriptionHelpFormatter allows custom formatting of Description text
	parser = argparse.ArgumentParser(formatter_class=CustomFormatter,description=description_string,usage=usage_short)
	parser.add_argument("-v","--version", action="store_true", help="Current version of convert_aid_to_vcf")

	# Read arguments
	parser.add_argument("-i",dest="input", help="Input variantstats file. Required.")
	parser.add_argument("-o", dest="out", help="Output VCF file. Required.")
	args = parser.parse_args()

	# if user asks for the version
	if args.version:
		print version_string
		sys.exit(0)

	# Print version always
	print "Converting AID output to VCF file..."
	print version_string

	# if user did not input all required arguments
	if (args.input == None) or (args.out == None):
		parser.print_help()
		sys.exit(1)

	main(args.input, args.out)
