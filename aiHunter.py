#!/usr/bin/env python

import argparse
import os
import sys
import time
import amplicon_indel_hunter
import amplicon_indel_diagnoser
import convert_aid_to_vcf

###########################
# Amplicon Indel Hunter (aiHunter) v1.1.0
###########################
# @ authors: skadri, jsegal
# copyright The University of Chicago
# Amplicon Indel Hunter is a reference genome-independent large (>5-bp) indel detection method
# used for the identification of somatic indels in amplicon-based, paired-end,
# NGS data. It must be used with professional care and judgement.

# This program takes as input paired end fastq files from amplicon assays and
# a tab-delimited information file of amplicon data in the following format:
# AmpliconName  Primer1Sequence  Primer2Sequence  AmpliconLength Chr  GenomicAmpliconStart  GenomicAmpliconEnd
# The input is not required to be quality filtered or primer matched, but can be done by the users separately.

def main(R1fastq, R2fastq, infofile, inserts,outdir, cushion, sig_threshold):

	# Error checks for files existence
	if (not os.path.exists(R1fastq)):
		print R1fastq + " does not exist!"
		sys.exit(1)

	if (not os.path.exists(args.R2fastq)):
		print R2fastq + " does not exist!"
		sys.exit(1)

	if (not os.path.exists(infofile)):
		print args.infofile + " does not exist!"
		sys.exit(1)

	if (not os.path.exists(inserts)):
		print args.infofile + " does not exist!"
		sys.exit(1)

	if (not os.path.exists(outdir)):
		print outdir + " does not exist!"
		sys.exit(1)

	# Call AIH
	print ("##################################################"+"\n"+"CALLING AMPLICON INDEL HUNTER"+"\n")
	amplicon_indel_hunter.main(R1fastq, R2fastq, infofile, outdir, cushion,sig_threshold)

	# Call AID
	print ("##################################################"+"\n"+"CALLING AMPLICON INDEL DIAGNOSER"+"\n")
	amplicon_indel_diagnoser.main(R1fastq, R2fastq, infofile, outdir, cushion, inserts,outdir+"/"+os.path.basename(R1fastq)+".R2fastq.indelcalls.significant.txt")

	# Make VCF file of final output
	print ("##################################################"+"\n"+"CONVERT AID OUTPUT TO VCF"+"\n")
	convert_aid_to_vcf.main(outdir+"/"+os.path.basename(R1fastq)+".R2fastq.finalindelstats",outdir+"/"+os.path.basename(R1fastq)+".R2fastq.finalindelstats.vcf")

# Returns the directory of the script
def getScriptPath():
	return os.path.dirname(os.path.realpath(sys.argv[0]))

# Adds two classes to the formatter_class
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
	pass

if __name__ == "__main__":
	version_string="##################\nAmplicon Indel Hunter (aiHunter) v1.1.0"

	description_string="##########################################################\n| Amplicon Indel Hunter is a reference genome-independent |\n| large (>5-bp) indel detection method used for the	   |\n| identification of somatic indels in amplicon-based,	 |\n| paired-end, NGS data. It must be used with professional |\n| care and judgement.					 |\n##########################################################"

	usage_short="aiHunter --read1 <Read1 Fastq. Required.> --read2 <Read2 Fastq. Required.> --amp <Amplicon information sheet. Required.> --inserts <Inserts Fasta file. Required> [options]"

	# RawDescriptionHelpFormatter allows custom formatting of Description text
	parser = argparse.ArgumentParser(formatter_class=CustomFormatter,description=description_string,usage=usage_short)
	parser.add_argument("-v","--version", action="store_true", help="Current version of aiHunter")

	# if user asks for information on info file
	parser.add_argument("-i","--info", action="store_true")

	# Read arguments
	parser.add_argument("--read1",dest="R1fastq", help="Read1 Fastq file path. Required.")
	parser.add_argument("--read2", dest="R2fastq", help="Read2 Fastq file path. Required.")
	parser.add_argument("--amp", dest="infofile",help="Amplicon information file. Use aiHunter --info for more details")
	parser.add_argument("--inserts",dest="inserts",help="Inserts fasta file. Required.")
	parser.add_argument("--out", dest="outdir",default=os.getcwd(),help="Output directory")
	parser.add_argument("--cushion", dest="cushion",type=int,help="Cushion Value",default=5)
	parser.add_argument("--maf", dest="sig_threshold",type=float,help="MAF threshold",default=0.05)
	args = parser.parse_args()

	# if user asks for the version
	if args.version:
		print version_string
		sys.exit(0)

	# Print version always
	print version_string

	if args.info:
		print "\n###########################\n\nThe amplicon info file is a tab-delimited file with the following columns:\nAmpliconName\tPrimer1Sequence\tPrimer2Sequence\tAmpliconLength\tChr\tGenomicAmpliconStart\tGenomicAmpliconEnd\n\nNote: Primer2Sequence should be the reverse complement of the genomic sequence (as used in the assay)\nNote: GenomicAmpliconStart = Primer1 Start\tNote:  GenomicAmpliconEnd = Primer2 End"
		sys.exit(0)

	# if user did not input all required arguments
	if (args.R1fastq == None) or (args.R2fastq == None) or (args.infofile == None) or (args.inserts == None):
		parser.print_help()
		sys.exit(1)

	# Output the input parameters:
	print "\nINPUT PARAMETERS: "
	print "Read1 Fastq file: "+args.R1fastq
	print "Read2 Fastq file: "+args.R2fastq
	print "Amplicon info file: "+args.infofile
	print "Inserts Fasta file: "+args.inserts
	print "Output directory: "+args.outdir
	print "Cushion value: "+str(args.cushion)
	print "MAF threshold: "+str(args.sig_threshold)

	main(args.R1fastq, args.R2fastq, args.infofile, args.inserts,args.outdir,args.cushion,args.sig_threshold)
