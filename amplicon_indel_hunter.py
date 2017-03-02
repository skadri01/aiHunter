#!/usr/bin/env python

from Bio.Seq import Seq
import argparse
import os
import sys

###########################
# Amplicon Indel Hunter v0.1.0
###########################
# @ authors: skadri, jsegal

# This program is supposed to take in 2 paired end fastq files from amplicon assays to detect indels 
#It also needs to take in a file of amplicon data in the following format:
# Amplicon_Name  Primer1  qualmod1  Primer2  qualmod2  amplicon_length
# The program should take in the two fastqs and the info and then iterate read by read through both fastq's together. 

# For each pair that matches a primer pair, 
#       - Find the primer set that matches the read pair
#       - find the expected amplicon length from the info sheet
#       - use this length to try a test hybridization between reads at that position
#       - if hyb is good, give read pair an "OK", if hyb is bad give it a "FLAG"
#       - count up total statistics across all primer pairs in the info sheet
# Functionality to handle amplicons shorter than read length

def main(R1fastq, R2fastq, infofile, outdir, cushion, sig_threshold):

	print("*** Starting Indel Hunter ***")
	# Read the fastq files:
        R1reads=[]
        R1qscores=[]
        R2reads=[]
        R2qscores=[]
        
	print("Reading Read1 fastq file")
        R1fastqfile=open(R1fastq, 'rU')
        #put reads and qscores and titles into memory for R1
        counter=0
        for line in R1fastqfile:
                counter=counter+1
                linemod=counter % 4
                if (linemod == 2):
                        R1reads.append(line.strip())
                if (linemod == 0):
                        R1qscores.append(line.strip())

        R1fastqfile.close()     
	print("***Done***")

	print("Reading Read2 fastq file")
        R2fastqfile=open(R2fastq, 'rU')
        #put reads and qscores and titles into memory for R2
        counter=0
        for line in R2fastqfile:
                counter=counter+1
                linemod=counter %4
                if (linemod == 2):
                        R2reads.append(line.strip())
                if (linemod == 0):
                        R2qscores.append(line.strip())

        R2fastqfile.close()
	print("***Done***")

	print("Reading Amplicon Infosheet.")
	# read the info sheet
	# Format is:
	# amplicon_name  primer1  primer2  length chr start end
        inputfile=open(infofile, 'rU')
        header=inputfile.readline()
	#Brings all amplicon info into ampliconinfo and initializes amplicon_indel_data with amplicon_name -> 0
        ampliconinfo={}
        amplicon_indel_data={}
        for lin in inputfile:
                line=lin.strip()
                ampliconinfo[line.split('\t')[0]]=line
                amplicon_indel_data[line.split('\t')[0]]="0:0"
        inputfile.close()
	print("***Done***\n")

	acceptablemismatch = 3
	cnt=0
	useful=0
	# 10/3/15 - write in a quality check
	# IF FILTERED NOT INPUT, THEN WRITE THE QUALITY CHECK IN
	for i in range(0,len(R1reads)):
		cnt=cnt+1
		# Which amplicon does the read pair represent - use read 1 and read 2 first
		ampliconIDhit=searchforamplicon(R1reads[i], R2reads[i], ampliconinfo)
		# if the input to this is filtered reads, then we will not have False
                if not (ampliconIDhit == "False"):
			useful=useful+1
			ampliconData=ampliconinfo[ampliconIDhit]
			# read1 is upstream of read2
			pairedreadresult=hybridize(R2reads[i], R1reads[i], R2qscores[i], R1qscores[i], ampliconinfo[ampliconIDhit], cushion)
			amplicon_indel_data[ampliconIDhit]=resultadder(amplicon_indel_data[ampliconIDhit], pairedreadresult)
		else:
			ampliconIDhit=searchforamplicon(R2reads[i], R1reads[i], ampliconinfo)
			if not(ampliconIDhit == "False"):
				useful=useful+1
				ampliconData=ampliconinfo[ampliconIDhit]
				# read2 is upstream of read1 - like Illumina
				pairedreadresult=hybridize(R1reads[i], R2reads[i], R1qscores[i], R2qscores[i], ampliconinfo[ampliconIDhit], cushion)
				amplicon_indel_data[ampliconIDhit]=resultadder(amplicon_indel_data[ampliconIDhit], pairedreadresult)

	print("Total reads processed = "+str(cnt))
	print("Reads matching primers in infosheet = "+str(useful))

	outputfile=open(outdir + "/" + os.path.basename(R1fastq) + ".R2fastq.indelcalls.txt", 'w')
	# significant file
	sigfile=open(outdir + "/" + os.path.basename(R1fastq) + ".R2fastq.indelcalls.significant.txt", 'w')

        outputfile.write("Amplicon\tReads_w_indel>" + str(cushion) + "bp\t#readpairspassingfilter\t%indel\n")
	results_num=0
        for k,v in amplicon_indel_data.iteritems():
                pairs_with_indel = v.split(":")[0]
                pairs_passing_filter = v.split(":")[1]
		#Each pair represents 2 reads, so depths and pos counts need to be doubled
		depth_okqual = 2*int(pairs_passing_filter)
                reads_withindel = 2*int(pairs_with_indel)
		try:
                        indelpercent=str(float(pairs_with_indel)/float(pairs_passing_filter))
			if float(indelpercent)>sig_threshold:
				sigfile.write(k+'\n')
				results_num+=1
                except:
                        indelpercent="N/A"
		
		outputfile.write(k + '\t' + str(reads_withindel) + '\t' + str(depth_okqual) + '\t' + indelpercent + '\n')

        outputfile.close()
	sigfile.close()

	print("\nNumber of amplicons with indel = "+str(results_num)+"\n")

	print("***Amplicon Indel Hunter finished***")


############################
# FUNCTIONS
############################
# Returns the amplicon that this read pair matches. If no primer pair match the read pair, returns "False"
# Amplicon info will have following information
# AMPLICON_NAME      PRIMER1_SEQUENCE       PRIMER2_SEQUENCE(Reverse complement)
def searchforamplicon(read1, read2, ampliconinfo):
        acceptablemismatch = 3
	for ampl, value in ampliconinfo.iteritems():
                # No reverse complement needed
                primer2=value.split('\t')[2]
                primer1=value.split('\t')[1]
                # if read1 starts with primer 1 and read2 starts with primer 2
                # OR
                # if read2 starts with primer1 and read1 starts with primer 2
                if((loosematch(primer1, read1, acceptablemismatch) and loosematch(primer2, read2, acceptablemismatch))): 
                        return ampl
        return "False"

# This function hybridizes the overlapping region of one read mate with the other mate
# It allows for an cushion region. 
# Returns:
# 0:0 - No hybridization
# 0:1
# 1:1
def hybridize(read1, read2, read1q, read2q,  ampliconData, cushion):
	allowedMismatches = (int)(0.1*len(read1))

        if len(read1) == len(read2):
		ampliconlength=int(ampliconData.split('\t')[3])
		
		# Function assumes that Read1 is always downstream and RC strand - make sure the function call is correct 
		read1rev=Seq(read1).reverse_complement()
		read1qrev=read1q[::-1]
		# If the amplicon is longer than the read length - normal case
		if(ampliconlength>len(read1)):
			offset=ampliconlength-len(read1)
			for i in range(0,int(cushion)+1):
				if((offset-i)>=0):
					if (loosematchq(read2[offset-i:], read1rev,read2q[offset-i:],read1qrev)):
						return "0:1"
                                if (loosematchq(read2[offset+i:], read1rev, read2q[offset+i:],read1qrev)):
                                        return "0:1"
			return "1:1"
		# Amplicon is shorter than read - we are reading adapter sequence
		else:
                        offset=len(read1)-ampliconlength
                        for i in range(0,int(cushion)+1):
				if((offset-i)>=0):
					if (loosematchq(read1rev[offset-i:], read2, read1qrev[offset-i:],read2q)):
						return "0:1"
                                if (loosematchq(read1rev[offset+i:], read2, read1qrev[offset+i:],read2q)):
                                        return "0:1"
                        return "1:1"
			
                                
        print "Hyb FAILURE"
        return "0:0"

#Returns true if the read start matches the querystring with atmost "acceptablemismatch" number of mismatches 
def loosematch(querystring, read, acceptablemismatch):
        mismatch=0
	
        for i in range(0,len(querystring)):
                if not(querystring[i] == read[i]):
                        mismatch=mismatch+1
	if mismatch > acceptablemismatch:
                return False
        else:
                return True

def loosematchq(querystring, read, queryqscores ,readqscores):
        mismatch=0
	okbases=0

	minqual=25
        qualchecker = minqual + 33

        for i in range(0,len(querystring)):
                if(ord(queryqscores[i])>qualchecker and ord(readqscores[i])>qualchecker):
			okbases=okbases+1
			if not(querystring[i] == read[i]):
				mismatch=mismatch+1
	
	# calculate the number of acceptable mismatches
	acceptablemismatch=0.1*(okbases+1)

        if mismatch >= acceptablemismatch:
                return False
        else:
                return True



def resultadder(current, addition):
        currentcomponents=current.split(':')
        additioncomponents=addition.split(':')
        result=str(int(currentcomponents[0]) + int(additioncomponents[0])) + ':' + str(int(currentcomponents[1]) + int(additioncomponents[1]))
        return result


if __name__ == "__main__":
        parser = argparse.ArgumentParser()
        parser.add_argument("-f", dest="R1fastq")
        parser.add_argument("-r", dest="R2fastq")
        parser.add_argument("-o", dest="outdir")
        parser.add_argument("-i", dest="infofile")
        parser.add_argument("-c", dest="cushion")       
	parser.add_argument("-s", dest="sigthresh")
        args = parser.parse_args()

	if (args.R1fastq == None) or (args.R2fastq == None) or (args.infofile == None):
            print 'Failed: Usage amplicon_indel_hunter.py -f <R1fastq> -r <R2fastq> -i <amplicon_infosheet> -o <outdir, def=R1fastq.dir> -c <cushion, def=5>.'
            print
            sys.exit(1)

        if (args.cushion == None):
            cushionvalue=5
        else:
            cushionvalue=args.cushion

	if(args.sigthresh == None):
		args.sigthresh=0.05
	    

        if (args.outdir == None):            
            outdir=os.path.dirname(os.path.realpath(args.R1fastq))
        else:
            outdir=args.outdir

        print outdir + " is entered as outdir"

        if (not os.path.exists(args.R1fastq)):
            print args.R1fastq + " does not exist!"
            sys.exit(1)
        
        if (not os.path.exists(args.R2fastq)):
             print args.R2fastq + " does not exist!"
             sys.exit(1)

        if (not os.path.exists(args.infofile)):
             print args.infofile + " does not exist!"
             sys.exit(1)

        if (not os.path.exists(outdir)):
             print outdir + " does not exist!"
             sys.exit(1)

        main(args.R1fastq, args.R2fastq, args.infofile, outdir, cushionvalue,args.sigthresh)
