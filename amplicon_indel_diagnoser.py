#!/usr/bin/env python

from Bio.Seq import Seq
from Bio import SeqIO
import argparse
import os
import sys
import operator

#This program is supposed to take in 2 paired end fastq files from amplicon assays and an amplicon ID that was flagged as containing an indel
#It also needs to take in a file of amplicon info sheet
# AmpliconID  Primer1  QualMod  Primer2 QualMod  AmpliconLength

#The program needs to do a few things:
# 1. it needs to make sure the amplicon ID is actually in the info sheet
# 2. it reads in the reads from R1fastq and R2fastq
# 3. it extracts the upstream oligo, downstream oligo and intervening sequence from the datasheet
# 4. it finds all hopefully mutant reads and sends them to an output file and to a dict file (using hybridize with the "correct" amplicon size)
# 5. it pulls out the "mutant" (non hybridizing) read pair that is present the most frequently in the fastq files.  This is the assumed "consensus mutant pair"
# 6. it slides this mutant pair against each other to find out the "true" offset for this read pair
# 7. from the expected amplicon length and expected offset and the "true" offset, the actual amplicon length can be determined (and thus the indel size)
# 8. we need to take the reference sequence (given in the infosheet file) and compare that with the new sequence that comes from the reads
# 9. If it's an insertion of X bases, we will try matching the reference with X number of N's inserted at every possible position with the new amplicon (looking for highest match %)
# 10. If it's a deletion, we will try matching the reference with X number of bases deleted from every position with the new amplicon (again looking for highest match%)
#     This will give us the best guess of the expected position of the insertion or deletion
#     Knowing this, we can extract the added or deleted sequence either from the read or the reference (depending on whether it's an insertion or deletion
#     Then we can determine the genomic coordinates of the amplicon from the infosheet and put together the final position and base changes and %mutant
#     for the final output.

def main(R1fastq, R2fastq, infofile, outdir, cushionvalue, fastafile,amplicons):

        print "\n***Initializing Amplicon Indel Diagnoser (AID)***\n"
        

	outputfilename=outdir + "/" + os.path.basename(R1fastq) + ".R2fastq.finalindelstats"
	outfile=open(outputfilename,"w")

	R1reads=[]
        R2reads=[]

	print("Reading Read1 file: "+R1fastq)
        R1fastqfile=open(R1fastq, 'r')
        #put reads into memory for R1
        counter=0
        for line in R1fastqfile:
                counter=counter+1
                if (counter%4 == 2):
                        R1reads.append(line.strip())
        R1fastqfile.close()     
	print("***Done***")

	print("Reading Read2 file: "+R2fastq)
        R2fastqfile=open(R2fastq, 'r')
        #put reads (not qscores) into memory for R2
        counter=0
        for line in R2fastqfile:
                counter=counter+1
                if (counter%4 == 2):
                        R2reads.append(line.strip())
        R2fastqfile.close()
	print("***Done***")
        
        print R1fastq + " has " + str(len(R1reads)) + " reads."
        print R2fastq + " has " + str(len(R2reads)) + " reads."

	# Read the list
	ampliconIDs=[]
	ampFile=open(amplicons,'rU')
	for lin in ampFile:
		ampliconIDs.append(lin.strip())

	# read amplicon data in dictionary ONLY for the IDs in the IDs list
        inputfile=open(infofile, 'rU')
        header=inputfile.readline()
	amplicon_info={}
	for line in inputfile:
		id=line.split('\t')[0]
		if(id in ampliconIDs):
			amplicon_info[id]=line.strip()

	inputfile.close()
	print(str(len(ampliconIDs))+" in significant list")
	print(str(len(amplicon_info))+" found in info sheet")

	if(len(ampliconIDs)!=len(amplicon_info)):
		print("ERROR: Significant amplicon list has amplicons not provided in "+ampFile)

        #Search through the amplicon infosheet for an entry that matches the ampliconID given.  When found, extract the ULSO and DLSO and the amplicon length.
        ampliconfound=False
	for ampliconID,data in amplicon_info.iteritems():
                print "*** Testing for *** "+ampliconID
		short_amplicon=False;
		#The info file lists a primer1 and primer2 
		# Primer1 is read in the correct phase by the FIRST/SECOND read
                # Primer2 is already RC of reference and read by the SECOND/FIRST read
                #This program is written to (for each line in the amplicon data), try matching primers with read1/read2 and then read2/read1
#                print(" In amplicon info :"+ampliconID+" "+lin); 
		primer1=data.split('\t')[1]
                primer2=data.split('\t')[2]
		ampliconlength=int(data.split('\t')[3])
		
		print "Length of reads = "+str(len(R1reads[0]))
		
		if(ampliconlength<len(R1reads[0])):
			short_amplicon=True
			print "#######################################"
			print "AMPLICON LENGTH IS LESS THAN READ LENGTH"

		print(" Reading from Info Sheet: Amplicon length = "+str(ampliconlength))
		# TODO: Add functionality for when amplicon is shorter than read length
                normaloffset=abs(ampliconlength-len(R1reads[0]))
		print "NORMAL OFFSET IS :"+str(normaloffset)
                chrom=data.split('\t')[4]
                chromstartpos=data.split('\t')[5]
                chromendpos=data.split('\t')[6]
        
		ampliconFound=False
		handle=open(fastafile,"rU")
		for record in SeqIO.parse(handle,"fasta"):
			if(record.id==ampliconID):
				# Read the fasta file and get the middle sequence:
				#The middle intervening sequencing is in sense phase - that is same as primer1 and the read primer 1 overlaps                                                           
                                #The usable reference sequence to compare with the read that overlaps with primer 1 and revcom of the read that overlaps with primer2
				ampliconReferenceSequence = primer1 + record.seq + Seq(primer2).reverse_complement()
				print "Amplicon Reference Sequence: "+ampliconReferenceSequence
#				print "Length of reference "+str(len(ampliconReferenceSequence))
				ampliconFound=True
				break;
		if not ampliconFound:
			print "No sequence found for "+ampliconID
			break;
		handle.close()

		slideresults={}

		# HERE READ HAVING PRIMER 1 IS STORED IN R1READS_W_PRIMERMATCH AND SIMILARLY FOR PRIMER2...
		R1reads_w_primermatch=[]
		R2reads_w_primermatch=[]
		reads_w_primermatch=[]
	        #Now go through the data read pair by read pair, and if they match the amplicon ID...
		for i in range(0,len(R1reads)):
			#First find out if the reads are from the correct amplicon      
			#then if the reads are both from the right amplicon
		        # READ1 IS UPSTREAM OF READ2
			if loosematch(primer2, R2reads[i], 2) and loosematch(primer1, R1reads[i], 2):
				R1reads_w_primermatch.append(R1reads[i])
				R2reads_w_primermatch.append(R2reads[i])
				reads_w_primermatch.append(i)
				# the read that starts downstream is the first argument
				# Check for amplicon length to determine which one this is
				if(short_amplicon):
#					slidevalue=slide(Seq(R2reads[i]).reverse_complement(), R1reads[i])
					slidevalue=slide(R1reads[i],Seq(R2reads[i]).reverse_complement())
				else:
					slidevalue=slide(Seq(R2reads[i]).reverse_complement(), R1reads[i])
				if slidevalue in slideresults:
					slideresults[slidevalue]=slideresults[slidevalue]+1
				else:
					slideresults[slidevalue]=1
			else:
				# READ 2 IS UPSTREAM OF READ1
				if loosematch(primer1, R2reads[i], 2) and loosematch(primer2, R1reads[i], 2):
					R1reads_w_primermatch.append(R2reads[i])
					R2reads_w_primermatch.append(R1reads[i])
					reads_w_primermatch.append(i)
					if(short_amplicon):
						slidevalue=slide(R2reads[i],Seq(R1reads[i]).reverse_complement())
					else:
						slidevalue=slide(Seq(R1reads[i]).reverse_complement(), R2reads[i])
					if slidevalue in slideresults:
						slideresults[slidevalue]=slideresults[slidevalue]+1
					else:
						slideresults[slidevalue]=1

		print "Number of reads with primermatch " + ampliconID +" is " + str(len(reads_w_primermatch))
		print "Offset_size\tCount"
        
		for k,v in slideresults.iteritems():
			print str(k) + '\t' + str(v)    

   	        #get number with the correct offset (is this necessary?) and then erase that from the slideresults so we can get the next highest value
                #also erase the entries within the cushion range of the normal offset.  I believe this is necessary because in an amplicon
                #that is baseline bad like NPM1, if there is a low % NPM1 mutation, the indel returned would actually be the one-offset
		#instead of the true mutant 4-offset because of the homopolymer in that amplicon.
		try:
			NumReadPairswithCorrectOffset=slideresults[normaloffset]
		except: 
			NumReadPairswithCorrectOffset=0

	        #Here we purge the results that are within 5 bases of the normal offset
		slideresults[normaloffset]=0
		for i in range(1,(cushionvalue + 1)):
			slideresults[(normaloffset+i)]=0
			slideresults[(normaloffset-i)]=0
        
		print "Slide results after purging results withint 5bp of normal offset:"
		for k,v in slideresults.iteritems():
                        print str(k) + '\t' + str(v)

		mutantoffset=max(slideresults.iteritems(), key=operator.itemgetter(1))[0]
		NumReadPairswithMaxMutantIndel=slideresults[mutantoffset]

		try:
			mutantoffset=int(mutantoffset)
		except:
                        #converting mutantoffset to int could fail if the most common mutant offset is actually "Fail", meaning the reads didn't match.
			#This would happen if there is an extreme insertion, or if the caller is called by accident on an ok amplicon.
			print "Mutantoffset can not be converted to int because mutantoffset = " + str(mutantoffset)
			final_variant_outputline="Failed indel identification - manual review required for " + ampliconID + "!"
			outfile.write(final_variant_outputline + '\n')
			continue;

		print "The mutant offset is " + str(mutantoffset)
		print "The normal offset is " + str(normaloffset)
		# Mutant offset is the actual offset by which the reads align.
		# If there is a deletion, the mutant offset will be smaller than the normal offset
		# If there is an insertion, the mutant offset will be larger than the normal offset
		if(short_amplicon):
			indelsize=normaloffset-mutantoffset
		else:
			indelsize=mutantoffset-normaloffset
		print "The mutant indel size is " + str(indelsize)

                #Here we go back to the reads that had a primer match and keep all the ones that have the mutant indel size in a new pair of lists                                                       
		reads_w_primermatchANDmutantindel=[]
		mutantampliconlength=ampliconlength + indelsize
                print "Mutant amplicon length = "+str(mutantampliconlength)
		# Reads that match with the mutant amplicon length 
		for i in range(0,len(R1reads_w_primermatch)):
			# NOTE THAT R1READS_W_PRIMERMATCH HAVE THE READS WITH PRIMER1 AND R2... HAS THAT WITH PRIMER2 - SO R2 HAVE TO BE RC - OPP. OF ONCOSCREEN PANELS
			# hybridize function takes into account the shorter amplicon functionality
			mutantindelsizeflag=hybridize(R2reads_w_primermatch[i],R1reads_w_primermatch[i], mutantampliconlength, 0)
			if mutantindelsizeflag:
				reads_w_primermatchANDmutantindel.append(i)
				
	        print "Correct Reference is..."
		print ampliconReferenceSequence
		print "Length of correctphasereference is... " + str(len(ampliconReferenceSequence))
		print "ampliconlength is... " + str(ampliconlength)
		print "Found "+str(len(reads_w_primermatchANDmutantindel))+ " reads with mutant indel"
		print "Mutant offset is "+str(mutantoffset)
		print "Mutant amplicon length = "+str(mutantampliconlength)

		#print("Len(reads with indel:"+str(len(reads_w_primermatchANDmutantindel)))
	        
		correctphaseamplicons_for_consensus={}
		FindIndelResults={}
		for i in range(0,len(reads_w_primermatchANDmutantindel)):
			index=reads_w_primermatchANDmutantindel[i]
			#Iterating through the mutant indel reads, we need to assemble a full amplicon (read 2 up to the mutant offset plus the reverse of read 1)
			if mutantoffset > 0:
				if(mutantampliconlength<len(R1reads_w_primermatch[0])):
					mutantamplicon=R1reads_w_primermatch[index][0:mutantampliconlength ]
				else:
					mutantamplicon=R1reads_w_primermatch[index][0:mutantoffset] + Seq(R2reads_w_primermatch[index]).reverse_complement()
			# Mutation has resulted in the amplicon length being smaller than the read length
			else:
				if(mutantampliconlength<len(R1reads_w_primermatch[0])):
					mutantamplicon=R1reads_w_primermatch[index][0:mutantampliconlength ]
				else:
					mutantamplicon=R1reads_w_primermatch[index]+Seq(R2reads_w_primermatch[index]).reverse_complement()[mutantoffset:len(Seq(R2reads_w_primermatch[index]))]
			
			correctphasemutantamplicon=mutantamplicon

	                #Add it to a dictionary of amplicons that have the correct mutant indel in order to find the most common indel pair sequence (amplicon)                                          
			try:
				correctphaseamplicons_for_consensus[correctphasemutantamplicon]=correctphaseamplicons_for_consensus[correctphasemutantamplicon] + 1
			except:
				correctphaseamplicons_for_consensus[correctphasemutantamplicon]=1

			IndelResult=FindIndel(ampliconReferenceSequence, correctphasemutantamplicon, indelsize)
			try:
				FindIndelResults[IndelResult]=FindIndelResults[IndelResult] + 1
			except:
				FindIndelResults[IndelResult]=1

		for k,v, in FindIndelResults.iteritems():
			print str(k) + " - this result is seen " + str(v) + " times."

		consensusindelposition=max(FindIndelResults.iteritems(), key=operator.itemgetter(1))[0]
		print "The consensus indel position within the reference segment is " + str(consensusindelposition)
		genomicIndelPosition=int(chromstartpos) + int(consensusindelposition) -1
		print "The genomic coordinate of the indel is " + str(chrom) + ":" + str(genomicIndelPosition)
                #This will get us the most commonly observed mutant amplicon with the correct mutant indel size                                                                                                                                    
		consensuscorrectphasemutantamplicon=max(correctphaseamplicons_for_consensus.iteritems(), key=operator.itemgetter(1))[0]
		print "The consensus correct phase mutant amplicon detected is... "
		print consensuscorrectphasemutantamplicon
		if indelsize>0:
		      indelstring="+" + str(indelsize) + consensuscorrectphasemutantamplicon[consensusindelposition:(consensusindelposition + indelsize)]
		else:
		      indelstring=str(indelsize) + ampliconReferenceSequence[consensusindelposition:(consensusindelposition - indelsize)]
		refbaseatindelstart=ampliconReferenceSequence[int(consensusindelposition)-1]
		print "Number of read pairs with primer match and mutant indel is " + str(len(reads_w_primermatchANDmutantindel))
		print "Number of read pairs with primer match is " + str(len(R1reads_w_primermatch))
		print("Len(reads with indel:"+str(len(reads_w_primermatchANDmutantindel)))
		print("Denominator = "+str(len(R1reads_w_primermatch)))
		mutantreadpercentage=float(len(reads_w_primermatchANDmutantindel))/float(len(R1reads_w_primermatch))
		referencereadpercentage=1-mutantreadpercentage
		depth = 2*len(R1reads_w_primermatch)

                #Data Output                                                                                                                                                                                     
		final_variant_outputline=chrom + '\t' + str(genomicIndelPosition) + '\t' + str(depth) + '\t' + str(depth) + '\t' + str(refbaseatindelstart) + '\t' + str(referencereadpercentage) + '\t' + str(indelstring) + '\t' + str(mutantreadpercentage) + '\n'
		#If the indel is "detected" at a location within the primer sequences, then it is likely garbage.  Output should reflect that.                                                                                    
		indellocationinprimer=indellocationchecker(ampliconReferenceSequence, primer1, primer2, consensusindelposition, indelsize)
                #This seems to be messing up for FLT3 ITDs, so if the primer is in the duplicated region it is failing.  Removing it now. 5-1-14.                                                                                         
		print final_variant_outputline

		outfile.write(final_variant_outputline)
		      
		# Remove reads that belong to these primers                                                                                                                                      
		for i in sorted(reads_w_primermatch, reverse=True):
			x=R1reads.pop(i)
			x=R2reads.pop(i)
		print "Remaining reads =  " + str(len(R1reads))
	outfile.close()
 
##########################################
# FUNCTIONS
##########################################
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


def slide(read1rev, read2):
#        slidedata={}
#        slidedata[0]=loosematch(read1rev,read2,len(read1rev)/10)
        for i in range(0,len(read1rev)-20):
                slidepos=loosematch(read1rev[0:len(read1rev)-i], read2[i:], len(read2[i:])/20)
                if slidepos:
                        return i
                slideneg=loosematch(read1rev[i:], read2[0:len(read2)-i], len(read1rev[i:])/20)
                if slideneg:
                        returnvalue = -i
                        return returnvalue
        return "Fail"   

def indellocationchecker(correctphasereference, primerforread2, primerforread1, consensusindelposition, indelsize):
        read2primerlength=len(primerforread2)
        read1primerlength=len(primerforread1)
        lengthofreference=len(correctphasereference)
        if indelsize>=0:
                furthestreach=consensusindelposition
        else:
                furthestreach=consensusindelposition - indelsize
		if consensusindelposition <= read2primerlength:
			return True
	if furthestreach <= (lengthofreference - read1primerlength):
                return True
        return False

def FindIndel(reference, mutant, indelsize):
        #Here we take in a reference sequence and a mutant (indel) sequence, and an indel size which can be positive or negative
	#print("FindIndel for reference:"+reference+" vs mutant: "+mutant+" and indel size: "+str(indelsize))
        matchdatabyposition={}
        if int(indelsize)>0:
                insertstring="X"*int(indelsize)
                for i in range(0, len(reference)+1):
                        testreference=reference[0:i] + insertstring + reference[i:]
                        matchdatabyposition[i]=matchcounter(testreference, mutant)
        else:
                for i in range(0, len(reference)+int(indelsize)+1):
                        testreference=reference[0:i] + reference[i-int(indelsize):]
                        matchdatabyposition[i]=matchcounter(testreference, mutant)
        import operator
        bestmatchposition=max(matchdatabyposition.iteritems(), key=operator.itemgetter(1))[0]
        return bestmatchposition


def matchcounter(sequenceA, sequenceB):
	if not (len(sequenceA) == len(sequenceB)):
                print "There seems to be a problem with the comparison lengths."
		print(sequenceA)
		print(sequenceB)
        matchcount=0
        for i in range(0,len(sequenceA)):
                if sequenceA[i]==sequenceB[i]:
                        matchcount=matchcount+1
        return matchcount

# This function hybridizes the overlapping region of one read mate with the other mate
# It allows for an cushion region. 
# Returns:
# 0:0 - No hybridization
# 0:1
# 1:1
def hybridize(read1, read2, ampliconlength, cushion):
        allowedMismatches = (int)(0.1*len(read1))

        if len(read1) == len(read2):                
                # Function assumes that Read1 is always downstream and RC strand - make sure the function call is correct 
                read1rev=Seq(read1).reverse_complement()
                # If the amplicon is longer than the read length - normal case
                if(ampliconlength>len(read1)):
                        offset=ampliconlength-len(read1)
                        for i in range(0,int(cushion)+1):
                                if (loosematch(read2[offset-i:], read1rev, allowedMismatches)):
                                        return True
                                if (loosematch(read2[offset+i:], read1rev, allowedMismatches)):
                                        return True                            
                        return False
                # Amplicon is shorter than read - we are reading adapter sequence
                else:
                        offset=len(read1)-ampliconlength
                        for i in range(0,int(cushion)+1):
                                if (loosematch(read1rev[offset-i:], read2, allowedMismatches)):
                                        return True
                                if (loosematch(read1rev[offset+i:], read2, allowedMismatches)):
                                        return True
                        return False
                                       
        print "Hyb FAILURE"
        return False


if __name__ == "__main__":
        parser = argparse.ArgumentParser()
        parser.add_argument("-f", dest="R1fastq")
        parser.add_argument("-r", dest="R2fastq")
        parser.add_argument("-i", dest="infofile")
        parser.add_argument("-o", dest="outdir")
        parser.add_argument("-c", dest="cushion")
	parser.add_argument("-s", dest="fastafile")
	parser.add_argument("-a", dest="amplicons")
        args = parser.parse_args()

        if (args.R1fastq == None) or (args.R2fastq == None) or (args.infofile == None)  or (args.fastafile == None) or (args.amplicons==None):
                print 'Failed: Usage amplicon_indel_diagnoser.py -f <R1fastq> -r <R2fastq> -i <amplicon_infosheet> -o <outdir, def=R1fastqdir> -c <cushionvalue_def=5> -s <ampliconfastafile> -a <amplicons_to_test_file>.'
                print
                sys.exit(1)
        
        if args.cushion == None:
                cushionvalue=5
        else:
                try:
                        cushionvalue=int(args.cushion)
                except:
                        print "Cushion value needs to be int."
                        sys.exit(1)

        if (not os.path.exists(args.R1fastq)):
                print args.R1fastq + " does not exist!"
                sys.exit(1)
        
        if (not os.path.exists(args.R2fastq)):
                print args.R2fastq + " does not exist!"
                sys.exit(1)

        if (not os.path.exists(args.infofile)):
                print args.infofile + " does not exist!"
                sys.exit(1)
        
        if (not os.path.exists(args.outdir)):
                print args.outdir + " does not exist!"
                sys.exit(1)

	if (not os.path.exists(args.fastafile)):
		print args.fastafile + "does not exist!"
		sys.exit(1)

	if (not os.path.exists(args.amplicons)):
		print args.amplicons + "does not exist!"
		sys.exit(1)

        main(args.R1fastq, args.R2fastq, args.infofile, args.outdir, cushionvalue,args.fastafile,args.amplicons)
