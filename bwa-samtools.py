#! /usr/bin/env python

import os, sys, math

#usage: bwa-samtools.py [-q trim quality] [-pb path to bwa] [-ps path to samtools] database_file.fa reads.fq 

mode = 0
index = ""
trimQual = "20"
pathBWA = "/share/apps/bwa/"
pathSam = "/share/apps/samtools/"
if "-pb" in sys.argv:
   pathBWA = sys.argv[sys.argv.index("-pb")+1]
   mode+=2
if "-ps" in sys.argv:
   pathSam = sys.argv[sys.argv.index("-ps")+1]
   mode+=2
if "-q" in sys.argv:
   trimQual = sys.argv[sys.argv.index("-q")+1]
   mode+=2

database = sys.argv[mode+1]
file = sys.argv[mode+2]
name = file.split('.')[0]

dsize = os.path.getsize(database)
if dsize < 1048576000:
   index = " is "
else:
   index = " bwtsw "   
   
os.system(pathBWA+"bwa index -a"+index+" "+database)
os.system(pathBWA+"bwa aln -t 1 -q "+trimQual+" "+database+" "+file+" > "+name+"_aln.sai")
os.system(pathBWA+"bwa samse "+database+" "+name+"_aln.sai "+file+" > "+name+"_aln.sam")
os.system(pathSam+"samtools view -bS "+name+"_aln.sam"+" > "+name+"_aln.bam")
os.system(pathSam+"samtools sort "+name+"_aln.bam"+" "+name+"_aln.sorted")
#New line below
#samtools index <aln.bam>
#samtools idxstats <aln.bam> 
os.system(pathSam+"samtools index "+name+"_aln.sorted.bam")
os.system(pathSam+"samtools idxstats "+name+"_aln.sorted.bam")
os.system(pathSam+"samtools pileup -f "+database+" -s "+name+"_aln.sorted.bam > "+name+"_pileup.txt")