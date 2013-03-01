#!/usr/bin/env python3
import sys; sys.path.append('../Packages/')
from GenomeTools import Pileup

#filename = "/home/mattchat/SuecicaDupSearch/Data/Sue/sue1/single_bp_preprocess/set5/libSUE1_set5_aln_Aa_Filtered.pileup.gz"
#Pileup.coverage(filename, "AtChr1")
#Pileup.coverage(filename, "AtChr2")
#Pileup.coverage(filename, "AtChr3")
#Pileup.coverage(filename, "AtChr4")
#Pileup.coverage(filename, "AtChr5")
#Pileup.coverage(filename, "AtChr1;AtChr2;AtChr3;AtChr4;AtChr5")

#filename = "/home/mattchat/SuecicaDupSearch/Data/weigel_col-0/SRR013327,SRR013328_aln.pileup.gz"
#filename = "/home/mattchat/SuecicaDupSearch/Data/weigel_bur-0/SRR013331,SRR013333_aln.pileup.gz"
#filename = "/home/mattchat/SuecicaDupSearch/Data/weigel_tsu-1/SRR013335,SRR013337_aln.pileup.gz"
#filename = "/home/mattchat/SuecicaDupSearch/Data/C24/C24_all_reads_aln.pileup.gz"
filename = "/home/mattchat/SuecicaDupSearch/Data/Ler-0/ERR031544,SRR279136_TAIR9_aln.pileup.gz"
Pileup.coverage(filename, "Chr1")
Pileup.coverage(filename, "Chr2")
Pileup.coverage(filename, "Chr3")
Pileup.coverage(filename, "Chr4")
Pileup.coverage(filename, "Chr5")
Pileup.coverage(filename, "Chr1;Chr2;Chr3;Chr4;Chr5")