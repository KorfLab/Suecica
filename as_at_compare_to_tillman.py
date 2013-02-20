#!/usr/bin/env python
import re, time; from random import sample; from scipy.stats import fisher_exact

def get_genes():
    filename = "Data/t9_sorted.gff"
    #Chr1	TAIR9	gene	3631	5899	.	+	.	ID=AT1G01010;Note=protein_coding_gene;Name=AT1G01010
    #Note: Do not include mitochondrial or chloroplast genes
    pattern = r"^Chr([1-5])\tTAIR9\tgene\t(\d*)\t(\d*)\t\.\t[+-]\t\.\tID=(AT.{1}G\d{5});.*$"
    recomp = re.compile(pattern)
    gene_dict = {}
    with open(filename) as gene_file:
        linelist = gene_file.readlines()
        for i in range(0,len(linelist)):
            match = recomp.match(linelist[i])
            if match:
                chr_num = str(match.group(1)); start_pos = int(match.group(2)); end_pos = int(match.group(3)); gene_name = str(match.group(4));
                gene_dict[gene_name] = [chr_num, start_pos, end_pos]
    return (gene_dict)

start = time.time()

gene_dict = get_genes()

mydata = "Data/Sue_SUN_Genes.txt"                    # Compare my top Sue candidate duplicates based on SUN counts to Tillman's June 2011 top duplicate candidate list (probabilities of duplication based on read counts)
#mydata = "Data/AsAt_ratios (Tillman).txt"           # Compare my top Sue candidate duplicates based on Tillman's May 2010 read counts to Tillman's June 2011 top duplicate candidate list (probabilities of duplication based on read counts)
#mydata = "Data/AsAt_ratios (IsabelleParents).txt"   # Compare my top Sue candidate duplicates based on Isabelle's F2 parents' read counts to Tillman's June 2011 top duplicate candidate list (probabilities of duplication based on read counts)
tillmandata = "Data/Tillman_significant_genes_June_2011/Significant_Thaliana_Genes_061511.txt"

#Gene	Chr	Start_pos	End_pos	SUN_Count
#AT1G50030	1	18522441	18539995	23
pattern = r"^(AT[1-5]G\d{5})\t\d{1}\t\d*\t\d*\t(\d*)$"
recomp = re.compile(pattern)
if mydata[0:16] == "Data/AsAt_ratios":
    #Gene	Chr	Start_pos	As_norm-At_norm	As_avg_norm	At_avg_norm	As_avg	At_avg	Col_avg	Tsu_avg	Bur_avg
    #AT1G40104	1	15081952	96.2686	101.0534	4.7848	4481.5	344.8333	231.5	527.0	276.0
    pattern = r"^(AT\d{1}G\d{5})\t\d{1}\t\d*\t([-]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?).*$"
    recomp = re.compile(pattern)
#Locus	Bur1	Bur2	Col1	Col2	Col3	Col4	Col5	Col6	Col7	Tsu1	Tsu2	Sue1	Sue2	Sue2_5	Sue4	Over.Dispersed.Indicator	Pvalue	AdjustPvalue	LogFoldChange
#AT4G10780	150	113	38	105	93	115	85	50	37	249	87	0	0	2	0	NotOverDispersed	2.49E-95	1.78E-91	-4.828965497
pattern2 = r"^(AT[1-5]G\d{5})\t\d*\t\d*\t\d*\t\d*\t\d*\t\d*\t\d*\t\d*\t\d*\t\d*\t\d*\t\d*\t\d*\t\d*\t\d*\t(?:Not)?OverDispersed\t[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?\t[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?\t([-]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?).*$"
recomp2 = re.compile(pattern2)

my_list = []; till_list = []

with open(mydata) as infile:
    for line in infile:
        match = recomp.match(line)
        if match:
            gene_name = str(match.group(1))
            if mydata == "Data/Sue_SUN_Genes.txt":
                SUN_counts = int(match.group(2))
                if SUN_counts > 2:
                    my_list.append(gene_name)
            else:
                asat_norm = float(match.group(2))
                if asat_norm >= 1.96:
                    my_list.append(gene_name)

my_dups = len(my_list)

with open(tillmandata) as infile:
    linelist = infile.readlines()
    for line in linelist:
        match = recomp2.match(line)
        if match:
            gene_name = str(match.group(1))
            foldchange = str(match.group(2))
            if foldchange[0:1] != '-':
                till_list.append(gene_name)
            else:
                break
till_dups = len(till_list)

common_num = len([gene for gene in my_list if gene in till_list])
myprob = common_num / my_dups
print("My duplicate list length: ", my_dups, "\nTillman duplicate list length: ", till_dups, "\nDuplicates in common: ", common_num, "\nProbability of overlap: ", round(myprob,2),sep="")

gene_list = list(gene_dict.keys())
fail = 0
repeat = 10000
for i in range (0,repeat):
    myrandlist = sample(gene_list,my_dups)
    common_num_t = len([gene for gene in myrandlist if gene in till_list])
    #				Random		Actual
    #   	Common	  a			  b
    #   NotInCommon	  c			  d
    oddsratio, pvalue = fisher_exact([[common_num_t, common_num], [my_dups-common_num_t, my_dups-common_num]])
    if pvalue > 0.05 or pvalue <= 0.05 and myprob < common_num_t / my_dups:
        fail += 1
    #print(common_num_t,"\t",common_num,"\n",my_dups-common_num_t,"\t",my_dups-common_num,"\n",pvalue,"\n\n",sep="")

print("\n", fail," out of ", repeat, " (", int(round(fail/repeat,2)*100), "%) trials showed that my data has a probability of overlap with Tillman's data that is indistinguishable from overlap of a comparably sized list of randomly selected genes", sep="")

total = (time.time() - start) / 60
print("\nComparison of Matt and Tillman's As-At gene duplication data completed (" + str(round(total,2)) + " minutes)")