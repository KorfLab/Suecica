#!/usr/bin/env python3.2
import re
import sys; sys.path.append('/home/mattchat/Packages/')
from GenomeTools import Pileup, GFF, Viterbi

def Output(cur_chr, obs, path, Suecica_posdict):
    pattern = r"D{5,}"
    recomp = re.compile(pattern)
    
    outfilename = "Output/" + str(cur_chr) + "_ViterbiObs.txt"
    outobs = re.sub(" ", "", obs)
    with open(outfilename, 'w') as outfile:
        outfile.write(outobs)
    
    outfilename = "Output/" + str(cur_chr) + "_ViterbiPath.txt"
    with open(outfilename, 'w') as outfile:
        outfile.write(path)
    
    outfilename = "Output/" + str(cur_chr) + "_Duplications.txt"
    with open(outfilename, 'w') as outfile:
        for m in recomp.finditer(path):
            s_chr, s_pos = Suecica_posdict[max(m.start()-10,0)]
            e_chr, e_pos = Suecica_posdict[min(m.end()+10,len(path)-1)]
            line = str(s_chr) + ", " + str(s_pos) + "\t" + path[m.start()-10:m.end()+10] + "\t" + str(e_chr) + ", " + str(e_pos) + "\n\n"
            outfile.write(line)
    

filename = "Data/SueAll/Sue-OurRESCAN/libsSue/ThalyrataAlignment/cleaned-RESCANlibsSue_sxn_pileup_simplified.txt"
Suecica = Pileup.Parse(filename, True, "AtChr3:3000000-4000000")
#Suecica = Pileup.Parse(filename, True, "AtChr3")#

cur_chr = ""
obs = ""
filename = "SuecicaViterbiParam - At.txt"
#filename = "SuecicaViterbiParam - Al.txt"
Suecica.pileup_dict[("END",0)] = [0,0]
Suecica_posdict = {}
i = -1

for pos_info in Suecica:
    (chr, pos), [ref, nuclist] = pos_info
    i += 1
    Suecica_posdict[i] = (chr, pos)
    if cur_chr == "": cur_chr = chr
    if cur_chr != chr:
        DupRegions = Viterbi.Viterbi_0thOrder(filename, obs.split(" "))
        (state, path, prob) = DupRegions.GetOptimalPath()
        path = ''.join(path)
        Output(cur_chr, obs, path, Suecica_posdict)
        cur_chr = chr
        obs = ""
        i = -1
    if chr == "END": break
    region = str(chr) + ":" + str(pos) + "-" + str(pos)
    cov = Suecica.Coverage(region)
    if cov <= 20: obs += str(cov) + " "