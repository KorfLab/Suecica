#!/usr/bin/env python3
# Written by Matthew Porter @ UC Davis, Korf and Comai Labs
import os, gzip

rootdir = "/home/mattchat/SuecicaDupSearch/Data/Sue/SxN_F2s"

for root, subfolders, files in os.walk(rootdir):
    for filename in files:
        normal_reads = 0
        corrupt_reads = 0
        line_num = 0
        line_num_list = []
        
        f_in = os.path.join(root,filename)
        gzipped = False
        if f_in[-7:] == ".sam.gz":
            infile = gzip.open(f_in, 'rb')
            gzipped = True
        elif f_in[-4:] == ".sam":
            infile = open(f_in)
        
        if f_in[-7:] == ".sam.gz" or f_in[-4:] == ".sam": # Make sure we're working with .sam, not .fq
            print("Parsing ", str(f_in), sep='')
            f_out = f_in[:-7] + "_fixed.sam.gz" # Assume .sam.gz first
            if f_in[-4:] == ".sam":
                f_out = f_in[:-4] + "_fixed.sam.gz" # Change filename since original is .sam
        
            with gzip.open(f_out, 'wb') as outfile:
                for line in infile:
                    line_num += 1
                    if gzipped == True: line = str(line, encoding='utf8')
                    if line[:3] != "@SQ":
                        if line[:6] == "SOLEXA":
                            normal_reads += 1
                            outfile.write(bytes(line,"UTF-8"))
                        else:
                            corrupt_reads += 1
                            line_entry = str(line_num) + ": " + str(line) + "\n"
                            line_num_list.append(line_entry)
                    else:
                        outfile.write(bytes(line,"UTF-8"))
            infile.close()
            
            if corrupt_reads == 0:
                try:
                    os.remove(f_out)
                except:
                    print("Could not remove ", str(f_out), sep='')
            else:
                total_reads = normal_reads + corrupt_reads
                corrupt_pct = round(corrupt_reads / total_reads * 100,2)
                f_out2 = f_out[:-7] + "_stats.txt"
                with open(f_out2, 'w') as outfile:
                    outline = str(corrupt_reads) + " out of " + str(total_reads) + " reads (" + str(corrupt_pct) + "%) were corrupted and removed.\n"
                    outfile.write(outline)
                    for entry in line_num_list:
                        outfile.write(entry)
        
        else:
            print("Skipping ", str(f_in), sep='')