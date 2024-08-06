#! ~/miniconda2/bin/python

#simple script to make Relate map recombination file from LDhelmet output, but still needs some manual editing: 
#e.g. for scaffold12 of Trichaptum:
#cat output_Circumboreal_NorAm_NB_scaffold12_LDhat_lk_post_to_text.txt | grep -v "#\|a"  > output_Circumboreal_NorAm_NB_scaffold12_LDhat_lk_post_to_text_no_header.txt
#get CM by multiplying second column with 100/(4*N_e) which is 100/(4*100000)= 0.00025
#cat output_Circumboreal_NorAm_NB_scaffold12_LDhat_lk_post_to_text_no_header.txt | awk '{FS=" "}{print $1}' > testc1
#cat output_Circumboreal_NorAm_NB_scaffold12_LDhat_lk_post_to_text_no_header.txt | awk '{FS=" "}{print $3}' > testc2
#cat testc2 | awk '{print $1*0.00025}' > testc2_r
#cat testc2_r | awk '{print $1/2.381207}' > testc2_r_mb #(length of scaffold12 is 2381207)
#python make_recombination.py testc1_c2_r_mb testc1c2c3
#edit in last line: 2346732 0	        7.6464006440790200
#edit in headers: Position.bp. Rate.cM.Mb. Map.cM.
#convert from tab to space delimitation

import os
import sys
import io

infile = sys.argv[1]
outfile = sys.argv[2]

my_data = io.open(infile, "rU")
#my_data.readline() #to skip first line with headers
my_output = io.open(outfile, "w")

#my_dict={}
pos_list=list()
r_list=list()

for line in my_data:
        line_list = line.strip().split("\t")
        pos = line_list[0]
        r = line_list[1]
        pos_list.append(pos)
        r_list.append(r)
        #my_dict[pos] = r

#print(my_dict)
#print(my_dict.keys())
#print(my_dict.values())

index=0
rdist=0

for i in pos_list:
        #print(i)
        pos=pos_list[index]
        r = r_list[index]
        index = index + 1
        diff =  int(pos_list[index]) - int(pos_list[index-1])
        print(pos + "\t" + r + "\t" + str(diff) + "\t" + str(rdist) + "\n")
        my_output.write(pos + "\t" + r  + "\t" + str(rdist) + "\n")
        rdist_next = (float(r) * int(diff)) + float(rdist)
        rdist = rdist_next
