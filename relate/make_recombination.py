#! ~/miniconda2/bin/python
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
