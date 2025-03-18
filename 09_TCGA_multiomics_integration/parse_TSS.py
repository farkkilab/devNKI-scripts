#! /usr/bin/python3

import sys
import getopt
import numpy as np


######## Taking arguments ########
try:
    opts, args = getopt.getopt(sys.argv[1:], 'i:o:g:m:h', ['input=', 'output=', 'genes=', 'minlen=', 'help'])
except getopt.GetoptError:
    usage()
    sys.exit(2)

for opt, arg in opts:
    if opt in ('-h', '--help'):
        usage()
        sys.exit(2)
    elif opt in ('-i', '--input'):
        inputf = arg
    elif opt in ('-o', '--output'):
        outputf = arg
    elif opt in ('-g' , '--genes'):
        genesf = arg
    elif opt in ('-m', '--minlen'):
        minlength = arg
    else:
        usage()
        sys.exit(2)

############ Defining functions ##########
def findcolumns (firstline):
    count = 0
    for col in firstline.split():
        if (col == 'Gene_Symbol'):
            genecol = count
        elif (col == 'Position_to_TSS'):
            tsscol = count
        count += 1
    return(genecol, tsscol)


#Reading input files
input_file = open(inputf, "r")
output_file = open(outputf, "w")
genes_file = open(genesf, "r")

print ("Input file:", inputf)
print ("Finding TSS in a position less than:", minlength)

genes_to_find = genes_file.readlines()
fl = input_file.readlines()
#print(fl[0])

target_columns = findcolumns(fl[0])
#print (*target_columns)

#Printing header in ouput file
output_file.write("\"input_gene\"\t\"TSS_distance\"\t" + fl[0])

for line in fl:
    line = line.rstrip()
    bins = line.split()[0]
    genes = line.split()[target_columns[0]]
    pos = line.split()[target_columns[1]]
    genes_listed = genes.split(';')
    TSSites = pos.split(';')
    for j in genes_to_find:
        input_gene = j.rstrip()
        if (input_gene in genes_listed):
            genepos = [i for i, n in enumerate(genes_listed) if n == input_gene]
            for position in genepos:
                TSS = TSSites[position]
                if (abs(int(TSS)) < int(minlength)):
                    output_file.write(input_gene + "\t" + TSS + "\t" + line + "\n")
                    break
output_file.close()
