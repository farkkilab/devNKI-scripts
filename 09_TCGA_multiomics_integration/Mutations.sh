#! /bin/bash

#This script will take importan columns from the file mc3.v0.2.8.PUBLIC.maf, ovarian samples and mutations in genes of interest
#This file was downloaded from the GDC portal in this link: https://gdc.cancer.gov/about-data/publications/pancanatlas

project_folder="/mnt/d/users/fperez/NKI_TMAs_AF/TCGA_methylation_mutations_loh/"
INPUT_MC3="Input/Somatic_mutations/mc3.v0.2.8.PUBLIC.maf"
OVAsamples="Input/clinical/TCGA-CDR-SupplementalTableS1_OVsamples.txt"
genes="Input/MHCII_genes.txt"
outputsuffix="PUBLIC_OVA-MHCIIgenes" #There will be as output a .maf file with all columns and a csv file with only columns of interest
#The outputprefix will be "mc3.v0.2.8."

nameinput=`basename ${INPUT_MC3}`

cd $project_folder

#Creating outputdirectory
mkdir -p Intermediate_files/Somatic_mutations/

### Select only OVA samples
head -n 1 \
	${INPUT_MC3} \
	> Intermediate_files/Somatic_mutations/${nameinput}

grep -f \
	${OVAsamples} \
	${INPUT_MC3} >>Intermediate_files/Somatic_mutations/${nameinput}

#List of the samples
cut -f 16 Intermediate_files/Somatic_mutations/${nameinput}  \
	| cut -d "-" -f 1-5 \
	| sort -u > Intermediate_files/Somatic_mutations/mc3.v0.2.8.PUBLIC_OVA_SAMPLES


### Select only HR genes from the OVA samples mutations
head -n 1 Intermediate_files/Somatic_mutations/${nameinput} \
	> Intermediate_files/Somatic_mutations/mc3.v0.2.8.${outputsuffix}".maf"


for gene in `cat ${genes}`
 do 
	 grep -P "^${gene}\t" Intermediate_files/Somatic_mutations/${nameinput}
 done >> Intermediate_files/Somatic_mutations/mc3.v0.2.8.${outputsuffix}".maf"


awk '{if ($11 != "-" && $47 != "-")  print "\""$16"\","$1",chr"$5","$6","$11","$47","$14","$10","$9","$18","$19","$72","$73","$93","$41","$42",chr"$5":"$6"-"$11"-"$47} 
     {if ($11 == "-")  print "\""$16"\","$1",chr"$5","$6","$11","$47","$14","$10","$9","$18","$19","$72","$73","$93","$41","$42",chr"$5":"$6"--"$47} 
     {if ($47=="-") print "\""$16"\","$1",chr"$5","$6","$11","$47","$14","$10","$9","$18","$19","$72","$73","$93","$41","$42",chr"$5":"$6"-"$11"-"}' \
	     Intermediate_files/Somatic_mutations/"mc3.v0.2.8."${outputsuffix}".maf" \
	| grep -v "Intron" \
	| grep -v "UTR" \
	| grep -v "Silent" \
	> Intermediate_files/Somatic_mutations/"mc3.v0.2.8."${outputsuffix}".csv" \
	&& echo "File created: intermediate_files/Somatic_mutations/mc3.v0.2.8.${outputsuffix}.csv"

# | sed -E "s/\-[0-9]+[A-Z]\-[0-9]+[A-Z]+\-[A-Z]?[0-9]+[A-Z]?\-[0-9]+\",/\",/" \ 
