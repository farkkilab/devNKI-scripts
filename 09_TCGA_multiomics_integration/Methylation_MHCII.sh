#! /bin/bash

#fperez
#This script will parse the methylation files downlades from GDC wiht betavalues
#Will merge them and and then using the python script will generate a matrix of betavalues per genes of interest

results_folder="/mnt/d/users/fperez/NKI_TMAs_AF/TCGA_methylation_mutations_loh/Intermediate_files/Methylation/"
methylation_folder="/mnt/d/users/fperez/NKI_TMAs_AF/TCGA_methylation_mutations_loh/input/Methylation/HumanMethylation27/"
basename=`basename ${methylation_folder}`
samplesheet="/mnt/d/users/fperez/NKI_TMAs_AF/TCGA_methylation_mutations_loh/input/Methylation/HumanMethylation27/gdc_sample_sheet.2020-03-27.tsv"
genesIDS="/mnt/d/users/fperez/NKI_TMAs_AF//TCGA_methylation_mutations_loh/input/MHCII_genes.txt"

mkdir -p $results_folder

####  Selection of ovarian samples, ignore normal ###
random_file=`ls ${methylation_folder}/*/*_hg38.txt | head -n 1`

cut -f 1,3-  ${random_file} \
	| grep -v "^Composite" > ${results_folder}"/"${basename}_methylation-annotation.txt

#Produce files with the beta_values
for file in `ls ${methylation_folder}/*/*_hg38.txt`
 do
	cut -f 2 ${file} \
		| grep -v "^Beta" > ${file}"_betavalue"
done

#Generate file with betavalues of each sample
header=`head -n 1 ${random_file} | cut -f 3-`
printf "probe\t${header}\t" > ${results_folder}"/"${basename}_methylationTEMP.txt

for ID in `ls ${methylation_folder}/*/*_hg38.txt | xargs -I {} basename {}`
do
	grep $ID ${samplesheet} \
		| cut -f 2 \
		| cut -d "." -f 6 \
		| cut -d "-" -f 1-5 \
		| sed -e "s/^/\"/" -e "s/$/\"/"
done | tr -s "\n" "\t" \
     | sed "s/\t$/\n/"  >> ${results_folder}"/"${basename}_methylationTEMP.txt

betavalues=`ls ${methylation_folder}/*/*_hg38.txt_betavalue`

paste ${results_folder}"/"${basename}_methylation-annotation.txt \
	${betavalues} \
	>> ${results_folder}"/"${basename}_methylationTEMP.txt \
&& mv ${results_folder}"/"${basename}_methylationTEMP.txt ${results_folder}"/"${basename}_methylation.txt


#Select only beans in a span no greather than the minlen (1500bp)
python3 parse_TSS.py \
	--input ${results_folder}"/"${basename}_methylation.txt \
	--output ${results_folder}"/"${basename}_methylation_byGene.txt \
	--genes $genesIDS \
	--minlen 1500
