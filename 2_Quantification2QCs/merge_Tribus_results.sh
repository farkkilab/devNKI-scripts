#! /bin/bash

#Just an script for concatenating each core results per slide
for j in `ls -d TMA_*`
	do cd ${j}/Tribus_celltype/
	echo "Core_Names,Cellid,GlobalCellType" > ${j}_celltypes.csv
       	cut -d "," -f 1,2,4 *_cellTypesNapari.csv \
		| sed "s/_Probabilities_cytomask2//" \
		| sed "s/core//" \
		| sed "s/\"//g" \
		| grep -v "^Core" >>${j}_celltypes.csv
       	cd ../../
done
