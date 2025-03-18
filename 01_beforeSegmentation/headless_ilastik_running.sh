#! /bin/bash

#Define where is ilastik
ilastik="/home/fperez/apps/ilastik/ilastik-1.3.2post1-Linux/run_ilastik.sh"

#ilastik project path
ilastikfile="/mnt/d/users/fperez/NKI_TMAs_AF/Ilastik_labeling/Cellring_nuclei_background_Selected-Channels_3.ilp"

#Project TMA path
projectfolder="/mnt/d/users/fperez/NKI_TMAs_AF/"

#Dearray folder inside Project
dearrayfolder="dearray/selected_Channels/"

#Output pixel probabilities folder
probfolder="Probmaps_cyto_nuc_backg4"

#Suffix for input images
inputsuffix="*.tif"

#Suffix for output images
outputsuffix="_Probabilities.tiff"


#Running for each sample
for samplename in `ls -d ${projectfolder}/TMA* | xargs -I {} basename {}`
do
	inputimages=`ls ${projectfolder}/${samplename}/${dearrayfolder}/${inputsuffix}`

	#Creating outputfolder
	mkdir -p ${projectfolder}/${samplename}/${probfolder}

	#Running ilastik
	${ilastik} --headless \
		--project=${ilastikfile} \
		--output_format=tiff \
		--export_dtype=uint16 \
		--output_filename_format=${projectfolder}/${samplename}/${probfolder}/{nickname}${outputsuffix} \
		${inputimages}

done
