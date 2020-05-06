#!/bin/bash

#How to:
#1.: fill in or check file paths of data files and script files
#2.: parameters are defined
#3.: make sure Rscript command is in your path
#3.: chomd +x
#4.: > ./runDTU.sh


#FILL IN YOUR FILE PATHS
##################################################################
#where to store outputs of this pipeline, holding DRIMSeq and DEXSeq result objects and processed result tables
out="./results/"
#specify your reference gtf file 
gtf="./referenceData/ucsc.hg19.gtf"
#specify the directory of the salmon output files for tximport
salmon="./rawData"
#specify sample metadata
csv="./metaData/phenoData.csv"


#PARAMETERS: levels of group variable, names of cohorts, covariates, expression filter threshold, tximport scaling method, postHoc filter (as described in: PMID: 30356428)

#################################################################
condition_names="Control:Case"
cohort_names="discovery:replication"
batch="rin:age_years:sex:Microglia_Genes:Oligo_Genes"
filter="10:10"
scalingMethod="scaledTPM"
posthocfilter="FALSE"


#SCRIPTS 
##################################################################
drim="./scripts/drimStage.R"
dex="./scripts/dexStage.R"
intersect="./scripts/bindToolRes.R"
pipes=($drim $dex) 


timestamp=$(date +%s)
timestamp="04-05-20"
LOG_FILE="$timestamp"
LOG_FILE="./logs/"$LOG_FILE".txt"

exec > >(tee -a ${LOG_FILE} )
exec 2> >(tee -a ${LOG_FILE} >&2)




" \e[92mrunning three DTU pipelines in sequence\e[0m"

if [ -d "$out" ]; then
	echo -e " \e[92mSaving results to "$out"\e[0m"
else
	echo -e "\e[31mDir "$out" not found"
fi
if [ -f "$gtf" ]; then
	echo -e " \e[92musing "$gtf" as annotation file\e[0m"
else
	echo -e "\e[31mfile "$gtf" not found\e[0m"
fi
if [ -d "$salmon" ]; then
	echo -e " \e[92musing transcript abundances found in "$salmon"\e[0m"
else
	echo -e "\e[31mDir "$salmon" not found\e[0m"
fi


pids=""
script_names=()
for script in "${pipes[@]}"
do
(	#split by slash, reverse and select first, which is filename with extension
	script_name=$(echo "$script" | rev | cut -d/ -f1 | rev)
	#strip filename of extension
	script_name="${script_name%.*}"
#save name for later when acessing results for intersection script
	script_names+=("$script_name")
	echo -e " \e[92mNow running $script_name\e[0m"

	mkdir -p ""$out"$script_name" 

	if Rscript --vanilla "$script" -o ""$out""$script_name"" -f "$filter" -g "$gtf" -m "$csv" -s "$salmon" -z "$scalingMethod" -n "$cohort_names" -c "$condition_names" -d "$posthocfilter" -t "$timestamp" -b "$batch" ; then
		echo -e " \e[92mDone!Check Results in "$out"\e[0m" 
	else 
		echo -e "\e[31m"$script_name" pipeline did not finish without error\e[0m"
	fi 
) &
pids="$pids $!"
done

wait $pids

echo -e "\e[92mMoving on to intersecting results\e[0m"
gene_paths=()
Ds_paths=()
Ss_paths=()
filtInfo_paths=()
Ds_unfilt_paths=()

for script in "${pipes[@]}"
do	

	toolResFolder=$(basename "$script" .R)
	Ds=""$out""$toolResFolder"/Ds/Ds_"$timestamp".rds"
	Ds_unfilt=""$out""$toolResFolder"/Ds/Ds_preFilt_"$timestamp".rds"
	filtInfo=""$out""$toolResFolder"/Ds/paramFilt_"$timestamp".rds"
	gene=""$out""$toolResFolder"/genes/genelist_"$timestamp".rds"


FileArr=("$Ds" "$Ds_unfilt" "$filtInfo" "$gene")
 

	for file in ${FileArr[*]}; do
	    if [ ! -f "$file" ] ;  then 
	    	echo ""$file" doesn't exist"
	    	exit 1
	    fi
	done



	Ds_paths+=("$Ds")
	Ds_unfilt_paths+=("$Ds_unfilt")
	filtInfo_paths+=("$filtInfo")	
	gene_paths+=("$gene")

done

echo "$Ds_paths"
Ds_paths=$(printf ":%s" "${Ds_paths[@]}")
Ds_paths=${Ds_paths:1}
Ds_paths=${Ds_paths#":"}
echo "$Ds_paths"
Ds_unfilt_paths=$(printf ":%s" "${Ds_unfilt_paths[@]}")
Ds_unfilt_paths=${Ds_unfilt_paths:1}
Ds_unfilt_paths=${Ds_unfilt_paths#":"}
echo "$Ds_unfilt_paths"
filtInfo_paths=$(printf ":%s" "${filtInfo_paths[@]}")
filtInfo_paths=${filtInfo_paths:1}
filtInfo_paths=${filtInfo_paths#":"}
echo "$filtInfo_paths"
gene_paths=$(printf ":%s" "${gene_paths[@]}")
gene_paths=${gene_paths:1}
echo "$gene_paths"



mkdir -p ""$out"rds/"
if Rscript  --vanilla "$intersect" -t "$timestamp" -g "$gene_paths" -d "$Ds_paths" -n "$cohort_names" -o ""$out"rds/" -f "$filtInfo_paths" -u "$Ds_unfilt_paths" ; then
		echo -e " \e[92mSucessfully intersected gene_lists, results can be found in "$out"\e[0m"
	else
		echo -e "\e[31mIntersecting geneLists failed\e[0m"
	fi

