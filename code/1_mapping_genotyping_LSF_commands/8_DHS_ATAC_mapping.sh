LUS=/lustre/scratch126/cellgen/team297/bs16
cd $LUS/kidney_glomTI_response/data/seq_data/10X_snATACseq

chmod +x ./map_ATAC_DHS.sh
chmod +x ./map_ATAC_DHS_500.sh


#run for DHS sites where we have merged closer than 500 bp
LUS=/lustre/scratch126/cellgen/team297/bs16
cd $LUS/current_projects/kidney_glomTI_response/data/seq_data/10X_snATACseq
DHS_IN=./DHS_input.txt
while read line
do
DIR=$(echo "$line" | awk '{print $1'})
CHEM=$(echo "$line" | awk '{print $2'})
IR=$(echo "$line" | awk '{print $3'})
echo ${IR}
bsub -G teichlab -J DHS_mapping_${DIR} -q "long" -o ./LSF_out/DHS_mapping_500_${IR}.txt  -R "select[mem>100000] rusage[mem=100000] span[hosts=1]" -M 100000 -n 2 ./map_ATAC_DHS_500.sh ${DIR} ${CHEM}
done < ${DHS_IN}

#run for all DHS sites
LUS=/lustre/scratch126/cellgen/team297/bs16
cd $LUS/current_projects/kidney_glomTI_response/data/seq_data/10X_snATACseq
DHS_IN=./DHS_input.txt
while read line
do
DIR=$(echo "$line" | awk '{print $1'})
CHEM=$(echo "$line" | awk '{print $2'})
IR=$(echo "$line" | awk '{print $3'})
echo ${IR}
bsub -G teichlab -J DHS_mapping_${DIR} -q "long" -o ./LSF_out/DHS_mapping_${IR}.txt  -R "select[mem>100000] rusage[mem=100000] span[hosts=1]" -M 100000 -n 2 ./map_ATAC_DHS.sh ${DIR} ${CHEM}
done < ${DHS_IN}


#now with multiome datasets
LUS=/lustre/scratch126/cellgen/team297/bs16
cd $LUS/current_projects/kidney_glomTI_response/data/seq_data/10X_multiome

chmod +x ./map_ATAC_DHS.sh
chmod +x ./map_ATAC_DHS_500.sh

#run for DHS sites where we have merged closer than 500 bp
LUS=/lustre/scratch126/cellgen/team297/bs16
cd $LUS/current_projects/kidney_glomTI_response/data/seq_data/10X_multiome
DHS_IN=./DHS_input.txt
while read line
do
DIR=$(echo "$line" | awk '{print $1'})
CHEM=$(echo "$line" | awk '{print $2'})
IR=$(echo "$line" | awk '{print $3'})
echo ${IR}
bsub -G teichlab -J DHS_mapping_${DIR} -q "long" -o ./LSF_out/DHS_mapping_500_${IR}.txt  -R "select[mem>100000] rusage[mem=100000] span[hosts=1]" -M 100000 -n 2 ./map_ATAC_DHS_500.sh ${DIR} ${CHEM}
done < ${DHS_IN}

#run for all DHS sites
LUS=/lustre/scratch126/cellgen/team297/bs16
cd $LUS/current_projects/kidney_glomTI_response/data/seq_data/10X_snATACseq
DHS_IN=./DHS_input.txt
while read line
do
DIR=$(echo "$line" | awk '{print $1'})
CHEM=$(echo "$line" | awk '{print $2'})
IR=$(echo "$line" | awk '{print $3'})
echo ${IR}
bsub -G teichlab -J DHS_mapping_${DIR} -q "long" -o ./LSF_out/DHS_mapping_${IR}.txt  -R "select[mem>100000] rusage[mem=100000] span[hosts=1]" -M 100000 -n 2 ./map_ATAC_DHS.sh ${DIR} ${CHEM}
done < ${DHS_IN}
