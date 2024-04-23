#glom project LSF commands
LUS=/lustre/scratch126/cellgen/team297/bs16
BEN=/nfs/team297/bs16

#Now look at the single cell data..
cd $LUS/current_projects/kidney_glomTI_response/data/seq_data/10X_multiome

#run amulet on multiome ATAC data  - not working very well. may have to do this in R
cd $LUS/current_projects/kidney_glomTI_response/data/seq_data/10X_multiome
GENOTYPE_COUNTS=./SOC_scripts/SOC_genotype_counts_multiome.txt
chmod +x ./amulet_atac_doublet.sh
#do SOC listwise.. 
while read line
do
SMP=$(echo "$line" | awk '{print $1'})
mkdir /lustre/scratch126/cellgen/team297/bs16/current_projects/kidney_glomTI_response/data/seq_data/10X_multiome/${SMP}/amulet_output
#now do bap
bsub -G teichlab -J amulet_$SMP -q "normal" -o $LUS/current_projects/kidney_glomTI_response/data/seq_data/10X_multiome/LSF_out/amulet_$SMP.txt   -R "select[mem>250000] rusage[mem=250000] span[hosts=1]" -M 250000 -n 5 ./amulet_atac_doublet.sh $LUS/current_projects/kidney_glomTI_response/data/seq_data/10X_multiome/${SMP} atac_fragments.tsv.gz single_cell.csv
done < ${GENOTYPE_COUNTS}

