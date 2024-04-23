#LSF commands
LUS=/lustre/scratch126/cellgen/team297/bs16
BEN=/nfs/team297/bs16

#Now look at the single cell data..
cd $LUS/current_projects/kidney_glomTI_response/data/seq_data
KG_PATH=$LUS/current_projects/kidney_glomTI_response/data/seq_data

#then grab our already aligned adult datasets.
cd 10X_3p
mkdir LSF_out
for smp in $(cat /lustre/scratch126/cellgen/team297/bs16/current_projects/kidney_glomTI_response/data/seq_data/10X_3p/3p_irods.txt)
do
	bsub -G teichlab -J ${smp} -q 'normal' -o ./LSF_out/iget_${smp}.txt -R "select[mem>100000] rusage[mem=100000] span[hosts=1]" -M 100000 -n 1 iget -r ${smp}
done

#run velocyto on 10X3p kidney
cd $LUS/current_projects/kidney_glomTI_response/data/seq_data/10X_3p
CHANNELS=./3p_channels.txt
while read line
do
SMP=$(echo "$line" | awk '{print $1'})
#now do velocyto
bsub -G teichlab -J ${SMP}_velocyto -q "long" -o ./LSF_out/velocyto_$SMP.txt -R "select[mem>500000] rusage[mem=500000] span[hosts=1]" -M 500000 -n 1 ./run_velocity_10x3p.sh ${SMP}
done < ${CHANNELS}


#do denoising using scAR
cd $LUS/current_projects/kidney_glomTI_response/data/seq_data/10X_3p
CHANNELS=./3p_channels.txt
while read line
do
SMP=$(echo "$line" | awk '{print $1'})
#now do scAR
bsub -G teichlab -J ${SMP} -q gpu-huge -o ./LSF_out/scAR_${SMP}.txt -R "select[mem>100000] rusage[mem=100000, ngpus_physical=1.00] span[hosts=1]" -gpu "mode=shared:j_exclusive=yes" -M 100000 -n 1 /nfs/team297/bs16/tools/conda_envs/scvi-env/bin/python ./scar_denoising.py ${SMP}
done < ${CHANNELS}

#grab our 10X ATAC datasets
cd 10X_snATACseq
mkdir LSF_out
for smp in $(cat /lustre/scratch126/cellgen/team297/bs16/current_projects/kidney_glomTI_response/data/seq_data/10X_snATACseq/ATAC_irods.txt)
do
	bsub -G teichlab -J ${smp} -q 'normal' -o ./LSF_out/iget_${smp}.txt -R "select[mem>100000] rusage[mem=100000] span[hosts=1]" -M 100000 -n 1 iget -r ${smp}
done

#then grab our 10X multiome datasets
cd $LUS/current_projects/kidney_glomTI_response/data/seq_data/10X_multiome
mkdir LSF_out
for smp in $(cat /lustre/scratch126/cellgen/team297/bs16/current_projects/kidney_glomTI_response/data/seq_data/10X_multiome/multiome_irods.txt)
do
	bsub -G teichlab -J ${smp} -q 'normal' -o ./LSF_out/iget_${smp}.txt -R "select[mem>100000] rusage[mem=100000] span[hosts=1]" -M 100000 -n 1 iget -r ${smp}
done

#run velocyto on 10X multiome data
cd $LUS/current_projects/kidney_glomTI_response/data/seq_data/10X_multiome
CHANNELS=./multiome_channels.txt
chmod +x ./run_velocity_10x_multiome.sh
while read line
do
SMP=$(echo "$line" | awk '{print $1'})
#now do velocyto
bsub -G teichlab -J ${SMP}_velocyto -q "long" -o ./LSF_out/velocyto_$SMP.txt -R "select[mem>500000] rusage[mem=500000] span[hosts=1]" -M 500000 -n 1 ./run_velocity_10x_multiome.sh ${SMP}
done < ${CHANNELS}
