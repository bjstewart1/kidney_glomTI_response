LUS=/lustre/scratch117/cellgen/team297/bs16
BEN=/nfs/team297/bs16
cd $LUS/current_projects/kidney_glomTI_response/code/3_public_data_preprocessing
mkdir LSF_out
for smp in $(cat /lustre/scratch117/cellgen/team297/bs16/10X_data/science_paper_data/channels_use.txt)
do
#we do this on a normal queue for now, but might be better to do on a GPU queue
bsub -G teichlab -J ${smp} -q gpu-huge -o ./LSF_out/${smp}_scAR.txt -R "select[mem>100000] rusage[mem=100000, ngpus_physical=1.00] span[hosts=1]" -gpu "mode=shared:j_exclusive=yes" -M 100000 -n 1 /nfs/team297/bs16/tools/conda_envs/scvi-env/bin/python ./scar_denoising.py ${smp}
done
