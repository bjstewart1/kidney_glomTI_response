LUS=/lustre/scratch126/cellgen/team297/bs16
BEN=/nfs/team297/bs16

cd $LUS/current_projects/kidney_glomTI_response/data/raw/10X_multiome/archR_files
chmod +x do_aggregate_arrow.R
conda activate r-4.2
bsub -G teichlab -J aggregate_arrow -q "yesterday" -o ./LSF_out/aggregate_arrow.txt  -R "select[mem>70000] rusage[mem=70000] span[hosts=1]" -M 70000 -n 1 Rscript do_aggregate_arrow.R
