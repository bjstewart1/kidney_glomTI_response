cd /lustre/scratch126/cellgen/team297/bs16/current_projects/kidney_glomTI_response/data/raw/10X_multiome/archR_files
#mkdir LSF_out
chmod +x ./make_h5_gact_aggregate.py
ARROW=/lustre/scratch126/cellgen/team297/bs16/current_projects/kidney_glomTI_response/data/raw/10X_multiome/archR_files/aggregated_fragments.arrow
bsub -G teichlab -J arrow_agg -q normal -o ./LSF_out/arrow_to_h5ad_aggregate.txt -R "select[mem>100000] rusage[mem=100000] span[hosts=1]" -M 100000 -n 1 -M 100000 -n 1 /nfs/team297/bs16/tools/conda_envs/archr_env/bin/python ./make_h5_gact_aggregate.py $ARROW
