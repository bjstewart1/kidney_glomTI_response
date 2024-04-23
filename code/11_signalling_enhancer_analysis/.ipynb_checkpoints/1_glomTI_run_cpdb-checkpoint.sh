cd /lustre/scratch126/cellgen/team297/bs16/current_projects/kidney_glomTI_response/code/16_signalling_analysis
source activate cpdb
bsub -G teichlab -J glomTI_cpdb -q "yesterday" -o ./cpdb_glomTI.txt -R "select[mem>100000] rusage[mem=100000] span[hosts=1]" -M 100000 -n 1 python glom_niche_cpdb.py 
