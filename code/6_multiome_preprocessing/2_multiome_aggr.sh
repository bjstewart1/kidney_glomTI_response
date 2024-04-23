CRA=/nfs/team297/bs16/tools/cellranger-arc-2.0.0/cellranger-arc
cd /lustre/scratch126/cellgen/team297/bs16/current_projects/kidney_glomTI_response/data/seq_data/10X_multiome

bsub -G teichlab -J aggregate_multiome -q "hugemem" -o ./LSF_out/multiome_aggr.txt -R "select[mem>500000] rusage[mem=500000] span[hosts=1]" -M 500000 -n 10 $CRA aggr --id=aggregate_multiome --csv=./aggregation.csv --reference=/nfs/srpipe_references/downloaded_from_10X/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 --localcores 10

#get this directly from pipeline

rsync -r /lustre/scratch125/ids/team117/npg/ygu/10X/bs16/cellranger-arc200_aggr_RT771977 /lustre/scratch126/cellgen/team297/bs16/current_projects/kidney_glomTI_response/data/seq_data/10X_multiome/

mv /lustre/scratch126/cellgen/team297/bs16/current_projects/kidney_glomTI_response/data/seq_data/10X_multiome/cellranger-arc200_aggr_RT771977 /lustre/scratch126/cellgen/team297/bs16/current_projects/kidney_glomTI_response/data/seq_data/10X_multiome/multiome_aggr