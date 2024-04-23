cd /lustre/scratch126/cellgen/team297/bs16/current_projects/kidney_glomTI_response/data/seq_data/10X_multiome
bgzip aggregated_fragments.tsv
tabix --preset=bed -0 aggregated_fragments.tsv.gz
