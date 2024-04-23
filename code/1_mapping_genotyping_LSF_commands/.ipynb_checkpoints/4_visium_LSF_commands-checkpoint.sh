#LSF commands
LUS=/lustre/scratch126/cellgen/team297/bs16
BEN=/nfs/team297/bs16

#Now look at the single cell data..
cd $LUS/current_projects/kidney_glomTI_response/data/seq_data
KG_PATH=$LUS/current_projects/kidney_glomTI_response/data/seq_data

#then grab our already aligned adult datasets.
cd visium
mkdir LSF_out
for smp in $(cat /lustre/scratch126/cellgen/team297/bs16/current_projects/kidney_glomTI_response/data/seq_data/visium/visium_irods.txt)
do
	bsub -G teichlab -J ${smp} -q 'yesterday' -o ./LSF_out/iget_${smp}.txt -R "select[mem>50000] rusage[mem=50000] span[hosts=1]" -M 50000 -n 1 iget -r ${smp}
done

#tweak the tissue positions file name
cd $LUS/current_projects/kidney_glomTI_response/data/seq_data/visium
for smp in $(cat /lustre/scratch126/cellgen/team297/bs16/current_projects/kidney_glomTI_response/data/seq_data/visium/visium_slides.txt)
do
rsync ./${smp}/spatial/tissue_positions.csv ./${smp}/spatial/tissue_positions_list.csv
done


#get fresh frozen visium
cd $LUS/current_projects/kidney_glomTI_response/data/seq_data/visium
for smp in $(cat /lustre/scratch126/cellgen/team297/bs16/current_projects/kidney_glomTI_response/data/seq_data/visium/ff_visium_irods.txt)
do
	bsub -G teichlab -J ${smp} -q 'yesterday' -o ./LSF_out/iget_${smp}.txt -R "select[mem>50000] rusage[mem=50000] span[hosts=1]" -M 50000 -n 1 iget -r ${smp}
done



#tweak the tissue positions file name
cd $LUS/current_projects/kidney_glomTI_response/data/seq_data/visium
for smp in $(cat /lustre/scratch126/cellgen/team297/bs16/current_projects/kidney_glomTI_response/data/seq_data/visium/visium_slides_extra.txt)
do
rsync ./${smp}/spatial/tissue_positions.csv ./${smp}/spatial/tissue_positions_list.csv
done
