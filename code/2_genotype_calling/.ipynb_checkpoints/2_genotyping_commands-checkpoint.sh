#genotyping LSF commands
LUS=/lustre/scratch126/cellgen/team297/bs16
BEN=/nfs/team297/bs16

cd $LUS/current_projects/kidney_glomTI_response/data/seq_data/10X_3p
mkdir 3p_reference
rsync -r /nfs/srpipe_references/downloaded_from_10X/refdata-gex-GRCh38-2020-A/fasta ./3p_reference/

##### 10X3P kidney GLOM TI RESPONSE 
cd $LUS/current_projects/kidney_glomTI_response/data/seq_data/10X_3p
chmod +x ./SOC_scripts/run_SOC.sh
GENOTYPE_COUNTS=./SOC_scripts/SOC_genotype_counts.txt
#get singularity up and running
prepend_path MODULEPATH /software/modules/
module load ISG/singularity/3.9.0
#do SOC listwise.. 
while read line
do
SMP=$(echo "$line" | awk '{print $1'})
KVAL=$(echo "$line" | awk '{print $2'})
#now do SOC
bsub -G teichlab -J SOC_10x3p_$SMP -q "long" -o ./LSF_out/souporcell_$SMP.txt  -R "select[mem>500000] rusage[mem=500000] span[hosts=1]" -M 500000 -n 2 ./SOC_scripts/run_SOC.sh ./${SMP} ${KVAL} ./filtered_feature_bc_matrix/barcodes.tsv.gz ./${SMP}/SOC_common_variants ./3p_reference/resort+chr.vcf ./3p_reference/fasta/genome.fa possorted_genome_bam.bam False
done < ${GENOTYPE_COUNTS}


##### 10X3P kidney GLOM TI RESPONSE 
cd $LUS/current_projects/kidney_glomTI_response/data/seq_data/10X_3p
chmod +x ./SOC_scripts/run_SOC.sh
cat ./SOC_scripts/SOC_genotype_counts.txt | tail -n 16 > ./SOC_scripts/SOC_genotype_counts_extra.txt
GENOTYPE_COUNTS=./SOC_scripts/SOC_genotype_counts_extra.txt
#get singularity up and running
prepend_path MODULEPATH /software/modules/
module load ISG/singularity/3.9.0
#do SOC listwise.. 
while read line
do
SMP=$(echo "$line" | awk '{print $1'})
KVAL=$(echo "$line" | awk '{print $2'})
#now do SOC
bsub -G teichlab -J SOC_10x3p_$SMP -q "long" -o ./LSF_out/souporcell_$SMP.txt  -R "select[mem>500000] rusage[mem=500000] span[hosts=1]" -M 500000 -n 2 ./SOC_scripts/run_SOC.sh ./${SMP} ${KVAL} ./filtered_feature_bc_matrix/barcodes.tsv.gz ./${SMP}/SOC_common_variants ./3p_reference/resort+chr.vcf ./3p_reference/fasta/genome.fa possorted_genome_bam.bam False
done < ${GENOTYPE_COUNTS}



##### ATAC kidney 
cd $LUS/current_projects/kidney_glomTI_response/data/seq_data/10X_snATACseq
chmod +x ./SOC_scripts/run_SOC.sh
GENOTYPE_COUNTS=./SOC_scripts/SOC_genotype_counts_ATAC.txt
#get singularity up and running
prepend_path MODULEPATH /software/modules/
module load ISG/singularity/3.9.0
#do SOC listwise.. 
while read line
do
SMP=$(echo "$line" | awk '{print $1'})
KVAL=$(echo "$line" | awk '{print $2'})
#now do SOC
bsub -G teichlab -J SOC_ATAC_$SMP -q "long" -o ./LSF_out/souporcell_$SMP.txt  -R "select[mem>500000] rusage[mem=500000] span[hosts=1]" -M 500000 -n 2 ./SOC_scripts/run_SOC.sh ./${SMP} ${KVAL} filtered_peak_bc_matrix/barcodes.tsv ./${SMP}/SOC_common_variants ./atac_reference/resort+chr.vcf ./atac_reference/fasta/genome.fa possorted_bam.bam True
done < ${GENOTYPE_COUNTS}


##### multiome kidney #we do this on the ATAC libraries
cd $LUS/current_projects/kidney_glomTI_response/data/seq_data/10X_multiome
chmod +x ./SOC_scripts/run_SOC.sh
GENOTYPE_COUNTS=./SOC_scripts/SOC_genotype_counts_multiome.txt
#get singularity up and running
prepend_path MODULEPATH /software/modules/
module load ISG/singularity/3.9.0
#do SOC listwise.. 
while read line
do
SMP=$(echo "$line" | awk '{print $1'})
KVAL=$(echo "$line" | awk '{print $2'})
#now do SOC
bsub -G teichlab -J SOC_multiome_$SMP -q "long" -o ./LSF_out/souporcell_$SMP.txt  -R "select[mem>500000] rusage[mem=500000] span[hosts=1]" -M 500000 -n 2 ./SOC_scripts/run_SOC.sh ./${SMP} ${KVAL} filtered_feature_bc_matrix/barcodes.tsv.gz ./${SMP}/SOC_common_variants ./multiome_reference/resort+chr.vcf ./multiome_reference/fasta/genome.fa atac_possorted_bam.bam True
done < ${GENOTYPE_COUNTS}
