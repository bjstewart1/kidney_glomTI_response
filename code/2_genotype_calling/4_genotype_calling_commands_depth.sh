LUS=/lustre/scratch126/cellgen/team297/bs16
BEN=/nfs/team297/bs16

cd $LUS/current_projects/kidney_glomTI_response/data/genotypes/depth_vcf
#paths to tools
PICARD=/nfs/team297/bs16/tools/picard.jar
BCTOOLS=$BEN/tools/bcftools-1.3.1/bin/bcftools 

DEPTH_SAMPLES=./depth_samples.txt
while read line
do
SMP=$(echo "$line" | awk '{print $1'})
echo $SMP
CLUSTER=$(echo "$line" | awk '{print $3'})
echo $CLUSTER
VCF_FILE_PATH=$(echo "$line" | awk '{print $4}')
#rsync the cluster genotype file
rsync ${VCF_FILE_PATH}/cluster_genotypes.vcf ${SMP}.vcf
java -jar ${PICARD} FixVcfHeader -I ${SMP}.vcf -O ${SMP}_reheader.vcf
rm ${SMP}.vcf 
mv ${SMP}_reheader.vcf ${SMP}.vcf 
#select just the sample we want
${BCTOOLS} view -s ${CLUSTER} ${SMP}.vcf > ${SMP}_filtered.vcf
rm ${SMP}.vcf
echo ${SMP} > ${SMP}_sample_name.txt
${BCTOOLS} reheader -s ${SMP}_sample_name.txt ${SMP}_filtered.vcf > ${SMP}_filtered_renamed.vcf
rm ${SMP}_sample_name.txt
rm ${SMP}_filtered.vcf 
done < ${DEPTH_SAMPLES}

#now merge
parallel bgzip {} ::: *filtered_renamed.vcf
for f in ./*filtered_renamed.vcf.gz; do tabix -p vcf -f $f;done
${BCTOOLS} merge *filtered_renamed.vcf.gz > merged_vcf.vcf
rm *filtered_renamed.vcf.gz
rm *idx
rm *tbi
$BEN/tools/bcftools-1.3.1/bin/bcftools view merged_vcf.vcf | head -n 200

#for the depth samples we rerun souporcell but in supervised mode.
rsync ./merged_vcf.vcf $LUS/current_projects/kidney_glomTI_response/data/seq_data/10X_3p/3p_reference/depth_samples.vcf

##### 10X3P kidney DEPTH SAMPLES
cd $LUS/current_projects/kidney_glomTI_response/data/seq_data/10X_3p
chmod +x ./SOC_scripts/run_SOC_supervised.sh
GENOTYPE_COUNTS=./SOC_scripts/SOC_genotype_counts_depth.txt
#get singularity up and running
prepend_path MODULEPATH /software/modules/
module load ISG/singularity/3.9.0
#do SOC listwise.. 
while read line
do
SMP=$(echo "$line" | awk '{print $1'})
KVAL=$(echo "$line" | awk '{print $2'})
#now do SOC but in a supervised mode
bsub -G teichlab -J SOC_10x3p_$SMP -q "yesterday" -o ./LSF_out/souporcell_supervised_$SMP.txt  -R "select[mem>100000] rusage[mem=100000] span[hosts=1]" -M 100000 -n 2 ./SOC_scripts/run_SOC_supervised.sh ./${SMP} ${KVAL} ./filtered_feature_bc_matrix/barcodes.tsv.gz ./${SMP}/SOC_common_variants ./3p_reference/depth_samples.vcf ./3p_reference/fasta/genome.fa possorted_genome_bam.bam False
done < ${GENOTYPE_COUNTS}
