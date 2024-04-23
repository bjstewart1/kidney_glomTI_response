 #grab the zenodo dataset
cd /lustre/scratch117/cellgen/team297/bs16/ENCODE_reference
zenodo_get 10.5281/zenodo.3838751

#unzip the TF associated DHS file
tar -xf TF_associated_DHSs_hg38.tar.gz
mkdir TF_associated_DHS
#move these around
mv ./*bed ./TF_associated_DHS/
mv ./TF_associated_DHS/DHS_ENCODE_ref.bed ./

#get a fasta from this bed file
FASTA=/lustre/scratch117/cellgen/team297/bs16/refdata-cellranger-atac-GRCh38-1.2.0/fasta/genome.fa
BT=/nfs/team297/bs16/tools/bedtools2/bin/bedtools

cd /lustre/scratch117/cellgen/team297/bs16/ENCODE_reference
mkdir LSF_out
bsub -G teichlab -J peak_to_bed -q "yesterday" -o ./LSF_out/bed_to_fasta.txt  -R "select[mem>100000] rusage[mem=100000] span[hosts=1]" -M 100000 -n 2 $BT getfasta -fi ${FASTA} -fo ./DHS_ENCODE_ref.fa -bed ./DHS_ENCODE_ref.bed


#take a look at this - this gives a fasta of all the regions we are interested in... 
cat ./DHS_ENCODE_ref.fa | head -n 50

#index the fasta
samtools faidx ./DHS_ENCODE_ref.fa

#now do kallisto index
bsub -G teichlab -J kallisto_indexing -q "yesterday" -o ./LSF_out/kallisto_indexing.txt  -R "select[mem>200000] rusage[mem=200000] span[hosts=1]" -M 200000 -n 4 /nfs/team297/bs16/tools/kallisto/kallisto index -i DHS_ENCODE_ref.ki DHS_ENCODE_ref.fa

#then make the 'map' which gives us a "transcript to genes mapping" equivalent
cat DHS_ENCODE_ref.fa | awk '{if($1~/>/)print $1"\t"$1"\t"$1}' > tr2g.txt
sed -i 's/>//g' tr2g.txt


#now make a reference where we merge regions closer than 500 bp
${BT} sort -i ./DHS_ENCODE_ref.bed | ${BT} merge -d 500 > ./DHS_ENCODE_ref_500.bed

bsub -G teichlab -J peak_to_bed -q "yesterday" -o ./LSF_out/bed_to_fasta_500.txt  -R "select[mem>100000] rusage[mem=100000] span[hosts=1]" -M 100000 -n 2 $BT getfasta -fi ${FASTA} -fo ./DHS_ENCODE_ref_500.fa -bed ./DHS_ENCODE_ref_500.bed

cat ./DHS_ENCODE_ref_500.fa | head -n 50

#index the fasta
samtools faidx ./DHS_ENCODE_ref_500.fa

#then make the 'map' which gives us a "transcript to genes mapping" equivalent
cat DHS_ENCODE_ref_500.fa | awk '{if($1~/>/)print $1"\t"$1"\t"$1}' > tr2g_500.txt
sed -i 's/>//g' tr2g_500.txt

#now do kallisto index
bsub -G teichlab -J kallisto_indexing -q "yesterday" -o ./LSF_out/kallisto_indexing_500.txt  -R "select[mem>200000] rusage[mem=200000] span[hosts=1]" -M 200000 -n 4 /nfs/team297/bs16/tools/kallisto/kallisto index -i DHS_ENCODE_ref_500.ki DHS_ENCODE_ref_500.fa


