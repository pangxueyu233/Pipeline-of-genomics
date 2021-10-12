# Step2 . somatic SNVs detected and annotated

This page showed the somatic SNV called by `Mutect implemented in GATK4` without paired normal sample. To better visualize the SNVs results, we used `vep` to annotate the SNVs info, and transfer `.vcf` files to `.maf` files. Then, we could used `maftools >= v2.8` to summary the somatic SNVs results.

- 2.1 somatic SNVs detected by Mutec2 *(without paired normal sample)*

~~~shell
####execute codes
cat samplelist | while read id ; do
arr=($id)
fq2=${arr[1]}'_R1.fq.gz'
fq1=${arr[1]}'_R2.fq.gz'
sample=${arr[0]}
echo ${fq2}
echo ${fq1}
echo $sample
$GATK --java-options "-Xmx100G -Djava.io.tmpdir=./" FixMateInformation \
-I $BAMOUT/${sample}_marked.bam \
-O $BAMOUT/${sample}_fix.bam \
-SO coordinate --CREATE_INDEX true 

$GATK --java-options "-Xmx100G -Djava.io.tmpdir=./" BaseRecalibrator \
-R $GENOME -I $BAMOUT/${sample}_fix.bam \
--known-sites $mm10/chr_mm10.INDELS.dbSNP142.vcf \
--known-sites $mm10/chr_mm10.dbSNP142.vcf \
-O $BAMOUT/$sample.wes.recal_data.table

$GATK --java-options "-Xmx100G -Djava.io.tmpdir=./" ApplyBQSR \
-R $GENOME -I $BAMOUT/${sample}_fix.bam \
-bqsr $BAMOUT/$sample.wes.recal_data.table \
-O $BAMOUT/${sample}_BQSR.bam 

$GATK --java-options "-Xmx100G -Djava.io.tmpdir=./" Mutect2 \
-L $EXON_BED \
-R $GENOME \
-L $EXON_BED # add the exon bed
-I $BAMOUT/${sample}_BQSR.bam \
-O $SNVOUT/$sample.Mutect2only_tumor.vcf.gz ;
done

rm -r *_fix.*
~~~

- 2.2. SNVs annotated by vep and maf files generated *(without paired normal sample)*

~~~shell
####execute codes
cat samplelist | while read id ; do
arr=($id)
tumor_name=${arr[0]}
echo $tumor_name
gunzip $tumor_name.Mutect2only_tumor.vcf.gz
perl /mnt/data/userdata/xiangyu/programme/vcf2maf-master/vcf2maf.pl \
--input-vcf $SNVOUT/$tumor_name.Mutect2only_tumor.vcf \
--output-maf $SNVOUT/$tumor_name.Mutect2only_tumor.vep.maf \
--tumor-id $tumor_name \
--vcf-tumor-id $tumor_name \
--ref-fasta /mnt/data/userdata/xiangyu/programme/genome_index/bwa_index/bwa_mm10_index/genome.fa \
--vep-path /mnt/data/userdata/xiangyu/programme/ensembl-vep-release-95.3/vep_all \
--vep-data /mnt/data/userdata/xiangyu/programme/ensembl-vep-release-95.3/vep_all/.vep \
--vep-forks 40 --species mus_musculus --species mus_musculus --ncbi-build GRCm38 ;
done

cd /mnt/data/user_data/zhanglu/01_WES_NSCLC_wml/BWA_out
scp -P 22 -r *_BQSR* xiangyu@172.16.1.20:/mnt/data/userdata/xiangyu/workshop/NSCLC_PD1_WES/bam_files
~~~

- 2.3 somatic SNVs detected by `samtools + varscan`*(with paired normal sample)*

~~~shell
vim sample_file
NSCLC_10
NSCLC_17
NSCLC_37
NSCLC_54
NSCLC_7
NSCLC_83
NSCLC_8

GENOME=/mnt/data/user_data/xiangyu/8922_server/programme/genome_index/bwa_index/bwa_mm10_index/genome.fa
mkdir mpileup
cat sample_file | while read id ; do
arr=($id)
sample=${arr[0]}
echo $sample
echo mpileupping...
samtools mpileup -d 1000 -q 1 -Q 15 -A -f $GENOME ${sample}_BQSR.bam > ./mpileup/$sample.mpileup
echo varscanning...
java -jar $VarScan mpileup2cns -Xmx100G  mpileup2cns \
./mpileup/$sample.mpileup --variants --strand-filter 0 --p-value 0.01 --min-var-freq 0.05 --output-vcf 1 > ./mpileup/$sample.varscan.vcf ;
done
~~~

- 2.4 SNVs annotated by vep and maf files generated *(with paired normal sample)*

~~~shell
scp -P 22 -r /mnt/data/user_data/zhanglu/01_WES_NSCLC_wml/BWA_out/mpileup/*.vcf xiangyu@172.16.1.20:/mnt/data/userdata/xiangyu/workshop/NSCLC_PD1_WES

cd /mnt/data/userdata/xiangyu/workshop/NSCLC_PD1_WES

sed 's/Sample1/NSCLC_10/g' NSCLC_10.vcf > NSCLC_10.new.vcf
sed 's/Sample1/NSCLC_17/g' NSCLC_17.vcf > NSCLC_17.new.vcf
sed 's/Sample1/NSCLC_37/g' NSCLC_37.vcf > NSCLC_37.new.vcf
sed 's/Sample1/NSCLC_54/g' NSCLC_54.vcf > NSCLC_54.new.vcf
sed 's/Sample1/NSCLC_7/g' NSCLC_7.vcf > NSCLC_7.new.vcf
sed 's/Sample1/NSCLC_83/g' NSCLC_83.vcf > NSCLC_83.new.vcf
sed 's/Sample1/NSCLC_8/g' NSCLC_8.vcf > NSCLC_8.new.vcf

cat sample_file | while read id ; do
arr=($id)
tumor_name=${arr[0]}
echo $tumor_name
perl /mnt/data/userdata/xiangyu/programme/vcf2maf-master/vcf2maf.pl \
--input-vcf $tumor_name.new.vcf \
--output-maf $tumor_name.vep.maf \
--tumor-id $tumor_name \
--vcf-tumor-id $tumor_name \
--ref-fasta /mnt/data/userdata/xiangyu/programme/genome_index/bwa_index/bwa_mm10_index/genome.fa \
--vep-path /mnt/data/userdata/xiangyu/programme/ensembl-vep-release-95.3/vep_all \
--vep-data /mnt/data/userdata/xiangyu/programme/ensembl-vep-release-95.3/vep_all/.vep \
--vep-forks 40 --species mus_musculus --species mus_musculus --ncbi-build GRCm38 ;
done
mkdir vep.maf
mv *vep.maf ./vep.maf/

scp -P 22 -r /mnt/data/userdata/xiangyu/workshop/NSCLC_PD1_WES/vep.maf/*vep.maf xiangyu@172.16.1.21:/mnt/data/user_data/zlu/01_job/WES_NSCLC_wml/snv_out/xy_workshop/mpileup
~~~

- 2.5 somatic SNVs detected by Mutec2 *(with paired normal sample)*



- 2.6 SNVs annotated by vep and maf files generated *(with paired normal sample)*

