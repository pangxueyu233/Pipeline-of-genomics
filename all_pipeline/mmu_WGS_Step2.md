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
--known-sites $mm10/mm10.INDELS.dbSNP142.vcf \
--known-sites $mm10/mm10.dbSNP142.vcf \
-O $BAMOUT/$sample.wes.recal_data.table

$GATK --java-options "-Xmx100G -Djava.io.tmpdir=./" ApplyBQSR \
-R $GENOME -I $BAMOUT/${sample}_fix.bam \
-bqsr $BAMOUT/$sample.wes.recal_data.table \
-O $BAMOUT/${sample}_BQSR.bam 

$GATK --java-options "-Xmx100G -Djava.io.tmpdir=./" Mutect2 \
-R $GENOME \
-I $BAMOUT/${sample}_BQSR.bam \
-O $SNVOUT/$sample.Mutect2only_tumor.vcf.gz ;
done

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
~~~
