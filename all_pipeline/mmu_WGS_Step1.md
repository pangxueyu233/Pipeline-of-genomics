

# Step1. pre-processing of genomics data

This page showed the pre-processing of WGS on mice data, including alignment, PCR duplicates removing, and BQSR.

~~~shell
#for reference,to match the mm10 version, we need to add the 'chr' in each chromosome name in reference files. 
cd /mnt/data/user_data/xiangyu/8922_server/programme/genome_index/GATK/GATK_REFER/bundle/mm10/dbsnp/dbsnp
wget -O ./mm10.INDELS.dbSNP142.vcf.gz.tbi ftp://ftp-mouse.sanger.ac.uk/current_indels/strain_specific_vcfs/C57BL_6NJ.mgp.v5.indels.dbSNP142.normed.vcf.gz.tbi
wget -O ./mm10.dbSNP142.vcf.gz.tbi ftp://ftp-mouse.sanger.ac.uk/current_indels/strain_specific_vcfs/C57BL_6NJ.mgp.v5.snps.dbSNP142.vcf.gz.tbi

awk -v OFS='\t' '{chromosome="chr"$1;gsub($1,chromosome,$1);print $0}' mm10.INDELS.dbSNP142.vcf > chr_mm10.INDELS.dbSNP142.vcf1
awk -v OFS='\t' '{chromosome="chr"$1;gsub($1,chromosome,$1);print $0}' mm10.dbSNP142.vcf > chr_mm10.dbSNP142.vcf1
sed '1,70d' chr_mm10.INDELS.dbSNP142.vcf1 > chr_mm10.INDELS.dbSNP142.vcf.tmp
head -n 70 mm10.INDELS.dbSNP142.vcf > head.tmp.INDEL
cat head.tmp.INDEL chr_mm10.INDELS.dbSNP142.vcf.tmp > chr_mm10.INDELS.dbSNP142.vcf.tmp1
sed '1,69d' chr_mm10.dbSNP142.vcf1 > chr_mm10.dbSNP142.vcf.tmp
head -n 69 mm10.dbSNP142.vcf > head.tmp.dbSNP142
cat head.tmp.dbSNP142 chr_mm10.dbSNP142.vcf.tmp > chr_mm10.dbSNP142.vcf.tmp1
sed 's/##contig=<ID=/##contig=<ID=chr/g' chr_mm10.dbSNP142.vcf.tmp1 > chr_mm10.dbSNP142.vcf
sed 's/##contig=<ID=/##contig=<ID=chr/g' chr_mm10.INDELS.dbSNP142.vcf.tmp1 > chr_mm10.INDELS.dbSNP142.vcf
rm -r chr_mm10.INDELS.dbSNP142.vcf1 chr_mm10.dbSNP142.vcf1 chr_mm10.INDELS.dbSNP142.vcf.tmp head.tmp.INDEL chr_mm10.dbSNP142.vcf.tmp head.tmp.dbSNP142 chr_mm10.dbSNP142.vcf.tmp1 chr_mm10.INDELS.dbSNP142.vcf.tmp1

####build the index file of reference
GATK=/mnt/data/user_data/xiangyu/programme/gatk-4.1.3.0/gatk
$GATK --java-options "-Xmx60G -Djava.io.tmpdir=./" IndexFeatureFile -F chr_mm10.dbSNP142.vcf
$GATK --java-options "-Xmx60G -Djava.io.tmpdir=./" IndexFeatureFile -F chr_mm10.INDELS.dbSNP142.vcf


####reference and software
mm10=/mnt/data/user_data/xiangyu/8922_server/programme/genome_index/GATK/GATK_REFER/bundle/mm10/dbsnp/dbsnp
GENOME=/mnt/data/user_data/xiangyu/8922_server/programme/genome_index/bwa_index/bwa_mm10_index/genome.fa
GATK=/mnt/data/user_data/xiangyu/programme/gatk-4.1.3.0/gatk

####outpath for genomics
BAMOUT=/mnt/data/userdata/abao/project/4_whole_genome_sequencing/1_Lishujun/other/
SNVOUT=$BAMOUT/snv_out
mkdir $SNVOUT

####sample list
vim samplelist
33_2R /mnt/data/sequencedata/whole-genome_sequencing/Mki67/ANYWSC180111_PM-YWSC180111-01_AH5MFMDSXX_2018-12-01/Cleandata/33-2R/33-2R
33_B /mnt/data/sequencedata/whole-genome_sequencing/Mki67/ANYWSC180111_PM-YWSC180111-01_AH5MFMDSXX_2018-12-01/Cleandata/33-B/33-B
33_L /mnt/data/sequencedata/whole-genome_sequencing/Mki67/ANYWSC180111_PM-YWSC180111-01_AH5MFMDSXX_2018-12-01/Cleandata/33-L/33-L
35_B /mnt/data/sequencedata/whole-genome_sequencing/Mki67/ANYWSC180111_PM-YWSC180111-01_AH5MFMDSXX_2018-12-01/Cleandata/35-B/35-B
35_L /mnt/data/sequencedata/whole-genome_sequencing/Mki67/ANYWSC180111_PM-YWSC180111-01_AH5MFMDSXX_2018-12-01/Cleandata/35-L/35-L
35_N /mnt/data/sequencedata/whole-genome_sequencing/Mki67/ANYWSC180111_PM-YWSC180111-01_AH5MFMDSXX_2018-12-01/Cleandata/35-N/35-N

####execute codes
cat samplelist | while read id ; do
arr=($id)
fq2=${arr[1]}'_R1.fq.gz'
fq1=${arr[1]}'_R2.fq.gz'
sample=${arr[0]}
echo ${fq2}
echo ${fq1}
echo $sample
bwa mem -t 30 -M  -R "@RG\tID:$sample\tSM:$sample\tLB:WGS\tPL:Illumina" $GENOME $fq1 $fq2 > $BAMOUT/$sample.sam
$GATK  --java-options "-Xmx60G -Djava.io.tmpdir=./"  SortSam -SO coordinate  -I $BAMOUT/$sample.sam -O $BAMOUT/$sample.sort.bam > $BAMOUT/$sample.log.sort	
$GATK  --java-options "-Xmx60G -Djava.io.tmpdir=./"  MarkDuplicates \
-I $BAMOUT/$sample.sort.bam \
-O $BAMOUT/${sample}_marked.bam \
-REMOVE_DUPLICATES=true \
-M $BAMOUT/$sample.metrics 1> $BAMOUT/$sample.log.mark  2>&1 
rm -r $BAMOUT/$sample.sam ;
done
~~~

