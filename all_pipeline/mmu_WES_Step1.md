

# Step1. pre-processing of genomics data

This page showed the pre-processing of WGS on mice data, including alignment, PCR duplicates removing, and BQSR.

- *Here, we had specialize the position of `exon bed files` of mm10 refernce.*

~~~shell
####reference and software
mm10=/mnt/data/user_data/xiangyu/8922_server/programme/genome_index/GATK/GATK_REFER/bundle/mm10/dbsnp/dbsnp
GENOME=/mnt/data/user_data/xiangyu/8922_server/programme/genome_index/bwa_index/bwa_mm10_index/genome.fa
GATK=/mnt/data/user_data/xiangyu/programme/gatk-4.1.3.0/gatk
EXON_BED=/mnt/data/user_data/xiangyu/programme/genome_index/reference_data/reference_data/mm10/bed/Exons.with_genes.bed

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

