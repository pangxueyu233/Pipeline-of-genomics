# step1 . Alignment

This page showed the alignment of WGS for mouse data.

~~~shell
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

# step2 . SNV detected by Mutect2

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

# step3. SNV annotated by vep and maf files generated

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

# step4. CNV called by cnvkit (for individual samples)

~~~shell
####reference and software
refFlat=/mnt/data/userdata/xiangyu/programme/genome_index/bwa_index/bwa_mm10_index/refFlat.txt
GENOME=/mnt/data/userdata/xiangyu/programme/genome_index/bwa_index/bwa_mm10_index/genome.fa
extend_tools=/mnt/data/userdata/xiangyu/programme/R_PACKAGES/my_code/CNV_pipeline/
CNVKIT=/mnt/data/userdata/xiangyu/programme/cnvkit-0.9.9/cnvkit.py
PROCESSES=50
parameter='-1.3,-0.4,0.3,0.9'
ratio=0.4

####outpath for genomics
BAMOUT=/mnt/data/userdata/abao/project/4_whole_genome_sequencing/1_Lishujun/other
CNVOUT=$BAMOUT/cnvkit_out
mkdir $CNVOUT
rm_XY_out=$CNVOUT/rm_XY
all_chr_out=$CNVOUT/all_chr
mkdir $rm_XY_out
mkdir $all_chr_out

cat samplelist | while read id ; do
arr=($id)
sample=${arr[0]}
echo $sample
####execute codes
##=============== call batch to WGS ===============##
$CNVKIT batch $BAMOUT/${sample}_marked.bam -m wgs -f $GENOME -n -p $PROCESSES --output-reference $CNVOUT/$sample.cnn \
--segment-method hmm --target-avg-size 1000 --drop-low-coverage --output-dir $CNVOUT --annotate $refFlat
# default bin >10kb
# cutoff probe >10
# chromosome      start   end     gene    log2    depth   probes  weight
perl -ne '@F=split/\t/; $probe=$F[6]; $len=$F[2]-$F[1]; if($.==1){print}elsif($probe>10 && $len>1000){print}' $CNVOUT/${sample}_marked.call.cns > $CNVOUT/$sample.filtered.cns ;
done
~~~

## 4.1 The coverage of X/Y always too slow to detect, consequently, following analysis had removed the XY chromosome for flowing analysis 

~~~shell
####execute codes
cat samplelist | while read id ; do
arr=($id)
sample=${arr[0]}
echo $sample
# chr1-22 (XY also removed because it had low accuracy)
perl -ne 'print if (/^chr\d+/ || /^chromosome/ || /^\d+/)' $CNVOUT/${sample}_marked.cnr > $rm_XY_out/$sample.rmXY.cnr
perl -ne 'print if (/^chr\d+/ || /^chromosome/ || /^\d+/)' $CNVOUT/$sample.filtered.cns > $rm_XY_out/$sample.rmXY.cns
##=============== Add absolute CN ===============##
$CNVKIT call $rm_XY_out/$sample.rmXY.cns -m threshold -t=${parameter} -o $rm_XY_out/$sample.rmXY.call.cns
##=============== Plot cn pattern ===============##
# *.rmXY.cns      # all of cnv is yellow color
# *.rmXY.call.cns   # only abnormal cnv is yellow color
$CNVKIT scatter $rm_XY_out/$sample.rmXY.cnr -s $rm_XY_out/$sample.rmXY.call.cns --y-min -4 --y-max 4 -w 10000000 -o $rm_XY_out/$sample.rmXY.call.cns.pdf
##=============== Gene-level ===============##
# call genelevel cnv
#https://cnvkit.readthedocs.io/en/stable/reports.html
# use *.call.rmXY.cns
$CNVKIT genemetrics $rm_XY_out/$sample.rmXY.cnr > $rm_XY_out/$sample.ratio_gene.rmXY.raw.tsv
$CNVKIT genemetrics $rm_XY_out/$sample.rmXY.cnr -s $rm_XY_out/$sample.rmXY.call.cns -t ${ratio} -m 5 > $rm_XY_out/$sample.segment_gene.rmXY.raw.tsv
sed '1d' $rm_XY_out/$sample.ratio_gene.rmXY.raw.tsv | cut -f 1 | sort -u > $rm_XY_out/$sample.ratio_gene.gene
sed '1d' $rm_XY_out/$sample.segment_gene.rmXY.raw.tsv | cut -f 1 | sort -u > $rm_XY_out/$sample.segment_gene.gene
comm -12 $rm_XY_out/$sample.ratio_gene.gene $rm_XY_out/$sample.segment_gene.gene > $rm_XY_out/$sample.trusted-genes.txt
Rscript $extend_tools/extract_and_filter_Scripts --l=$rm_XY_out/$sample.trusted-genes.txt --m=$rm_XY_out/$sample.segment_gene.rmXY.raw.tsv --c=gene --o=$rm_XY_out/$sample.segment_gene.trusted.gainloss.rmXY.tsv
# remove process files
rm -f $rm_XY_out/genome.*bed $rm_XY_out/$sample.*targetcoverage.cnn $rm_XY_out/$sample.ratio_gene.gene $rm_XY_out/$sample.segment_gene.gene $rm_XY_out/$sample.*.tmp
rm -f $rm_XY_out/$sample.*.filtered.cns $rm_XY_out/$sample.*.ratio_gene.trusted.gainloss.tsv $rm_XY_out/$sample.*trusted-genes.txt $rm_XY_out/$sample.segment_gene.rmXY.raw.tsv ;
done
~~~

## 4.2 This part includes all chromosome data for flowing analysis 

~~~shell
####execute codes
cat samplelist | while read id ; do
arr=($id)
sample=${arr[0]}
echo $sample
# chr1-22 (KEEP XY)
##=============== Add absolute CN ===============##
$CNVKIT call $CNVOUT/$sample.filtered.cns -m threshold -t=${parameter} -o $all_chr_out/$sample.allChr.call.cns
##=============== Plot cn pattern ===============##
# *.allChr.cns      # all of cnv is yellow color
# *.allChr.call.cns   # only abnormal cnv is yellow color
$CNVKIT scatter $CNVOUT/${sample}_marked.cnr -s $all_chr_out/$sample.allChr.call.cns --y-min -4 --y-max 4 -w 10000000 -o $all_chr_out/$sample.allChr.call.cns.pdf
##=============== Gene-level ===============##
# call genelevel cnv
#https://cnvkit.readthedocs.io/en/stable/reports.html
# use *.call.allChr.cns
$CNVKIT genemetrics $CNVOUT/${sample}_marked.cnr > $all_chr_out/$sample.ratio_gene.allChr.raw.tsv
$CNVKIT genemetrics $CNVOUT/${sample}_marked.cnr -s $all_chr_out/$sample.allChr.call.cns -t ${ratio} -m 5 > $all_chr_out/$sample.segment_gene.allChr.raw.tsv
sed '1d' $all_chr_out/$sample.ratio_gene.allChr.raw.tsv | cut -f 1 | sort -u > $all_chr_out/$sample.ratio_gene.gene
sed '1d' $all_chr_out/$sample.segment_gene.allChr.raw.tsv | cut -f 1 | sort -u > $all_chr_out/$sample.segment_gene.gene
comm -12 $all_chr_out/$sample.ratio_gene.gene $all_chr_out/$sample.segment_gene.gene > $all_chr_out/$sample.trusted-genes.txt
Rscript $extend_tools/extract_and_filter_Scripts --l=$all_chr_out/$sample.trusted-genes.txt --m=$all_chr_out/$sample.segment_gene.allChr.raw.tsv --c=gene --o=$all_chr_out/$sample.segment_gene.trusted.gainloss.allChr.tsv
# remove process files
rm -f $all_chr_out/genome.*bed $all_chr_out/$sample.*targetcoverage.cnn $all_chr_out/$sample.ratio_gene.gene $all_chr_out/$sample.segment_gene.gene $all_chr_out/$sample.*.tmp
rm -f $all_chr_out/$sample.*.filtered.cns $all_chr_out/$sample.*.ratio_gene.trusted.gainloss.tsv $all_chr_out/$sample.*trusted-genes.txt $all_chr_out/$sample.segment_gene.allChr.raw.tsv ;
done
~~~

## 4.3 merge all individual results

~~~
~~~

