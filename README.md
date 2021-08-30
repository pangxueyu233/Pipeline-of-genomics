# Pipeline-of-genomics
This workshop recorded the whole processing steps of genomics analysis in CC-LY Lab, including alignment, quality control, somatic and germline SNV detected and CNV detected from fasta files. There were two individual pipelines to be executed for mouse and human data, respectively.  Here, we referred the GATK workshop for SNV calling, CNVkit workshop for CNV detecting. 

The software we need: 

~~~shell
#bwa
Program: bwa (alignment via Burrows-Wheeler transformation)
Version: 0.7.17-r1188
Contact: Heng Li <lh3@sanger.ac.uk>

#gatk4.1
GATK=/mnt/data/user_data/xiangyu/programme/gatk-4.1.3.0/gatk

#cnvkit
cnvkit.py version
0.9.9

#vep
/mnt/data/userdata/xiangyu/programme/vcf2maf-master/vcf2maf.pl

#Rscript
Rscript --version
R scripting front-end version 3.5.1 (2018-07-02)

#perl
perl --version
This is perl 5, version 22, subversion 1 (v5.22.1) built for x86_64-linux-gnu-thread-multi
(with 77 registered patches, see perl -V for more detail)
~~~


