# Pipeline of genomics on mice data

​	Based on the current study, there are really barely analysis of `WGS` and `WES` on mice data. Although, [MoCaSeq](https://github.com/roland-rad-lab/MoCaSeq) had provided a systematic analysis pipelines on mice data, it almost based on the `GRcm38` version genome reference. If you used (or plan to use) the `mm10` as the references to process, that's  not easy  to make it on  [MoCaSeq](https://github.com/roland-rad-lab/MoCaSeq). Here, we provide a pipeline of genomics analysis on mice data based on the `mm10`, and we plan to make it much more flexible and easier to  understand each step in processing. 

Now, let's begin our pipeline of **whole genome sequencing (WGS)** data analysis on mice data

- ### Step1. pre-processing of  WGS data

  You could access the [Step1](all_pipeline/mmu_WGS_Step1.md) by clicking [here](all_pipeline/mmu_WGS_Step1.md)

- ### Step2. somatic SNVs detcted in each sample

  You could access the [Step2](all_pipeline/mmu_WGS_Step2.md) by clicking [here](all_pipeline/mmu_WGS_Step2.md)

- ### Step3. somatic CNV detcted in each sample

- You could access the [Step3](all_pipeline/mmu_WGS_Step3.md) by clicking [here](all_pipeline/mmu_WGS_Step3.md)

And we also share the visualization codes about the `Step2` results, and you also could access by clicking [here](

---

Now, let's begin our pipeline of **whole exon sequencing (WES/WXS)** data analysis on mice data

- There is no differences of `Step1` and `BQSR` in `WGS` and `WES` analysis, but you could specialize the position of `exon bed files` in [Step1](all_pipeline/mmu_WES_Step1.md), which you could follow by clicking [here](all_pipeline/mmu_WES_Step1.md)

- Besides, you just need to add the `-L` parameter into the pipelines of SNV and CNV detected to accelerate the processing. 

- ### somatic SNVs detcted in each sample

  You could access the [Step2](all_pipeline/mmu_WES_Step2.md) by clicking [here](all_pipeline/mmu_WES_Step2.md)

- ### somatic CNV detcted in each sample

- You could access the [Step3](all_pipeline/mmu_WES_Step3.md) by clicking [here](all_pipeline/mmu_WES_Step3.md)



