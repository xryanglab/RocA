# -*- coding:UTF-8 -*-
'''
Author: Li Fajin
Date: 2020-12-21 19:59:09
LastEditors: Li Fajin
LastEditTime: 2022-04-01 14:26:35
Description: snakemake pipeline for ribo-seq data analyses, just for basic control steps such fastqc, cutadapt, qfiltering, mapping, samtools sort and index, 3nt periodicity checking, et al.
'''

import os

## samples defination
SAMPLES=["SRR5937640",
         "SRR5937641",
         "SRR5937642",
         "SRR5937643",
         "SRR5937644",
	 "SRR5937645"]

## parameters defination
BOWTIE_NONCODING_INDEX="/workdata/home/lifj/lifj/data/Reference/human/rRNA_bowtieIndex/human_rRNA"
STAR_GENOME_INDEX="/workdata/home/lifj/lifj/data/Reference/human/STAR_Human_Ensembl_GRCh38_Ensembl"
RIBOCODE_ANNOT="/workdata/home/lifj/lifj/data/Reference/human/RiboCode_annot_Human/RiboCode_annot"
GTFFILE="/workdata/home/lifj/lifj/data/Reference/human/Homo_sapiens.GRCh38.88.gtf"
ADAPTER="CTGTAGGCACCATCAAT"


## softwares defination
# make sure the all softwares needed are all installed and in your PATH

with os.popen("which fastqc") as path:
    FASTQC = path.read().strip()
with os.popen("which cutadapt") as path:
    CUTADAPT = path.read().strip()
with os.popen("which fastq_quality_filter") as path:
    FILTER = path.read().strip()
with os.popen("which bowtie") as path:
    BOWTIE = path.read().strip()
with os.popen("which STAR") as path:
    STAR = path.read().strip()
with os.popen("which samtools") as path:
    SAMTOOLS = path.read().strip()
with os.popen("which metaplots") as path:
    METAPLOTS = path.read().strip()
with os.popen("which LengthDistribution") as path:
    LENGTHDISTRIBUTION = path.read().strip()
with os.popen("which ModifyHTseq") as path:
    MODIFYHTSEQ = path.read().strip()
with os.popen("which StatisticReadsOnDNAsContam") as path:
    STATISTICREADSONDNASCONTAM = path.read().strip()
with os.popen("which bamCoverage") as path:
    BAMCOVERAGE = path.read().strip()


## snakemakes pipeline

# ruleorder:statistic_contamination>summary>mergeLogs

rule all:
    input:
        expand("01.beforeQC/{sample}",sample=SAMPLES),
        expand("02.cutadapt/{sample}_trimmed.fastq",sample=SAMPLES),
        expand("03.filter/{sample}_trimmed_Qfilter.fastq",sample=SAMPLES),
        expand("04.afterQC/{sample}",sample=SAMPLES),
        expand("05.contam/noncontam_{sample}.fastq",sample=SAMPLES),
        expand("06.finalQC/{sample}",sample=SAMPLES),
        expand("07.STAR/{sample}_STAR/{sample}.Aligned.sortedByCoord.out.bam",sample=SAMPLES),
        expand("07.STAR/{sample}_STAR/{sample}.Aligned.toTranscriptome.out.bam",sample=SAMPLES),
        expand("07.STAR/{sample}_STAR/{sample}.Log.final.out",sample=SAMPLES),
        expand("07.STAR/{sample}_STAR/{sample}.Aligned.toTranscriptome.out.sorted.bam",sample=SAMPLES),
        expand("07.STAR/{sample}_STAR/{sample}.Aligned.{type}.bam.bai",sample=SAMPLES,type={"sortedByCoord.out","toTranscriptome.out.sorted"}),
        #expand("07.STAR/{sample}_STAR/{sample}.bw",sample=SAMPLES),
        expand("08.periodicity/{sample}{sample}.Aligned.toTranscriptome.out.pdf",sample=SAMPLES),
        expand("08.periodicity/{sample}_pre_config.txt",sample=SAMPLES),
        expand("09.LengthDistribution/{sample}_reads_length.pdf",sample=SAMPLES),
        expand("09.LengthDistribution/{sample}_reads_length.txt",sample=SAMPLES),
        expand("10.read_counts/{sample}.counts",sample=SAMPLES),
        expand("11.RPF_statistics/{sample}_DNA.pdf",sample=SAMPLES),
        expand("11.RPF_statistics/{sample}_Intron.pdf",sample=SAMPLES),
        expand("11.RPF_statistics/{sample}_RNA.pdf",sample=SAMPLES),
        expand("11.RPF_statistics/{sample}_reads_distribution.txt",sample=SAMPLES),
        expand("12.summary/{sample}_Cutadapt.txt",sample=SAMPLES),
        expand("12.summary/{sample}_Filtering.txt",sample=SAMPLES),
        expand("12.summary/{sample}_RRNA_contam.txt",sample=SAMPLES),
        expand("12.summary/{sample}_Star_mapping.txt",sample=SAMPLES),
        expand("12.summary/{sample}_Statistics.txt",sample=SAMPLES)

rule beforeQC:
    input:
        "00.rawdata/{sample}.fastq"
    output:
        directory("01.beforeQC/{sample}")
    threads: 1
    shell:
        """
        mkdir -p {output}
        {FASTQC} {input} -o {output}
        """

rule Cutadapt:
    input:
        "00.rawdata/{sample}.fastq"
    output:
        "02.cutadapt/{sample}_trimmed.fastq"
    log:
        "02.cutadapt/{sample}_trimmed.log"
    threads: 1
    shell:
        """
        {CUTADAPT} -m 15 -M 45 --match-read-wildcards -a {ADAPTER} -o {output} {input} > {log} 2>&1
        """

rule quality_filter:
    input:
        "02.cutadapt/{sample}_trimmed.fastq"
    output:
        "03.filter/{sample}_trimmed_Qfilter.fastq"
    log:
        "03.filter/{sample}_trimmed_Qfilter.log"
    threads: 1
    shell:
        """
        {FILTER} -Q33 -v -q 25 -p 75 -i {input} -o {output} > {log} 2>&1
        """

rule afterQC:
    input:
        "03.filter/{sample}_trimmed_Qfilter.fastq"
    output:
        directory("04.afterQC/{sample}")
    threads: 1
    shell:
        """
        mkdir -p {output}
        {FASTQC} {input} -o {output}
        """

rule remove_rRNA:
    input:
        "03.filter/{sample}_trimmed_Qfilter.fastq"
    output:
        "05.contam/{sample}.aln",
        "05.contam/noncontam_{sample}.fastq"
    log:
        "05.contam/{sample}.log"
    threads: 1
    shell:
        """
        {BOWTIE} -n 0 -y -a --norc --best --strata -S -p 4 -l 15 --un={output[1]} {BOWTIE_NONCODING_INDEX} -q {input} {output[0]} > {log} 2>&1
        """


rule finalQC:
    input:
        "05.contam/noncontam_{sample}.fastq"
    output:
        directory("06.finalQC/{sample}")
    threads: 1
    shell:
        """
        mkdir -p {output}
        {FASTQC} {input} -o {output}
        """

rule STAR:
    input:
        "05.contam/noncontam_{sample}.fastq"
    output:
        "07.STAR/{sample}_STAR/{sample}.Aligned.sortedByCoord.out.bam",
        "07.STAR/{sample}_STAR/{sample}.Aligned.toTranscriptome.out.bam",
        "07.STAR/{sample}_STAR/{sample}.Log.final.out",
        "07.STAR/{sample}_STAR/{sample}.Log.out",
        "07.STAR/{sample}_STAR/{sample}.Log.progress.out",
        "07.STAR/{sample}_STAR/{sample}.ReadsPerGene.out.tab",
        "07.STAR/{sample}_STAR/{sample}.Signal.UniqueMultiple.str1.out.wig",
        "07.STAR/{sample}_STAR/{sample}.Signal.UniqueMultiple.str2.out.wig",
        "07.STAR/{sample}_STAR/{sample}.Signal.Unique.str1.out.wig",
        "07.STAR/{sample}_STAR/{sample}.Signal.Unique.str2.out.wig",
        "07.STAR/{sample}_STAR/{sample}.SJ.out.tab"
    threads: 4
    shell:
        """
        ## part of parameters from GSE125114
        {STAR} --runThreadN {threads} --outFilterType Normal --outWigType wiggle --outWigStrand Stranded \
            --outWigNorm RPM --alignEndsType EndToEnd --outFilterMismatchNmax 5 \
            --outFilterMultimapNmax 1  --outFilterMatchNmin 16 --genomeDir {STAR_GENOME_INDEX} --readFilesIn {input} \
            --outFileNamePrefix  07.STAR/{wildcards.sample}_STAR/{wildcards.sample}. \
            --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts \
            --outSAMattributes All
        """

rule samtools_sort:
    input:
        "07.STAR/{sample}_STAR/{sample}.Aligned.toTranscriptome.out.bam"
    output:
        "07.STAR/{sample}_STAR/{sample}.Aligned.toTranscriptome.out.sorted.bam"
    threads: 1
    shell:
        "{SAMTOOLS} sort -T {output}.temp -o {output} {input}"


rule samtools_index:
    input:
        "07.STAR/{sample}_STAR/{sample}.Aligned.{type}.bam"
    output:
        "07.STAR/{sample}_STAR/{sample}.Aligned.{type}.bam.bai"
    threads: 1
    shell:
        "{SAMTOOLS} index {input}"

#rule bam2bigwig:
#    input:
#        "07.STAR/{sample}_STAR/{sample}.Aligned.sortedByCoord.out.bam"
#    output:
#        "07.STAR/{sample}_STAR/{sample}.bw"
#    log:
#        "07.STAR/{sample}_STAR/{sample}.bam2bw.log"
#    threads: 2
#    shell:
#        """
#        {BAMCOVERAGE} -b {input} -o {output} --normalizeUsing CPM --binSize 1 > {log} 2>&1
#        """

rule periodicity:
    input:
        "07.STAR/{sample}_STAR/{sample}.Aligned.toTranscriptome.out.bam"
    output:
        "08.periodicity/{sample}{sample}.Aligned.toTranscriptome.out.pdf",
        "08.periodicity/{sample}_pre_config.txt"
    threads: 1
    shell:
        "{METAPLOTS} -a {RIBOCODE_ANNOT} -r {input} -o 08.periodicity/{wildcards.sample} -m 15 -M 45"

rule length_distribution:
    input:
        "07.STAR/{sample}_STAR/{sample}.Aligned.toTranscriptome.out.sorted.bam"
    output:
        "09.LengthDistribution/{sample}_reads_length.pdf",
        "09.LengthDistribution/{sample}_reads_length.txt"
    log:
        "09.LengthDistribution/{sample}_reads_length.log"
    threads: 1
    shell:
        "{LENGTHDISTRIBUTION} -i {input} -o 09.LengthDistribution/{wildcards.sample} -f bam > {log} 2>&1"


# specifically for Ribo-seq data, if you have RNA-seq data, using htseq-count please.
rule htseq_count:
    input:
        "07.STAR/{sample}_STAR/{sample}.Aligned.sortedByCoord.out.bam"
    output:
        "10.read_counts/{sample}.counts"
    log:
        "10.read_counts/{sample}.counts.log"
    threads: 1
    shell:
        """
        {MODIFYHTSEQ} -i {input} -g {GTFFILE} -o {output} -t CDS -m union -q 10 --minLen 15 --maxLen 45 --exclude-first 45 --exclude-last 15 --id-type gene_id > {log} 2>&1
        """

rule statistic_contamination:
    input:
        "07.STAR/{sample}_STAR/{sample}.Aligned.sortedByCoord.out.bam"
    output:
        "11.RPF_statistics/{sample}_DNA.pdf",
        "11.RPF_statistics/{sample}_Intron.pdf",
        "11.RPF_statistics/{sample}_RNA.pdf",
        "11.RPF_statistics/{sample}_reads_distribution.txt"
    log:
        "11.RPF_statistics/{sample}.statistics.log"
    threads: 1
    shell:
        """
        {STATISTICREADSONDNASCONTAM} -i  {input}  -g {GTFFILE}  -o  11.RPF_statistics/{wildcards.sample} > {log} 2>&1
        """


rule summary:
    input:
        cutadapt="02.cutadapt/{sample}_trimmed.log",
        filtering="03.filter/{sample}_trimmed_Qfilter.log",
        removeRNAcontam="05.contam/{sample}.log",
        starMapping="07.STAR/{sample}_STAR/{sample}.Log.final.out",
        statistics="11.RPF_statistics/{sample}_reads_distribution.txt"
    output:
        cutadaptLog="12.summary/{sample}_Cutadapt.txt",
        filteringLog="12.summary/{sample}_Filtering.txt",
        rRNAContamLog="12.summary/{sample}_RRNA_contam.txt",
        starMappingLog="12.summary/{sample}_Star_mapping.txt",
        statisticsLog="12.summary/{sample}_Statistics.txt"
        # finalLog="RNA_seq_summary_statistics.txt"

    shell:
        """
        python summary.py -f {input.filtering}  --of {output.filteringLog} -c {input.cutadapt}  --oc {output.cutadaptLog} \
            -r {input.removeRNAcontam} --or {output.rRNAContamLog} -s {input.statistics} --os {output.statisticsLog}\
            -m {input.starMapping} --om {output.starMappingLog}
        """
