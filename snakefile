import os
from snakemake.logging import logger
import pandas as pd
import numpy as np
import re

#configfile: "/home/fengchr/staff/RNA_edit_pipeline/configs/config_E_colli.yaml"
mapping=pd.read_csv(config['mapping'])
test_controls=mapping[mapping['sample_type']=='test']
normals=mapping[mapping['use_as_pon']=='yes']

ruleorder: mutect2 > mutect2_tumor_only
ruleorder: STAR > STAR_alt

rule all:
    input:
        "success"

rule generate_STAR_ref:
    input:
        genome_fasta=config['reference'],
        annotation_gtf=config['annotation_gtf']
    output:
        "reference/STAR_ref/chrName.txt"
    threads: 64
    params:
        genomeSAindexNbases=config['STAR']['genomeSAindexNbases'],
        sjdbOverhang=config['STAR']['sjdbOverhang']
    shell:
        """
        source activate STAR
        STAR --runThreadN {threads} \
        --runMode genomeGenerate \
        --genomeDir reference/STAR_ref \
        --genomeFastaFiles {input.genome_fasta} \
        --sjdbGTFfile {input.annotation_gtf} \
        --sjdbOverhang {params.sjdbOverhang} \
        --genomeSAindexNbases {params.genomeSAindexNbases}
        """

def get_fastq(wildcards):
    try:
        fq1,fq2 = list(mapping['file'][mapping['sample_id']==wildcards.sample])[0].split(';')
        return({"fq1":fq1, "fq2":fq2})
    except:
        return({"fq1":'', "fq2":''})

rule STAR:
    input:
        config['STAR']["STARref"]+'/chrName.txt',
        unpack(get_fastq)
    output:
        "{sample}_STAR/Aligned.out.sam"
    log:
        "log/{sample}.STAR.log.txt"
    threads: 64
    params:
        genomeDir=config['STAR']["STARref"],
        outFilterMultimapNmax=config['STAR']['outFilterMultimapNmax'],
        sjdbOverhang=config['STAR']['sjdbOverhang'],
        readcommand=config['readcommand']
    shell:
        """
        source activate STAR
        STAR --runThreadN {threads} \
        --genomeDir {params.genomeDir} \
        --readFilesIn {input.fq1} {input.fq2} {params.readcommand}\
        --outSAMmapqUnique 60 \
        --outFilterMultimapNmax {params.outFilterMultimapNmax} \
        --sjdbOverhang {params.sjdbOverhang} \
        --outFileNamePrefix {wildcards.sample}_STAR/ \
        > {log} 2>&1
        conda deactivate
        """

rule STAR_alt:
    input:
        "reference/STAR_ref/chrName.txt",
        unpack(get_fastq)
    output:
        "{sample}_STAR/Aligned.out.sam"
    log:
        "log/{sample}.STAR.log.txt"
    threads: 64
    params:
        genomeDir="reference/STAR_ref",
        outFilterMultimapNmax=config['STAR']['outFilterMultimapNmax'],
        sjdbOverhang=config['STAR']['sjdbOverhang'],
        readcommand=config['readcommand']
    shell:
        """
        source activate STAR
        STAR --runThreadN {threads} \
        --genomeDir {params.genomeDir} \
        --readFilesIn {input.fq1} {input.fq2} {params.readcommand}\
        --outSAMmapqUnique 60 \
        --outFilterMultimapNmax {params.outFilterMultimapNmax} \
        --sjdbOverhang {params.sjdbOverhang} \
        --outFileNamePrefix {wildcards.sample}_STAR/ \
        > {log} 2>&1
        conda deactivate
        """

rule samtools_sort:
    input:
        "{sample}.sam"
    output:
        "{sample}.sorted.bam"
    threads: 64
    shell:
        "samtools sort {input} -o {output} -@ {threads}; rm {input}"

rule mark_dup:
    input:
        "{sample}_STAR/Aligned.out.sorted.bam"
    output:
        "{sample}_STAR/Aligned.out.sorted.dupmkd.bam"
    log:
        "log/{sample}.mark_dup.log.txt"
    threads: 32
    shell:
        """
        java -jar ~/software/picard/build/libs/picard.jar MarkDuplicates -I {input} -M {wildcards.sample}_STAR/picard.dup.metrics.txt -O {output} > {log} 2>&1
        """

rule add_GP:
    input:
        "{sample}_STAR/Aligned.out.sorted.dupmkd.bam"
    output:
        "{sample}_STAR/Aligned.out.sorted.dupmkd.addGP.bam"
    params:
        RGLB="{sample}",
        RGPL="ILLUMINA",
        RGPU="{sample}",
        RGSM="{sample}"
    log:
        "log/{sample}.addGP.log.txt"
    shell:
        """
        gatk AddOrReplaceReadGroups -I {input} \
        -O {output} \
        --RGLB {params.RGLB} \
        --RGPL {params.RGPL} \
        --RGPU {params.RGPU} \
        --RGSM {params.RGSM} > {log} 2>&1
        """

rule samtools_index:
    input:
        "{sample}_STAR/Aligned.out.sorted.dupmkd.addGP.bam"
    output:
        "{sample}_STAR/Aligned.out.sorted.dupmkd.addGP.bam.bai"
    threads: 32
    shell:
        "samtools index -@ {threads} {input}"

rule index_ref:
    input:
        config['reference']
    output:
        config['reference']+".fai",
        re.sub("\.[^\./]*$","",config['reference'])+".dict"
    shell:
        """
        samtools faidx {input}
        gatk CreateSequenceDictionary -R {input}
        """

rule pon_mutect2:
    input:
        bam="{sample}_STAR/Aligned.out.sorted.dupmkd.addGP.bam",
        index="{sample}_STAR/Aligned.out.sorted.dupmkd.addGP.bam.bai",
        ref_fai=config['reference']+".fai",
        ref_dict=re.sub("\.[^\./]*$","",config['reference'])+".dict"
    output:
        "pon/{sample}_normal.vcf.gz"
    log:
        "log/{sample}.pon_mutect2.log.txt"
    params:
        ref=config['reference']
    threads: 64
    shell:
        """
        gatk Mutect2 -R {params.ref} \
        -I {input.bam} \
        --tumor-sample {wildcards.sample} \
        -O {output} \
        --max-mnp-distance 0 \
        > {log} 2>&1
        """

rule GenomicsDBImport:
    input:
        expand("pon/{i}_normal.vcf.gz", i=list(normals['sample_id']))
    output:
        "pon_db/callset.json"
    log:
        "log/GenomicsDBImport.log.txt"
    params:
        ref=config['reference'],
        intervals_list=config['intervals_list']
    shell:
        """
        rm -rf pon_db
        gatk GenomicsDBImport -R {params.ref} \
        -V `echo "{input}"|sed 's/ / -V /g'` \
        --genomicsdb-workspace-path pon_db \
        -L {params.intervals_list} > {log} 2>&1
        """

rule create_pon:
    input:
        "pon_db/callset.json"
    output:
        "pon.vcf"
    log:
        "log/create_pon.log.txt"
    params:
        ref=config['reference'],
        min_sample_count=config['min_sample_count']
    shell:
        "gatk CreateSomaticPanelOfNormals -R {params.ref} -V gendb://pon_db --min-sample-count {params.min_sample_count} -O {output} > {log} 2>&1"

def get_control_s(sample):
    try:
        return(list(test_controls['control_sample'][test_controls['sample_id']==sample])[0].split(','))
    except:
        return('')

rule mutect2:
    input:
        test="{sample}_STAR/Aligned.out.sorted.dupmkd.addGP.bam",
        index="{sample}_STAR/Aligned.out.sorted.dupmkd.addGP.bam.bai",
        control=lambda wildcards: expand("{controls}_STAR/Aligned.out.sorted.dupmkd.addGP.bam", controls=get_control_s(f"{wildcards.sample}")),
        control_index=lambda wildcards: expand("{controls}_STAR/Aligned.out.sorted.dupmkd.addGP.bam.bai", controls=get_control_s(f"{wildcards.sample}")),
        pon="pon.vcf",
        ref_fai=config['reference']+".fai",
        ref_dict=re.sub("\.[^\./]*$","",config['reference'])+".dict"
    output:
        vcf="{sample}.mutect2.vcf",
        stats="{sample}.mutect2.vcf.stats"
    log:
        "log/{sample}.mutect2.log.txt"
    params:
        ref=config['reference'],
        control_s=lambda wildcards: get_control_s(f"{wildcards.sample}")
    threads: 64
    shell:
        """
        gatk Mutect2 -R {params.ref} \
        -I {input.test} \
        -I `echo "{input.control}"|sed 's/ / -I /g'` \
        -normal `echo "{params.control_s}"|sed 's/ / -normal /g'` \
        -tumor {wildcards.sample} \
        -pon {input.pon} \
        -O {output.vcf} \
        -bamout {wildcards.sample}.mutect2.bam \
        > {log} 2>&1
        """

rule mutect2_tumor_only:
    input:
        test="{sample}_STAR/Aligned.out.sorted.dupmkd.addGP.bam",
        index="{sample}_STAR/Aligned.out.sorted.dupmkd.addGP.bam.bai",
        pon="pon.vcf"
    output:
        vcf="{sample}.mutect2.vcf",
        stats="{sample}.mutect2.vcf.stats"
    log:
        "log/{sample}.mutect2.log.txt"
    params:
        ref=config['reference']
    threads: 64
    shell:
        """
        gatk Mutect2 -R {params.ref} \
        -I {input.test} \
        -pon {input.pon} \
        -O {output.vcf} \
        -bamout {wildcards.sample}.mutect2.bam \
        > {log} 2>&1
        """

rule filter:
    input:
        vcf="{sample}.mutect2.vcf",
        stats="{sample}.mutect2.vcf.stats"
    output:
        "{sample}.mutect2.filtered.vcf"
    log:
        "log/{sample}.mutect2.filter.log.txt"
    params:
        ref=config['reference']
    shell:
        "gatk FilterMutectCalls -R {params.ref} -V {input.vcf} --stats {input.stats} -O {output} > {log} 2>&1"

rule select_snp:
    input:
        "{sample}.mutect2.filtered.vcf"
    output:
        "{sample}.mutect2.filtered.subs.vcf"
    params:
        ref=config['reference']
    shell:
        "gatk SelectVariants -V {input} -O {output} -select-type SNP -R {params.ref}"

rule extract_snp_site:
    input:
        "{sample}.mutect2.filtered.subs.vcf"
    output:
        "{sample}.mutect2.filtered.subs.vcf.site"
    shell:
        """awk -v OFS="\t" '/^[^#]/ {{print $1,$2,$2}}' {input} > {output}"""

rule bam_readcountm:
    input:
        bam="{sample}_STAR/Aligned.out.sorted.dupmkd.addGP.bam",
        site="{site_sample}.mutect2.filtered.subs.vcf.site"
    output:
        "{sample}.site_{site_sample}.bam-readcount.txt"
    log:
        "log/{sample}.site_{site_sample}.bam_readcount.log.txt"
    params:
        ref=config['reference'],
        min_mapping_quality=config['min_mapping_quality'],
        min_base_quality=config['min_base_quality']
    threads: 64
    shell:
        """
        bam-readcountm.sh -q {params.min_mapping_quality} -b {params.min_base_quality} -f {params.ref} -i {input.bam} -l {input.site} -o {output} -m `expr {threads} / 2` >{log} 2>&1
        """

rule readcount_annotate:
    input:
        vcf="{sample}.mutect2.filtered.subs{annotation}.vcf",
        txt="{sample_annotate}.site_{sample}.bam-readcount.txt"
    output:
        "{sample}.mutect2.filtered.subs{annotation,.*}.annotated_{sample_annotate}.vcf"
    log:
        "log/{sample}{annotation}.annotate_{sample_annotate}.log.txt"
    shell:
        """
        source activate bam-readcount
        vcf-readcount-annotator {input.vcf} {input.txt} DNA -t snv -o {output} -s {wildcards.sample_annotate} > {log} 2>&1
        conda deactivate
        """

rule vcf2bed:
    input:
        "{sample}.mutect2.filtered.subs.vcf"
    output:
        "{sample}.mutect2.filtered.subs.vcf.bed"
    shell:
        """awk -v OFS="\t" '/^[^#]/ {{ss=$2 - 1;print $1,ss,$2}}' {input} > {output}"""

rule gtf2bed:
    input:
        "{annotation}.gtf"
    output:
        "{annotation}.gtf.bed"
    params:
        field=config["gtf2bed_field"]
    shell:
        """
        awk -v FS="\t" -v OFS="\t" '/^[^#]/{{if ($3=="{params.field}") {{m=$4-1;match($9,/gene_id *"([^;]*)"[;$]/,n);match($9,/gene_name *"([^;]*)"[;$]/,o);match($9,/transcript_id *"([^;]*)"[;$]/,p);match($9,/transcript_name *"([^;]*)"[;$]/,q);match($9,/gene_type *"([^;]*)"[;$]/,r);match($9,/transcript_type *"([^;]*)"[;$]/,s);print $1,m,$5,n[1]";"o[1]";"p[1]";"q[1]";"r[1]";"s[1],100,$7}}}}' {input} > {output}
        """

rule bedtool_intersect:
    input:
        vcf_bed="{sample}.mutect2.filtered.subs.vcf.bed",
        annotation_bed=config["annotation_gtf"]+".bed"
    output:
        "{sample}.mutect2.filtered.subs.vcf.bedtools.txt"
    log:
        "log/{sample}.bedtool_intersect.log.txt"
    shell:
        "bedtools intersect -loj -a {input.vcf_bed} -b {input.annotation_bed} > {output} 2>{log}"

rule extract_site_info:
    input:
        "{sample}.mutect2.filtered.subs.vcf.bedtools.txt"
    output:
        "{sample}.annotation"
    shell:
        """awk -v OFS="\\t" '{{print $1,$3,$5,$6,$7,$9}}' {input} | awk -F'[\\t;]' -v OFS="\\t" '{{if ($3 != -1) {{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}}}}' > {output}"""

rule extract_vcf_info:
    input:
        lambda wildcards: f"{wildcards.sample}.mutect2.filtered.subs.annotated_{wildcards.sample}"+''.join(['.annotated_'+i for i in get_control_s(f"{wildcards.sample}")])+".vcf"
    output:
        "{sample}.res"
    params:
        target=lambda wildcards: f"{wildcards.sample} "+' '.join(get_control_s(f"{wildcards.sample}"))
    shell:
        """
        awk -v OFS="\\t" '/^[^#]/ {{print $1,$2,$4,$5,$7}}' {input} > {wildcards.sample}.res1
        ori=`grep "^#CHROM" {input}|head -1|awk '{{for (i=10;i<=NF;i++) {{printf "%s ",$i}};printf "\\n"}}'`
        target="{params.target}"
        matches=`python -c 'import sys;a=sys.argv[1].split(" ");b=sys.argv[2].split(" ");order=list(map(lambda x:[i+1 for i in range(len(a)) if a[i]==x][0],b));print("".join(["match($"+str(order[i]+9)+",/^[^:]*:([^:,]*),([^:,]*):.*/,m"+str(i)+");" for i in range(len(order))]))' "$ori" "$target"`
        inds=`python -c 'import sys;l=len(sys.argv[1].split(" "));print(",".join(["m"+str(i)+"[1],m"+str(i)+"[2]" for i in range(l)]))' "$target"`
        bash -c "awk -v OFS=\\\"\\t\\\" '/^[^#]/ {{"${{matches}}"print "${{inds}}"}}' "{input}" > "{wildcards.sample}".res2"
        paste {wildcards.sample}.res1 {wildcards.sample}.res2 > {output}
        rm {wildcards.sample}.res1 {wildcards.sample}.res2
        """

rule draw_all_test_sample:
    input:
        res=expand("{sample}.res",sample=list(test_controls['sample_id'])),
        annotations=expand("{sample}.annotation",sample=list(test_controls['sample_id']))
    output:
        "success"
    log:
        log="log/draw_all_test_sample.stdout.txt"
    params:
        ref=config["reference"],
        samples=' '.join(list(test_controls['sample_id'])),
        maximum_alt_percent_in_control=config['maximum_alt_percent_in_control'],
        minimum_editing_site_coverage_in_test=config['minimum_editing_site_coverage_in_test'],
        sequencelogo_span_len=config['sequencelogo_span_len']
    conda:
        "envs/draw_all_test_sample.yaml"
    script:
        "scripts/draw_all_test_samples.R"

