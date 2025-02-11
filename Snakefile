with open(config['SAMPLES']) as fp:
    samples = fp.read().splitlines()
TREATMENT_SAMPLES = ["ZF2","ZF4", "ZF6", "ZF8"]
CONTROL_SAMPLE = "ZF7"

rule all:
         input:
            #Prepare samples
            #================
            expand("galore/{sample}_R1_001_val_1.fq.gz", sample = samples),
            expand("galore/{sample}_R2_001_val_2.fq.gz", sample = samples),
            expand("{sample}.sam", sample = samples),
            expand("{sample}.bam", sample = samples), 
            expand("{sample}.sorted.bam", sample =samples),
            expand("{sample}.sorted.rmDup.bam", sample =samples),
            expand("{sample}.bigwig", sample = samples),
            expand("macs/{sample}_peaks.narrowPeak", sample=TREATMENT_SAMPLES),
            expand( "macs/{sample}_summits.bed",  sample=TREATMENT_SAMPLES),
            expand("{sample}.annotatedPeak", sample=TREATMENT_SAMPLES),
            expand("Motif_{sample}/seq.autonorm.tsv", sample=TREATMENT_SAMPLES),
            expand("{sample}.bw", sample=TREATMENT_SAMPLES),
            expand("{sample}.bed",sample = CONTROL_SAMPLE),
            expand("{sample}_QC.pdf",sample = CONTROL_SAMPLE),
            expand("{sample}_QC.pdf",sample=TREATMENT_SAMPLES),
            "merged_peaks.narrowPeak",
            "merged_KEGGpathways.pdf" 
rule trim: 
       input:  
            r1 = "{sample}_R1_001.fastq.gz",
            r2 = "{sample}_R2_001.fastq.gz"
       output: 
          "galore/{sample}_R1_001_val_1.fq.gz",
          "galore/{sample}_R2_001_val_2.fq.gz"
       shell: 
           """
           mkdir -p galore 
           mkdir -p fastqc 
           trim_galore --cores 4 --gzip --retain_unpaired --fastqc --fastqc_args "--outdir fastqc" -o galore --paired {input.r1} {input.r2}
           """ 


rule align:
              input:               
                  "galore/{sample}_R1_001_val_1.fq.gz",
                  "galore/{sample}_R2_001_val_2.fq.gz"
              params:
                   index=config['INDEX'],
                   mem = config['MEMORY'],
                   cores = config['CORES']
              output:
                   "{sample}.sam",
                   "{sample}_hist.txt" 
              shell:
                   """
                   bowtie2 --local --very-sensitive-local --no-mixed --no-discordant --phred33 -I 10 -X 700 -p {params.cores} -x {params.index} -1 {input[0]} -2 {input[1]} -S {output[0]}  &> {output[1]}
                   """




rule samTobam:
             input: 
                 "{sample}.sam",
             output: 
                 "{sample}.bam"
             shell: 
                   """
                   samtools view -bS {input} > {output}
                   """


rule sort:
             input:
                  "{sample}.bam"
             output:
                  "{sample}.sorted.bam"
             shell:
                  """
                  picard SortSam I={input}  O={output} SORT_ORDER=coordinate
                  """


rule remove_duplicates:
       input: 
        "{sample}.sorted.bam"
       output: 
         "{sample}.sorted.rmDup.bam",
         "{sample}.rmDup.txt"
       shell: 
           """
            picard MarkDuplicates I={input} O={output[0]} REMOVE_DUPLICATES=true METRICS_FILE={output[1]} 
           """

rule index: 
      input: 
         "{sample}.sorted.rmDup.bam"
      output: 
         "{sample}.sorted.rmDup.bam.bai"
      shell: 
          """
          samtools index {input} 
          """ 
rule bamCoverage: 
       input: 
        "{sample}.sorted.rmDup.bam",
        "{sample}.sorted.rmDup.bam.bai" 
       output: 
        "{sample}.bigwig" 
       params: 
         genome_size = config['Genome_Size'], 
         binsize = config['BINSIZE'], 
         num_processors = config['Num_Processors'] 
       shell: 
          """ 
          bamCoverage -b {input[0]} -p {params.num_processors}  --normalizeUsing RPGC --effectiveGenomeSize {params.genome_size} --binSize {params.binsize} -o {output} 
          """ 


rule macs_bed:
    input:
        "{sample}.sorted.rmDup.bam",
        expand("{sample}.sorted.rmDup.bam", sample =CONTROL_SAMPLE) 
    params: 
        "{sample}"
    output:
        "macs/{sample}_peaks.narrowPeak",
         "macs/{sample}_summits.bed"
    shell:
        """
        macs2 callpeak -t {input[0]} -c {input[1]} -f BAMPE -g mm --outdir macs -n {params}
        """


rule annotateNarrowPeaks:
    input:
        "macs/{sample}_peaks.narrowPeak"
    params:
        genome=config['GENOME'],
        gtf=config['GTF']
    output:
        "{sample}.annotatedPeak"
    run:
        for sample in TREATMENT_SAMPLES:
            shell("annotatePeaks.pl macs/{sample}_peaks.narrowPeak {params.genome} -gtf {params.gtf} > {sample}.annotatedPeak")



rule findMotifs:
    input:
        "macs/{sample}_summits.bed"
    params:
        genome=config['GENOME']
    output:
        "Motif_{sample}/seq.autonorm.tsv"
    shell:
         """
          findMotifsGenome.pl {input} {params.genome} {output} -size 200 -mask
         """



rule prepareQC_bw:
    input:
       "{sample}.sorted.rmDup.bam", 
    output:
       "{sample}.bw"
    shell:
        """
            bamCoverage -b {input} -o {output}  --normalizeUsing BPM
        """

rule Control_QCheatmap: 
      input: 
         expand("{sample}.sorted.rmDup.bam",sample = CONTROL_SAMPLE), 
         expand("{sample}.bw",sample = CONTROL_SAMPLE)
      output: 
           expand("{sample}.bed",sample = CONTROL_SAMPLE),
           expand("{sample}_QC.pdf",sample = CONTROL_SAMPLE) 
      shell: 
          """
            #bedtools bamtobed -i {input} > {output[0]}
            Rscript hm.R {output[0]} {input[1]} 
          """          
rule QCheatmap: 
     input: 
        "macs/{sample}_summits.bed",
        "{sample}.bw" 
     output: 
        "{sample}_QC.pdf"
     shell: 
        """
         Rscript hm.R {input[0]} {input[1]} 
        """

rule merge: 
     input: 
         expand("macs/{sample}_peaks.narrowPeak", sample =TREATMENT_SAMPLES) 
     output: 
         "merged_peaks.narrowPeak"
     shell: 
        """
        cat {input} > {output} 
        """
 
rule pathway: 
    input: 
       "merged_peaks.narrowPeak"
    output: 
       "merged_KEGGpathways.pdf" 
    shell: 
        "Rscript pathways.R "
