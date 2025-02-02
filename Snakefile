with open(config['SAMPLES']) as fp:
    samples = fp.read().splitlines()


rule all:
         input:
            #Prepare samples
            #================
            #expand("galore/{sample}_R1_001_val_1.fq.gz", sample = samples),
            #expand("galore/{sample}_R2_001_val_2.fq.gz", sample = samples),
            #expand("{sample}.sam", sample = samples),
            #expand("{sample}.bam", sample = samples), 
            #expand("{sample}.sorted.bam", sample =samples),
            #expand("{sample}.sorted.rmDup.bam", sample =samples),
            #expand("{sample}.bigwig", sample = samples),
            "macs/ZF1_peaks.narrowPeak",
            "macs/ZF9_peaks.narrowPeak",
            "macs/ZF10_peaks.narrowPeak",
            "macs/ZF8_peaks.narrowPeak",
            "macs/ZF11_peaks.narrowPeak",
            "macs/ZF12_peaks.narrowPeak",
            "macs/ZF5_peaks.narrowPeak",
            "macs/ZF3_peaks.narrowPeak",
            "Motif_ZF1/seq.autonorm.tsv",
            "Motif_ZF9/seq.autonorm.tsv",
            "Motif_ZF10/seq.autonorm.tsv",
            "Motif_ZF8/seq.autonorm.tsv",
            "Motif_ZF11/seq.autonorm.tsv",
            "Motif_ZF12/seq.autonorm.tsv",
            "Motif_ZF5/seq.autonorm.tsv",
            "Motif_ZF3/seq.autonorm.tsv",
            "shared_peaks.narrowPeak"
            #"ZFs_KEGGpathways.pdf"  
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
         expand("{sample}.sorted.rmDup.bam", sample =samples)
      output: 
          "macs/ZF1_peaks.narrowPeak",
          "macs/ZF9_peaks.narrowPeak",
          "macs/ZF10_peaks.narrowPeak",
          "macs/ZF8_peaks.narrowPeak",
          "macs/ZF11_peaks.narrowPeak",
          "macs/ZF12_peaks.narrowPeak",
          "macs/ZF5_peaks.narrowPeak",
          "macs/ZF3_peaks.narrowPeak",
      shell: 
           """
	   macs2 callpeak -t 12701-ZF-1_S1.sorted.rmDup.bam -c 12701-ZF-7_S7.sorted.rmDup.bam -f BAMPE -g mm --outdir macs -n ZF1
	   macs2 callpeak -t 12701-ZF-9_S9.sorted.rmDup.bam -c 12701-ZF-7_S7.sorted.rmDup.bam -f BAMPE -g mm --outdir macs -n ZF9
	   macs2 callpeak -t 12701-ZF-10_S10.sorted.rmDup.bam -c 12701-ZF-7_S7.sorted.rmDup.bam -f BAMPE -g mm --outdir macs -n ZF10
	   macs2 callpeak -t 12701-ZF-8_S8.sorted.rmDup.bam -c 12701-ZF-7_S7.sorted.rmDup.bam -f BAMPE -g mm --outdir macs -n ZF8
	   macs2 callpeak -t 12701-ZF-11_S11.sorted.rmDup.bam -c 12701-ZF-7_S7.sorted.rmDup.bam -f BAMPE -g mm --outdir macs -n ZF11
	   macs2 callpeak -t 12701-ZF-12_S12.sorted.rmDup.bam -c 12701-ZF-7_S7.sorted.rmDup.bam -f BAMPE -g mm --outdir macs -n ZF12
	   macs2 callpeak -t 12701-ZF-5_S5.sorted.rmDup.bam -c 12701-ZF-7_S7.sorted.rmDup.bam -f BAMPE -g mm --outdir macs -n ZF5 
	   macs2 callpeak -t 12701-ZF-3_S3.sorted.rmDup.bam -c 12701-ZF-7_S7.sorted.rmDup.bam -f BAMPE -g mm --outdir macs -n ZF3
           """


rule annotateNarrowPeaks: 
      input: 
          "macs/ZF1_peaks.narrowPeak",
          "macs/ZF9_peaks.narrowPeak",
          "macs/ZF10_peaks.narrowPeak",
          "macs/ZF8_peaks.narrowPeak",
          "macs/ZF11_peaks.narrowPeak",
          "macs/ZF12_peaks.narrowPeak",
          "macs/ZF5_peaks.narrowPeak",
          "macs/ZF3_peaks.narrowPeak",
      params: 
           genome= config['GENOME'], 
           gtf = config['GTF']  
      output:
          "ZF1.annotatedPeak",
          "ZF9.annotatedPeak",
          "ZF10.annotatedPeak",
          "ZF8.annotatedPeak",
          "ZF11.annotatedPeak",
          "ZF12.annotatedPeak",
          "ZF5.annotatedPeak",
          "ZF3.annotatedPeak"
      shell: 
          """
          annotatePeaks.pl {input[[0]} {params.genome} -gtf {params.gtf}  > {output[0]} 
          annotatePeaks.pl {input[[1]} {params.genome} -gtf {params.gtf}  > {output[2]} 
          annotatePeaks.pl {input[[0]} {params.genome} -gtf {params.gtf}  > {output[3]}
          annotatePeaks.pl {input[[1]} {params.genome} -gtf {params.gtf}  > {output[4]} 
          annotatePeaks.pl {input[[0]} {params.genome} -gtf {params.gtf}  > {output[5]}
          annotatePeaks.pl {input[[1]} {params.genome} -gtf {params.gtf}  > {output[6]} 
          annotatePeaks.pl {input[[0]} {params.genome} -gtf {params.gtf}  > {output[7]}
          """ 


rule findMotifs:
      input:
          "macs/ZF1_summits.bed",
          "macs/ZF9_summits.bed",
          "macs/ZF10_summits.bed",
          "macs/ZF8_summits.bed",
          "macs/ZF11_summits.bed",
          "macs/ZF12_summits.bed",
          "macs/ZF5_summits.bed",
          "macs/ZF3_summits.bed",

      params:
          genome = config['GENOME'],
          output_dir0 = "Motif_ZF1",
          output_dir1 = "Motif_ZF9",
          output_dir2 = "Motif_ZF10",
          output_dir3 = "Motif_ZF8",
          output_dir4 = "Motif_ZF11",
          output_dir5 = "Motif_ZF12",
          output_dir6 = "Motif_ZF5",
          output_dir7 = "Motif_ZF3"
      output:
          "Motif_ZF1/seq.autonorm.tsv",
          "Motif_ZF9/seq.autonorm.tsv",
          "Motif_ZF10/seq.autonorm.tsv",
          "Motif_ZF8/seq.autonorm.tsv",
          "Motif_ZF11/seq.autonorm.tsv",
          "Motif_ZF12/seq.autonorm.tsv",
          "Motif_ZF5/seq.autonorm.tsv",
          "Motif_ZF3/seq.autonorm.tsv", 
      shell:
         """
          findMotifsGenome.pl {input[0]} {params.genome} {params.output_dir0} -size 200 -mask
          findMotifsGenome.pl {input[1]} {params.genome} {params.output_dir1} -size 200 -mask
          findMotifsGenome.pl {input[2]} {params.genome} {params.output_dir2} -size 200 -mask
          findMotifsGenome.pl {input[3]} {params.genome} {params.output_dir3} -size 200 -mask
          findMotifsGenome.pl {input[4]} {params.genome} {params.output_dir4} -size 200 -mask
          findMotifsGenome.pl {input[5]} {params.genome} {params.output_dir5} -size 200 -mask
          findMotifsGenome.pl {input[6]} {params.genome} {params.output_dir6} -size 200 -mask
          findMotifsGenome.pl {input[7]} {params.genome} {params.output_dir7} -size 200 -mask
         """


rule sharedPeaks:
     input:
           "macs/ZF9_peaks.narrowPeak",
           "macs/ZF5_peaks.narrowPeak", 
           "macs/ZF3_peaks.narrowPeak",
           "macs/ZF10_peaks.narrowPeak",
     output:
           "shared_peaks.narrowPeak"
     shell:
       """
       bedtools intersect -a macs/ZF9_peaks.narrowPeak -b macs/ZF5_peaks.narrowPeak -u | \
       bedtools intersect -a stdin -b macs/ZF3_peaks.narrowPeak -u | \
       bedtools intersect -a stdin -b macs/ZF10_peaks.narrowPeak > shared_peaks.narrowPeak
       """

rule plotQC: 
      input:
          "{sample}.sorted.rmDup.bam",
      output: 
          "{sample}.bw"   
      shell: 
        """
        bamCoverage -b {input} -o {output} --binSize 20 --normalizeUsing BPM --smoothLength 60 --extendReads 150 --centerReads 
        computeMatrix reference-point --referencePoint TSS -b 1000 -a 1000 -R mm10.bed -S 12701-ZF-7_S7.bw 12701-ZF-1_S1.bw 12701-ZF-3_S3.bw 12701-ZF-5_S5.bw 12701-ZF-8_S8.bw 12701-ZF-9_S9.bw \
        12701-ZF-10_S10.bw 12701-ZF-11_S11.bw --skipZeros -o matrix.gz
        plotHeatmap -m matrix.gz -out ZF_heatmap.png --colorMap coolwarm --heatmapHeight 10 --heatmapWidth 7 --sortUsing mean --zMin -4 --zMax 4 --samplesLabel "IgG" "Nfi_1" "Nfi_2" "IgG"

        """
         

rule pathway: 
    input: 
       "shared_peaks.narrowPeak"
    output: 
       "shared_KEGGpathways.pdf" 
    shell: 
        "Rscript pathways.R "
