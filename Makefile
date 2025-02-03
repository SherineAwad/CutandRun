genome=  "/nfs/turbo/umms-thahoang/sherine/REFERENCES/Mus_musculus/Ensembl/GRCm38/Sequence/WholeGenomeFasta/genome.fa"
gtf= "/nfs/turbo/umms-thahoang/sherine/REFERENCES/Mus_musculus/Ensembl/GRCm38/Annotation/Genes/genes.gtf" 

SRR15632729: 
	#GSM5542341
	wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR15632726/SRR15632726
	#GSM5542342
	wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR15632727/SRR15632727
	#GSM5542343
	wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR15632728/SRR15632728
	#GSM5542344
	wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR15632729/SRR15632729


SRR15632728_2.fastq:
	fastq-dump --split-files SRR15632727
	fastq-dump --split-files SRR15632729
	fastq-dump --split-files SRR15632726
	fastq-dump --split-files SRR15632728


# Makefile rule to rename files
rename_files:
	for file in 12701-ZF-*_*_R2_001.fastq.gz; do \
        new_name=$$(echo $$file | sed -E 's/^12701-ZF-([0-9]+)_S[0-9]+_R2_001.fastq.gz/ZF\1_R2_001.fastq.gz/'); \
        mv $$file $$new_name; \
	done
	# Makefile rule to rename files
	for file in 12701-ZF-*_*_R1_001.fastq.gz; do \
        new_name=$$(echo $$file | sed -E 's/^12701-ZF-([0-9]+)_S[0-9]+_R1_001.fastq.gz/ZF\1_R1_001.fastq.gz/'); \
        mv $$file $$new_name; \
    done


