# Makefile rule to rename files
rename_files:
	for file in 12701-ZF-*_*_R2_001.fastq.gz; do \
        new_name=$$(echo $$file | sed -E 's/^12701-ZF-([0-9]+)_S[0-9]+_R2_001.fastq.gz/ZF\1_R2_001.fastq.gz/'); \
        mv $$file $$new_name; \
	done
	for file in 12701-ZF-*_*_R1_001.fastq.gz; do \
        new_name=$$(echo $$file | sed -E 's/^12701-ZF-([0-9]+)_S[0-9]+_R1_001.fastq.gz/ZF\1_R1_001.fastq.gz/'); \
        mv $$file $$new_name; \
    done



