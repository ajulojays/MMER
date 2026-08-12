.PHONY: metadata test download fastqc dada2 resistance all
metadata:
	python scripts/00_build_metadata.py

test: metadata
	python tests/test_metadata.py

download: test
	bash scripts/01_download_fastq.sh

fastqc:
	mkdir -p results/qc/fastqc
	fastqc -t 8 -o results/qc/fastqc data/raw/*.fastq.gz
	multiqc -o results/qc/multiqc results/qc/fastqc

dada2:
	Rscript scripts/02_dada2_asv.R

resistance:
	Rscript scripts/03_resistance_metrics.R
	Rscript scripts/04_within_cow_contrasts.R

all: metadata test download fastqc dada2 resistance
