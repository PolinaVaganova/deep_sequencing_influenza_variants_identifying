URL_roommate = "http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/001/SRR1705851/SRR1705851.fastq.gz"
URL_control_1 = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/008/SRR1705858/SRR1705858.fastq.gz"
URL_control_2 = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/009/SRR1705859/SRR1705859.fastq.gz"
URL_control_3 = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/000/SRR1705860/SRR1705860.fastq.gz"

rule roommate_data_download:
	output:
		"roommate.fastq.gz"
	shell:
		"wget -O {output} {URL_roommate}"

rule control_1_data_download:
	output:
		"control_1.fastq.gz"
	shell:
		"wget -O {output} {URL_control_1}"

rule control_2_data_download:
	output:
		"control_2.fastq.gz"
	shell:
		"wget -O {output} {URL_control_2}"

rule control_3_data_download:
	output:
		"control_3.fastq.gz"
	shell:
		"wget -O {output} {URL_control_3}"

rule data_unzip:
	input:
		"{sample}.fastq.gz"
	output:
		"{sample}.fastq"
	shell:
		"gunzip -c {input} > {output}"

rule bwa_index:
	input:
		"{reference}.fasta"
	output:
		"{reference}.fasta.amb",
		"{reference}.fasta.ann",
		"{reference}.fasta.bwt",
		"{reference}.fasta.pac",
		"{reference}.fasta.sa"
	shell:
		"bwa index {input}"

rule bwa_alignment:
	input:
		"{reference}.fasta.amb",
		"{reference}.fasta.ann",
		"{reference}.fasta.bwt",
		"{reference}.fasta.pac",
		"{reference}.fasta.sa",
		ref="{reference}.fasta",
		reads="{sample}.fastq" 
	log:
		"logs/bwa.{reference}.{sample}.log"
	output:
		temporary("{reference}.{sample}.unsorted.bam")
	shell:
		"bwa mem {input.ref} {input.reads} 2>{log} | samtools view -S -b > {output}"

rule bam_sort:
	input:
		"{reference}.{sample}.unsorted.bam"
	output:
		protected("{reference}.{sample}.sorted.bam")
	threads: 8
	shell:
		"samtools sort --threads {threads} {input} > {output}"

rule mpileup:
	input:
		ref="{reference}.fasta",
		align="{reference}.{sample}.sorted.bam"
	output:
		"{reference}.{sample}.mpileup"
	shell:
		"samtools mpileup -f {input.ref} {input.align} > {output} -d=0"

rule varscan_sample_rare:
	input:
		"{reference}.{sample}.mpileup"
	output:
		"{reference}.{sample}.rare.vcf"
	shell:
		"varscan mpileup2snp {input} --min-var-freq 0.001 --variants --output-vcf 1 > {output}"

rule collect_variants_characteristics:
	input:
		"{reference}.{sample}.rare.vcf"
	output:
		"{reference}.{sample}.variants.csv"
	shell:
		"awk 'NR>24 {{split($10, arr, \":\"); print $4, $2, $5, arr[7]}}' {input} > {output}"


rule filter_variants_based_on_statistics:
	input:
		"{reference}.{sample}.variants.csv"
	output:
		"filtered.{reference}.{sample}.variants.csv"
	shell:
		"python filter_variants.py"


