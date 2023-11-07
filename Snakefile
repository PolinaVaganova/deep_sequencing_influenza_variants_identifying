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
		"{refernce}.fasta"
	output:
		"{refernce}.fasta.amb",
		"{refernce}.fasta.ann",
		"{refernce}.fasta.bwt",
		"{refernce}.fasta.pac",
		"{refernce}.fasta.sa"
	shell:
		"bwa index {input}"

rule bwa_alignment:
	input:
		"{refernce}.fasta.amb",
		"{refernce}.fasta.ann",
		"{refernce}.fasta.bwt",
		"{refernce}.fasta.pac",
		"{refernce}.fasta.sa",
		ref="{refernce}.fasta",
		reads="{sample}.fastq" 
	log:
		"logs/bwa.{refernce}.{sample}.log"
	output:
		temporary("{refernce}.{sample}.unsorted.bam")
	shell:
		"bwa mem {input.ref} {input.reads} 2>{log} | samtools view -S -b > {output}"

rule bam_sort:
	input:
		"{refernce}.{sample}.unsorted.bam"
	output:
		protected("{refernce}.{sample}.sorted.bam")
	threads: 8
	shell:
		"samtools sort --threads {threads} {input} > {output}"

rule mpileup:
	input:
		ref="{refernce}.fasta",
		align="{refernce}.{sample}.sorted.bam"
	output:
		"{refernce}.{sample}.mpileup"
	shell:
		"samtools mpileup -f {input.ref} {input.align} > {output} -d=0"

rule varscan_roommate_common:
	input:
		"{refernce}.{roommate}.mpileup"
	output:
		"{refernce}.{roommate}.common.vcf"
	shell:
		"varscan mpileup2snp {input} --min-var-freq 0.95 --variants --output-vcf 1 > {output}"

rule varscan_sample_rare:
	input:
		"{refernce}.{sample}.mpileup"
	output:
		"{refernce}.{sample}.rare.vcf"
	shell:
		"varscan mpileup2snp {input} --min-var-freq 0.001 --variants --output-vcf 1 > {output}"

rule collect_variants_characteristics:
	input:
		"{refernce}.{sample}.rare.vcf"
	output:
		"{refernce}.{sample}.csv"
	shell:
		"awk 'NR>24 {{split($10, arr, \":\"); print $4, $2, $5, arr[7]}}' {input} > {output}"

rule filter_variants_based_on_statistics:
	input:
		"{refernce}.{sample}.variants.csv"
	output:
		"controls_statistics.csv"
		"filtered_roommate_variants.csv"
	shell:
		"python3 get_true_mutations.py"
