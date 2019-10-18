IDS, = glob_wildcards("raw_reads/{id}_R1.fastq.gz")

report: "report/workflow.rst"

rule all:
	input:
		"abricate_summary_out/abricate_summary.tsv",
		"multiqc_fastp_out/fastp.html",
		"multiqc_quast_out/quast.html",
		"summary/summary.csv",
		"plots/tree.svg",
		"poppunk_export.zip",
		"snp-dists_out/snpmat.tsv"

rule fastp:
	input:
		fw = "raw_reads/{sample}_R1.fastq.gz",
		rv = "raw_reads/{sample}_R2.fastq.gz"
	output:
		fw = "ATQT_reads/{sample}_ATQT_R1.fastq.gz",
		rv = "ATQT_reads/{sample}_ATQT_R2.fastq.gz",
		json = "fastp_out/{sample}_fastp.json",
		html = "fastp_out/{sample}_fastp.html"
	conda:
		"envs/fastp.yaml"
	params:
		general = "--disable_length_filtering",
		compression_level = "9"
	log:	
		"logs/fastp/fastp_{sample}.log"
	threads: 6
	shell:
		"""
		fastp -w {threads} -z {params.compression_level} -i {input.fw} -o {output.fw} -I {input.rv} -O {output.rv} {params.general} --html {output.html} --json {output.json} 2>&1>{log}
		"""

rule kraken2:
	input:
		fw = "ATQT_reads/{sample}_ATQT_R1.fastq.gz",
		rv = "ATQT_reads/{sample}_ATQT_R2.fastq.gz",
		db = "db/kraken_db_nohuman"
	output:
		report = "kraken_out/{sample}_kraken2_report.txt"
	conda:
		"envs/kraken2.yaml"
	params:
		general = "--output - --gzip-compressed --paired"
	log:
		"logs/kraken2/kraken2_{sample}.log"
	threads: 6
	shell:
		"""
		kraken2 --db {input.db} {params.general} --threads {threads} --report {output.report} {input.fw} {input.rv} 2>&1>{log}
		"""

rule skesa:
	input:
		fw = "ATQT_reads/{sample}_ATQT_R1.fastq.gz",
		rv = "ATQT_reads/{sample}_ATQT_R2.fastq.gz"
	output:
		assembly = "skesa_out/{sample}_AIGHD.fasta"
	conda:
		"envs/skesa.yaml"
	params:
		min_length = "500"
	log:
		"logs/skesa/skesa_{sample}.log"
	threads: 6
	shell:
		"""
		skesa --fastq {input.fw},{input.rv} --use_paired_ends --min_contig {params.min_length} 1> {output.assembly} 2>{log}
		"""

rule mlst:
	input:
		assembly = "skesa_out/{sample}_AIGHD.fasta"
	output:
		mlst = "mlst/{sample}.tsv"
	conda:
		"envs/mlst.yaml"
	log:
		"logs/mlst/mlst_{sample}.log"
	threads: 6
	shell:
		"""
		mlst {input.assembly} 1> {output.mlst} 2>{log}
		"""

rule coverage:
	input:
		fw = "ATQT_reads/{sample}_ATQT_R1.fastq.gz",
		rv = "ATQT_reads/{sample}_ATQT_R2.fastq.gz",
		assembly = "skesa_out/{sample}_AIGHD.fasta"
	output:
		"coverage_out/{sample}.txt"
	conda:
		"envs/coverage.yaml"
	params:
		minimap_x = "sr"
	log:
		"logs/coverage/{sample}.log"
	threads: 6
	shell:
		"""
		minimap2 -a -x {params.minimap_x} -t {threads} {input.assembly} {input.fw} {input.rv} | samtools sort -l 0 --threads {threads} | bedtools genomecov -d -ibam stdin | awk '{{t += $3}} END {{print t/NR}}' 1>{output} 2>{log}
		"""

rule ska_fasta:
	input:
		assembly = "skesa_out/{sample}_AIGHD.fasta"
	output:
		skf = "ska_fasta_out/{sample}.skf"
	conda:
		"envs/ska.yaml"
	params:
		name = "ska_fasta_out/{sample}"
	log:
		"logs/ska_fasta/ska_fasta_{sample}.log"
	shell:
		"""
		ska fasta {input.assembly} -o {params.name} 2>&1>{log}
		"""

rule ska_align:
	input:
		expand("ska_fasta_out/{sample}.skf", sample=IDS)
	output:
		"ska_align_out/alignment_variants.aln"
	params:
		percentage = "0",
		variants = "-v",
		name = "ska_align_out/alignment"
	log:
		"logs/ska_align/ska_align.log"
	shell:
		"""
		ska align {params.variants} -p {params.percentage} -o {params.name} {input} 2>&1>{log}
		"""

rule abricate:
	input:
		assembly = "skesa_out/{sample}_AIGHD.fasta"
	output:
		"abricate_out/{sample}.tsv"
	conda:
		"envs/abricate.yaml"
	params:
		minid = "95",
		mincov = "60",
		db = "ncbi"
	log:
		"logs/abricate_out/abricate_{sample}.log"
	shell:
		"""
		abricate --db {params.db} --mincov {params.mincov} --minid {params.minid} {input.assembly} 1> {output} 2>{log}
		"""

rule abricate_summary:
	input:
		expand("abricate_out/{sample}.tsv", sample=IDS)
	output:
		"abricate_summary_out/abricate_summary.tsv"
	conda:
		"envs/abricate.yaml"
	log:
		"logs/abricate_summary/abricate_summary.log"
	shell:
		"""
		abricate --summary {input} 1>{output} 2>{log}
		"""

rule quast:
	input:
		assembly = "skesa_out/{sample}_AIGHD.fasta"
	output:
		directory("quast_out/{sample}")
	conda:
		"envs/quast.yaml"
	log:
		"logs/quast/quast_{sample}.log"
	shell:
		"""
		quast -o {output} {input.assembly}
		"""

rule multiqc_quast:
	input:
		expand("quast_out/{sample}", sample=IDS)
	output:
		html = report("multiqc_quast_out/quast.html", caption="report/multiqc_quast.rst", category="Quality control"),
		data = directory("multiqc_quast_out/quast_data")
	conda:
		"envs/multiqc.yaml"
	log:
		"logs/multiqc_quast/multiqc_quast.log"
	shell:
		"""
		OUTPUT={output.html}
		DIR=${{OUTPUT%/*}}
		NAME=${{OUTPUT##*/}}
		INPUT=$(for i in {input}; do echo ${{i%/*}}; done | sort | uniq)
		multiqc -n ${{NAME}} -o ${{DIR}} ${{INPUT}} 2>&1>{log}
		"""

rule multiqc_fastp:
	input:
		expand("fastp_out/{sample}_fastp.json", sample=IDS)
	output:
		html = report("multiqc_fastp_out/fastp.html", caption="report/multiqc_fastp.rst", category="Quality control"),
		data = directory("multiqc_fastp_out/fastp_data")
	conda:
		"envs/multiqc.yaml"
	log:
		"logs/multiqc_quast/multiqc_fastp.log"
	shell:
		"""
		OUTPUT={output.html}
		DIR=${{OUTPUT%/*}}
		NAME=${{OUTPUT##*/}}
		INPUT=$(for i in {input}; do echo ${{i%/*}}; done | sort | uniq)
		multiqc -n ${{NAME}} -o ${{DIR}} ${{INPUT}} 2>&1>{log}
		"""

rule poppunk_K2:
	input:
		expand("skesa_out/{sample}_AIGHD.fasta", sample=IDS)
	output:
		K2 = directory("poppunk_out/K2")
	conda:
		"envs/poppunk.yaml"
	params:
		K = "2",
		mink = "15",
		kstep = "2",
		name = "K2",
		perplexity = "20",
		maxadist = "1"
	log:
		"../logs/poppunk_K2.log"
	threads: 6
	shell:
		"""
		cd poppunk_out
		bash ../scripts/prepare_poppunk.sh {input}
		ls poppunk_temp/* > list_poppunk.txt
		
		poppunk --create-db --overwrite --r-files list_poppunk.txt --output {params.name} --threads {threads} --plot-fit 5 --min-k {params.mink} --k-step {params.kstep} 2>&1>{log}
		poppunk --fit-model --distances {params.name}/{params.name}.dists --max-a-dist {params.maxadist} --ref-db {params.name} --output {params.name} --full-db --K {params.K} --microreact --cytoscape --phandango --grapetree --perplexity {params.perplexity} 2>&1>{log}
		
		rm -rf poppunk_temp list_poppunk.txt
		"""

rule poppunk_K3:
	input:
		expand("skesa_out/{sample}_AIGHD.fasta", sample=IDS)
	output:
		K3 = directory("poppunk_out/K3")
	conda:
		"envs/poppunk.yaml"
	params:
		K = "3",
		mink = "15",
		kstep = "2",
		name = "K3",
		perplexity = "20",
		maxadist = "1"
	log:
		"../logs/poppunk_K3.log"
	threads: 6
	shell:
		"""
		cd poppunk_out
		bash ../scripts/prepare_poppunk.sh {input}
		ls poppunk_temp/* > list_poppunk.txt
		
		poppunk --create-db --overwrite --r-files list_poppunk.txt --output {params.name} --threads {threads} --plot-fit 5 --min-k {params.mink} --k-step {params.kstep} 2>&1>{log}
		poppunk --fit-model --distances {params.name}/{params.name}.dists --max-a-dist {params.maxadist} --ref-db {params.name} --output {params.name} --full-db --K {params.K} --microreact --cytoscape --phandango --grapetree --perplexity {params.perplexity} 2>&1>{log}
		
		rm -rf poppunk_temp list_poppunk.txt
		"""

rule poppunk_K4:
	input:
		expand("skesa_out/{sample}_AIGHD.fasta", sample=IDS)
	output:
		K4 = directory("poppunk_out/K4")
	conda:
		"envs/poppunk.yaml"
	params:
		K = "4",
		mink = "15",
		kstep = "2",
		name = "K4",
		perplexity = "20",
		maxadist = "1"
	log:
		"../logs/poppunk_K4.log"
	threads: 6
	shell:
		"""
		cd poppunk_out
		bash ../scripts/prepare_poppunk.sh {input}
		ls poppunk_temp/* > list_poppunk.txt
		
		poppunk --create-db --overwrite --r-files list_poppunk.txt --output {params.name} --threads {threads} --plot-fit 5 --min-k {params.mink} --k-step {params.kstep} 2>&1>{log}
		poppunk --fit-model --distances {params.name}/{params.name}.dists --max-a-dist {params.maxadist} --ref-db {params.name} --output {params.name} --full-db --K {params.K} --microreact --cytoscape --phandango --grapetree --perplexity {params.perplexity} 2>&1>{log}
		
		rm -rf poppunk_temp list_poppunk.txt
		"""

rule iqtree:
	input:
		"ska_align_out/alignment_variants.aln"
	output:
		directory("iqtree_out")
	conda:
		"envs/iqtree.yaml"
	params:
		model = "GTR+G+ASC",
		prefix = "iqtree"
	log:
		"logs/iqtree/iqtree.log"
	shell:
		"""
		mkdir -p {output} && cd {output}
		iqtree -fast -s ../{input} -m {params.model} -nt AUTO -pre {params.prefix} 
		"""

rule summary:
	input:
		fastp = expand("fastp_out/{sample}_fastp.json", sample=IDS),
		kraken = expand("kraken_out/{sample}_kraken2_report.txt", sample=IDS),
		quast = expand("quast_out/{sample}", sample=IDS),
		abricate = expand("abricate_out/{sample}.tsv", sample=IDS),
		mlst = expand("mlst/{sample}.tsv", sample=IDS),
		coverage = expand("coverage_out/{sample}.txt", sample=IDS)
	output:
		report("summary/summary.csv", caption="report/summary.rst", category="Summary")
	shell:
		"""
		bash scripts/summary_isolates.sh {input.quast} > {output}
		"""

rule plottree:
	input:
		"iqtree_out"
	output:
		report("plots/tree.svg", caption="report/plottree.rst", category="Plots")
	conda:
		"envs/plottree.yaml"
	log:
		"logs/plottree.log"
	shell:
		"""
		Rscript scripts/plottree.R {input}/iqtree.treefile {output}
		"""

rule zip:
	input:
		K2 = "poppunk_out/K2",
		K3 = "poppunk_out/K3",
		K4 = "poppunk_out/K4" 
	output:
		report("poppunk_export.zip", caption="report/poppunk_zip.rst", category="Poppunk export")
	log:
		"logs/poppunk_export_zip.log"
	shell:
		"""
		bash scripts/zip_poppunk.sh {input.K2} {input.K3} {input.K4} 1>&2>{log}
		"""

rule snpdists:
	input:
		"ska_align_out/alignment_variants.aln"
	output:
		"snp-dists_out/snpmat.tsv"
	conda:
		"envs/snpdists.yaml"
	log:
		"logs/snpdists.log"
	shell:
		"""
		snp-dists {input} 1> {output} 2>{log}
		"""
