_cat_threads= 2

rule scatac_preprocess:
	input:
		r1 = lambda wildcards: FILES[wildcards.sample]["R1"],
		r2 = lambda wildcards: FILES[wildcards.sample]["R2"],
		r3 = lambda wildcards: FILES[wildcards.sample]["R3"]
	output:
		r1cat = "Result/fastq/{sample}/{sample}_R1.fastq",
		r2cat = "Result/fastq/{sample}/{sample}_R2.fastq",
		r3cat = "Result/fastq/{sample}/{sample}_R3.fastq"
	log:
		"Result/Log/{sample}_preprocess.log"
	benchmark:
		"Result/Benchmark/{sample}_Preprocess.benchmark" 
	threads: _cat_threads
	shell:
		"""
		gunzip -c {input.r1} > {output.r1cat} 2> {log}
		gunzip -c {input.r2} > {output.r2cat} 2>> {log}
		gunzip -c {input.r1} > {output.r3cat} 2>> {log}
		
		"""
