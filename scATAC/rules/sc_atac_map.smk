_minmap_threads = 8
_picard_threads = 6

rule scatac_map:
	input:
		fasta = config["genome"]["fasta"],
		r1 = "Result/fastq/{sample}/{sample}_R1.fastq",
		r3 = "Result/fastq/{sample}/{sample}_R3.fastq"
	output:
			bam = temp("Result/minimap2/{sample}/{sample}.sortedByPos.bam")
	threads:
		_minmap_threads
	benchmark:
			"Result/Benchmark/{sample}_Minimap2.benchmark"
	shell:
		"""
		minimap2 -ax sr -t {threads} {input.fasta} {input.r1} {input.r3} \
		| samtools view --threads {threads} -b \
		| samtools sort --threads {threads} -o {output.bam}
		"""


rule scatac_rmdp:
    input:
        bam = "Result/minimap2/{sample}/{sample}.sortedByPos.bam",
    output:
        bam = "Result/minimap2/{sample}/{sample}.sortedByPos.rmdp.bam",
        metric = "Result/minimap2/{sample}/{sample}.rmdp.txt",
        fragbed = "Result/QC/{sample}/{sample}_frag.bed",
        #tmp = temp(directory("Result/{sample}"))
    params:
        sam = "Result/minimap2/{sample}/{sample}.sortedByPos.rmdp.sample.sam"
    threads:
        _picard_threads
    benchmark:
        "Result/Benchmark/{sample}_Rmdp.benchmark" 
    shell:
    	"""
        picard MarkDuplicates INPUT={input.bam} OUTPUT={output.bam} METRICS_FILE={output.metric} TMP_DIR=Results/tmp

        samtools view -@ {threads} -s 0.01 -o {params.sam} {input.bam}

        awk '{{if ($9>0) print $9}}' {params.sam} > {output.fragbed}
        """


