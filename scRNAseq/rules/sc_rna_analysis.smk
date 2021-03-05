SCRIPT_PATH = config["SCRIPT_PATH"]
RSCRIPT_PATH = config["RSCRIPT_PATH"]

rule scrna_analysis:
    input:
        expression = "Result/QC/{sample}/{sample}_filtered_gene_count.h5"
    output:
        specificgene = "Result/Analysis/{sample}/{sample}_DiffGenes.tsv",
        clusterplot = "Result/Analysis/{sample}/{sample}_cluster.png",
        annotateplot = "Result/Analysis/{sample}/{sample}_annotated.png",
        tflist = "Result/Analysis/{sample}/{sample}.PredictedTFTop10.txt"    
    params:
        expression = "../QC/{sample}/{sample}_filtered_gene_count.h5",
        species = config["species"],
        outpre = "{sample}/{sample}",
        outdir = "Result/Analysis",
        lisadir = config["lisadir"],
        signature = config["signature"]
    benchmark:
        "Result/Benchmark/{sample}/{sample}_Analysis.benchmark"    
    threads:
        config["cores"]
    shell:
	    """
        python {SCRIPT_PATH}/lisa_path.py --species {params.species} --input {params.lisadir};
        Rscript {RSCRIPT_PATH}/scRNAseq_pipe.R --expression {params.expression} --species {params.species} \
        --prefix {params.outpre} --signature {params.signature} \
        --outdir {params.outdir} --thread {threads}
        """