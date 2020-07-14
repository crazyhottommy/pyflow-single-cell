
## when one processes multiple samples from the same experiment, he/she 
## will want to merge all the peaks across samples together and then get a new count matrix
## with the new peak set.

_countpeak_threads = 4


rule scatac_merge_peaks_batch:
	input: expand("Result/Analysis/{sample}/{sample}_final_peaks.bed", sample = ALL_SAMPLES)
	output: "Result/Analysis/Batch/all_samples_peaks.bed"
	params: catpeaksort = "Result/Analysis/Batch/all_samples_peaks_cat.bed"
	log: "Result/Log/merge_peaks_batch.log"
	shell:
		"""	
		cat {input}  \
		| sort -k1,1 -k2,2n | cut -f 1-4 > {params.catpeaksort}

		mergeBed -i {params.catpeaksort} | grep -v '_' | grep -v 'chrEBV' > {output}

		rm {params.catpeaksort}

		"""

rule scatac_countpeak_batch:
	input:
		finalpeak = "Result/Analysis/Batch/all_samples_peaks.bed",
		validbarcode = "Result/QC/{sample}/{sample}_scATAC_validcells.txt",
		frag = "Result/minimap2/{sample}/fragments_corrected_count.tsv"
	output:
		count = "Result/Analysis/Batch/{sample}/{sample}_peak_count.h5"
	params:
		species = config["species"],
		outdir = "Result/Analysis/Batch/{sample}",
		outpre = "{sample}"
	threads:
		_countpeak_threads
	benchmark:
		"Result/Benchmark/{sample}_PeakCount_batch.benchmark" 
	shell:
		"""
		MAESTRO scatac-peakcount --peak {input.finalpeak} --fragment {input.frag} --barcode {input.validbarcode} \
		--species {params.species} --cores {threads} --directory {params.outdir} --outprefix {params.outpre}
		"""


