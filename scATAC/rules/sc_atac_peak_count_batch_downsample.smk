### for multiple scATACseq samples, some samples maybe much more deeply sequenced
### let's downsample every CB tag added bam file to a certain number, and call peaks
### using all bam files by macs2 to get a peak set. Then go back to the original bam 
### file to get the counts in the peak set.

_downsample_threads = 4

if config.get("downsample"):
    target_reads = config['target_reads']


rule scatac_downsample_bam_batch:
    input: 
        bam = "Result/minimap2/{sample}/{sample}.sortedByPos.rmdp.CBadded.bam",
        bulk_stat = "Result/QC/{sample}/flagstat.txt"
    output:
        bam = "Result/miniap2/{sample}/{sample}.sortedByPos.rmdp.CBadded.downsample.bam"
    threads:
        _downsample_threads
    message:
        "downsampling {input}"
    log: 
        "Result/Log/{sample}_downsample_batch.log"
    benchmark:
        "Result/Benchmark/{sample}_downsample_batch.benchmark" 
    run:
        import re 
        with open(input.bulk_stat, "r") as f:
            # fifth line contains the number of mapped reads
            line = f.readlines()[4]
            match_number = re.match(r'(\d.+) \+.+', line)
            total_reads = int(match_number.group(1))
        if total_reads > target_reads:
            down_rate = target_reads/total_reads
        else:
            down_rate = 1
        shell("sambamba view -f bam -t {threads} --subsampling-seed=3 -s {rate} {inbam} \
            | samtools sort -m 2G -@ {threads} -T {outbam}.tmp > {outbam} 2> {log}".format(rate = down_rate, inbam = input[0], threads = threads, outbam = output[0], log = log))


rule scatac_downsample_peak_call:
    input: 
        bams = expand("Result/miniap2/{sample}/{sample}.sortedByPos.rmdp.CBadded.downsample.bam", sample = ALL_SAMPLES)
    output:
        peak = "Result/Analysis/Batch/all_samples_peaks.narrowPeak"
    params:
        name = "all_samples",
        genome = macs2_genome
    log:
        "Result/Log/batch_downsample_macs2_allpeak.log" 
    benchmark:
        "Result/Benchmark/batch_downsample_AllPeakCall.benchmark" 
    shell:
        "macs2 callpeak -f BAMPE -g {params.genome} --outdir Result/Analysis/Batch -n {params.name} -B -q 0.05 --nomodel --extsize=50 --keep-dup all -t {input.bams}"


rule scatac_countpeak_batch:
    input:
        finalpeak = "Result/Analysis/Batch/all_samples_peaks.narrowPeak",
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



