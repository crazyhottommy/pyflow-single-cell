_shortfragment_threads = 8

macs2_genome = "hs" if config["species"] == "GRCh38" else "mm"

rule scatac_allpeakcall:
    input:
        bam = "Result/minimap2/{sample}/{sample}.sortedByPos.rmdp.CBadded.bam" 
    output:
        peak = "Result/Analysis/{sample}/{sample}_all_peaks.narrowPeak",
        bdg = "Result/Analysis/{sample}/{sample}_all_treat_pileup.bdg"
    params:
        name = "{sample}_all",
        genome = macs2_genome
    log:
        "Result/Log/{sample}_macs2_allpeak.log" 
    benchmark:
        "Result/Benchmark/{sample}_AllPeakCall.benchmark" 
    shell:
        "macs2 callpeak -f BAMPE -g {params.genome} --outdir Result/Analysis/{wildcards.sample} -n {params.name} -B -q 0.05 --nomodel --extsize=50 --SPMR --keep-dup all -t {input.bam}"

if config["shortpeaks"]:
    rule scatac_shortfragment:
        input:
            bam = "Result/minimap2/{sample}/{sample}.sortedByPos.rmdp.CBadded.bam" 
        output:
            shortbam = "Result/minimap2/{sample}/{sample}.sortedByPos.rmdp.CBadded.150bp.bam" 
        threads:
            _shortfragment_threads
        benchmark:
            "Result/Benchmark/{sample}_ShortFrag.benchmark" 
        shell:
            "samtools view -@ {threads} -h {input.bam} | "
            "awk -F'\\t' 'function abs(x){{return ((x < 0.0) ? -x : x)}} {{if (abs($9)<=150) print}}' | "
            "samtools view -@ {threads} -b -o {output.shortbam}"
    
    rule scatac_shortpeakcall:
        input:
            shortbam = "Result/minimap2/{sample}/{sample}.sortedByPos.rmdp.CBadded.150bp.bam" 
        output:
            bed = "Result/Analysis/{sample}/{sample}_150bp_peaks.narrowPeak"
        params:
            name = "{sample}_150bp",
            genome = macs2_genome 
        log:
            "Result/Log/{sample}_macs2_shortpeak.log" 
        benchmark:
            "Result/Benchmark/{sample}_ShortPeakCall.benchmark" 
        shell:
            "macs2 callpeak -f BAMPE -g {params.genome} --outdir Result/Analysis/{wildcards.sample} -n {params.name} -B -q 0.05 --nomodel --extsize=50 --SPMR -t --keep-dup all {input.shortbam}"

if config["custompeaks"] and config["shortpeaks"]:
    rule scatac_mergepeak:
        input:
            allpeak = "Result/Analysis/{sample}/{sample}_all_peaks.narrowPeak" ,
            shortpeak = "Result/Analysis/{sample}/{sample}_150bp_peaks.narrowPeak",
            custompeak = config["custompeaksloc"]
        output:
            finalpeak = "Result/Analysis/{sample}/{sample}_final_peaks.bed" 
        params:
            catpeaksort = "Result/Analysis/{sample}/{sample}_cat_peaks.bed" 
        benchmark:
            "Result/Benchmark/{sample}_PeakMerge.benchmark" 
        shell:
            """
            cat {input.allpeak} {input.shortpeak} {input.custompeak} \
            | sort -k1,1 -k2,2n | cut -f 1-4 > {params.catpeaksort}

            mergeBed -i {params.catpeaksort} | grep -v '_' | grep -v 'chrEBV' > {output.finalpeak}

            rm {params.catpeaksort}
            """

elif config["custompeaks"]:
    rule scatac_mergepeak:
        input:
            allpeak = "Result/Analysis/{sample}/{sample}_all_peaks.narrowPeak",
            custompeaks = config["custompeaksloc"]
        output:
            finalpeak = "Result/Analysis/{sample}/{sample}_final_peaks.bed" 
        params:
            catpeaksort = "Result/Analysis/{sample}/{sample}_cat_peaks.bed"  
        benchmark:
            "Result/Benchmark/{sample}_PeakMerge.benchmark"          
        shell:
            """
            cat {input.allpeak} {input.custompeaks} \
            | sort -k1,1 -k2,2n | cut -f 1-4 > {params.catpeaksort}

            mergeBed -i {params.catpeaksort} | grep -v '_' | grep -v 'chrEBV' > {output.finalpeak}

            rm {params.catpeaksort}
            """

elif config["shortpeaks"]:
    rule scatac_mergepeak:
        input:
            allpeak = "Result/Analysis/{sample}/{sample}_all_peaks.narrowPeak",
            shortpeak = "Result/Analysis/{sample}/{sample}_150bp_peaks.narrowPeak" 
        output:
            finalpeak = "Result/Analysis/{sample}/{sample}_final_peaks.bed" 
        params:
            catpeaksort = "Result/Analysis/{sample}/{sample}_cat_peaks.bed" 
        benchmark:
            "Result/Benchmark/{sample}_PeakMerge.benchmark" 
        shell:
            """
            cat {input.allpeak} {input.shortpeak} \
            | sort -k1,1 -k2,2n | cut -f 1-4 > {params.catpeaksort}

            mergeBed -i {params.catpeaksort} | grep -v '_' | grep -v 'chrEBV' > {output.finalpeak}

            rm {params.catpeaksort}
            """

else:
    rule scatac_mergepeak:
        input:
            allpeak = "Result/Analysis/{sample}/{sample}_all_peaks.narrowPeak" 
        output:
            finalpeak = "Result/Analysis/{sample}/{sample}_final_peaks.bed" 
        params:
            catpeaksort = "Result/Analysis/{sample}/{sample}_cat_peaks.bed" 
        benchmark:
            "Result/Benchmark/{sample}_PeakMerge.benchmark" 
        shell:
            """
            cat {input.allpeak} \
            | sort -k1,1 -k2,2n | cut -f 1-4 > {params.catpeaksort}

            mergeBed -i {params.catpeaksort} | grep -v '_' | grep -v 'chrEBV' > {output.finalpeak}

            rm {params.catpeaksort}
            """





