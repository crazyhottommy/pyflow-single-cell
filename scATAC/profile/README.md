# submit MAESTRO snakemake workflow to HPC 

read https://snakemake.readthedocs.io/en/v5.1.4/executable.html#profiles and 
https://github.com/Snakemake-Profiles/slurm for more details


```bash
snakemake --profile absolute/path/to/profiles/slurm

```

You may want to change the memory or time request in `slurm/cluster_config.yaml` file.