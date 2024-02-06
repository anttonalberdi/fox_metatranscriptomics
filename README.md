# Fox metatranscriptomics
Repository for the fox metatranscriptomic analysis

## Prepare working environment

```sh
# Load dependencies mamba and snakemake (not needed if they are already installed in root)
module load mamba/1.5.6 snakemake/7.20.0

# Clone this repository to the base directory
cd greenland_shotgun
git clone https://github.com/anttonalberdi/fox_metatranscriptomics.git
cd fox_metatranscriptomics

# Create screen session 
screen -S fox_metatranscriptomics

```

## Prepare reads
Store reads in the 'resources/reference/host' directory.

## Prepare host genome
Fox genome: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_003160815.1/
Store the files in the 'resources/reference/host' directory.

```sh
cd resources/reference/host
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/160/815/GCF_003160815.1_VulVul2.2/GCF_003160815.1_VulVul2.2_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/160/815/GCF_003160815.1_VulVul2.2/GCF_003160815.1_VulVul2.2_genomic.gtf.gz
cd ../../../
```

## Prepare microbial genomes
Store the microbial genomes in a single file called: 'resources/reference/microbiome/MAGs.fa'.

## Run snakemake
```sh
# Reopen screen session if needed
screen -r fox_metatranscriptomics

#Run snakemake
snakemake \
    --use-conda \
    --conda-frontend mamba \
    --rerun-incomplete \
    --jobs 50 \
    --cluster 'sbatch -o logs/{params.jobname}-slurm-%j.out --mem {resources.mem_gb}G --time {resources.time} -c {threads} --job-name={params.jobname} -v'
    --notemp \
    --slurm \
    --latency-wait 600
```
