# fox_metatranscriptomics
Repository for the fox metatranscriptomic analysis

## Prepare working environment

```sh
# Load dependencies mamba and snakemake (not needed if they are already installed in root)
module load mamba/1.5.6 snakemake/7.20.0

# Clone this repository to the base directory
cd greenland_shotgun
git clone https://github.com/3d-omics/mg_assembly.git
cd fox_metatranscriptomics

# Create screen session 
screen -S fox_metatranscriptomics

```

## Get reference genome and annotation

Fox genome: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_003160815.1/
Store the files in the 'resources/reference/host' directory.

```sh
cd resources/reference/host
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/160/815/GCF_003160815.1_VulVul2.2/GCF_003160815.1_VulVul2.2_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/160/815/GCF_003160815.1_VulVul2.2/GCF_003160815.1_VulVul2.2_genomic.gtf.gz
cd ../../../
```

