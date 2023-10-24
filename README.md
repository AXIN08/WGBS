# A simple and brutal WGBS workflow

## Start
### 1. Prepare reference genome using Bismark
You can choose hisat2 or bowtie2 for preparing, and for the next step, the alignment tool you used in Bismark should be the same.

`bismark_genome_preparation --hisat2/bowtie2 --verbose bisulfite_genome/`

### 2. Prepare config file
Prepare the config file for Snakemake, better to choose set the config path in the workflow/config/ directory. After generating the config file, you can change the parameters or paths as needed.

`python scripts/generate_config.py -f sample.txt -o output.path -C workflow/config/ -r bisulfite_genome/`

### 3. Run Snakemake

`snakemake -s workflow/Snakefile --cores [Int]`

### Noted
The file "sample.txt" should have "id  fq1  fq2" as the first line

