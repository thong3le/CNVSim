# CNVSim

This is a simulation that can generate high-throughput sequencing (HTS) data (BAM file format) for copy number variations such as deletion, inversion, tandem duplication, inverted duplication, interspersed direct duplication. 

## Prerequisites

* [wgsim](https://github.com/lh3/wgsim) 
* [bwa](https://github.com/lh3/bwa) 
* [samtools](https://github.com/samtools/) 

## Usage

```
python simulate_cnvs.py -h


usage: simulate_cnvs.py [-h] [--out OUT] [--chr CHR] ref n c

positional arguments:
  ref         reference genome (fasta file), e.g. hg38.fa
  n           number of variations, e.g. 1200
  c           sequencing coverage, e.g. 10

optional arguments:
  -h, --help  show this help message and exit
  --out OUT   output directory, default is the current directory
  --chr CHR   name of the chromosome that will contain the CNVs, default =
              'chr1'        
```

## Example

```
python simulate_cnvs.py --chr 1 --out out human_g1k_v37_gatk.fasta 1200 10
```

