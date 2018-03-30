# CNVSim

This is a simulation that can generate high-throughput sequencing (HTS) data (BAM file format) for copy number variations such as deletion, inversion, tandem duplication, inverted duplication, interspersed direct duplication. 

## Prerequisites

* [wgsim](https://github.com/lh3/wgsim) simulation 
* [bwa](https://github.com/lh3/bwa) 
* [samtools](https://github.com/samtools/) 

## Usage

```
python simulate_cnvs.py -h
```

## Example

```
python simulate_cnvs.py --chr 1 --out out human_g1k_v37_gatk.fasta 1200 10
```

