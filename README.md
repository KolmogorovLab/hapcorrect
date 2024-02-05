# hapcorrect

###### Note: This repository is under extensive updates.

hapcorrect takes long-read alignment and phased heterozygous variants as input, and corrects the haplotypes phase-switching errors in BAMs around phased blocks as well as inside phase-blocks.

## Usage
```
python main.py --threads <4> --reference <ref.fa>  --target-bam <data.haplotagged.bam>  --out-dir-plots <coverage_plots>  --phased-vcf <data.phased.vcf.gz>  --genome-name <cellline/dataset name> --cut-threshold <150> 
```