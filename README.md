# hapcorrect

###### Note: This repository is under extensive updates.

hapcorrect takes long-read alignment and phased heterozygous variants as input, and corrects the haplotypes phase-switching errors in BAMs around phased blocks as well as inside phase-blocks.

## Usage 

## For phase-switch errors correction
Input tumor BAM should be haplotagged and use `--phaseblock-flipping-enable True` enabled (default: disabled).

### Tumor-Normal Mode (requires normal phased VCF)
```
python main.py --threads <4> --reference <ref.fa>  --target-bam <data.tumor.bam>  --tumor-vcf <data.tumor.vcf.gz>  --normal-phased-vcf <data.normal_phased.vcf.gz>  --unphased-reads-coverage-enable True --phaseblock-flipping-enable True  --genome-name <cellline/dataset name> --cut-threshold <150> --out-dir-plots <genome_abc_output>
```

### Tumor-only (requires tumor phased/haplotagged BAM and phased VCF)
```
python main.py --threads <4> --reference <ref.fa>  --target-bam <data.tumor_haplotagged.bam>  --tumor-vcf <data.tumor_phased.vcf.gz>  --unphased-reads-coverage-enable True --phaseblock-flipping-enable True  --genome-name <cellline/dataset name> --cut-threshold <150> --out-dir-plots <genome_abc_output>
```
