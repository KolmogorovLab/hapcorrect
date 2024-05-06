# hapcorrect

###### Note: This repository is under extensive updates.

hapcorrect takes long-read alignment and phased heterozygous variants as input, and corrects the haplotypes phase-switching errors in BAMs around phased blocks as well as inside phase-blocks.


#### Phasing errors:
<img width="1373" alt="plots_example" src="src/images/1.png">

#### Phasing errors correction:
<img width="1373" alt="plots_example" src="src/images/2.png">

## Usage 

### Tumor-Normal Mode (requires normal phased VCF)
```
python main.py --threads <4> --reference <ref.fa>  --target-bam <data.tumor.bam>  --tumor-vcf <data.tumor.vcf.gz>  --normal-phased-vcf <data.normal_phased.vcf.gz>  --genome-name <cellline/dataset name> --cut-threshold <150> --out-dir-plots <genome_abc_output>
```
## Optional parameters
* `--rephase-normal-vcf` Set true to rephase normal VCF with hapcorrect phasing correction (default: disabled)
* `--rehaplotag-tumor-bam` Set true to rehaplotag the tumor BAM with new rephased VCF from above step (default: disabled)


### Tumor-only (requires tumor phased/haplotagged BAM and phased VCF)
```
python main.py --threads <4> --reference <ref.fa>  --target-bam <data.tumor_haplotagged.bam>  --tumor-vcf <data.tumor_phased.vcf.gz>  --genome-name <cellline/dataset name> --cut-threshold <150> --out-dir-plots <genome_abc_output>
```
## Optional parameters
* `--rephase-tumor-vcf` Set true to rephase tumor VCF with hapcorrect phasing correction (default: disabled)

## Note:
In some cases, phaseblocks are not good enough for `hapcorrect` to correct phase-switch errors (too small phaseblocks), in that scenrio, user can use `enable-simple-heuristics` parameter to apply simple heuristics which will assign higher coverage bins values to HP-1 and lower bins coverage values to HP-2.  