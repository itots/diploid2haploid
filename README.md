# diploid2haploid.py
The script converts diploid genotypes to haploid genotypes for the sex
 chromosomes of males in a vcf file. The scripts modifies only GT fields
  with others remain. Heterozygous genotypes in the sex chromosomes of males
  , if any, are converted to be missing.  
    
  The script was tested in Python 3.7 with the vcf file produced by stacks 2.5.

## Dependencies
Python 3
```
argparse
```
## Inputs
- vcf file
- sex assignment file, which is tab-delimited text file with each line
 consisting of sample name and sex (female or male), like this: 
  ```
  lib1_1_A	female
  lib1_1_B	male
  lib1_1_C	female
  lib1_1_D	male
  ```
  It is required when a input vcf file contains not only males but also
   females.
   
## Outputs
The script outputs modified vcf file and summary report as follows. 
- output.vcf
- output_summary.txt

These are provided in a current directory by default.  
Output directory and file name can be changed via `-o` option.
## Usage
```
usage: diploid2haploid.py [-h] -i INPUT_VCF -y YCHROM -x XCHROM
                          [-o OUTPUT_VCF] [-s SEX]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_VCF, --input_vcf INPUT_VCF
                        Name of input vcf file
  -y YCHROM, --ychrom YCHROM
                        Name of Y chromosome in reference sequence
  -x XCHROM, --xchrom XCHROM
                        Name of X chromosome in reference sequence
  -o OUTPUT_VCF, --output_vcf OUTPUT_VCF
                        Name of output vcf file. By default, it is
                        ./output.vcf
  -s SEX, --sex SEX     File of sex assignment, required when input vcf file
                        contains females
```

## Examples
```
./diploid2haploid.py -i input.vcf -y NC_027914.1 -x NC_041774.1
```
This assumes that all the samples in an input vcf file are males. If an input
 vcf file contains not only males but also females, use `-s` option.
## Author
ITO Tsuyoshi

## License
This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

