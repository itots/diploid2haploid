#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The script converts diploid genotypes to haploid genotypes for the sex
 chromosomes of males in a vcf file. The scripts modifies only GT fields
  with others remain. Heterozygous genotypes in the sex chromosomes of males
  , if any, are converted to be missing.
"""

import argparse

__author__ = "ITO Tsuyoshi"
__version__ = "1.1"
__email__ = "ito.tsuyoshi.3a@kyoto-u.ac.jp"
__date__ = "2020-04-01"


def main():

    args = get_arguments()

    with open(args.input_vcf, 'r') as input_vcf, \
            open(args.output_vcf, 'w') as output_vcf:

        y_het = 0
        y_hom = 0
        y_mis = 0
        x_het = 0
        x_hom = 0
        x_mis = 0

        y_het_loci = 0
        x_het_loci = 0

        n_loci = 0
        n_y_loci = 0
        n_x_loci = 0

        for line in input_vcf:
            if not line.startswith("#"):
                n_loci += 1

            if "#CHROM" in line:
                output_vcf.write(line)

                columns = line.split('\n')[0].split('\t')
                n_samples = len(columns[9:])

                if args.sex is not None:

                    f_ids = []  # The id of females
                    m_ids = []  # The id of males
                    with open(args.sex) as sex_lines:
                        for sex_line in sex_lines:
                            id_sex = sex_line.split('\n')[0].split('\t')
                            if "female" in id_sex[1]:
                                f_ids += [id_sex[0]]
                            elif "male" in id_sex[1]:
                                m_ids += [id_sex[0]]

                    # This checks whether the number of samples matches between
                    # vcf file and sex assignment file

                    if len(f_ids + m_ids) != n_samples:
                        raise Exception("The name and number of samples must "
                                        "be consistent between input vcf "
                                        "file and sex assignment file")

                    # This checks whether the name of samples matches between
                    # vcf file and sex assignment file
                    elif set(f_ids + m_ids) != set(columns[9:]):
                        raise Exception("The name and number of samples must "
                                        "be consistent between input vcf "
                                        "file and sex assignment file")

                    else:
                        m_indexes = list(
                            map(lambda x: columns.index(
                                x) if x in columns else None,
                                m_ids))

                else:
                    # The samples begins at tenth columns in vcf
                    m_indexes = list(range(len(columns)))[9:]

                n_males = len(m_indexes)

            # Y chromosome
            elif args.ychrom in line:

                n_y_loci += 1

                new_sample = []
                y_het_sample = 0
                for i in range(len(line.split('\n')[0].split('\t'))):

                    if i in m_indexes:
                        gt = line.split('\n')[0].split('\t')[i].split(':')[
                            0].split('/')

                        if gt[0] != gt[1]:
                            new_gt = '.'
                            y_het += 1
                            y_het_sample += 1

                        else:
                            new_gt = gt[0]
                            if gt[0] == '.':
                                y_mis += 1
                            else:
                                y_hom += 1

                        new_sample += [
                            ':'.join(
                                [new_gt] + line.split('\n')[0].split('\t')[
                                               i].split(':')[1:]
                            )
                        ]

                    else:
                        new_sample += [line.split('\n')[0].split('\t')[i]]

                new_line = '\t'.join(new_sample)
                output_vcf.write(new_line + '\n')

                if y_het_sample != 0:
                    y_het_loci += 1

            # X chromosome
            elif args.xchrom in line:

                n_x_loci += 1

                new_sample = []
                x_het_sample = 0
                for i in range(len(line.split('\n')[0].split('\t'))):
                    if i in m_indexes:
                        gt = line.split('\n')[0].split('\t')[i].split(':')[
                            0].split('/')

                        if gt[0] != gt[1]:
                            new_gt = '.'
                            x_het += 1
                            x_het_sample += 1

                        else:
                            new_gt = gt[0]
                            if gt[0] == '.':
                                x_mis += 1
                            else:
                                x_hom += 1

                        new_sample += [
                            ':'.join(
                                [new_gt] + line.split('\n')[0].split('\t')[
                                               i].split(':')[1:]
                            )
                        ]

                    else:
                        new_sample += [line.split('\n')[0].split('\t')[i]]

                if x_het_sample != 0:
                    x_het_loci += 1

                new_line = '\t'.join(new_sample)
                output_vcf.write(new_line + '\n')

            else:
                output_vcf.write(line)

    # Output summary
    output_prefix = args.output_vcf.split(".vcf")[0]
    with open(output_prefix + "_summary.txt", "w") as out_sum:
        out_sum.write(str(n_males)
                      + " males in "
                      + str(n_samples)
                      + " total samples\n\n")

        out_sum.write(str(n_loci) + " loci in total\n\n")

        out_sum.write("Y chromosome:\n")
        out_sum.write("\t" + str(n_y_loci) + " loci in total\n")
        out_sum.write("\t" + str(y_het_loci) + " loci contain hetero snps\n\n")
        out_sum.write("\t" + str(y_het) + " snps are hetero\n")
        out_sum.write("\t" + str(y_hom) + " snps are homo\n")
        out_sum.write("\t" + str(y_mis) + " snps are missing\n\n")

        out_sum.write("X chromosome:\n")
        out_sum.write("\t" + str(n_x_loci) + " loci in total\n")
        out_sum.write("\t" + str(x_het_loci) + " loci contain hetero snps\n\n")
        out_sum.write("\t" + str(x_het) + " snps are hetero\n")
        out_sum.write("\t" + str(x_hom) + " snps are homo\n")
        out_sum.write("\t" + str(x_mis) + " snps are missing\n")


def get_arguments():

    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("-i",
                        "--input_vcf",
                        action="store",
                        dest="input_vcf",
                        required=True,
                        help="Name of input vcf file")

    parser.add_argument("-y",
                        "--ychrom",
                        action="store",
                        dest="ychrom",
                        required=True,
                        help="Name of Y chromosome in reference sequence")

    parser.add_argument("-x",
                        "--xchrom",
                        action="store",
                        dest="xchrom",
                        required=True,
                        help="Name of X chromosome in reference sequence")

    parser.add_argument("-o",
                        "--output_vcf",
                        action="store",
                        dest="output_vcf",
                        required=False,
                        default="./output.vcf",
                        help="Name of output vcf file. By default, it is "
                             "./output.vcf")

    parser.add_argument("-s",
                        "--sex",
                        action="store",
                        dest="sex",
                        required=False,
                        help="File of sex assignment, required when input "
                             "vcf file contains females")

    args = parser.parse_args()

    if args.sex is not None:
        print("Sex assignment file was specified. ")
    else:
        print("All samples are assumed to be males. Please specify sex "
              "assignment file if input vcf file contains not only "
              "males but also females.")

    return args


if __name__ == "__main__":
    main()
