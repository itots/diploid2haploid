#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The script converts diploid genotypes to haploid genotypes for the sex
chromosomes of males in vcf file.
"""

import argparse

__author__ = "ITO Tsuyoshi"
__version__ = "0.1.0"
__email__ = "ito.tsuyoshi.3a@kyoto-u.ac.jp"
__date__ = "2020-03-25"


def main():
    sex = "../../list/sexlist.txt"
    input_vcf = "../../populations_y_vcf/populations.snps.vcf"
    output_vcf = "../../populations_y_vcf/haploid.vcf"
    ychrom = "NC_027914.1"
    xchrom = "NC_041774.1"
    female = False

    args = get_arguments()

    if args.female:
        f_ids = []  # The id of females
        m_ids = []  # The id of males
        with open(args.sex) as lines:
            for line in lines:
                id_sex = line.split('\n')[0].split('\t')
                if "female" in id_sex[1]:
                    f_ids += [id_sex[0]]
                elif "male" in id_sex[1]:
                    m_ids += [id_sex[0]]

    with open(args.input_vcf, 'r') as input_vcf, \
            open(args.output_vcf, 'w') as output_vcf:

        het_y = 0
        hom_y = 0
        mis_y = 0
        het_x = 0
        hom_x = 0
        mis_x = 0

        het_loci_y = 0
        het_loci_x = 0

        for line in input_vcf:
            if '#CHROM' in line:
                columns = line.split('\n')[0].split('\t')

                if args.female:
                    # This checks whether the number of samples matches between
                    # vcf file and sex assignment file
                    if len(f_ids + m_ids) != len(columns[9:]):
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

            # Y chromosome
            elif args.ychrom in line:

                new_sample = []
                het_y_per_sample = 0
                for i in range(len(line.split('\n')[0].split('\t'))):

                    if i in m_indexes:
                        gt = line.split('\n')[0].split('\t')[i].split(':')[
                            0].split('/')

                        if gt[0] != gt[1]:
                            new_gt = '.'
                            het_y += 1
                            het_y_per_sample += 1

                        else:
                            new_gt = gt[0]
                            if gt[0] == '.':
                                mis_y += 1
                            else:
                                hom_y += 1

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

                if het_y_per_sample != 0:
                    het_loci_y += 1

            elif args.xchrom in line:

                new_sample = []
                het_x_per_sample = 0
                for i in range(len(line.split('\n')[0].split('\t'))):
                    if i in m_indexes:
                        gt = line.split('\n')[0].split('\t')[i].split(':')[
                            0].split('/')

                        if gt[0] != gt[1]:
                            new_gt = '.'
                            het_x += 1
                            het_x_per_sample += 1

                        else:
                            new_gt = gt[0]
                            if gt[0] == '.':
                                mis_x += 1
                            else:
                                hom_x += 1

                        new_sample += [
                            ':'.join(
                                [new_gt] + line.split('\n')[0].split('\t')[
                                               i].split(':')[1:]
                            )
                        ]

                    else:
                        new_sample += [line.split('\n')[0].split('\t')[i]]

                if het_x_per_sample != 0:
                    het_loci_x += 1

                new_line = '\t'.join(new_sample)
                output_vcf.write(new_line + '\n')

            else:
                output_vcf.write(line)

    # This shows the number of snps with hetero, homo, and missing
    print("# hetero snps Y\t" + str(het_y) + " in " + str(het_loci_y) + " "
                                                                        "loci")
    print("# homo snps Y\t" + str(hom_y))
    print("# missing snps Y\t" + str(mis_y))
    print("# hetero snps X\t" + str(het_x) + " in " + str(het_loci_x) + " "
                                                                        "loci")
    print("# homo snps X\t" + str(hom_x))
    print("# missing snps X\t" + str(mis_x))


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
                        help="The name of output vcf file. By default, "
                             "./output.vcf")

    parser.add_argument("-f",
                        "--female",
                        action="store_true",
                        dest="female",
                        required=False,
                        default=False,
                        help="Input vcf file contains females, disabled "
                             "by default")

    parser.add_argument("-s",
                        "--sex",
                        action="store",
                        dest="sex",
                        required=False,
                        help="File of sex assignment, required when input "
                             "file contains females")

    args = parser.parse_args()

    if args.female:
        if args.sex is None:
            raise Exception("Sex assignment file must be specified when "
                            "input vcf file contains females")

    return args


def show_counts():



if __name__ == "__main__":
    main()
