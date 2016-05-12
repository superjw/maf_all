#!/usr/bin/python3
"""
this is a rewrite version of maf calculation script(08/Apr/2016)
# exmaple usage: python3 v2sum_of_maf.py 22 1000
output file: v2.avg.maf.each.gene.chr22.tsv
"""
import os
import sys
import re
from decimal import *
# getcontext().prec = 5


def gene_length_dict(chromosome, flank):
    """
    generate {gene, length} dictionary based on chromosome
    flank length was considered
    :param chromosome: chromosome number
    :param flank: flank length
    :return: dictionary {gene, length}
    """
    file_name = os.path.abspath('/home/nis/jingwei/work/1k_analysis/gene_location/pos.chr' +
                                str(chromosome) + '.gene.tsv')
    d = {}
    with open(file_name, 'r') as f:
        for line in f:
            line_lst = line.strip().split('\t')
            g_id = line_lst[0]
            g_length = int(line_lst[2]) - int(line_lst[1]) + 2*int(flank)
            d[g_id] = g_length
    return d


def set_of_gene(infile):
    """
    generate a list of genes located on one chromosome based \
    on the inputted file object
    :param: string, infile name
    :return: list, a list of mapped genes on the current chromosome
    """
    f = open(infile, 'r')
    all_gene_id = []
    myre = re.compile(r"ENSG\d+")
    # gene = myre.findall(f)
    for line in f:
        # print(line)
        gene = myre.findall(line)
        # print(gene)
        for g in gene:
                # print(g)
                # if not g:
            all_gene_id.append(g)
    all_gene_id = filter(None, set(all_gene_id))
    f.close()
    # print(chrome + 'list building finished!')
    return all_gene_id


def deal_biallelic(l):
    """
    combine all the comma separated values into one value
    :param l:list, input list
    :return:list, list in the same order
    """
    bia = []
    for a in l:
        # print(a)
        the_sum = Decimal(0)
        if ',' in a:
            a = a .split(',')
            for i in a:
                the_sum += Decimal(i)
                # bia.append(the_sum)
        else:
            the_sum += Decimal(a)
            # print(type(a))
        bia.append(the_sum)
    return bia


def maf(gene_id, infile_name):
    """
   check if gene in the line, and sum up the maf
    :param gene_id: ensemble gene_id
    :param infile_name: mapped file with MAF in, generated from zcat | tt3.py
    :return: decimal, maf of each superpopulation
    """
    with open(infile_name, 'r') as f:
        eas = amr = afr = eur = sas = Decimal(0)
        # i = Decimal(0)
        for line in f:
            if gene_id in line:
                line_lst = line.strip().split('\t')
                new_maf = deal_biallelic(line_lst[3:8])
                eas += new_maf[0]
                amr += new_maf[1]
                afr += new_maf[2]
                eur += new_maf[3]
                sas += new_maf[4]
                # i += Decimal(1)
            else:
                pass
        # print((eas/i).quantize(Decimal('0.0001')))
        # return [(eas/length).quantize(Decimal('0.00001')),
        #         (amr/length).quantize(Decimal('0.00001')),
        #         (afr/length).quantize(Decimal('0.00001')),
        #         (eur/length).quantize(Decimal('0.00001')),
        #         (sas/length).quantize(Decimal('0.00001'))]
        return eas, amr, afr, eur, sas


def main():
    # exmaple usage: python3 v2sum_of_maf.py 22 1000
    chromosome = sys.argv[1]
    flank = sys.argv[2]
    infile = str(chromosome) + '.rewrite.script.mapping.tsv'
    # the infile is the mapping result file, generated from the zcat | tt3.py
    outf = open('v2.avg.maf.each.gene.chr' + str(chromosome) + '.tsv', 'w')
    gene_list = set_of_gene(infile)
    # print(gene_list)
    length_dict = gene_length_dict(chromosome, flank)
    for gene_id in gene_list:
        length = Decimal(length_dict.get(gene_id))
        eas, amr, afr, eur, sas = maf(gene_id, infile)
        avg_eas = (eas / length).quantize(Decimal('0.00001'))
        avg_amr = (amr / length).quantize(Decimal('0.00001'))
        avg_afr = (afr / length).quantize(Decimal('0.00001'))
        avg_eur = (eur / length).quantize(Decimal('0.00001'))
        avg_sas = (sas / length).quantize(Decimal('0.00001'))
        # print(gene_id + '\t' +
        #       str(avg_eas) + '\t' +
        #       str(avg_amr) + '\t' +
        #       str(avg_afr) + '\t' +
        #       str(avg_eur) + '\t' +
        #       str(avg_sas))
        outf.write(gene_id + '\t' +
              str(avg_eas) + '\t' +
              str(avg_amr) + '\t' +
              str(avg_afr) + '\t' +
              str(avg_eur) + '\t' +
              str(avg_sas) + '\n')
    outf.close()
    print('job done! ' + chromosome)


main()
