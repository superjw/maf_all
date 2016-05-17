# This file is used to extract allele frequency of all the SNPs discovered by 1k genome project
# outline
# 1. get AF from 1kgp file
# 2. add AF to 19.rewrite.script.mapping.tsv
# 3. recalculate sum of maf/kb for AF(modify v2sum_of_maf.py)
import os
import re
import sys


# the below section is used for generating dict files
# chromosome = sys.argv[1]
# errorfile = open('errorlines.txt', 'w')
# dict_file = open('./tmp/pos_af_dict_file' + str(chromosome) + '.tmp','w')
# expr = r'^\d{1,2}\t(\d+)\t\S+\t.*;AF=(.+);AN.*EAS_AF=(.*);AMR_AF=(.*);AFR_AF=(.*);EUR_AF=(.*);SAS_AF=(.*\d);'
# for line in sys.stdin:
#     if not line.startswith('#'):
#         # print(line)
#         try:
#             m = re.search(expr, line)
#             id = m.group(1)
#             af = m.group(2)
#             # d[rs] = af
#             # print(id, end='\t')
#             # print(af)
#             dict_file.write(id + '\t' + af + '\n')
#         except AttributeError:
#             print(line)
#             errorfile.write(line)
# errorfile.close()
# dict_file.close()
# print(str(chromosome) + 'done!')

# build dictionary

# chromosome = str(sys.argv[1])
# dict_source_file = './tmp/pos_af_dict_file' + str(chromosome) + '.tmp'
# dict = {}
# with open(dict_source_file, 'r') as f:
#     for line in f:
#         line = line.strip().split('\t')
#         dict[line[0]] = line[1]

# add af to the 19.rewrite.script.mapping.tsv filed
def build_dict(chromosome):
    d = {}
    with open('./tmp/pos_af_dict_file' + str(chromosome) + '.tmp', 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            d[line[0]] = line[1]
    return d


chromosome = sys.argv[1]
outfile = open('./' + chromosome + 'mapped_snp_af_1k_flank_.tsv', 'a')
header = 'chr\tpos\trsid\teas\tamr\tafr\teur\sas\tgid\tgname\taf\n'
dictionary = build_dict(chromosome)
with open('../gwas/' + str(chromosome) + '19.rewrite.script.mapping.tsv', 'r') as f:
    for line in f:
        line = line.replace('\n', '\t')
        pos = line.split('\t', 2)[1]
        af = dictionary.get(pos)
        outfile.write(line + af +'\n')
outfile.close()

