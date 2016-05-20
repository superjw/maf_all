# scripts that used to generate correlation analysis in R
# 1. build dictionary from `mart_export_gid_gname_37.txt` {ENSG: g_name}
# 2.build dictionary from 'no_of_mapped_traits.tsv' {g_name: no_of_m_traits}
# 3. read '\*.gid.maf.af.tsv' file: use the 1st column get g_name from dictionary in step 1
# 4. use values returned from step 4(g_name)as key, get values from dictionary in step 2
# input files: `mart_export_gid_gname_37.txt`, 'no_of_mapped_traits.tsv',  '\*.gid.maf.af.tsv'
import sys


def build_dict_gid_name(gid_gname_file_name):
    """
    build a dictionary of {ensg_id: g_name} from mart_exort_gid_gname_37.txt
    :param file_name: the full path to file: mart_export_gid_gname_37.txt
    :return d: a dictionary {ensg-id: g_name}
    """
    with open(gid_gname_file_name, 'r') as f:
        next(f)
        d = {}
        for line in f:
            g_id = line.strip().split('\t')[0]
            g_name = line.strip().split('\t')[1]
            d[g_id] = g_name
        return d


def build_dict_name_mtratis(no_of_m_traits_file):
    """

    :param no_of_m_traits_file:
    :return d:
    """
    with open(no_of_m_traits_file, 'r') as f:
        next(f)
        d = {}
        for line in f:
            g_name = line.strip().split('\t')[0]
            no_of_mapped_traits = line.strip().split('\t')[2]
            d[g_name] = no_of_mapped_traits
        return d



def main(chromosome):
    gene_id_name_37 = 'mart_export_gid_gname_37.txt'
    no_of_mapped_trait_file = 'no_of_mapped_trait.tsv'
    id_name_dict = build_dict_gid_name(gene_id_name_37)
    name_m_traits = build_dict_name_mtratis(no_of_mapped_trait_file)
    af_file_name = str(chromosome) + '.gid.maf.af.tsv'
    outfile = open(chromosome + '_id_af_m_traits.tsv', 'w')
    with open(af_file_name, 'r') as f:
        header = next(f).strip() + '\tgene_name\tno_of_m_traits\n'
        outfile.write(header)
        for line in f:
            gid = line.strip().split('\t')[0]
            g_name = id_name_dict.get(gid)
            no_m_traits = name_m_traits.get(g_name)
            outfile.write(line.strip() + '\t' + g_name + '\t' + no_m_traits + '\n')
    outfile.close()

chromosome = sys.argv[1]
main(chromosome)
