# Read me  
## files  
* `\*.gid.maf.af.tsv`: are files contain gene_id, MAF/kb in 5 major population in 1kG project and AF/kb for each gene(in ENSG_ids)  
* `\*mapped_snp_af_1k_flank_.tsv`: are files that contain MAF, AF for each snps in each one of the five major populations; and the mapped genes based on 1kb flank length. They are source files for use in `sum_af_maf.py` to generate `\*gid.amf.af.tsv` files    
* `gwas_catalog_trait-mappings_r2016-04-24.tsv`: exported form GWAS-Catalog. It contains the info how the reported traits are mapped to the mapped traits.   
* `gwas_catalog_v1.0.1-associations_e84_r2016-04-24.tsv`: association file downloaded from GWAS-Catalog website.  
* `mart_export_gid_gname_37.txt`: ENSG_ids and corresponded gene name  
* `no_of_mapped_trait_sorted.tsv`:  
* `no_of_mapped_trait.tsv`: number of mapped traits for each gene name  

## scripts  
* `sum_af_maf.py`:  is a modified version of `v2sum_of_maf.py`, which is used in calculating MAF/kb and AF/kb. 
* `af.bash` : is used to run `sum_af_maf.py` in a batch for chromosome 1\-22  
* `af_extraction.py`:  
* `generate_dictfile.bash`:
* `demo.py`:  
* `add_af_to_file.bash`: 
