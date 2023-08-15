#!/usr/bin/env python3

"""
calculate the genebody coverage for RNAseq data

alternative:
1. geneBody_coverage.py (RSeQC)
2. ngsplot

Why this method:
- to speed up (RSeQC takes half hour, even hours for single BAM file)

How to:
- use `bamtocov` (https://github.com/telatin/BamToCov) to calculate the coverage first.
- parse coverage for each transcript (exon-exon) [using multiple threads, each chromosomes] 
- calculate 100 percentiles for each transcript 
- generate plot

"""

# public modules
import os
import pathlib
import argparse
import logging
import numpy as np
from xopen import xopen
from datetime import datetime


logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    stream=sys.stdout)
log = logging.getLogger(__name__)
log.setLevel('INFO')


def cur_time():
    return datetime.now().strftime('%Y-%m-%d %H:%M:%S')


def calculate_percentile_averages(data, num=100):
    if not isinstance(data, np.ndarray):
        log.error('require np.ndarray, failed.')
        return None

    # Sort the data array
    sorted_data = np.sort(data)
    
    # Calculate indices for percentiles
    indices = np.linspace(0, len(sorted_data) - 1, num=num, dtype=int)
    
    # Calculate percentiles and their average values
    percentiles = np.percentile(sorted_data, np.linspace(0, 100, num=num))
    averages = []
    
    for i in range(len(indices) - 1):
        start_idx = indices[i]
        end_idx = indices[i + 1] + 1
        segment = sorted_data[start_idx:end_idx]
        average = np.mean(segment)
        averages.append(average)
    
    return percentiles, averages


## load transcript from GTF file
def parse_gtf(fh, **kwargs):
    """
    Parse genome annotation file (GTF), Ensembl, Gencode
    optional: BED

    Filtering
    3. by biotype (protein_coding) 
    4. by gene_size 

    # for Ensembl-release 102; Gencode release-M25
    # 110 genes were duplicated
    # 
    # keep the smallest gene_version:
    Ptp4a1 Arhgef4 Gm24826 Gm26457 Septin2 Zfp813-ps Gm23370 Gm28040 Zc3h11a 
    Gm28724 Snord80 Gm16701 Ndor1 Ndor1 Snora43 Nron Gm24350 Olfr1073-ps1 
    Gm22813 A530058N18Rik Gm4430 Nkx2-2os Gm23925 1600017P15Rik Gm7270 Gm20690
    Gm25820 4933434E20Rik Gm18433 Terc Tmigd3 Gm13301 Rmrp Gm25053 Pakap Pakap
    Gm26379 Gm26047 Snora16a Gm22897 Gm28710 Jakmip1 Gm27680 Fam220a Gm24105 
    Dlx6os1 Gm16499 Rprl1 Gm27013 Tmem147os Lim2 Mir1839 Olfr290 Aldoa Gm23128
    Gm23604 Gm18645 Gm26265 Dpep2 Gm23377 Gm38642 Zkscan7 Gm16364 Rnu3a
    C730027H18Rik Gm6729 Ddit3 Rnu3b1 Rnu3b3 Rnu3b4 Rnu3b2 Gm22711 Gm26413 
    St6galnac2 1700030C10Rik Gm22149 Gm25203 Gm27528 Gm35558 Gm35558 Gm24022 
    Ighv5-8 Ighv1-13 Vmn1r216 Gm36638 Gm9025 Nnt Vmn2r-ps111 Gm6740 Gm5089 
    4930594M22Rik Mirt2 Gcat Gm27825 Gm28023 Gm41392 Gm38619 Gm25617 Atp5o Gm29719
    Gm5966 Snhg4 Mir1949 Pcdha11 Arhgap26 Gm23927 Zfp91 Gm35438 Gm17522 Gm23786

    # number of genes
    source           ensembl    gencode
    release          102        M25
    date             2020-11    2020-03-24
    genes            55401      55401
    protein_coding   21859      21859

    filtered:
    genes            55367      55291 (remove duplicated gene_names)
    protein_coding   21898      21831 (only protein_coding)
    """
    pass


def read_gtf(gtf, **kwargs):
    """
    Read Ensembl GTF file

    Parameters:
    -----------
    feature : str
        gene|exon|CDS|transcript, default: None
    gene_biotype : str
        TEC|protein_coding, default: None
    Convert GTF to BED format
    filt by: gene_biotype (ENSEMBL), gene_type (GENCODE)

    duplicate gene records found in Ensembl.
    eg: Zfp91, Zkscan7, ...
    filtered by: "gene_version" !!!
    """
    features = kwargs.get('features', ['transcript', 'exon']) # UTR !!
    line = 0 # line number
    with xopen(gtf) as r:
        for l in r:
            line += 1
            p = l.strip().split('\t')
            if len(p) < 9: continue # 
            if l.startswith('#'): continue # comment lines
            # add chr to chromosome name
            # chr = p[0] # if p[0].startswith('chr') else 'chr'+p[0]
            tx_id = parse_gtf_desc(p[8], 'transcript_id')
            gene_name = parse_gtf_desc(p[8], 'gene_name') # column-9
            # gene version in "GENCODE"
            gene_id = parse_gtf_desc(p[8], 'gene_id')
            if '.' in gene_id: # gencode gene_id, ENSG00000237613.2
                gene_id, ver = gene_id.split('.', 1)
                # ver = int(ver)
            else: # Ensembl
                ver = parse_gtf_desc(p[8], 'gene_version')
                ver = 1 if ver is None else int(ver)
            gene_name = gene_id if gene_name is None else gene_name
            # warning: if n is None
            if gene_name is None:
                log.warning(f'line-{line} failed, could not found')
            biotype = parse_gtf_desc(p[8], 'gene_biotype')
            # 6-columns + gene_biotype + gene_version
            yield [p[0], p[3], p[4], tx_id, '255', p[6], gene_id, gene_name, biotype, ver]


def parse_gtf(fh, **kwargs):
    """
    filt GTF records by [feature] and [gene_biotype]

    Parameters:
    -----------
    feature : list
        gene|exon|CDS|transcript, default: None (all)
    gene_biotype : str
        TEC|protein_coding, default: None (all)
    """
    # if not isinstance(fh, io.TextIOWrapper):
    #     log.error('not a io.TextIOWrapper, [fh]')
    #     return None
    feature = kwargs.get('feature', 'transcript')
    gene_biotype = kwargs.get('gene_biotype', 'protein_coding')
    for l in fh:
        p = l.strip().split('\t')
        # filtering
        if len(p) < 9:
            continue
        # filt by feature
        if not feature is None:
            if isinstance(feature, str) and not p[2] == feature:
                continue
            elif isinstance(feature, list) and not p[2] in feature:
                continue
            else:
                pass
        # filt by gene_biotype
        gb = parse_gtf_desc(p[8], 'gene_biotype')
        if not gene_biotype is None:
            if isinstance(gene_biotype, str) and not gb == gene_biotype:
                continue
            elif isinstance(gene_biotype, list) and not gb in gene_biotype:
                continue
            else:
                pass
        yield p


def parse_gtf_desc(s, key='gene_name'):
    """
    read Description column of GTF (column-9), gene_id "ENSMUSG00000118491"; ...
    gene:
    ENSEMBL:
    gene_id "ENSMUSG00000102693"; gene_version "1"; gene_name "4933401J01Rik"; 
    gene_source "havana"; gene_biotype "TEC"; havana_gene "OTTMUSG00000049935";
    havana_gene_version "1";
    GENCODE:
    gene_id "ENSG00000227232.4"; transcript_id "ENSG00000227232.4"; 
    gene_type "pseudogene"; gene_status "KNOWN"; gene_name "WASH7P"; 
    transcript_type "pseudogene"; transcript_status "KNOWN"; 
    transcript_name "WASH7P"; level 2; havana_gene "OTTHUMG00000000958.1";
    """
    d = {}
    if isinstance(s, str):
        for p in s.strip().split(';'):
            if " " in p:
                a, b = p.split()[:2]
                b = b.replace('"', '')
                d.update({a:b})
    # fix gene_type (GENCODE), gene_biotype (ENSEMBL)
    if key == 'gene_type' and 'gene_biotype' in s:
        key = 'gene_biotype'
    elif key == 'gene_biotype' and 'gene_type' in s:
        key = 'gene_type'
    return d.get(key, None)


## calculate the coverages using bamtocov


## calculate the coverages for each transcript (gene_body), strandness


## summary


## plot


# # Example usage
# # data_array = np.array(list(range(1, 1000)))
# data_array = np.array([1] * 1000)
# percentiles, averages = calculate_percentile_averages(data_array)

# # print(len(percentiles))
# for p, avg in zip(percentiles, averages):
#     print(f"Percentile {p:.2f}%: Average = {avg:.2f}")





