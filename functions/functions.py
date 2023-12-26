import pandas as pd
from pyliftover import LiftOver
import numpy as np

# Regions to annotate:
par1 = [1,2781479]
cen = [10072350,11686750]
dyz19 = [20054914,20351054]
ppr = [26637971,26673210]
par2 = [56887903,57217415]

def check_region(position):
    if (par1[0] <= position <= par1[1]):
        i_reg = 'par1'
    elif (cen[0] <= position <= cen[1]):
        i_reg = 'cen'
    elif (dyz19[0] <= position <= dyz19[1]):
        i_reg = 'dyz19'
    elif (ppr[0] <= position <= ppr[1]):
        i_reg = 'ppr'
    elif (par2[0] <= position <= par2[1]):
        i_reg = 'par2'
    else:
        i_reg = ''
    return i_reg


def check_mutation(reference, genotype):
    i_mut = ''
    if (reference == 'C') & (genotype == 'T'):
        i_mut = 'C to T deamination';
    elif (reference == 'G') & (genotype == 'A'):
        i_mut = 'G to A deamination'
    return i_mut

def change_to_hg38_and_annotate_positions(vcf):
    vc_list = []
    lo = LiftOver('hg19', 'hg38')
    for rec in vcf.fetch():
        i_pos_0 = lo.convert_coordinate('chrY', int(rec.pos))
        if i_pos_0 != []:
            i_pos = int(i_pos_0[0][1])        
            i_reg = check_region(i_pos)
            i_gen_0 = rec.alts
            i_ref = rec.ref
            i_mut = ''
            if i_gen_0 == None:
                i_gen = i_ref
            else:
                i_gen = i_gen_0[0]
                i_mut = check_mutation(i_ref, i_gen)
            vc_list.append([i_pos, rec.ref, rec.alts, i_gen, i_mut, i_reg, rec.qual])
        else:
            vc_list.append([int(0), rec.ref, rec.alts, rec.alts, '', '', rec.qual])
    return vc_list

def snps_aggregation(group):
    annotations = ['Name_yleaf', 'Name_ybrowse', 'haplogroup_yleaf', 'ISOGG_haplogroup', 'allele_anc', 'allele_der', 'mutation_yleaf', 'mutation_ybrowse']
    result = {'start': group['start'].iloc[0]}
    result.update({col: '/'.join(map(str, np.unique(group[col].dropna().values))) for col in annotations})
    return pd.Series(result)

def load_yleaf(start_positions):
    df_yleaf = pd.read_csv('./SNP_annotations_hg38/yleaf_new_positions.txt', sep="\t", header=None)
    df_yleaf.columns = ['chr', 'Name', 'haplogroup_yleaf', 'start', 'mutation', 'allele_anc', 'allele_der']
    df_yleaf = df_yleaf.drop(columns=['chr'])
    df_yleaf['start'] = df_yleaf['start'].astype(int)
    df_yleaf = df_yleaf[df_yleaf['start'].isin(start_positions)]
    return df_yleaf

def loaf_ybrowse(start_positions):
    df_ybrowse = pd.read_csv('./SNP_annotations_hg38/ybrowse_snps_hg38.csv')
    df_ybrowse = df_ybrowse[['Name', 'ISOGG_haplogroup', 'start', 'allele_anc', 'allele_der', 'mutation']]
    df_ybrowse['start'] = df_ybrowse['start'].astype(int)
    df_ybrowse = df_ybrowse[df_ybrowse['start'].isin(start_positions)]
    return df_ybrowse

def generate_annotation_df(df, df_yleaf=None, df_yleaf_old=None, df_ybrowse=None):
    # Get the unique start positions from the input df
    start_positions = df['start'].astype(int).unique()

    # Read Yleaf
    if df_yleaf is None:
        df_yleaf = load_yleaf(start_positions)
    # Read Ybrowse
    if df_ybrowse is None:
        df_ybrowse = loaf_ybrowse(start_positions)

    # Merge the two dfs
    df_annotations = pd.merge(df_yleaf, df_ybrowse, on=['start', 'allele_anc', 'allele_der'], how='outer', suffixes=('_yleaf', '_ybrowse'))
    # Merge on 'start' to get all possible matches
    df_merged = pd.merge(df, df_annotations, on='start', how='left')
    del df_annotations

    # If there is a mutation in the sample (non-empty "Alt"), filter to matching alleles only
    matching_alleles = (df_merged['Ref'] == df_merged['allele_anc']) & (df_merged['Genotype'] == df_merged['allele_der'])
    non_matching_mask = df_merged['Alt'].notna() & ~matching_alleles
    columns_to_update = ['haplogroup_yleaf', 'Name_yleaf', 'Name_ybrowse', 'ISOGG_haplogroup', 'mutation_yleaf', 'mutation_ybrowse']
    df_merged.loc[non_matching_mask, columns_to_update] = None

    # Find the counts of each start position and split the df into unique and repeated start positions
    start_counts = df_merged['start'].value_counts()
    df_unique_starts = df_merged[df_merged['start'].isin(start_counts.index[start_counts == 1])]
    df_repeated_starts = df_merged[df_merged['start'].isin(start_counts.index[start_counts > 1])]

    # Apply the groupby operation only to the df with repeated start positions
    df_repeated_agg = df_repeated_starts.groupby('start').apply(snps_aggregation).reset_index(drop=True)

    # Concatenate the df with unique start positions and the aggregated df with repeated start positions
    annotations = ['Name_yleaf', 'Name_ybrowse', 'haplogroup_yleaf', 'ISOGG_haplogroup', 'allele_anc', 'allele_der', 'mutation_yleaf', 'mutation_ybrowse']
    df_result = pd.concat([df_unique_starts[['start'] + annotations], df_repeated_agg], ignore_index=True)

    # Merge the original df with the result, keeping all original rows
    df_final = pd.merge(df, df_result, on='start', how='left')
    print("Final merge done")
    
    return df_final

if __name__ == '__main__':
    print('This is a module')