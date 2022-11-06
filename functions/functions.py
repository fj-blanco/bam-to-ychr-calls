import pandas as pd
from pyliftover import LiftOver

# Regions to annotate:
par1 = [1,2781479]
cen = [10072350,11686750]
dyz19 = [20054914,20351054]
ppr = [26637971,26673210]
par2 = [56887903,57217415]

# mutations to annotate:
mutations_to_annotate = ["C to T", "A to T", "G to C", "A to G"]

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
    elif (reference == 'A') & (genotype == 'T'):
        i_mut = 'A to T risk of strand misidentification'
    elif (reference == 'G') & (genotype == 'C'):
        i_mut = 'G to C risk of strand misidentification'
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

def generate_annotation_df(df):
    # Yleaf:
    df_yleaf = pd.read_csv('SNP_annotations_hg38/yleaf_new_positions.txt', sep="\t", header=None)
    df_yleaf.columns = ['chr', 'Name', 'haplogroup', 'start', 'mutation', 'allele_anc', 'allele_der'];
    #Ybrowse:
    df_ybrowse = pd.read_csv('SNP_annotations_hg38/ybrowse_snps_hg38.csv')[(['seqid', 'Name', 'ISOGG_haplogroup', 'start', 'mutation', 'allele_anc', 'allele_der'])];
    df_ybrowse = df_ybrowse.rename(columns={'seqid': 'chr'});
    # We  select all the Yleaf annotated positions found in the sample 
    # and then we add the Ybrowse annotated positions that are in the sample but not in Yleaf file:
    # giving thus priority to Yleaf annotations:
    df_int_yleaf = df_yleaf[df_yleaf["start"].isin(df["start"].tolist())]
    df_int_ybrowse_no_yleaf = df_ybrowse[(df_ybrowse["start"].isin(df["start"].tolist())) & ~(df_ybrowse["start"].isin(df_int_yleaf["start"].tolist()))]
    df_anotation = pd.concat([df_int_yleaf, df_int_ybrowse_no_yleaf])
    return df_anotation