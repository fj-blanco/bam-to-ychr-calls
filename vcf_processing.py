import pysam
from pysam import VariantFile
import pandas as pd
from pyliftover import LiftOver

sample = "I13518"
#bcf_route = "./arslantepe/" + sample + "/" + sample + "_Y"
bcf_route = "/home/javi/Documents/genomics/samples/southern_arc/" + sample + "/" + sample + "_Y"

bcf_in = VariantFile(bcf_route + ".vcf")

par1 = [1,2781479]
cen = [10072350,11686750]
dyz19 = [20054914,20351054]
ppr = [26637971,26673210]
par2 = [56887903,57217415]
mutations_to_exclude = ["C to T", "A to T", "G to C", "A to G"]


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
        i_mut = 'C to T deamination'
    elif (reference == 'A') & (genotype == 'T'):
        i_mut = 'A to T risk of strand misidentification'
    elif (reference == 'G') & (genotype == 'C'):
        i_mut = 'G to C risk of strand misidentification'
    elif (reference == 'G') & (genotype == 'A'):
        i_mut = 'G to A deamination'
    return i_mut


convertToHG38 = True
regionsToFilter = ['par1', 'cen', 'dyz19', 'ppr']
mutationsToFilter = ['C to T deamination', 'A to T risk of strand misidentification', 'G to C risk of strand misidentification', 'G to A deamination']
regions2 = ["tal"]
# Cambiando a hg38 (sin excluir y anotando):
vc_list = []
lo = LiftOver('hg19', 'hg38')
for rec in bcf_in.fetch():
    i_pos_0 = lo.convert_coordinate('chrY', int(rec.pos))
    if i_pos_0 != []:
        # Anotamos las posiciones no adecuadas para filogenia:
        i_pos = int(i_pos_0[0][1])
        i_reg = check_region(i_pos)
        # if i_reg not in regions2 list print i_reg:
        if not (i_reg in regionsToFilter):
            # Extraemos el genotipo:
            i_gen_0 = rec.alts
            i_ref = rec.ref
            i_mut = ''
            if i_gen_0 == None:
                i_gen = i_ref
            else:
                i_gen = i_gen_0[0]
                i_mut = check_mutation(i_ref, i_gen)
            if not (i_mut in mutationsToFilter):
                # Creamos la lista:
                vc_list.append([i_pos, rec.ref, rec.alts, i_gen, i_mut, i_reg, rec.qual])
                # Order:
                #['Position', 'Ref', 'Alt', 'Genotype', 'problematic_mutation', 'problematic_region', 'Qual']
    else:
        vc_list.append([int(0), rec.ref, rec.alts, rec.alts, '', '', rec.qual])

# Cambiando a hg38 excluyendo:
vc_list = [];
lo = LiftOver('hg19', 'hg38');
for rec in bcf_in.fetch():
    i_pos_0 = lo.convert_coordinate('chrY', int(rec.pos));
    if i_pos_0 != []:    
        i_pos = int(i_pos_0[0][1]);
        if ~((par1[0] <= i_pos <= par1[1]) | (cen[0] <= i_pos <= cen[1]) | (dyz19[0] <= i_pos <= dyz19[1]) |
           (ppr[0] <= i_pos <= ppr[1]) | (par2[0] <= i_pos <= par2[1])):
            vc_list.append([i_pos, rec.ref, rec.alts, rec.qual]);
    else:
        vc_list.append([int(0), rec.ref, rec.alts, rec.qual]);


# Sin convertir excluyendo:
vc_list = [];
for rec in bcf_in.fetch():
    if rec.chrom:
        i_pos = int(rec.pos)
        if ~((par1[0] <= i_pos <= par1[1]) | (cen[0] <= i_pos <= cen[1]) | (dyz19[0] <= i_pos <= dyz19[1]) |
        (ppr[0] <= i_pos <= ppr[1]) | (par2[0] <= i_pos <= par2[1])):
            vc_list.append([i_pos, rec.ref, rec.alts, rec.qual]);


# Sin convertir:
vc_list = [];
for rec in bcf_in.fetch():
    if rec.chrom:
        i_pos = int(rec.pos);
        vc_list.append([i_pos, rec.ref, rec.alts, rec.qual]);