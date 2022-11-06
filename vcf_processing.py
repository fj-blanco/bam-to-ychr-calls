import os
from datetime import datetime
import pandas as pd
from pysam import VariantFile
from functions.functions import *

sample = "sample"
vcf_route = "./genomes/" + sample + "_Y"
vcf = VariantFile(vcf_route + ".vcf")

# Changing to hg38 coordinates and annotating the problematic regions and positions:
print(datetime.now(), "Changing to hg38 coordinates and annotating the problematic regions and positions (1/3)...", flush=True)
vc_list = change_to_hg38_and_annotate_positions(vcf)

# Generating a dataframe from the .vcf file with the annotated positions:
df = pd.DataFrame(vc_list, columns =['Position', 'Ref', 'Alt', 'Genotype', 'problematic_mutation', 'problematic_region', 'Qual'], dtype = str)
df = df.rename(columns={'Position': 'start'})
df = df.astype({"start": int})

# Generating the annotation dataframe from the sample positions:
print(datetime.now(), "2/3. Generating the annotation dataframe from the sample positions (2/3)...", flush=True)
df_anotation = generate_annotation_df(df)

# Mixing and saving the output:
print(datetime.now(), "Mixing (3/3)...", flush=True)
df_final = pd.merge(df, df_anotation, on='start')
# Ybrowse file contains some positions with a different reference allele (perhaps bacmutations?)
# We remove those positions:
df_final = df_final[df_final['Ref'] == df_final['allele_anc']]
df_final = df_final.sort_values(by=['Qual'], ascending=False)
df_final.to_csv(vcf_route + "_annotated_" + ".csv", encoding='utf-8', index=False)