import argparse
from datetime import datetime
import pandas as pd
from pysam import VariantFile
from functions.functions import *

def main(vcf_route, all_calls, quality_threshold):

    print(datetime.now(), "Loading vcf file (1/3)...", flush=True)
    vcf = VariantFile(vcf_route)

    # Changing to hg38 coordinates and annotating the problematic regions and positions:
    print(datetime.now(), "Changing to hg38 coordinates and annotating the problematic regions and positions (2/3)...", flush=True)
    vc_list = change_to_hg38_and_annotate_positions(vcf)

    # Generating a dataframe from the .vcf file with the annotated positions:
    print(datetime.now(), "Generating the annotation dataframe from the sample positions (3/3)...", flush=True)
    df = pd.DataFrame(vc_list, columns =['Position', 'Ref', 'Alt', 'Genotype', 'problematic_mutation', 'problematic_region', 'Qual'], dtype = str)
    df = df.rename(columns={'Position': 'start'})
    df = df.astype({"start": int})

    if not all_calls:
        # Filtering rows where 'Alt' is not null
        df = df[df['Alt'].notna()]

    # Generating the annotation dataframe from the sample positions:
    df_final = generate_annotation_df(df)

    # Mixing and saving the output:
    print(datetime.now(), "Finalizing (3/3)...", flush=True)

    # Filtering by 'Qual' column and sorting:
    df_final['Qual'] = pd.to_numeric(df_final['Qual'], errors='coerce')
    if quality_threshold is not None:
        df_final = df_final[df_final['Qual'] >= quality_threshold]
    df_final = df_final.sort_values(by=['Qual'], ascending=False)

    if all_calls:
        df_final.to_csv(vcf_route.replace(".vcf", "") + "_annotated_all_variants" + ".csv", encoding='utf-8', index=False)
        print(f'Annotated variants saved in {vcf_route.replace(".vcf", "") + "_annotated_all_variants" + ".csv"}')
        # Filtering rows where 'Alt' is not null
        df_final = df_final[df_final['Alt'].notna()]
    
    # Saving annotated SNPs
    df_final.to_csv(vcf_route.replace(".vcf", "") + "_annotated_SNPs.csv", encoding='utf-8', index=False)
    print(f'Annotated SNPs saved in {vcf_route.replace(".vcf", "") + "_annotated_SNPs.csv"}')
    print("Annotation finished!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process VCF and Annotate variants.')
    parser.add_argument('-v', '--vcf_route', required=True, help='VCF route (path without .vcf extension)')
    parser.add_argument('-a', '--all_calls', action='store_true', help='Annotate all calls (default is False)')
    parser.add_argument('-q', '--quality_threshold', type=float, default=None, help='Quality threshold for filtering (default is None)')

    args = parser.parse_args()

    main(args.vcf_route, args.all_calls, args.quality_threshold)