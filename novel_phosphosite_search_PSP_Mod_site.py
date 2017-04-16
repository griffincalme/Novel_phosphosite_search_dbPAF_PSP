import pandas as pd
import time


# Need to open original file, filter out non class1
phospho_file = input('Enter phospho filepath: (default: Phospho (STY)Sites.txt) ') or 'Phospho (STY)Sites.txt'
PSP_dataset_file = input('Enter PhosphoSite Plus dataset: (default: Phosphorylation_site_dataset.xlsx) ') or 'Phosphorylation_site_dataset.xlsx'
localiation_cutoff = float(input('Enter Localization prob cutoff: (default: .75) ') or .75)


if phospho_file.endswith('.txt'):
    phospho_df = pd.read_table(phospho_file, dtype=object)
elif phospho_file.endswith('.xlsx'):
    phospho_df = pd.read_excel(phospho_file)
elif phospho_file.endswith('.csv'):
    phospho_df = pd.read_csv(phospho_file)
else:
    raise Exception('Please use tab-delimited (.txt), .xlsx or .csv')


if PSP_dataset_file.endswith('.txt'):
    PSP_df = pd.read_table(PSP_dataset_file, dtype=object)
elif PSP_dataset_file.endswith('.xlsx'):
    PSP_df = pd.read_excel(PSP_dataset_file)
elif PSP_dataset_file.endswith('.csv'):
    PSP_df = pd.read_csv(PSP_dataset_file)
else:
    PSP_df = pd.read_table(PSP_dataset_file, dtype=object)


start_secs = time.time()
print('\nAll reversed "REV" peptides will be removed\n')


def remove_REV(phospho_df):
    # Remove any reversed peptide
    phospho_df = phospho_df[phospho_df['Protein'].str.contains('REV__') == False]
    phospho_df = phospho_df.reset_index(drop=True)

    return phospho_df


phospho_df = remove_REV(phospho_df)

PSP_df = PSP_df[['ACC_ID','MOD_RSD']].copy()
phospho_df = phospho_df[['Fasta headers', 'Proteins', 'Localization prob', 'Amino acid', 'Position']].copy()


# Remove phosphos below localization threshold
phospho_df = phospho_df[phospho_df['Localization prob'].astype(float) >= localiation_cutoff]
original_phospho_df = phospho_df  # Hold a copy of original df

# Label phosphosite position like PSP's MOD column
phospho_df['MOD_RSD'] = phospho_df['Amino acid'] + phospho_df['Position'] + '-p'

# Drop redundant letter and mod sites
phospho_df = phospho_df.drop('Amino acid', 1)
phospho_df = phospho_df.drop('Position', 1)



# Compare phospho_df sequence windows to dbPSP
phospho_comparison_df = pd.concat([phospho_df['Proteins'], phospho_df['MOD_RSD']], axis=1, keys=['ACC_ID', 'MOD_RSD'])
PSP_comparison_df = pd.concat([PSP_df['ACC_ID'], PSP_df['MOD_RSD']], axis=1, keys=['ACC_ID', 'MOD_RSD'])

phospho_flat_df = phospho_comparison_df[['MOD_RSD']].join(phospho_comparison_df.ACC_ID.str.split(';', expand=True))
phospho_flat_df = pd.melt(phospho_flat_df.reset_index(), ['index', 'MOD_RSD'], value_name='ACC_ID').dropna()

common = phospho_flat_df.merge(PSP_comparison_df)
phospho_comparison_df['Novel'] = 1
phospho_comparison_df.loc[common['index'].values, 'Novel'] = 0

phospho_df = phospho_comparison_df
phospho_df = phospho_df.rename(columns={'ACC_ID': 'Proteins'})



# Save merged dataframe
current_time = time.localtime()
time_string = time.strftime("%m-%d-%y %H_%M", current_time)

phospho_df.to_excel('novel_phospho_sites_PSP_by_site' + time_string + '.xlsx', sheet_name='Sheet1', index=False)
print('\nFile saved as novel_phospho_sites_PSP_by_site' + time_string + '.xlsx')

end_secs = time.time()
runsecs = end_secs - start_secs
print('\n Took ' + str(runsecs) + ' seconds')


