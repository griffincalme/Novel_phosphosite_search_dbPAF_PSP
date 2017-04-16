import pandas as pd
import time


# Need to open original file, filter out non class1
phospho_file = input('Enter phospho filepath: (default: Phospho (STY)Sites.txt) ') or 'Phospho (STY)Sites.txt'
PAF_dataset_file = input('Enter dbPAF phosphosite dataset: (default: RAT.elm) ') or 'RAT.elm'
localiation_cutoff = float(input('Enter Localization prob cutoff: (default: .75) ') or .75)


if phospho_file.endswith('.txt'):
    phospho_df = pd.read_table(phospho_file, dtype=object)
elif phospho_file.endswith('.xlsx'):
    phospho_df = pd.read_excel(phospho_file)
elif phospho_file.endswith('.csv'):
    phospho_df = pd.read_csv(phospho_file)
else:
    raise Exception('Please use tab-delimited (.txt), .xlsx or .csv')


if PAF_dataset_file.endswith('.txt'):
    PAF_df = pd.read_table(PAF_dataset_file, dtype=object)
elif PAF_dataset_file.endswith('.xlsx'):
    PAF_df = pd.read_excel(PAF_dataset_file)
elif PAF_dataset_file.endswith('.csv'):
    PAF_df = pd.read_csv(PAF_dataset_file)
else:
    PAF_df = pd.read_table(PAF_dataset_file, dtype=object)


start_secs = time.time()
print('\nAll reversed "REV" peptides will be removed\n')


def remove_REV(phospho_df):
    # Remove any reversed peptide
    phospho_df = phospho_df[phospho_df['Protein'].str.contains('REV__') == False]
    phospho_df = phospho_df.reset_index(drop=True)

    return phospho_df


phospho_df = remove_REV(phospho_df)


PAF_df = PAF_df[['Uniprot','Type', 'Position']].copy()
phospho_df = phospho_df[['Fasta headers', 'Proteins', 'Localization prob', 'Amino acid', 'Position', 'Modified sequence']].copy()

# Remove phosphos below localization threshold
phospho_df = phospho_df[phospho_df['Localization prob'].astype(float) >= localiation_cutoff]


# Combine amino acid letter and modification site for easier cross-comparison
PAF_df['Modification site'] = PAF_df['Type'] + PAF_df['Position'].astype(str)
phospho_df['Modification site'] = phospho_df['Amino acid'] + phospho_df['Position'].astype(str)

# Drop redundant letter and mod sites
PAF_df = PAF_df.drop('Type', 1)
PAF_df = PAF_df.drop('Position', 1)
phospho_df = phospho_df.drop('Amino acid', 1)
phospho_df = phospho_df.drop('Position', 1)



#total_phospho_rows = str(len(phospho_df.index))

phospho_comparison_df = pd.concat([phospho_df['Proteins'], phospho_df['Modification site']], axis=1, keys=['Uniprot', 'Modification site'])
PAF_comparison_df = pd.concat([PAF_df['Uniprot'], PAF_df['Modification site']], axis=1, keys=['Uniprot', 'Modification site'])

phospho_flat_df = phospho_comparison_df[['Modification site']].join(phospho_comparison_df.Uniprot.str.split(';', expand=True))
phospho_flat_df = pd.melt(phospho_flat_df.reset_index(), ['index', 'Modification site'], value_name='Uniprot').dropna()

common = phospho_flat_df.merge(PAF_comparison_df)
phospho_comparison_df['Novel'] = 1
phospho_comparison_df.loc[common['index'].values, 'Novel'] = 0

phospho_df = phospho_comparison_df
phospho_df = phospho_df.rename(columns={'Uniprot': 'Proteins'})


# Save merged dataframe
current_time = time.localtime()
time_string = time.strftime("%m-%d-%y %H_%M", current_time)

phospho_df.to_excel('novel_phospho_sites_by_site_index' + time_string + '.xlsx', sheet_name='Sheet1', index=False)
print('File saved as novel_phospho_sites_by_site_index' + time_string + '.xlsx')



end_secs = time.time()
runsecs = end_secs - start_secs
print('\n Took ' + str(runsecs) + ' seconds')

