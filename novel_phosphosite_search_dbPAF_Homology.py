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

PAF_df = PAF_df[['Uniprot','Type', 'Position', 'Sequence']].copy()
phospho_df = phospho_df[['Fasta headers', 'Proteins', 'Localization prob', 'Amino acid', 'Position', 'Sequence window']].copy()

# Remove phosphos below localization threshold
phospho_df = phospho_df[phospho_df['Localization prob'].astype(float) >= localiation_cutoff]
original_phospho_df = phospho_df  # Hold a copy of original df

# Rename "type" to "Amino acid"
PAF_df = PAF_df.rename(columns={'Type': 'Amino acid'})

#Give 15-underscore buffer for sequence window slices
PAF_df['Sequence underscores'] = PAF_df['Sequence']
PAF_df['Sequence underscores'] = 15*'_' + PAF_df['Sequence'].astype(str) + 15*'_'


PAF_df['Sequence window'] = 'NaN'
for index_PAF, row_PAF in PAF_df.iterrows():
    position = int(row_PAF['Position']) + 15  # account for the underscore buffer
    position -= 1  # convert from 1 to 0-based numbering

    sequence_with_spacers = row_PAF['Sequence underscores']
    sequence_window = sequence_with_spacers[position-15:position+16]

    PAF_df.loc[index_PAF, 'Sequence window'] = sequence_window


# Compare phospho_df sequence windows to dbPAF
phospho_comparison_df = pd.concat([phospho_df['Proteins'], phospho_df['Sequence window']], axis=1, keys=['Uniprot', 'Sequence window'])
PAF_comparison_df = pd.concat([PAF_df['Uniprot'], PAF_df['Sequence window']], axis=1, keys=['Uniprot', 'Sequence window'])

phospho_flat_df = phospho_comparison_df[['Sequence window']].join(phospho_comparison_df.Uniprot.str.split(';', expand=True))
phospho_flat_df = pd.melt(phospho_flat_df.reset_index(), ['index', 'Sequence window'], value_name='Uniprot').dropna()

common = phospho_flat_df.merge(PAF_comparison_df)
phospho_comparison_df['Novel'] = 1
phospho_comparison_df.loc[common['index'].values, 'Novel'] = 0

phospho_df = phospho_comparison_df
phospho_df = phospho_df.rename(columns={'Uniprot': 'Proteins'})


# Add back in missing columns
phospho_df['Amino acid'] = original_phospho_df['Amino acid']
phospho_df['Position'] = original_phospho_df['Position']

# Save merged dataframe
current_time = time.localtime()
time_string = time.strftime("%m-%d-%y %H_%M", current_time)

phospho_df.to_excel('novel_phospho_sites_dbPAF_by_homology' + time_string + '.xlsx', sheet_name='Sheet1', index=False)
print('\nFile saved as novel_phospho_sites_dbPAF_by_homology' + time_string + '.xlsx')

end_secs = time.time()
runsecs = end_secs - start_secs
print('\n Took ' + str(runsecs) + ' seconds')


