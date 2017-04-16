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

PSP_df = PSP_df[['SITE_+/-7_AA']].copy()
phospho_df = phospho_df[['Fasta headers', 'Proteins', 'Localization prob', 'Amino acid', 'Position', 'Sequence window']].copy()


# Remove phosphos below localization threshold
phospho_df = phospho_df[phospho_df['Localization prob'].astype(float) >= localiation_cutoff]
original_phospho_df = phospho_df  # Hold a copy of original df


# Capitalize AAs in "SITE_7+/-_AA"
PSP_df['Sequence_window_small'] = PSP_df['SITE_+/-7_AA'].str.upper()

# Truncate MaxQuant +/-15 window to +/-7
phospho_df['Sequence_window_small'] = phospho_df['Sequence window'].str[8:-8]

phospho_df['is_in_PSP'] = phospho_df.Sequence_window_small.isin(PSP_df.Sequence_window_small)
phospho_df['Novel'] = -phospho_df['is_in_PSP'].astype(bool)
phospho_df['Novel'] = phospho_df['Novel'].astype(int)


# Save merged dataframe
current_time = time.localtime()
time_string = time.strftime("%m-%d-%y %H_%M", current_time)

phospho_df.to_excel('novel_phospho_sites_PSP_by_homology' + time_string + '.xlsx', sheet_name='Sheet1', index=False)
print('\nFile saved as novel_phospho_sites_PSP_by_homology' + time_string + '.xlsx')

end_secs = time.time()
runsecs = end_secs - start_secs
print('\n Took ' + str(runsecs) + ' seconds')


