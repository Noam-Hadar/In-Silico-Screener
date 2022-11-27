import gzip
from glob import glob
import os
import warnings
import shutil
warnings.filterwarnings('ignore')
import sys
sys.stderr = object

temp = 'In-Silico Screener TEMP directory'
ref = None

def skipHeader(vcf):
    skip = 0
    if vcf.endswith('.gz'):
        with gzip.open(vcf, 'r') as f:        
            for line in f:
                if line.decode('utf-8-sig').startswith('##'):
                    skip += 1
                else:
                    return skip
    else:
        with open(vcf) as f:
            for line in f:
                if line.startswith('##'):
                    skip += 1
                else:
                    return skip

vcfs = glob('**/*.vcf', recursive = True)
gz_vcfs = glob('**/*.vcf.gz', recursive = True)
vcfs = vcfs + gz_vcfs

total = len(vcfs)
if total == 0:
    input('No vcf files found in this directory nor subdirectories\nPress the enter key to exit')
    exit()
else:
    print('Found ' + str(total) + ' vcf files.\n')
    print('What reference genome are you using?')
    while ref not in ['19','38']:
        ref = input('Type 19 for hg19 or 38 for hg38 and press enter: ').lower().strip().replace('hg','')
        if ref not in ['19','38']:
            print('Wrong input')

try:
    os.mkdir(temp)
except:
    False

try:
    import pandas as pd
except:
    import subprocess
    print('\nInstalling needed modules for data analysis...')
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'numpy==1.19.3', '--no-warn-script-location'])
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'pandas', '--no-warn-script-location'])
    import pandas as pd

def GenerateByAssembly(assembly):
    df = pd.read_csv('https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz', sep = '\t')
    df = df[['Assembly', 'Type', 'GeneSymbol', 'ClinicalSignificance', 'ReviewStatus', 'PhenotypeList', 'Chromosome', 'Start', 'VariationID', 'ReferenceAlleleVCF', 'AlternateAlleleVCF']]
    df = df[df['Assembly'] == assembly]
    del df['Assembly']
    df = df[df['ClinicalSignificance'].str.lower().str.replace('conflicting interpretations of pathogenicity','').str.contains('pathogenic')]
    df['Type'] = df['Type'].str.replace('single nucleotide variant', 'SNV')
    df.columns = ['Mutation Type', 'Gene', 'Significance', 'ReviewStatus', 'Phenotypes', 'Chromosome', 'Position', 'VariationID', 'Reference', 'Variant']
    df = df[['VariationID', 'Gene', 'Mutation Type', 'Significance', 'ReviewStatus', 'Phenotypes', 'Chromosome', 'Position', 'Reference', 'Variant']]
    df.to_csv('Variants_' + assembly + '.tsv', sep = '\t', index = False)


print('\nStarting analysis')

if ref == '38':
    try:
        cleanerVar = pd.read_csv('Variants_GRCh38.tsv', dtype = str, sep = '\t')
    except:
        print('variants file not found, downloading...')
        GenerateByAssembly('GRCh38')
        print('Finished downloading')
        cleanerVar = pd.read_csv('Variants_GRCh38.tsv', dtype = str, sep = '\t')
if ref == '19':
    try:
        cleanerVar = pd.read_csv('Variants_GRCh37.tsv', dtype = str, sep = '\t')
    except:
        print('variants file not found, downloading...')
        GenerateByAssembly('GRCh37')
        print('Finished downloading')
        cleanerVar = pd.read_csv('Variants_GRCh37.tsv', dtype = str, sep = '\t')
cleanerVar['variants_coordinates'] = cleanerVar['Chromosome'] + ':' + cleanerVar['Position']
variants_coordinates = cleanerVar['variants_coordinates'].tolist()

for n, vcf in enumerate(vcfs):
    print(str(n + 1) + '/' + str(total) + ' - Processing ' + vcf + '...')
    try:
        df = pd.read_csv(vcf, skiprows = skipHeader(vcf), sep = '\t', dtype = str)
        for col in ['ID', 'QUAL', 'FILTER', 'INFO', 'FORMAT']:
            del df[col]
        df['#CHROM'] = df['#CHROM'].astype(str).str.lower()
        df['#CHROM'] = df['#CHROM'].str.replace('chr','').str.upper().str.replace('MT','M')
        chromosomes = [str(v) for v in range(1,23)] + ['X','Y','M']
        df = df[df['#CHROM'].isin(chromosomes)]
        df['POS'] = df['POS'].astype(str)
        df['#CHROM'] = df['#CHROM'] + ':' + df['POS']
        del df['POS']
        df.rename(columns = {'#CHROM' : 'Coordinates - hg' + ref}, inplace = True)
        samples = df.columns[3:]
        for sample in samples:
            df[sample] = df[sample].apply(lambda x : x.split(':')[0].replace('|','/').replace('0/1','1').replace('1/2','1').replace('1/1','2').replace('2/2','2'))
            cols_for_sample = list(df.columns[:3]) + [sample]
            sample_df = df[cols_for_sample]
            sample_df = df[df[sample].isin(['1','2'])]
            sample_df = df[df['Coordinates - hg' + ref].isin(variants_coordinates)]
            sample_df['REF'] = sample_df['REF'].str.replace(',','|')
            sample_df['ALT'] = sample_df['ALT'].str.replace(',','|')
            sample_df.to_csv(temp + '/' + sample + '_single.YYF.gz', index = None, compression='gzip')
            print(vcf + ' processed successfully')
    except:
        open(temp + '/' + 'Errors.txt', 'a').write('Error processing ' + vcf + '\n')
        print('Error processing ' + vcf)

single_YYFs = glob(temp + '/*_single.YYF.gz')
if len(single_YYFs) == 0:
    input('None of the VCF files were processed successfully')
    exit()

print('\nMerging into a single file...')
df = pd.read_csv(single_YYFs[0], dtype = 'str')
for n, YYF in enumerate(single_YYFs[1:]):
    print(str(n+2) + '/' + str(total))
    df = pd.merge(df, pd.read_csv(YYF, dtype = 'str'), on = ['Coordinates - hg' + ref,'REF','ALT'], how='outer')
print('Mapping variants...')

if ref == '38':
    ClinVarDF = pd.read_csv('Variants_GRCh38.tsv', dtype = str, sep = '\t')
    ClinVarDF['Coordinates - hg38'] = ClinVarDF['Chromosome'] + ':' + ClinVarDF['Position']
elif ref == '19':
    ClinVarDF = pd.read_csv('Variants_GRCh37.tsv', dtype = str, sep = '\t')
    ClinVarDF['Coordinates - hg19'] = ClinVarDF['Chromosome'] + ':' + ClinVarDF['Position']
ClinVarDF.rename(columns = {'Reference' : 'REF', 'Variant' : 'ALT'}, inplace = True)

ClinVarDf2 = ClinVarDF.copy()
ClinVarDf2['REF'] = ClinVarDf2['REF'].apply(lambda x : x.replace('A','t').replace('T','a').replace('C','g').replace('G','c').upper())
ClinVarDf2['ALT'] = ClinVarDf2['ALT'].apply(lambda x : x.replace('A','t').replace('T','a').replace('C','g').replace('G','c').upper())
ClinVarDF = pd.concat([ClinVarDF,ClinVarDf2]).reset_index(drop=True)
ClinVarDF.rename(columns={'VariationID' : 'ClinVar ID', 'ReviewStatus' : 'Review Status', 'Significance': 'Severity'}, inplace = True)

df = pd.merge(df, ClinVarDF, how = 'inner', on = ['Coordinates - hg' + ref,'REF','ALT'])
df['Coordinates - hg' + ref] = df['Coordinates - hg' + ref] + ' ' + df['REF'] + '>' + df['ALT']
del df['REF']
del df['ALT']
first_columns = ['Phenotypes','Gene','Coordinates - hg' + ref,'ClinVar ID','Severity','Review Status']
df = df[['Phenotypes','Gene','Coordinates - hg' + ref,'ClinVar ID','Severity','Review Status'] + [col for col in df.columns if col not in first_columns]]
custom_ID = 0
def customizeClinVarID(ClinVarID):
    global custom_ID
    ClinVarID = str(ClinVarID).replace('.0','')
    if ClinVarID.isdigit():
        return ClinVarID
    else:
        custom_ID -= 1
        return custom_ID

df['ClinVar ID'] = df['ClinVar ID'].apply(customizeClinVarID)
df.drop_duplicates(subset = ['ClinVar ID'], keep = 'first', inplace = True)
df.drop_duplicates(subset = list(df.columns[:3]), keep = 'first', inplace = True)
df.drop_duplicates(subset = list(df.columns[1:4]), keep = 'first', inplace = True)

print('Generating the output file...')
df.fillna('', inplace = True)
for col in df.columns:
    df[col] = df[col].astype(str)
    df[col] = df[col].str.replace(',','~~')
del df['Chromosome']
del df['Position']
del df['Mutation Type']

df['Alleles'] = df.apply(lambda x : sum([int(v) for v in x[6:].values if v in ('1','2')]), axis = 1)
df = df[df['Alleles'].isin([0,'None']) == False]
df.sort_values('Alleles', ascending = False, inplace = True)
first_columns = ['Alleles','Phenotypes','Gene','Coordinates - hg' + ref,'ClinVar ID','Severity','Review Status']
df = df[['Alleles', 'Phenotypes','Gene','Coordinates - hg' + ref,'ClinVar ID','Severity','Review Status'] + [col for col in df.columns if col not in first_columns]]
df.to_csv("merged.YYF.gz", index = None, compression='gzip')

groups_file = 'Sample,Group\n' + '\n'.join([sample.replace(',','').strip() + ',' for sample in df.columns[7:]])
if len(glob('groups.csv')) > 0:
    input('\nFile groups.csv is already present,\nMove or delete it and then press enter,\notherwise the program can crush.')
open('groups.csv','w').write(groups_file)
shutil.rmtree(temp)
print('\nDone!')
input('You can optionally tag each sample in the groups.csv file\nor start analysing right away\nGood Luck')
