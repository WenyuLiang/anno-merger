import pandas as pd
import glob
import itertools
import os
import sys
import psutil
import time
import subprocess
from multiprocessing import Pool, cpu_count
COLUMNS = ['sampleId','readDepth','alleleFrequency','coverage','chrom','inputPos','REF','ALT','rsId','transcript','nucChange','cNomen','pNomen','cosmicIds','clinVarIds','gene','codingEffect','varLocation','ExAC_ALL','ExAC_AFR','ExAC_AMR','ExAC_EAS','ExAC_FIN','ExAC_NFE','ExAC_OTH','ExAC_SAS','gnomAD_exome_ALL','gnomAD_exome_AFR','gnomAD_exome_AMR','gnomAD_exome_ASJ','gnomAD_exome_EAS','gnomAD_exome_FIN','gnomAD_exome_NFE','gnomAD_exome_OTH','gnomAD_exome_SAS','SIFT_score','Polyphen2_HDIV_score','CADD_phred','CADD_raw','CLINSIG','1000g_AF','1000g_EAS_AF','1000g_AMR_AF','transcriptIds','cosmics','chrom_pos_ref_alt_gene','#Chr_Start_Ref_Alt_Ref.Gene','Consequence','varHGVSc','varHGVSp','EXON','INTRON','DOMAINS','1000g_AFR_AF','1000g_EUR_AF','1000g_SAS_AF','AA_AF','EA_AF','MAX_AF','MAX_AF_POPS','SOMATIC','PHENO','PUBMED','MOTIF_NAME','MOTIF_POS','HIGH_INF_POS','MOTIF_SCORE_CHANGE','CADD_PHRED','CADD_RAW','CANONICAL','CLINSIG_PRIORITY','CLINSIG_FINAL','hasClinicalSynopsis','lossOfFunction','inputPosInt','gnomAD_exome_ALL_Int','gnomAD_exome_AFR_Int','gnomAD_exome_AMR_Int','CDS_position','selected_variant','HGNC_SYMONYMS','HGNC_PRE_SYMBOL','VariantMatching','withdrawn_gene','SIFT','Polyphen2','gnomAD_genome_ALL','gnomAD_genome_AFR','gnomAD_genome_AMR','gnomAD_genome_ASJ','gnomAD_genome_EAS','gnomAD_genome_FIN','gnomAD_genome_NFE','gnomAD_genome_OTH','gnomADe_ALL','gnomADe_AFR','gnomADe_AMR','gnomADe_ASJ','gnomADe_EAS','gnomADe_FIN','gnomADe_NFE','gnomADe_OTH','gnomADe_SAS','Clinvar_VARIANT_ID','masterMind_MMID3','masterMind_MMCNT3','masterMind_GENE','GeneSplicer','IMPACT','STRAND','VARIANT_CLASS','VAR_GENE','VAR_SCORE','QUAL','FILTER','GT','Trimmed_variant','AF_hom','AF_het','pop_AF_hom','pop_AF_het','hg38_position','Clinvar_CLNSIG','CLNACC','CLNSIG_BTG','review_status','last_evaluated','gold_stars','consensus_score','curation','BTG_Concensus','HGMD','InterVar','ACMG_classification','ACMG_priority','inheritance','inheritance_score','is_compounds_hets','index','is_ROH']
COLUMNS_SET = set(COLUMNS)
COLUMNS_TYPE = dict(zip(range(len(COLUMNS)), [str] * len(COLUMNS)))

def print_memory_usage():
    process = psutil.Process(os.getpid())
    print(f'Memory usage: {process.memory_info().rss / 1024 ** 2} MB')  # rss: Resident Set Size

def process_file(datafile):
    # Load data
    df = pd.read_csv(datafile, parse_dates=True, na_values='.', sep='\t', on_bad_lines='skip', quoting=3, usecols= [0,4,5,6,7,15], dtype={0: int, 4: str, 5: str, 6: str, 7: str, 15: str})
    sampleId = df['sampleId'][0]
    # Create unique_id column
    df['unique_id'] = df['chrom'].astype(str) + '_' + df['inputPos'].astype(str) + '_' + df['REF'] + '_' + df['ALT'] + '_' + df['gene']
    # Create a local dictionary to hold unique IDs and corresponding sample sets for this file
    local_dict = {}
    anno_index = {}
    # Add sample_id to the set corresponding to each unique_id in the dictionary
    row_index = 0
    for unique_id in df['unique_id']:
        local_dict[unique_id] = {sampleId}
        anno_index[unique_id] = (str(datafile), row_index)
        row_index += 1
    return local_dict, anno_index

def merge_dicts(dict_list):
    # Function to merge multiple dictionaries
    final_dict = {}
    for single_dict in dict_list:
        for key, value in single_dict.items():
            if key in final_dict:
                final_dict[key].update(value)
            else:
                final_dict[key] = value
    return final_dict

def merge_index(anno_index):
    final_index = {}
    final_list = []
    for single_index in anno_index:
        for key, value in single_index.items():
            if key in final_index: # each unique_id only associated with 1 anno file because the information is the same
                continue
            else:
                final_index[key] = value
                final_list.append((key, value))
    with open('sorted_index.txt', 'w') as f:
        for item in final_list:
            f.write("%s\n" % str(item))
    return final_list

def column_alignment(column):
    # pairwise alignment between two columns
    if len(column) == len(COLUMNS):
        return column
    elif len(column) > len(COLUMNS):
        exit('Error: column length is greater than COLUMNS length')
    else:
        for i in range(len(COLUMNS)):
            if COLUMNS[i] not in column:
                column.insert(i, COLUMNS[i])
    assert COLUMNS == column
    return column

def read_anno(datafile):
    #--------------------------------------------------------------pad header--------------------------------------------------------------#
    data = pd.read_csv(datafile, parse_dates=True, na_values='.', sep='\t', on_bad_lines='skip', quoting=3, dtype=COLUMNS_TYPE, nrows=1)
    file_header = data.columns.tolist()
    padded_header = column_alignment(file_header)
    #print(f'padded_header: {padded_header}')
    #--------------------------------------------------------------pad header--------------------------------------------------------------#
    data = pd.read_csv(datafile, parse_dates=True, na_values='.', sep='\t', on_bad_lines='skip', quoting=3, dtype=COLUMNS_TYPE, names=padded_header)
    return data

def merge_anno(dict, group):
    print_memory_usage()  # Monitor memory before function starts
    group = list(group)
    file = group[0]
    group = group[1]
    data = read_anno(file)
    print(f'Total memory usage of DataFrame in MB: {data.memory_usage(deep=True).sum() / 1024 ** 2}')
    with open(f'temp/{file}_temp.csv', 'w') as f:
        for id, anno in group:
            if anno[1] == 0:
                continue
            row = data.iloc[anno[1]].tolist()
            try:
                row[0] = list(dict[id])
            except:
                continue
            pd.DataFrame([row]).to_csv(f, index=False, sep='\t', header=False, mode='a')
    del data
    print_memory_usage()  # Monitor memory after function ends

# def merge_anno_parallel(dict, anno_index):
#     # Separate annotation index into groups of rows by file
#     grouped_anno_index = itertools.groupby(anno_index, key=lambda x: x[1][0])
#     print('grouped_anno_index done')
#     sample = [(dict, (k, list(group))) for k, group in grouped_anno_index]
#     chunknum = len(sample)//4
#     rest = True if len(sample)%4 != 0 else False
#     for i in range(chunknum):     
#         try:
#             with Pool(processes=cpu_count()) as pool:
#                 pool.starmap(merge_anno, sample[i*4:(i+1)*4])
#         except Exception as e:
#             print(e)
#     if rest:
#         try:
#             with Pool(processes=cpu_count()) as pool:
#                 pool.starmap(merge_anno, sample[chunknum*4:])
#         except Exception as e:
#             print(e)
    
#     print('merge_anno done')
#     # write COLUMNS to a file, separate by tab
#     with open('temp/COLUMNS.txt', 'w') as f:
#         for item in COLUMNS:
#             f.write("%s\t" % item)
#     # merge all temp files with subprocess
#     subprocess.run('cat temp/COLUMNS.txt temp/*.csv > merged_data', shell=True)

    # Concatenate all dataframes
    #merged_data = pd.concat(dfs)

    # merged_data.to_csv('merged_data.csv', index=False, sep='\t')
    # return merged_data
    
# def merge_anno_parallel(dict, anno_index):
#     current_file = anno_index[0][1][0]
#     data = read_anno(current_file)
#     with open(f'merged_data.csv', 'w') as f:
#         for id, anno in anno_index:
#             file, row = anno
#             if row == 0:
#                 continue
#             if file != current_file:
#                 data = read_anno(file)
#                 row = data.iloc[row].tolist()
#                 try:
#                     row[0] = list(dict[id])
#                 except:
#                     continue
#                 pd.DataFrame([row]).to_csv(f, index=False, sep='\t', header=False, mode='a')
#                 current_file = file
            # else:
            #     row = data.iloc[row].tolist()
            #     try:
            #         row[0] = list(dict[id])
            #     except:
            #         continue
            #     pd.DataFrame([row]).to_csv(f, index=False, sep='\t', header=False, mode='a')
#     print('merge_anno done')
#     # write COLUMNS to a file, separate by tab
#     with open('temp/COLUMNS.txt', 'w') as f:
#         for item in COLUMNS:
#             f.write("%s\t" % item)
#     # merge all temp files with subprocess
#     subprocess.run('cat temp/COLUMNS.txt merged_data.csv > merged_data', shell=True)

def merge_anno_parallel(dict, anno_index):
    current_file = anno_index[0][1][0]
    data = read_anno(current_file)
    rows = []  # Initiate an empty list to hold rows
    count = 0
    for id, anno in anno_index:
        file, row = anno
        if row == 0:
            continue
        
        if file != current_file:
            # If rows list has more than 10000 rows, write to file and clear list
            if len(rows) % 50000 == 0 and len(rows) != 0:
                count += 1
                pd.DataFrame(rows, columns=COLUMNS).to_csv(f'temp/merged_data{count}.csv', index=False, sep='\t', header=False, mode='a')
                print(f'{count*50000} rows written to file')
                rows.clear()  # Clear the list
            data = read_anno(file)
            current_file = file
        else:
            if len(rows) % 50000 == 0 and len(rows) != 0:
                count += 1
                pd.DataFrame(rows, columns=COLUMNS).to_csv(f'temp/merged_data{count}.csv', index=False, sep='\t', header=False, mode='a')
                print(f'{count*50000} rows written to file')
                rows.clear()  # Clear the list

        row = data.iloc[row].tolist()
        try:
            row[0] = list(dict[id])
        except:
            continue
        rows.append(row)  # Append the row to the list of rows
    # Write remaining rows to file
    if rows:
        pd.DataFrame(rows, columns=COLUMNS).to_csv('temp/merged_data_remain.csv', index=False, sep='\t', header=False, mode='a')
    print('merge_anno done')
    # write COLUMNS to a file, separate by tab
    with open('temp/COLUMNS.txt', 'w') as f:
        for item in COLUMNS:
            f.write("%s\t" % item)
    # merge all temp files with subprocess
    subprocess.run('cat temp/COLUMNS.txt temp/*.csv > merged_data', shell=True)


if __name__ == '__main__':
    # Get a list of all datafiles
    #datafiles = ['79/20230523211720604581.anno','68/202304060856259201461.anno', '82/202305291685355720.anno', '99/202306081686194590.anno']
    #datafiles = ['202305091339456206492.anno', '202305261348457776023.anno', '202305311557043118420.anno', '202306081523306812195.anno', '20230614130630582438.anno', '202305091615071721756.anno', '202305261508326386455.anno', '202305311557072326180.anno', '20230608152332253135.anno', '202306141306316491415.anno', '202305101011006054528.anno', '202305261523296634444.anno', '20230531155710991252.anno', '202306081523337884526.anno', '202306141328449649269.anno']
    #datafiles = ['202305231439012277422.anno', '202305301100151694198.anno', '202306061635129158014.anno', '202306121718165174341.anno']
    #datafiles = ['202305091339456206492.anno', '202305261348457776023.anno', '202305311557043118420.anno', '202306081523306812195.anno']
    datafiles = glob.glob('*.anno')
    # Start the timer
    start = time.time()
    # Create a pool of processes
    with Pool(processes=cpu_count()) as pool:
        # Process the files in parallel
        result_list = pool.map(process_file, datafiles)
    # Stop the timer
    print('generating hash taken = {} seconds'.format(time.time() - start))

    # Start the timer
    start = time.time()
    # only merge_dict the first element of result_list
    data_dict = merge_dicts([r[0] for r in result_list])
    print('Size of data_dict: ', sys.getsizeof(data_dict)/1024/1024, 'MB')  # Prints size of dict in MB
    # Stop the timer
    print('merge_dicts taken = {} seconds'.format(time.time() - start))

    # Start the timer
    start = time.time()
    # merge anno_index
    anno_index = merge_index([r[1] for r in result_list])
    print('Size of anno_index : ', sys.getsizeof(anno_index)/1024/1024, 'MB')  # Prints size of anno_index  in MB
    # Stop the timer
    print('merge_index taken = {} seconds'.format(time.time() - start))
    # Start the timer
    start = time.time()
    merge_anno_parallel(data_dict, anno_index)
    # Stop the timer
    print('merge_anno taken = {} seconds'.format(time.time() - start))