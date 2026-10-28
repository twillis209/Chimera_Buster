## Chimera_buster
## By: Jessica L Albert
## Last edited : 05/17/24 by JLA

import csv
import math
import time
import edlib
import pandas as pd
from concurrent.futures import ProcessPoolExecutor

##############################################################################

# Take cluster concensus fasta and create list for each entry
## [name, fwd umi, rev umi, size]
def readFastq(filename):
    #Reads FASTA file into [name, seq] pairs; sequences may be wrapped over several lines
    samples = []
    with open(filename) as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                samples.append([line, ""])
            elif samples:
                samples[-1][1] += line
    return [sample for sample in samples if sample[1]]

##############################################################################

#changes a list of lists into a sigle list
def flatten_list(data):
    flat_list = []
    for element in data:
        if type(element) == list:
            for item in element:
                flat_list.append(item)
        else:
            flat_list.append(element)
    return flat_list

##############################################################################

#compares two sequences and determines if the edit distances is within the set mismatch_tolerance
def check_mismatch(ref_seq, test_seq, mismatch_tolerance):
    alignment = edlib.align(ref_seq,test_seq,k=mismatch_tolerance)
    mismatching_num = alignment['editDistance']
    mismatching = True
    if mismatching_num < 0 : #if the edit distance is over the k value in edlib, it returns -1 meaning more mismatches than should be tolerated
        #so the sample should be considered distinct and not chimeric
        mismatching = False
    return mismatching

##############################################################################

#cycle through fwd umis and check for any matches with the allowed mismatch_tolerence
def check_FWD_UMI(fwd_sample_num_list, fwd_umi, mismatch_tolerance):
    fwd_lap_time = time.time()
    fwd_check_count = 0
    for seq in fwd_umi:
        fwd_check_count = fwd_check_count + 1
        fwd_sample_num = fwd_umi.index(seq)
        while fwd_sample_num+1 < len(fwd_umi):
            mismatching_fwd = check_mismatch(seq, fwd_umi[fwd_sample_num+1], mismatch_tolerance)
            # removes any sequences deterimined to be chimeras
            if mismatching_fwd:
                del fwd_sample_num_list[fwd_sample_num+1]
                del fwd_umi[fwd_sample_num+1]
            else:
                fwd_sample_num = fwd_sample_num + 1
        if fwd_check_count % 1000 == 0:
            print(str(round((fwd_check_count/len(fwd_sample_num_list)*100), 2)) + " percent of UMIs processed. (%s seconds)              " % (time.time() - fwd_lap_time), end = '\r' )
            fwd_lap_time = time.time()
    return fwd_sample_num_list
    

##############################################################################

#cycle through rev umis and check for any matches with the allowed mismatch_tolerence
def check_REV_UMI(rev_sample_num_list, rev_umi, mismatch_tolerance):
    rev_lap_time = time.time()
    rev_check_count = 0
    for seq in rev_umi:
        rev_check_count = rev_check_count + 1
        rev_sample_num = rev_umi.index(seq)
        while rev_sample_num+1 < len(rev_umi):
            mismatching_rev = check_mismatch(seq, rev_umi[rev_sample_num+1], mismatch_tolerance)
            # removes any sequences deterimined to be chimeras
            if mismatching_rev:
                del rev_sample_num_list[rev_sample_num+1]
                del rev_umi[rev_sample_num+1]
            else:
                rev_sample_num = rev_sample_num + 1
    return rev_sample_num_list

##############################################################################

def chimera_buster(sample_file, size_file, output_name, mismatch_tolerance, check_clusters_status):
    start_time = time.time()

    print ("Loading and preparing samples...")
    sample_list = readFastq(sample_file) 
    sample_size_list = readFastq(size_file)

    entry_num = 0

    for entry in sample_list:
        sample_list[entry_num] = [entry[0].split(";"), entry[1]]
        sample_list[entry_num] = flatten_list(sample_list[entry_num])
        entry_num = entry_num + 1
        
    entry_num = 0

    for entry in sample_size_list:
        sample_size_list[entry_num] = entry[0].split(";")
        sample_size_list[entry_num] = flatten_list(sample_size_list[entry_num])
        entry_num = entry_num + 1
        
    size_dict = {}

    for sublist in sample_size_list:
        size_dict[sublist[-1][10:]] = sublist[-2][5:]
        
    working_list = []

    for sublist in sample_list:
        sample = []
        #add cluster name
        name = sublist[0][10:]
        sample.append(name)
        #add fwd and rev umi seqs from header
        umi_fwd = sublist[4][12:]
        umi_rev = sublist[5][12:]
        sample.append(str(umi_fwd))
        sample.append(str(umi_rev))
        # add final cluster size
        size = sublist[-3][5:]
        sample.append(int(size))
        # add earlier cluster size
        if name[:-2] not in size_dict:
            raise ValueError("Cluster " + name + " has no cluster " + name[:-2] + " in " + size_file)
        pre_size = size_dict[name[:-2]]
        sample.append(int(pre_size))
        if len(sample) != 5:
            print ("Error")
        working_list.append(sample)    
        
    #set up pandas df and put in decending order of concesus size 1st then prelim consensus size 2nd
    df = pd.DataFrame(working_list, columns = ['Name', 'UMI_fwd', 'UMI_rev', 'Consensus_size', 'prelim_size'])
    df = df.sort_values(by=['Consensus_size', 'prelim_size'], ascending = [False, False], ignore_index=True)
    #df.to_csv(("test_outputs.csv"), index=False)
    

    lap_time = time.time()
    
    ##########################################################################
    
    print ("Beginning filtering...")

    #drop any exact matches in the UMIs and keep only top result (which will have the highest prevelence)
    df_no_f_dups = df.drop_duplicates(subset='UMI_fwd', keep='first')
    df_no_dups = df_no_f_dups.drop_duplicates(subset='UMI_rev', keep='first')
    #df_no_dups.to_csv(("test_no_dups_outputs.csv"), index=False)
    

    
    fwd_sample_num_list = list(df_no_dups['Name'])
    rev_sample_num_list = list(df_no_dups['Name'])
    fwd_umi = list(df_no_dups['UMI_fwd'])
    rev_umi = list(df_no_dups['UMI_rev'])
    nonchimera_list =[]

    print ("Checking for duplicate UMIs...")
    #check for any non-exact matches using mismatch_tolerance as the cut off
    if mismatch_tolerance != 0:

        #the fwd and rev checks run in parallel; a failure in either is raised here
        with ProcessPoolExecutor(max_workers=2) as pool:
            fwd_job = pool.submit(check_FWD_UMI, fwd_sample_num_list, fwd_umi, mismatch_tolerance)
            rev_job = pool.submit(check_REV_UMI, rev_sample_num_list, rev_umi, mismatch_tolerance)
            fwd_sample_list = fwd_job.result()
            rev_sample_list = rev_job.result()

        for sample in fwd_sample_list:
            if sample in rev_sample_list:
                nonchimera_list.append(sample)
                rev_sample_list.remove(sample)

    #allows for skipping the lengthy process above when not needed
    else:
        nonchimera_list = fwd_sample_num_list

    #a removed cluster whose fwd AND rev UMIs both match a kept cluster is the same
    #molecule split by the clustering, not a chimera: report it as a duplicate
    duplicate_list = []
    if check_clusters_status == "True":
        df_non_chimeras_prelim = df.loc[df['Name'].isin(nonchimera_list)]
        df_chimeras_prelim = df.loc[~df['Name'].isin(nonchimera_list)]
        print ("Checking for missorted clusters...")

        #check for clustering issues
        cluster_lap_time = time.time()
        num1 = 0
        while num1 < len(df_non_chimeras_prelim["Name"]):
            num2 = 0
            while num2 < len(df_chimeras_prelim["Name"]):
                fwd_chimera = check_mismatch(str(df_non_chimeras_prelim.iloc[num1,1]), str(df_chimeras_prelim.iloc[num2,1]), mismatch_tolerance)
                rev_chimera = check_mismatch(str(df_non_chimeras_prelim.iloc[num1,2]), str(df_chimeras_prelim.iloc[num2,2]), mismatch_tolerance)
                if fwd_chimera == True and rev_chimera == True:
                        duplicate_list.append(df_chimeras_prelim.iloc[num2,0])
                num2 = num2 + 1
            num1 = num1 + 1
            if num1 % 100 == 0:
                print(str(round((num1/len(df_non_chimeras_prelim["Name"])*100), 2)) + " percent of clusters checked. (%s seconds)              " % (time.time() - cluster_lap_time), end = '\r' )
                cluster_lap_time = time.time()

    #make df of final chimeras and non-chimeras
    df_non_chimeras = df.loc[df['Name'].isin(nonchimera_list)]
    df_duplicates = df.loc[df['Name'].isin(duplicate_list) & ~df['Name'].isin(nonchimera_list)]
    df_chimeras = df.loc[~df['Name'].isin(nonchimera_list) & ~df['Name'].isin(duplicate_list)]
    
    print("All samples processed. (%s seconds)                                              " % (time.time() - lap_time))
    lap_time = time.time()
    print("Writing outputs...")

    #outputs csv of chimeras and non-chimeras
    df_non_chimeras.to_csv((output_name + "_non_chimeras.csv"), index=False)
    df_chimeras.to_csv((output_name + "_chimeras.csv"), index=False)
    df_duplicates.to_csv((output_name + "_duplicates.csv"), index=False)

    #writes text file of chimeric transcript IDs        
    with open(output_name+"_chimera_list.txt",'w') as tfile:
        for sample in list(df_chimeras['Name']):
            tfile.write(sample + "\n")

    #writes text file of duplicate transcript IDs (empty unless check_clusters is True)
    with open(output_name+"_duplicate_list.txt",'w') as tfile:
        for sample in list(df_duplicates['Name']):
            tfile.write(sample + "\n")

    #writes text file of non-chimeric transcript IDs         
    with open(output_name+"_nonchimera_list.txt",'w') as tfile:
        for sample in list(df_non_chimeras['Name']):
            tfile.write(sample + "\n")
            
    print("# of Chimeras: " + str(len(list(df_chimeras['Name']))))
    print("# of Non-Chimeras: " + str(len(list(df_non_chimeras['Name']))))
    print("# of Duplicates: " + str(len(list(df_duplicates['Name']))))
    print ("Complete.\n--- %s minutes total ---" % ((time.time() - start_time)/60))
    
