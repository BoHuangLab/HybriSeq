#!/usr/bin/env python
# coding: utf-8

#HybriSeq analysis functions 

import os
import sys
import numpy as np
import pandas as pd
import csv
import dask.dataframe as dd
from Bio import pairwise2

def reverse_complement(dna):
    complement = {'A': 'T','a': 'T', 'C': 'G','c': 'G', 'G':'C','g': 'C','T': 'A','t': 'A','N':'N','n':'N'}
    return ''.join([complement[base] for base in dna[::-1]])

def count_sam_headers(file_path):
    header_count = 0
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('@'):
                header_count += 1
            else:
                break  # Stop once we reach the alignment data
    return header_count

def getKmer(seq, k):
    if k>len(seq):
         return([])
    # catch case len seq = k
    if len(seq)==k:
        return([seq])
    kmer = []
    for i in range(len(seq)-k):
        kmer.append(seq[i:i+k])
    return(kmer)

def lookUPalg(BC_lookup, 
              BC_lookup_1,
              seq_line,
              start = -9,
              stop = 0,
              k = 7):
    '''
    function to search a sequences for a barcode 
    
    BC_lookup: dict, {'barcode seq': 'well'}, non exact matches, all posible barcodes, barcode seq length = k
    BC_lookup_1: dict, {'barcode seq': 'well'}, exact matches, all posible barcodes, barcode seq length = k
    seq_line: str, sequence to do search of 
    start: int, where to start search 
    stop: int, where to stop search 
    k: int, lenth of barcode sequence in BC_lookup and BC_lookup_1
    
    returns: (exact,close,exactP,closeP) or None, 
            exact: list of barcode ids that match exactly (in BC_lookup_1)
            close: list of barcode ids that are close (in BC_lookup)
            exactP: list of position on read of exact 
            closeP: list of position on read of close
    
    '''
    kmer = getKmer(seq_line[start:stop], k = k)
    exact = []
    exactP = []
    close = []
    closeP = []
    for sub in kmer:
        if BC_lookup_1.get(sub)!= None:
            exact.append(BC_lookup_1.get(sub))
            
            exactP.append(seq_line.find(sub,start,stop))
    if len(exact) != 0:
        return(exact,close,exactP,closeP)
    
    else:
        for sub in kmer:
            if BC_lookup.get(sub)!= None:
                close.append(BC_lookup.get(sub))
                closeP.append(seq_line.find(sub,start,stop))
    if len(close) != 0:
        return(exact,close,exactP,closeP)
    return(None)
  
def read_HybriSeq_reads(FASTQ_R1, # path to FASTQ read 1
                        FASTQ_R2, # path to FASTQ read 2
                        sam_path, # path to sam file, output of bowtie2 alingment of probe region
                        BC_lookupR3, # dict of barcode 3 seq and well id that are close in sequence
                        BC_lookupR3_1, # dict of barcode 3 seq and well id that are exact in sequence
                        BC_lookupR1_R2,# dict of barcode 1 & 2 seq and well id that are close in sequence
                        BC_lookupR1_R2_1,# dict of barcode 1 & 2 seq and well id that are exact in sequence
                        exact_only = True, # only exact match to BC_lookupR3_1 and BC_lookupR1_R2_1
                        MaxReads = None,
                        output_file = ''):
    
    '''
    function to call barcodes for HybriSeq reads
    
    FASTQ_R1: str, path to FASTQ read 1
    FASTQ_R2: str, path to FASTQ read 2
    sam_path: str, path to sam file, output of bowtie2 alingment of probe region
    BC_lookupR3: dict, dict of barcode 3 seq and well id that are close in sequence {barcode_seq : well_ID}
    BC_lookupR3_1: dict, dict of barcode 3 seq and well id that are exact in sequence {barcode_seq : well_ID}
            for exact matching only as in the publication, BC_lookupR3 = BC_lookupR3_1
    
    BC_lookupR1_R2: dict, dict of barcode 1 & 2 seq and well id that are close in sequence {barcode_seq : well_ID}
    BC_lookupR1_R2_1: dict, dict of barcode 1 & 2 seq and well id that are exact in sequence {barcode_seq : well_ID}
            for exact matching only as in the publication, BC_lookupR1_R2 = BC_lookupR1_R2_1
    
    exact_only: bool (True), only exact match to BC_lookupR3_1 and BC_lookupR1_R2_1
    
    MaxReads: int or None (None), number of reads to assess 
    output_file: str, path to output csv 
    
    '''
    
    #Read 1 structure 
    #_CCAGAGCATTCGNNNNNNNCATCGGCGTACGACTATCCACGTGCTTGAGNNNNNNNGTGGCCGATGTTTCGNNNNNNNN___

    if exact_only == True:
        BC_lookupR3 = BC_lookupR3_1
        BC_lookupR1_R2 = BC_lookupR1_R2_1

    if True:
        # read through header of sam file 
        header_len = count_sam_headers(os.path.normpath(sam_path))
        
        # open fastqs and sam
        fastq_fileR1 = open(os.path.normpath(FASTQ_R1), 'r')
        fastq_fileR2 = open(os.path.normpath(FASTQ_R2), 'r')
        sam = open(os.path.normpath(sam_path), 'r')
        
        # create and open output file
        file_out = open(output_file, 'w')
        file_out.write('BC1,BC2,BC3,UMI,call\n')
        
        # read all the comment lines in sam
        for x in range(header_len):
            sam_line = sam.readline()
        
        count = 0
        
        while True:
            lineA = fastq_fileR1.readline().replace('\n','').replace(' ','_') # read 1
            lineB = fastq_fileR2.readline().replace('\n','').replace(' ','_') # read 2
            
            if lineA == '':
                continue
            if lineA == 'end':
                break
            if lineA[0] == '@':
                lineA = lineA[1:lineA.find('_')]
                line = fastq_fileR1.readline().replace('\n','')
                Plus_line1 = fastq_fileR1.readline().replace('\n','')
                quality_line1 = fastq_fileR1.readline().replace('\n','')
                
                lineB = lineB[1:lineB.find('_')]
                line2 = fastq_fileR2.readline().replace('\n','')
                Plus_line2 = fastq_fileR2.readline().replace('\n','')
                quality_line2 = fastq_fileR2.readline().replace('\n','')
                
                sam_line = sam.readline()
                sam_line_ID = sam_line.split('\t')[0]
                call = sam_line.split('\t')[2]
                
            
            # check if all line IDs are the same
            if lineA != lineB or lineA != sam_line_ID or lineB != sam_line_ID:
                print('___________________________________')
                print(lineA)
                print(lineB)
                print(sam_line_ID)
                print(sam_line)
                break
                continue

            count += 1
            
            # catch unalinged reads in sam 
            if call == '*': # from sam file
                continue
            
            if MaxReads == None:
                MaxReads = sys.maxsize
                
            if count >= MaxReads:
                break
                
            seq_line2 = line2 # read 2
            seq_line = reverse_complement(line) # read 1
            
            
            # where to start search for sub sequence 
            R1_offset = 0
            R2_offset = 35
            R3_offset = 57
            
            # these are start location relative to offset
            R1_start = seq_line[R1_offset:25].find('ATTCG') 
            R2_start = seq_line[R2_offset:56].find('TGCTTGAG')
            R3_start = seq_line[R3_offset:].find('GTTTCG')
                    
            R1_seq = None
            if R1_start != -1:
                # absolute start location on read
                R1_s = R1_offset+R1_start+len('ATTCG')
                R1_seq = lookUPalg(BC_lookupR1_R2,
                            BC_lookupR1_R2_1,
                            seq_line = seq_line,
                            start = R1_s,
                            stop = R1_s + 7+2,
                            k = 7)
            
            #R2
            R2_seq = None
            if R2_start != -1:
                # absolute start location on read
                R2_s = R2_offset+R2_start+len('TGCTTGAG')
                R2_seq = lookUPalg(BC_lookupR1_R2,
                            BC_lookupR1_R2_1,
                            seq_line = seq_line,
                            start = R2_s,
                            stop = R2_s + 7+2,
                            k = 7)
                    
            #R3 
            R3_seq = None
            if R3_start != -1:
                # absolute start location on read
                R3_s =  R3_offset + R3_start+ len('GTTTCG')
                R3_seq = lookUPalg(BC_lookupR3,
                                   BC_lookupR3_1,
                                   seq_line = seq_line,
                                   start = R3_s,
                                   stop = R3_s+7,
                                   k = 7)  
            
            #R1
            if R1_seq == None:
                R1_seq = 'none'
            else:
                exact = R1_seq[0]
                colose = R1_seq[1]
                
                #check for ambugious identification 
                if len(exact) == 1:
                    R1_seq = exact[0]
                elif len(colose) == 1:
                    R1_seq = colose[0]
                else:
                    R1_seq = 'none'    
            
            #R2
            if R2_seq == None:
                R2_seq = 'none'
            else:
                exact = R2_seq[0]
                colose = R2_seq[1]
                
                #check for ambugious identification 
                if len(exact) == 1:
                    R2_seq = exact[0]
                elif len(colose) == 1:
                    R2_seq = colose[0]
                else:
                    R2_seq = 'none'    

            
            #R3       
            if R3_seq == None:
                R3_seq = 'none'
            else:
                exact = R3_seq[0]
                colose = R3_seq[1]
                
                #check for ambugious identification 
                if len(exact) == 1:
                    R3_seq = exact[0]
                elif len(colose) == 1:
                    R3_seq = colose[0]
                else:
                    R3_seq = 'none'    

               
            BC1 = R1_seq
            BC2 = R2_seq
            BC3 = R3_seq
        
            UMI = seq_line2[60:68]
            
            # exclude un IDed reads
            if 'none' in [BC1,BC2,BC3]:
                continue
            
            # write to output file
            file_out.write(BC1+','+BC2+','+BC3+','+UMI+','+call+'\n')
                 
        # close files
        fastq_fileR1.close()
        fastq_fileR2.close()  
        sam.close()  
        file_out.close()  

def remove_duplicates_dask(input_file, output_file):
    
    '''
    function to remove duplicates (UMI)
    
    '''
    
    # Read the large dataset into a Dask DataFrame
    ddf = dd.read_csv(input_file)

    # Drop duplicates
    unique_ddf = ddf.drop_duplicates(keep = 'first')
    # Write the result to a new file
    unique_ddf.to_csv(output_file, index=False, single_file=True)
        
def make_cell_count_matrix(csv_path,
                           cell_count_save_path):
    
    '''
    function to take barcode-probe CSV and make cell-count matrix

    csv_path: str, path to csv with columns in the order: barcode1,  barcode2,  barcode3, UMI, probeID 
    cell_count_save_path: str, path to save cell-count matrix csv 
    
    
    '''
    
    all_good = {}
    all_probes = {}
    # read CSV and count probes 
    with open(csv_path, mode='r', newline='') as file:
        csv_reader = csv.reader(file)
        c = 0
        for row in csv_reader:
            if c==0:
                c+=1
                continue
            c+=1
            probe = row[4]
            
            cell = row[0]+row[1]+row[2]
            if all_good.get(cell) == None:
                 all_good[cell]={}
    
            if all_good[cell].get(probe) == None:
                 all_good[cell][probe] = 1
                 all_probes[probe]=1
            else:
                all_good[cell][probe]+=1
                all_probes[probe]+=1
    

    all_probes = list(all_probes.keys())
         
    # create dataframe/ cell cound matrix 
    probe_count_list = {}
    for probe in all_probes:
        probe_count_list[probe] = []
    
    cell_col = []    
    for cell in all_good:
        cell_col.append(cell)
        counts = all_good[cell]
        for probe in all_probes:
            probe_count_list[probe].append(counts.get(probe,0))
        
    probe_count_list['cell'] = cell_col
    cell_count = pd.DataFrame(data = probe_count_list,columns=probe_count_list.keys())
    
    # save as CSV
    cell_count.to_csv(cell_count_save_path,index = False)
    
def LR_splitpipe_get_BC(read1):
    
    '''
    
    adapted barcode calling from LR-splitpipe (Github https://doi.org/10.5281/zenodo.5168059) 
    for HybriSeq data
    
    read1: tuple (str, str), (read_name (str), read1 sequence (str))
    
    returns: (read_name ,[Barcode1,Barcode2,Barcode3]), 
    
            Barcode1: str of well position for barcode 1 A1-H12
            Barcode2: str of well position for barcode 2 A1-H12
            Barcode3: str of well position for barcode 3 A1-H12
            
            returns (read_name, ['none','none','none']) if no barcode is found in any position
            
    '''
    
    name, read1 = read1
    
    BC3_dict = {'GACCAAA': 'A1',
             'CAGTTGC': 'A2',
             'TTCATGG': 'A3',
             'TGTTACC': 'A4',
             'ACATTGC': 'A5',
             'ACCACGA': 'A6',
             'GAGTAGT': 'A7',
             'GCCATAG': 'A8',
             'GCACATG': 'A9',
             'GTAGGGT': 'A10',
             'ATGCCTA': 'A11',
             'CTTTAGG': 'A12',
             'ACCTCAG': 'B1',
             'CAGTGTT': 'B2',
             'AGATTCG': 'B3',
             'GACTCTA': 'B4',
             'ACGCTGA': 'B5',
             'AGGCATT': 'B6',
             'ACCAAAC': 'B7',
             'CTCCTTG': 'B8',
             'CGAGAGA': 'B9',
             'GGTTGTC': 'B10',
             'GAAGCTT': 'B11',
             'CAGCTAT': 'B12',
             'ACGGAAT': 'C1',
             'TCTATGC': 'C2',
             'CCCTAGA': 'C3',
             'CCTTGCA': 'C4',
             'GATGGAG': 'C5',
             'CGACAAT': 'C6',
             'CCCAATT': 'C7',
             'CCCTGAT': 'C8',
             'GTCAATG': 'C9',
             'TACAGAC': 'C10',
             'GGTGATT': 'C11',
             'GGCTGAA': 'C12',
             'TGGTCCA': 'D1',
             'TAGGTGT': 'D2',
             'CTGAGCT': 'D3',
             'TCTCTCG': 'D4',
             'AATGCTC': 'D5',
             'CCATCGT': 'D6',
             'CACTTCT': 'D7',
             'AGCAACA': 'D8',
             'CGGTATG': 'D9',
             'ACATCCA': 'D10',
             'TTCACTC': 'D11',
             'AAAACCC': 'D12',
             'CGTGGTA': 'E1',
             'GAGCTTA': 'E2',
             'AAGACGT': 'E3',
             'AAAGGCA': 'E4',
             'CTACGAG': 'E5',
             'GTTTTGC': 'E6',
             'TGGCGAA': 'E7',
             'TCGCTTC': 'E8',
             'TATCCCC': 'E9',
             'CAAGTCG': 'E10',
             'GTTATCG': 'E11',
             'CATACGG': 'E12',
             'AAACACG': 'F1',
             'TCAGCTG': 'F2',
             'CCACGTT': 'F3',
             'TCTGGAT': 'F4',
             'GTATACC': 'F5',
             'TCTTCGA': 'F6',
             'GCGTTTG': 'F7',
             'ATGGTCA': 'F8',
             'TGAAGGT': 'F9',
             'TCGTTCT': 'F10',
             'CAGGACT': 'F11',
             'AGCGAGT': 'F12',
             'AAGAGTC': 'G1',
             'TAGCCTT': 'G2',
             'GTCGTTC': 'G3',
             'AACACTG': 'G4',
             'CTAAAGC': 'G5',
             'CACCGTA': 'G6',
             'CGTGAAG': 'G7',
             'GTGGATA': 'G8',
             'TGAGGAG': 'G9',
             'TGCCATA': 'G10',
             'AACTAGG': 'G11',
             'AGCTTTC': 'G12',
             'TTTCGAC': 'H1',
             'TATGAGG': 'H2',
             'TAGACAG': 'H3',
             'CTTAGTC': 'H4',
             'TCACCAA': 'H5',
             'CTCATCC': 'H6',
             'ATACTCC': 'H7',
             'GTCTTCA': 'H8',
             'GACGTAT': 'H9',
             'CAATGCC': 'H10',
             'ACGAGCA': 'H11',
             'AACCCAT': 'H12'}
    
    BC1_dict = {'TCTTCTG': 'A1',
                 'ACGAACT': 'A2',
                 'GGTTATG': 'A3',
                 'GACTTGA': 'A4',
                 'CCTAAAG': 'A5',
                 'GAACGAC': 'A6',
                 'GGTATAC': 'A7',
                 'AGTAGGT': 'A8',
                 'GAAAGCT': 'A9',
                 'CTTATGG': 'A10',
                 'TCCGATA': 'A11',
                 'AGTGAAG': 'A12',
                 'CAGTGTT': 'B1',
                 'TAGTAGG': 'B2',
                 'TGGATGT': 'B3',
                 'TGCGAAT': 'B4',
                 'TCTGAGC': 'B5',
                 'CCGATAT': 'B6',
                 'GTCACAG': 'B7',
                 'CTGGAAT': 'B8',
                 'AGCTCAG': 'B9',
                 'GTAGATC': 'B10',
                 'TTCGCCA': 'B11',
                 'ATTGGGA': 'B12',
                 'GTACCTT': 'C1',
                 'ATCTTGC': 'C2',
                 'CGACTTG': 'C3',
                 'ACAGACA': 'C4',
                 'ATAACCC': 'C5',
                 'CATACCG': 'C6',
                 'TACACGA': 'C7',
                 'TAAGCTC': 'C8',
                 'CTCTGAT': 'C9',
                 'TTTGGTC': 'C10',
                 'GCTGACT': 'C11',
                 'AAAGGGT': 'C12',
                 'TATCCAC': 'D1',
                 'CGCATAA': 'D2',
                 'TAGGCAT': 'D3',
                 'TACTGGC': 'D4',
                 'AGGGAGT': 'D5',
                 'CAGATTG': 'D6',
                 'TGCTTAC': 'D7',
                 'TTGGATG': 'D8',
                 'ACTCGCT': 'D9',
                 'TCAAGCA': 'D10',
                 'CCCTTTT': 'D11',
                 'CGGTTAG': 'D12',
                 'TGCCTTT': 'E1',
                 'CCTTCAT': 'E2',
                 'CGCTATA': 'E3',
                 'CAAGGAG': 'E4',
                 'CCGTATG': 'E5',
                 'ATGTGAG': 'E6',
                 'ACATCCG': 'E7',
                 'CGTCGAT': 'E8',
                 'GGGGATA': 'E9',
                 'CCCAACA': 'E10',
                 'TGTTGCT': 'E11',
                 'TAGAGTC': 'E12',
                 'ATAGCTG': 'F1',
                 'TACGGTG': 'F2',
                 'GTTTGGT': 'F3',
                 'AGGGTTG': 'F4',
                 'TGGCTAA': 'F5',
                 'GGATGAG': 'F6',
                 'GGTCGTA': 'F7',
                 'AACACTG': 'F8',
                 'AATCACC': 'F9',
                 'GAATCCC': 'F10',
                 'GGTAACA': 'F11',
                 'CAACAGT': 'F12',
                 'CACTCAA': 'G1',
                 'GGGTTTT': 'G2',
                 'CGAGAGA': 'G3',
                 'GTGCAAG': 'G4',
                 'TTCTCTC': 'G5',
                 'GATTAGC': 'G6',
                 'GTTTCCG': 'G7',
                 'GACAACC': 'G8',
                 'ACAGTTC': 'G9',
                 'GCATGTT': 'G10',
                 'TCTCTCG': 'G11',
                 'TGCTCGT': 'G12',
                 'CTACCGA': 'H1',
                 'TAATGCG': 'H2',
                 'CGTTAGT': 'H3',
                 'GCGATTA': 'H4',
                 'TTACGCT': 'H5',
                 'CCTAGTT': 'H6',
                 'TATCGGT': 'H7',
                 'AGGTTGA': 'H8',
                 'GCCAATT': 'H9',
                 'GAAGCGA': 'H10',
                 'CCATGAA': 'H11',
                 'GACCTAT': 'H12'}
    BC2_dict =BC1_dict
    
    L2 = 'CTCAAGCACGTGGATAGTCGTACGCCGATG' # min score 22 (4 miss matched bases)
    P3 = 'CGAAACATCGGCCAC'# min score 11 (2 miss matched bases)

    L2a = pairwise2.align.localms( read1, L2,1, -1, -1, -1,
                 one_alignment_only=True)
    scoreL2 = L2a
    entryL2 = {}
    prefL2 = 'L2'
    entryL2['{}_score'.format(prefL2)] = scoreL2
    alin_scoreL2 = entryL2['L2_score'][0].score
    
    P3a = pairwise2.align.localms( read1, P3,
                  1, -1, -1, -1,
                 one_alignment_only=True)
    scoreP3 = P3a
    entryP3 = {}
    prefP3 = 'P3'
    entryP3['{}_score'.format(prefP3)] = scoreP3
    alin_scoreP3 = entryP3['P3_score'][0].score
    
    # catch poor alingment
    if alin_scoreL2<22 or alin_scoreP3<11:
        return(name,['none','none','none'])
    
    startL2 = entryL2['L2_score'][0].start
    stopL2 = entryL2['L2_score'][0].end
    
    if (stopL2+7) > len(read1) or (startL2-7) < 0:
        return(name,['none','none','none'])
    
    BC1_seq = read1[stopL2:stopL2+7]
    BC2_seq = read1[startL2-7:startL2]
    
    startP3 = entryP3['P3_score'][0].start
    if (startP3 -7) <0: # catch barcods outside of seq 
        return(name,['none','none','none'])
    
    BC3_seq = read1[startP3 -7:startP3]
    
    BC1 = BC1_dict.get(BC1_seq, 'none')
    BC2 = BC2_dict.get(BC2_seq, 'none')
    BC3 = BC3_dict.get(BC3_seq, 'none')
    if 'none' in [BC1,BC2,BC3]:
        return(name,['none','none','none'])
    
    BC_full = [BC1,BC2,BC3]
    
    return(name, BC_full)

