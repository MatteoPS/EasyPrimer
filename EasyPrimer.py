#! /usr/bin/env python3

# python3.6 EasyPrimer.py -in Example.fas

# Developed by Matteo Perini 2019
# SkyNet UNIMI
# Pediatric Clinical Research Center
# Romeo ed Enrica Invernizzi
# Universita degli Studi di Milano 

# Please cite our work
# 

from Bio import SeqIO
import os
import csv
import pandas as pd
import numpy as np
import logging
from datetime import date
import time
import subprocess
import argparse

start_time = time.time()


def primerHRM(out_folder, tmp_folder, inputfile, conslimit, prilen, amplen, prim_thr, amp_thr, n_adj, w_adj, alg, jobname, snp):

    cwd = os.getcwd() + '/'
    path_scripts = cwd + "Scripts_and_tools/"

    # the following will create the folders needed if they don't exist already
    if not os.path.exists(cwd + tmp_folder):
        os.makedirs(cwd + tmp_folder)
        os.makedirs(cwd + tmp_folder + '/log')
        os.makedirs(cwd + tmp_folder + '/tmp')
    if not os.path.exists(cwd + out_folder):
        os.makedirs(cwd + out_folder)

    os.chdir(cwd + tmp_folder)

    # ########setting log file###################################################
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)

    # create a file handler
    date1 = "{:%m-%d-%Y}".format(date.today())
    handler = logging.FileHandler('log/MAIN_LOG_EasyPrimer_' + date1 + '.txt')
    handler.setLevel(logging.INFO)

    # create a logging format
    formatter = logging.Formatter('%(asctime)s - %(message)s')
    handler.setFormatter(formatter)

    # add the handlers to the logger
    logger.addHandler(handler)
    ########log file set########################################################

    logger.info(jobname + " " + os.path.basename(inputfile) + ' analysis starts'+"======================================")
    logger.info('Consensus limit set at ' + str(conslimit))
    logger.info('Primer lenght ranges from ' + str(min(prilen)) + ' to ' + str(max(prilen)) + ' nts')
    logger.info('Amplicon lengh ranges from ' + str(min(amplen)) + ' to ' + str(max(amplen)) + ' nts')

    
    # to get just the gene name   
    genename_ext = os.path.basename(inputfile)     #mdh.fas
    genename = os.path.splitext(genename_ext)[0]   #mdh
    
    tmp_path = cwd + tmp_folder + "/tmp/"  



    #++++++++++++++++muscle alignment++++++++++++++++++++++++++++++++++++++
    
    if alg != 'n':
        check_algn = False
        mcomand = path_scripts + 'muscle3.8.31_i86linux64 -in '+ cwd + inputfile + ' -out '+ cwd + tmp_folder +'/tmp/'+ genename_ext + ' -loga '+ cwd + tmp_folder +'/log/muscle_log.txt -quiet'
        logger.info('muscle called with the command: ' + mcomand)
        subprocess.run(mcomand, shell=True)
        #change the inputfile path to perform the analysis on the aligned sequences inyo the tmp folder
        inputpath = cwd + tmp_folder + "/tmp/"
    else:
        check_algn = True
        inputpath = cwd
        
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    logger.info('Anlysis STARTS on ' + genename)


    # apre il file e crea una lista delle sequenze
    allseq = []
    tocheck = []
    tocheck1 = []
    k = 0
    
    with open(inputpath + genename_ext, "rU") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            myseq = record.seq
            allseq.append(myseq)
            tocheck.append(str(myseq))
            tocheck1.append(record)
            k = k+1
    # check to see if there are idententical sequences with different names
    tocheck = set(tocheck)
    if len(tocheck) != len(tocheck1):
        logger.info('There is one or more duplicate sequences in the variants list of ' + genename)
        logger.info('PROGRAM HAS STOPPED ###ERROR######ERROR######ERROR######ERROR######ERROR###')
        quit('###ERROR### PROGRAM HAS STOPPED' + '\n' + '---> There is one or more duplicate sequences in the variants list of '+ genename)
    else:
        pass
    
    #check sequences lenght
    if check_algn == True:
        logger.info('MUSCLE alignment avoided, checking if all sequences have the same lenght....')
        the_len = len(allseq[0])
        if not all(len(l) == the_len for l in allseq):
            logger.info('PROGRAM HAS STOPPED ###ERROR######ERROR######ERROR######ERROR######ERROR###')
            logger.info('The alignement provided contains sequences of different lenghts')
            quit('###ERROR### PROGRAM HAS STOPPED' + '\n' + '---> the alignement provided contains sequences of different lenghts')
        else:
            logger.info('Check passed: all sequences have the same lenght')

    snptab = np.zeros((len(allseq[0]), 15), dtype=object)

    consensus = ''
    conscode = ''
    consdic = {             # dictionary of consensus IUPAC ambiguous characters
            '1000': 'A',
            '0100': 'T',
            '0010': 'G',
            '0001': 'C',
            '1010': 'R',  # A or G
            '0101': 'Y',  # C or T  
            '0011': 'S',  # C or G
            '1100': 'W',  # A or T
            '0110': 'K',  # T or G
            '1001': 'M',  # A or C
            '0111': 'B',  # not A (C or G or T)
            '1110': 'D',  # not C (A or G or T)
            '1101': 'H',  # not G (A or C or T)
            '1011': 'V',  # not T (A or C or G)
            '1111': 'N',  # any base
            '0000': '-',  # gap
            '*': '*'      # base insertion in few genes in the alinment
    }
    logger.info('Building the table with residue frequencies...')
    # handler = string with all the bases for the position k
    for k in range(len(allseq[0])): 
        handler = ''
        # selection of all the bases in position k in the alignment
        for currSeq in allseq:        
            handler += (currSeq[k])
        snptab[k][0] = handler.count("A")
        snptab[k][1] = handler.count("T")
        snptab[k][2] = handler.count("G")
        snptab[k][3] = handler.count("C")
        snptab[k][4] = handler.count("-")
        totalleles = sum(snptab[k][0:4])
        snptab[k][5] = float(snptab[k][0]) / totalleles  # freqA
        snptab[k][6] = float(snptab[k][1]) / totalleles  # freqT
        snptab[k][7] = float(snptab[k][2]) / totalleles  # freqG
        snptab[k][8] = float(snptab[k][3]) / totalleles  # freqC
        snptab[k][9] = float(snptab[k][4]) / totalleles  # freq-
        

        # consensus
        conscode = ''

        # A
        if snptab[k][5] >= conslimit:
            conscode += '1' 
        else:
            conscode += '0'    
        # T
        if snptab[k][6] >= conslimit:
            conscode += '1' 
        else:
            conscode += '0'          
        # G
        if snptab[k][7] >= conslimit:
            conscode += '1' 
        else:
            conscode += '0'          
        # C
        if snptab[k][8] >= conslimit:
            conscode += '1' 
        else:
            conscode += '0'
        
        # gaps with frequencies > than 1 - thr get the '*'
        if snptab[k][9] >= (1-conslimit):
            conscode = '*'

        consensus += consdic.get(conscode)
        snptab[k][14] = consdic.get(conscode)

        if snptab[k][13] != '*':
            # fHRM+ (frequency of HRM detectable SNPs)
            
            snptab[k][10] = np.multiply(snptab[k][5], snptab[k][7]) + \
            np.multiply(snptab[k][5], snptab[k][8]) + \
            np.multiply(snptab[k][6], snptab[k][7]) + \
            np.multiply(snptab[k][6], snptab[k][8])
            
            snptab[k][11] = np.multiply(snptab[k][5], snptab[k][7]) + \
            np.multiply(snptab[k][5], snptab[k][8]) + \
            np.multiply(snptab[k][6], snptab[k][7]) + \
            np.multiply(snptab[k][6], snptab[k][8]) + \
            np.multiply(snptab[k][6], snptab[k][5]) + \
                np.multiply(snptab[k][7], snptab[k][8])


            # Shannon index
            shannon = 0
            for j in range(5, 10):
                if snptab[k][j] != 0:
                    shannon += snptab[k][j] * (np.log(snptab[k][j]))
                
            snptab[k][12] = -(shannon)

            # Simpson index
            simpson = 0
            N = sum(snptab[k][0:5])
            for j in range(0, 5):
                simpson += (snptab[k][j]*(snptab[k][j]-1))/(N*(N-1))
            snptab[k][13] = 1 - simpson
        

    # add headers
    header1 = ['Position', 'A', 'T', 'G', 'C', 'gap', 'fA', 'fT', 'fG',
                'fC', 'fgap', 'fHRM+', 'fSNP', 'Shannon', 'Simpson', 'consensus']
    header2 = np.arange(1, len(allseq[0]) + 1)
    snptab = np.insert(snptab, 0, header2, axis=1)
    snptab = np.insert(snptab, 0, header1, axis=0)
    
    arrayFile = str(tmp_path+ jobname + genename + '_TAB.csv')
    with open(arrayFile, "w") as output:
        writer = csv.writer(output, lineterminator='\n')
        writer.writerows(snptab)

    logger.info('The table of residue frequencies is in: ' + arrayFile)

    # .txt file with consensus seq of the alignment
    consensusFile = str(tmp_path + jobname + genename + '_consensus.txt')
    hand = open(consensusFile, 'w')
    hand.write(consensus + '\n')
    hand.close()
    # .txt file with consensus threshold
    consensusTHR = str(tmp_path+ jobname + genename + '_consensus_threshold.txt')
    hand = open(consensusTHR, 'w')
    hand.write('Consenus_Threshold\t' + str(conslimit) + '\n')
    hand.close()

    logger.info('consensus sequence file: '+ consensusFile)

    # Building the Dataframe of all the possible combinations of primers
    logger.info('Building the Dataframe of all the possible combinations of primers and their scores...')
    scoreset = []
    scoreset.append(['ID', 'p1', 'p2', 'p3', 'p4',
    'Score_primer1', 'Score_primer2', 'Score_MaxPrimer', 'Score_amplicon', 
    'Primer_Lenght', 'Amplicon_Lenght', 
    'Primer1_consensus', 'Primer2_consensus', 'Pri1_amb ', 'Pri2_amb', 'amp_HRMamb'])

    ID = 0
    pri_count = 0
    for pri in prilen:
	# list of weight for mean primer shannon score caluclation
        w_list = []
        # this loops will weight the mean shannon index in the primer according to the distance from the amplicon
        for xx in range(prilen[pri_count]-n_adj):
            w_list.append((1-w_adj)/(prilen[pri_count]-n_adj))
        for xx in range(n_adj):
            w_list.append(w_adj/n_adj)
        pri_count += 1
        for amp in amplen:        
            setlen = pri + amp + pri  # total lenght of the current combination of primer1 + amplicon + primer2
            for k in range(1, len(snptab) - setlen):
                ID += 1                 # 'ID': univocal ascending identifier 
                pri1score = float(0)    # 'Score_primer1' mean shannon index for primer1 region
                pri2score = float(0)    # 'Score_primer2'mean shannon index for primer2 region
                maxpri = float(0)       # 'Score_MaxPrimer': max between 'Score_primer1', 'Score_primer2' (the wrose score between them)
                ampscore = float(0)     # 'Score_amplicon': mean FHRM+ in amplicon region 
                pri1cons = ''           # 'Primer1 consensus': consensus seq for primer 1
                pri2cons = ''           # 'Primer2 consensus': consensus seq for primer 2
                pri1_amb = int          # 'Pri1_amb ': number of consensus IUPAC ambiguous characters in primer1
                pri2_amb = int          # 'Pri2_amb': number of consensus IUPAC ambiguous characters in primer2
                amp_HRMamb = int        # 'ampHRMamb': number of consensus IUPAC ambiguous characters that represent a HRM-detectable substitution
                handamp = ''   # handler
                currset = []   # list of values for current combination

                # key positions:
                # combination scheme:  'p1'__primer1__'p2'__amplicon__'p3'__primer2__'p4'
                p1 = k          # p1: primer1 start position
                p2 = p1 + pri   # p2: primer1 end position
                p3 = p2 + amp   # p3: primer2 start position 
                p4 = p3 + pri   # p4: primer2 end position
                # 'Score_primer1
                w_count = 0
                for j in range(p1, p2):
                    pri1score += (float(snptab[j][13]) * w_list[w_count])  # Shannon index
                    pri1cons = pri1cons + snptab[j][15]
                    w_count += 1
                #pri1score = pri1score/pri
                # 'Score_primer2'
                w_count = 0
                rw_list=list(reversed(w_list))
                rw_count = 0
                for jj in range(p3, p4):
                    pri2score += (float(snptab[j][13]) * rw_list[rw_count])  # Shannon index    
                    pri2cons = pri2cons + snptab[jj][15]
                    rw_count += 1
                #pri2score = pri2score/pri
                # 'Score_amplicon'
                for h in range(p2, p3):
                    ampscore += (float(snptab[h][11]))  # FHRM+ (summation of frequencies of the HRM detectable SNPs inside the amplicon)
                    if snptab[h][14] == ('R' or 'M' or 'K' or 'Y'):
                        handamp = handamp + snptab[j][14]
                ampscore = ampscore/amp
                maxpri = max(pri1score, pri2score)
                pri1_amb = pri1cons.count("R") + pri1cons.count("Y") + \
                    pri1cons.count("S") + pri1cons.count("W") + \
                    pri1cons.count("K") + pri1cons.count("M") + \
                    pri1cons.count("B") + pri1cons.count("D") + \
                    pri1cons.count("H") + pri1cons.count("V") + pri1cons.count("N")
                pri2_amb = pri2cons.count("R") + pri2cons.count("Y") + \
                    pri2cons.count("S") + pri2cons.count("W") + \
                    pri2cons.count("K") + pri2cons.count("M") + \
                    pri2cons.count("B") + pri2cons.count("D") + \
                    pri2cons.count("H") + pri2cons.count("V") + pri2cons.count("N")
                amp_HRMamb = pri2cons.count("R") + pri2cons.count("M") + \
                    pri2cons.count("K") + pri2cons.count("Y")
                currset = [ID, p1, p2, p3, p4,
                        pri1score, pri2score, maxpri, ampscore, 
                        pri, amp,
                        pri1cons, pri2cons,
                        pri1_amb, pri2_amb, amp_HRMamb]

                scoreset.append(currset)

    arrayFile = str(tmp_path+ jobname + genename + '_SCORES.csv')
    with open(arrayFile, "w") as output:
        writer = csv.writer(output, lineterminator='\n')
        writer.writerows(scoreset)

    logger.info('The Dataframe with all the combinations of scores is in: '+ arrayFile)

    if not jobname:
        jobname_flag = 'empty'
    else:
        jobname_flag = jobname

    #++++++++++++++++Plotting in R++++++++++++++++++++++++++++++++++++++
    rcomand = "Rscript " + path_scripts + "Plot_OUTPUT.R "+tmp_folder+" "+ genename +" "+ out_folder +" " + cwd +" "+ snp +" "+ str(prim_thr) +" "+ str(amp_thr) +" "+ jobname_flag+" "+ str(conslimit)

    try:
        subprocess.check_call(rcomand, shell=True)
        logger.info('Rscript called with the command: ' + rcomand)
        logger.info('Selecting and plotting data...')
        logger.info('final PDF graph file written in '+ cwd+out_folder)
    except subprocess.CalledProcessError:
        logger.info('PROGRAM HAS STOPPED ###ERROR######ERROR######ERROR######ERROR######ERROR###')
        logger.info('Failing to write the final PDF graph, problems found calling the command: ' + rcomand)
        quit('###ERROR### PROGRAM HAS STOPPED' + '\n' + '---> Failing to write the final PDF graph')
    

    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    finaltime = (time.time() - start_time)
    runtime = str(round(finaltime, 2)) + ' s'
    developed_by = "MP"
    footer = {'Runtime':runtime,
            'Run on':date1,
            'Version':'1.0_linux',
            'Developed by':developed_by}

    df_footer = pd.DataFrame.from_dict(footer, columns=[''], orient='index')
    logger.info(df_footer)
    logger.info(jobname + " " + os.path.basename(inputfile) + ' analysis ends'+"========================================"+ '\n')
    print('==================DONE!=================\nYour Analysis on '+ genename +' is in '+ out_folder+' folder.'+'\n========================================')
    print(df_footer)

def run(args):
    out_folder = args.out_folder
    tmp_folder = args.tmp_folder
    alg = args.aln
    inputfile= args.inputfile
    conslimit = args.consensus_limit
    minpri = args.minpri
    maxpri = args.maxpri
    minamp = args.minamp
    maxamp = args.maxamp
    prim_thr = args.prim_thr
    amp_thr = args.amp_thr
    n_adj = args.n_adj
    w_adj = args.w_adj 
    jobname = args.jobname
    snp = args.snp

    alg = alg.lower()

    prilen = list(range(minpri, maxpri + 1))   # range of possible primer lenghts
    amplen = list(range(minamp, maxamp + 1))   # range of possible amplicon lenghts

    if not jobname:
        jobname = jobname
    else:
        jobname = jobname + '_'
    ########## MAIN FUNCTION ##########
    primerHRM(out_folder, tmp_folder, inputfile, conslimit, prilen, amplen, prim_thr, amp_thr, n_adj, w_adj, alg, jobname, snp)
    ###################################

def main():    
    parser=argparse.ArgumentParser(description="This tool assists the primer design procedure for typing \n if you use it in your work please cite: ????")

    parser.add_argument("-in",help="fasta file of the gene to be analyzed" ,dest="inputfile" , type=str, required=True)

    parser.add_argument("-out",help="Output Folder. Default = Out" ,dest="out_folder" , type=str, default='Out')
    parser.add_argument("-tmp",help="folder with tmp and log files, created automatically if it doesn't exist. Default = 'tmp'" ,dest="tmp_folder" , type=str, default='tmp')
    parser.add_argument("-aln",help="write 'n' to skip sequence alignment. In this case the sequences must be already aligned and the file in the Current Working Directory" ,dest="aln" , type=str, default='y')
    parser.add_argument("-consthr",help="Consensus sequence limit, positive real number between 0 and 1. Default = 0.05" ,dest="consensus_limit" , type=float, default=0.05)
    parser.add_argument("-minpri",help="Minimum primer lenght. Default = 15" ,dest="minpri" , type=int, default=15)
    parser.add_argument("-maxpri",help="Maxium primer lenght. Default = 25" ,dest="maxpri" , type=int, default=25)
    parser.add_argument("-minamp",help="Minimum amplicon lenght. Default = 70" ,dest="minamp" , type=int, default=70)
    parser.add_argument("-maxamp",help="Maximum amplicon lenght. Default = 100" ,dest="maxamp" , type=int, default=100)
    parser.add_argument("-prithr",help="Variability threshold for the selection of the primers, positive real number between 0 and 1. Default = 0.2" ,dest="prim_thr" , type=float, default=0.2)
    parser.add_argument("-ampthr ",help="Variability threshold for the selection of the amplicon, positive real number between 0 and 1. Default = 0.95" ,dest="amp_thr" , type=float, default=0.95)
    parser.add_argument("-npri",help="Number of primer residues next to the aplicon considered fundamental for primer anneling. default = 5" ,dest="n_adj" , type=int, default=5)
    parser.add_argument("-wnpri",help="Weight for the adjustment in the '-npri' residues, positive real number between 0 and 1. default = 0.7" ,dest="w_adj" , type=float, default=0.7)
    parser.add_argument("-prefix",help="Job name to be added as a prefix in final PDF and in tmp files as well" ,dest="jobname" , type=str, default="")
    parser.add_argument("-snp",help="Either 'HRM' or 'ALL', dafault is 'HRM'" ,dest="snp" , choices=['HRM', 'ALL'], default="HRM")
    
    parser.set_defaults(func=run)
    args=parser.parse_args()
    args.func(args)


if __name__=="__main__":
	main()

 
#python3.6 EasyPrimer.py -in Example.fas
