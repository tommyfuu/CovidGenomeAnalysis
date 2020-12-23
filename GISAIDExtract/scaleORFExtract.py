"""
Author      : Tom Fu, Lucy Paddock
Date        : 2020 Oct 30
Description : scaleORFExtract.py for covid biomakerspace project
"""

import pandas as pd
import numpy as np
from ORFExtracter import *
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2
# from Bio.Alphabet.IUPAC import ambiguous_dna

BigFastaAddress = "sequences_2020-10-29_07-32.fasta"
refFasta = 'REFEPI_ISL_402124.fasta'


def scaleOrfAlignScore(BigFastaAddress, scoreDict={}):
    """Go through all the sequences of the giant fasta file, get the orf for each sequnece if any,
    align those orfs with reference orfs, then for each orf, save them in a specific csv file"""
    # get the reference seq
    refDF = wrapOrfClean(refFasta)
    refDict = {}
    for index, row in refDF.iterrows():
        refDict.update({row['orf']: row['sequence']})
    # starting to deal with the giant seq file
    raw = open(BigFastaAddress, 'r')
    lines = raw.readlines()
    start = 0
    end = 300000
    interval = 50
#     for lineIdx in range(len(lines)):
    for lineIdx in range(start, end, interval):
        # info line: even
        if lineIdx % 2 == 0:
            currentInfo = lines[lineIdx][1:].split('/')
            currentLocation = currentInfo[0]
            currentAscensionNum = currentInfo[1]
            currentYear = currentInfo[2][:-1]
            # go to the next line to get the sequence
            currentSeq = Seq(lines[lineIdx+1][:-1])
            # get orfs for the current case
            currentOrfDF = wrapOrfSeqClean(
                currentSeq, currentAscensionNum, currentLocation, currentYear)
            # get the scoreDict so far
            scoreDict = orfAlignScoreDF(refDict, currentOrfDF, scoreDict)
            print(lineIdx)
    print("scoreDict DONE")
    print(scoreDict)
    # go through the scoreDict and build alignment score dfs
    for orf, orfInfo in scoreDict.items():
        currentDF = pd.DataFrame(
            orfInfo, columns=['AscensionNum', 'AlignmentScore', 'Location'])
        currentCSVName = orf + 'AlignmentScore.csv'
        currentDF.to_csv(currentCSVName)
    return scoreDict


def scaleWholeGenomeAlignScore(BigFastaAddress):
    """Go through all the sequences of the giant fasta file, align each sequence
    for the reference sequence, save them in a specific csv file"""
    # get the reference seq
    refSeq = fastaSeqExtract(refFasta)[3]
    # starting to deal with the giant seq file
    raw = open(BigFastaAddress, 'r')
    lines = raw.readlines()
    start = 0
    end = 300000
    interval = 50
    ascensionNumL = []
    scoreL = []
    locationL = []
#     for lineIdx in range(len(lines)):
    for lineIdx in range(start, end, interval):
        # info line: even
        if lineIdx % 2 == 0:
            print(lineIdx)
            currentInfo = lines[lineIdx][1:].split('/')
            currentLocation = currentInfo[0]
            currentAscensionNum = currentInfo[1]
            # go to the next line to get the sequence
            currentSeq = Seq(lines[lineIdx+1][:-1])
            # align currentSeq
            currentWGAlignScore = pairwise2.align.globalxx(refSeq, currentSeq)
            # save into lists
            ascensionNumL.append(currentAscensionNum)
            scoreL.append(currentWGAlignScore)
            locationL.append(currentLocation)
    # go through the scoreDict and build alignment score dfs
    currentDF = pd.DataFrame(
        {'AscensionNum': ascensionNumL, 'AlignmentScore': scoreL, 'Location': locationL})
    currentCSVName = 'WGAlignmentScore.csv'
    currentDF.to_csv(currentCSVName)
    return


def ZScore(scoreDict):
    '''idea: make scoreDict a dictionary of dictionaries. But besides that, this
    function takes in the scoreDict produced from our GISAID file and adds
    an entry to each sequence info tuple of that sequence's Z score relative
    to the rest of its ORF'''
    # loop through orfs
    for orf in scoreDict:

        orfL = scoreDict[orf]
        scoreL = []

        # loop through each sequence sample of one orf
        for seqT in orfL:
            # extract the alignment score from dict
            scoreL.append(seqT[1])
        mean = np.mean(scoreL)
        std = np.std(scoreL)
        print(mean)
        print(std)

        # update sequnce tuples with Z score
        newOrfL = []
        for seqT in orfL:
            ZScore = (seqT[1]-mean)/std
            newOrfL.append((seqT[0], seqT[1], ZScore, seqT[2]))
        scoreDict[orf] = newOrfL

    return scoreDict


'''
Tests
scoreDict = {'ORF1a': [('IVDC-HB-01', 13230.0, 'Wuhan'), ('IVDC-HB-04', 13230.0, 'Wuhan'), ('IVDC-HB-05', 13230.0, 'Wuhan'), ('IPBCAMS-WH-01', 13227.0, 'Wuhan'), ('WIV04', 13230.0, 'Wuhan')],
'ORF1b': [('IVDC-HB-01', 8106.0, 'Wuhan'), ('IVDC-HB-04', 8106.0, 'Wuhan'), ('IVDC-HB-05', 8104.0, 'Wuhan'), ('IPBCAMS-WH-01', 8106.0, 'Wuhan'), ('WIV04', 8106.0, 'Wuhan')],
'ORFS': [('IVDC-HB-01', 3882.0, 'Wuhan'), ('IVDC-HB-04', 3882.0, 'Wuhan'), ('IVDC-HB-05', 3882.0, 'Wuhan'), ('IPBCAMS-WH-01', 3885.0, 'Wuhan'), ('WIV04', 3882.0, 'Wuhan')],
'ORF3a': [('IVDC-HB-01', 873.0, 'Wuhan'), ('IVDC-HB-04', 873.0, 'Wuhan'), ('IVDC-HB-05', 873.0, 'Wuhan'), ('IPBCAMS-WH-01', 873.0, 'Wuhan'), ('WIV04', 874.0, 'Wuhan')],
'ORFN': [('IVDC-HB-01', 1302.0, 'Wuhan'), ('IVDC-HB-04', 1301.0, 'Wuhan'), ('IVDC-HB-05', 1302.0, 'Wuhan'), ('IPBCAMS-WH-01', 1302.0, 'Wuhan'), ('WIV04', 1302.0, 'Wuhan')]}

ZScore(scoreDict) returns
{'ORFS': [('IVDC-HB-01', 3882.0, -0.4999999999999242, 'Wuhan'), ('IVDC-HB-04', 3882.0, -0.4999999999999242, 'Wuhan'), ('IVDC-HB-05', 3882.0, -0.4999999999999242, 'Wuhan'), ('IPBCAMS-WH-01', 3885.0, 2.000000000000076, 'Wuhan'), ('WIV04', 3882.0, -0.4999999999999242, 'Wuhan')],
'ORF1b': [('IVDC-HB-01', 8106.0, 0.4999999999995453, 'Wuhan'), ('IVDC-HB-04', 8106.0, 0.4999999999995453, 'Wuhan'), ('IVDC-HB-05', 8104.0, -2.0000000000004547, 'Wuhan'), ('IPBCAMS-WH-01', 8106.0, 0.4999999999995453, 'Wuhan'), ('WIV04', 8106.0, 0.4999999999995453, 'Wuhan')],
'ORF3a': [('IVDC-HB-01', 873.0, -0.5000000000001137, 'Wuhan'), ('IVDC-HB-04', 873.0, -0.5000000000001137, 'Wuhan'), ('IVDC-HB-05', 873.0, -0.5000000000001137, 'Wuhan'), ('IPBCAMS-WH-01', 873.0, -0.5000000000001137, 'Wuhan'), ('WIV04', 874.0, 1.9999999999998863, 'Wuhan')],
'ORF1a': [('IVDC-HB-01', 13230.0, 0.5000000000003032, 'Wuhan'), ('IVDC-HB-04', 13230.0, 0.5000000000003032, 'Wuhan'), ('IVDC-HB-05', 13230.0, 0.5000000000003032, 'Wuhan'), ('IPBCAMS-WH-01', 13227.0, -1.999999999999697, 'Wuhan'), ('WIV04', 13230.0, 0.5000000000003032, 'Wuhan')],
'ORFN': [('IVDC-HB-01', 1302.0, 0.5000000000001136, 'Wuhan'), ('IVDC-HB-04', 1301.0, -1.9999999999998859, 'Wuhan'), ('IVDC-HB-05', 1302.0, 0.5000000000001136, 'Wuhan'), ('IPBCAMS-WH-01', 1302.0, 0.5000000000001136, 'Wuhan'), ('WIV04', 1302.0, 0.5000000000001136, 'Wuhan')]}
'''


# def scaleExtract(fastaAddress):
#     """extract location, ascensionNum, date, sequence info from gisaid fasta file"""
#     raw = open(fastaAddress, 'r')
#     lines = raw.readlines()
#     orfDF = pd.DataFrame(
#         columns=['seqID', 'location', 'date', 'orf', 'sequence'])
#     # for lineIndex in range(len(lines)):
#     for lineIndex in range(0, 5):
#         if (lineIndex+1) % 2 == 1:
#             currentInfo = lines[lineIndex]
#             infoL = currentInfo.split("/")
#             location = infoL[0][1:]
#             seqID = infoL[1]
#             date = infoL[2][:-1]
#             rawseq = lines[lineIndex+1]
#             if len(rawseq) % 3 == 1:
#                 rawseq += "NN"
#             elif len(rawseq) % 3 == 2:
#                 rawseq += "N"
#             print(seqID)
#             seq = Seq(rawseq)
#             orf_list = find_orfs_with_trans(seq)
#             seqCounter = 0
#             for start, end, strand, pro in orf_list:
#                 for orf, startEndPairs in ORFDict.items():
#                     if orf == "ORF3a-8" and ((start in range(startEndPairs[0], startEndPairs[1])) or (end in range(startEndPairs[2], startEndPairs[3]))):
#                         orfDF.loc[seqCounter] = [seqID, location,
#                                                  date, orf, str(seq[start:end])]
#                         seqCounter += 1
#                     elif (start in range(startEndPairs[0], startEndPairs[1])) and (end in range(startEndPairs[2], startEndPairs[3])):
#                         orfDF.loc[seqCounter] = [seqID, location,
#                                                  date, orf, str(seq[start:end])]
#                         seqCounter += 1
#     return orfDF
