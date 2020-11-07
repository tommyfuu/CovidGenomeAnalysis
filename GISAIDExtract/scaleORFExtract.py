"""
Author      : Tom Fu
Date        : 2020 Oct 30
Description : scaleORFExtract.py for covid biomakerspace project
"""

import pandas as pd
from ORFExtracter import *
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
# from Bio.Alphabet.IUPAC import ambiguous_dna

BigFastaAddress = "/Users/chenlianfu/Documents/Github/CovidGenomeAnalysis/GISAIDExtract/sequences_2020-10-29_07-32.fasta"


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
    # for lineIdx in range(len(lines)):
    for lineIdx in range(0, 10):
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
    return

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
