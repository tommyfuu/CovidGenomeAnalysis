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

fastaAddress = "/Users/chenlianfu/Documents/Github/CovidGenomeAnalysis/GISAIDExtract/sequences_2020-10-29_07-32.fasta"


def scaleExtract(fastaAddress):
    """extract location, ascensionNum, date, sequence info from gisaid fasta file"""
    raw = open(fastaAddress, 'r')
    lines = raw.readlines()
    orfDF = pd.DataFrame(
        columns=['seqID', 'location', 'date', 'orf', 'sequence'])
    # for lineIndex in range(len(lines)):
    for lineIndex in range(0, 5):
        if (lineIndex+1) % 2 == 1:
            currentInfo = lines[lineIndex]
            infoL = currentInfo.split("/")
            location = infoL[0][1:]
            seqID = infoL[1]
            date = infoL[2][:-1]
            rawseq = lines[lineIndex+1]
            if len(rawseq) % 3 == 1:
                rawseq += "NN"
            elif len(rawseq) % 3 == 2:
                rawseq += "N"
            print(seqID)
            seq = Seq(rawseq)
            orf_list = find_orfs_with_trans(seq)
            seqCounter = 0
            for start, end, strand, pro in orf_list:
                for orf, startEndPairs in ORFDict.items():
                    if orf == "ORF3a-8" and ((start in range(startEndPairs[0], startEndPairs[1])) or (end in range(startEndPairs[2], startEndPairs[3]))):
                        orfDF.loc[seqCounter] = [seqID, location,
                                                 date, orf, str(seq[start:end])]
                        seqCounter += 1
                    elif (start in range(startEndPairs[0], startEndPairs[1])) and (end in range(startEndPairs[2], startEndPairs[3])):
                        orfDF.loc[seqCounter] = [seqID, location,
                                                 date, orf, str(seq[start:end])]
                        seqCounter += 1
    return orfDF
