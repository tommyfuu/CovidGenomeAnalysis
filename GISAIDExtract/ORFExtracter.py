"""
Author      : Tom Fu
Date        : 2020 Oct 15
Description : ORFExtracter.py for covid biomakerspace project
Reference   : http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec380
ORF Ref     : https://www.ncbi.nlm.nih.gov/labs/gquery/all/?term=covid+19&utm_source=Datasets
ORFFutureRef: https://virologyj.biomedcentral.com/articles/10.1186/s12985-020-01402-1#Sec6
"""

import pandas as pd
from Bio import SeqIO
from Bio import pairwise2

ORFDict = {"ORF1a": (0, 400, 13000, 13600),
           "ORF1b": (13000, 13600, 21000, 21700),
           "ORFS": (21000, 21700, 25000, 25700),
           "ORF3a": (25000, 25700, 25800, 26500),
           # "ORF3a-8": (25200, 25500, 28000, 28500),
           "ORFN": (28000, 28700, 28700, 30000)}
refFasta = "./EPI_ISL_402124.fasta"


def fastaSeqExtract(fastaAddress):
    """extract location, ascensionNum, date, sequence info from gisaid fasta file"""
    record = SeqIO.read(fastaAddress, "fasta")
    infoL = record.description.split('/')
    yearAcNTime = infoL[-1].split('|')
    location = infoL[1]
    ascensionNum = yearAcNTime[1]
    date = yearAcNTime[2]
    return location, ascensionNum, date, record.seq


def find_orfs_with_trans(seq, trans_table=1, min_protein_length=100):
    """extract orf from biopython seq from gisaid fasta"""
    answer = []
    seq_len = len(seq)
    for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
        for frame in range(3):
            trans = str(nuc[frame:].translate(trans_table))
            trans_len = len(trans)
            aa_start = 0
            aa_end = 0
            while aa_start < trans_len:
                aa_end = trans.find("*", aa_start)
                if aa_end == -1:
                    aa_end = trans_len
                if aa_end - aa_start >= min_protein_length:
                    if strand == 1:
                        start = frame + aa_start * 3
                        end = min(seq_len, frame + aa_end * 3 + 3)
                    else:
                        start = seq_len - frame - aa_end * 3 - 3
                        end = seq_len - frame - aa_start * 3
                    answer.append((start, end, strand, trans[aa_start:aa_end]))
                aa_start = aa_end + 1
    answer.sort()
    return answer


def wrapOrfPrint(fastaAddress):
    """print approximate orf info of the current covid sequence and save into a pandas df"""
    location, ascensionNum, date, seq = fastaSeqExtract(fastaAddress)
    currentDF = pd.DataFrame(
        columns=['ascensionNum', 'location', 'date', 'orf', 'sequence'])
    print("Info of the current strain:")
    print("Ascension #: ", ascensionNum)
    print("Location: ", location)
    print("Submission date: ", date)
    orf_list = find_orfs_with_trans(seq)
    seqCounter = 0
    for start, end, strand, pro in orf_list:
        for orf, startEndPairs in ORFDict.items():
            if orf == "ORF3a-8" and ((start in range(startEndPairs[0], startEndPairs[1])) or (end in range(startEndPairs[2], startEndPairs[3]))):
                print("Current ORF", orf)
                print(
                    "%s...%s - length %i, strand %i, %i:%i"
                    % (pro[:30], pro[-3:], len(pro), strand, start, end)
                )
                currentDF.loc[seqCounter] = [ascensionNum, location,
                                             date, orf, str(seq[start:end])]
                seqCounter += 1
            elif (start in range(startEndPairs[0], startEndPairs[1])) and (end in range(startEndPairs[2], startEndPairs[3])):
                print("Current ORF", orf)
                print(
                    "%s...%s - length %i, strand %i, %i:%i"
                    % (pro[:30], pro[-3:], len(pro), strand, start, end)
                )
                currentDF.loc[seqCounter] = [ascensionNum, location,
                                             date, orf, str(seq[start:end])]
                seqCounter += 1

    return currentDF


def wrapOrfClean(fastaAddress):
    """print approximate orf info of the current covid sequence and save into a pandas df"""
    location, ascensionNum, date, seq = fastaSeqExtract(fastaAddress)
    currentDF = pd.DataFrame(
        columns=['ascensionNum', 'location', 'date', 'orf', 'sequence'])
    orf_list = find_orfs_with_trans(seq)
    seqCounter = 0
    for start, end, strand, pro in orf_list:
        for orf, startEndPairs in ORFDict.items():
            if orf == "ORF3a-8" and ((start in range(startEndPairs[0], startEndPairs[1])) or (end in range(startEndPairs[2], startEndPairs[3]))):
                currentDF.loc[seqCounter] = [ascensionNum, location,
                                             date, orf, str(seq[start:end])]
                seqCounter += 1
            elif (start in range(startEndPairs[0], startEndPairs[1])) and (end in range(startEndPairs[2], startEndPairs[3])):
                currentDF.loc[seqCounter] = [ascensionNum, location,
                                             date, orf, str(seq[start:end])]
                seqCounter += 1

    return currentDF


def wrapOrfPrintAll(fastaAddress):
    """print approximate orf info of the current covid sequence and save into a pandas df"""
    location, ascensionNum, date, seq = fastaSeqExtract(fastaAddress)
    currentDF = pd.DataFrame(
        columns=['ascensionNum', 'location', 'date', 'orf', 'sequence'])
    print("Info of the current strain:")
    print("Ascension #: ", ascensionNum)
    print("Location: ", location)
    print("Submission date: ", date)
    orf_list = find_orfs_with_trans(seq)
    seqCounter = 0
    for start, end, strand, pro in orf_list:
        for orf, startEndPairs in ORFDict.items():
            if orf == "ORF3a-8" and ((start in range(startEndPairs[0], startEndPairs[1])) or (end in range(startEndPairs[2], startEndPairs[3]))):
                print("Current ORF", orf)
                print(
                    "%s...%s - length %i, strand %i, %i:%i"
                    % (pro[:30], pro[-3:], len(pro), strand, start, end)
                )
            elif (start in range(startEndPairs[0], startEndPairs[1])) and (end in range(startEndPairs[2], startEndPairs[3])):
                print("Current ORF", orf)
                print(
                    "%s...%s - length %i, strand %i, %i:%i"
                    % (pro[:30], pro[-3:], len(pro), strand, start, end)
                )
            else:
                orf = "-"
            currentDF.loc[seqCounter] = [ascensionNum, location,
                                         date, orf, str(seq[start:end])]
            seqCounter += 1
    return currentDF


def orfAlignScore(fastaAddress, scoreDict={}):
    refDF = wrapOrfClean(refFasta)
    refDict = {}
    for index, row in refDF.iterrows():
        refDict.update({row['orf']: row['sequence']})
    currentDF = wrapOrfClean(fastaAddress)
    for index, row in currentDF.iterrows():
        if row['orf'] in refDict:
            currentOrf = row['orf']
            currentAscensionNum = row['ascensionNum']
            currentLocation = row['location']
            currentSeq = row['sequence']
            currentRefSeq = refDict[currentOrf]
            currentScore = pairwise2.align.globalxx(
                currentSeq, currentRefSeq, score_only=True)
            scoreDict.update(
                {currentOrf: (currentScore, currentLocation, currentAscensionNum)})
    return scoreDict
