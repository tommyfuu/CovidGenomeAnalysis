"""
Author      : Tom Fu
Date        : 2020 Oct 15
Description : ORFExtracter.py for covid biomakerspace project
Reference   : http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec380
"""

import pandas as pd
from Bio import SeqIO

ORFDict = {"ORF1a": (100, 400, 13200, 13600),
           "ORF1b": (13200, 13600, 21200, 21700),
           "ORFS": (21200, 21700, 25200, 25500),
           "ORF3a": (25200, 25500, 26000, 26500),
           # "ORF3a-8": (25200, 25500, 28000, 28500),
           "ORFN": (28000, 28500, 29200, 30000)}


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
