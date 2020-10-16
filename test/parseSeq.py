"""
Author      : Tom Fu
Date        : 2020 Oct 15
Description : parseSeq.py for covid biomakerspace project
"""

import pandas as pd


def parseSeqFile(address, destination):
    """parse the input tsv file and generate a pandas dataframe

    Args:
        address: str address of the input tsv file
        destination: address of the output file

    Returns:
        df: a pandas df with categorised columns
        a destination file that summarises the info from the current file
    """
    currentDF = pd.DataFrame(
        columns=['gene', 'protein', 'proteinID', 'location', 'sequence'])
    seqCounter = 0
    raw = open(address, 'r')
    lines = raw.readlines()
    for lineIdx in range(len(lines)):
        line = lines[lineIdx]
        if line[0] == ">":
            meetingNext = False
            lineSeqIdx = lineIdx+1
            currentSeq = ''
            currentGeneStartIndex = line.find('[gene=')
            currentProteinStartIndex = line.find('[gprotein=')
            currentExceptionStartIndex = line.find('[exception=')
            currentProteinIDStartIndex = line.find('[protein_id=')
            currentLocationStartIndex = line.find('[location=')
            currentgbKeyStartIndex = line.find('[gbkey==')
            currentGene = line[(currentGeneStartIndex + 1):(currentProteinStartIndex-3)]
            currentProtein = line[(currentProteinStartIndex + 1):(currentExceptionStartIndex-3)]
            currentProteinID = line[(
                currentProteinIDStartIndex + 1):(currentLocationStartIndex-3)]
            currentLocation = line[(
                currentLocationStartIndex + 1):(currentgbKeyStartIndex-3)]
            while meetingNext == False:
                if lineSeqIdx+1 < len(lines) and lines[lineSeqIdx][0] != ">":
                    currentSeq += lines[lineSeqIdx][:-2]
                    lineSeqIdx += 1
                else:
                    meetingNext = True
            currentDF.loc[seqCounter] = [currentGene, currentProtein,
                                         currentProteinID, currentLocation, currentSeq]
            seqCounter += 1
    currentDF.to_csv(destination)
    return currentDF


if __name__ == "__main__":
    address = './20CovidSequences.txt'
    destination = './20CovidSequencesParsed.csv'
    df = parseSeqFile(address, '20CovidSequencesParsed.csv')
