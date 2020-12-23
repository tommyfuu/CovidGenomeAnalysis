"""
Author      : Lucy Paddock, Tom Fu
Date        : 2020 dec 20
Description : clustering.py for covid biomakerspace project
"""

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.preprocessing import StandardScaler


from sklearn.cluster import KMeans
import matplotlib.pyplot as plt

address1 = 'outputs1107/ORF1a0300000AlignmentScore.csv'
address2 = 'outputs1107/ORF1b0300000AlignmentScore.csv'
address3 = 'outputs1107/ORF3a0300000AlignmentScore.csv'
address4 = 'outputs1107/ORFN0300000AlignmentScore.csv'
address5 = 'outputs1107/ORFS0300000AlignmentScore.csv'
addressL = [address1, address2, address3, address4, address5]

# Planning
# 1. Convert csv files to dicts and do z score normalization
# 2. Convert data into matrices to use sklearn
# 3. Models to be used: PCA, https://scikit-learn.org/stable/modules/clustering.html#k-means

ORFL = ['ORF1a', 'ORF1b', 'ORF3a', 'ORFN', 'ORFS']
NUMOFORFs = 5


def csvToScoreDict(addressL):
    '''Turn the orf csv files into scoreDict format where keys 
    are ascensionNums'''
    scoreDict = {}
    for addressIndex in range(len(addressL)):
        address = addressL[addressIndex]
        currentDF = pd.read_csv(address)
        # for determining orf name
        for orfName in ORFL:
            if orfName in address:
                currentORF = orfName
        scoreL = currentDF['AlignmentScore'].tolist()
        ascensionNumL = currentDF['AscensionNum'].tolist()
        for i in range(len(ascensionNumL)):
            if not ascensionNumL[i] in scoreDict:
                scoreDict[ascensionNumL[i]] = (scoreL[i],)
            else:
                scoreDict[ascensionNumL[i]] += (scoreL[i],)
        # mean = np.mean(scoreL)
        # std = np.std(scoreL)
        # if addressIndex == 0:
        #     ZScoreMatrix = [[] for x in scoreL]
        # print(len(scoreL))
        # for index in range(len(scoreL)):
        #     ZScore = (scoreL[index]-mean)/std
        #     # ascensionNum = ascensionNumL[index]
        #     ZScoreMatrix[index].append(ZScore)
    #X = np.matrix(ZScoreMatrix)
    return scoreDict


def scoreDictToZMatrix(scoreDict):
    '''Converts a dictionary of alignment scores (where keys are 
    ascension numbers) to a matrix of Z scores we can input into
    sklearn. Columns are ORfs, rows are ascension numbers.'''
    scoreMatL = []
    for scoreTup in scoreDict.values():
        if len(scoreTup) == NUMOFORFs:
            scoreMatL.append(list(scoreTup))
    scoreMat = np.array(scoreMatL)
    # return stats.zscore(scoreMat, axis=0)
    scaler = StandardScaler()
    scaler.fit(scoreMat)

    return scaler.transform(scoreMat)


def kMeans(X):
    ''''''
    kmeans = KMeans(n_clusters=5, random_state=10).fit(X)
    return kmeans.labels_


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
