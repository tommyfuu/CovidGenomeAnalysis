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
from sklearn.decomposition import PCA

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

# def PCAZMatrix(ZMatrix):


def PCAAnalysis(ZMatrix):

    kmeansLabels = kMeans(ZMatrix, 4)
    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(ZMatrix)
    principalComponents = principalComponents.tolist()
    for i in range(len(principalComponents)):
        principalComponents[i].append(kmeansLabels[i])
    print(principalComponents[:2])
    finalDf = pd.DataFrame(data=principalComponents, columns=[
        'principal component 1', 'principal component 2', 'Batch'])

    # finalDf = pd.concat([principalDf, df[['Batch']]], axis=1)
    # print(finalDf)

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlabel('PC 1', fontsize=15)
    ax.set_ylabel('PC 2', fontsize=15)
    ax.set_title('Gusdon/2 component PCA', fontsize=20)
    targets = [0, 1, 2, 3]
    colors = ['r', 'b', 'y', 'green']
    for target, color in zip(targets, colors):
        indicesToKeep = finalDf['Batch'] == target
        ax.scatter(finalDf.loc[indicesToKeep, 'principal component 1'],
                   finalDf.loc[indicesToKeep, 'principal component 2'], c=color, s=50)
    ax.legend(targets)
    ax.text(7.5, 7.5, "Colors/Batches")
    ax.grid()
    fig.show()
    return


def kMeans(X, components):
    ''''''
    kmeans = KMeans(n_clusters=components, random_state=10).fit(X)
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
