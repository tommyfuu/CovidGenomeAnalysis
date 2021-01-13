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

distAddress1 = '/Users/chenlianfu/Documents/Github/CovidGenomeAnalysis/GISAIDExtract/relativeDistAlloutputs1107/ORF1a0300000AlignmentScore.csv'
distAddress2 = '/Users/chenlianfu/Documents/Github/CovidGenomeAnalysis/GISAIDExtract/relativeDistAlloutputs1107/ORF1b0300000AlignmentScore.csv'
distAddress3 = '/Users/chenlianfu/Documents/Github/CovidGenomeAnalysis/GISAIDExtract/relativeDistAlloutputs1107/ORF3a0300000AlignmentScore.csv'
distAddress4 = '/Users/chenlianfu/Documents/Github/CovidGenomeAnalysis/GISAIDExtract/relativeDistAlloutputs1107/ORFN0300000AlignmentScore.csv'
distAddress5 = '/Users/chenlianfu/Documents/Github/CovidGenomeAnalysis/GISAIDExtract/relativeDistAlloutputs1107/ORFS0300000AlignmentScore.csv'
distAddressL = [distAddress1, distAddress2,
                distAddress3, distAddress4, distAddress5]
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


def relativeDistFilesToDict(distAddressL):
    """turn relative distance files into a python dictionary where 
    keys are ascension numbers and values are relative distance to Wuhan.
    Duplicates were removed by default."""
    relativeDistDict = {}
    for addressIndex in range(len(distAddressL)):
        address = distAddressL[addressIndex]
        currentDF = pd.read_csv(address)
        currentRelativeDistDict = dict(
            zip(currentDF['AscensionNum'].tolist(),
                currentDF['relativeDistance'].tolist()))
        relativeDistDict.update(currentRelativeDistDict)
    return relativeDistDict


def relativeDistFilesToContinentDict(distAddressL):
    """turn relative distance files into a python dictionary where 
    keys are ascension numbers and values are relative distance to Wuhan.
    Duplicates were removed by default."""
    relativeDistContDict = {}
    for addressIndex in range(len(distAddressL)):
        address = distAddressL[addressIndex]
        currentDF = pd.read_csv(address)
        currentRelativeDistDict = dict(
            zip(currentDF['AscensionNum'].tolist(),
                currentDF['continentCode'].tolist()))
        relativeDistContDict.update(currentRelativeDistDict)
    return relativeDistContDict


def scoreDictToZMatrix(scoreDict, relativeDistDict, relativeDistContDict):
    '''Converts a dictionary of alignment scores (where keys are
    ascension numbers) to a matrix of Z scores we can input into
    sklearn. Columns are ORfs, rows are ascension numbers. Also returns
    a relativeDistL with the distances corresponding to each row
    (indexes are the same)'''
    scoreMatL = []
    relativeDistL = []
    continentCodeL = []
    CONTINENTCONVERT = {"AS": 0, "NorthA": 1, "SA": 2,
                        "AF": 3, "OC": 4, "EU": 5, "AC": 6}

    # if we have all ORF scores and location for an ascNum
    for ascNum in scoreDict.keys():
        if len(scoreDict[ascNum]) == NUMOFORFs:
            if ascNum in relativeDistDict.keys():
                scoreMatL.append(list(scoreDict[ascNum]))
                relativeDistL.append(relativeDistDict[ascNum])
                print(relativeDistContDict[ascNum])
                continentCodeL.append(
                    CONTINENTCONVERT[relativeDistContDict[ascNum]])
    scoreMat = np.array(scoreMatL)

    # normalize scores in matrix (make them Z scores)
    scaler = StandardScaler()
    scaler.fit(scoreMat)

    return scaler.transform(scoreMat), relativeDistL, continentCodeL

def ZMatrixToCsv():
    scoreDict = csvToScoreDict(addressL)
    relativeDistDict = relativeDistFilesToDict(distAddressL)
    relativeDistContDict = relativeDistFilesToContinentDict(distAddressL)
    ZMatrix = scoreDictToZMatrix(scoreDict, relativeDistDict, relativeDistContDict)

    labelArr = ['Order', 'ZScore1', 'ZScore2', 'ZScore3', 'ZScore4', 'ZScore5', 'Distance', 'Continent']
    arr1 = np.column_stack((list(range(0, 5606)), ZMatrix[0], ZMatrix[1], ZMatrix[2]))
    arr = np.vstack((labelArr, arr1))
    np.savetxt("ZMatrixOutput2.csv", arr, delimiter=',', fmt = '%s')

def PCAAnalysis(ZMatrix):
    '''Runs a PCA analysis on our ZMatrix, plotting the PCA
    results in two dimensions, and coloring the plot according
    to clusters determined  by kMeans'''

    kmeansLabels = kMeans(ZMatrix, 4)
    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(ZMatrix)
    principalComponents = principalComponents.tolist()
    for i in range(len(principalComponents)):
        principalComponents[i].append(kmeansLabels[i])
    print(principalComponents[:2])
    finalDf = pd.DataFrame(data=principalComponents, columns=[
        'principal component 1', 'principal component 2', 'Batch'])

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlabel('PC 1', fontsize=15)
    ax.set_ylabel('PC 2', fontsize=15)
    ax.set_title('Covid Meeting', fontsize=20)
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


def gradientPCA(ZMatrix, relativeDistL):

    # run PCA on ZMatrix
    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(ZMatrix)
    principalComponents = principalComponents.tolist()

    finalPrincipalComponents = []
    finalrelativeDistL = []
    # get rid of outliers
    for pcaIndex in range(len(principalComponents)):
        currentX = principalComponents[pcaIndex][0]
        currentY = principalComponents[pcaIndex][1]
        if currentX < 8 and currentY > -2:
            finalPrincipalComponents.append([currentX, currentY])
            finalrelativeDistL.append(relativeDistL[pcaIndex])

    principalComponents = finalPrincipalComponents
    relativeDistL = finalrelativeDistL

    # normalize distance values to between 0 and 1
    maxDist = max(relativeDistL)
    normRelativeDistL = []
    for dist in relativeDistL:
        normRelativeDistL.append((dist/maxDist, dist/maxDist, dist/maxDist))

    for i in range(len(principalComponents)):
        principalComponents[i].append(normRelativeDistL)
    # print(principalComponents[:2])
    finalDf = pd.DataFrame(data=principalComponents, columns=[
        'principal component 1', 'principal component 2', 'Colors'])

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlabel('PC 1', fontsize=15)
    ax.set_ylabel('PC 2', fontsize=15)
    ax.set_title('Covid Meeting', fontsize=20)
    ax.set_facecolor("red")
    ax.scatter(finalDf['principal component 1'],
               finalDf['principal component 2'], c=normRelativeDistL, s=50)
    # ax.legend(targets)
    ax.text(7.5, 7.5, "Colors/Batches")
    ax.grid()
    fig.show()
    return


def continentPCA(ZMatrix, continentCodeL):
    '''Runs a PCA analysis on our ZMatrix, plotting the PCA
    results in two dimensions, and coloring the plot according
    to clusters determined  by kMeans'''

    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(ZMatrix)
    principalComponents = principalComponents.tolist()

    # get rid of outliers
    finalPrincipalComponents = []
    finalcontinentCodeL = []
    for pcaIndex in range(len(principalComponents)):
        currentX = principalComponents[pcaIndex][0]
        currentY = principalComponents[pcaIndex][1]
        if currentX < 8 and currentY > -2:
            finalPrincipalComponents.append([currentX, currentY])
            finalcontinentCodeL.append(continentCodeL[pcaIndex])

    principalComponents = finalPrincipalComponents
    continentCodeL = finalcontinentCodeL

    for i in range(len(principalComponents)):
        principalComponents[i].append(continentCodeL[i])
    finalDf = pd.DataFrame(data=principalComponents, columns=[
        'principal component 1', 'principal component 2', 'Batch'])

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlabel('PC 1', fontsize=15)
    ax.set_ylabel('PC 2', fontsize=15)
    ax.set_title('Covid Meeting', fontsize=20)
    targets = [0, 1, 2, 3, 4, 5, 6]
    colors = ['r', 'b', 'y', 'green', 'black', 'purple', 'cyan']
    for target, color in zip(targets, colors):
        indicesToKeep = finalDf['Batch'] == target
        ax.scatter(finalDf.loc[indicesToKeep, 'principal component 1'],
                   finalDf.loc[indicesToKeep, 'principal component 2'], c=color, s=50)
    ax.legend(targets)
    ax.text(7.5, 7.5, "Colors/Batches")
    ax.grid()
    fig.show()
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
