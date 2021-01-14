'''
'''
import numpy as np
import pandas as pd
from clustering import *

import matplotlib.pyplot as plt

address1 = 'outputs1107/ORF1a0300000AlignmentScore.csv'
address2 = 'outputs1107/ORF1b0300000AlignmentScore.csv'
address3 = 'outputs1107/ORF3a0300000AlignmentScore.csv'
address4 = 'outputs1107/ORFN0300000AlignmentScore.csv'
address5 = 'outputs1107/ORFS0300000AlignmentScore.csv'
addressL = [address1, address2, address3, address4, address5]

distAddress1 = "relativeDistOutputs1107/ORF1a0300000AlignmentScore.csv"
distAddress2 = "relativeDistOutputs1107/ORF1b0300000AlignmentScore.csv"
distAddress3 = "relativeDistOutputs1107/ORF3a0300000AlignmentScore.csv"
distAddress4 = "relativeDistOutputs1107/ORFN0300000AlignmentScore.csv"
distAddress5 = "relativeDistOutputs1107/ORFS0300000AlignmentScore.csv"
distAddressL = [distAddress1, distAddress2, distAddress3, distAddress4, distAddress5]

ORFL = ['ORF1a', 'ORF1b', 'ORF3a', 'ORFN', 'ORFS']
NUMOFORFs = 5

scoreDict = csvToScoreDict(addressL)
relativeDistDict = relativeDistFilesToDict(distAddressL)
relativeDistContDict = relativeDistFilesToContinentDict(distAddressL)
zscoreMatrix, relativeDistL, continentCodeL = scoreDictToZMatrix(scoreDict, relativeDistDict, relativeDistContDict)
avgZScoreL = np.mean(zscoreMatrix, axis = 1)

def plottingCorrelation():
    '''Plots the graph of relative distance vs alignment score for each ORF'''
    for addressIndex in range(len(addressL)):
        currentScoreDict = {}
        address = addressL[addressIndex]
        currentScoreDF = pd.read_csv(address)
        distAddress = distAddressL[addressIndex]
        currentDistDF = pd.read_csv(distAddress)
        
        currentORF = ORFL[addressIndex]

        scoreL = currentScoreDF["AlignmentScore"].tolist()
        ascensionNumL = currentScoreDF["AscensionNum"].tolist()
        distL = currentDistDF["relativeDistance"].tolist()
        distNumL = currentDistDF["AscensionNum"].tolist()

        for i in range(len(distNumL)):
            currentScoreDict[distNumL[i]] = (distL[i],)
        for i in range(len(ascensionNumL)):
            if ascensionNumL[i] in currentScoreDict:
                currentScoreDict[ascensionNumL[i]] += (scoreL[i],)
    
        relativeDistanceL = [x[0] for x in list(currentScoreDict.values())]
        alignmentScoreL = [x[1] for x in list(currentScoreDict.values())]

        outputFileName = "./plots/" + currentORF + ".png"

        plt.plot(relativeDistanceL, alignmentScoreL, "ro")
        plt.suptitle(currentORF)
        plt.xlabel("Relative Distance (km)")
        plt.ylabel("Alignment Score")
        plt.savefig(outputFileName)
        plt.close()

def plotByAvgZScore():
    plt.plot(relativeDistL, avgZScoreL, "ro")
    plt.suptitle("Plot by Average Z Score")
    plt.xlabel("Relative Distance (km)")
    plt.ylabel("Normalized Average Alignment Score")
    plt.savefig("./avgZScore.png")
    plt.close()

def plotWithColors():
    plt.scatter(relativeDistL, avgZScoreL, s = 20, c = continentCodeL)
    # plt.legend()

    plt.suptitle("Plot by Average Z Score")
    plt.xlabel("Relative Distance (km)")
    plt.ylabel("Normalized Average Alignment Score")
    plt.savefig("./plots/colors.png")
    plt.close()

def plotByTimeSeries():
    timeStampL = list(range(1, len(avgZScoreL)+1))

    plt.plot(timeStampL, avgZScoreL, "ro")
    plt.ylim(-6, 1)
    plt.savefig("./Plot By Time Series.png")

def outputCSVFile():
    df = pd.DataFrame()
    df["Continent"] = continentCodeL
    df["relativeDistance"] = relativeDistL
    df["avgZScore"] = avgZScoreL
    df.to_csv("zscores.csv")