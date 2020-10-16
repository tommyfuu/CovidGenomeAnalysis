
# -*- coding: utf-8 -*-

"""
Author      : Tom Fu
Date        : 2020 June 29
FileName    : rfModel.py (for Coarfa Lab)
Description : Vectorised implementation
"""

import pickle
import numpy as np
import pandas as pd
from boruta import BorutaPy
from sklearn.model_selection import GridSearchCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectPercentile, chi2
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
from sklearn.metrics import f1_score
from sklearn.metrics import accuracy_score
from sklearn.metrics import roc_auc_score

import preprocessing

################################## global variables ##################################

VARIABLENAMES = ['Depression_Remission', 'Anxiety_Remission', 'DepOutMod', 'AnxOutMod', 'SSRI',
                    'SNRI', 'ANTIDEPRESSANTOTHER', 'MOODSTABILIZER_OR_ANTICONVULSANT', 'ATYPICALANTIPSYCHOTIC', 
                    'BENZODIAZEPINE', 'ALTERNATIVE_OR_COMPLEMENTARY']




################################## model codes ##################################


def evaluate(addressX, addressY, model, modelName, parameter_grid = []):
    # initialise
    numOfRepetitions = 100
    gridBool = False
    # normalise and arrayise
    X, binaryY, categoricalY, numericY, binaryLabelDict, categoricLabelDict, numericLabelDict, allFeatureList = preprocessing.BINUnormalise(addressX, addressY)
    
    
    #### prepare to optimise w gridSearchCV to get the best model
    if parameter_grid != []:
        grid_search = GridSearchCV(estimator = model, param_grid = parameter_grid, cv = 3, n_jobs = -1, verbose = 2)
        gridBool = True
    else: grid_search = model
    
    # split, optimise to find the best parameters 
    numOfYVars = len(binaryY[0]) # will be equal to the length of RFParaList
    for index in range(0, numOfYVars):
        ## get the current Y
        currentY = [eachY[index] for eachY in binaryY]
        ## variableName
        currentVarName = VARIABLENAMES[index]
        ## drop out missing values
        currentRealX, currentRealY = preprocessing.sampleDropper(X, currentY)
        ## check homogeneity
        currentRealYList = currentRealY.tolist()
        zeroRatio = (currentRealYList.count(0))/len(currentRealYList)
        
        counter = 0
        F1List = []
        AccuracyScoreList = []
        AUCList = []
        ### feature selector
        currentRealX = SelectPercentile(chi2, percentile=25).fit_transform(currentRealX, currentRealY)
        #### repeat 100 times to test robustness
        for repeat in range(0, numOfRepetitions):
            ## count
            counter +=1 
            print(counter)
            ## split
            trainX, testX, trainYs, testYs = train_test_split(currentRealX, currentRealY, test_size = 0.3)
            ## fit
            grid_search.fit(trainX, trainYs)
            if gridBool == True:
                currentModel = grid_search.best_estimator_
                currentLegend = grid_search.best_params_
            else: 
                currentModel = grid_search
                currentLegend = "No necessary parameters needed"
            ## predict
            currentPredictions = currentModel.predict(testX)
            ## compare and evaluate
            ### AUC scores
            try:
                currentAUCscore = roc_auc_score(testYs, currentPredictions)
                #### make sure there is no value smaller than 0.5
                # if currentAUCscore < 0.5:
                #     for prediction in currentPredictions:
                #         if prediction == 0: prediction = 1
                #         else: prediction = 0
                currentAUCscore = roc_auc_score(testYs, currentPredictions)
                AUCList.append(currentAUCscore)
                print(currentAUCscore)
            except ValueError:
                currentAUCscore = 0
                print("currentAUCscore1 = 0")
                pass
            ### F1 scores
            currentF1Score = f1_score(testYs, currentPredictions, average="weighted")
            F1List.append(currentF1Score)
            print(currentF1Score)
            ### accuracy scores
            currentAccuracyScore = accuracy_score(testYs, currentPredictions)
            AccuracyScoreList.append(currentAccuracyScore)
            print(currentAccuracyScore)


        # plot histograms
        # stacking plots
        fig, axs = plt.subplots(3)
        
        #fig.suptitle('Performance Consistency - binary '+str(index), fontsize=16)
        ## F1 scores
        axs[0].hist(F1List)
        axs[0].set_title('F1 scores')

        ## accuracy scores
        axs[1].hist(AccuracyScoreList)
        axs[1].set_title('Accuracy scores')

        # AUC scores
        axs[2].hist(AUCList)
        axs[2].set_title('AUC scores')

        ## show
        fig.tight_layout()
        fig.subplots_adjust(bottom=0.2)
        fig.text(.5,0.06,'Performance Consistency - binary'+str(index+1), fontsize=18, ha='center')
        fig.text(.5,0.02,'Parameters: '+str(currentLegend),fontsize=10,ha='center')

        plt.savefig('./outputs' + modelName + '/testNoAUChistogramsSel' + str(index+1) + ".png")
        plt.show(block=False)
        plt.pause(3)
        plt.close()

        # form pds and csvs
        dictData = {'F1 scores': F1List, 'Accuracy Scores': AccuracyScoreList, 'AUC Scores': AUCList}
        performanceDF1 = pd.DataFrame(dictData)
        performanceDF1.to_csv('./outputs' + modelName + '/testNOAUCSelPerformance'+str(index+1)+'-'+str(numOfRepetitions)+'.csv')
        
    return
                
