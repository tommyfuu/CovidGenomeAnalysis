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
from sklearn.model_selection import ParameterGrid
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

EXCLUDEDLIST = ['DepOutMod', 'AnxOutMod','SNRI', 'ATYPICALANTIPSYCHOTIC']

# dict to find index for the member lists of MODELDICT, used in predict
MODELIDEXDICT = {'Depression_Remission':0, 'Anxiety_Remission':1, 'DepOutMod':2, 'AnxOutMod':3, 
                 'SSRI':4,'SNRI':5,'ANTIDEPRESSANTOTHER': 6, 'MOODSTABILIZER_OR_ANTICONVULSANT':7, 
                 'ATYPICALANTIPSYCHOTIC':8, 'BENZODIAZEPINE':9, 'ALTERNATIVE_OR_COMPLEMENTARY':10}

# TRAINED MODELS
RFFILELIST = ['./RFModels/RFModel1Depression_Remission.pkl', 
                './RFModels/RFModel2Anxiety_Remission.pkl', 
                './RFModels/RFModel3DepOutMod.pkl', 
                './RFModels/RFModel4AnxOutMod.pkl', 
                './RFModels/RFModel5SSRI.pkl', 
                './RFModels/RFModel6SNRI.pkl', 
                './RFModels/RFModel7ANTIDEPRESSANTOTHER.pkl', 
                './RFModels/RFModel8MOODSTABILIZER_OR_ANTICONVULSANT.pkl', 
                './RFModels/RFModel9ATYPICALANTIPSYCHOTIC.pkl', 
                './RFModels/RFModel10BENZODIAZEPINE.pkl', 
                './RFModels/RFModel11ALTERNATIVE_OR_COMPLEMENTARY.pkl']


FEATUREDICT =  {'Depression_Remission': ['PANTO-PWY: phosphopantothenate biosynthesis I|unclassified', 'PWY-3841: folate transformations II|unclassified', 
  'PWY-6122: 5-aminoimidazole ribonucleotide biosynthesis II|unclassified'], 'Anxiety_Remission': ['PANTO-PWY: phosphopantothenate biosynthesis I|unclassified'],
   'DepOutMod': [], 'AnxOutMod': [], 'SSRI': [], 'SNRI': ['PWY-3781: aerobic respiration I (cytochrome c)'], 
   'ANTIDEPRESSANTOTHER': ['PWY-6703: preQ0 biosynthesis|g__Bacteroides.s__Bacteroides_thetaiotaomicron'], 'MOODSTABILIZER_OR_ANTICONVULSANT': [], 
   'ATYPICALANTIPSYCHOTIC': ['PWY-5384: sucrose degradation IV (sucrose phosphorylase)|g__Escherichia.s__Escherichia_coli', 
   'PWY-621: sucrose degradation III (sucrose invertase)|g__Escherichia.s__Escherichia_coli', 'PWY-6305: putrescine biosynthesis IV|unclassified'], 
   'BENZODIAZEPINE': ['GLUCONEO-PWY: gluconeogenesis I'], 'ALTERNATIVE_OR_COMPLEMENTARY': ['PWY-5897: superpathway of menaquinol-11 biosynthesis']}

################################## check model ##################################


def evaluatePipeline(addressX, addressY):
    print("DULU")
    # initialise
    numOfRepetitions = 100

    # normalise and arrayise
    X, binaryY, categoricalY, numericY, binaryLabelDict, categoricLabelDict, numericLabelDict, allFeatureList  = preprocessing.BINUnormalise(addressX, addressY)
    
    # Create the random grid to prep for gridSearchCV
    criterions = ['gini', 'entropy']
    n_estimators = [10, 50, 100, 500, 1000] 
    max_depth = [1, 5, 10, 20, 30]

    parameter_grid = {'criterion': criterions,
                'n_estimators': n_estimators,
                'max_depth': max_depth,
                }
    # paramGrid = ParameterGrid(parameter_grid)
    #### prepare to optimise w gridSearchCV to get the best model
    grid_search = GridSearchCV(estimator = RandomForestClassifier(), param_grid = parameter_grid, cv = 3, n_jobs = -1, verbose = 2)
           
    
    # split, optimise to find the best parameters 
    numOfYVars = len(binaryY[0]) # will be equal to the length of RFParaList
    for index in range(0, numOfYVars):
        ## get the current Y
        currentY = [eachY[index] for eachY in binaryY]
        ## variableName
        currentVarName = VARIABLENAMES[index]
        ## drop out missing values
        currentRealX, currentRealY = preprocessing.sampleDropper(X, currentY)
        print(len(currentRealY))
        
        #### initialise things
        counter = 0
        F1List = []
        AccuracyScoreList = []
        AUCList = []
        #### also optimise by feature selection - boruta
        feat_selector = BorutaPy(estimator = RandomForestClassifier(), max_iter = 50, verbose=2)
        feat_selector.fit(currentRealX, currentRealY)
        X_filtered = feat_selector.transform(currentRealX)
        
        #### repeat 100 times to test robustness
        for repeat in range(0, numOfRepetitions):
            ## count
            counter +=1 
            print(counter)
            if len(X_filtered[0]) == 0:
                print("All features are not relevant to the current label")
                X_filtered = X
            ## split
            trainX, testX, trainYs, testYs = train_test_split(X_filtered, currentRealY, test_size = 0.3)
            ## fit
            grid_search.fit(trainX, trainYs)
            currentModel = grid_search.best_estimator_
            currentPara = grid_search.best_params_
            ## predict
            currentPredictions = currentModel.predict(testX)
            print(testYs)
            print(currentPredictions)
            ## compare and evaluate
            ### AUC scores
            try:
                currentAUCscore = roc_auc_score(testYs, currentPredictions)
                #### make sure there is no value smaller than 0.5
                # if currentAUCscore < 0.5:
                #     for prediction in currentPredictions:
                #         if prediction == 0: prediction = 1
                #         else: prediction = 0
                # currentAUCscore = roc_auc_score(testYs, currentPredictions)
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
        fig.text(.5,0.02,'Parameters: '+str(currentPara),fontsize=10,ha='center')
        # plt.title('Parameters: '+str(currentPara), fontsize=16)
        plt.savefig('./outputsRF/NOAUChistogramsRF' + str(index+1) + ".png")
        plt.show(block=False)
        plt.pause(3)
        plt.close()

        # form pds and csvs
        dictData = {'F1 scores': F1List, 'Accuracy Scores': AccuracyScoreList, 'AUC Scores': AUCList}
        performanceDF1 = pd.DataFrame(dictData)
        performanceDF1.to_csv('./outputsRF/NOAUCrfPerformance'+str(index+1)+'-'+str(numOfRepetitions)+'.csv')
    else: 
        ### heterogeneous enough, we are good
        print("heterogeneous enough, we are good")
    return


################################## model codes ##################################
    

def scaleRFModels(addressX, addressY):
    """
    Intro: get the optimised rf models for all existing binary ys
    Para: addressX, a string address of a import file for X
          addressY, a string address of a import file for Y
    Output: svmModelList, a list of optimised svm models 
    """
    # initialise things
    binaryParaDict = {}
    rfFileList = []
    featureDict = {}
    
    # Create the random grid to prep for gridSearchCV
    criterions = ['gini', 'entropy']
    n_estimators = [10, 50, 100, 500, 1000] 
    max_depth = [1, 5, 10, 20, 30]

    parameter_grid = {'criterion': criterions,
                'n_estimators': n_estimators,
                'max_depth': max_depth,
                }

    # normalise and arrayise
    X, binaryY, categoricalY, numericY, binaryLabelDict, categoricLabelDict, numericLabelDict, allFeatureList = preprocessing.BINUnormalise(addressX, addressY)

    # split, optimise to find the best parameters 
    numOfYVars = len(binaryY[0])
    for index in range(0, numOfYVars):
        currentVarName = VARIABLENAMES[index]
        ## get the current Y
        currentY = [eachY[index] for eachY in binaryY]
        
        ## drop out missing values
        currentRealX, currentRealY = preprocessing.sampleDropper(X, currentY)
        
        #### also optimise by feature selection - boruta
        currentFeatureList = []
        feat_selector = BorutaPy(estimator = RandomForestClassifier(), max_iter = 50, verbose=2)
        feat_selector.fit(currentRealX, currentRealY)
        X_filtered = feat_selector.transform(currentRealX)
        if len(X_filtered[0]) == 0:
            print("All features are not relevant to the current label")
            X_filtered = X
        else:
            featureIndexes = np.where(feat_selector.support_ == True)
            for x in np.nditer(featureIndexes):
                currentFeatureList.append(allFeatureList[x])
        featureDict.update({currentVarName: currentFeatureList})
        ## split
        trainX, testX, trainYs, testYs = train_test_split(X_filtered, currentRealY, test_size = 0.3)
        ## optimise w gridSearchCV to get the best model
        grid_search = GridSearchCV(estimator = RandomForestClassifier(), param_grid = parameter_grid, cv = 3, n_jobs = -1, verbose = 2)
        grid_search.fit(trainX, trainYs)
        best_grid_Model = grid_search.best_estimator_
        best_grid_para = grid_search.best_params_
        
        ## append to binary dictionary
        binaryParaDict.update( {index : best_grid_para} )
        
        currentRFFileName = './RFModels/RFModel' + str(index+1) + currentVarName +'.pkl'
        with open(currentRFFileName, 'wb') as f:
            pickle.dump(best_grid_Model, f)
        rfFileList.append(currentRFFileName)

    return binaryParaDict, binaryLabelDict, rfFileList, featureDict
            
def predict(addressX, featureDict = FEATUREDICT, rfFileList = RFFILELIST, variableNames = 'All'):
    """
    Intro: with the best models, predict ys for certain variables the user defined
    Para: modelDict, output of scaleFindParameters, the optimised models found, BY DEFAULT = MODELDICT, models found on 06252020
          X, a normalised version of feature values for a set of samples in need of prediction
          variableName, a user input string that decides the variable we are predicting here. BY DEFAULT
                          variableNames = "All", meaning predict for all the variables the model was
                          trained for. It can also be a string or a list of strings,
                          representing the variables we want to predict for
    Output: predictionDict, a python dict where the keys are the variableNames and 
                          the values are predictions (a list ordered by the order
                          samples were given in X)
    """
    X, allFeatureList = preprocessing.XGetter(addressX)
    allFeatureDict = {val : idx for idx, val in enumerate(allFeatureList)} # features as keys and indices as values
    def variablePredict(rfFileList, currentFeatureList, X, variableName):
        """ 
        helper function to predict for ONE variable
        """
        # find X of matching features
        currentX = X
        if currentFeatureList != []:
            featureIndices = [allFeatureDict[feature] for feature in currentFeatureList]
            effectiveFeatureData = []
            for featureIndex in featureIndices:
                effectiveFeatureData.append([sample[featureIndex] for sample in X])
            newXsArray = np.array(effectiveFeatureData)
            currentX = newXsArray.T
        # get model
        ## get the index for the modelList based on the variableName
        index = MODELIDEXDICT[variableName]
        currentFile = rfFileList[index]
        with open(currentFile, 'rb') as f:
            currentrfModel = pickle.load(f)
        
        # get predictions 
        predictions = currentrfModel.predict(currentX)
        return predictions
    
    # intialise
    predictionDict = {}
    
    # 1. if variableNames is default/a string "All":
    if variableNames == "All":
        for variable in VARIABLENAMES:
            predictionDict.update({variable: variablePredict(rfFileList, featureDict[variable], X, variable)})
    # 2. if variableNames is an individual string/variable
    elif isinstance(variableNames, str): 
        predictionDict.update({variableNames: variablePredict(rfFileList, featureDict[variable], X, variableNames)})
    # 3. if variableNames is a list of strings/variables
    else:
        for variable in variableNames:
            predictionDict.update({variable: variablePredict(rfFileList, featureDict[variable], X, variable)})
    
    return predictionDict