# MetagenomicsMLClean

Machine learning deliverables for Metagenomics data to provide insights on the association between microbiome pathway activities and psychiatric symptom severity and treatments.

**author**: Tom Fu\
**affiliation**: Coarfa Lab, Baylor College of Medicine (BCM SMART Program 2020)\
**dataSource**: Dominique Thompson, Petrosino Lab, Baylor College of Medicine

## To run this code and check performance

- Clone the repo and enter the directory.

### 1 To check you have the right libraries

This repo uses sklearn and keras for learning, BorutaPy for feature selection, pandas and numpy for data and arithmetics, and matplotlib for graphing. To make sure you have the right versions of the mentioned libraries installed, in the terminal, enter `pip install -r requirements.txt` to install the right packages and check. BorutaPy might need to be manually installed by `pip install Boruta`. Note that you need to make sure the version of python you are using - if you are using python3, then you might want to use `pip3` instead of `pip` in the above command lines.

### 2 The binary prediction/classification algorithm

#### 2.1 What this algorithm can do

(1) Take as input an X file with all the patient samples you have with their WGS pathway feature data. An example for an X file can be seen in `/sourceFiles/WGS_pathway_Menninger.txt`\
(2) The model selected with the optimal performance - the random forest model - will be able to predict for you for all those variables you desire with a dictionary. Steps for this can be seen in part 2.2.\
(3) You can train the rfModel further. Steps for this can be seen in part 2.3.\
(4) Pipeline designs and usage for additional individual models in part 2.4.\
(5) Pipeline designs and usage for mode consensus/voting ensembles in part 2.5.\
(6) Pipeline designs and usage for multi-label random forest model in part 2.6.\
(7) Pipeline designs and usage for DBSCAN + random forest model in part 2.7.\
(8) Performance analysis for binary variables in part 2.8.

#### 2.2 Predict with the current rf model

1. Enter the rfModel's directory `binaryModels/rf_bestModel`. Enter into python mode with `python` or `python3`.
2. Import the model by

```
import rfModel
```

3. Prepare an X file formatted like [this file](https://github.com/tommyfuu/MetagenomicsMLClean/blob/master/sourceFiles/WGS_pathway_Menninger.txt). In the terminal, set the variable `addressX` to the path of this X file. Based on what you want, you can also decide one or some or all labels in the following list that you would like to predict for

```
binary_cols_Y = ['Depression_Remission', 'Anxiety_Remission', 'DepOutMod', 'AnxOutMod', 'SSRI',
                 'SNRI', 'ANTIDEPRESSANTOTHER', 'MOODSTABILIZER_OR_ANTICONVULSANT', 'ATYPICALANTIPSYCHOTIC',
                 'BENZODIAZEPINE', 'ALTERNATIVE_OR_COMPLEMENTARY']
```

4. Once you decide which variables/labels you desire to predict for, check which one of the following scenarios you are doing:

- You would like to predict all possible variables, then directly enter `rfModel.predict(addressX)` to predict.
- You would like to predict one of the variables only, then enter `rfModel.predict(addressX, variableNames = 'variableNameAsAString')`.
- You would like to predict some of the variables, then make them into a list of variables such as `variableList = ['Depression_Remission', 'Anxiety_Remission']`, then enter `rfModel.predict(addressX, variableNames = variableList`.
  The model selected with the optimal performance - the random forest model - will be able to predict for you for all those variables you desire with a dictionary. For example, for the following input:

```
variableList = ['Depression_Remission', 'Anxiety_Remission']
rfModel.predict(addressX, variableNames = variableList
```

the output would be something that looks like this:\\

```
{'Depression_Remission': array([1., 0., 0., 0., 1., 1., 0., 1., 1., 0., 1., 0., 0., 1., 0., 1., 0.,
        1., 0., 0., 0., 0., 0., 1., 0., 0., 1., 1., 1., 0., 0., 1., 0., 0.,
        0., 0., 1., 0., 0., 0., 1., 1., 1., 1., 1., 1., 1., 1., 0., 1., 1.,
        0., 0., 1., 1., 0., 0., 1., 1., 0., 0., 0., 1., 0., 1., 0., 1., 1.,
        1., 0., 0., 0., 0., 0., 1., 0., 1., 1., 1., 1., 1., 1., 1., 1., 0.,
        0., 1., 1., 1., 0., 1., 0., 0., 0., 1., 0., 1., 1., 0., 1.]),
 'Anxiety_Remission': array([1., 0., 1., 1., 1., 1., 0., 1., 0., 0., 1., 0., 0., 1., 0., 1., 0.,
        1., 0., 0., 0., 0., 0., 1., 0., 0., 1., 1., 0., 0., 1., 1., 0., 1.,
        0., 0., 0., 0., 0., 0., 1., 0., 1., 0., 1., 1., 1., 1., 1., 1., 1.,
        0., 1., 1., 1., 0., 0., 1., 1., 1., 0., 0., 1., 1., 1., 0., 1., 1.,
        1., 0., 1., 0., 0., 0., 0., 0., 0., 1., 0., 0., 1., 1., 1., 1., 1.,
        0., 1., 1., 1., 0., 0., 0., 0., 0., 0., 0., 1., 1., 1., 1.])}
```

If you would like to evaluate the random forest model, check out 2.4.1.

#### 2.3 Further train the rfModel

If you would like to train the model further by providing more samples for the same set of binary variables, please check out the format of the source files for [Xs](https://github.com/tommyfuu/MetagenomicsMLClean/blob/master/sourceFiles/WGS_pathway_Menninger.txt) and [ys](https://github.com/tommyfuu/MetagenomicsMLClean/blob/master/sourceFiles/WGS_Metadata_remission_Menninger_DataOnly.xlsx). Combine your dataset with the already existing Xs and ys seamlessly. Then you are ready to go to train the rf model again and evaluate its performance. You can do so by running the following codes

```
import rfModel
binaryParaDict, binaryLabelDict, rfFileList, featureDict = rfModel.scaleRFModels(addressX, addressY) # after this replace the global variables RFFILELIST and FEATUREDICT with the new rfFileList, featureDict
```

And then, open the python file `rfModel.py` in the programming environment you desire. Update the global variables `RFFILELIST` and `FEATUREDICT` with the outputs `rfFileList` and `featureDict` you got from the command line above. Then you can predict as usual.

If you would like to evaluate the random forest model, check out 2.4.1.

#### 2.4 Pipeline designs and usage for additional individual models

A summary of pipeline designs including optimisations and nuances for all individual models (including random forest, support vector machine, k-neighbours, linear discriminant analysis, neural networks, and decision trees) can be seen in the graph below:
![optimisation flowchart](https://github.com/tommyfuu/MetagenomicsMLClean/blob/master/performanceAnalysis/classification.png)

To preprocess the data, we

- vectorised data into to train models for individual variables
- implemented a text conversion function that covers more binary data: we can do both yes/no and 1/0 for training. Manual addition to this text conversion can be done quickly if needed.
- implemented a sample dropper function that removes samples that have NAN values for each specific variables. This way we keep all of the useful samples while removing the NANs.
- implemented a normalising function that normalise the X data into the range of (0,1) for each feature.
- implemented a categorising function that distinguishes labels into binary/categorical/numeric labels. Currently it's implemented according to the given instruction, we can quickly implement a function that distinguishes based on text scanning if needed.

To optimise the performances of individual models, we

- used feature selection functions. We used BorutaPy for random forest according to past studies, and sklearn's selectPercentile with 25% threshold for LDA/svm/CART/Logistic/Kneighbors. We did not use a feature selection for the neural networks model because keras is not compatible with Boruta and sklearn packages.
- implemented an AUC filter. Whenever a model gets AUC score smaller than 0.5 in testing, we reverse the prediction of the model. The filter is currently commented out - just decomment it if you want to use it.

##### 2.4.1 Evaluate individual models (including random forest)

In order to evaluate them correctly, you will first need to enter the directories where the python file of the model you are looking at is located and create the output folders. The dictionary below shows the names of the output folders you should have for your individual models.

```
outputDirName = {'randomForest': 'outputsRF', 'supportVectorMachine': 'outputsSVM', 'cart(DecisionTree)': 'outputscart',
                 'kneighbors': 'outputskn', 'linearDiscriminantAnalysis': 'outputsLDA', 'neuralNetworks': 'outputsNN',
                 'logistic': 'outputsLogistic'}
```

The way to evaluate the models are shown below (we assume that you are already in the respective directory and have the right output directory made, and have your addressX and addressY predefined):

```
# random forest
import rfModel
rfModel.evaluatePipeline(addressX, addressY)

# support vector machine
import svmModel
svmModel.svmEvaluate(addressX, addressY)

# k neighbours
import knModel
knModel.knEvaluate(addressX, addressY)

# LDA
import LDAModel
LDAModel.LDAEvaluate(addressX, addressY)

# Logistic
import LogisticModel
LogisticModel.LogisticEvaluate(addressX, addressY)

# cart (decision tree)
import cartModel
cartModel.cartEvaluate(addressX, addressY)

# neural networks
import NNModel
NNModel.evaluateNewNN(addressX, addressY)
```

The predict functions for the models other than random forest are not implemented, so we cannot predict with these models yet. I did not implement them because they are not high AUC models and so not really worth implementing.

To use these outputs and analyse the performance, please go to part 2.8 for reference.

#### 2.5 Pipeline designs and usage for mode consensus/voting ensembles

I implemented two ensemble models. These two models were implemented in a way that is tailored to cater to each individual binary variables. Only those with median AUC above certain threshold (here 0.55) for the variable we are looking at are selected. Feature selection was also used for all of the models involved. There are two such ensemble models, generated from tie breaking. Since to get a majority vote all the time without a tie, we have to have odd number of models. Therefore, when there are even number of models above the 0.55 AUC threshold, there are two ways to break ties:

- Rise to odd model (`voting.py`): we add the next best individual model with AUC score lower than 0.55 to the list of models used for the current variable
- Reduce to odd model (`votingR.py`): we remove the last best individual model with AUC score higher than 0.55 to the list of models used for the current variable

The implementation details are shown in the graph here. Note that the current voting model does not include neural networks anymore due to compatibility issues between sklearn and keras.

![votingModels](https://github.com/tommyfuu/MetagenomicsMLClean/blob/master/performanceAnalysis/voting.png)

##### 2.5.1 predict with ensemble model

Similar to the random forest model, I also wrote a predict function for both ensemble models. You enter the folder `binaryModels/otherModels/votingModels`.

Note that you might need to change your global variables `OLDADDRESSX` and `OLDADDRESSY` to whereever `WGS_pathway_Menninger.txt` and `WGS_Metadata_remission_Menninger_DataOnly.xlsx` are. They should be in the directory `sourceFiles`

For the rise-to-odd model, you can do

```
import voting
voting.predict(addressX, variableNames)
```

to predict.\
As for the reduce-to odd model, you can do

```
import votingR
votingR.predict(addressX, variableNames)
```

##### 2.5.2 further training for the ensemble models

Similar to the random forest model, you can further train with external X and y datasets if you have any. Generally follow the instructions in 2.3, except that you change the training codes into

```
import voting
votingModelDict = voting.scaleVotingModels(addressX, addressY) # after this replace the global variable VOTINGMODELDICT with the new votingModelDict
```

OR

```
import votingR
votingModelDict = votingR.scaleVotingModels(addressX, addressY) # after this replace the global variable VOTINGRMODELDICT with the new votingModelDict
```

And then, you change the global variables `VOTINGMODELDICT` into the outputs you got from the command lines above.

#### 2.6 Pipeline designs and usage for multi-label random forest model

There are two kinds written with two multi-label methods in sklearn - one versus all AND multioutput. A summary of implementation details for the multi-label random forest model is shown here:

![multi-label](https://github.com/tommyfuu/MetagenomicsMLClean/blob/master/performanceAnalysis/multi-labelRF.png)

##### 2.6.1 Predict with multi-label random forest

I wrote predict functions for both multi-label random forest methods. Enter the directory `binaryModels/otherMOdels/multi_label_rf/`. Assuming that you have your addressX predefined, the way to call them is shown below:

```
# one versus all method
import multiLabelRF1VR
multiLabelRF1VR.predict(addressX)

# multi-output method
import multiLabelRFModel
multiLabelRFModel.predict(addressX)
```

##### 2.6.2 Further train the multi-label random forest models

```
# one versus all method
import multiLabelRF1VR
currentModel, OVAFileName = multiLabelRF1VR.scaleMulti1VRModels(addressX, addressY) # after this, replace the global variable OVAFILENAME with the new multiOutputFileName

# multi-output method
import multiLabelRFModel
currentModel, multiOutputFileName = multiLabelRFModel.scaleMultiModels(addressX, addressY) # after this, replace the global variable MULTIOUTPUTFILENAME with the new multiOutputFileName
```

#### 2.7 Pipeline designs and usage for DBSCAN + random forest model

DBSCAN + random forest model essentially just replaces the feature selection method (BorutaPy) to DBSCAN clustering. The clusters we get are used as the training and testing sets. This method might make sense in some cases because inside a cluster, features might share similarities and so it's easier to train the models. In our case, this is not obvious.

##### 2.7.1 Evaluating the performance of DBSCAN + random forest

Enter the directory `binaryModels/otherModels/DBSCAN_rf/`. Make sure to make an output folder called Then do the following to call the evaluating methods

```
import DBSCAN
DBSCAN.scaleClusterEval(addressX, addressY)
```

Note that I did not write a predict function for this case, since its performance is not ideal.

#### 2.8 To evaluate the performance of different models

I chose random forest over others by this evaluation pipeline. Analysis can be seen [here](https://github.com/tommyfuu/MetagenomicsMLClean/blob/master/performanceAnalysis/binaryClassificationAnalysisFinal.pdf).

Since DBSCAN produce a set of outputs slightly different from the way other models do - DBSCAN produce more than one way of clustering for each variable, each clustering produces more than one sets of evaluation results since one cluster will have its own evaluation, we made a different evaluation file, which can be seen [here](https://github.com/tommyfuu/MetagenomicsMLClean/blob/master/performanceAnalysis/DBSCAN-Analysis.pdf).

### 3. Numeric Models

#### 3.1 What this algorithm can do

(1) Take as input an X file with all the patient samples you have with their WGS pathway feature data. An example for an X file can be seen in `/sourceFiles/WGS_pathway_Menninger.txt`\
(2) Pipeline designs and usage for individual models (random forest, support vector machine k-neighbours, and elastic net) in part 3.2.\
(3) Performance analysis for binary variables in part 3.3.

#### 3.2 Pipeline designs and usage for individual models

A summary of implementation details can be seen in the picture below:

![regression](https://github.com/tommyfuu/MetagenomicsMLClean/blob/master/performanceAnalysis/regression.png)

##### 3.2.1 Evaluating performance

Enter the directory `numericModels`. Make sure to make a directory for the numeric model you are using. The dictionary below shows the names of the output folders you should have for your individual models.

```
outputDirName = {'randomForest': 'outputsRF', 'supportVectorMachine': 'outputssvm',
                 'kneighbors': 'outputskn', 'elastic net": 'outputElastic'}
```

Then, assuming that you have addressX and addressY predefined, to evaluate performance, you can use the following:

```
# random forest
import rfModel
rfModel.evaluatePipeline(addressX, addressY)

# k neighbours
import knModel
knModel.knEvaluate(addressX, addressY)

# elastic net
import elasticNetModel
elasticNetModel.evaluatePipeline(addressX, addressY)

# support vector machine
import svmModel
svmModel.svmEvaluate(addressX, addressY)
```

No predict functions were written for the numeric regression problem since the R^2 values were not decent enough.

##### 3.3 Performance analysis

Analysis can be seen [here](https://github.com/tommyfuu/MetagenomicsMLClean/blob/master/performanceAnalysis/numericAnalysis.pdf).

## Acknowledgement

Dr. Cristian Coarfa’s Lab: Tanmay Gandhi, Matthew Robertson, Sandra Grimm.
Dr. Joe Petrosino’s Lab: Dominique S. Thompson
Suggestions and more: research discussion group (Joshua Sepulveda, Zachary Foulks, James Conde, Phyllis Wei)
SMART Program Coordinators: Dr Laurie Connor and Kerri Mejia
Funding support: Harvey Mudd College Office of Community Engagement
Programming Languages: Python, R, and JS.
Programming Environments: Anaconda Spyder and Visual Studio Code.
Main libraries: scikit-learn, keras, BorutaPy, numpy, and pandas.

Cheers,\
TF\
July 27, 2020
