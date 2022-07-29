#!/usr/bin/python

# manojmw
# 29 June, 2022

import argparse, sys, os
import logging
import pandas as pd
from imblearn.over_sampling import SMOTE
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix, classification_report
from joblib import dump


###########################################################

# Parses the canonical cohort file (produced by the grexome-TIMC-secondary 
# pipeline) for which the model is to be built
# Also parses the causalVariant_dict dictionary returned by CausalVariantParser
# function
#
# Preprocessing of the data, selection of features, adjusting 
# for imbalanced classes, training the model.
#
# The trained model will be saved to a file (RandomForest_Model.joblib) 
# in the current working directory. 
def BuildModel(args):

    logging.info("starting to run")

    # Secondary analysis cohort file
    CohortFile = args.incohort
    
    # Getting the exact file name excluding path
    try:
        CohortFileName_fields = CohortFile.split('/')
        Cohort_FName = CohortFileName_fields[-1]
    except:
        Cohort_FName = CohortFile

    # Getting the pathology name from the file name Ex: MMAF.grexT1297.HG35FT.GATK
    FName_fields = Cohort_FName.split('.')
    
    # "cohort" as global variable so that it can be accessed later and also in other 
    # functions
    global cohort
    cohort = FName_fields[0]

    # calling the function
    causalVariant_dict = CausalVariantParser(cohort, args.insample, args.indir)

    # columns of interest
    Required_cols = ['POSITION', 'REF', 'ALT', 'SYMBOL', 'INTERACTORS_PVALUE', 'COUNT_HR', 'COUNT_'+cohort+'_HV', 'COUNT_'+cohort+'_HET', 'COUNT_'+cohort+'_OTHERCAUSE_HV', 'COUNT_'+cohort+'_OTHERCAUSE_HET', 'COUNT_COMPAT_HV', 'COUNT_COMPAT_HET', 'COUNT_NEGCTRL_HV', 'COUNT_NEGCTRL_HET', 'COUNT_OTHERGENO', 'IMPACT', 'gnomADe_AF', 'gnomADg_AF', 'GTEX_testis_RATIO', 'GTEX_ovary_RATIO', 'GTEX_testis', 'GTEX_ovary', 'GTEX_blood', 'GTEX_cerebellar_hemisphere', 'GTEX_liver', 'GTEX_lung', 'CADD_PHRED']

    # Reading the Cohort file
    Features = pd.read_csv(CohortFile, usecols = Required_cols, sep='\t', low_memory = False)

    # Combining the variant info columns into a single column
    # and inserting at the first position
    Features.insert(0, 'VARIANT_ID', Features['POSITION'] + '_' + Features['REF'] + '_' +  Features['ALT'])

    Features['POTENTIALLY_CAUSAL'] = " "

    # Marking causal variants
    for sampleID in causalVariant_dict:
        Features.loc[(causalVariant_dict[sampleID][0] == Features['VARIANT_ID']) & (causalVariant_dict[sampleID][1] == Features['SYMBOL']) & (causalVariant_dict[sampleID][2] == Features['IMPACT']), 'POTENTIALLY_CAUSAL'] = 'YES'
            
    Features["POTENTIALLY_CAUSAL"] = Features["POTENTIALLY_CAUSAL"].fillna('NO')

    # Removing the variant info columns ('POSITION', 'REF', 'ALT') because the data
    # for these 3 columns is stored in the new "VARIANT_ID" column
    # Also dropping 'SYMBOL' and 'IMPACT' columns: causal variants have been marked
    Features.drop(Features.columns[[1,2,3,4,16]], axis='columns', inplace = True)

    # If gnomAD values missing: we assign the value as 0
    Features["gnomADe_AF"] = Features["gnomADe_AF"].fillna(0)
    Features["gnomADg_AF"] = Features["gnomADg_AF"].fillna(0)

    # If INTERACTORS_PVALUE is missing: we assign p-value as 1
    Features["INTERACTORS_PVALUE"] = Features["INTERACTORS_PVALUE"].fillna(1)

    # If the variant is potentially causal in the patient, 
    # but there is no pathogenecity prediction score, then we assign "CADD_PHRED" score as 20
    Features["CADD_PHRED"] = Features["CADD_PHRED"].fillna(20)

    # Drop rows with empty values
    Features = Features.dropna()

    # # print the dataset
    # print(Features.to_csv(sys.stdout, index=False))

    # Feature and target objects
    y = Features['POTENTIALLY_CAUSAL']
    X = Features.drop(['POTENTIALLY_CAUSAL', 'VARIANT_ID'], axis=1)

    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state = 42, test_size=0.33)

    # In our data, the classes are imbalanced 
    # Handling Imbalanced Data by oversampling using SMOTE (Synthetic Minority Oversampling Technique)
    sm = SMOTE(random_state = 42)
    X_train_res, y_train_res = sm.fit_resample(X_train, y_train)

    # Instantiate and fit the RandomForest Classifier
    # Random forest classifier object
    forest = RandomForestClassifier(n_estimators = 100, random_state = 42, bootstrap = True)
    
    # Train the model on resampled train datasets
    forest.fit(X_train_res, y_train_res)

    # # Saving feature names for later use
    # feature_list = list(X.columns)

    # # Get numerical feature importances
    # importances = list(forest.feature_importances_)
    # # List of tuples with variable and importance
    # feature_importances = [(feature, round(importance, 2)) for feature, importance in zip(feature_list, importances)]
    # # Sort the feature importances by most important first
    # feature_importances = sorted(feature_importances, key = lambda x: x[1], reverse = True)
    # # Print out the feature and importances 
    # [print('Variable: {:20} Importance: {}'.format(*pair)) for pair in feature_importances]

    # # Make predictions
    # y_pred_test = forest.predict(X_test)

    # # Cross validation (CV) - 10 fold
    # forest_cv_score = cross_val_score(forest, X, y, cv=10, scoring='roc_auc')

    # print('--------------------------')
    # print("Confusion Matrix")
    # print('--------------------------')
    # print(confusion_matrix(y_test, y_pred_test), "\n")
    # print('--------------------------')
    # print("Classification Report")
    # print('--------------------------')
    # print(classification_report(y_test, y_pred_test), "\n")
    # print('--------------------------------------------------')
    # print("ROC_AUC scores computed at each CV iteration")
    # print('--------------------------------------------------')
    # print(forest_cv_score, "\n")
    # print('--------------------------')
    # print("Mean ROC_AUC score")
    # print('--------------------------')
    # print(forest_cv_score.mean(),"\n")

    # # save the model
    dump(forest, 'RandomForest_Model.joblib')

    logging.info("All done, completed successfully!")

    return 

###########################################################

# Parses the sample metadata file in .xlsx format
# Required columns are: 'sampleID', 'Causal gene' & 'pathologyID' 
# (can be in any order, but they MUST exist)
# Also, takes an input - path to the directory containing sample
# sample results file produced by the grexome-TIMC-secondary pipeline
#
# Extracts the causal variants for each sample (where gene is known)
#
# Returns a dictionary
# - Key: sampleID
# - Value: A list containing variant information such as:
#         - Chromosome Position, REF and ALT allels seperated by an '_'
#         - Gene Symbol
#         - Variant Impact
def CausalVariantParser(cohort, insample, indir):

    # Dictinary to store causal variants
    # key: POSITION, REF & ALT seperated by an Underscore
    # value: List containing SYMBOL, GENOTYPE
    causalVariant_dict = {}

    # Input sample metdata file
    sampleMetaFile = insample

     # columns of interest in sample metadata file
    Required_cols = ['sampleID', 'pathologyID', 'Causal gene']

    # Reading the samples metadata file into a dataframe
    sampleMetaData = pd.read_excel(sampleMetaFile, usecols = Required_cols, engine='openpyxl')
    
    # eliminating rows where sampleID = 0
    sampleMetaData = sampleMetaData[sampleMetaData['sampleID'] != 0]

    # eliminating rows where pathologyID does not match with the cohort
    # for which ML model is to be built
    sampleMetaData = sampleMetaData[sampleMetaData['pathologyID'] == cohort]

    # remove spaces in columns name
    sampleMetaData.columns = sampleMetaData.columns.str.replace(' ','_')

    # Dictionary to store sample and cohort info
    # key = sample ID, value = cohort this sample belongs to (ie pathologyID)
    sample2cohortR = pd.Series(sampleMetaData.pathologyID.values, index = sampleMetaData.sampleID).to_dict()

    # Dictionary to store sample and causal gene info (if known)
    # key = sample id, value = HGNC gene name of the known causal gene
    # also, do not store samples where causal gene is not known
    sample2causalR = pd.Series(sampleMetaData.Causal_gene.values, index = sampleMetaData.sampleID).dropna().to_dict()

    # parsing each sample result file in the indir that belongs to the 
    # cohort for which the ML model is to be built
    for path, dir, sampleFiles in os.walk(indir):
        for sampleF in sampleFiles:

            FName_fields = sampleF.split(".")

            # Getting the cohort of the sample result file
            sampleF_cohort = FName_fields[0]
            # Getting the sampleID of the sample result file
            sampleF_ID = FName_fields[1]

            if sampleF_ID in sample2cohortR:
                # Parse only those files where cohort matches the cohort
                # for which the model is to be built
                if sampleF_cohort == cohort:
                    sampleF_path = os.path.join(path, sampleF)
                    sampleFdata = open(sampleF_path)
                    # If this sample has a causal gene: see if and how it is hit
                    if sampleF_ID in sample2causalR:
                        
                        # get the causal gene
                        causalGene = sample2causalR[sampleF_ID]

                        # indexes of columns of interest
                        (symbolCol, ImpactCol, canonCol, positionCol, refcol, altcol, genotypecol) = (-1, -1, -1, -1, -1, -1, -1)

                        # grab header line
                        sampleF_header = sampleFdata.readline()
                        sampleF_headerFieds = sampleF_header.split('\t')

                        # grab column indexes of interest
                        for i in range(len(sampleF_headerFieds)):
                            if sampleF_headerFieds[i] == 'SYMBOL':
                                symbolCol = i
                            elif sampleF_headerFieds[i] == 'IMPACT':
                                ImpactCol = i
                            elif sampleF_headerFieds[i] == 'CANONICAL':
                                canonCol = i
                            elif sampleF_headerFieds[i] == 'POSITION':
                                positionCol = i
                            elif sampleF_headerFieds[i] == 'REF':
                                refcol = i
                            elif sampleF_headerFieds[i] == 'ALT':
                                altcol = i
                            elif sampleF_headerFieds[i] == 'GENOTYPE':
                                genotypecol = i
                        if not (symbolCol >= 0 and ImpactCol >= 0 and canonCol >= 0 and positionCol >= 0 and refcol >=0 and altcol >= 0 and genotypecol >= 0):
                            logging.error("Cannot find required column headers in the file: %s \n" % sampleF, "\n", sampleF_header, "\n")
                            sys.exit()
                        # else grabbed the required column indexes -> PROCEED

                        # Data lines
                        # find lines affecting a CANONICAL transcript of the causal gene in the sample result file
                        for line in sampleFdata:
                            line = line.rstrip("\n\r")
                            line_fields = line.split("\t")

                            # in the symbol column the gene name is in the format: ' Gene
                            # Ex: ' MIB2
                            # And We want only canonical transcripts
                            if not ((line_fields[symbolCol] == "' "+ causalGene) and (line_fields[canonCol] == 'YES')):
                                continue 
                            else: 
                                # Get the most severly hit variant, 
                                variant_info = [line_fields[positionCol]+'_'+ line_fields[refcol]+'_'+line_fields[altcol], line_fields[symbolCol], line_fields[ImpactCol]]
                                if line_fields[ImpactCol] == 'HIGH':
                                    causalVariant_dict[sampleF_ID] = variant_info
                                    break # We found the high impact variant
                                elif line_fields[ImpactCol] == 'MODHIGH':
                                    if sampleF_ID in causalVariant_dict:
                                        impact = causalVariant_dict[sampleF_ID][2]
                                        if not impact == 'HIGH':
                                            causalVariant_dict[sampleF_ID] = variant_info
                                    else:
                                        causalVariant_dict[sampleF_ID] = variant_info
                                elif line_fields[ImpactCol] == 'MODERATE':
                                    if sampleF_ID in causalVariant_dict:
                                        impact = causalVariant_dict[sampleF_ID][2]
                                        if not (impact == 'HIGH' or impact == 'MODHIGH'):
                                            causalVariant_dict[sampleF_ID] = variant_info
                                    else:
                                        causalVariant_dict[sampleF_ID] = variant_info
                                elif line_fields[ImpactCol] == 'LOW':
                                    if sampleF_ID in causalVariant_dict:
                                        impact = causalVariant_dict[sampleF_ID][2]
                                        if not (impact == 'HIGH' or impact == 'MODHIGH' or impact == 'MODERATE'):
                                            causalVariant_dict[sampleF_ID] = variant_info
                                    else:
                                        causalVariant_dict[sampleF_ID] = variant_info
                #else: sample doesn't have a causal gene
                
    return causalVariant_dict

###########################################################

# Taking and handling command-line arguments
def main():
    file_parser = argparse.ArgumentParser(description =
    """
---------------------------------------------------------------------------------------------------------------------
Program: Parses the canonical cohort result file (produced by the grexome-TIMC-Secondary pipeline). Also parses the 
         sample metadata file and sample results file to extract causal variant for each sample. Finally, builds a 
         Random Forest Model using different features. 
         
         The trained model will be saved to a file (RandomForest_Model.joblib) in the current working directory.
---------------------------------------------------------------------------------------------------------------------

Arguments [defaults] -> Can be abbreviated to shortest unambiguous prefixes
    """,
    formatter_class = argparse.RawDescriptionHelpFormatter)

    required = file_parser.add_argument_group('Required arguments')
    optional = file_parser.add_argument_group('Optional arguments')

    required.add_argument('--insample', metavar = "Input File", dest = "insample", help = 'Sample metadata file', required=True)
    required.add_argument('--incohort', metavar = "Input File", dest = "incohort", help = 'Single Cohort result file CSV (The cohort for which you would want to build the Machine learning model)', required=True)
    required.add_argument('--indir', metavar = "Input Directory", dest = "indir", help='Path to the directory containing sample result files CSV produced by the grexome-TIMC-Secondary pipeline')
 
    args = file_parser.parse_args()
    BuildModel(args)

if __name__ == "__main__":
    # Logging to Standard Error
    logging.basicConfig(format = "%(levelname)s %(asctime)s: %(filename)s - %(message)s", datefmt='%Y-%m-%d %H:%M:%S', stream = sys.stderr, level = logging.DEBUG)
    logging.addLevelName(logging.INFO, 'I' )
    logging.addLevelName(logging.ERROR, 'E')
    logging.addLevelName(logging.WARNING, 'W')
    main()
