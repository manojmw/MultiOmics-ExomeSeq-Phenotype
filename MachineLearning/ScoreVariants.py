#!/usr/bin/python

# manojmw
# 17 July, 2022

import argparse, sys
import logging
import pandas as pd
from joblib import load

###########################################################

# Parses the sample result file (produced by the grexome-TIMC-secondary 
# pipeline) 
# Also loads the model built using Build_RFModel.py script.
#
# Preprocessing of the data and predicting causal variants
#
# Prints to STDOUT the same data present in sample result file 
# but with an additional column "MODEL_CAUSAL_SCORE" inserted at the beginning
# The data is also sorted (descending) by this column
def ScoreVariants(args):

    logging.info("starting to run")

    # Load the model
    forest = load(args.inmodel)

    # sample result file
    sampleresultF = args.insampleresult
    
    # Getting the exact file name excluding path
    try:
        SampleFileName_fields = sampleresultF.split('/')
        Sample_FName = SampleFileName_fields[-1]
    except:
        Sample_FName = sampleresultF

    # Getting the pathology name from the file name Ex: MMAF.grexT1297.HG35FT.GATK
    FName_fields = Sample_FName.split('.')

    # Get the cohort name
    cohort = FName_fields[0]

    if not cohort:
        logging.error("Could not identify the cohort for the file: %s" % sampleresultF, "Please fix the name of the file!")
        sys.exit()

    # Storing the data in dataframe
    Data = pd.read_csv(sampleresultF, sep='\t', low_memory = False)

    # Dropping last 2 columns
    # Ex: max_ctrl_hv=3 max_ctrl_het=10 min_hr=117 no_mod max_af_gnomad=0.01 max_af_1kg=0.03
    Data = Data.iloc[:, :-2]

    # If PVALUE is missing: we assign p-value as 1
    Data["INTERACTORS_PVALUE"] = Data["INTERACTORS_PVALUE"].fillna(1)
    Data["ENRICHED_CLUSTER_PVALUE"] = Data["ENRICHED_CLUSTER_PVALUE"].fillna(1)

    # Unwanted columns for prediction
    uwcolumns = ['POSITION', 'REF', 'ALT', 'SYMBOL', 'KNOWN_CANDIDATE_GENE', 'GENOTYPE', 'DP:AF/BF:RR', 'BIALLELIC', 'INTERACTORS_COUNT', 'INTERACTORS', 'ENRICHED_CLUSTER', 'ENRICHED_CLUSTER_ID', 'ENRICHED_CLUSTER_SIZE', 'ENRICHED_CLUSTER_PVALUE', 'IMPACT', 'Consequence', 'HGVSc', 'HGVSp', 'Protein_position', 'GTEX_ovary_RATIO', 'GTEX_ovary', 'GTEX_blood', 'GTEX_cerebellar_hemisphere', 'GTEX_liver', 'Gene', 'Feature', 'CANONICAL', 'BIOTYPE', 'VARIANT_CLASS', 'RefSeq', 'MANE_SELECT', 'MANE_PLUS_CLINICAL', 'ALLELE_NUM', 'EXON', 'INTRON', 'cDNA_position', 'CDS_position', 'SIFT', 'PolyPhen', 'MetaRNN_pred', 'MetaRNN_rankscore', 'CADD_raw_rankscore', 'MutationTaster_pred', 'REVEL_rankscore', 'ada_score', 'rf_score', 'SpliceAI_pred_DS_AG', 'SpliceAI_pred_DP_AG', 'SpliceAI_pred_DS_AL', 'SpliceAI_pred_DP_AL', 'SpliceAI_pred_DS_DG', 'SpliceAI_pred_DP_DG', 'SpliceAI_pred_DS_DL', 'SpliceAI_pred_DP_DL', 'AF', 'MOTIF_NAME', 'MOTIF_POS', 'HIGH_INF_POS', 'MOTIF_SCORE_CHANGE', 'TRANSCRIPTION_FACTORS', 'Existing_variation', 'CLIN_SIG', 'SOMATIC', 'PHENO', 'GTEX_Brodmann_(1909)_area_24', 'GTEX_Brodmann_(1909)_area_9', 'GTEX_C1_segment_of_cervical_spinal_cord', 'GTEX_EBV-transformed_lymphocyte', 'GTEX_adrenal_gland', 'GTEX_amygdala', 'GTEX_aorta', 'GTEX_atrium_auricular_region', 'GTEX_breast', 'GTEX_caudate_nucleus', 'GTEX_cerebellum', 'GTEX_cerebral_cortex', 'GTEX_coronary_artery', 'GTEX_cortex_of_kidney', 'GTEX_ectocervix', 'GTEX_endocervix', 'GTEX_esophagogastric_junction', 'GTEX_esophagus_mucosa', 'GTEX_esophagus_muscularis_mucosa', 'GTEX_fallopian_tube', 'GTEX_greater_omentum', 'GTEX_heart_left_ventricle', 'GTEX_hippocampus_proper', 'GTEX_hypothalamus', 'GTEX_lower_leg_skin', 'GTEX_minor_salivary_gland', 'GTEX_nucleus_accumbens', 'GTEX_pancreas', 'GTEX_pituitary_gland', 'GTEX_prostate_gland', 'GTEX_putamen', 'GTEX_sigmoid_colon', 'GTEX_skeletal_muscle_tissue', "GTEX_small_intestine_Peyer's_patch", 'GTEX_spleen', 'GTEX_stomach', 'GTEX_subcutaneous_adipose_tissue', 'GTEX_substantia_nigra', 'GTEX_suprapubic_skin', 'GTEX_thyroid_gland', 'GTEX_tibial_artery', 'GTEX_tibial_nerve', 'GTEX_transformed_skin_fibroblast', 'GTEX_transverse_colon', 'GTEX_urinary_bladder', 'GTEX_uterus', 'GTEX_vagina']

    # Fill missing values with 0
    Data = Data.fillna(0)

    Test_dataset = Data.drop(uwcolumns, axis=1)

    # Make predictions
    Model_Causal_Predprob = forest.predict_proba(Test_dataset)[:,1]

    Data.insert(0, 'MODEL_CAUSAL_SCORE', Model_Causal_Predprob)

    Data = Data.sort_values(by = 'MODEL_CAUSAL_SCORE', ascending = False)

    # Replace 0 with NA
    Data = Data.replace(0, " ")

    # now Assign missing scores as 0
    Data['MODEL_CAUSAL_SCORE'] = Data['MODEL_CAUSAL_SCORE'].replace(" ", 0)

    # Printing to STDOUT
    print(Data.to_csv(sys.stdout, index = False, line_terminator = '\r'))

    logging.info("All done, completed successfully!")

    return

###########################################################

# Taking and handling command-line arguments
def main():
    file_parser = argparse.ArgumentParser(description =
    """
---------------------------------------------------------------------------------------------------------------------
Program: Parses the sample result file (produced by the grexome-TIMC-Secondary pipeline). Loads the model built with 
         Build_RFModel.py script. Predicts the causal variants and prints to STDOUT in .csv format.
---------------------------------------------------------------------------------------------------------------------
- The output consists of same data present in sample result file, but with an additional column "MODEL_CAUSAL_SCORE" 
  added at the the beginning.
- A single integrative score is assigned between 0 and 1 indicating how likely the variant is potentially causal 
  in the patient. 
- Higher score: more likely that the variant is causal
- Lower score: less likely that the variant is causal
- The data is also sorted by the "MODEL_CAUSAL_SCORE" column in the descending order
         
---------------------------------------------------------------------------------------------------------------------

Arguments [defaults] -> Can be abbreviated to shortest unambiguous prefixes
    """,
    formatter_class = argparse.RawDescriptionHelpFormatter)

    required = file_parser.add_argument_group('Required arguments')
    optional = file_parser.add_argument_group('Optional arguments')

    required.add_argument('--inmodel', metavar = "Input File", dest = "inmodel", help = 'The trained model file (.joblib) generated by Build_RFModel.py script', required=True)
    required.add_argument('--insampleresult', metavar = "Input File", dest = "insampleresult", help = 'Sample result file (.csv - the cohort should match the trained model cohort, else you might get incorrect predictions!!!) for which you want to score the variants', required=True)
 
    args = file_parser.parse_args()
    ScoreVariants(args)

if __name__ == "__main__":
    # Logging to Standard Error
    logging.basicConfig(format = "%(levelname)s %(asctime)s: %(filename)s - %(message)s", datefmt='%Y-%m-%d %H:%M:%S', stream = sys.stderr, level = logging.DEBUG)
    logging.addLevelName(logging.INFO, 'I' )
    logging.addLevelName(logging.ERROR, 'E')
    logging.addLevelName(logging.WARNING, 'W')
    main()
