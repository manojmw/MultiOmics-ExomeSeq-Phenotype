## Introduction

This is the main repository containing all the scripts for the project: MultiOmics-ExomeSeq-Phenotype. Here, I am developing a Multi-Omics method that scores genomic variants might be responsible for a particular phenotype in the patient (**The project is still in progress and I will add more scripts, such as those that provide scoring component for the Machine Learning step**). 
</br></br>
This repository contains individual scripts. I have integrated these into current **Exome-Seq Secondary Analysis Pipeline** ([click here](https://github.com/manojmw/grexome-TIMC-Secondary))
</br>

- [Example Usage](#example-usage-and-details-of-the-scripts)
   - [UniProt Parser](#uniprotparser)
   - [Protein-Protein Interaction Parser](#ppiparser)
   - [PPI Experiment Count](#ppiexpcount)
   - [Build Interactome](#interactome)
   - [Module Input File Generator](#modulefile)
   - [Uniprot2ENSG Mapper](#uniprotensgmapper)
   - [Machine Learning ScoreComponent-1](#MLScoreComp1)
- [Arguments](#Arguments)
- [Metadata files](#metadata-files)
- [Dependencies](#dependencies)
- [License](#license)

## Example Usage and Details of the Scripts

<a name="uniprotparser"></a>**UniProt Parser**

- Parses a UniProt file and extracts the required data from each record
- Populates the OUTPUT files with processed data

-> Grab the latest UniProt data with:

```console
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz
```

-> Parse UniProt data to produce output files with:

```console
gunzip -c uniprot_sprot.dat.gz | python3 1_Uniprot_parser.py --outPrimAC uniprot_main.tsv --outSecAC uniprot_secondary.tsv ---outGeneID geneID.tsv
```    

</br>

<a name="ppiparser"></a>**Protein-Protein Interaction Parser**

- Parses a Protein-Protein Interaction (PPI) File (miTAB 2.5 or 2.7)
- Maps to UniProt using the output files produced by `1_Uniprot_parser.py` and prints to STDOUT in .tsv format

-> Grab the latest PPI data (Ex: BioGRID) with:
```console
wget https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-ORGANISM-LATEST.mitab.zip
```  
-> Unzip with:
```console
unzip BIOGRID-ORGANISM-LATEST.mitab
```
-> This will produce one miTAB file per Organism (Use BIOGRID-ORGANISM-Homo_sapiens*.mitab.txt for human data)
-> Parse PPI data with:
```console
python3 2_Interaction_parser.py --inInteraction BIOGRID-ORGANISM-Homo_sapiens-4.4.207.mitab.txt --inPrimAC uniprot_main.tsv --inSecAC uniprot_secondary.tsv --inGeneID geneID.tsv > Exp_Biogrid.tsv
```

</br>

<a name="ppiexpcount"></a>**PPI Experiment Count**

- Parses a Protein-Protein Interaction File (miTAB 2.5 or 2.7)
- Prints the count of Human-Human Protein Interaction experiments to STDOUT

-> Provide a STDIN miTAB 2.5 or 2.7 file with:
```console
python3 3_Count_HumanPPIExp.py < BIOGRID-ORGANISM-Homo_sapiens-4.4.207.mitab.txt
```                   

</br>

<a name="interactome"></a>**Build Interactome**

- High-Quality Interactome Criteria:

    1] Filtering Interactions based on Interaction Detection Method:
      - We filter out pull down (MI:0096), genetic interference (MI:0254) & unspecified method (MI:0686)

    2] Filtering Interactions based on Interaction Type:
      - We keep only direct interaction (MI:0407) & physical association (MI:0915)

    3] We try to eliminate most of the SPOKE EXPANSION DATA, and consider only TRUE BINARY INTERACTIONS

    4] Each Interaction has ≥ 2 experiments, of which at least one of them should be proved by any BINARY METHOD

-> Build High-Quality Human Interactome with:      
```console
python3 4_BuildInteractome.py --inExpFile Exp_Biogrid.tsv --inPrimAC uniprot_main.tsv --inCanonicalFile canonicalTranscripts_*.tsv.gz > Interactome_human.tsv
```                      
-> For getting `canonical transcripts file`, please refer to [grexome-TIMC-Secondary](https://github.com/ntm/grexome-TIMC-Secondary/tree/master/Transcripts_Data)

</br>

<a name="modulefile"></a>**Module Input File Generator**

- Parses the output produced by `4_BuildInteractome.py`
- Assigns a weight to each interaction and prints to STDOUT in .tsv format
- This can be used as INPUT for most of the module identification/clustering methods

-> Generate Module Input File with:
```console
python3 5_ModuleInputFile.py < Interactome_human.tsv
```                      
</br>

<a name="uniprotensgmapper"></a>**UniProt2ENSG Mapper**

- Parses the output files produced by `1_Uniprot_parser.py` and the `canonical transcripts file` (Ex: canonicalTranscripts_220221.tsv)
- Maps UniProt accession to ENSG and prints to STDOUT

-> Run UniProt2ENSG Mapper with:
```console
python3 6_Uniprot2ENSG.py --inPrimAC uniprot_main.tsv --inCanonicalFile canonicalTranscripts_220221.tsv
```                      

</br>

<a name="MLScoreComp1"></a>**Machine Learning ScoreComponent-1**
- Parses the UniProt Primary Accession file (Ex:uniprot_main.tsv), canonical transcripts file (Ex: canonicalTranscripts_220221.tsv), Candidate Genes file (Ex: candidateGenes.xlsx) & Interactome file generated by `4_BuildInteractome.py`
- Checks the number of Interactors for each gene
- Checks the number of Interactors that are known candidate genes and calculates P-values
- Further, computes Benjamini-Hochberg adjusted P-value and prints to STDOUT in .tsv format
- This script provides the first scoring component for the Machine Learning step

-> Run Machine Learning ScoreComponent-1 script with:
```console
python3 7_Score1ML_candidateGenes.py --inPrimAC uniprot_main.tsv --inCandidateFile candidateGenes.xlsx --inCanonicalFile canonicalTranscripts_220221.tsv --inInteractome Interactome_human.tsv
```                      

## Arguments

**Arguments [defaults] -> Can be abbreviated to shortest unambiguous prefixes**

```shell
# UniProt Files
   --outPrimAC                          A tab-seperated UniProt Primary Accession Output File name containing UniProt Primary Accession, Taxonomy Identifier, ENST(s) & ENSG(s)
   --outSecAC                           A tab-seperated UniProt Secondary Accession Output File name containing UniProt Secondary Accession & Corresponding UniProt Primary Accession
   --outGeneID                          A tab-seperated GeneID Output File name (.tsv) containing GeneID & Corresponding UniProt Primary Accession

# Protein-Protein Interaction File(s)                                     
   --inInteraction                      miTAB 2.5 or 2.7 Input File name (Protein-Protein Interaction File)
   --inPrimAC                           Uniprot Primary Accession Input File name produced by 1_Uniprot_parser.py
   --inSecAC                            Uniprot Secondary Accession Input File name produced by 1_Uniprot_parser.py
   --inGeneID                           GeneID Input File name produced by 1_Uniprot_parser.py

# Interactome File(s)   
   --inExpFile                          PPI Experiments Input File name (produced by 2_Interaction_parser.py)

# Canonical Transcripts File
   --inCanonicalFile                    Canonical Transcripts Input File name (.gz or non .gz)

# Candidate Gene File(s)
   --inCandidateFile                    Candidate Gene Input Files(s) name (.xlsx)

# Interactome File
   --inInteractome                      High-Quality Interactome Input File name (produced by 4_BuildInteractome.py)

# Help
   -h, --help                           Show the help message and exit
```


## Metadata files

- Currently, the scripts use only one metadata file i.e. `candidateGenes.xlsx`. Later, description of other files will be added.

1. candidateGenes.xlsx: </br>
  * Lists known candidate genes. This eases the identification of a patient's likely causal variant: variants impacting a known candidate gene can be easily selected.  
  * Required columns: </br>
    - Gene: name of gene (should be the HGNC name, see www.genenames.org).
    - pathologyID: pathology/phenotype
    - Confidence score: indicates how confident you are that LOF variants in this gene are causal for this pathology. Value: integers from 1 and 5 (5 meaning the gene is definitely causal, while 1 is a lower-confidence candidate).


## Dependencies

* Requires **Python version ≥ 3**
* External dependencies are kept to minimum in all the scripts. The only required python modules are listed below: </br>
  - Pandas ([Installation Guide](https://pandas.pydata.org/docs/getting_started/install.html))
  - SciPy ([Installation Guide](https://scipy.org/install/))
* Most other standard core modules should already be available on your system

## License

Licensed under GNU General Public License v3.0 (Refer LICENSE file for more details)
