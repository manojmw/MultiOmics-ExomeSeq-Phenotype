## Introduction

This is the main repository containing all the scripts for my Master's thesis project. Here, I am developing a Multi-Omics method to score genomic variants that might be causal for a particular phenotype in the patient. (**Currently in progress, I will add more scripts - such as those that provide scoring component for the Machine Learning step**). 
</br></br>
This repository contains individual scripts. I have integrated these into the current **Exome-Seq Secondary Analysis Pipeline** ([click here](https://github.com/manojmw/grexome-TIMC-Secondary))

</br>

- [Example Usage](#example-usage)
   - [UniProt Parser](#uniprotparser)
   - [Protein-Protein Interaction Parser](#ppiparser)
   - [PPI Experiment Count](#ppiexpcount)
   - [Build Interactome (True Binary Interactions only)](#interactomebinary)
   - [Build Interactome (True Binary Interactions with Expansion)](#interactomebinarywithexpansion)
   - [Module Input File Generator](#modulefile)
   - [Uniprot2ENSG Mapper](#uniprotensgmapper)
   - [Naïve Approach](#NaïveApproach)
   - [DREAM Challenge: Cluster File Processing](#ProcessMonet)
   - [Naïve with Clustering Approach](#NaïvewithClusteringApproach)
- [Output](#output)
- [Interactome Clustering Methods](#Interactome-Clustering-Methods)
- [Arguments](#Arguments)
- [Metadata files](#metadata-files)
- [Dependencies](#dependencies)
- [License](#license)

## Example Usage

<a name="uniprotparser"></a>**UniProt Parser**

- Parses on STDIN a UniProt file and extracts the required data from each record
- Prints to STDOUT in .tsv format

-> Grab the latest UniProt data with:

```console
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz
```

-> Parse UniProt data to produce output with:

```console
gunzip -c uniprot_sprot.dat.gz | python3 Uniprot_parser.py > Uniprot_output.tsv
```    

</br>

<a name="ppiparser"></a>**Protein-Protein Interaction Parser**

- Parses a Protein-Protein Interaction (PPI) File (miTAB 2.5 or 2.7)
- Maps to UniProt using the output files produced by `1_Uniprot_parser.py` and prints to STDOUT in .tsv format

1] Grab the latest **BioGRID data** with:
   ```console
   wget https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-ORGANISM-LATEST.mitab.zip
   ```  
   -> Unzip with:
   ```console
   unzip BIOGRID-ORGANISM-LATEST.mitab.zip
   ```
   -> This will produce one miTAB file per Organism (Use BIOGRID-ORGANISM-Homo_sapiens*.mitab.txt for human data)

   </br>

2] Grab the latest **IntAct data** with:
   ```console
   wget ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.zip
   ```  
   -> Unzip with:
   ```console
   unzip intact.zip
   ```
   -> This will produce 2 files (intact.txt & intact_negative.txt). Use intact.txt for further steps

</br>

3] Parse PPI data with:
   ```console
   python3 Interaction_parser.py --inInteraction BIOGRID-ORGANISM-Homo_sapiens*.mitab.txt --inUniprot Uniprot_output.tsv > Exp_Biogrid.tsv
   ```
   ```console
   python3 Interaction_parser.py --inInteraction intact.txt --inUniprot Uniprot_output.tsv > Exp_Intact.tsv
   ```

-> The above example is for the Protein-Protein Interaction data from BioGRID and IntAct. But, you can retrieve PPI data (in miTAB format) from any database and feed it to the script to produce output file.

</br>

<a name="ppiexpcount"></a>**PPI Experiment Count**

- Parses a Protein-Protein Interaction File (miTAB 2.5 or 2.7)
- Prints the count of Human-Human Protein Interaction experiments to STDOUT

-> Provide a STDIN miTAB 2.5 or 2.7 file with:
```console
python3 3_Count_HumanPPIExp.py < miTAB File
```                   

</br>

<a name="interactomebinary"></a>**Build Interactome (True Binary Interactions only)**

- High-Quality Interactome Criteria:

    1] Filtering Interactions based on Interaction Detection Method:
      - We filter out pull down (MI:0096), genetic interference (MI:0254) & unspecified method (MI:0686)

    2] Filtering Interactions based on Interaction Type:
      - We keep only direct interaction (MI:0407) & physical association (MI:0915)

    3] Here, we try to eliminate most of the EXPANSION DATA, and consider only TRUE BINARY INTERACTIONS

    4] Each Interaction has ≥ 2 experiments, of which at least one of them should be proved by any BINARY METHOD

-> Build High-Quality Human Interactome with:      
```console
python3 4_BuildInteractome_BinaryPPIonly.py --inExpFile Exp_Biogrid.tsv Exp_Intact.tsv --inUniprot Uniprot_output.tsv --inCanonicalFile canonicalTranscripts_*.tsv.gz > Interactome_human_binaryonly.tsv
```                      
-> For getting `canonical transcripts file`, please refer to [grexome-TIMC-Secondary](https://github.com/ntm/grexome-TIMC-Secondary/tree/master/Transcripts_Data)
-> Build Interactome scripts accepts multiple processed Protein-Protein Interaction experiment file (--inExpFile)


</br>

<a name="interactomebinarywithexpansion"></a>**Build Interactome (True Binary Interactions with Expansion)**

- High-Quality Interactome Criteria:

    1] Filtering Interactions based on Interaction Detection Method:
      - We filter out genetic interference (MI:0254) & unspecified method (MI:0686)

    2] Here, we consider both TRUE BINARY INTERACTIONS and PPIs derived from EXPANSION

    4] Each Interaction should be proven by ≥ 2 experiments

-> Build High-Quality Human Interactome with:      
```console
python3 5_BuildInteractome_BinaryPPIwithExpansion.py --inExpFile Exp_Biogrid.tsv Exp_Intact.tsv --inUniprot Uniprot_output.tsv --inCanonicalFile canonicalTranscripts_*.tsv.gz > Interactome_human_binarywithexpansion.tsv
```                      

</br>

<a name="modulefile"></a>**Module Input File Generator**

- Parses the output produced by `4_BuildInteractome_BinaryPPIonly.py` or `5_BuildInteractome_BinaryPPIwithExpansion.py`
- Assigns a default edge weight = 1 for each interaction and prints to STDOUT in .tsv format
- This can be used as INPUT for most of the module identification/clustering methods

-> Generate Module Input File with:
```console
python3 6_ModuleInputFile.py < Interactome_human_binaryonly.tsv
```                      
</br>

<a name="uniprotensgmapper"></a>**UniProt2ENSG Mapper**

- Parses the output files produced by `1_Uniprot_parser.py` and the `canonical transcripts file` (Ex: canonicalTranscripts_220221.tsv)
- Maps UniProt accession to ENSG and prints to STDOUT

-> Run UniProt2ENSG Mapper with:
```console
python3 7_Uniprot2ENSG.py --inUniprot Uniprot_output.tsv --inCanonicalFile canonicalTranscripts_220221.tsv
```                      
</br>

<a name="NaïveApproach"></a>**Naïve Approach**
- Parses the Sample metadata file (.xlsx), UniProt File, Canonical transcripts file, Candidate Gene file(s), Interactome file and GTEX File
- For a given gene:
   * Checks if the gene is already a known candidate
   * Checks the number of Interactors 
   * Checks the number of Interactors that are known candidates
   * Eliminates Hub/Sticky proteins
   * Next, applies Fisher's Exact test to compute P-values
   * Additionally adds GTEX data
   * Prints to STDOUT in .tsv format
- This script provides the first scoring component for the Machine Learning step

-> Run 8_NaiveApproach.py script with:
```console
python3 8_NaiveApproach.py --inSampleFile sample.xlsx --inUniprot Uniprot_output.tsv --inCandidateFile candidateGenes.xlsx --inCanonicalFile canonicalTranscripts_220221.tsv --inInteractome Interactome_human.tsv --inGTEXFile E-MTAB-5214-query-results.tpms.tsv
```      
-> You can use the GTEX file provided in this repository. (**Note**: The GTEX file provided in the repository might not be the latest. If you want to retrieve the latest GTEX file, please visit https://www.ebi.ac.uk/gxa/home).  

</br>

<a name="ProcessMonet"></a>**DREAM Challenge: Cluster File Processing**
- This script should be run only if you are using the cluster file produced by the MONET TOOL (DREAM Challenge)
- Parses on STDIN a .tsv file produced by the MONET tool, processes it and prints to STDOUT in .cls format
- The output can be used as the input Cluster File for `10_Naive_withClusteringApproach.py`
- **Note:** For producing Interactome Clustering File and using this script, please refer to the "Interactome Clustering Methods" section

</br>

<a name="NaïvewithClusteringApproach"></a>**Naïve with Clustering Approach**
- This script is similar to `8_NaiveApproach.py`, but the output additionally contains Interactome Clustering data

-> Run 10_Naive_withClusteringApproach.py script with:
```console
python 10_Naive_withClusteringApproach.py --inSampleFile sample.xlsx --inUniprot Uniprot_out.tsv --inCandidateFile candidateGenes_*.xlsx --inCanonicalFile canonicalTranscripts_220221.tsv --inInteractome Interactome_human --inClusterFile K1Clustering_clusterFile.cls --inGTEXFile E-MTAB-5214-query-results.tpms.tsv
```

## Output

- For detailed description on the output generated by the scripts, please use the --help, -h option.
- You can also view the sample output files provided in the Sample_Output_Files directory of the repository


## Interactome Clustering Methods

- We consider Clusters with a size of >= 2 and 100 (max)
- If the cluster size exceeds 100, the methods are applied recursively (MONET tool automatically does this) to obtain clusters of desired size.
- I have tested mainly three types of clustering methods:

   1] Kernel clustering approach (top ranking method K1 from DREAM Challenge) (Choobdar, Sarvenaz et al. “Assessment of network module identification across complex diseases.” Nature methods vol. 16,9 (2019): 843-852. doi:10.1038/s41592-019-0509-5)
   
   
   2] Modularity Optimization method (top ranking method M1 from DREAM Challenge) (Choobdar, Sarvenaz et al. “Assessment of network module identification across complex diseases.” Nature methods vol. 16,9 (2019): 843-852. doi:10.1038/s41592-019-0509-5)

   </br>
   
       - To run this clustering method on the Interactome file generated by Build_Interactome.py, please use the MONET tool described at: https://github.com/BergmannLab/MONET

       - If you will be using the cluster file produced by these 2 methods, then please process it using ProcessClusterFile_MONET.py script using the command:

         % cat cluster_outputFile.tsv | python3 9_ProcessClusterFile_MONET.py > newClusterFile.cls

       - Input File description:

            cluster_outputFile.tsv:    Clustering Output File produced by the MONET tool (DREAM Challenge)



   3] Randomized optimization of modularity (Didier, Gilles et al. “Identifying communities from multiplex biological networks by randomized optimization of modularity.” F1000Research vol. 7 1042. 10 Jul. 2018, doi:10.12688/f1000research.15486.2)
      
   </br>
         
      - To run this clustering method on the Interactome file generated by Build_Interactome.py, please use the MolTi-DREAM tool described at: https://github.com/gilles-didier/MolTi-DREAM

      - This might generate some large clusters (i.e size > 100), in such cases, please run the tool recurscively as described at MolTi-DREAM GitHub page

      - The output produced by this tool need not be processed further and can be directly used for the 5.2_addInteractome.py script

  
- You can use one of the above or any other clustering methods, but the Cluster File (please refer to the sample_clusterFile.cls file in the Sample_Input_Files directory of the repository) should be of the format:

   </br>

      - Header: (Ex: #ClustnSee analysis export)
      - Followed by ClusterID (Ex: ClusterID:1||)
      - Followed by Name(ENSG) of the Cluster(s) (Ex: ENSG00000162819)
      - End of a given Cluster is indicated by an empty line


## Arguments

**Arguments [defaults] -> Can be abbreviated to shortest unambiguous prefixes**

```shell
# UniProt Files
   --inUniProt                          A tab-seperated Input File name (produced by 1_Uniprot_parser.py) containing UniProt Primary Accession, Taxonomy Identifier, ENST(s), ENSG(s), UniProt Secondary Accession(s), Gene ID(s) & Gene name(s)

# Protein-Protein Interaction File(s)                                     
   --inInteraction                      miTAB 2.5 or 2.7 Input File name (Protein-Protein Interaction File)

# Interactome File(s)   
   --inExpFile                          PPI Experiments Input File name (produced by 2_Interaction_parser.py)

# Canonical Transcripts File
   --inCanonicalFile                    Canonical Transcripts Input File name (.gz or non .gz)
   
# Sample File
   --inSampleFile                       Sample Metadata Input File name (.xlsx)   

# Candidate Gene File(s)
   --inCandidateFile                    Candidate Gene Input Files(s) name (.xlsx)

# Interactome File
   --inInteractome                      High-Quality Interactome Input File name (produced by 4_BuildInteractome_BinaryPPIonly.py/5_BuildInteractome_BinaryPPIwithExpansion.py)
   
# GTEX File
   --inGTEXFile                         GTEX Input File name (.tsv)

# Help
   -h, --help                           Show the help message and exit
```


## Metadata files

- Currently, the scripts use only one metadata file i.e. `candidateGenes.xlsx`. Later, description of other files will be added.

1. samples.xlsx: </br>
   * This metadata file describes the samples. 
   * Required column: </br>
      - pathologyID: the phenotype of each patient/sample, used to define the "cohorts".
   * Optional columns (currently not used by the scripts) such as:
      - sampleID: unique identifier for each sample
      - specimenID: external identifier for each sample, typically related to the BAM or FASTQ filenames.
      - patientID: a more user-friendly identifier for each sample
      - Sex: 'F' or 'M'
      
2. candidateGenes.xlsx: </br>
   * Lists known candidate genes. This eases the identification of a patient's likely causal variant: variants impacting a known candidate gene can be easily selected.  
   * Required columns: </br>
      - Gene: name of gene (should be the HGNC name, see www.genenames.org).
      - pathologyID: pathology/phenotype
   * Optional column (currently not used by the scripts) such as:
      - Confidence score: indicates how confident you are that LOF variants in this gene are causal for this pathology. Value: integers from 1 and 5 (5 meaning the gene is definitely causal, while 1 is a lower-confidence candidate).


## Dependencies

* Requires **Python version >= 3**
* External dependencies are kept to minimum in all the scripts. The only required python modules are listed below: </br>
  - OpenPyXl 
  - SciPy
* You can install these with pip/conda Ex: pip3 install openpyxl scipy OR conda install openpyxl scipy
* Most other standard core modules should already be available on your system

## License

Licensed under GNU General Public License v3.0 (Refer LICENSE file for more details)
