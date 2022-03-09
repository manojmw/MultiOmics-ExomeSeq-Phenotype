## Introduction

This is the main repository containing all the scripts for the project: MultiOmics-ExomeSeq-Phenotype. All the scripts are written in python. (**The project is still in progress**)
</br>
- [Example Usage](#example-usage-of-the-scripts)
   - [Uniprot Parser](#uniprotparser)
   - [Protein-Protein Interaction Parser](#ppiparser)
   - [PPI Experiment Count](#ppiexpcount) 
   - [Interactome generator](#interactome)
   - [Module Input File Generator](#modulefile)
   - [Uniprot2ENSG Mapper](#uniprotensgmapper)
- [Detailed Description(Arguments, Input Files and Output)](#detailed-description)
- [Metadata files](#metadata-files)
- [Dependencies](#dependencies)

## Example Usage of the Scripts

<a name="uniprotparser"></a>**Uniprot Parser**

```console
python 1_uniprot_parser.py --inUniprot uniprot_sprot.dat --outPrimAC uniprot_main.tsv --outSecAC uniprot_secondary.tsv --outGeneID geneID.tsv
```    
</br>

- This will basically parse a uniprot file (Ex: uniprot_sprot.dat)
- Extracts the required data from each record 
- Populates the output files with the processed data.

</br>

<a name="ppiparser"></a>**Protein-Protein Interaction Parser**

```console
python 2_interaction_parser.py --inInteraction intact.txt --inPrimAC uniprot_main.tsv --inSecAC uniprot_secondary.tsv --inGeneID geneID.tsv
```  
</br>

- Parses a miTAB 2.5 or 2.7 file (Ex: intact.txt)
- Maps the data for each protein-protein interaction experiment to the output files produced by `1_uniprot_parser.py`
- Extracts the uniprot accessions of the interacting proteins as well as other associated data
- Prints to STDOUT in tab-separated format

</br>

<a name="ppiexpcount"></a>**PPI Experiment Count**

```console
python 3_check_HumanPPIExp.py < intact.txt
```                      

</br>

- Parses a miTAB 2.5 or 2.7 file (Ex: intact.txt)
- Prints the count of Human-Human Protein Interaction experiments to STDOUT 

</br>

<a name="interactome"></a>**Interactome generator**

```console
python 4_Interactome.py --inCuratedFile curatedPPI_BioGRID.tsv curatedPPI_Intact.tsv --inPrimAC uniprot_main.tsv --inCanonicalFile canonicalTranscripts_220221.tsv
```                      

</br>

- Parses the output files (Ex: curatedPPI_BioGRID.tsv and curatedPPI_Intact.tsv) produced by `2_interaction_parser.py` and generates a high-quality human interactome
- Further, maps the uniprot accessions to ENSG using the canonical transcripts file and prints to STDOUT
- To produce the `canonical transcripts file`, please refer to [grexome-TIMC-Secondary](https://github.com/ntm/grexome-TIMC-Secondary/tree/master/Transcripts_Data)

</br>

<a name="modulefile"></a>**Module Input File Generator**

```console
python 5_ModuleInputFile.py < Interactome.tsv
```                      

</br>

- Parses the output (Ex: Interactome.tsv) produced by `4_Interactome.py` 
- Produces the input files that can be used for most of the module identification/clustering methods

</br>

<a name="uniprotensgmapper"></a>**Uniprot2ENSG Mapper**

```console
python 6_Uniprot2ENSG.py --inPrimAC uniprot_main.tsv --inCanonicalFile canonicalTranscripts_220221.tsv
```                      

</br>

- Parses the output files produced by `1_uniprot_parser.py` and the `canonical transcripts file` (Ex: canonicalTranscripts_220221.tsv)
- Maps uniprot accession to ENSG and prints to STDOUT

## Detailed Description

- For more detailed description on arguments, input files and the output generated, please use the help option with the scripts - `python script.py --help` OR `python script.py -h`

## Metadata files

- Currently, the scripts uses only one metadata file i.e. `candidateGenes.xlsx`. Later, description of other files will be added.

1. candidateGenes.xlsx: </br>
  * Lists known candidate genes. This eases the identification of a patient's likely causal variant: variants impacting a known candidate gene can be easily selected.  
  * Required columns: </br>
    - Gene: name of gene (should be the HGNC name, see www.genenames.org).
    - pathologyID: pathology/phenotype, as in the previous metadata files.
    - Confidence score: indicates how confident you are that LOF variants in this gene are causal for this pathology. We recommend using integers between 1 and 5, 5 meaning the gene is definitely causal while 1 is a lower-confidence candidate.


## Dependencies

* Requires **python ≥ 3**
* External dependencies are kept to minimum in all the scripts. The only required python modules are listed below: </br>
  - Pandas ([Installation Guide](https://pandas.pydata.org/docs/getting_started/install.html))
* Most other standard core modules should already be available on your system
</br>
