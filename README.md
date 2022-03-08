## Introduction

This is the main repository containing all the scripts for the project: 

"Developing a Multi-Omics Integrative Analysis method to score Exome-Seq derived Genomic Variants that are causal for the Patient's Phenotype"

</br>

## Example Usage of the Scripts:

> 1_uniprot_parser.py

</br>

```console
python 1_uniprot_parser.py --inUniprot uniprot_sprot.dat --outPrimAC uniprot_main.tsv --outSecAC uniprot_secondary.tsv --outGeneID geneID.tsv
```    
</br>

- This will basically parse a uniprot file (Ex: uniprot_sprot.dat)
- Extracts the required data from each record
- Populates the output files with the processed data
- For more detailed description on the output files generated, please use the help option with the script: </br> `python 1_uniprot_parser.py --help`

</br>

> 2_interaction_parser.py

</br>

```console
python 2_interaction_parser.py --inInteraction intact.txt --inPrimAC uniprot_main.tsv --inSecAC uniprot_secondary.tsv --inGeneID geneID.tsv
```  
</br>

- Parses a miTAB 2.5 or 2.7 file (Ex: intact.txt)
- Maps the data for each protein-protein interaction experiment to the output files produced by `1_uniprot_parser.py`
- Extracts the uniprot accessions of the interacting proteins as well as other associated data
- Prints to STDOUT in tab-separated format
- For more detailed description on the arguments and the output, please use the help option with the script: </br> `python 2_interaction_parser.py --help`

</br>

> 3_check_HumanPPIExp.py

</br>

```console
python 3_check_HumanPPIExp.py < intact.txt
```                      

</br>

- Parses a miTAB 2.5 or 2.7 file (Ex: intact.txt)
- Prints the count of Human-Human Protein Interaction experiments to STDOUT 

</br>

> 4_Interactome.py

</br>

```console
python 4_Interactome.py --inCuratedFile curatedPPI_BioGRID.tsv curatedPPI_Intact.tsv --inPrimAC uniprot_main.tsv --inCanonicalFile canonicalTranscripts_220221.tsv
```                      

</br>

- Parses the output files (Ex: curatedPPI_BioGRID.tsv and curatedPPI_Intact.tsv) produced by 2_interaction_parser.py and generates a high-quality human interactome
- Further, maps the uniprot accessions to ENSG using the canonical transcripts file and prints to STDOUT
- For more detailed description on the arguments and the output, please use the help option with the script: </br> `python 4_Interactome.py --help`

</br>

> 5_ModuleInputFile.py

</br>

```console
python 5_ModuleInputFile.py < Interactome.tsv
```                      

</br>

- Parses the output (Ex: Interactome.tsv) produced by 4_Interactome.py 
- For more detailed description on the arguments and the output, please use the help option with the script: </br> `python 5_ModuleInputFile.py --help`

</br
