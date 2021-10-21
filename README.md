################################################################################
README
################################################################################

*Title:*
IUCN Red List protects avian genetic diversity

*Authors:*
**Elisabetta Canteri**, Damien A. Fordham, Sen Li, Peter A. Hosner, Carsten Rahbek, and David Nogues-Bravo

*Contacts:*
Elisabetta Canteri -
The Environment Institute and School of Biological Sciences, The University of Adelaide, 5005 Adelaide (AU)
elisabetta.canteri@adelaide.edu.au

David Nogues-Bravo -
GLOBE Institute, University of Copenhagen, 2100 Copenhagen (DK)
dnogues@sund.ku.dk

################################################################################

The dataset contains the following tables and code:

1. **Master_db.csv** -
This is a table containing the information of the cytochrome-b sequences of birds available in GenBank. The table has the following columns:
* ACC_NUM = Accession number - sequence identifier in GenBank.
* ORDER = Order which the individual (sequence) belongs to.
* OscSubOsc = whether the individual (sequence) belongs to the Oscines or the Suboscines clade.
* FAMILY = Family which the individual (sequence) belongs to.
* GB_NAME = Name of the organism associated to the sequence, as reported in GenBank. This can be at the species or subspecies level.
* GB_SPECIES = Name of the organism cut to the species level. This column is the result of processing on "GB_NAME".
* IOC_SPECIES = Taxonomic standardisation of "GB_SPECIES" to the IOC v2.2 taxonomy (Gill et al. 2009).
* JETZ_NEW = Taxonomic standardisation of "GB_SPECIES" to the BirdLife International Checklist v3 (2010) (BirdLife International 2010).
* NEW_NAME_PHYLO = For some sequences the taxonomy was corrected using phylogenetic trees. New names were assigned based on the position of the sequences within the phylogenetic tree.
* FINAL_NAME = Species names after the taxonomic reconciliation.
* RED_LIST_JETZ.2017. = IUCN Red List categories assigned to the new names, according to the HBW-BirdLife Checklist v2 (2017) (HBW and BirdLife International 2017).
* NOTES = Notes about sequences. Sequences belonging to domesticated individuals are noted as "DOM", hybrids as "HYB", extinct as "EX", and duplicates as "IDENTICAL".
* GENE = mtDNA gene.
* NUM_BP = length of the sequence (number of base-pairs)

2. **dataset_for_matlab.csv** - subset of the Master_db.csv table, which contains only the sequences to be used in the calculation of nucleotide diversity. Columns specifications are reported above.

3. **matlab_output.csv** - output table from Matlab, with nucleotide diversity values for each species. Column names were added manually to the output table from Matlab. The columns are:
* SP = species name
* GD = nucleotide diversity value
* NUM_MUT = number of mutations
* DROPPED_PAIRS = number of paired sequences not included in the calculation
* SEQS = total number of sequences for the species

4. **dataset_analysis.csv** - table that is used in the analyses. This table associates the IUCN Red List categories to the matlab_output.csv table. Species with < 5 sequences are removed.

5. **SpeciesSequenceLengths.csv** - table with maximum, minimum and average length of sequences for each species. The unit for the lengths is number of base-pairs.

6. **Final_Alignments.fasta** - Checked and corrected species-specific alignments used in the analyses.

7. **AvianGD.Rproj** - R project linked to the R code used for data processing and analysing.

8. **dataset_for_matlab.R** - code to generate the dataset that is used in Matlab for the calculation of nucleotide diversity.

9. **prepare_matlab_imput.sh** - code to process sequence alignments to the format required by Matlab. Alignments are saved as single files in a new folder.

10. **FunPairSegSites.m** and **nuc_div.m** - code to calculate nucleotide diversity in Matlab.

11. **dataset_paper_analysis.R** - code to generate the data used in the genetic diversity analyses.

12. **resampling_function.R** and **phylANOVA_modified.R** - code to calculate perform phylANOVAs and test for differences in intra-specific genetic diversity between threatened and non-threatened species.

13. **paper_analysis.R** - code to perform the analyses on genetic diversity and for producing the plots. This is the main code used to obtain the results presented in the manuscript.

14. **Realms.R** - code used to calculate the % of species belonging to each Zoogeographic realm.

################################################################################

To replicate the analyses, please follow these steps:
1. Open **AvianGD.Rproj** in R and run **dataset_for_matlab.R**.
2. Run **prepare_matlab_imput.sh** with **dataset_for_matlab.csv** and **Final_Alignments.fasta**
3. Run **FunPairSegSites.m** and **nuc_div.m** in Matlab.
4. Run **dataset_paper_analysis.R** in R.
5. Run **paper_analysis.R** in R, sourcing **resampling_function.R** and **phylANOVA_modified.R**.
6. Run **Realms.R** in R.

The phylogenetic tree used in the analyses is the MC_Hackett2_Full.tre, which can be found in Jetz et al. (2012) - https://www.nature.com/articles/nature11631

################################################################################
*References*:

BirdLife International 2010. The BirdLife checklist of the birds of the world, with conservation status and taxonomic sources. Version 3.

Gill, F. et al. 2009. IOC World Bird List (v2.2). - IOC

HBW and BirdLife International 2017. Handbook of the Birds of the World and BirdLife International digital checklist of the birds of the world. Version 2.

Jetz, W. et al. 2012. The global diversity of birds in space and time. - Nature 491: 444–448.

################################################################################
