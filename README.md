Microbial Fatty Acid Compositions
=================================

Overview
--------

This repository contains microbial fatty acid compositions from the literature. Most of the data is linked to a species-level or below numeric taxonomy identifier from the NCBI (https://www.ncbi.nlm.nih.gov/taxonomy) and cross-literature lipid nomenclatural differences have been normalized to generate two summary tables: `lo_qual.tsv` and `md_qual.tsv`.


Fatty Acid Nomenclature
-----------------------

All raw and summary data uses a abbreviated fatty acid labelling nomenclature. Briefly, the length of the longest carbon chain is prefixed by C and followed by a colon and the number of unsaturated positions to form the core (e.g. myristic acid -> `C14:0`). If there are double bonds and their positions are known, a lower-case omega (ω) follows along with the positions, optionally ending with c/t (cis/trans) modifiers (e.g. palmitoleic acid -> `C16:1ω7c`). Before these descriptors any modifying functional groups are listed (hydroxy -> `OH`, keto -> `Oxo`, methyl -> `Me`) prefixed by their positions.


Summarization Procedure
-----------------------

Compositions marked as below the limit of quantification in the raw data were copied into the summary tables as 1/3 that limit of quantification (to ensure all values are numeric).

Some aldehydes, dicarboxylic acids, etc were quantified along with monocarboxylic in the original literature. These are included in the raw data, but were grouped into the Other category for the summary tables.


Files
-----

`reference_data/*.tsv` - Raw data copied from the literature
`downsample.py` - Script used to generate summary tables from the raw data
`references.tsv` - The citations and locations in the references where all of hte raw data came from
`lo_qual.tsv` - A quality filtered and merged version of all of the raw data with poorer-quality lipid descriptions (neither constitutional isomers nor not stereochemical isomers seperated)
`md_qual.tsv` - A quality filtered and merged version of all of the raw data with higher-quality lipid descriptions (constitutional isomers separate but not stereochemical isomers)

Applications
------------

Some ideas for projects with this dataset:

 1. _Evolutionary Biology_ What correlations are there between taxonomic trees derived from marker genes or whole genome analysis and trees derived from metabolic data like fatty acid abundances? How conserved is membrane lipid composition phylogenetically and how far up the tree of life does this conservation go?

 2. _Environmental Microbiology & Geology_ Can you make up roughly "independent" vectors of lipid compositions at some taxonomic level and then predict the taxonomic composition of a mixed environmental sample solely through analysing lipid abundances? Could you apply this same treatment to a geological sample where this isn't a known truth to determine the taxonomic composition of extinct microbial ecosystems?

 3. _Systems Biology_ What correlations are there between lipid-synthesis genes present in sequenced genomes and actual fatty acid abundances? How much of an organisms membrane composition is set solely through synthetic genes versus regulatory elements, etc.?

License
-------
This compilation is released under the same license as Wikipedia:
https://creativecommons.org/licenses/by-sa/4.0/

There is no associated paper to cite (at this time), but you can cite the GitHub repository if you wish.
