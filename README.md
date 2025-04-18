# The diversity and spread of cytoplasmic incompatibility genes among maternally inherited symbionts: *Wolbachia* do not walk alone - GitHub page

Amoros J, Buysse M, Floriano AM, Moumen B, Vavre F, Bouchon D, Duron D

This GitHub page shares the command lines and additional data used in this study. It is organized according to the typical workflow used to screen *cifA-cifB* genes from `.fasta` bacterial genome files retrieved from new assemblies or public databases, verify the presence of *cif* genes, determine their flanking genetic environment, visualize the results, and perform phylogenetic analyses to place each element in an evolutionary context.

The page is divided as follows: 

#### 1. Assembling of new genomes and collection of publicly available genomes
This section provides details about the methods used to assemble the eight *Rickettsia* genomes from metagenomes of tick organs (assembly, correction, binning, taxonomic assignation, genome representations, etc.).

#### 2. Genome annotation
This step analyzes gene boundaries and provides structural context for each genome using the Prokka annotation tool.

#### 3. Screening of *cifA* and *cifB* and their predicted protein domains
This part presents all the tools and command lines used to detect *cifA* and *cifB* genes and their predicted protein domains in bacterial genomes. The final *cifA* and *cifB* sequences, along with their predicted domains, used in this study are shared and can be used as queries to screen for *cif* presence in future genome analyses. 

#### 4. Determining *cif* flanking region
Here are the guidelines to identify *cif* flanking regions based on Prokka `.faa` protein files, RAGE and WO prophage databases, and BlastP analyses.

#### 5. Drawing of *cif* operons and genetic environment
This section proposes a detailed method to represent *cifA* and *cifB* structures and their flanking regions using `.gbk` files and homemade R scripts. Example files are provided to test the pipeline on a reduced dataset.

#### 6. Phylogenies
All the scripts used for phylogenetic analyses, from single-gene phylogenies to whole-genome phylogenies. Example files are provided to test the workflow on a reduced dataset.
