# Drawing of *cif* operons and genetic environment

## 1. Drawing of CifA and CifB protein domains

First, an information matrix is a prerequisited for drawing cifA and cifB. This matrix must include the total gene length, the start and end positions of each domain, and positions of ORF-disturbing mutations. The cifA and cifB matrices constructed for this study are attached to this GitHub repository in `cifA.txt` and `cifB.txt`. In `cifA.txt`, 'Apopt' is the abbreviation for the Apoptosis regulator-like domain and 'RNA' stands for the RNA-binding-like domain. 'STOP' correspond to the ORF-disturbing mutation codons. Same for the `cifB.txt` where 'OTU' stands for the OTU-like cysteine protease, 'AAA' for the AAA-ATPases-like, 'IPPDE' for the PD-(D/E)XK nuclease N-terminal, 'IIPDDE' for the PD-(D/E)XK nuclease C-terminal, 'TPR' for the TPR repeats, 'DUB' for the Deubiquitinase DUB domain, 'TOXIN' for the Pore forming toxin TcdA/TcdB domain, 'DUF' for the DUF3491 domain, 'RTX' for the RTX toxin, 'SAL' for the salivary-gland toxin, 'ANK' for the Ankyrin repeat domains. and 'LATRO' for the Latrotoxin domain.

Then, these two matrices are imported into RStudio v4.3.1 (http://www.rstudio.com/):
```
CIFA <- read.table("cifA.txt",sep="\t",dec=",",header =T,row.names=NULL)
CIFB <- read.table("cifB.txt",sep="\t",dec=",",header =T,row.names=NULL)
```

The following libraries are loaded : ggplot2 v.3.4.4 (doi: <10.1007/978-3-319-24277-4>), cowplot v1.1.1 (doi: <10.32614/CRAN.package.cowplot>), and gridExtra v.2.3 (doi: <10.32614/CRAN.package.gridExtra>).
```
library(ggplot2)
library(cowplot)
library(gridExtra)
```

Here is the command to create a plot specifically for the different cifA ORF. Briefly, we first draw a line corresponding to the length of the cifA ORF for each symbiont, then we add the specific domains for each cifA ORF as well as the ORF-disturbing positions, and finally, we set the titles and theme for better representation:
```
plot_CIFA<-ggplot(CIFA) +
  geom_segment(aes(x = start, xend = end, y = fragment, yend = fragment), color = "black", size = 1) +
  geom_segment(aes(x = STOP10_START, xend = STOP10_END, y = fragment, yend = fragment), color = "black", size = 10) +
  geom_segment(aes(x = RNA_START, xend = RNA_END, y = fragment, yend = fragment), color = "purple", size = 10) +
  geom_segment(aes(x = Apopt_START, xend = Apopt_END, y = fragment, yend = fragment), color = "orange", size = 10) +
  geom_segment(aes(x = STOP1_START, xend = STOP1_END, y = fragment, yend = fragment), color = "black", size = 10) +
  geom_segment(aes(x = STOP2_START, xend = STOP2_END, y = fragment, yend = fragment), color = "black", size = 10) +
  geom_segment(aes(x = STOP3_START, xend = STOP3_END, y = fragment, yend = fragment), color = "black", size = 10) +
  geom_segment(aes(x = STOP4_START, xend = STOP4_END, y = fragment, yend = fragment), color = "black", size = 10) +
  geom_segment(aes(x = STOP5_START, xend = STOP5_END, y = fragment, yend = fragment), color = "black", size = 10) +
  geom_segment(aes(x = STOP6_START, xend = STOP6_END, y = fragment, yend = fragment), color = "black", size = 10) +
  geom_segment(aes(x = STOP7_START, xend = STOP7_END, y = fragment, yend = fragment), color = "black", size = 10) +
  geom_segment(aes(x = STOP8_START, xend = STOP8_END, y = fragment, yend = fragment), color = "black", size = 10) +
  geom_segment(aes(x = STOP9_START, xend = STOP9_END, y = fragment, yend = fragment), color = "black", size = 10) +
  labs(x = "Size in amino acids", y = NULL, title = "cifA-like") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line.x = element_line(color = "black", size = 0.5))
```

The cifB commands are similars:
```
plot_CIFB<-ggplot(CIFB) +
  geom_segment(aes(x = start, xend = end, y = fragment, yend = fragment), color = "black", size = 1) +
  geom_segment(aes(x = OTU_START, xend = OTU_END, y = fragment, yend = fragment), color = "#FF6600", size = 10) +
  geom_segment(aes(x = AAA_START, xend = AAA_END, y = fragment, yend = fragment), color = "green4", size = 10) +
  geom_segment(aes(x = IPDDE_START, xend = IPDDE_END, y = fragment, yend = fragment), color = "#CC0000", size = 10) +
  geom_segment(aes(x = IIPDDE_START, xend = IIPDDE_END, y = fragment, yend = fragment), color = "#CC0000", size = 10) +
  geom_segment(aes(x = ITPR_START, xend = ITPR_END, y = fragment, yend = fragment), color = "slateblue4", size = 10) +
  geom_segment(aes(x = IITPR_START, xend = IITPR_END, y = fragment, yend = fragment), color = "slateblue4", size = 10) +
  geom_segment(aes(x = DUB_START, xend = DUB_END, y = fragment, yend = fragment), color = "#0099FF", size = 10) +
  geom_segment(aes(x = TOXIN_START, xend = TOXIN_END, y = fragment, yend = fragment), color = "pink", size = 10) +
  geom_segment(aes(x = IDUF_START, xend = IDUF_END, y = fragment, yend = fragment), color = "grey", size = 10) +
  geom_segment(aes(x = IIDUF_START, xend = IIDUF_END, y = fragment, yend = fragment), color = "grey", size = 10) +
  geom_segment(aes(x = RTX_START, xend = RTX_END, y = fragment, yend = fragment), color = "blue", size = 10) +
  geom_segment(aes(x = SAL_START, xend = SAL_END, y = fragment, yend = fragment), color = "#003300", size = 10) +
  geom_segment(aes(x = IANK_START, xend = IANK_END, y = fragment, yend = fragment), color = "#33FFCC", size = 10) +
  geom_segment(aes(x = IIANK_START, xend = IIANK_END, y = fragment, yend = fragment), color = "#33FFCC", size = 10) +
  geom_segment(aes(x = IIIANK_START, xend = IIIANK_END, y = fragment, yend = fragment), color = "#33FFCC", size = 10) +
  geom_segment(aes(x = LATRO_START, xend = LATRO_END, y = fragment, yend = fragment), color = "deeppink3", size = 10) +
  geom_segment(aes(x = STOP1_START, xend = STOP1_END, y = fragment, yend = fragment), color = "black", size = 10) +
  geom_segment(aes(x = STOP2_START, xend = STOP2_END, y = fragment, yend = fragment), color = "black", size = 10) +
  geom_segment(aes(x = STOP3_START, xend = STOP3_END, y = fragment, yend = fragment), color = "black", size = 10) +
  geom_segment(aes(x = STOP4_START, xend = STOP4_END, y = fragment, yend = fragment), color = "black", size = 10) +
  geom_segment(aes(x = STOP5_START, xend = STOP5_END, y = fragment, yend = fragment), color = "black", size = 10) +
  geom_segment(aes(x = STOP6_START, xend = STOP6_END, y = fragment, yend = fragment), color = "black", size = 10) +
  geom_segment(aes(x = STOP7_START, xend = STOP7_END, y = fragment, yend = fragment), color = "black", size = 10) +
  geom_segment(aes(x = STOP8_START, xend = STOP8_END, y = fragment, yend = fragment), color = "black", size = 10) +
  geom_segment(aes(x = STOP9_START, xend = STOP9_END, y = fragment, yend = fragment), color = "black", size = 10) +
  geom_segment(aes(x = STOP10_START, xend = STOP10_END, y = fragment, yend = fragment), color = "black", size = 10) +
  geom_segment(aes(x = STOP11_START, xend = STOP11_END, y = fragment, yend = fragment), color = "black", size = 10) +
  geom_segment(aes(x = STOP12_START, xend = STOP12_END, y = fragment, yend = fragment), color = "black", size = 10) +
  labs(x = "Size in amino acids", y = NULL, title = "cifB-like") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line.x = element_line(color = "black", size = 0.5))
```

Finally, we combine the cifA and cifB of each symbiont on the same plot and add a scale to visualize the length of the two genes:
```
max_length <- max(max(CIFA$end), max(CIFB$end))
plot_CIFA <- plot_CIFA + 
  geom_vline(xintercept = seq(0, max_length, by = 500), color = "lightgray", linetype = "solid", size = 0.1) + 
  xlim(0, max_length)
plot_CIFB <- plot_CIFB +
  geom_vline(xintercept = seq(0, max_length, by = 500), color = "lightgray", linetype = "solid", size = 0.1) +
  xlim(0, max_length) +
  scale_y_discrete(labels = NULL)

combined_plot <- plot_grid(plot_CIFA, plot_CIFB, ncol = 2, align = "v")
print(combined_plot)
```

For the final plot, we can merge the cif operon phylogeny (see the directory part.5) and arrange *cifA-cifB* according to their phylogenetic placement.


## 2. Drawing of RAGE and WO prophage environments
First, using the Prokka annotation, database searches in RAGE and WO with OrthoFinder and blastP (see the directory part.3), we created a table grouping key information about syntenic gene environments bordering *cifA-cifB*. An example of this table, named `Environment.txt`, is provided for two typical *cifA-cifB* environments: RAGE and WO prophages. eg. (i) *Rickettsia hoogstraalii* strain CS carry *cif* operon in a RAGE module on contig29 (see `RhooCS.gbk` for full genome annotated raw data Prokka file, or `RhooCS_29.txt` for contig29 only), (ii) *Mesenet longicola* strain GL2 carry *cif* operon in a WO prophage or WO-like island on contig9 (see `MesenetGL2.gbk` for full genome annotated raw data Prokka file, or `MesenetGL2_9.txt` for contig9 only). 
The `Environment.txt` table indicates the distinct ORFs identified by Prokka ('Prokka_name') and their correspondance in the RAGE/WO databases ('Gene'). Regarding the two examples, contig29 of *R. hoogstraalii* CS contains 92 ORFs (only 56 are presented in the `Environment.txt` table); contig9 of *M. lonficola* GL2 hold 22 ORFs (all presented in `Environment.txt`).

We develop a specific R script capable of (i) reading the `RhooCS_29.txt` and `MesenetGL2_9.txt` files (derived from `RhooCS.gbk` and `MesenetGL2.gbk`),  extracting key information such as the number of ORFs, ORF lengths, distances between ORFs, and coding strand orientation, and (ii) reading the `Environment.txt` table to compare Prokka annotation names from `RhooCS_29.txt` and `MesenetGL2_9.txt` with their corresponding entries in the RAGE/WO databases. Finally, this script generates a genomic plot of the selected contigs (see `plot_examples.png`), highlighting *cif* operons and their neighboring RAGE or WO prophage environments.

This R script is as follows:

First, load the libraries and navigate to a working directory where your files are located.
```
library(openxlsx)
library(dplyr)
library(stringr)
library(ggplot2)

directory_path <- "indicate your working directory here"
```

Then, the script processes each `.txt` files, extracts key information for each ORF, and stores it in a synthenic matrix.
```
files <- list.files(path = directory_path, pattern = "\\.txt$", full.names = TRUE)
all_gene_data <- list() 

for (file_path in files) {
  lines <- readLines(file_path)
  fragment <- sub("LOCUS\\s+([A-Za-z0-9_]+)\\s+.*", "\\1", lines[grep("^LOCUS", lines)])
  fragment_size <- as.numeric(sub(".*\\s(\\d+)\\s+bp.*", "\\1", lines[grep("^LOCUS", lines)]))
  fragment_position <- paste0("1..", fragment_size)
  gene_data <- list()
  record_gene <- FALSE
  for (i in seq_along(lines)) {
    line <- lines[i]
    if (grepl("^\\s+gene\\s", line) || grepl("^\\s+CDS\\s", line)) {
      if (grepl("^\\s+CDS\\s", line)) {
        record_gene <- TRUE
        position <- sub(".*?(\\d+..\\d+|complement\\(\\d+..\\d+\\)).*", "\\1", line)
        sense <- ifelse(grepl("complement", line), "anti_sens", "standard")
        pos_clean <- gsub("complement\\(|\\)", "", position)
        start_end <- as.numeric(unlist(strsplit(pos_clean, "\\.\\.")))
        gene_size <- abs(start_end[2] - start_end[1]) + 1
      }
    }
    if (record_gene && grepl("/locus_tag=", line)) {
      gene <- sub('.*/locus_tag="([^"]+)".*', "\\1", line)
      record_gene <- FALSE
      gene_data[[length(gene_data) + 1]] <- c(
        fragment,
        fragment_size,
        gene,
        position,
        gene_size,
        sense
      )
    }
  }
  gene_matrix <- as.data.frame(do.call(rbind, gene_data))
  colnames(gene_matrix) <- c("fragment", "fragment_size", "gene", "gene_position", "gene_size", "gene_sens")
  gene_matrix$fragment_size <- as.numeric(gene_matrix$fragment_size)
  gene_matrix$gene_size <- as.numeric(gene_matrix$gene_size)
  file_name <- basename(file_path)
  gene_matrix$file <- file_name 
  all_gene_data[[file_name]] <- gene_matrix
}

synthetic_gene_matrix <- do.call(rbind, all_gene_data)
synthetic_gene_matrix <- synthetic_gene_matrix[, c("fragment", "fragment_size", "gene", "gene_position", "gene_size", "gene_sens")]
print(synthetic_gene_matrix)
```

Then, for each ORF in the matrix, the script adds its corresponding RAGE gene or WO prophage by comparing it with the information in the 'environment.txt' table. A new synthenic matrix with RAGE or WO gene informations is provided.
```
environment_table <- read.delim(".../environment.txt", header = TRUE, sep = "\t")
environment_table <- environment_table[, c("Gene", "Prokka_name")]
environment_table$Prokka_name <- sub(" .*", "", environment_table$Prokka_name)
environment_table$gene <- environment_table$Prokka_name
environment_table$Fonction <- sub(" .*", "", environment_table$Gene)
environment_table <- environment_table[, c("gene", "Fonction")]
synthetic_gene_matrix <- merge(synthetic_gene_matrix, environment_table, by.x = "gene", by.y = "gene", all.x = TRUE)
synthetic_gene_matrix_merged <- synthetic_gene_matrix[, c("fragment", "fragment_size", "gene", "gene_position", "gene_size", "gene_sens", "Fonction")]
synthetic_gene_matrix_merged$Fonction[synthetic_gene_matrix_merged$Fonction == "" | is.na(synthetic_gene_matrix_merged$Fonction)] <- "other_or_unspecific"
print(synthetic_gene_matrix_merged)
```

Prepare different parameters to create a graphical representation from the synthetic_gene_matrix_merged: Indicate and apply the colors of your choice to the genes and adjust the genes on each contig.
```
color_map <- c(
  "TraA" = "#009933",
  "TraL" = "#009933",
  "TrakN" = "#009933",
  "TrakC" = "#009933",
  "TraV" = "#009933",
  "TraE" = "#009933",
  "TraB" = "#009933",
  "TraF" = "#009933",
  "TraH" = "#009933",
  "TraN" = "#009933",
  "TrbC" = "#009933",
  "TraU" = "#009933",
  "TraW" = "#009933",
  "TraG" = "#009933",
  "TraC" = "#009933", 
  "TraD_Flike" = "#009933", 
  "TraAI" = "#99CC33",
  "TraD_Ti_like" = "#99CC33",
  "Integrase" = "#FF0000",
  "Transposase" = "#FFFF00",
  "SpoT_hydrolase" = "#FF9900",
  "Dna_methyl" = "#FF9900",
  "DnaI" = "#FF9900",
  "TPR" = "#FF9900",
  "HisKinase" = "#0033FF",
  "cifA" = "#FF99FF",
  "cifB" = "#FF99FF",
  "Clp" = "#333399",
  "CK" = "#669999",
  "Ribosomal" = "#000000",
  "other_or_unspecific" = "#CCCCCC",
  "ToxinTAT" = "#FF9999",
  "Phage_like" = "#996633",
  "ANK" = "#00FFFF",
  "Reverse" = "#FF9966",
  "ABC" = "#660000",
  "membrane" = "#660000",
  "Tail" = "#FF9999",
  "Tail_Fiber" = "#FFCCCC",
  "Patatin" = "#996633",
  "MutL" = "#000066"
)

synthetic_gene_matrix_merged <- synthetic_gene_matrix_merged %>%
  mutate(color = color_map[as.character(Fonction)])
create_gene_data <- function(fragment_name, gene_matrix) {
  gene_matrix %>%
    filter(fragment == fragment_name) %>%
    mutate(
      start_position = as.numeric(sapply(strsplit(as.character(gene_position), "\\.\\."), `[`, 1)),
      end_position = as.numeric(sapply(strsplit(as.character(gene_position), "\\.\\."), `[`, 2)),
      y_position = which(unique(gene_matrix$fragment) == fragment_name),
      # Calculer la position y pour le texte en fonction de gene_sens
      text_y_position = ifelse(gene_sens == "standard", y_position + 0.3, y_position - 0.3)
    )
}
gene_data_list <- lapply(unique(synthetic_gene_matrix_merged$fragment), 
                         create_gene_data, 
                         gene_matrix = synthetic_gene_matrix_merged)

all_gene_data <- bind_rows(gene_data_list)

fragment_positions <- synthetic_gene_matrix_merged %>%
  select(fragment, fragment_size) %>%
  distinct() %>%
  mutate(y_position = row_number())
```

A final plot can then be created with ggplot (see `plot_examples.png`). The parameters can be adjusted to your preference. For example, here we have specifically indicated the *cifA* and *cifB* genes, as well as the *Tra* genes, which primarily make up the RAGE.
```
ggplot() +
  geom_segment(data = fragment_positions, 
               aes(x = 0, xend = fragment_size, 
                   y = y_position, yend = y_position), 
               color = "black", size = 1.5) +
  geom_rect(data = all_gene_data,
            aes(xmin = start_position, xmax = end_position, 
                ymin = ifelse(gene_sens == "standard", y_position, y_position - 0.175), 
                ymax = ifelse(gene_sens == "standard", y_position + 0.175, y_position), 
                fill = Fonction),
            color = "black", alpha = 0.7) + 
  geom_text(data = all_gene_data %>% 
              filter((grepl("^Tra", Fonction) & !grepl("Transposase", Fonction)) | 
                       grepl("^TrbC", Fonction) | 
                       grepl("^cif", Fonction)), 
            aes(x = (start_position + end_position) / 2, 
                y = ifelse(gene_sens == "standard", y_position + 0.2, y_position - 0.2),  
                label = Fonction),  
            angle = 45, hjust = 0.5, size = 3.5) +
   geom_text(data = fragment_positions, 
            aes(x = -500, y = y_position, label = fragment),  
            hjust = 1, size = 3.5) +
  theme_minimal() +
  labs(x = "Taille du Fragment", 
       y = "Fragments",
       fill = "Type of genes")) +
  theme(axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(), 
        plot.title = element_text(hjust = 0.5)) +
  xlim(-600, max(fragment_positions$fragment_size)) + 
  scale_fill_manual(values = color_map)
```

We then manually adjust the plots to best fit the figures.
