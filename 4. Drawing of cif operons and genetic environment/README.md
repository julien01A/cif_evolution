# Drawing of *cif* operons and genetic environment

## 1. Drawing of CifA and CifB protein domains

First, an information matrix is a prerequisited for cifA and cifB. This matrix must include the total gene length, the start and end positions of each domain, and positions of ORF-disturbing mutations. The cifA and cifB matrices constructed for this study are attached to this GitHub repository in `cifA.txt` and `cifB.txt`. In `cifA.txt`, 'Apopt' is the abbreviation for the Apoptosis regulator-like domain and 'RNA' stands for the RNA-binding-like domain. 'STOP' correspond to the ORF-disturbing mutation codons. Same for the `cifB.txt` where 'OTU' stands for the OTU-like cysteine protease, 'AAA' for the AAA-ATPases-like, 'IPPDE' for the PD-(D/E)XK nuclease N-terminal, 'IIPDDE' for the PD-(D/E)XK nuclease C-terminal, 'TPR' for the TPR repeats, 'DUB' for the Deubiquitinase DUB domain, 'TOXIN' for the Pore forming toxin TcdA/TcdB domain, 'DUF' for the DUF3491 domain, 'RTX' for the RTX toxin, 'SAL' for the salivary-gland toxin, 'ANK' for the Ankyrin repeat domains. and 'LATRO' for the Latrotoxin domain.

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

For the final plot, we can merge the cif operon phylogeny (see the directory part.5) and arrange cifA-cifB according to their phylogenetic placement.


## 2. Drawing of RAGE, WO prophage and SMGE-rich environments
First, using the Prokka annotation, database searches in RAGE and WO with OrthoFinder and blastP (see the directory part.4), we created a table grouping key information about syntenic gene environments bordering cifA-cifB. An example of this table, named `Environment.txt`, is provided for three typical cifA-cifB environments: RAGE, WO prophages, and SMGE-rich islands.
