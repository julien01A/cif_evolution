# Genome annotation

We get the `.fasta` of the 786 genomes (786 $ech) on a repertory. Then, each genome were annotaded using Prokka (Seemann, 2014, <doi:10.1093/bioinformatics/btu153>) with the following command:

```
prokka $ech-genome.fasta --locustag $ech --prefix $ech --outdir Prokka-$ech --rfam --compliant --cpus 6
```


# Drawing of *cif* operons and genetic environment

## 1. Drawing of CifA and CifB protein domains

First, an information matrix is a prerequisited for cifA and cifB. This matrix must include the total gene length, the start and end positions of each domain, and positions of ORF-disturbing mutations. The cifA and cifB matrices constructed for this study are attached to this GitHub repository in `cifA.txt` and `cifB.txt`. In `cifA.txt`, 'Apopt' is the abbreviation for the Apoptosis regulator-like domain and 'RNA' stands for the RNA-binding-like domain. Same for the `cifB.txt` where 'OTU' stands for the OTU-like cysteine protease, 'AAA' for the AAA-ATPases-like, 'IPPDE' for the PD-(D/E)XK nuclease N-terminal, 'IIPDDE' for the PD-(D/E)XK nuclease C-terminal, 'TPR' for the TPR repeats, 'DUB' for the Deubiquitinase DUB domain, 'TOXIN' for the Pore forming toxin TcdA/TcdB domain, 'DUF' for the DUF3491 domain, 'RTX' for the RTX toxin, 'SAL' for the salivary-gland toxin, 'ANK' for the Ankyrin repeat domains. and 'LATRO' for the Latrotoxin domain.

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
  # Personnalisation des axes et du titre
  labs(x = "Size in amino acids", y = NULL, title = "cifA-like") +
  # Spécifier un fond blanc
  theme(panel.background = element_rect(fill = "white"),
        axis.line.x = element_line(color = "black", size = 0.5))
