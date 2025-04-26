# Genome annotation

We get the `.fasta` of the 786 genomes (786 $ech) on a repertory. Then, each genome were annotaded using `Prokka (v1.14.6)` (doi: <https://doi.org/10.1093/bioinformatics/btu153>) with the following command:

```
prokka $ech-genome.fasta --locustag $ech --prefix $ech --outdir Prokka-$ech --rfam --compliant --cpus 6
```

Note that for the *Spiroplasma* genomes (Mollicutes), we used genetic code 4 instead of the default code 11. For *Hodgkinia*, each genome was annotated twice: once with code 4 and once with code 11. Here is the command for code 4:

```
prokka --gcode 4 $ech-genome.fasta --locustag $ech --prefix $ech --outdir Prokka-$ech --rfam --compliant --cpus 6
```

For each genome, Prokka give the list of all the protein annotated in a `.faa` file, which are necessary for Orthofinder analysis (Protein screening against specific databases, see part.4, or whole genome phylogenies, see part.6). The `.gbk` file will also be used to detect the position of *cif* genes in the bacterial genome and trace their neighboring genetic environments.


# Determining the Average Nucleotide Identity (ANI) for each genus

To determine the ANI of genomes for each genus, we perform an "all-vs-all" genome comparison using `fastANI (v.1.33)` ([GitHub link](https://github.com/ParBLiSS/FastANI)). For each genus, all the genomes in `.fna` format were placed in the same folder, and we run the following command to generate a list of genome files:

```
ls /Genus1_genomes/*.fna > Genus1_genomes_list.txt
```

Then, we run the comparison using fastANI:

```
fastANI --ql Genus1_genomes_list.txt --rl Genus1_genomes_list.txt -o Genus1_genomes_result_ani.txt -t 4
```

We removed lines comparing the same genomes (100% identity).
