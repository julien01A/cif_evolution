# Phylogenies

## Phylogeny of one or a few genes

For each gene, operon, or concatenated sequence of interest targeted for phylogenetic analysis, an alignment file was required. All the alignment files were generated by aligning the desired sequences using `Clustal Omega v1.2.2` (doi: <https://doi.org/10.1002/pro.3290>) implemented on `UGENE v52.0` (<https://ugene.net/>, doi: <https://doi.org/10.1093/bioinformatics/bts091>), and then the positions with gaps '-' were removed. All the inital sequence files (before aligment and concatenation, if any) and the final alignment files were joined: `cifA_RNA_binding_like.faa`, `cifB_AAA.faa`, `cifB_PDDEXK_N_terminal.faa`, `cifB_PDDEXK_C_terminal.faa` and the final `cif_operon_CONCATENATED_ALIGNMENT.faa` for *cif* operon, `RAGE_5_Tra_genes.faa` and `RAGE_CONCATENANED_ALIGNMENT.faa` for RAGE, `WOrecomb.faa` and `WOrecomb_ALIGNMENT.faa` for WO prophage, `PDDEXK2.faa` and `PDDEXK2_ALIGNMENT.faa` for PDDEXK2 transposon.

Then, for each ALIGNMENT.faa files, substitution models were evaluated using `modeltest v0.1.7` (doi: <https://doi.org/10.1093/molbev/msz189>) to determine the most appropriate ML substitution model (based on the AICc criterion), followed by phylogenetic tree construction with `raxml-ng v1.1.0` (doi: <https://doi.org/10.1093/bioinformatics/btz305>):
```
modeltest-ng -i ALIGNMENT.faa -p 12 -T raxml -d aa # for example: FLU+G4

raxml-ng --all --msa ALIGNMENT.faa --model FLU+G4 --prefix Your_Tree-raxmlng --seed 5 --threads 4 --bs-trees 1000
raxml-ng --support --tree Your_Tree-raxmlng.raxml.bestTree --bs-trees 1000 --prefix Your_Tree-boot --threads 2
```

Finally, the phylogenetic tree was visualized and adapted using `figtree` (<https://github.com/rambaut/figtree/>) and `MEGA11` (<https://megasoftware.net/>)


## Whole genome phylogenies

We made the whole genome phylogenies of the Rickettsiaceae family and the *Rickettsiella* genus. The methodology used is very similar to that described for basic phylogenies, but it requires a few additional steps upstream for the preparation of the `ALIGNMENT.faa` file.

First, single-copy orthologs (SCOs) were identified using `OrthoFinder v2.5.4` (<https://github.com/davidemms/OrthoFinder>, doi: <https://doi.org/10.1186/s13059-019-1832-y>):
```
orthofinder -f ./OrthoFinder_genomes/ -t 4 -S blast 
```
with `OrthoFinder_genomes` a directory containing corresponding `.faa` files obtained from Prokka. 

For each SCO, sequences were individually aligned using `MAFFT v7.505` (<https://github.com/GSLBiotech/mafft>, doi: <https://doi.org/10.1093/molbev/mst010>):
```
for file in /Single_Copy_Orthologue_Sequences/*
do mafft "$file" > "$file"
done
```

For each SCO, ambigious hypervariable regions were removed using `trimAl v1.2rev59` (<https://github.com/inab/trimal>, doi: <https://doi.org/10.1093/bioinformatics/btp348>):
```
cp ./Single_Copy_Orthologue_Sequences/*_align.fasta ./Single_Copy_Orthologue_Sequences_trimal/
for file in /Single_Copy_Orthologue_Sequences_trimal/*
do trimal -in "$file" -out "$file" -fasta -gt 1 -cons 50
done
```

Then, all SCO sequences were concatenated using `Amas v1.0` (<https://github.com/marekborowiec/AMAS>, doi: <https://doi.org/10.7717/peerj.1660>) in a single file:
```
for file in /Single_Copy_Orthologue_Sequences_trimal/*
do awk '/^>/{print ">organism" ++i; next}{print}' < "$file" > "${file%_align.fasta}_rename.fasta"
done
cp ./Single_Copy_Orthologue_Sequences_trimal/*_rename.fasta ./Single_Copy_Orthologue_Sequences_AMAS/
AMAS.py concat -f fasta -d aa --in-files ./Single_Copy_Orthologue_Sequences_AMAS/*.fasta
```

The final file produced is used as ALIGNMENT.faa, which is submitted to `modeltest` and `raxml-ng` like the basic phylogenies (see above). We provide the final alignment file used for the Rickettsiaceae (`Whole_Rickettsiaceae_ALIGNMENT.faa`) and the Rickettsiella (`Whole_Rickettsiella_ALIGNMENT.faa`) whole genome phylogenies.
