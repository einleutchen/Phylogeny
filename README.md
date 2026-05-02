# Phylogeny
A project from Comparative Genomics course
Species belong to the same genus, like Clostridium, can be much diverged in
genomes, however, proteins might be conserved. Therefore, to select homologous
prote, I set the threshold for percentage of identity at 70% and query coverage
90% to obtain diverged-but-whole sequences. The cutoffs for selection of
homologous genes, that are 70% identity and 50% coverage. When targetting higher
cutoffs for genes, I could not obtain any set of sequences for further analysis,
plus to avoid collecting sets of fragmented sequences (query coverage < 40) I
think these cutoffs are reasonable. The E-value for protein is 1e-05 while for
genes 1e-06 was applied as two nucleotide sequence has a higher probability to
match, thus, stricter threshold is essential to filter out random matches.
To choose the best genome among 11 Clostridium genomes for reference, I
investigates available metrics that come along with the genome itself:
1. complete genome with 1 contigous sequence
2. checkM anaylsis shows high completeness percentage
3. contamination is much low
4. total length covers most of the other genomes
Eventually, I define Clostridium diolis as reference genome as it meets all the
criteria above.
IQTree2 iterates over different models, such as Jukes-Cantor, LG, WAG, Blosum62,
Dayhoff with different combinations. The ModelFinder will then return the
optimal model for evolution according to BIC (Bayesian Information Criterion by
default). IQTree2 eventually returns one of the files with the best tree that
has the highest likelihood of representing data (Maximum Likelihood).
The evolutionary distances in 5 out of 6 my protein sets (including concatenated
one) are calculated by LG+G4 while one is estimated by LG+I+G4. While the
evolutionary changes in nucleotide sequences are estimated with GTR+F+G4 (5
sets), GTR+F+I+G4 (1 set).
Annotation retrieval from NCBI Blast of each representative sequence
Tree
Protein
Nucleotide
1
Glycine-tRNA ligase
Elongation factor Tu
2
RNA polymerase
Elongation factor Tu
sporulattion factor
3
50S ribosomal protein L11 DNA polymerase subunit B
4
50S ribosomal protein L1 x
5
50S ribosomal L7/L12
x
Concat
Glycine-tRNA ligase (cov: DNA polymerase (cov: 44%,
39%, ident: 97%)
ident: 98%
Manual batch-blasting of all sequences in the set returns the same annotations.
This confirms the reliability of representative sequence and the effectiveness
of the pipeline in finding homologous sequences.
Protein and gene trees in general shows C. aceticum and C. formicaceticum as
distantly relatated sister groups to the others which form a seperate group.
This demonstrates that the selected genes and proteins are conserved among thesespecies and probably theses two species live in another environment compared to
the others. C. butyricum, C. beijerinckii and C. diolis consistently co-appear
in a cluster, indicating their close relatedness. C. beijerinckii and C. diolis
are indeed sister taxa while C. butyricum exhibited longer branch lengths/higher
rate of substitution. so that it is distinguished from the other two though they
share a common ancestor. Addtionally, the branch length in gene tree Similarly,
C. drakei and C. scatologenes are also a sister group. On the other hands, C.
cochlearium, C. perfringens, C. septicum and C. sporogenes shows less consistent
in grouping across different trees. C. perfringens and C. septicum are grouped
together in trees of RNA polymerase sporulation factor and 50S ribosomal
protein. In other cases, they form in a paraphyletic group. In majority of gene
trees, C. cochlearium does not form an apparent sister group with other species.
This means that the selected genes and proteins in C. cochlearium have much more
evolutionary changes. Protein trees, 3 out of 6 trees (including concatenated
one), reveal the opposite in which C. cochlearium occasionally shares the recent
ancestor with C. drakei and C. scatologenes.
C. sporogenes also groups with C. drakei and C. scarologenes and happens to be a
basal taxon in some trees, probably some proteins/genes in C. sporogenes retain
ancestral characteristics.
The concatenated alignment based tree captures the close relatedness of C.
beijerinckii and C. diolis as well as C. drakei and C. scatologenes. In
individual gene or protein tree, the relationship between C. perfringens and C.
sporogene remains unclear. However, in concatenated trees for both genes and
proteins they are immediate descendants of a common ancestor, indicating that
the more information (multiple genes or proteins) that is fed to tree-building
model, the higher possibility of obtaining the more accurate tree topology. In
some cases, the branch length of individual gene or protein trees are shorter
than concatenated tree. This can be explained by the difference in rate of
genetic change across genes or proteins. The branch length of concatenated tree
does not account for genetic change of single gene or protein, but the average
amount of change across the number of genes or proteins, that can also mean the
list of genes or proteins with different evolutionary history or mutation rate
can make the branch length either longer or shorter. A single gene or protein
trees reflect the true evolutionary history of a single locus, while
concatenated tree neutralizes this. Depending on the research questions and
concerns, one can choose one method over another.
