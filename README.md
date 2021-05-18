# On the identifiability of the isoform deconvolution problem: application to select the proper fragment length in an RNAseq library

**Abstract** <br>
**Motivation**: Isoform deconvolution is an NP-hard problem. The accuracy of the proposed solutions are far from perfect. At present, it is not known if gene structure and isoform concentration can be uniquely inferred given paired-end reads, and there is no objective method to select the fragment length to improve the number of identifi-able genes. Different pieces of evidence suggest that the optimal fragment length is gene-dependent, stressing the need for a method that selects the fragment length according to a reasonable trade-off across all the genes in the whole genome. <br>
**Results**: A gene is considered to be identifiable if it is possible to get both the structure and concentration of its transcripts univocally. Here, we present a method to state the identifiability of this deconvolution problem. Assum-ing a given transcriptome and that the coverage is sufficient to interrogate all junction reads of the transcripts, this method states whether or not a gene is identifiable given the read length and fragment length distribution. <br>
Applying this method using different read and fragment length combinations, the optimal average fragment length for the human transcriptome is around 400-600nt for coding genes and 150-200nt for long non-coding RNAs. The optimal read length is the largest one that fits in the fragment length. It is also discussed the potential profit of combining several libraries to reconstruct the transcriptome. Combining two libraries of very different fragment lengths results in a significant improvement in gene identifiability. 

<br><br>

**Code** <br>

* C1orf174_gene_example.R: main code to replicate the example depicted in the associated article.

* gencode24_chr22.R: main code to apply the method over the chr22 of the GenCode24 reference Transcriptome.

* twolibraries_simulation.R: Code to replicate the bootstrap analysis performed in the associate article. The aim of this analysis is to discuss the potential profit of combining several libraries to reconstruct the transcriptome.





