#The definition of a high-confidence orthologue gene set from separate de novo RNA-seq assemblies

This is a step-by-step guide to the bioinformatics analyses that generated a high-confidence orthologue gene set and input for coalescence analyses in the following study:

Nürnberger, Lohse, Fijarczyk, Szymura and Blaxter. Para-allopatry in hybridising fire-bellied toads (Bombina bombina and B. variegata): inference from transcriptome-wide coalescence analyses. [bioRxiv manuscript](http://biorxiv.org/content/early/2015/10/28/030056)

This repository is provided to document the analyses. The pipeline is highly customised and would need extensive modification in order to be applied to other datasets. At present, there are no plans to turn this pipeline into a software tool. The rationale for the analyses is presented in the manuscript. 


####Input data:
Illumina RNAseq reads, assembled with Trinity<sup>1</sup> (fasta file)\n
Roche 454 RNA-seq reads, assembled Newbler/Mira<sup> 2</sup>/CAP3<sup>3</sup> (fasta file)
 MySQL, 
####Additional Dependencies
BLAST, EMBOSS revseq and water, Clustal-omega<sup>4</sup>, Gblocks<sup>5</sup>, Muscle<sup>6</sup>, MySQL, OrthoMCL<sup>7</sup>,  PAL2NAL<sup>8</sup>, RSEM<sup>9</sup>
