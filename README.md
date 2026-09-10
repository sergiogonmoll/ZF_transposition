# ZF_transposition
The genomic data for this chapter can be found in the European Nucleotide Archive (ENA) under accession PRJNA1152924. These files were aligned to the zebra finch reference genome "GCA_003957565.4_bTaeGut1.4.pri_genomic.fna" using the "BWA_script.sh". Aligned BAM files were then used to run RetroSeq.

RetroSeq was run calling "retroseq.pl" and using a custom TE library for the zebra finch genome "ZF_reduced.fasta" on the aligned BAM files. These jobs were initiated by "ZF_array_script.sh". Comparison between parents and offspring are coded in "RetroSeq_output.R", filtering out TEs that are reference TEs and are within 300 base pairs of parental TEs to be defined as private TEs. This parent-offspring comparisons were run with the script "ZF_pedigree_RS.sh" by piping the data pedigree into the script. The pedigree is found in "ZF_pedigree_SRR.xlsx". The output of this batch script is the main data for the analysis found in "TEs_ZF.csv". Metadata used for the analysis is found in "TEs_ZF_metadata.csv".

Additionally, the TE landscape and divergence data for the reference genome and the chromosomal classifications are found in the files "TE_landscape.csv", "TE_landscape_macro.csv", "TE_landscape_micro.csv" and "TE_landscape_W.csv".

The analysis for the chapter and the figures presented within it are produced using the script "Ch4_analysis.R". In this script, TEs outside of the reliable mapapble areas from the genome are excluded from the analysis using the "mappable.bed.gz" BED file.

A complete repository of the PhD thesis has been created elsewhere since file size limitations make it hard to upload all data, and can be found here: https://drive.google.com/file/d/1W5KDuDYTOoX9cbLOTBuyc-lM4j5ymqMH/view?usp=sharing

For any additional questions: sergiogonmoll[at]gmail.com or h.l.dugdale[at]rug.nl
