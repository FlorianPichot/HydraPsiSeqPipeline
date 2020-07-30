# HydraPsiSeqPipeline
You will find here the different HydraPsiSeq pipelines used by the NGS Core Facility "EpiRNA-Seq" from Nancy, France (https://umsibslor.univ-lorraine.fr/en/facility/epitranscriptomics-sequencing-epirna-seq).

!!! Be aware !!! All the scripts has been used in a local computer, with relatives paths for the inputs and outputs. If you want to use these scripts on your own, please be careful of the inputs and directories names along the script.

Every manual input that you should modify will be indicated by this symbol: //!\\\

# Requirements

An Unix shell and the R software are needed to use this scripts.

The "Unix_scripts" needs Trimmomatic, Bowtie2 and Bedtools. You can naturally use your own trimmming and alignment tools if you prefer, but be sure that the trimming output is in FASTQ format, and the alignment output in SAM format. Bedtools is mandatory.

# Pipeline workflow

General workflow :

| Trimming (Unix script) |    -->      | Alignement (Unix script) |     -->     | NormUCount generation (R script) |    -->     | Analysis with list (R script) |
                                                                                                                                   
Transcriptome workflow :             

| Trimming (Unix script) |    -->      | Alignement (Unix script) |     -->    | NormUCount transcriptome (R script) |   -->     | Analysis with list (R script) |



Analysis with list (R script) requires a list of known pseudouridine sites or putative pseudouridylated candidates.
The file "Hs_rRNA_modList.csv" is a template of such list which regroups known modified sites in human rRNAs.



