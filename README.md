# HydraPsiSeqPipeline
You will find here the different RiboMethSeq pipelines used by the NGS Core Facility "EpiRNA-Seq" from Nancy, France (https://umsibslor.univ-lorraine.fr/en/facility/epitranscriptomics-sequencing-epirna-seq).

!!! Be aware !!! All the scripts has been used in a local computer, with relatives paths for the inputs and outputs. If you want to use these scripts on your own, please be careful of the inputs and directories names along the script.

Every manual input that you should modify will be indicated by this symbol: //!\\\

# Pipeline workflow

General workflow :

| Trimming (Unix script) |    -->      | Alignement (Unix script) |     -->     | NormUCount generation (R script) |    -->     | Analysis with list (R script ) |
                                                                                                                                   
Transcriptome workflow :             

| Trimming (Unix script) |    -->      | Alignement (Unix script) |     -->    | NormUCount transcriptome (R script) |   -->     | Analysis with list (R script ) |
