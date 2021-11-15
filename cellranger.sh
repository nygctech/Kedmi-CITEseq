#please note that for this analysis, we used cellranger v5.0.0
#reference here is derived from the cellranger mm10 2020A reference with three additional sequences appended to the reference

#lane 1 GEX data
cellranger count --sample=cDNA1A-AAGACGGA,cDNA1C-CGAGGCTC,cDNA1G-GTCCTTCT,cDNA1T-TCTTAAAG --transcriptome=mouse_Ranit_cre_GFP_tdTomato --id=cDNA1_GEX --fastqs=./fastq/

#lane 2 GEX data
cellranger count --sample=cDNA2A-AGCGAAAG,cDNA2C-CCGTGTGA,cDNA2G-GAAACCCT,cDNA2T-TTTCTGTC --transcriptome=mouse_Ranit_cre_GFP_tdTomato --id=cDNA2_GEX --fastqs=./fastq/

