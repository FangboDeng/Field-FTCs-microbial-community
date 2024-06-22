setwd('D:/冻融/冻融田间实验/FTC field 稿/To STOEN/修回/Rcode')
##bacterial communites analysis as an example
	library(edgeR)
    library(indicspecies)
    library(tidyr)
    library(dplyr)
    library(stringr)

#####FTC 
     OTUs <- read.csv('OTUs.csv', stringsAsFactors = FALSE, row.names=1)
     OTUs <- OTUs[rowSums(OTUs) > 0,]

     keep_OTUs <- which(rowSums(OTUs >= 2) >= 16) 
     OTU_final <- OTUs[keep_OTUs,]
	 
	 taxon_otu <- read.csv('OTUs_taxonomy.csv',stringsAsFactors = FALSE, row.names=1) 
	 tax_otu <- taxon_otu[rownames(OTU_final),]

	 design <- read.csv('design_name.csv') %>% select(-Treatment,-R)
	 edgeR_otu <- DGEList(counts=OTU_final, 
                          group=design$FTC,
                          genes=tax_otu)
						  
     edgeR_otu <- calcNormFactors(edgeR_otu)

	 ## Get TMM normalized counts expressed as relative abundance counts per million 
     otu_norm <- cpm(edgeR_otu, normalized.lib.sizes=T, log=F)
	 
     ##### indicator species analysis #####
      indic <- as.data.frame(t(otu_norm))
      indic_groups <- design$FTC
      length(unique(indic_groups))
      set.seed(168168)
      indicatorsp <- multipatt(indic,indic_groups,func = "r.g",control=how(nperm=9999))
      summary(indicatorsp,alpha=1,indvalcomp=T)
	  summary(indicatorsp)
      indic_df <- indicatorsp$sign
      write.table(indic_df,paste0("indic_df.txt"),sep="\t",quote=F)
	  
     ## Import data frame of indicator species to save time
      indic_df <- read.table("indic_df.txt", header=T, sep="\t")
	 
	BF<- as.matrix(indic_df[which(indic_df$s.BF == 1 & indic_df$p.value < 0.05),])
	DF <- as.matrix(indic_df[which(indic_df$s.DF == 1 & indic_df$p.value < 0.05),])
	TF <- as.matrix(indic_df[which(indic_df$s.TF == 1 & indic_df$p.value < 0.05),])
	     

	indic_sign <- rbind(BF,DF,TF)
	colnames(indic_sign)[1:3] <-c("BF","DF","TF")
	
	## Range of correlation coefficients
	range(indic_sign[,"stat"])

    ## Total number of indicator OTUS
    length(unique(rownames(indic_sign)))

    ## Proportion of soil bacteria OTUs responding to treatment
    length(unique(rownames(indic_sign)))/nrow(indic_sign)

    ## Proportion 
    final_relative <- t(t(OTU_final)/colSums(OTU_final)) * 100
    sum(colSums(final_relative[unique(rownames(final_relative)),]))/sum(colSums(final_relative))

     write.table(indic_sign,paste0("indic_sign.txt"),sep="\t",quote=F)
    View <- filter(as.data.frame(indic_sign), index=="1")
	length(rownames(View))

	 tax_indic_sign <- taxon_otu[rownames(indic_sign),]
	 t <- as.data.frame(tax_indic_sign)
	 t_indic <-as.data.frame(indic_sign)
	 indi_otu_taxon <-  bind_cols(t,t_indic)
     write.csv(indi_otu_taxon,file="indi_otu_taxon.csv")

	 indic_str <- indicatorsp$str
	 indi_otu_taxon_str <- indic_str[rownames(indic_sign),]
     strt<- as.data.frame(indi_otu_taxon_str)
	 taxon<- as.data.frame(indi_otu_taxon)
	 Indi_otu_taxon_str <-  bind_cols(strt,taxon)
     write.csv(Indi_otu_taxon_str,file="Indi_otu_taxon_str.csv")
    
#####Autumn 
     OTUs_BF <- read.csv('BF_OTUs.csv', stringsAsFactors = FALSE, row.names=1)
     OTUs_BF <- OTUs_BF[rowSums(OTUs_BF) > 0,]

     keep_OTUs_BF <- which(rowSums(OTUs_BF >= 2) >= 4) 
     OTU_final_BF <- OTUs_BF[keep_OTUs_BF,]
	 taxon_otu <- read.csv('OTUs_taxonomy.csv',stringsAsFactors = FALSE, row.names=1) 
	 tax_otu_BF <- taxon_otu[rownames(OTU_final_BF),]

	 design_BF <- read.csv('design_name.csv') %>% filter(FTC=="BF")
	 edgeR_otu_BF <- DGEList(counts=OTU_final_BF, 
                          group=design_BF$Treatment,
                          genes=tax_otu_BF)
						  
     edgeR_otu_BF <- calcNormFactors(edgeR_otu_BF)

	 ## Get TMM normalized counts expressed as relative abundance counts per million 
     otu_norm_BF <- cpm(edgeR_otu_BF, normalized.lib.sizes=T, log=F)
	 
     ##### indicator species analysis #####
      indic_BF <- as.data.frame(t(otu_norm_BF))
      indic_groups_BF <- design_BF$Treatment
      length(unique(indic_groups_BF))
      set.seed(168168)
      indicatorsp_BF <- multipatt(indic_BF,indic_groups_BF,func = "r.g",control=how(nperm=9999))
      summary(indicatorsp_BF,alpha=1,indvalcomp=T)
	  summary(indicatorsp_BF)
      BF_indic_df <- indicatorsp_BF$sign
      write.table(BF_indic_df,paste0("BF_indic_df.txt"),sep="\t",quote=F)
	  
     ## Import data frame of indicator species to save time
      indic_df_BF <- read.table("BF_indic_df.txt", header=T, sep="\t")
	 
	S0F1_BF <- as.matrix(indic_df_BF[which(indic_df_BF$s.S0F135 == 1 & indic_df_BF$p.value < 0.05),])
	S0F2_BF <- as.matrix(indic_df_BF[which(indic_df_BF$s.S0F240 == 1 & indic_df_BF$p.value < 0.05),])
	S1F1_BF <- as.matrix(indic_df_BF[which(indic_df_BF$s.S100F135 == 1 & indic_df_BF$p.value < 0.05),])
	S1F2_BF <- as.matrix(indic_df_BF[which(indic_df_BF$s.S100F240 == 1 & indic_df_BF$p.value < 0.05),])
	
	indic_sign_BF <- rbind(S0F1_BF,S0F2_BF,S1F1_BF,S1F2_BF)
	colnames(indic_sign_BF)[1:4] <-c("S0F135","S0F240","S100F135","S100F240")
	## Range of correlation coefficients
	range(indic_sign_BF[,"stat"])

    ## Total number of indicator OTUS
    length(unique(rownames(indic_sign_BF)))

    ## Proportion of soil bacteria OTUs responding to treatment
    length(unique(rownames(indic_sign_BF)))/nrow(indic_sign_BF)

    ## Proportion 
    BF_final_relative <- t(t(OTU_final_BF)/colSums(OTU_final_BF)) * 100#求每一个OTU在整个处理中的百分比
    sum(colSums(BF_final_relative[unique(rownames(BF_final_relative)),]))/sum(colSums(BF_final_relative))

    write.table(indic_sign_BF,paste0("indic_sign_BF.txt"),sep="\t",quote=F)
    
	 tax_indic_sign_BF <- taxon_otu[rownames(indic_sign_BF),]
	 t_BF <- as.data.frame(tax_indic_sign_BF)
	 t_indic_BF <-as.data.frame(indic_sign_BF)
	 BF_indi_otu_taxon <-  bind_cols(t_BF,t_indic_BF)
     write.csv(BF_indi_otu_taxon,file="BF_indi_otu_taxon.csv")

	 BF_indic_str <- indicatorsp_BF$str
	 BF_indi_otu_taxon_str <- BF_indic_str[rownames(indic_sign_BF),]
     str_BF <- as.data.frame(BF_indi_otu_taxon_str)
	 taxon_BF <- as.data.frame(BF_indi_otu_taxon)
	 BF_indi_otu_taxon_str <-  bind_cols(str_BF,taxon_BF)
     write.csv(BF_indi_otu_taxon_str,file="BF_indi_otu_taxon_str.csv")
	 
#####Winter 
     
     OTUs_DF <- read.csv('DF_OTUs.csv', stringsAsFactors = FALSE, row.names=1)
     OTUs_DF <- OTUs_DF[rowSums(OTUs_DF) > 0,]

     keep_OTUs_DF <- which(rowSums(OTUs_DF >= 2) >= 4) 
     OTU_final_DF <- OTUs_DF[keep_OTUs_DF,]
	 
	 taxon_otu <- read.csv('OTUs_taxonomy.csv',stringsAsFactors = FALSE, row.names=1) 
	 tax_otu_DF <- taxon_otu[rownames(OTU_final_DF),]

	 design_DF <- read.csv('design_name.csv') %>% filter(FTC=="DF")
	 edgeR_otu_DF <- DGEList(counts=OTU_final_DF, 
                          group=design_DF$Treatment,
                          genes=tax_otu_DF)
						  
     edgeR_otu_DF <- calcNormFactors(edgeR_otu_DF)

	 ## Get TMM normalized counts expressed as relative abundance counts per million 
     otu_norm_DF <- cpm(edgeR_otu_DF, normalized.lib.sizes=T, log=F)
	 
     ##### indicator species analysis #####
      indic_DF <- as.data.frame(t(otu_norm_DF))
      indic_groups_DF <- design_DF$Treatment
      length(unique(indic_groups_DF))
      set.seed(168168)
      indicatorsp_DF <- multipatt(indic_DF,indic_groups_DF,func = "r.g",control=how(nperm=9999))
      summary(indicatorsp_DF,alpha=1,indvalcomp=T)
	  summary(indicatorsp_DF)
      DF_indic_df <- indicatorsp_DF$sign
      write.table(DF_indic_df,paste0("DF_indic_df.txt"),sep="\t",quote=F)

	  
     ## Import data frame of indicator species to save time
      indic_df_DF <- read.table("DF_indic_df.txt", header=T, sep="\t")
	 
	S0F1_DF <- as.matrix(indic_df_DF[which(indic_df_DF$s.S0F135 == 1 & indic_df_DF$p.value < 0.05),])
	S0F2_DF <- as.matrix(indic_df_DF[which(indic_df_DF$s.S0F240 == 1 & indic_df_DF$p.value < 0.05),])
	S1F1_DF <- as.matrix(indic_df_DF[which(indic_df_DF$s.S100F135 == 1 & indic_df_DF$p.value < 0.05),])
	S1F2_DF <- as.matrix(indic_df_DF[which(indic_df_DF$s.S100F240 == 1 & indic_df_DF$p.value < 0.05),])
	
	indic_sign_DF <- rbind(S0F1_DF,S0F2_DF,S1F1_DF,S1F2_DF)
	colnames(indic_sign_DF)[1:4] <-c("S0F135","S0F240","S100F135","S100F240")
	
	## Range of correlation coefficients
	range(indic_sign_DF[,"stat"])

    ## Total number of indicator OTUS
    length(unique(rownames(indic_sign_DF)))

    ## Proportion of soil bacteria OTUs responding to treatment
    length(unique(rownames(indic_sign_DF)))/nrow(indic_sign_DF)

    ## Proportion 
    DF_final_relative <- t(t(OTU_final_DF)/colSums(OTU_final_DF)) * 100
    sum(colSums(DF_final_relative[unique(rownames(DF_final_relative)),]))/sum(colSums(DF_final_relative))
	
    write.table(indic_sign_DF,paste0("indic_sign_DF.txt"),sep="\t",quote=F)

	 tax_indic_sign_DF <- taxon_otu[rownames(indic_sign_DF),]
	 t_DF <- as.data.frame(tax_indic_sign_DF)
	 t_indic_DF <-as.data.frame(indic_sign_DF)
	 DF_indi_otu_taxon <-  bind_cols(t_DF,t_indic_DF)
     write.csv(DF_indi_otu_taxon,file="DF_indi_otu_taxon.csv")
	 write.csv(DF_final_relative,file="DF_final_relative.csv")

	 DF_indic_str <- indicatorsp_DF$str
	 DF_indi_otu_taxon_str <- DF_indic_str[rownames(indic_sign_DF),]
     str_DF <- as.data.frame(DF_indi_otu_taxon_str)
	 taxon_DF <- as.data.frame(DF_indi_otu_taxon)
	 DF_indi_otu_taxon_str <-  bind_cols(str_DF,taxon_DF)
     write.csv(DF_indi_otu_taxon_str,file="DF_indi_otu_taxon_str.csv")


#####Spring 
     OTUs_TF <- read.csv('TF_OTUs.csv', stringsAsFactors = FALSE, row.names=1)
     OTUs_TF <- OTUs_TF[rowSums(OTUs_TF) > 0,]

     keep_OTUs_TF <- which(rowSums(OTUs_TF >= 2) >= 4) 
     OTU_final_TF <- OTUs_TF[keep_OTUs_TF,]
	 
	 taxon_otu <- read.csv('OTUs_taxonomy.csv',stringsAsFactors = FALSE, row.names=1) 
	 tax_otu_TF <- taxon_otu[rownames(OTU_final_TF),]

	 design_TF <- read.csv('design_name.csv') %>% filter(FTC=="TF")
	 edgeR_otu_TF <- DGEList(counts=OTU_final_TF, 
                          group=design_TF$Treatment,
                          genes=tax_otu_TF)
						  
     edgeR_otu_TF <- calcNormFactors(edgeR_otu_TF)

	 ## Get TMM normalized counts expressed as relative abundance counts per million 
     otu_norm_TF <- cpm(edgeR_otu_TF, normalized.lib.sizes=T, log=F)
	 
     ##### indicator species analysis #####
      indic_TF <- as.data.frame(t(otu_norm_TF))
      indic_groups_TF <- design_TF$Treatment
      length(unique(indic_groups_TF))
      set.seed(168168)
      indicatorsp_TF <- multipatt(indic_TF,indic_groups_TF,func = "r.g",control=how(nperm=9999))
      summary(indicatorsp_TF,alpha=1,indvalcomp=T)
	  summary(indicatorsp_TF)
      TF_indic_df <- indicatorsp_TF$sign
      write.table(TF_indic_df,paste0("TF_indic_df.txt"),sep="\t",quote=F)
	  
     ## Import data frame of indicator species to save time
      indic_df_TF <- read.table("TF_indic_df.txt", header=T, sep="\t")
	 
	S0F1_TF <- as.matrix(indic_df_TF[which(indic_df_TF$s.S0F135 == 1 & indic_df_TF$p.value < 0.05),])
	S0F2_TF <- as.matrix(indic_df_TF[which(indic_df_TF$s.S0F240 == 1 & indic_df_TF$p.value < 0.05),])
	S1F1_TF <- as.matrix(indic_df_TF[which(indic_df_TF$s.S100F135 == 1 & indic_df_TF$p.value < 0.05),])
	S1F2_TF <- as.matrix(indic_df_TF[which(indic_df_TF$s.S100F240 == 1 & indic_df_TF$p.value < 0.05),])
	
	BF_taxon <- as.data.frame(S1F2_TF)
	BF_taxon$OTU <- rownames(S1F2_TF)
	taxon_otu$OTU <- rownames(taxon_otu)
	BF_taxon <- left_join(BF_taxon, taxon_otu)
	write.csv(BF_taxon, "S1F2_TF_taxon.csv")

	indic_sign_TF <- rbind(S0F1_TF,S0F2_TF,S1F1_TF,S1F2_TF)
	colnames(indic_sign_TF)[1:4] <-c("S0F135","S0F240","S100F135","S100F240")
	
	## Range of correlation coefficients
	range(indic_sign_TF[,"stat"])

    ## Total number of indicator OTUS
    length(unique(rownames(indic_sign_TF)))

    ## Proportion of soil bacteria OTUs responding to treatment
    length(unique(rownames(indic_sign_TF)))/nrow(indic_sign_TF)

    ## Proportion 
    TF_final_relative <- t(t(OTU_final_TF)/colSums(OTU_final_TF)) * 100
    sum(colSums(TF_final_relative[unique(rownames(TF_final_relative)),]))/sum(colSums(TF_final_relative))
    
	write.table(indic_sign_TF,paste0("indic_sign_TF.txt"),sep="\t",quote=F)
 
    View <- filter(as.data.frame(indic_sign_TF), index=="4")
	length(rownames(View))

	 tax_indic_sign_TF <- taxon_otu[rownames(indic_sign_TF),]
	 t_TF <- as.data.frame(tax_indic_sign_TF)
	 t_indic_TF <-as.data.frame(indic_sign_TF)
	 TF_indi_otu_taxon <-  bind_cols(t_TF,t_indic_TF)
     write.csv(TF_indi_otu_taxon,file="TF_indi_otu_taxon.csv")

	 TF_indic_str <- indicatorsp_TF$str
	 TF_indi_otu_taxon_str <- TF_indic_str[rownames(indic_sign_TF),]
     str_TF <- as.data.frame(TF_indi_otu_taxon_str)
	 taxon_TF <- as.data.frame(TF_indi_otu_taxon)
	 TF_indi_otu_taxon_str <-  bind_cols(str_TF,taxon_TF)
     write.csv(TF_indi_otu_taxon_str,file="TF_indi_otu_taxon_str.csv")

