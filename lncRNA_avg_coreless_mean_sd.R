library(openxlsx)
setwd("/Users/uemurakohei/Desktop/ncRNA/lncRNA/16_categories_coreless_mean_sd")
reptss <- read.table("/Users/uemurakohei/Desktop/ncRNA/lncRNA/representative_tss/representative_tss_from_TIE_score.txt", stringsAsFactors=F)
reptss_gene_id <- sapply(reptss[,4] , function(x) unlist(strsplit(x, split="\\|"))[1])
lncRNA_id <- read.csv("/Users/uemurakohei/Desktop/ncRNA/lncRNA/16_categories/genes.csv")

### element
data_elements <- c()
for(i_chr in paste("chr", c( 1:22, "X", "Y"), sep="")){
  data_elements <- rbind(data_elements, read.table(sprintf("/Users/uemurakohei/Desktop/ncRNA/lncRNA/fig1_core_promoter_element/Result_pls5/%s.core_motifs.txt", i_chr), header=T))
}

DCE_I_II_III <-  apply(data_elements[,c("DCE_I", "DCE_II", "DCE_III")], 1, sum)
DCE_I_II_III[which(DCE_I_II_III != 3)] <- 0
DCE_I_II_III[which(DCE_I_II_III == 3)] <- 1
data_elements <- cbind(data_elements, DCE_I_II_III)

data_ex <- data_elements
data_ex <- data_ex[,-which(colnames(data_ex) == "DCE_III")]
data_ex <- data_ex[,-which(colnames(data_ex) == "DCE_II")]
data_ex <- data_ex[,-which(colnames(data_ex) == "DCE_I")]
data_ex <- data_ex[,-which(colnames(data_ex) == "TATA_7mer")]
data_ex <- data_ex[,-which(colnames(data_ex) == "Inr_2er")]

data_ex <- cbind(data_ex, as.numeric((data_ex[,which(colnames(data_ex) == "DPE_1")]==1)|(data_ex[,which(colnames(data_ex) == "DPE_2")]==1)))
colnames(data_ex) <- c(colnames(data_ex)[1:(ncol(data_ex)-1)], "DPE")
data_ex <- data_ex[,-which(colnames(data_ex) == "DPE_1")]
data_ex <- data_ex[,-which(colnames(data_ex) == "DPE_2")]

data_ex <- data_ex[reptss[,4],]
data_ex_count <- apply(data_ex, 1, sum)

### core less
gene_id_coreless <- names(which(data_ex_count == 0))


for(i_data in list.files("/Users/uemurakohei/Desktop/ncRNA/lncRNA/Result", full=T)){

  property_name <- unlist(strsplit(i_data, split="/Users/uemurakohei/Desktop/ncRNA/lncRNA/Result/"))[2]
  data_each <- read.table(paste(i_data, "/", property_name, "_all.txt", sep=""), stringsAsFactors=F, row.names=1)


  if(property_name == "Bruckner"){
  y_lim <- c(-0.2 ,0.1)
  }

  if(property_name == "DNA_bending_stiffness"){
  y_lim <- c(40 ,85)
  }

  if(property_name == "DNA_denaturation"){
  y_lim <- c(-1.7 ,-1.0)
  }

  if(property_name == "Duplex_disrupt_energy"){
  y_lim <- c(1.6 ,2.6)
  }

  if(property_name == "Duplex_free_energy"){
  y_lim <- c(-2.0 ,-1.2)
  }

  if(property_name == "flexibility"){
  y_lim <- c(10.5 ,20)
  }

  if(property_name == "Protein_induced_deformability_Bp"){
  y_lim <- c(1.2 ,3.3)
  }

  if(property_name == "Stabilizing_energy_of_Z_DNA_AS"){
  y_lim <- c(1.5 ,4.5)
  }

  if(property_name == "Stabilizing_energy_of_Z_DNA_SA"){
  y_lim <- c(2 ,5)
  }

  if(property_name == "Stacking_energy"){
  y_lim <- c(-10 ,-6) #c(-12 ,-3.5)
  }


  # data
  data_xlsx <- matrix(-150:150)
  colnames(data_xlsx) <- "position"

  for(DHS.Support in unique(lncRNA_id[,"DHS.Support"])){


    for(gene_class in unique(lncRNA_id[,"Gene.Class"])){



      lncRNA_id.i <- lncRNA_id[which((lncRNA_id[,"DHS.Support"]==DHS.Support)&(lncRNA_id[,"Gene.Class"]==gene_class)),"GeneID"]
      pos_rep_tss_coding_geneid <- which(reptss_gene_id %in% lncRNA_id.i)
      gene_id_rep_coding_tss <- reptss[pos_rep_tss_coding_geneid, 4]

      pos_res <- which(row.names(data_each) %in% intersect(gene_id_coreless, gene_id_rep_coding_tss))

      data_mRNA_averagesd <- apply(data_each[pos_res,], 2, mean)
      data_mRNA_averagesd <- cbind(-150:150, data_mRNA_averagesd)
      data_mRNA_averagesd <- cbind(data_mRNA_averagesd, apply(data_each[pos_res,], 2, sd))
      colnames(data_mRNA_averagesd) <- c("position", "mean","sd")


      x<- data_mRNA_averagesd[,"position"]
      mean_force<- data_mRNA_averagesd[,"mean"]
      mean_shuffle <- apply(t(apply(data_each[pos_res,], 1, function(x) sample(x))) ,2, mean)
      #assume constant standard deviation across the
      sd<-data_mRNA_averagesd[,"sd"]
      #determine error band
      psd<-mean_force+sd
      nsd<-mean_force-sd

      pdf(sprintf("%s_%s_%s.pdf", property_name, DHS.Support, gene_class),width=7.5, height=3)
      #dev.new(width=7.5, height=3)
      par(mar = c(4, 5, 2, 1))
      par(family = "Helvetica")


      plot(x, mean_force, type="l", col="black",
           ylab="",
           xlab="",
           las=1,
           cex.axis = 1.8,
           #main=property_name,
           lty=1,lwd=2,

           ylim=c(min(c(data_mRNA_averagesd[,"mean"]-data_mRNA_averagesd[,"sd"])),max(c(data_mRNA_averagesd[,"mean"]+data_mRNA_averagesd[,"sd"])))
           #ylim=c(y_lim)
         )
      #draw boundary and fill
      lines(x, psd, col="gray")
      lines(x, nsd, col="gray")

      polygon(x=c(x, rev(x)), y=c(psd, rev(nsd)), col="gray", density = 100, angle=90)
      #redraw line on top

      #shuffle
      #lines(x, mean_shuffle, col=rgb(1,0,0, alpha=0.5), lwd=1,lty=1)
      #lines(x, mean_shuffle, col="red", lwd=3)
      lines(x, mean_force, col="black",lwd=4)
      #points(-27, mean_force[which(x == -27)], col=rgb(1,0,0),cex = 1.5, lwd=2)
      dev.off()

      colnames(data_mRNA_averagesd) <- paste(colnames(data_mRNA_averagesd), sprintf("%s_%s",DHS.Support, gene_class))

      data_xlsx <- cbind(data_xlsx, data_mRNA_averagesd[,2:3])

    }
  }

  write.xlsx(data.frame(data_xlsx), file=sprintf("%s.xlsx", property_name), Rownamse=F)

}
