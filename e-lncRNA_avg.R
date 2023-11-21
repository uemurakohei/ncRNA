library(openxlsx)
setwd("/Users/uemurakohei/Desktop/paper3/lncRNA/average_physical/avg_profile/e_lncRNA")
reptss <- read.table("/Users/uemurakohei/Desktop/paper3/lncRNA/representative_tss/representative_tss_from_TIE_score.txt", stringsAsFactors=F)
lncRNA_id <- read.table("/Users/uemurakohei/Desktop/paper3/lncRNA/divergent_e_lncRNA/FANTOM_CAT.lv3_robust.only_e_lncRNA.bed")
lncRNA_id <- unique(sapply(lncRNA_id[,4], function(x) unlist(strsplit(x, "\\|"))[1]))
reptss_gene_id <- sapply(reptss[,4] , function(x) unlist(strsplit(x, split="\\|"))[1])

pos_rep_tss_coding_geneid <- which(reptss_gene_id %in% lncRNA_id)
gene_id_rep_coding_tss <- reptss[pos_rep_tss_coding_geneid, 4]

properties <- list.files("/Users/uemurakohei/Desktop/paper3/lncRNA/Result")
data_features <- matrix(rep("zzzz", length(properties)*12), ncol=length(properties), nrow=12)
colnames(data_features) <- properties
rownames(data_features) <- c("pos_max", "value_max", "pos_min", "value_min", "pos_max_delta", "value_max_delta", "pos_min_delta", "value_min_delta", "pos_up_max", "value_up_max", "pos_up_min", "value_up_min")


for(i_data in list.files("/Users/uemurakohei/Desktop/paper3/lncRNA/Result", full=T)){

  property_name <- unlist(strsplit(i_data, split="/Users/uemurakohei/Desktop/paper3/lncRNA/Result/"))[2]
  data_each <- read.table(paste(i_data, "/", property_name, "_all.txt", sep=""), stringsAsFactors=F, row.names=1)
  pos_res <- which(row.names(data_each) %in% gene_id_rep_coding_tss)

  data_mRNA_averagesd <- apply(data_each[pos_res,], 2, mean)
  data_mRNA_averagesd <- cbind(-150:150, data_mRNA_averagesd)
  data_mRNA_averagesd <- cbind(data_mRNA_averagesd, apply(data_each[pos_res,], 2, sd))
  colnames(data_mRNA_averagesd) <- c("position", "value","sd")

  pos_max <- data_mRNA_averagesd[which(data_mRNA_averagesd[,"value"] == max(data_mRNA_averagesd[,"value"])), "position"]
  value_max <- max(data_mRNA_averagesd[,"value"])

  pos_min <- data_mRNA_averagesd[which(data_mRNA_averagesd[,"value"] == min(data_mRNA_averagesd[,"value"])), "position"]
  value_min <- min(data_mRNA_averagesd[,"value"])

  #delta#
  data_mRNA_delta <- data_mRNA_averagesd[,"value"][1:(length(data_mRNA_averagesd[,"value"])-1)] - data_mRNA_averagesd[,"value"][2:length(data_mRNA_averagesd[,"value"])]
  data_mRNA_delta <- cbind(-50:49, data_mRNA_delta)
  colnames(data_mRNA_delta) <- c("position", "value")

  pos_max_delta <- data_mRNA_delta[which(data_mRNA_delta[,"value"] == max(data_mRNA_delta[,"value"])), "position"]
  value_max_delta <- max(data_mRNA_delta[,"value"])

  pos_min_delta <- data_mRNA_delta[which(data_mRNA_delta[,"value"] == min(data_mRNA_delta[,"value"])), "position"]
  value_min_delta <- min(data_mRNA_delta[,"value"])

  data_mRNA_up_averaged <- data_mRNA_averagesd[which(data_mRNA_averagesd[,"position"] %in% (-15):(-50)),]
  pos_up_max <- data_mRNA_up_averaged[which(data_mRNA_up_averaged[,"value"] == max(data_mRNA_up_averaged[,"value"])), "position"]
  value_up_max <- max(data_mRNA_up_averaged[,"value"])

  pos_up_min <- data_mRNA_up_averaged[which(data_mRNA_up_averaged[,"value"] == min(data_mRNA_up_averaged[,"value"])), "position"]
  value_up_min <- min(data_mRNA_up_averaged[,"value"])

  data_features[,property_name] <- sapply(rownames(data_features), function(x){data_features[x,property_name] <- get(x)})

  print(property_name)
  #pdf(sprintf("%s.pdf", property_name))


  x<- data_mRNA_averagesd[,"position"]
  mean_force<- data_mRNA_averagesd[,"value"]
  mean_shuffle <- apply(t(apply(data_each[pos_res,], 1, function(x) sample(x))) ,2, mean)
  #assume constant standard deviation across the
  sd<-data_mRNA_averagesd[,"sd"]
  #determine error band
  psd<-mean_force+sd
  nsd<-mean_force-sd

  if(property_name == "Bruckner"){
  y_lim <- c(-0.3 ,0.3)
  }

  if(property_name == "DNA_bending_stiffness"){
  y_lim <- c(10 ,120)
  }

  if(property_name == "DNA_denaturation"){
  y_lim <- c(-2.4 ,-0.5)
  }

  if(property_name == "Duplex_disrupt_energy"){
  y_lim <- c(0.8 ,3.5)
  }

  if(property_name == "Duplex_free_energy"){
  y_lim <- c(-2.4 ,-0.75)
  }

  if(property_name == "flexibility"){
  y_lim <- c(2 ,25)
  }

  if(property_name == "Protein_induced_deformability_Bp"){
  y_lim <- c(0 ,4)
  }

  if(property_name == "Stabilizing_energy_of_Z_DNA_AS"){
  y_lim <- c(0.8 ,6)
  }

  if(property_name == "Stabilizing_energy_of_Z_DNA_SA"){
  y_lim <- c(1 ,6)
  }

  if(property_name == "Stacking_energy"){
  y_lim <- c(-12 ,-3.5)
  }


  pdf(sprintf("%s_average.pdf", property_name),width=7.5, height=3.5)
  #dev.new(width=5, height=4)
  par(mar = c(4, 5, 2, 1))
  par(family = "Helvetica")
  plot(x, mean_force, ty="l", col="black",
       ylab="",
       xlab="",
       las=1,
       cex.axis = 1.8,
       #main=property_name,
       lty=1,lwd=2,

       #ylim=c(min(c(data_mRNA_averagesd[,"value"]-data_mRNA_averagesd[,"sd"])),max(c(data_mRNA_averagesd[,"value"]+data_mRNA_averagesd[,"sd"]))))
       ylim=c(y_lim))
  #draw boundary and fill
  lines(x, psd, col="gray")
  lines(x, nsd, col="gray")

  polygon(x=c(x, rev(x)), y=c(psd, rev(nsd)), col="gray", density = 100, angle=90)
  #redraw line on top

  #shuffle
  #lines(x, mean_shuffle, col=rgb(1,0,0, alpha=0.5), lwd=1,lty=1)
  #lines(x, mean_shuffle, col="red", lwd=3)
  lines(x, mean_force, col="black",lwd=4)


  dev.off()

  write.xlsx(data.frame(cbind(x, mean_force, sd)), file=sprintf("%s_average.xlsx", property_name), colNames=T)

}


data_elements <- c()
for(i_chr in paste("chr", c( 1:22, "X", "Y"), sep="")){
  data_elements <- rbind(data_elements, read.table(sprintf("/Users/uemurakohei/Desktop/program_paper/core_element/Result_pls5/%s.core_motifs.txt", i_chr), header=T))
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


data_ex <- data_ex[gene_id_rep_coding_tss,]
data_ex_count <- apply(data_ex, 1, sum)

gene_id_coreless <- names(which(data_ex_count == 0))
for(i_id in colnames(data_ex)){
  eval(parse(text=sprintf("gene_id_%s <- rownames(data_ex)[which(data_ex[,\"%s\"]==1)]", i_id, i_id)))
}

#only plot for coreless#
##[1] i_data <- "/Users/uemurakohei/Desktop/program_paper/average_physical/all_data/Bruckner_chr_all.txt"
##[2] i_data <- "/Users/uemurakohei/Desktop/program_paper/average_physical/all_data/DNA_bending_stiffness_chr_all.txt"
##[3] i_data <- "/Users/uemurakohei/Desktop/program_paper/average_physical/all_data/DNA_denaturation_chr_all.txt"
##[4] i_data <- "/Users/uemurakohei/Desktop/program_paper/average_physical/all_data/Duplex_disrupt_energy_chr_all.txt"
##[5] i_data <-"/Users/uemurakohei/Desktop/program_paper/average_physical/all_data/Duplex_free_energy_chr_all.txt"
##[6] i_data <- "/Users/uemurakohei/Desktop/program_paper/average_physical/all_data/flexibility_chr_all.txt"
##[7] i_data <- "/Users/uemurakohei/Desktop/program_paper/average_physical/all_data/Protein_induced_deformability_Bp_chr_all.txt"
##[8] "/Users/uemurakohei/Desktop/program_paper/average_physical/all_data/Stabilizing_energy_of_Z_DNA_AS_chr_all.txt"
##[9] "/Users/uemurakohei/Desktop/program_paper/average_physical/all_data/Stabilizing_energy_of_Z_DNA_SA_chr_all.txt"
##[10] i_data <- "/Users/uemurakohei/Desktop/program_paper/average_physical/all_data/Stacking_energy_chr_all.txt"

#setwd("")
for(i_data in list.files("/Users/uemurakohei/Desktop/program_paper/average_physical/all_data", full=T)){

  data_each <- read.table(i_data, stringsAsFactors=F, row.names=1)
  property_name <- unlist(strsplit(i_data, split="/Users/uemurakohei/Desktop/program_paper/average_physical/all_data/|_chr_all.txt"))[2]

  for(i_core in c("coreless", colnames(data_ex))){
    #i_core <- "BREd"
    #i_core <- "BREu"
    #i_core <- "TATA_8mer"
    #i_core <- "TCT"
    #i_core <- "Inr_6mer"
    #i_core <- "XCPE1"
    #i_core <- "DRE"
    #i_core <- "DCE_I_II_III"

    pos_res <- which(row.names(data_each) %in% get(sprintf("gene_id_%s", i_core)))

    if(!length(pos_res)){
      next
    }

    print(property_name)
    print(i_core)

    data_mRNA_averagesd <- apply(data_each[pos_res,], 2, mean)
    data_mRNA_averagesd <- cbind(-50:50, data_mRNA_averagesd)
    data_mRNA_averagesd <- cbind(data_mRNA_averagesd, apply(data_each[pos_res,], 2, sd))

    colnames(data_mRNA_averagesd) <- c("position", "value","sd")

    x<- data_mRNA_averagesd[,"position"]
    mean_force<- data_mRNA_averagesd[,"value"]
    #assume constant standard deviation across the
    sd<-data_mRNA_averagesd[,"sd"]
    #determine error band
    psd<-mean_force+sd
    nsd<-mean_force-sd
    mean_shuffle <- apply(t(apply(data_each[pos_res,], 1, function(x) sample(x))) ,2, mean)

    if(property_name == "Bruckner"){
    y_lim <- c(-0.3 ,0.3)
    }

    if(property_name == "DNA_bending_stiffness"){
    y_lim <- c(10 ,120)
    }

    if(property_name == "DNA_denaturation"){
    y_lim <- c(-2.4 ,-0.5)
    }

    if(property_name == "Duplex_disrupt_energy"){
    y_lim <- c(0.8 ,3.5)
    }

    if(property_name == "Duplex_free_energy"){
    y_lim <- c(-2.4 ,-0.75)
    }

    if(property_name == "flexibility"){
    y_lim <- c(2 ,25)
    }

    if(property_name == "Protein_induced_deformability_Bp"){
    y_lim <- c(0 ,4)
    }

    if(property_name == "Stabilizing_energy_of_Z_DNA_AS"){
    y_lim <- c(0.8 ,6)
    }

    if(property_name == "Stabilizing_energy_of_Z_DNA_SA"){
    y_lim <- c(1 ,6)
    }

    if(property_name == "Stacking_energy"){
    y_lim <- c(-12 ,-3.5)
    }

    pdf(sprintf("%s_%s_average.pdf", property_name, i_core),width=6, height=4)
    #dev.new(width=5.2, height=4)
    par(mar = c(4, 5, 2, 1))
    par(family = "Helvetica")
    plot(x, mean_force, ty="l", col="black",
         ylab="",
         xlab="",
         las=1,
         cex.axis = 1.8,
         #main=property_name,
         lty=1,lwd=2,
         #ylim=c(min(c(data_mRNA_averagesd[,"value"]-data_mRNA_averagesd[,"sd"])), max(c(data_mRNA_averagesd[,"value"]+data_mRNA_averagesd[,"sd"])))
         ylim=c(y_lim)
       )
    #draw boundary and fill
    lines(x, psd, col="gray")
    lines(x, nsd, col="gray")
    polygon(x=c(x, rev(x)), y=c(psd, rev(nsd)), col="gray", density = 100, angle=90)
    #redraw line on top
    #lines(x, mean_shuffle, col="red",lwd=3)
    lines(x, mean_force, col="black",lwd=4)

    dev.off()

    write.xlsx(data.frame(cbind(x, mean_force, sd)), file=sprintf("%s_%s_average.xlsx", property_name, i_core), colNames=T)

  }

}
