library(openxlsx)
setwd("/Users/uemurakohei/Desktop/ncRNA/lncRNA/16_categories_inr")
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
gene_id_inr <- rownames(data_ex[which(data_ex[,"Inr_6mer"] == 1),])


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


 #


  for(DHS.Support in unique(lncRNA_id[,"DHS.Support"])){

    color.i.j <- 1
    count.i.j <-  0

    for(gene_class in unique(lncRNA_id[,"Gene.Class"])){



      lncRNA_id.i <- lncRNA_id[which((lncRNA_id[,"DHS.Support"]==DHS.Support)&(lncRNA_id[,"Gene.Class"]==gene_class)),"GeneID"]
      pos_rep_tss_coding_geneid <- which(reptss_gene_id %in% lncRNA_id.i)
      gene_id_rep_coding_tss <- reptss[pos_rep_tss_coding_geneid, 4]

      pos_res <- which(row.names(data_each) %in% intersect(gene_id_inr, gene_id_rep_coding_tss))

      data_mRNA_averagesd <- apply(data_each[pos_res,], 2, mean)
      data_mRNA_averagesd <- cbind(-150:150, data_mRNA_averagesd)
      data_mRNA_averagesd <- cbind(data_mRNA_averagesd, apply(data_each[pos_res,], 2, sd))
      colnames(data_mRNA_averagesd) <- c("position", "value","sd")


      x<- data_mRNA_averagesd[,"position"]
      mean_force<- data_mRNA_averagesd[,"value"]
      mean_shuffle <- apply(t(apply(data_each[pos_res,], 1, function(x) sample(x))) ,2, mean)
      #assume constant standard deviation across the
      #sd<-data_mRNA_averagesd[,"sd"]
      #determine error band
      #psd<-mean_force+sd
      #nsd<-mean_force-sd

      if(count.i.j == 0){

        pdf(sprintf("%s_%s_average.pdf", property_name, DHS.Support),width=10, height=3)
        #dev.new(width=7.5, height=3)
        par(mar = c(4, 5, 2, 1))
        par(family = "Helvetica")
        plot(x, mean_force, type="l", col="black",
             ylab="",
             xlab="",
             las=1,
             cex.axis = 1.8,
             #main=property_name,
             #lty=color.i.j,
             lwd=0.5,


             #ylim=c(min(c(data_mRNA_averagesd[,"value"]-data_mRNA_averagesd[,"sd"])),max(c(data_mRNA_averagesd[,"value"]+data_mRNA_averagesd[,"sd"]))))
             ylim=c(y_lim))

             #縦線
             abline(v=0, col="red", lty="dotted", lwd=0.5)
             abline(v=-27, col="blue", lty="dotted", lwd=0.5)

             points(x, mean_force,cex=0.3,
               #lty="dotted",
               pch=16,
               col=color.i.j)

               box(lwd = 2)

        #lines(-150:(-27-5), mean_force[(-150:(-27-5))+151], lwd=2, col=colors[color.i.j])
        #lines((-27-5):(-27+5), mean_force[((-27-5):(-27+5))+151], lty="dotted",lwd=2, col=colors[color.i.j])
        #lines((-27+5):(0-5), mean_force[((-27+5):(0-5))+151], lwd=2, col=colors[color.i.j])

        #draw boundary and fill
        #lines(x, psd, col="gray")
        #lines(x, nsd, col="gray")
        #points(x, mean_force,col=color.i.j, pch=color.i.j)

      }else if(count.i.j == 1){

        lines(x, mean_force,lwd=0.5,
          #lty="dotted",
          col="black")

          points(x, mean_force,cex=0.3,
            pch=16,
            #lty="dotted",
            col=color.i.j)

        #points(x, mean_force,col=color.i.j, pch=color.i.j)

      }else if(count.i.j == 2){

        lines(x, mean_force,lwd=0.5,
          #lty="dotted",
          col="black")

          points(x, mean_force,cex=0.3,
            pch=16,
            #lty="dotted",
            col=color.i.j)

      }else{

        lines(x, mean_force,lwd=0.5,
          #lty="dotted",
          col="black")

          points(x, mean_force,cex=0.3,
            pch=16,
            #lty="dotted",
            col=color.i.j)

      }

      color.i.j <- color.i.j + 1
      #polygon(x=c(x, rev(x)), y=c(psd, rev(nsd)), col="gray", density = 100, angle=90)
      #redraw line on top

      #shuffle
      #lines(x, mean_shuffle, col=rgb(1,0,0, alpha=0.5), lwd=1,lty=1)
      #lines(x, mean_shuffle, col="red", lwd=3)

      count.i.j <- count.i.j + 1

    }
      dev.off()
  }

}



for(i_data in list.files("/Users/uemurakohei/Desktop/ncRNA/lncRNA/Result", full=T)){

  property_name <- unlist(strsplit(i_data, split="/Users/uemurakohei/Desktop/ncRNA/lncRNA/Result/"))[2]
  data_each <- read.table(paste(i_data, "/", property_name, "_all.txt", sep=""), stringsAsFactors=F, row.names=1)
  colnames(data_each) <- -150:150


  print(property_name)

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


 #

  for(DHS.Support in unique(lncRNA_id[,"DHS.Support"])){

    color.i.j <- 1
    count.i.j <-  0

    for(gene_class in unique(lncRNA_id[,"Gene.Class"])){



      lncRNA_id.i <- lncRNA_id[which((lncRNA_id[,"DHS.Support"]==DHS.Support)&(lncRNA_id[,"Gene.Class"]==gene_class)),"GeneID"]
      pos_rep_tss_coding_geneid <- which(reptss_gene_id %in% lncRNA_id.i)
      gene_id_rep_coding_tss <- reptss[pos_rep_tss_coding_geneid, 4]

      pos_res <- which(row.names(data_each) %in% intersect(gene_id_inr, gene_id_rep_coding_tss))


      data_mRNA_averagesd <- apply(data_each[pos_res,as.character(-50:5)], 2, mean)
      data_mRNA_averagesd <- cbind(-50:5, data_mRNA_averagesd)
      data_mRNA_averagesd <- cbind(data_mRNA_averagesd, apply(data_each[pos_res,as.character(-50:5)], 2, sd))
      colnames(data_mRNA_averagesd) <- c("position", "value","sd")


      x<- data_mRNA_averagesd[,"position"]
      mean_force<- data_mRNA_averagesd[,"value"]
      mean_shuffle <- apply(t(apply(data_each[pos_res,], 1, function(x) sample(x))) ,2, mean)
      #assume constant standard deviation across the
      #sd<-data_mRNA_averagesd[,"sd"]
      #determine error band
      #psd<-mean_force+sd
      #nsd<-mean_force-sd

      if(count.i.j == 0){

        pdf(sprintf("%s_%s_average_50_20.pdf", property_name, DHS.Support),width=2.7, height=2)
        #dev.new(width=7.5, height=3
        par(mar = c(3, 4, 1, 1))
        par(oma = c(0, 0, 0, 0))
        par(family = "Helvetica")
        plot(x, mean_force, type="l", col="black",
             ylab="",
             xlab="",
             #xaxp=c(-50, 5, 10),
             las=1,
             cex.axis = 1.4,
             #main=property_name,
             #lty=color.i.j,
             lwd=0.5,


             #ylim=c(min(c(data_mRNA_averagesd[,"value"]-data_mRNA_averagesd[,"sd"])),max(c(data_mRNA_averagesd[,"value"]+data_mRNA_averagesd[,"sd"]))))
             ylim=c(y_lim))

             points(x, mean_force,cex=0.3,
               #lty="dotted",
               pch=16,
               col=color.i.j)

               box(lwd = 2)

        #lines(-150:(-27-5), mean_force[(-150:(-27-5))+151], lwd=2, col=colors[color.i.j])
        #lines((-27-5):(-27+5), mean_force[((-27-5):(-27+5))+151], lty="dotted",lwd=2, col=colors[color.i.j])
        #lines((-27+5):(0-5), mean_force[((-27+5):(0-5))+151], lwd=2, col=colors[color.i.j])

        #draw boundary and fill
        #lines(x, psd, col="gray")
        #lines(x, nsd, col="gray")
        #points(x, mean_force,col=color.i.j, pch=color.i.j)

      }else if(count.i.j == 1){

        lines(x, mean_force,lwd=0.5,
          #lty="dotted",
          col="black")

          points(x, mean_force,cex=0.3,
            pch=16,
            #lty="dotted",
            col=color.i.j)

        #points(x, mean_force,col=color.i.j, pch=color.i.j)

      }else if(count.i.j == 2){

        lines(x, mean_force,lwd=0.5,
          #lty="dotted",
          col="black")

          points(x, mean_force,cex=0.3,
            pch=16,
            #lty="dotted",
            col=color.i.j)

      }else{

        lines(x, mean_force,lwd=0.5,
          #lty="dotted",
          col="black")

          points(x, mean_force,cex=0.3,
            pch=16,
            #lty="dotted",
            col=color.i.j)

      }

      color.i.j <- color.i.j + 1
      #polygon(x=c(x, rev(x)), y=c(psd, rev(nsd)), col="gray", density = 100, angle=90)
      #redraw line on top

      #shuffle
      #lines(x, mean_shuffle, col=rgb(1,0,0, alpha=0.5), lwd=1,lty=1)
      #lines(x, mean_shuffle, col="red", lwd=3)

      count.i.j <- count.i.j + 1

    }
      dev.off()
  }

}




 # 50から20まで
for(i_data in list.files("/Users/uemurakohei/Desktop/ncRNA/lncRNA/Result", full=T)){

  property_name <- unlist(strsplit(i_data, split="/Users/uemurakohei/Desktop/ncRNA/lncRNA/Result/"))[2]
  data_each <- read.table(paste(i_data, "/", property_name, "_all.txt", sep=""), stringsAsFactors=F, row.names=1)
  colnames(data_each) <- -150:150


  print(property_name)

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


 #

  for(DHS.Support in unique(lncRNA_id[,"DHS.Support"])){

    color.i.j <- 1
    count.i.j <-  0

    for(gene_class in unique(lncRNA_id[,"Gene.Class"])){



      lncRNA_id.i <- lncRNA_id[which((lncRNA_id[,"DHS.Support"]==DHS.Support)&(lncRNA_id[,"Gene.Class"]==gene_class)),"GeneID"]
      pos_rep_tss_coding_geneid <- which(reptss_gene_id %in% lncRNA_id.i)
      gene_id_rep_coding_tss <- reptss[pos_rep_tss_coding_geneid, 4]

      pos_res <- which(row.names(data_each) %in% intersect(gene_id_inr, gene_id_rep_coding_tss))


      data_mRNA_averagesd <- apply(data_each[pos_res,as.character(-50:20)], 2, mean)
      data_mRNA_averagesd <- cbind(-50:20, data_mRNA_averagesd)
      data_mRNA_averagesd <- cbind(data_mRNA_averagesd, apply(data_each[pos_res,as.character(-50:20)], 2, sd))
      colnames(data_mRNA_averagesd) <- c("position", "value","sd")


      x<- data_mRNA_averagesd[,"position"]
      mean_force<- data_mRNA_averagesd[,"value"]
      mean_shuffle <- apply(t(apply(data_each[pos_res,], 1, function(x) sample(x))) ,2, mean)
      #assume constant standard deviation across the
      #sd<-data_mRNA_averagesd[,"sd"]
      #determine error band
      #psd<-mean_force+sd
      #nsd<-mean_force-sd

      if(count.i.j == 0){

        pdf(sprintf("%s_%s_average_50_20.pdf", property_name, DHS.Support),width=3, height=2)
        #dev.new(width=7.5, height=3
        par(mar = c(3, 4, 1, 1))
        par(oma = c(0, 0, 0, 0))
        par(family = "Helvetica")
        plot(x, mean_force, type="l", col="black",
             ylab="",
             xlab="",
             #xaxp=c(-50, 5, 10),
             las=1,
             cex.axis = 1.4,
             #main=property_name,
             #lty=color.i.j,
             lwd=0.5,


             #ylim=c(min(c(data_mRNA_averagesd[,"value"]-data_mRNA_averagesd[,"sd"])),max(c(data_mRNA_averagesd[,"value"]+data_mRNA_averagesd[,"sd"]))))
             ylim=c(y_lim))

             #縦線

             abline(v=0, col="red", lty="dotted", lwd=0.5)
             abline(v=-27, col="blue", lty="dotted", lwd=0.5)

             points(x, mean_force,cex=0.3,
               #lty="dotted",
               pch=16,
               col=color.i.j)

               box(lwd = 2)

        #lines(-150:(-27-5), mean_force[(-150:(-27-5))+151], lwd=2, col=colors[color.i.j])
        #lines((-27-5):(-27+5), mean_force[((-27-5):(-27+5))+151], lty="dotted",lwd=2, col=colors[color.i.j])
        #lines((-27+5):(0-5), mean_force[((-27+5):(0-5))+151], lwd=2, col=colors[color.i.j])

        #draw boundary and fill
        #lines(x, psd, col="gray")
        #lines(x, nsd, col="gray")
        #points(x, mean_force,col=color.i.j, pch=color.i.j)

      }else if(count.i.j == 1){

        lines(x, mean_force,lwd=0.5,
          #lty="dotted",
          col="black")

          points(x, mean_force,cex=0.3,
            pch=16,
            #lty="dotted",
            col=color.i.j)

        #points(x, mean_force,col=color.i.j, pch=color.i.j)

      }else if(count.i.j == 2){

        lines(x, mean_force,lwd=0.5,
          #lty="dotted",
          col="black")

          points(x, mean_force,cex=0.3,
            pch=16,
            #lty="dotted",
            col=color.i.j)

      }else{

        lines(x, mean_force,lwd=0.5,
          #lty="dotted",
          col="black")

          points(x, mean_force,cex=0.3,
            pch=16,
            #lty="dotted",
            col=color.i.j)

      }

      color.i.j <- color.i.j + 1
      #polygon(x=c(x, rev(x)), y=c(psd, rev(nsd)), col="gray", density = 100, angle=90)
      #redraw line on top

      #shuffle
      #lines(x, mean_shuffle, col=rgb(1,0,0, alpha=0.5), lwd=1,lty=1)
      #lines(x, mean_shuffle, col="red", lwd=3)

      count.i.j <- count.i.j + 1

    }
      dev.off()
  }

}




#　フレームを取得
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
  y_lim <- c(1.0 ,3.5)
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


 #


  for(DHS.Support in unique(lncRNA_id[,"DHS.Support"])){

    color.i.j <- 1
    count.i.j <-  0

    for(gene_class in unique(lncRNA_id[,"Gene.Class"])){



      lncRNA_id.i <- lncRNA_id[which((lncRNA_id[,"DHS.Support"]==DHS.Support)&(lncRNA_id[,"Gene.Class"]==gene_class)),"GeneID"]
      pos_rep_tss_coding_geneid <- which(reptss_gene_id %in% lncRNA_id.i)
      gene_id_rep_coding_tss <- reptss[pos_rep_tss_coding_geneid, 4]

      pos_res <- which(row.names(data_each) %in% gene_id_rep_coding_tss)

      data_mRNA_averagesd <- apply(data_each[pos_res,], 2, mean)
      data_mRNA_averagesd <- cbind(-150:150, data_mRNA_averagesd)
      data_mRNA_averagesd <- cbind(data_mRNA_averagesd, apply(data_each[pos_res,], 2, sd))
      colnames(data_mRNA_averagesd) <- c("position", "value","sd")


      x<- data_mRNA_averagesd[,"position"]
      mean_force<- data_mRNA_averagesd[,"value"]
      mean_shuffle <- apply(t(apply(data_each[pos_res,], 1, function(x) sample(x))) ,2, mean)
      #assume constant standard deviation across the
      #sd<-data_mRNA_averagesd[,"sd"]
      #determine error band
      #psd<-mean_force+sd
      #nsd<-mean_force-sd

      #if(count.i.j == 0){

        pdf(sprintf("%s_%s_%s_average.pdf", property_name, DHS.Support, gene_class),width=7.5, height=3)
        #dev.new(width=7.5, height=3)
        par(mar = c(4, 5, 2, 1))
        par(family = "Helvetica")
        plot(x, mean_force, type="l", col=color.i.j,
             ylab="",
             xlab="",
             las=1,
             cex.axis = 1.8,
             #main=property_name,
             lty=1,lwd=2,

             #ylim=c(min(c(data_mRNA_averagesd[,"value"]-data_mRNA_averagesd[,"sd"])),max(c(data_mRNA_averagesd[,"value"]+data_mRNA_averagesd[,"sd"]))))
             ylim=c(y_lim))

             dev.off()
        #draw boundary and fill
        #lines(x, psd, col="gray")
        #lines(x, nsd, col="gray")
      #}else{
      #  lines(x, mean_force, col=color.i.j,lwd=2)
      #}

      color.i.j <- color.i.j + 1
      #polygon(x=c(x, rev(x)), y=c(psd, rev(nsd)), col="gray", density = 100, angle=90)
      #redraw line on top

      #shuffle
      #lines(x, mean_shuffle, col=rgb(1,0,0, alpha=0.5), lwd=1,lty=1)
      #lines(x, mean_shuffle, col="red", lwd=3)

      count.i.j <- count.i.j + 1

    }

  }

}
