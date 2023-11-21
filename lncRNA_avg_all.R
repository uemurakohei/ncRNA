library(openxlsx)
setwd("/Users/uemurakohei/Desktop/ncRNA/lncRNA/16_categories")
reptss <- read.table("/Users/uemurakohei/Desktop/ncRNA/lncRNA/representative_tss/representative_tss_from_TIE_score.txt", stringsAsFactors=F)
reptss_gene_id <- sapply(reptss[,4] , function(x) unlist(strsplit(x, split="\\|"))[1])
lncRNA_id <- read.csv("/Users/uemurakohei/Desktop/ncRNA/lncRNA/16_categories/genes.csv")

#colors <- rainbow(4, alpha =1)

for(i_data in list.files("/Users/uemurakohei/Desktop/ncRNA/lncRNA/Result", full=T)){

  property_name <- unlist(strsplit(i_data, split="/Users/uemurakohei/Desktop/ncRNA/lncRNA/Result/"))[2]
  data_each <- read.table(paste(i_data, "/", property_name, "_all.txt", sep=""), stringsAsFactors=F, row.names=1)
  colnames(data_each) <- -150:150


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





# -50から10まで
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

      pos_res <- which(row.names(data_each) %in% gene_id_rep_coding_tss)


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

        pdf(sprintf("%s_%s_average_50_10.pdf", property_name, DHS.Support),width=2.7, height=2)
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




 # 40から20まで
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

      pos_res <- which(row.names(data_each) %in% gene_id_rep_coding_tss)


      data_mRNA_averagesd <- apply(data_each[pos_res,as.character(-40:20)], 2, mean)
      data_mRNA_averagesd <- cbind(-40:20, data_mRNA_averagesd)
      data_mRNA_averagesd <- cbind(data_mRNA_averagesd, apply(data_each[pos_res,as.character(-40:20)], 2, sd))
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

        pdf(sprintf("%s_%s_average_40_20.pdf", property_name, DHS.Support),width=3, height=2)
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
             abline(v=0, col="red", lty="dotted", lwd=1.5)
             abline(v=-27, col="blue", lty="dotted", lwd=1.5)

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







#　フレームを取得 tss
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
    slide_length <- 0

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
      if(gene_class == unique(lncRNA_id[,"Gene.Class"])[1]){
        pdf(sprintf("%s_%s_average_tss.pdf", property_name, DHS.Support),width=7.5, height=3)
        par(mar = c(4, 5, 2, 1))
        par(family = "Helvetica")
        #par(mfrow = c(1, 4))
        plot(x, mean_force,  col=color.i.j,
          type="n",
             ylab="",
             xlab="",
             las=1,
             cex.axis = 1.8,
             #main=property_name,
             lty=1,lwd=2,

             #ylim=c(min(c(data_mRNA_averagesd[,"value"]-data_mRNA_averagesd[,"sd"])),max(c(data_mRNA_averagesd[,"value"]+data_mRNA_averagesd[,"sd"]))))
             ylim=c(y_lim),
             bty="n",
             xaxt="n",yaxt="n"
           )
      }

        #dev.new(width=7.5, height=3)


        lines((-100-4+slide_length):(-100+4+slide_length), mean_force[which(x == -4):which(x == +4)], col=color.i.j,cex.axis = 1.8,lty=1,lwd=2)

             #dev.off()
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

      slide_length <- slide_length +20

    }

    dev.off()


  }

}





#　フレームを取得 -27
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
      if(gene_class == unique(lncRNA_id[,"Gene.Class"])[1]){
        pdf(sprintf("%s_%s_each_m27.pdf", property_name, DHS.Support),width=8, height=3)
        par(mar = c(1, 1, 1, 1))
        par(family = "Helvetica")
        par(mfrow = c(1, 4))
      }

        #dev.new(width=7.5, height=3)

        plot((-27-5):(-27+5), mean_force[which(x == (-27-5)):which(x == (-27+5))],  col=color.i.j,
          type="l",
             ylab="",
             xlab="",
             las=1,
             cex.axis = 1.8,
             #main=property_name,
             lty=1,lwd=20,

             #ylim=c(min(c(data_mRNA_averagesd[,"value"]-data_mRNA_averagesd[,"sd"])),max(c(data_mRNA_averagesd[,"value"]+data_mRNA_averagesd[,"sd"]))))
             #ylim=c(y_lim),
             bty="n",
             xaxt="n",yaxt="n"
           )
        #lines((-5):(+5), mean_force[which(x == -5):which(x == +5)], col=color.i.j,cex.axis = 1.8,lty=1,lwd=2)

             #dev.off()
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

    dev.off()


  }

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
