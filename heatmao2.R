# heatmap using heatmap function
# if (!require("gplots")) {
#   install.packages("gplots", dependencies = TRUE)
#   library(gplots)
# }
# if (!require("RColorBrewer")) {
#   install.packages("RColorBrewer", dependencies = TRUE)
#   library(RColorBrewer)
# }

write.csv(GO_Biological_Process, file="GO_Biological_Process.csv")
write.csv(GO_Cellular_Component, file="GO_Cellular_Component.csv")
write.csv(GO_Molecular_Function, file="GO_Molecular_Function.csv")

#load all the files
library(readr)
CC <- read_csv("GO_Cellular_Component_new.csv")
MF <- read_csv("GO_Molecular_Function_new.csv")
BP <- read_csv("GO_Biological_Process_new.csv")

#########################################################
### B) Reading in data and transform it into matrix format
#########################################################

#CC
library(gplots)
CCnew<- unite(CC,GO_Term_Description, c(GO_Term, Description),sep="  ", remove = FALSE) 
CCrnames <- as.matrix(CCnew[,1])                       # assign labels in column 1 to "rnames"
CCmat_data <- data.matrix(CCnew[,4:6])  # transform column 2-5 into a matrix
rownames(CCmat_data) <- CCrnames 
CCmat_data <- as.matrix(CCmat_data)

# assign row names
BPnew<- unite(BP,GO_Term_Description, c(GO_Term, Description),sep="   ", remove = FALSE) 
BPrnames <- as.matrix(BPnew[,1])                       # assign labels in column 1 to "rnames"
BPmat_data <- data.matrix(BPnew[,4:6])  # transform column 2-5 into a matrix
rownames(BPmat_data) <- BPrnames
library(tidyr)
#create seperate graphs for 
#Biological Process

MFnew<- unite(MF,GO_Term_Description, c(GO_Term, Description),sep="   ", remove = FALSE) 
MFrnames <- as.matrix(MFnew[,1])                       # assign labels in column 1 to "rnames"
MFmat_data <- data.matrix(MFnew[,4:6])  # transform column 2-5 into a matrix
rownames(MFmat_data) <- MFrnames  
MFmat_data <- as.matrix(MFmat_data)
#


#########################################################
### C) Customizing and plotting the heat map
#########################################################

# creates a own color palette from red to green
my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)

# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(-1,0,length=100),  
               seq(0.01,0.8,length=100),           
               seq(0.81,1,length=100)) 

BPmat_data[is.na(BPmat_data)] <-0
CCmat_data[is.na(CCmat_data)] <-0
MFmat_data[is.na(MFmat_data)] <-0


# creates a 5 x 5 inch image
dev.off()
while (!is.null(dev.list()))  dev.off()

png("BPheatmaps_in_r.png",    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

heatmap.2(BPmat_data,
           col = c(colorpanel(1,"white","white"), colorpanel(700,"darkblue","blue","white")),
           margins = c(4, 28),
          dendrogram = 'none',
           trace = 'none', 
           lhei = c(2, 8), #plot layout
           scale = c("none"),
           symbreaks = min(BPmat_data, na.rm=TRUE),
           cexRow = 0.5, cexCol = 0.7,
           density.info = "none",
           sepcolor = colorpanel(1,"black","black"),
           sepwidth=c(0.0001,0.0001),
          key.xlab = "p-value",
          colsep=0:ncol(BPmat_data),
          rowsep=0:nrow(BPmat_data),
          srtCol=45)
# for (i in 1:nrow(BP)) {
# mtext(BP$GO_Term[i], side=2, at=c(5,i), cex=0.5) }

dev.off
     #      get.uni <- !duplicated(BP$GO_Term)
     # text(x =5, y = seq(0.2,1.2, by=(1/nrow(BPrnames))),
     # labels = BP$GO_Term[get.uni],
     # las = 2, col = "#000000", cex = 0.4, xpd = TRUE)



# creates a 5 x 5 inch image
dev.off()
while (!is.null(dev.list()))  dev.off()

png("CCheatmaps_in_r.png",    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)# smaller font size

get.uni <- !duplicated(CC$GO_Term)

heatmap.2(CCmat_data,
          col = c(colorpanel(1,"white","white"), colorpanel(700,"darkblue","blue","white")),
          margins = c(34, 28),
          dendrogram = 'none',
          trace = 'none', 
          lhei = c(2, 8), #plot layout
          scale = c("none"),
          symbreaks = min(CCmat_data, na.rm=TRUE),
          cexRow = 0.5, cexCol = 0.7,
          density.info = "none",
          sepcolor = colorpanel(1,"black","black"),
          sepwidth=c(0.0001,0.0001),
          key.xlab = "p-value",
          colsep=0:ncol(BPmat_data),
          rowsep=0:nrow(BPmat_data),
          srtCol=45)
#           text(x =5, y = seq(0.4,1.2, by=(1/nrow(CCrnames))),
# labels = CC$GO_Term[get.uni],
# las = 2, col = "#000000", cex = 0.4, xpd = TRUE)

dev.off

# creates a 5 x 5 inch image
dev.off()
while (!is.null(dev.list()))  dev.off()

png("MFheatmaps_in_r.png",    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

heatmap.2(MFmat_data,
          col = c(colorpanel(1,"white","white"), colorpanel(700,"darkblue","blue","white")),
          margins = c(18, 28),
          dendrogram = 'none',
          trace = 'none', 
          lhei = c(2, 8), #plot layout
          scale = c("none"),
          symbreaks = min(MFmat_data, na.rm=TRUE),
          cexRow = 0.5, cexCol = 0.7,
          density.info = "none",
          sepcolor = colorpanel(1,"black","black"),
          sepwidth=c(0.0001,0.0001),
          key.xlab = "p-value",
          colsep=0:ncol(MFmat_data),
          rowsep=0:nrow(MFmat_data),
          srtCol=45)

#           get.uni <- !duplicated(BP$GO_Term)
# text(x =5, y = seq(0.2,1.2, by=(1/nrow(BPrnames))),
#      labels = BP$GO_Term[get.uni],
#      las = 2, col = "#000000", cex = 0.4, xpd = TRUE)

dev.off


