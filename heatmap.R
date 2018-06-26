##load the results
#library(readxl)
# SEAcompare_results_all <- read_excel("SEAcompare_results_all.xlsx", skip = 1)
# View(SEAcompare_results_all)

# head(SEAcompare_results_all)

##sort out the headings and put NAs in
change <- c("FDR","het_Number","FDR","pan1_Number","FDR","b73_Number")


for (i in 1:5) {
colnames(SEAcompare_results_all)[i+4] <- change[i]
SEAcompare_results_all[i+4][SEAcompare_results_all[i+4] == "0.5"] <- "NA"
}

#str(SEAcompare_results_all)
#change the structure to make new table of GOterm, description, sample (het, pan1, b73), FDR Value
pan1 <- cbind(SEAcompare_results_all[2:4],SEAcompare_results_all[7], rep('pan1',108))
colnames(pan1)[5] <- "sample_type"
het <-cbind(SEAcompare_results_all[2:4],SEAcompare_results_all[5], rep('het',108))
colnames(het)[5] <- "sample_type"
b73 <-cbind(SEAcompare_results_all[2:4],SEAcompare_results_all[9], rep('b73',108))

#column names
colnames(b73)[5] <- "sample_type"

#form new database
heatmap_all <-rbind(pan1, het, b73)
colnames(heatmap_all)[1] <- "GO_Term"
heatmap_all$FDR <-as.numeric(heatmap_all$FDR)
heatmap_all$`sample_type` <-as.character(heatmap_all$`sample_type`)
heatmap_all[2][heatmap_all[2] == "P"] <- "GO Biological Process"
heatmap_all[2][heatmap_all[2] == "C"] <- "GO Cellular Component"
heatmap_all[2][heatmap_all[2] == "F"] <- "GO Molecular Function"


#bringing together the GO term and the description
library(tidyr)


#create a heat map that
# Create color palette
library("RColorBrewer", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("ggplot2", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library(colorspace)

# myColors <- choose_palette(pal=diverge_hcl, n=8L)
myColors2 <-myColors(11)
myColors <- brewer.pal(8, "Reds")


# Build the heat map from scratch
par(mar=c(5,6,4,2)+0.1,mgp=c(5,1,0))
 # Adjust colors

#create seperate graphs for 
#Biological Process
GO_Biological_Process <- subset(heatmap_all, heatmap_all$Onto=="GO Biological Process")
# GO_Biological_Process<- unite(GO_Biological_Process,GO_Term_Description, c(GO_Term, Description), remove = FALSE) 
# GO_Biological_Process<- subset(GO_Biological_Process, select = c(GO_Term_Description, Onto, FDR, sample_type))

png(file="Biological_Process.png",width=600,height=1980)
ggplot(GO_Biological_Process, aes(x =sample_type, y = GO_Term_Description, fill = FDR)) +
  geom_tile() + # Geom layer 
  scale_fill_gradientn(colors = myColors2, na.value = myColors[1]) +
  labs(title="",x="", y="") +
          theme_bw()
dev.off()

#create seperate graphs for 
#Molecular_Function
GO_Molecular_Function <- subset(heatmap_all, heatmap_all$Onto=="GO Molecular Function")
# GO_Molecular_Function<- unite(GO_Molecular_Function,GO_Term_Description, c(GO_Term, Description), remove = FALSE)
# GO_Molecular_Function <- GO_Molecular_Function[order(subset(GO_Molecular_Function, sample_type==GO_Molecular_Function$FDR), ]

# GO_Molecular_Function<- subset(GO_Molecular_Function, select = c(GO_Term_Description, Onto, FDR, sample_type))

black<- "000000"
png(file="Molecular_Function.png",width=700,height=1000)
ggplot(GO_Molecular_Function, aes(x =sample_type, y = Description, fill = FDR), colour = colors()[265]) +
  geom_tile() + # Geom layer 
  scale_fill_gradientn(colors = myColors2, na.value = myColors[1])+
  labs(title="",x="", y="")+
  theme_bw()
dev.off()

#GO Cellular_Component"
GO_Cellular_Component <- subset(heatmap_all, heatmap_all$Onto=="GO Cellular Component")
# GO_Cellular_Component <- unite(GO_Cellular_Component,GO_Term_Description, c(GO_Term, Description), remove = FALSE) 
# GO_Cellular_Component <- GO_Cellular_Component[order(GO_Cellular_Component$FDR, GO_Cellular_Component$sample_type), ]
# #subset(GO_Cellular_Component, select = c(GO_Term_Description, Onto, FDR, sample_type))

GO_Cellular_Component$FDR <- as.numeric(GO_Cellular_Component$FDR)
GO_Cellular_Component[is.na(GO_Cellular_Component)]<-0.5

png(file="Cellular_Component.png",width=800,height=170)
ggplot(GO_Cellular_Component, aes(x =sample_type, y = GO_Term, fill = FDR)) +
  geom_tile() + # Geom layer 
  scale_fill_gradientn(colors = myColors2, na.value = myColors[1], limits=c(0,0.05) )+
  labs(title="",x="", y="")+
  theme_bw()

dev.off()

