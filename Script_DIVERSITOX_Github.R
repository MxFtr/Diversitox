#######################################################
#                 DIVERSITOX 20203                    #
#                     R script                        #
#   Maxime Fuster https://github.com/MxFtr/Diversitox #
#######################################################

# Library -----------------------------------------------------------------
library(gdata)
library(ggplot2)
library(vegan)
library(phyloseq)
library(paletteer)
library(dplyr)
library(knitr)
library("VennDiagram")
library(tidyr)
library(car)
library(rstatix)
library(UpSetR)
library(scales)

BiocManager::install("MicrobiotaProcess") 
library(MicrobiotaProcess) # WARNING : MicrobiotaProcess and phyloseq can't be charged together



# Creation of phyloseq object = treatment of OTU table --------------------

# Before that you need to have in your working directory three .csv files 
# corresponding to your OTUs table : OTU= number of OTUs in your samples ;
# TAX= taxonomy affiliation of your OTUs ; META = metadata of your samples
# You can find exemples of these files on Github : https://github.com/MxFtr/Diversitox

dataOTU = read.table("OTU.csv", header=TRUE, sep=";")
dataOTU[,1] -> rownames(dataOTU)
dataOTU<-dataOTU[,2:37]#replace by your numbers of columnns 
dataTAX = read.table("TAX.csv", header=TRUE, sep=";")
dataTAX[,1] -> rownames(dataTAX)
dataTAX<-dataTAX[,2:8]#replace by your numbers of columnns
dataMETA = read.table("META.csv", header=TRUE, sep=";")
dataMETA[,1] -> rownames(dataMETA)
dataMETA<-dataMETA[,2:3]#replace by your numbers of columnns

rownames(dataMETA) -> colnames(dataOTU)

OTU = otu_table(dataOTU, taxa_are_rows = TRUE)
dataTAX<-as.matrix(dataTAX, rownames.force = TRUE)
TAX = tax_table(dataTAX)
SAM = sample_data(dataMETA)
physeq = phyloseq(OTU, TAX, SAM) #creaton of phyloseq object



# DIVERSITY INDEX ---------------------------------------------------------
#data preparation
indiv <- estimate_richness(physeq, split = TRUE, measures = NULL)
indiv$sample<-rownames(indiv)
indiv$sample[1:18]="Free"
indiv$sample[19:36]="Attached"
indiv <- select(indiv,one_of(c("Observed","Chao1","Shannon","sample")))
indiv<-pivot_longer(indiv,cols = 1:3,names_to = "Type",values_to = "Values")
indiv$sample<-as.factor(indiv$sample)

#creation of graphic
p<- ggplot(indiv,aes(x=sample,y=Values)) +
  geom_boxplot(aes(fill=sample))+
  facet_wrap("Type",scales = "free") +
  theme_bw() +
  labs(fill="Fraction",x="Fraction",y="Index") +
  scale_fill_manual(values = c("#8cae77","#758bb4"))
p
ggsave("Diversity_Index.svg",device = "svg")

#statiscal tests
chao1<-indiv %>%
  filter(Type=="Chao1") 
wilcox.test(Values~sample,data = chao1,paired=T)

Observed<-indiv %>%
  filter(Type=="Observed") 
wilcox.test(Values~sample,data = Observed,paired=T)

Shannon<-indiv %>%
  filter(Type=="Shannon") 
wilcox.test(Values~sample,data = Shannon,paired=T)


# Venn Diagram ------------------------------------------------------------
vennlist <- get_vennlist(obj=physeq, factorNames="Fraction")
vennp <- venn.diagram(vennlist, 
                      category.names = c("Attached" , "Free"),
                      filename = NULL,
                      compression = "lzw",
                      lwd = 1,
                      col=c("#8cae77","#758bb4"),
                      fill = c(alpha("#8cae77",0.3), alpha('#758bb4',0.3)),
                      cex = 0.8,
                      main.fontfamily = "sans",
                      main.fontface = "bold",
                      fontfamily = "sans",
                      cat.cex = 1,
                      cat.default.pos = "outer",
                      #cat.pos = c(-35,35),
                      cat.fontfamily = "sans",
                      cat.fontface = "bold",
                      cat.col = c("#8cae77", '#758bb4'))
grid::grid.draw(vennp)


# PCoA : Principal Coordonate Analysis ------------------------------------
#Create a list to organize the dates in the desired order
desired_order_Date <- list("September_14th","September_21th",
                           "October_10th","October_19th",
                           "October_26th","November_02nd")

#creation of ordination
Phy.ord.pcoa <- ordinate(physeq, "PCoA", "bray")

#creation of plot
pcoa<-plot_ordination(physeq,Phy.ord.pcoa,color = "Fraction",
                               shape = "Date") +
  geom_point(size=3)+
  theme_linedraw()+ 
  stat_ellipse(level = 0.95, linetype = 2)+
  scale_colour_manual(values = c("#8cae77","#758bb4"))

#date reorganization
pcoa$data$Date<- factor(pcoa$data$Date, levels = desired_order_Date)
pcoa




# MixMC pipeline for sPLS-DA applicated to OTU Data  ----------------------
#based on this good example: http://mixomics.org/mixmc/koren-bodysites-case-study/
#Read this article before using this part of script

#Import DATAs
#You can find exemples of these files on Github : https://github.com/MxFtr/Diversitox
library(mixOmics)
x=read.csv2("x_sPLSDA.csv",header = T,sep = ";")
x[,1] -> rownames(x)
x<-x[,2:688]
dim(x)
y=read.csv2("y_sPLSDA.csv",header = T,sep = ";")
y[,1] -> rownames(y)
y[,2] <- as.factor(y[,2])
y<-y[,2]
summary(y)

# Step 1 : Applying the offset 
x <- x+1 #essential : to remove 0
sum(which(x == 0))

# Step 2 : Pre-filtering 
#This step involves removing OTUs (features) for which the sum of counts are 
#below a certain threshold compared to the total sum of all counts. 
#The function is given below and was adapted from Arumugam et al., (2011).
#Function creation
low.count.removal <- function(
    data, # OTU count df of size n (sample) x p (OTU)
    percent=0.01 # cutoff chosen
) 
{
  keep.otu = which(colSums(data)*100/(sum(colSums(data))) > percent)
  data.filter = data[,keep.otu]
  return(list(data.filter = data.filter, keep.otu = keep.otu))
}
dim(x)
# call the function then apply on the  data
result.filter <- low.count.removal(x, percent=0.01)
data.filter <- result.filter$data.filter
length(result.filter$keep.otu)
x=data.filter

# Step 3 : INITIAL ANALYSIS
#Preliminary Analysis with PCA
cyanosphere.pca = pca(x, ncomp = 10, logratio = 'CLR') # undergo PCA with 10 comps
plot(cyanosphere.pca) # plot explained variance

cyanosphere.pca = pca(x, ncomp = 2, logratio = 'CLR') # undergo PCA with 10 comps
plotIndiv(cyanosphere.pca, # plot samples projected onto PCs
          ind.names = FALSE, # not showing sample names
          group = y, # color according to Y
          legend = TRUE,
          title = 'OTUs of cyanosphère, PCA Comps 1&2')

#Initial sPLS-DA model
basic.cyano.plsda = plsda(x, y, logratio = 'CLR', 
                          ncomp = nlevels(y))


# Step 4 : Tuning sPLS-DA 
#The Ncomp Parameter
# assess the performance of the sPLS-DA model using repeated CV
basic.cyano.perf.plsda = perf(basic.cyano.plsda,  
                              validation = 'Mfold', 
                              folds = 5, nrepeat = 50, 
                              progressBar = FALSE)
plot(basic.cyano.perf.plsda, overlay = 'measure', sd=TRUE) # plot this tuning

# extract the optimal component number
optimal.ncomp <- basic.cyano.perf.plsda$choice.ncomp["BER", "max.dist"] 

#The keepX Parameter
grid.keepX = c(seq(5,150, 5))

cyano.tune.splsda = tune.splsda(x, y,
                                ncomp = optimal.ncomp, # use optimal component number
                                logratio = 'CLR', # transform data to euclidean space
                                test.keepX = grid.keepX,
                                validation = c('Mfold'),
                                folds = 5, nrepeat = 50, # use repeated CV
                                dist = 'max.dist', # maximum distance as metric
                                progressBar = FALSE)

# extract the optimal component number and optimal feature count per component
optimal.keepX = cyano.tune.splsda$choice.keepX
optimal.ncomp = cyano.tune.splsda$choice.ncomp$ncomp 

plot(cyano.tune.splsda) # plot this tuning



# Step 5 : Final Model
cyano.splsda = splsda(x,  y, logratio= "CLR", # form final sPLS-DA model
                      ncomp = optimal.ncomp, 
                      keepX = optimal.keepX)


# Step 6 : Plots 
plotIndiv(cyano.splsda,
          comp = c(1,2),
          ind.names = FALSE,
          ellipse = TRUE, # include confidence ellipses
          legend = TRUE,
          legend.title = "Days",
          title = 'sPLS-DA')

cim<-cim(cyano.splsda,
         comp = c(1,2),
         row.sideColors = color.mixo(y), # colour rows based on bodysite
         legend = list(legend = c(levels(y))),
         title = 'Clustered Image Map of Cyanosphere bloom data')


# Step 7 : Keep Important OTUs 
library(dplyr)

write.csv2((first<-plotLoadings(cyano.splsda, comp = 1, 
                                method = 'mean', contrib = 'max',  
                                size.name = 0.8, legend = FALSE,  
                                ndisplay = 30,
                                title = "(a) Loadings of first component")),"first.csv")

write.csv2((second<-plotLoadings(cyano.splsda, comp = 2, 
                                 method = 'mean', contrib = 'max',   
                                 size.name = 0.7,
                                 ndisplay = 150,
                                 title = "(b) Loadings of second comp.")),"second.csv")

write.csv2((third<-plotLoadings(cyano.splsda, comp = 3, 
                                method = 'mean', contrib = 'max',   
                                size.name = 0.7,
                                ndisplay = 150,
                                title = "(c) Loadings of third comp.")),"third.csv")

write.csv2((fourth<-plotLoadings(cyano.splsda, comp = 4, 
                                 method = 'mean', contrib = 'max',   
                                 size.name = 0.7,
                                 ndisplay = 20,
                                 title = "(d) Loadings of fourth comp.")),"fourth.csv")

write.csv2((fifth<-plotLoadings(cyano.splsda, comp = 5, 
                                method = 'mean', contrib = 'max',   
                                size.name = 0.7,
                                ndisplay = 130,
                                title = "(e) Loadings of 5th comp.")),"fifth.csv")



# Relative Abundance of important OTUS determined with the sPLS-DA --------
#import data
taxo <- read.csv2("TAX_Attached.csv",header = T,sep = ";")#Taxonomy of your interesting OTUs
keepotu <- read.csv2("OTUs_intérêts.csv")#Number of your interesting OTUs
table <- dplyr::inner_join(taxo,keepotu,by="Observation_name")
Sept21<-table %>%
  filter(GroupContrib=="September_21th") %>%
  group_by(Order,Phylum) %>%
  count()
Sept14<-table %>%
  filter(GroupContrib=="September_14th") %>%
  group_by(Order) %>%
  count()
Oct12<-table %>%
  filter(GroupContrib=="October_12th") %>%
  group_by(Order) %>%
  count()
Oct19<-table %>%
  filter(GroupContrib=="October_19th") %>%
  group_by(Order) %>%
  count()
Oct26<-table %>%
  filter(GroupContrib=="October_26th") %>%
  group_by(Order) %>%
  count()
Nov02<-table %>%
  filter(GroupContrib=="November_02nd") %>%
  group_by(Order) %>%
  count()

test<-full_join(Sept14,Sept21,by="Order")
test<-full_join(test,Oct12,by="Order")
test<-full_join(test,Oct19,by="Order")
test<-full_join(test,Oct26,by="Order")
test<-full_join(test,Nov02,by="Order")
test[is.na(test)]<-0
test <- as.table(test)
colnames(test)<-c("Order","14_09","21_09","12_10","19_10",
                  "26_10","02_11")
write.csv2(test,"table_histo_splsda.csv")

#converting the last files in percent
#re-import data
data<-read.csv2("table_histo_splsda_%.csv",header = T,sep=";")
colnames(data)<-c("Phylum","Order","14_09","21_09","12_10","19_10","26_10","02_11")
test2 <- data[,2:8]%>% pivot_longer(-Order)
desired_order_genus <- list("Abditibacteriales",	"Blastocatellales",
                            "Pyrinomonadales",	"Vicinamibacterales",	"Corynebacteriales",	"Frankiales",
                            "Gaiellales",	"Micrococcales",	"Microtrichales",	"Propionibacteriales",
                            "Pseudonocardiales",	"Solirubrobacterales",	"Streptomycetales",	
                            "Fimbriimonadales",	"Chitinophagales",	"Cytophagales",	"Flavobacteriales",
                            "Kapabacteriales",	"Rhodothermales",	"Sphingobacteriales",	"Anaerolineales",
                            "Ardenticatenales",	"Kallotenuales",	"Thermomicrobiales",	"Deinococcales",
                            "Desulfobacterales",	"Desulfobulbales",	"Syntrophales",	"Bacillales",	
                            "Brevibacillales",	"Clostridiales",	"Erysipelotrichales",	"Lactobacillales",
                            "Peptostreptococcales-Tissierellales",	"Staphylococcales",	
                            "Thermoactinomycetales",	"Veillonellales-Selenomonadales",	"Fusobacteriales",
                            "Gemmatimonadales",	"Nannocystales",	"Polyangiales",	"Gemmatales",	
                            "Isosphaerales",	"Phycisphaerales",	"Pirellulales",	"Planctomycetales",
                            "Acetobacterales",	"Azospirillales",	"Burkholderiales",	"Caulobacterales",
                            "Holosporales","Methylococcales","Paracaedibacterales",	"Pseudomonadales",
                            "Reyranellales",	"Rhizobiales",	"Rhodobacterales",	"Rhodospirillales",
                            "Rickettsiales",	"Salinisphaerales",	"Sphingomonadales",	"Steroidobacterales",
                            "Xanthomonadales",	"Zavarziniales",	"unknown order")
test2$name_f <- factor(test2$name, levels=c('14_09','21_09','12_10',"19_10",
                                            "26_10",'02_11'))#Ordonate for facet 

#Creation of plot
plot <- ggplot(test2,aes(x=2, y=value, fill=Order)) +
  geom_bar(stat="identity", width=1 )+
  facet_wrap("name_f")+
  theme_minimal()+
  theme(legend.position = "none")+
  labs(x="",y="Relative Abundance %")+
  scale_fill_manual(values = c("#000000",	"#EFBFBF",	"#C48686",	"#7A4A4A",
                               "#CEF0BC",	"#A9C893",	"#7D9A69",	"#5B714B",	"#394B2A",	"#517432",	
                               "#7AAE4E",	"#C7FF8B",	"#B3FF58",	"#7A7A7A",	"#FFD8AA",	"#CCAE88",
                               "#92795B",	"#DCAC64",	"#FFC872",	"#FFBA3F",	"#EFFF99",	"#C2CF7C",
                               "#C4D858",	"#E1FF71",	"#9F9F9F",	"#A8FFCC",	"#61FFA9",	"#2E9E5F",
                               "#A3FFF7",	"#80C4BF",	"#5C9490",	"#397270",	"#329C9C",	"#42D2D2",
                               "#40FFFF",	"#008F8F",	"#004D4D",	"#E2E2E1",	"#543B3B",	"#003C4E",	
                               "#0088AE",	"#B4D7FF",	"#7997BB",	"#4C6681",	"#2C516E",	"#2C74A2",
                               "#F3E8FF",	"#EBCEFF",	"#B39CC8",	"#866F98",	"#573F6A",	"#3C234A",
                               "#2D032D",	"#540540",	"#801661",	"#7A2E62",	"#924A75",	"#A86A8A",
                               "#C96F9F",	"#F885C2",	"#FFAFD9",	"#FFD3E7",	"#FFBCBC",	"#D88F8F",
                               "#5B1414"))+
  coord_polar(theta = "y", start = 0) +
  xlim(0.5, 2.5)  

plot$data$Order<-factor(plot$data$Order, levels = desired_order_genus)
plot 
ggsave("Aundance_histo.svg",device = svg())

#Statistical tests
data<-read.csv2("table_histo_splsda_%.csv",header = T,sep=";")
colnames(data)<-c("Phylum","Order","14_09","X21_09","X12_10","19_10","26_10","X02_11")
data=data[,1:8]
data=data[,-1]
options("digits"=3)
row.names(data)=data[,1]
data=data[,-1]
data=data[,-1]
data=data[,-3]
data=data[,-3] #data formating


interest<-data %>% 
  mutate(sumcol= X21_09+X12_10+X02_11) %>%
  filter(sumcol>0) #suppression of orders missing during the 3 dates

interest=interest[,-4]
chisq.test(interest) -> khideux
khideux #Khi2 test

khideux$residuals ->resid.interest

dataframe_data  %>%
  filter(X21_09 > 1) #Identification of significant different orders: 21_09

dataframe_data  %>%
  filter(X12_10 > 1) #Identification of significant different orders: 12_10

dataframe_data  %>%
  filter(X02_11 > 1) #Identification of significant different orders: 02_11



# Identification of different functions -----------------------------------
# file = A table with column : Otu Name/ Date contribution (one the three important dates) / Affiliation Tax / Features / % relative abondance
read.csv2("yourtable.csv")->tabfun
tabfun=tabfun[,-13]
split.data<-split(tabfun,tabfun$GroupContrib)
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )
for (i in 1:length(split.data)){
  temp<-as.data.frame(split.data[i])
  colnames(temp)[12]<-"X.Abun"
  colnames(temp)[10]<-"Function"
  date<-temp[1,2]
  aggregate(X.Abun~Function,data = temp,FUN = sum) ->temp
  
  
  camembR<- ggplot(temp, aes(x="", y=X.Abun, fill=Function))+
    geom_bar(width = 1, stat = "identity") + coord_polar("y", start=0) +
    blank_theme +
    ggtitle(date) +
    theme(axis.text.x=element_blank()) +
    scale_fill_manual(values = c("#fed976","#ffeca0",
                                 "#d9f0a3","#9ecae0","#c6dbef",
                                 "#deeaf7","#d3b9da","#bababa"))
  
  print(camembR)
}
