#Janis Corona
#Lewis University Summer 2019 Capstone Project
#GEO, gene expression omnibus data on UL and ubiquitous genes associated with UL risk

# uses the following packageNames:
# dplyr, ggplot2, lattice, plotly, UsingR, heatmaply, 
# ComplexHeatmap (BiocManager::install)
# Gviz (bioconductor install), also GenomicRanges to test vignette demos of Gviz
# uncomment the install.packages('packageName'), or 
# install.packages("BiocManager")
# BiocManager::install("Gviz") #select 'a' to update all packages when asked
# BiocManager::install('ComplexHeatmap')

# Analytics with algorithms for machine learning used caret, randomForest, rpart,
# MaSS, and gbm
# Install first if you don't have them, may have to restart R, doesn't take much time to 
# install any of these packages:
# install.packages('caret')
# install.packages('randomForest')
# install.packages('MASS')
# install.packages('gbm')


# Series files use SOFT file downloads:
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68295 
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE13319
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE593
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE23112
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE2724

# Platform files, download full table:
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL570
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL96
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL6480


##############################################################################
#The microarray data is scaled differently for each output of Affymetrix
#this script reads in each microarray of UL/nonUL samples
#then scales two data set outputs to reverse log2, 
#and one data set reverse log2 then reverse avg diff normalization

# set as wd the setwd() in RStudio the original text files to turn into csv


#MAS 5.0 Algorithm (mean scaled to 100)
GSE2724 <- read.csv2("GSE2724_series_matrix.txt", sep='', comment.char='!',
                     na.strings = c('NA',''))#22283X19

#label the UL samples in GSE2724
colnames(GSE2724)[2:8] <- paste(colnames(GSE2724)[2:8], 'UL', sep='')

#RMA normalized results, raw results in CEL (ASCII) txt file encoding for download
GSE23112 <- read.csv2("GSE23112_series_matrix.txt", sep='', comment.char='!',
                      na.strings = c('NA',''))#22283X11

#label the UL in GSE23112
colnames(GSE23112)[7:11] <- paste(colnames(GSE23112)[7:11], 'UL', sep='')

# Log2-transformed normalized signal intensities (quantile normalization),
# raw data is txt supplementary
GSE68295 <- read.csv2("GSE68295_series_matrix.txt", sep='', comment.char='!',
                      na.strings = c('NA',''))#41078X13

#transform GSE68295 into inverse log2 values for same scale
names <- as.matrix(GSE68295$ID_REF, colnames='ID_REF')#41078X1
names2 <- colnames(GSE68295)
GSE68295m <- as.matrix(GSE68295[2:13]) #41078X12
nvrse <- 2^(as.numeric(GSE68295m)) #calculates as one list down each col to next col
nvrse <- matrix(nvrse, nrow=41078, ncol=12, byrow=FALSE) #creates the matrix to bind names
GSE68295 <- cbind(names,nvrse) #41078X13
GSE68295 <- as.data.frame(GSE68295) #41078X13
colnames(GSE68295) <- names2
rm(nvrse);rm(names);rm(names2);rm(GSE68295m) #removes supplementary objects
#now GSE68295 has same scale as other 5 data frames

#label the UL in GSE68295 and remove the sarcoma samples
GSE68295 <- GSE68295[,-c(8:13)] #41078X7
colnames(GSE68295)[5:7] <- paste(colnames(GSE68295)[5:7], 'UL', sep='')

#Microarray Suite 5.0 (MAS 5.0) expression value, not dChip calculated, no supplementary provided
GSE593 <-read.csv2("GSE593_series_matrix.txt", sep='', comment.char='!',
                   skip = 53, na.strings = c('NA',''))#22284X57

#label the UL in GSE593
colnames(GSE593)[2:6] <- paste(colnames(GSE593)[2:6], 'UL', sep='')

# GCOS Gene Chip Operating Software 1.0 signal intensity, 
# raw data supplementary as CEL file download
GSE13319 <- read.table("GSE13319-GPL570_series_matrix.txt", sep='', comment.char='!',
                       header=TRUE, skip = 65, na.strings = c('NA',''))#54675X78

#label the UL in GSE13319
colnames(GSE13319)[2:51] <- paste(colnames(GSE13319)[2:51], 'UL', sep='')

#write the data frame tables to csv export to folder
write.csv(GSE2724, file='GSE2724.csv')
write.csv(GSE764, file='GSE764.csv')
write.csv(GSE23112, file='GSE23112.csv')
write.csv(GSE68295, file='GSE68295.csv')
write.csv(GSE593, file='GSE593.csv')
write.csv(GSE13319, file='GSE13319.csv')

#######################################################################################
# All the ID_REF fields are the same for four tables, but the GSE68295 and GSE764
# tables have different ID_REF values

# Use the GPL-Gene Expression Omnibus platform container of the GSE-GEO series to align
# the ID_Ref by, GPL80:GSE764, GPL6480:GSE68295, GPL96:(GSE593,2724,23112), and
# GPL570:GSE13319, download the full table for each platform, save to original array folder
# the comment.char is '#' in the GPL tables

# there are errors opening the text files in csv,
# so use read.delim with quote="" to fix, because of the '///' and '//'

GPL570 <- read.delim('GPL570-55999.txt', sep='\t', quote="", header=TRUE,
                     comment.char='#', na.strings=c('','NA')) #54675X16
colnames(GPL570)
# [1] "ID"                               "GB_ACC"                          
# [3] "SPOT_ID"                          "Species.Scientific.Name"         
# [5] "Annotation.Date"                  "Sequence.Type"                   
# [7] "Sequence.Source"                  "Target.Description"              
# [9] "Representative.Public.ID"         "Gene.Title"                      
# [11] "Gene.Symbol"                      "ENTREZ_GENE_ID"                  
# [13] "RefSeq.Transcript.ID"             "Gene.Ontology.Biological.Process"
# [15] "Gene.Ontology.Cellular.Component" "Gene.Ontology.Molecular.Function"

GPL6480 <- read.delim('GPL6480-9577.txt', sep='\t', quote="", header=TRUE,
                      comment.char='#', na.strings=c('','NA')) #41108X17
colnames(GPL6480)
# [1] "ID"                   "SPOT_ID"              "CONTROL_TYPE"        
# [4] "REFSEQ"               "GB_ACC"               "GENE"                
# [7] "GENE_SYMBOL"          "GENE_NAME"            "UNIGENE_ID"          
# [10] "ENSEMBL_ID"           "TIGR_ID"              "ACCESSION_STRING"    
# [13] "CHROMOSOMAL_LOCATION" "CYTOBAND"             "DESCRIPTION"         
# [16] "GO_ID"                "SEQUENCE"   

GPL96 <- read.delim('GPL96-57554.txt', sep='\t', quote="", header=TRUE,
                    comment.char='#', na.strings=c('','NA')) #22283X16
colnames(GPL96)
# [1] "ID"                               "GB_ACC"                          
# [3] "SPOT_ID"                          "Species.Scientific.Name"         
# [5] "Annotation.Date"                  "Sequence.Type"                   
# [7] "Sequence.Source"                  "Target.Description"              
# [9] "Representative.Public.ID"         "Gene.Title"                      
# [11] "Gene.Symbol"                      "ENTREZ_GENE_ID"                  
# [13] "RefSeq.Transcript.ID"             "Gene.Ontology.Biological.Process"
# [15] "Gene.Ontology.Cellular.Component" "Gene.Ontology.Molecular.Function"

#########################################################################################

# merge the platforms to their series

#Microarray Suite 5.0 (MAS 5.0) expression value, not dChip calculated, no supplementary provided
GSE593_platform <- merge(GPL96, GSE593, by.x='ID', by.y="ID_REF")#22283X26

#RMA normalized results, raw results in CEL (ASCII) txt file encoding for download
GSE23112_platform <- merge(GPL96, GSE23112, by.x='ID', by.y="ID_REF")#22283X26

#MAS 5.0 Algorithm (mean scaled to 100)
GSE2724_platform <- merge(GPL96, GSE2724, by.x="ID", by.y="ID_REF")#22283X34

# GCOS Gene Chip Operating Software 1.0 signal intensity, 
# raw data supplementary as CEL file download
GSE13319_platform <- merge(GPL570, GSE13319, by.x="ID", by.y="ID_REF")#54675X93

# Log2-transformed GSE68295 normalized signal intensities were inversed earlier
GSE68295_platform <- merge(GPL6480, GSE68295, by.x='ID', by.y="ID_REF")#41078X23

#GSE68295_platform is single value for GENE_SYMBOL field
#The ENTREZ_GENE_ID field is the NCBI gene ID genenames.org for GSE764_platform
#The GENE field of GSE68295_platform is the NCBI gene ID
#The ENTREZ_GENE_ID field is the NCBI gene ID for GSE2724_platform, multiple entries
#The same with GSE23112_platform, GSE593_platform, and GSE13319_platform

Table2724 <- GSE2724_platform[,-c(1:11,13:16)] #22283X19
Table13319 <- GSE13319_platform[,-c(1:11,13:16)]#54675X78
Table23112 <- GSE23112_platform[,-c(1:11,13:16)]#22283X11
Table593 <- GSE593_platform[,-c(1:11,13:16)]#22283X11
Table68295 <- GSE68295_platform[,-c(1:5,7:17)] #41078X7

T2724 <- na.omit(Table2724) #20973X19
T13319 <- na.omit(Table13319) #44134X78
T23112 <- na.omit(Table23112)#201973X11
T593 <- na.omit(Table593) # 20973X11
T68295 <- na.omit(Table68295) #30396X7

#split the factors into first listed and drop the other list NCBI genes of Entrez field
Symbol <- T2724$ENTREZ_GENE_ID
Symbol <- strsplit(as.character(Symbol), '[///]' ) #splits field by '///'
Symbol <- lapply(Symbol, '[',1) #this takes the 1st item in each field row
length(Symbol) #20973
length(unique(Symbol)) # 13160
T2724$ENTREZ_GENE_ID <- as.character(Symbol) #this changes field to one value

Symbol <- T13319$ENTREZ_GENE_ID
Symbol <- strsplit(as.character(Symbol), '[///]' ) #splits field by '///'
Symbol <- lapply(Symbol, '[',1) #this takes the 1st item in each field row
length(Symbol) #44134
length(unique(Symbol)) # 21712 unique genes
T13319$ENTREZ_GENE_ID <- as.character(Symbol) #this changes field to one value

Symbol <- T23112$ENTREZ_GENE_ID
Symbol <- strsplit(as.character(Symbol), '[///]' ) #splits field by '///'
Symbol <- lapply(Symbol, '[',1) #this takes the 1st item in each field row
length(Symbol) #20973 genes
length(unique(Symbol)) # 13160 unique genes
T23112$ENTREZ_GENE_ID <- as.character(Symbol) #this changes field to one value

Symbol <- T593$ENTREZ_GENE_ID
Symbol <- strsplit(as.character(Symbol), '[///]' ) #splits field by '///'
Symbol <- lapply(Symbol, '[',1) #this takes the 1st item in each field row
length(Symbol) #20973
length(unique(Symbol)) # 13160 unique genes
T593$ENTREZ_GENE_ID <- as.character(Symbol) #this changes field to one value

# T68295 only has one entry for GENE field, no need to split and keep 1st in list
# Checking all the NCBI genes for TNRC6B:23112, BET1L:51272, CYTH4:27128, FASN:2194,
# HMGA2:8091, and CCDC57:284001,
# all the codes for those genes were in those data sets T593:T23112 to merge on,
# excluding GSE764 for now, makes 121 UL and non UL samples, of 70 UL and 51 non-UL

#Save meta fields of GSE68295
GSE_meta <- GSE68295_platform[,c(1:17)]#41078X17

#write file out for meta fields to use if needed from largest array GSE68295_platform
write.csv(GSE_meta, file='GSE_array_meta.csv',row.names=FALSE)

# Start merging the five data sets, do this the first time. Once folder has files, then
# just read in the table with csv file name, 
# takes a long while to do these merges df4 and df5

df1 <- merge(META, T68295,  by.x='GENE', by.y='GENE', no.dups=TRUE) #31137X8
write.csv(df1, 'mrg1.csv', row.names=FALSE)#3.9 MB
d <- read.csv('mrg1.csv', sep=',', header=TRUE, na.strings=c('','NA'))
d[3:8] <- round(d[3:8],1)

df2 <- merge(d, T13319, by.x='GENE', by.y='ENTREZ_GENE_ID', no.dups=TRUE) #71725X85
write.csv(df2, 'mrg2.csv', row.names=FALSE) #50.5 MB
d2 <- read.csv('mrg2.csv', sep=',', header=TRUE, na.strings=c('','NA'))
d2[3:85] <- round(d2[3:85],1)

df3 <- merge(d2, T2724, by.x='GENE', by.y='ENTREZ_GENE_ID', no.dups=TRUE) #122243X103
write.csv(df3,'mrg3.csv', row.names=FALSE)#66.2 MB
d3 <- read.csv('mrg3.csv', sep=',', header=TRUE, na.strings=c('','NA'))
d3[3:103] <- round(d3[3:103],1)

df4 <- merge(d3, T593, by.x='GENE', by.y='ENTREZ_GENE_ID', no.dups=TRUE) #405837X113
write.csv(df4, 'mrg4.csv', row.names=FALSE)#233.3 MB
d4 <- read.csv('mrg4.csv', sep=',', header=TRUE, na.strings=c('','NA'))
d4[3:113] <- round(d4[3:113],1)
write.csv(d4, 'mrg4.csv', row.names=FALSE)# 225.6 MB


df5 <- merge(d4, T23112, by.x='GENE', by.y='ENTREZ_GENE_ID', 
             no.dups=TRUE) #1954853X123
write.csv(df5, 'merged_121_samples.csv', row.names=FALSE) #1.3 GB
d5 <- read.csv('merged_121_samples.csv', sep=',', header=TRUE, na.strings=c('','NA'))
d5[3:123] <- round(d5[3:123],1)
write.csv(d5, 'mrg5.csv', row.names=FALSE)# 1.1 GB

#####################################################################################
#####################################################################################
# This code used dplyr package to group by GENEs in the combined data set and 
# summarise each group by mean for the group in each sample
install.packages('dplyr', suppressMessages=TRUE, suppressWarnings=TRUE)
library(dplyr)

# Read in the last data frame of all merged genes
df <- read.csv('mrg5.csv', sep=',', header=TRUE, na.strings=c('','NA')) # 1,954,853 X 123
META <- read.csv('META.csv',sep=',', header=TRUE, na.strings=c('','NA')) #19,634X3
META <- META[,2:3] # Gene and CYTOBAND

#Create a grouped table by NCBI genes and the mean of each sample for that one gene
samples <- colnames(df[3:123]) #GSE sample IDs
samples <- as.vector(samples)
gene <- df %>% group_by(GENE) %>% summarise_at(samples, mean, na.rm = TRUE)# 12,173 X 122
gene1 <- df %>% group_by(GENE) %>% summarise(Counts=n()) # 12,173 X 2

Genes <- merge(gene1,gene, by.x='GENE', by.y='GENE') # 12,173 X123

#use previous META data of genes and cytoband
GenesChromosMeans <- merge(META, Genes, by.x='GENE', by.y='GENE') #12238X124
meta <- read.csv('GSE_array_meta.csv', sep=',', 
                 na.strings=c('','NA'), header=TRUE) #41078X17

#now add the gene symbol from the meta file made earlier
meta1 <- meta[,c(6:7)] #41078X2
meta2 <- na.omit(meta1) #30936X2
GeneSymbolMeans <- merge(meta2, GenesChromosMeans, by.x='GENE', 
                         by.y='GENE') #20377X125
samples2 <- colnames(GeneSymbolMeans[5:125])#1:121
samples2 <- as.vector(samples)#1:121
gene2 <- GeneSymbolMeans %>% group_by(GENE) %>% summarise_at(samples2, mean, 
                                                             na.rm = TRUE)#12173X122
gene3 <- GeneSymbolMeans %>% group_by(GENE) %>% summarise(Counts=n()) #12173X2
gene4 <- GeneSymbolMeans[,1:2] #20377X2
gene5 <- gene4[!duplicated(gene4$GENE),] #12173X2
Genes2 <- merge(gene5, gene3, by.x='GENE', by.y='GENE') #12173X3
Genes3 <- merge(Genes2,gene2, by.x='GENE', by.y='GENE') #12173X124
gene6 <- GenesChromosMeans[,1:2] #12238X2
gene7 <- gene6[!duplicated(gene6$GENE),] #12173X2

Genes4 <- merge(gene7, Genes3, by.x='GENE', by.y='GENE') #12173X125
#check all six top genes are in this combined data of 121 gene sample means per gene
grep('23112', as.character(Genes4$GENE))#TNRC6B
#[1] 7382

grep('51272', as.character(Genes4$GENE))#BET1L
#[1] 8759

grep('27128', as.character(Genes4$GENE))#CYTH4
#[1] 8226

grep('2194', as.character(Genes4$GENE))#FASN
#[1]  1442 12090

grep('8091', as.character(Genes4$GENE))#HMGA2
#[1] 4953

grep('284001', as.character(Genes4$GENE))#CCDC57
#[1] 11965

# There are duplicate entries for FASN and shouldn't be, when checking, this is because
# two genes have '2194' in its GENE ID, '442194' and '2194', so ignore, its is valid

write.csv(Genes4, 'DE_means_Per_Gene_Chr.csv', row.names=FALSE) #9.9 MB

#######################################################################################
#####################################################################################
# This code grabs the cytoband chromosome location of the top ubiquitous genes from
# the larger merged data sets (except for GSE764)

# Read in the last data frame of all merged genes,
# these data sets are means per gene with cytoband location and counts per gene
# when generating the means per gene

# this table is the five data sets that have all six genes
df <- read.csv('DE_means_Per_Gene_Chr.csv', sep=',', 
               header=TRUE, na.strings=c('','NA')) 
#12173X125, 
#five data sets excluding GSE764

head(colnames(df),8)
# [1] "GENE"         "CYTOBAND"     "GENE_SYMBOL"  "Counts"       "GSM1667144"   "GSM1667145"  
# [7] "GSM1667146"   "GSM1667147UL"

unique(df$GENE_SYMBOL) # 12,173 Levels
###################################################################
#                                                                 #
# This code extracts the cytobands of only top six gene locations #
#                                                                 #
###################################################################

# https://www.genenames.org/
T <- grep('22q13.1', df$CYTOBAND) #56 gene findings for chromosome location of TNRC6B gene

# The other ubiquitous genes are FASN, CYTH4, BET1L, HGMA2, and CCDC57
# the chromosomal locations:
# FASN: 17q25.3
# CYTH4: 22q13.1, same as for TNRC6B
# BET1L: 11p15.5
# HMGA2: 12q14.3
# CCDC57: 17q25.3

F <- grep('17q25.3', df$CYTOBAND) #69 FASN results, same location as CCDC57
CY4 <- grep('22q13.1', df$CYTOBAND) #56 CYTH4 results, same location as TNRC6B
B <- grep('11p15.5', df$CYTOBAND) #49 BET1L results
H <- grep('12q14.3', df$CYTOBAND) #9 HMGA2 results
C57 <- grep('17q25.3', df$CYTOBAND) #69 CCDC57 results, same location as FASN 

# two of the loci are the same, so use one of either T or CY4, and C57 or F
DATA <- df[c(T,F,B,H),] #183X125, 
# compare to 112X128 of only genes not gene loci of ubiq3 data or UBIQ df

# The unique genes are 42 for unigene, gene, and gene symbol

unique(DATA$GENE_SYMBOL)
# [1] ADSL      ATF4      CSNK1E    H1F0      KCNJ4     LGALS1    LGALS2    MFNG      MGAT3     PDGFB    
# [11] POLR2F    RAC2      RPL3      SOX10     SSTR3     PLA2G6    GALR3     CACNA1I   SYNGR1    GRAP2    
# [21] PICK1     GTPBP1    APOBEC3B  JOSD1     DNAL4     TAB1      DDX17     KDELR3    TRIOBP    CDC42EP1 
# [31] DMC1      TNRC6B    GCAT      CBX6      NPTXR     CBX7      SLC16A8   SH3BP1    MAFF      CBY1     
# [41] SUN2      TMEM184B  GGA1      CYTH4     APOBEC3C  SGSM3     CBX2    EIF3L     SMCR7L    TOMM22   
# [51] MKL1      APOBEC3G  NOL12     BAIAP2L2  MICALL1   APOBEC3F  BIRC5     ARHGDIA   CD7       CSNK1D   
# [61] SLC25A10  FASN      GAA       GCGR      GPS1      FOXK2     LGALS3BP  MAFG      NPTX1     P4HB     
# [71] PDE6G     PYCR1     PCYT2     RAC3      RFNG      MRPL12    SECTM1    SGSH      TBCD      TIMP2    
# [81] TK1       CBX4      DNAH17    SOCS3     SYNGR2    HGS       CYTH1     PGS1      AATK      EIF4A3   
# [91] BAIAP2    SEPT9     RAB40B    TMC6      AZI1      FSCN2     NARF      DCXR      SIRT7     CCDC40   
# [101] NPLOC4    WDR45L    CBX8      BAHCC1    USP36     DUS1L     FN3K      ENGASE    ASPSCR1   CARD14   
# [111] C17orf62  CHMP6     FN3KRP    C17orf101 ZNF750    C17orf70  CBX2      C1QTNF1   SLC38A10  CANT1    
# [121] TBC1D16   B3GNTL1   STRA13    CCDC57    C17orf90  AP2A2     ASCL2     CD81      CD151     CTSD     
# [131] DRD4      HRAS      INS       IRF7      LSP1      MUC2      MUC5AC    MUC6      POLR2L    PSMD13   
# [141] ASPSCR1      MRPL23    RPLP2     SCT       TALDO1    TH        TSPAN4    TNNI2     TNNT3     RASSF7   
# [151] IFITM1    BRSK2     TSPAN32   TSSC4     IFITM3    DEAF1     IFITM2    PKP3      SIRT3     C11orf21 
# [161] IGF2-AS   BET1L     CEND1     CDHR5     TOLLIP    PIDD      PNPLA2    SIGIRR    EPS8L2    SLC25A22 
# [171] ATHL1     PTDSS2    TMEM80    MUC5B     GNS       HMGA2     WIF1      IRAK3     TBC1D30   SOCS3    
# [181] LEMD3     TMBIM4    LLPH     
# 12173 Levels: A1CF A2M A4GALT A4GNT AAAS AACS AADAC AAGAB AAK1 AAMP AANAT AARS AASDHPPT AASS AATF ... ZZZ3

write.csv(DATA, 'chr_loci_top_genes.csv', row.names=FALSE)#151.4 kb and 183X125

##########################################################################
# This section code will organize the samples so that the UL samples are at the end,
# and also so there are separate samples by GENE for UL and nonUL

big_df <- read.csv('chr_loci_top_genes.csv', sep=',', header=TRUE)#183X125
ul <- grep('UL', colnames(big_df))#1:70
UL <- big_df[,ul] # only ul no meta fields, 183X70
nonUL_df <- big_df[,-ul] #nonUL with meta data table #183X55
non_df <- nonUL_df[,-c(4)] #nonUL w/o original counts field, 183X54
ul_df <- big_df[,c(1:3,ul)] #ul data table w/o original counts field, 183X73
data_5 <- cbind(non_df, UL) #183X124, no orig. counts of each sample
#write out the csv files for gene expression with UL, nonUL, and Meta_UL_non


write.csv(data_5, 'Data_5_ready.csv', row.names=FALSE) #151 kb, 183X124
write.csv(non_df, 'Data_5_nonUL_only.csv', row.names=FALSE) #65.7 kb, 183X54
write.csv(ul_df, 'Data_5_ul_only.csv', row.names=FALSE) #90.4 kb, 183X73
############################################################################################

# Explore the means counts of the 5 data sets again for counts, maybe with development using
# bootstrapping

DE_5 <- read.csv('DE_means_Per_Gene_Chr.csv', sep=',', header=TRUE, na.strings=c('','NA'))
T <- grep('TNRC6B', DE_5$GENE_SYMBOL)
B <- grep('BET1L', DE_5$GENE_SYMBOL)
CY <- grep('CYTH4', DE_5$GENE_SYMBOL)
CC <- grep('CCDC57', DE_5$GENE_SYMBOL)
F <- grep('FASN', DE_5$GENE_SYMBOL)
H <- grep('HMGA2', DE_5$GENE_SYMBOL)

DE <- DE_5[c(T,B,CY,CC,F,H),]
DE[,c(3,4)]
# GENE_SYMBOL Counts
# 7382       TNRC6B      4
# 8759        BET1L      2
# 8226        CYTH4      2
# 11965      CCDC57      6
# 1442         FASN      1
# 4953        HMGA2      2
sum(DE$Counts) #17 combined observations, and some were combined means in making data set

All <- read.csv('Data_5_ready.csv', sep=',', header=TRUE, na.strings=c('','NA'))
UL <- read.csv('Data_5_ul_only.csv', sep=',', header=TRUE, na.strings=c('','NA'))
nonUL <- read.csv('Data_5_nonul_only.csv', sep=',', header=TRUE, na.strings=c('','NA'))

#keep only the GENE_SYMBOL field to map Genes to Samples

All1 <- All[,-c(1:2)] #183X122
UL1 <- UL[,-c(1:2)] #183X71
nonUL1 <- nonUL[,-c(1:2)] #183X52

#####################################################################################

# Extract meta fields from different GPL file, GPL6480, bc GSE68295 has all genes
meta <- read.delim('GPL6480-9577.txt', sep='\t', quote="", header=TRUE,
                   comment.char='#', na.strings=c('','NA')) #41108X17

colnames(meta)
# [1] "GENE_SYMBOL"          "ID"                   "SPOT_ID"             
# [4] "CONTROL_TYPE"         "REFSEQ"               "GB_ACC"              
# [7] "GENE"                 "GENE_NAME"            "UNIGENE_ID"          
# [10] "ENSEMBL_ID"           "TIGR_ID"              "ACCESSION_STRING"    
# [13] "CHROMOSOMAL_LOCATION" "CYTOBAND"             "DESCRIPTION"         
# [16] "GO_ID"                "SEQUENCE"    

#########################################################################################
# This code will add 3 fields to the meta table and drop fields not needed
# Split the CHROMOSOMAL_LOCATION field in the meta table by chr, and mb,
# then within the mb by start and end 

class(meta$CHROMOSOMAL_LOCATION)#factor

chr <- strsplit(as.character(meta$CHROMOSOMAL_LOCATION), ':') #list of 41108
chromosome <- lapply(chr, '[',1) #this takes the 1st item in each field row
mb <- lapply(chr, '[',2) #this takes the 2nd item in each field row

st <- strsplit(as.character(mb), '-')
start <- lapply(st, '[',1)
end <- lapply(st,'[',2)

CHROMOSOME <- as.matrix(as.character(chromosome))
START <- as.matrix(as.character(start))
END <- as.matrix(as.character(end))

meta <- cbind(meta,CHROMOSOME,START,END) #41108X20
colnames(meta)
# [1] "ID"                   "SPOT_ID"              "CONTROL_TYPE"        
# [4] "REFSEQ"               "GB_ACC"               "GENE"                
# [7] "GENE_SYMBOL"          "GENE_NAME"            "UNIGENE_ID"          
# [10] "ENSEMBL_ID"           "TIGR_ID"              "ACCESSION_STRING"    
# [13] "CHROMOSOMAL_LOCATION" "CYTOBAND"             "DESCRIPTION"         
# [16] "GO_ID"                "SEQUENCE"             "CHROMOSOME"          
# [19] "START"                "END"    

meta <- meta[,-c(1:6, 8:9,11:16)]
colnames(meta)
# [1] "GENE_SYMBOL" "ENSEMBL_ID"  "SEQUENCE"    "CHROMOSOME"  "START"       "END" 

META <- meta[unique(meta$GENE_SYMBOL),] #19596X6

#check the ubiquitous genes are in this meta file:
# grep('CYTH4', META$GENE_SYMBOL)
# [1] 15976
# > grep('BET1L', META$GENE_SYMBOL)
# [1] 16415
# > grep('TNRC6B', META$GENE_SYMBOL)
# [1] 2619 3251
# > grep('CCDC57', META$GENE_SYMBOL)
# [1]  6125 12335 19343
# > grep('FASN', META$GENE_SYMBOL)
# [1] 17302
# > grep('HMGA2', META$GENE_SYMBOL)
# [1] 5760

# They are all in the META table

# write out to csv file to add to DE data for Gviz
write.csv(META, 'meta_gviz.csv', row.names=FALSE) #2.2 MB, 19596X6

#########################################################################################
# this code merges meta fields for Gviz package to the All data table of UL and nonUL
# genes by chromosome region of ubiquitous genes in population studies

META <- read.csv('meta_gviz.csv', sep=',', header=TRUE, na.strings=c('','NA'))#19596X6
meta <- META[,-3] # remove sequence field for merge
All <- read.csv('Data_5_ready.csv', sep=',', header=TRUE, na.strings=c('','NA'))

#check all ubiquitous genes are in this file:
# grep('CYTH4', All$GENE_SYMBOL)
# [1] 44
# > grep('CCDC57', All$GENE_SYMBOL)
# [1] 124
# > grep('TNRC6B', All$GENE_SYMBOL)
# [1] 32
# > grep('BET1L', All$GENE_SYMBOL)
# [1] 162
# > grep('FASN', All$GENE_SYMBOL)
# [1] 62
# > grep('HMGA2', All$GENE_SYMBOL)
# [1] 176

All <- All[,-c(1:2)] # 183X122, removes fields not GENE_SYMBOL or sample of UL or nonUL


All_gviz <- merge(All, meta,  by.x='GENE_SYMBOL', by.y='GENE_SYMBOL',
                  no.dups=TRUE) #196X126
# Check all ubiquitous genes are in this data:

grep('CYTH4', All_gviz$GENE_SYMBOL)
grep('CCDC57', All_gviz$GENE_SYMBOL)
grep('TNRC6B', All_gviz$GENE_SYMBOL)
grep('BET1L', All_gviz$GENE_SYMBOL)
grep('FASN', All_gviz$GENE_SYMBOL)
grep('HMGA2', All_gviz$GENE_SYMBOL)

# All genes are in this data set of ubiquitous genes and loci for those genes

colnames(All_gviz) <- tolower(colnames(All_gviz)) #change colnames all lowercase
All_gviz <- All_gviz[,c(1,123:126,2:122)] #reorder columns by meta then samples
all <- na.omit(All_gviz) #176X126, removes ensemble IDs with NAs

#check all ubiquitous genes in table:
grep('CYTH4', all$gene_symbol) #[1] 46
grep('CCDC57', all$gene_symbol)# [1] 33 34
grep('TNRC6B', all$gene_symbol)# [1] 167 168
grep('BET1L', all$gene_symbol)# [1] 20
grep('FASN', all$gene_symbol)# [1] 61
grep('HMGA2', all$gene_symbol)# [1] 81

# all genes are in this data, write out to csv

write.csv(all, 'all_gviz.csv', row.names=FALSE) #152.9 kb, 176X126, 5 meta, 121 samples
#######################################################################################
# this code will change the ensemble id to 'transcript' and gene_symbol to 'symbol'
# also adds a 'width' field for end-start+1 to include start 
# search ucsc or ensembl for strand direction of transcript field

all <- read.csv('all_gviz.csv', sep=',', header=TRUE, 
                na.strings=c('','NA')) #176X126
colnames(all)
#  [1] "gene_symbol"  "ensembl_id"   "chromosome"   "start"        "end"...

colnames(all)[c(1,2)] <- c('symbol','transcript')
colnames(all)
#   [1] "symbol"       "transcript"   "chromosome"   "start"        "end"   ...

width <- abs(all$start-all$end)+1 #num[1:176]

all <- cbind(all, width) #176X127
all <- all[,c(1:5,127, 6:126)]
head(colnames(all))
# [1] "symbol"     "transcript" "chromosome" "start"      "end"        "width"

all1 <- all[order(all$width, decreasing=TRUE),] #176X127, 
# ordered by width largest to smallest

write.csv(all1, 'ub_genes_gviz.csv', row.names=FALSE) #153.5 kb, 176X127

###############################################################################################
# to search ensembl.org, go to BioMart -> New -> Ensembl 96 -> Human Genes (GRCh38.12)
# select all to display, copy the data from new tab, ctrl+a, ctrl+c, then paste into
# notepad as csv file, shows the ensembl transcript and ID next to each other,
# save as: ensembl_generated_id.csv

ensembl <- read.csv('ensembl_generated_id.csv', sep=',', 
                    header=TRUE, na.strings=c('','NA')) #has the gene = transcript,229428X2

###############################################################################################

genes <- read.csv('ub_genes_gviz.csv', sep=',', header=TRUE, 
                  na.strings=c('','NA')) #176X127
head(colnames(genes)) #this is the ensembl gene transcript ID, chromosome location,
# all ubiquitous genes listed, and loci of those genes, not the strand direction

# [1] "symbol"     "transcript" "chromosome" "start"      "end"        "width" 


#########################################################################################
#########################################################################################

genes <- read.csv('ub_genes_gviz.csv', sep=',', header=TRUE, 
                  na.strings=c('','NA')) #176X127

CYTH4 <- grep('CYTH4', genes$symbol) #[1] 75
CCDC57 <- grep('CCDC57', genes$symbol)# [1] 64 65
TNRC6B <- grep('TNRC6B', genes$symbol)# [1] 168 169
BET1L <- grep('BET1L', genes$symbol)# [1] 51
FASN <- grep('FASN', genes$symbol)# [1] 20
HMGA2 <- grep('HMGA2', genes$symbol)# [1] 101

#check out the duplicate gene entries
genes[64,] == genes[65,] #TRUE for all fields except transcript, start, and end
genes[168,] == genes[169,]#TRUE for all fields except transcript, start, and end

# read in the trascript ID from ENSEMBL downloaded earlier from ENSEMBLE
# to search ensembl.org, go to BioMart -> New -> Ensembl 96 -> Human Genes (GRCh38.12)
# select all to display, copy the data from new tab, ctrl+a, ctrl+c, then paste into
# notepad as csv file, shows the ensembl transcript and ID next to each other,
# save as: ensembl_generated_id.csv

ensembl <- read.csv('ensembl_generated_id.csv', sep=',', 
                    header=TRUE, na.strings=c('','NA')) #has the gene = transcript,229428X2
colnames(ensembl)
# [1] "Gene.stable.ID"       "Transcript.stable.ID"

#add the gene ID for ensemble by transcript
df <- merge(genes, ensembl, by.x='transcript', by.y='Transcript.stable.ID') #149X128

df$symbol #all ubiquitous genes are included

# reorder the fields to have the meta at front and the samples at the end
df <- df[,c(1,128,2:6,7:127)]
colnames(df)
# [1] "transcript"     "Gene.stable.ID" "symbol"         "chromosome"     "start"         
# [6] "end"            "width"          "gsm1667144"     "gsm1667145"     "gsm1667146"    
# [11] "gsm336252"      "gsm336253"      "gsm336254"      "gsm336255"      "gsm336256"     
# [16] "gsm336257"      "gsm336258"      "gsm336259"      "gsm336260"      "gsm336261"     
# [21] "gsm336262"      "gsm336263"      "gsm336264"      "gsm336265"      "gsm336266"     
# [26] "gsm336267"      "gsm336268"      "gsm336269"      "gsm336270"      "gsm336271"     
# [31] "gsm336272"      "gsm336273"      "gsm336274"      "gsm336275"      "gsm336276"     
# [36] "gsm336277"      "gsm336278"      "gsm52661"       "gsm52662"       "gsm52663"      
# [41] "gsm52664"       "gsm52665"       "gsm52666"       "gsm52667"       "gsm52668"      
# [46] "gsm52669"       "gsm52670"       "gsm52671"       "gsm9098"        "gsm9099"       
# [51] "gsm9100"        "gsm9101"        "gsm9102"        "gsm569424"      "gsm569425"     
# [56] "gsm569426"      "gsm569427"      "gsm569428"      "gsm1667147ul"   "gsm1667148ul"  
# [61] "gsm1667149ul"   "gsm336202ul"    "gsm336203ul"    "gsm336204ul"    "gsm336205ul"   
# [66] "gsm336206ul"    "gsm336207ul"    "gsm336208ul"    "gsm336209ul"    "gsm336210ul"   
# [71] "gsm336211ul"    "gsm336212ul"    "gsm336213ul"    "gsm336214ul"    "gsm336215ul"   
# [76] "gsm336216ul"    "gsm336217ul"    "gsm336218ul"    "gsm336219ul"    "gsm336220ul"   
# [81] "gsm336221ul"    "gsm336222ul"    "gsm336223ul"    "gsm336224ul"    "gsm336225ul"   
# [86] "gsm336226ul"    "gsm336227ul"    "gsm336228ul"    "gsm336229ul"    "gsm336230ul"   
# [91] "gsm336231ul"    "gsm336232ul"    "gsm336233ul"    "gsm336234ul"    "gsm336235ul"   
# [96] "gsm336236ul"    "gsm336237ul"    "gsm336238ul"    "gsm336239ul"    "gsm336240ul"   
# [101] "gsm336241ul"    "gsm336242ul"    "gsm336243ul"    "gsm336244ul"    "gsm336245ul"   
# [106] "gsm336246ul"    "gsm336247ul"    "gsm336248ul"    "gsm336249ul"    "gsm336250ul"   
# [111] "gsm336251ul"    "gsm38689ul"     "gsm38690ul"     "gsm38691ul"     "gsm38692ul"    
# [116] "gsm38693ul"     "gsm38694ul"     "gsm38695ul"     "gsm9093ul"      "gsm9094ul"     
# [121] "gsm9095ul"      "gsm9096ul"      "gsm9097ul"      "gsm569429ul"    "gsm569430ul"   
# [126] "gsm569431ul"    "gsm569432ul"    "gsm569433ul"   

colnames(df)[2] <- 'ensembl'

write.csv(df, 'ub_genes_ensembl.csv', row.names=FALSE) #134 kb

#clear environment and read in the latest data with ubiquitous genes and ensembl ID
df <- read.csv('ub_genes_ensembl.csv', sep=',', header=TRUE, 
               na.strings=c('','NA'))#149X128

#############################################################################################

df <- read.csv('ub_genes_ensembl.csv', sep=',', header=TRUE, 
               na.strings=c('','NA'))#149X128

#transform all the samples to log2 scale
df[,c(8:128)] <- log2(df[,c(8:128)])

#round to two decimal places
df[8:128] <- round(df[8:128],2) #149X128

# there is a file to download from ensembl.org -> Biomart -> 
# -> Ensembl Genes 96 -> Human genes (GRCh38.p12) -> select Structures -> Gene Stable ID, 
# Transcript Stable ID, Strand,# ,Chromosome/Scaffold name, 
# Gene Start (bp). Gene end (bp), Gene Name
# export 'Results' as csv format for 'all' entries
# saved as 'mart_export.txt' 14.1 MB in size

Strand <- read.csv('mart_export.txt', sep=',', header=TRUE, 
                   na.strings=c('','NA')) #229428X7
colnames(Strand)
# [1] "Gene.stable.ID"           "Transcript.stable.ID"     "Strand"                  
# [4] "Gene.end..bp."            "Gene.start..bp."          "Chromosome.scaffold.name"
# [7] "Gene.name"

#change the strand from 1 and -1 to + and - for gviz package
Strand$Strand <- gsub('-1','-', as.character(Strand$Strand))
Strand$Strand <- gsub('1','+', as.character(Strand$Strand))

DF <- merge(Strand,df, by.x='Transcript.stable.ID', by.y='transcript') #149X134
colnames(DF)[1:13]
# [1] "Transcript.stable.ID"     "Gene.stable.ID"           "Strand"                  
# [4] "Gene.end..bp."            "Gene.start..bp."          "Chromosome.scaffold.name"
# [7] "Gene.name"                "ensembl"                  "symbol"                  
# [10] "chromosome"               "start"                    "end"                     
# [13] "width"   ...   

DF <- DF[,-c(2,6,7,11:13)] #149X128
colnames(DF)[1:7]
# [1] "Transcript.stable.ID" "Strand"               "Gene.end..bp."        "Gene.start..bp."     
# [5] "ensembl"              "symbol"               "chromosome"    

colnames(DF)[1:7] <- c('transcript', 'strand','end', 'start', 'gene','symbol', 'chromosome')
#missing 'feature' and 'exon' fields from download

width <- abs(DF$start-DF$end)+1 #149 in length
width <- as.data.frame(width)
colnames(width) <- 'width'

df <- cbind(DF,width)#149X129
DF <- df[,c(1:7,129,8:128)]
df <- DF[,c(7,4,3,8,2,5,1,6,9:129)] #149X129

write.csv(df, 'ub_genes_ensembl_gviz.csv', row.names=FALSE) #100.2 kb

all <- df[,c(9:129)]#149X121
nonUL <- df[,c(9:59)] #149X51
UL <- df[,c(60:129)] #149X70
names <- as.character(df$symbol)

nonUL <- as.matrix(nonUL)
row.names(nonUL) <- names

UL <- as.matrix(UL)
row.names(UL) <- names

all <- as.matrix(all)
row.names(all) <- names


write.csv(all, 'ub_genes_ensembl_matrix.csv', row.names=names) #89.5 kb
write.csv(UL, 'ub_UL_genes_ensembl_matrix.csv', row.names=names) #52.4 kb
write.csv(nonUL, 'ub_nonUL_genes_ensembl_matrix.csv', row.names=names) #38.4 kb

# the row names do not get saved as such when saving as matrix or data frame
# use rownames=1 when reading in as csv to keep row names instead of adding 'X' field
##################################################################################

# This script uses plotting of ggplot2 to analyze data features
# Also, adds categorical and differential gene expression data between nonUL-UL means
# shows pairwise comparison plots in lattice
# simple plots, boxplots of genes in data are also shown

Data <- read.csv('ub_genes_ensembl_gviz.csv', sep=',', header=TRUE,
                 na.strings=c('','NA')) #149X129

# make subset of data to analyze

# The '!duplicated()' doesn't always omit duplicates, has to be done twice
d <- Data[!duplicated(Data),] #139X129
d <- d[!duplicated(d),c(8:129)] # only gene symbol and samples subset, 139X122

library(dplyr)

# counts each gene to make sure only one of each gene in table
Counts <- d %>% group_by(symbol) %>% summarise(Counts=n()) #130X2, 122:130 have 2 counts each

# This omits the duplicates, so counts are 1 for each gene symbol:
nonul <- d[!duplicated(d),c(1:52)] #51 non ul samples and symbol, 130X52
ul <- d[!duplicated(d),c(1,53:122)] #70 ul samples and gene, 130X71

names <- as.character(nonul[,1])

nonUL <- as.matrix(nonul)#130X52 matrix
nonUL <- nonUL[,2:52]#130X51 matrix
row.names(nonUL) <- names

UL <- as.matrix(ul)#130X71
UL <- UL[,2:71] #130X70
row.names(UL) <- names

df <- data.frame(ul)#keeps row names of symbols, no symbol field
df <- df[,2:71] #130X70 problems with data being factor not numeric when df->matrix

# each field has to be named, mutate() won't work, tried with vector names, all warning() errors
# This gives the true min and max, not the min of first or last... use word to grab
# by colnames(df), then remove the ',', and quotes and enter paragraphs to paste into each 
# min(), max(), mean(), this will add the true row values at the end as fields

UL_GeneStats <- df %>% rowwise() %>% mutate(Min = min(gsm1667147ul,gsm1667148ul,gsm1667149ul,gsm336202ul,gsm336203ul,gsm336204ul,
                                                      gsm336205ul,gsm336206ul,gsm336207ul,gsm336208ul,gsm336209ul,gsm336210ul,
                                                      gsm336211ul,gsm336212ul,gsm336213ul,gsm336214ul,gsm336215ul,gsm336216ul,
                                                      gsm336217ul,gsm336218ul,gsm336219ul,gsm336220ul,gsm336221ul,gsm336222ul,
                                                      gsm336223ul,gsm336224ul,gsm336225ul,gsm336226ul,gsm336227ul,gsm336228ul,
                                                      gsm336229ul,gsm336230ul,gsm336231ul,gsm336232ul,gsm336233ul,gsm336234ul,
                                                      gsm336235ul,gsm336236ul,gsm336237ul,gsm336238ul,gsm336239ul,gsm336240ul,
                                                      gsm336241ul,gsm336242ul,gsm336243ul,gsm336244ul,gsm336245ul,gsm336246ul,
                                                      gsm336247ul,gsm336248ul,gsm336249ul,gsm336250ul,gsm336251ul,gsm38689ul,
                                                      gsm38690ul,gsm38691ul,gsm38692ul,gsm38693ul,gsm38694ul,gsm38695ul,gsm9093ul,
                                                      gsm9094ul,gsm9095ul,gsm9096ul,gsm9097ul,gsm569429ul,gsm569430ul,gsm569431ul,
                                                      gsm569432ul,gsm569433ul), 
                                            
                                            Max = max(gsm1667147ul,gsm1667148ul,gsm1667149ul,gsm336202ul,gsm336203ul,gsm336204ul,
                                                      gsm336205ul,gsm336206ul,gsm336207ul,gsm336208ul,gsm336209ul,gsm336210ul,
                                                      gsm336211ul,gsm336212ul,gsm336213ul,gsm336214ul,gsm336215ul,gsm336216ul,
                                                      gsm336217ul,gsm336218ul,gsm336219ul,gsm336220ul,gsm336221ul,gsm336222ul,
                                                      gsm336223ul,gsm336224ul,gsm336225ul,gsm336226ul,gsm336227ul,gsm336228ul,
                                                      gsm336229ul,gsm336230ul,gsm336231ul,gsm336232ul,gsm336233ul,gsm336234ul,
                                                      gsm336235ul,gsm336236ul,gsm336237ul,gsm336238ul,gsm336239ul,gsm336240ul,
                                                      gsm336241ul,gsm336242ul,gsm336243ul,gsm336244ul,gsm336245ul,gsm336246ul,
                                                      gsm336247ul,gsm336248ul,gsm336249ul,gsm336250ul,gsm336251ul,gsm38689ul,
                                                      gsm38690ul,gsm38691ul,gsm38692ul,gsm38693ul,gsm38694ul,gsm38695ul,gsm9093ul,
                                                      gsm9094ul,gsm9095ul,gsm9096ul,gsm9097ul,gsm569429ul,gsm569430ul,gsm569431ul,
                                                      gsm569432ul,gsm569433ul))

UL_GeneStats <- data.frame(UL_GeneStats) #130X72
row.names(UL_GeneStats) <- names

Mean <- data.frame(rowMeans(UL_GeneStats)) #130X1, dplyr doesn't produce mean() 
#for all samples above, the dplyr() for mean on rowsums() gave values closer to max()
#use base mean() and cbind()

colnames(Mean) <- 'Mean'
UL_GeneStats <- cbind(UL_GeneStats,Mean) #130X73, added the Mean column last

dfn <- data.frame(nonul)#130X52, numeric values with symbol
dfn <- dfn[,2:52]#130X51

#use word after grabbing colnames(dfn), remove all ',' ,'"', and enter commands in view,
#then place in each min(), max(), mean()
nonUL_GeneStats <- dfn %>% rowwise() %>% mutate(Min=min(gsm1667144,gsm1667145,gsm1667146,gsm336252,gsm336253,
                                                        gsm336254,gsm336255,gsm336256,gsm336257,gsm336258,gsm336259,
                                                        gsm336260,gsm336261,gsm336262,gsm336263,gsm336264,gsm336265,
                                                        gsm336266,gsm336267,gsm336268,gsm336269,gsm336270,gsm336271,
                                                        gsm336272,gsm336273,gsm336274,gsm336275,gsm336276,gsm336277,
                                                        gsm336278,gsm52661,gsm52662,gsm52663, gsm52664,gsm52665,gsm52666,
                                                        gsm52667,gsm52668,gsm52669,gsm52670, gsm52671,gsm9098,gsm9099,gsm9100,
                                                        gsm9101,gsm9102,gsm569424,gsm569425,gsm569426,gsm569427,gsm569428),
                                                Max = max(gsm1667144,gsm1667145,gsm1667146,gsm336252,gsm336253,
                                                          gsm336254,gsm336255,gsm336256,gsm336257,gsm336258,gsm336259,
                                                          gsm336260,gsm336261,gsm336262,gsm336263,gsm336264,gsm336265,
                                                          gsm336266,gsm336267,gsm336268,gsm336269,gsm336270,gsm336271,
                                                          gsm336272,gsm336273,gsm336274,gsm336275,gsm336276,gsm336277,
                                                          gsm336278,gsm52661,gsm52662,gsm52663, gsm52664,gsm52665,gsm52666,
                                                          gsm52667,gsm52668,gsm52669,gsm52670, gsm52671,gsm9098,gsm9099,gsm9100,
                                                          gsm9101,gsm9102,gsm569424,gsm569425,gsm569426,gsm569427,gsm569428)) #130X53

nonUL_GeneStats <- data.frame(nonUL_GeneStats)#130X53
row.names(nonUL_GeneStats) <- names #adds the symbol names to nonUL stats table

Mean <- data.frame(rowMeans(nonUL_GeneStats)) #130X1
colnames(Mean) <- 'Mean'
nonUL_GeneStats <- cbind(nonUL_GeneStats,Mean) #130X54, added the Mean column last


nonul_T <- t(nonUL_GeneStats) #54X130
nonul_T <- data.frame(nonul_T) #makes data numeric data frame to plot
plot(nonUL_GeneStats$Mean)#plots scatterplot of first gene Mean
plot(nonul_T$DRD4)

UL_T <- t(UL_GeneStats) #73X130
UL_T <- data.frame(UL_T) #73X130
plot(UL_GeneStats$Mean)
plot(UL_T$DRD4)

#least and most expressed ul sample genes by mean
ul_Mean_lowest <- data.frame(UL_GeneStats[order(UL_GeneStats$Mean, decreasing=FALSE),])
leastExpressed_UL_mean <- ul_Mean_lowest[1:10,]
ul_Mean_highest <- data.frame(UL_GeneStats[order(UL_GeneStats$Mean, decreasing=TRUE),])
mostExpressed_UL_mean <- ul_Mean_highest[1:10,]

plot(mostExpressed_UL_mean$Mean)#decreasing as it should show
boxplot(t(mostExpressed_UL_mean))#shows by genes there are a few genes not decreasing overall
# indices 1,3,8, for RPLP2, DRD4, and GALR3

plot(leastExpressed_UL_mean$Mean)#increasing as it should show
boxplot(t(leastExpressed_UL_mean)) # some genes in boxplot, CBX2, HMGA2, and DMC1
# show as max of a subset, then a drop that climbs collectively per gene, 

#least and most expressed nonul sample genes by mean 
nonul_Mean_lowest <- data.frame(nonUL_GeneStats[order(nonUL_GeneStats$Mean, decreasing=FALSE),])
leastExpress_nonUL_mean <- nonul_Mean_lowest[1:10,]
nonul_Mean_highest <- data.frame(nonUL_GeneStats[order(nonUL_GeneStats$Mean, decreasing=TRUE),])
mostExpressed_nonUL_mean <- nonul_Mean_highest[1:10,]

plot(leastExpress_nonUL_mean$Mean)# increasing as it should show
boxplot(t(leastExpress_nonUL_mean))#increasing except drops at indices 2, 6, and 7
#those genes are CDHR5, GCGR, ZNF750

plot(mostExpressed_nonUL_mean$Mean) #decreasing as it should show
boxplot(t(mostExpressed_nonUL_mean))# RPLP2 and DRD4 drop before the next subsequent ordered
# gene overall

nonul_most <- data.frame(nonUL_GeneStats[order(nonUL_GeneStats$Max, decreasing=TRUE),])
#130X54, includes Mean field

nonul_least <- data.frame(nonUL_GeneStats[order(nonUL_GeneStats$Max, decreasing=FALSE),])#130X54
ul_most <- data.frame(UL_GeneStats[order(UL_GeneStats$Max, decreasing=TRUE),])#130X73
ul_least <- data.frame(UL_GeneStats[order(UL_GeneStats$Max, decreasing=FALSE),])#130X73

write.csv(ul_Mean_highest, 'ULgenesMostExpressedByMean.csv', row.names=TRUE) #by mean
write.csv(nonul_Mean_highest, 'mostExpressedNonULGeneByMean.csv', row.names=TRUE) #by mean
write.csv(nonul_most, 'nonul_most.csv', row.names=TRUE) #by Max
write.csv(nonul_least, 'nonul_least.csv', row.names=TRUE) #by Max
write.csv(ul_most, 'ul_most.csv', row.names=TRUE) #by max
write.csv(ul_least, 'ul_least.csv', row.names=TRUE)#by max

#############################################################################################
nonul_most <- read.csv('nonul_most.csv', sep=',', header=TRUE)#130X55, Min/Mean/Max,rm X field
ul_most <- read.csv('ul_most.csv', sep=',', header=TRUE) #130X74, Min/Mean/Max, rm x field
names <- ul_most[,1]
row.names(ul_most) <- names
ul_most <- ul_most[,2:74] #130X73
ul_most_t <- t(ul_most)#73X130
ul_most_t <- data.frame(ul_most_t) #73X130

names <- nonul_most[,1]
row.names(nonul_most) <- names
nonul_most <- nonul_most[,2:55] #130X54
nonul_most_t <- t(nonul_most)#54X130
nonul_most_t <- data.frame(nonul_most_t) #54X130

# ggplot of most expressed against least expressed with a fitted line
# by Max not mean values of each gene
g <- ggplot(ul_most_t, aes(x=ul_most_t$RPLP2, y=ul_most_t$ZNF750))
g= g+xlab('UL RPLP2 gene most expressed')
g= g+ xlim(-3,20)+ ylim(-3,20) 
g= g+ ylab('UL ZNF750 gene least expressed')
g= g+ geom_point(size=6, colour='black', alpha=0.2)
g= g+ geom_smooth(method='lm', colour='black')

g

boxplot(ul_most_t[,c(1,130)]) #boxplot of most expressed (RPLP2) and least (ZNF750) in UL
boxplot(nonul_most_t[,c(1,130)]) #boxplot of most expressed (RPLP2 and least (FSCN2)in nonUL
boxplot(nonul_most_t[,c(1,129)])
#when comparing the same two most and least expressed genes, the boxplots of UL and nonUL
# are similar, only the nonul shows the least expressed is now FSCN2. 

#adding fields that label each sample as which group they are derived from is helpful
#there are five data sets from GEO these 51 nonUL and 70 UL samples came from

#merge the two data sets of UL and nonUL together
nonUL <- data.frame(t(nonul_most_t)) #130X54, genes as obs with Min/Max/Mean
UL <- data.frame(t(ul_most_t)) #130X73, genes as obs, with Min/Max/Mean

colnames(nonUL)
# [1] "gsm1667144" "gsm1667145" "gsm1667146" "gsm336252"  "gsm336253"  "gsm336254"  "gsm336255" 
# [8] "gsm336256"  "gsm336257"  "gsm336258"  "gsm336259"  "gsm336260"  "gsm336261"  "gsm336262" 
# [15] "gsm336263"  "gsm336264"  "gsm336265"  "gsm336266"  "gsm336267"  "gsm336268"  "gsm336269" 
# [22] "gsm336270"  "gsm336271"  "gsm336272"  "gsm336273"  "gsm336274"  "gsm336275"  "gsm336276" 
# [29] "gsm336277"  "gsm336278"  "gsm52661"   "gsm52662"   "gsm52663"   "gsm52664"   "gsm52665"  
# [36] "gsm52666"   "gsm52667"   "gsm52668"   "gsm52669"   "gsm52670"   "gsm52671"   "gsm9098"   
# [43] "gsm9099"    "gsm9100"    "gsm9101"    "gsm9102"    "gsm569424"  "gsm569425"  "gsm569426" 
# [50] "gsm569427"  "gsm569428"  "Min"        "Max"        "Mean" 

colnames(UL)
# [1] "gsm1667147ul" "gsm1667148ul" "gsm1667149ul" "gsm336202ul"  "gsm336203ul"  "gsm336204ul" 
# [7] "gsm336205ul"  "gsm336206ul"  "gsm336207ul"  "gsm336208ul"  "gsm336209ul"  "gsm336210ul" 
# [13] "gsm336211ul"  "gsm336212ul"  "gsm336213ul"  "gsm336214ul"  "gsm336215ul"  "gsm336216ul" 
# [19] "gsm336217ul"  "gsm336218ul"  "gsm336219ul"  "gsm336220ul"  "gsm336221ul"  "gsm336222ul" 
# [25] "gsm336223ul"  "gsm336224ul"  "gsm336225ul"  "gsm336226ul"  "gsm336227ul"  "gsm336228ul" 
# [31] "gsm336229ul"  "gsm336230ul"  "gsm336231ul"  "gsm336232ul"  "gsm336233ul"  "gsm336234ul" 
# [37] "gsm336235ul"  "gsm336236ul"  "gsm336237ul"  "gsm336238ul"  "gsm336239ul"  "gsm336240ul" 
# [43] "gsm336241ul"  "gsm336242ul"  "gsm336243ul"  "gsm336244ul"  "gsm336245ul"  "gsm336246ul" 
# [49] "gsm336247ul"  "gsm336248ul"  "gsm336249ul"  "gsm336250ul"  "gsm336251ul"  "gsm38689ul"  
# [55] "gsm38690ul"   "gsm38691ul"   "gsm38692ul"   "gsm38693ul"   "gsm38694ul"   "gsm38695ul"  
# [61] "gsm9093ul"    "gsm9094ul"    "gsm9095ul"    "gsm9096ul"    "gsm9097ul"    "gsm569429ul" 
# [67] "gsm569430ul"  "gsm569431ul"  "gsm569432ul"  "gsm569433ul"  "Min"          "Max"         
# [73] "Mean" 

# remove the min/max/mean fields
nonUL <- nonUL[,-c(52:54)]
UL <- UL[,-c(71:73)]

Names <- data.frame(row.names(nonUL)) #130 X 1
Names2 <- data.frame(row.names(UL))#130 X 1
nonUL <- cbind(Names, nonUL) #130X52
UL <- cbind(Names2, UL) #130X71

colnames(nonUL)[1] <- 'genes'
colnames(UL)[1] <- 'genes'

all <- merge(nonUL, UL, by.x='genes', by.y='genes')#130X122, genes and samples
names <- all[,1]
all_t <- data.frame(t(all[,2:122]))
colnames(all_t) <- names # matrix of samplesXgenes in common bw UL/nonUL, 121X130

#save file out to use, this has all genes in common and all samples 121X130 matrix
write.csv(all_t, 'all-common-genes-samples.csv', row.names=TRUE)

# all_t data frame is 121 observations by 130 genes

#create a subset/grouping data frame to name which set the samples came from:
row.names(all_t)
# [1] "gsm1667144"   "gsm1667145"   "gsm1667146"   "gsm336252"    "gsm336253"    "gsm336254"   
# [7] "gsm336255"    "gsm336256"    "gsm336257"    "gsm336258"    "gsm336259"    "gsm336260"   
# [13] "gsm336261"    "gsm336262"    "gsm336263"    "gsm336264"    "gsm336265"    "gsm336266"   
# [19] "gsm336267"    "gsm336268"    "gsm336269"    "gsm336270"    "gsm336271"    "gsm336272"   
# [25] "gsm336273"    "gsm336274"    "gsm336275"    "gsm336276"    "gsm336277"    "gsm336278"   
# [31] "gsm52661"     "gsm52662"     "gsm52663"     "gsm52664"     "gsm52665"     "gsm52666"    
# [37] "gsm52667"     "gsm52668"     "gsm52669"     "gsm52670"     "gsm52671"     "gsm9098"     
# [43] "gsm9099"      "gsm9100"      "gsm9101"      "gsm9102"      "gsm569424"    "gsm569425"   
# [49] "gsm569426"    "gsm569427"    "gsm569428"    "gsm1667147ul" "gsm1667148ul" "gsm1667149ul"
# [55] "gsm336202ul"  "gsm336203ul"  "gsm336204ul"  "gsm336205ul"  "gsm336206ul"  "gsm336207ul" 
# [61] "gsm336208ul"  "gsm336209ul"  "gsm336210ul"  "gsm336211ul"  "gsm336212ul"  "gsm336213ul" 
# [67] "gsm336214ul"  "gsm336215ul"  "gsm336216ul"  "gsm336217ul"  "gsm336218ul"  "gsm336219ul" 
# [73] "gsm336220ul"  "gsm336221ul"  "gsm336222ul"  "gsm336223ul"  "gsm336224ul"  "gsm336225ul" 
# [79] "gsm336226ul"  "gsm336227ul"  "gsm336228ul"  "gsm336229ul"  "gsm336230ul"  "gsm336231ul" 
# [85] "gsm336232ul"  "gsm336233ul"  "gsm336234ul"  "gsm336235ul"  "gsm336236ul"  "gsm336237ul" 
# [91] "gsm336238ul"  "gsm336239ul"  "gsm336240ul"  "gsm336241ul"  "gsm336242ul"  "gsm336243ul" 
# [97] "gsm336244ul"  "gsm336245ul"  "gsm336246ul"  "gsm336247ul"  "gsm336248ul"  "gsm336249ul" 
# [103] "gsm336250ul"  "gsm336251ul"  "gsm38689ul"   "gsm38690ul"   "gsm38691ul"   "gsm38692ul"  
# [109] "gsm38693ul"   "gsm38694ul"   "gsm38695ul"   "gsm9093ul"    "gsm9094ul"    "gsm9095ul"   
# [115] "gsm9096ul"    "gsm9097ul"    "gsm569429ul"  "gsm569430ul"  "gsm569431ul"  "gsm569432ul" 
# [121] "gsm569433ul" 

samples <- data.frame(c(rep('GSE68295',3), rep('GSE13319', 27),rep('GSE2724',11),
                        rep('GSE593',5), rep('GSE23112',5), rep('GSE68295_UL',3),
                        rep('GSE13319_UL',50), rep('GSE2724_UL',7), rep('GSE593_UL',5),
                        rep('GSE23112_UL',5))) #121X1

# samples <- data.frame(c(rep('A',3), rep('B', 27),rep('C',11),
#              rep('D',5), rep('E',5), rep('A',3),
#              rep('B',50), rep('C',7), rep('D',5),
#              rep('E',5))) #121X1

#create a ggplot2 sample derived field to use color differentiation of sample source
colnames(samples) <- 'samples'
All <- cbind(samples,all_t) #matrix of samples X genes with meta field GSEsample derived

#write out this data table, samples added to matrix 121samples by 130 genes
write.csv(All, 'all_ggplot2_samples_derived.csv', row.names=TRUE)


#now categorical values can be added as a scatter showing color by sample group derived
g <- ggplot(All, aes(x=RPLP2, y=ZNF750))
g= g+xlab('UL RPLP2 gene most expressed')
g= g+ xlim(-3,20)+ ylim(-3,20) 
g= g+ ylab('UL ZNF750 gene least expressed')
g= g+ geom_point(aes(colour=samples), size=6, alpha=0.5)
g= g+ geom_smooth(method='lm', colour='black')

g

#zoom in on scales for same plot above with categorical samples colored

g <- ggplot(All, aes(x=RPLP2, y=ZNF750))
g= g+xlab('UL RPLP2 gene most expressed')
g= g+ xlim(min(All$RPLP2)-1,max(All$RPLP2)+1) + ylim(min(All$ZNF750)-1,max(All$ZNF750)+1)
g= g+ ylab('UL ZNF750 gene least expressed')
g= g+ geom_point(aes(colour=samples), size=6, alpha=0.5)
g= g+ geom_smooth(method='lm', colour='black')

g

# add a group that identifies as ul or not, to samplesXgenes & samples matrix
UL_nonUL <- data.frame(c(rep('nonUL',51), rep('UL',70)))#121X1
colnames(UL_nonUL) <- 'UL_nonUL'
All <- cbind(UL_nonUL,All)

#write this matrix out to csv, 121 samples X 132 (2 meta, 130 genes)
write.csv(All, 'All-ggplot2-type-sample-derived.csv', row.names=TRUE)

# plot most to least expressed genes RPLP2 to ZFN750 using UL or nonUL color categories
g <- ggplot(All, aes(x=RPLP2, y=ZNF750))
g= g+xlab('UL RPLP2 gene most expressed')
g= g+ xlim(min(All$RPLP2)-1,max(All$RPLP2)+1) + ylim(min(All$ZNF750)-1,max(All$ZNF750)+1)
g= g+ ylab('UL ZNF750 gene least expressed')
g= g+ geom_point(aes(colour=UL_nonUL), size=6, alpha=0.5)
g= g+ geom_smooth(method='lm', colour='black')

g

# plot TNRC6B, FASN, CYTH4, HMGA2, CCDC57, BET1L

# plot shows the FASN to HGMA2 expression in all samples
g <- ggplot(All, aes(x=FASN, y=HMGA2))
g= g+xlab('FASN')
g= g+ xlim(min(All$FASN)-1,max(All$FASN)+1) + ylim(min(All$HMGA2)-1,max(All$HMGA2)+1) 
g= g+ ylab('HMGA2')
g= g+ geom_point(aes(colour=samples), size=6, alpha=0.5)
g= g+ geom_smooth(method='lm', colour='black')

g

# plot shows UL and non UL expression levels of FASN to HMGA2
g <- ggplot(All, aes(x=FASN, y=HMGA2))
g= g+xlab('FASN')
g= g+ xlim(min(All$FASN)-1, max(All$FASN)+1)+ ylim(min(All$HMGA2)-1,max(All$HMGA2)+1) 
g= g+ ylab('HMGA2')
g= g+ geom_point(aes(colour=UL_nonUL), size=6, alpha=0.5)
g= g+ geom_smooth(method='lm', colour='black')

g

# plot shows UL_nonUL of TNRC6B to CYTH4
g <- ggplot(All, aes(x=TNRC6B, y=CYTH4))
g= g+xlab('TNRC6B')
g= g+ xlim(min(All$TNRC6B)-1,max(All$TNRC6B)+1) + ylim(min(All$CYTH4)-1,max(All$CYTH4)+1) 
g= g+ ylab('CYTH4')
g= g+ geom_point(aes(colour=UL_nonUL), size=6, alpha=0.5)
g= g+ geom_smooth(method='lm', colour='black')

g


#plot shows samples of TNRC6B to CYTH4

g <- ggplot(All, aes(x=TNRC6B, y=CYTH4))
g= g+xlab('TNRC6B')
g= g+ xlim(min(All$TNRC6B)-1,max(All$TNRC6B)+1) + ylim(min(All$CYTH4)-1,max(All$CYTH4)+1)  
g= g+ ylab('CYTH4')
g= g+ geom_point(aes(colour=samples), size=6, alpha=0.5)
g= g+ geom_smooth(method='lm', colour='black')

g

#plot shows UL_nonUL BET1L to TNRC6B

g <- ggplot(All, aes(x=BET1L, y=TNRC6B))
g= g+xlab('BET1L')
g= g+ xlim(min(All$BET1L)-1,max(All$BET1L)+1) + ylim(min(All$TNRC6B)-1,max(All$TNRC6B)+1)  
g= g+ ylab('TNRC6B')
g= g+ geom_point(aes(colour=UL_nonUL), size=6, alpha=0.5)
g= g+ geom_smooth(method='lm', colour='black')

g

#plot shows the samples of BET1L to TNRC6B
g <- ggplot(All, aes(x=BET1L, y=TNRC6B))
g= g+xlab('BET1L')
g= g+ xlim(min(All$BET1L)-1,max(All$BET1L)+1) + ylim(min(All$TNRC6B)-1,max(All$TNRC6B)+1)
g= g+ ylab('TNRC6B')
g= g+ geom_point(aes(colour=samples), size=6, alpha=0.5)
g= g+ geom_smooth(method='lm', colour='black')

g


########################################################################################
# The file to work with is now 'All-ggplot2-type-sample-derived.csv', 
# with 2 meta fields describing observation as the GSE sample derived, or if a UL or not
# this file is a matrix of genes as headers and samples as observations for all
# UL and nonUL samples from the five GSE samples used

all <- read.csv('All-ggplot2-type-sample-derived.csv', sep=',', header=TRUE,
                na.strings=c('','NA')) #121X133, an 'X' field of samples

names <- all[,1] #factor with 121 levels or gene samples
all <- all[,-1] #121X132
row.names(all) <- names

#####
#install.packages('lattice')
library(lattice)

splom(all[c(7:10,88:91)],
      main="Three nonUL and UL Genes")
#####

info <- all[,c(1,3:132)] #121X131, removes the GSE 'samples' field

library(dplyr)

# table of 2X2 for number of groups and count in each group
Means <- info %>% group_by(UL_nonUL) %>% summarise(n =n())
# # A tibble: 2 x 2
#    UL_nonUL     n
#     <fct>    <int>
#   1 nonUL       51
#   2 UL          70

names <- row.names(info)[1:51] #51 sample IDs for nonUL
names2 <- row.names(info)[52:121] #70 sample IDs for UL

UL <- filter(info, UL_nonUL=='UL') #70X131
nonUL <- filter(info, UL_nonUL=='nonUL') #51X131

UL <- data.frame(t(UL[,2:131])) #130X70, removes non-numeric field 'UL_nonUL'
colnames(UL) <- names2

nonUL <- data.frame(t(nonUL[,2:131])) #130X51, removes non-numeric field 'UL_nonUL'
colnames(nonUL) <- names

UL <- mutate(UL, UL_Mean = rowMeans(UL)) #130X71, a gene mean field added for UL
nonUL <- mutate(nonUL, nonUL_Mean = rowMeans(nonUL)) #130X52, a gene mene for nonUL added

MeanData <- cbind(nonUL,UL) #130X123
Means <- MeanData[,c(1:51,53:122,52,123)]

MeanData <- mutate(Means, Difference_UL_minus_non_means = UL_Mean - nonUL_Mean) #130X124
row.names(MeanData) <- row.names(info)

MeanData <- round(MeanData,2)#rounds the values, 130X124, 121 samples, 3 mean fields

write.csv(MeanData,'DE_data_unordered.csv', row.names=TRUE)

# The data is now ordered from least to most differential gene expression across samples
MeanData <- MeanData[order(MeanData$Difference_UL_minus_non_means),]#orders by DE least to grtst

DE_no_change <- filter(MeanData, Difference_UL_minus_non_means==0)#drops the gene name, not good


avg_change <- mean(MeanData$Difference_UL_minus_non_means) #0.018

# this data table shows genes that didn't change from UL-nonUL expression levels
DE_no_change <- MeanData[(MeanData$Difference_UL_minus_non_means==0), ]#GCGR gene originally,
# but now it is SIRT 1X124 after correcting row name assignment above

# this data table shows UL-nonUL expression levels giving decreased expression in UL tissue
Decr_UL_DE <- MeanData[(MeanData$Difference_UL_minus_non_means < 0.0),] #59X124

# this data table shows UL-nonUL expression level increases in UL tissue
Incr_UL_DE <- MeanData[(MeanData$Difference_UL_minus_non_means > 0.0),] #70X124

write.csv(DE_no_change, 'DE_noChange.csv', row.names=TRUE)
write.csv(Decr_UL_DE, 'DE_UL_decrease.csv', row.names=TRUE)
write.csv(Incr_UL_DE,'DE_UL_increase.csv', row.names=TRUE)
# Use lattice for splom() of pairwise layouts of genes

All <- read.csv('DE_data_unordered.csv', sep=',', header=TRUE)#130X125

names <- All[,1]
All <- All[,-1] #130X124, removes the row names as field 'X'

All <- data.frame(t(All)) 
colnames(All) <- names

DE <-All[124,] #difference UL-nonUL gene expression change in means, 1X130

non <- All[1:51,] #51X130
non_DE <- rbind(non,DE) #52X130

ul <- All[52:121,] #70X130
ul_DE <- rbind(ul,DE) #71X130

library(lattice)

# pairwise by genes in non-UL cases:
png('six_nonUL_lattice.png', width=600, height=600)
splom(non[c(10:15)],
      main="Six non-UL Genes")
dev.off()

#order non in decreasing differential expression
non_decr <- data.frame(t(non_DE)) #130X51, all row&col names included
non_most <- non_decr[order(non_decr$Difference_UL_minus_non_means,decreasing=TRUE),]#130X52

# the 10 genes with the most increase in gene expression from nonUL to UL sample means
non_10_most <- non_most[1:10,]
row.names(non_10_most)
# [1] "CBX2"   "CANT1"    "IRF7"     "ARHGDIA"  "NOL12"    "SLC25A10" "ASPSCR1"     "SYNGR1"  
# [9] "LLPH"     "BET1L"   

#when correcting the name assignment using the row.names(info) earlier and rerunning script:
# [1] "CBX2"    "RAC3"    "KDELR3"  "SOCS3"   "PYCR1"   "TH"      "ASPSCR1" "MICALL1" "EIF4A3" 
# [10] "NPTX1" 

# the 10 genes with the least increase or most decrease in gene expression from nonUL to UL
# sample means:
non_10_least <- non_most[121:130,] #10X52
non_10_least <- non_10_least[order(non_10_least$Difference_UL_minus_non_means),] #10X52
row.names(non_10_least)
# [1] "FSCN2"     "SOCS3"     "SLC38A10"  "TALDO1"    "FN3K"      "POLR2F"    "SMCR7L"   
# [8] "ASCL2"     "C17orf101" "HMGA2" 

#those results no longer are accurate:
# [1] "ZNF750"   "SOCS3"    "CBX7"     "C1QTNF1"  "SLC38A10" "CBX2"   "GRAP2"    "RAC2"    
# [9] "EPS8L2"   "TNNI2" 

non_10_least_t <- data.frame(t(non_10_least)) #52X10
non_10_most_t <- data.frame(t(non_10_most)) #52X10
plot(non_10_most$Difference_UL_minus_non_means)


##REPORT##
png('Pairwise_10_Most_DE_in_UL.png', width=800, height=800)
splom(non_10_most_t[1:10],
      main="Ten Genes Most DE in UL")
dev.off()

plot(non_10_least$Difference_UL_minus_non_means)
png('Pairwise_10_Least_DE_in_UL.png', width=800, height=800)
splom(non_10_least_t[1:10], main='Ten Genes Least DE in UL')
dev.off()

#order ul in decreasing differential expression
ul_decr <- data.frame(t(ul_DE)) #130X71, all row&col names included
ul_most <- ul_decr[order(ul_decr$Difference_UL_minus_non_means,decreasing=TRUE),]#130X71

ul_10_most <-ul_most[1:10,] #10X71
row.names(ul_10_most)
# [1] "CBX2"   "CANT1"    "IRF7"     "ARHGDIA"  "NOL12"    "SLC25A10" "ASPSCR1"     "SYNGR1"  
# [9] "LLPH"     "BET1L"

ul_10_most_t <- data.frame(t(ul_10_most))

ul_10_least <- ul_most[121:130,] #10X71
ul_10_least <- ul_10_least[order(ul_10_least$Difference_UL_minus_non_means),]
row.names(ul_10_least)
# [1] "FSCN2"     "SOCS3"     "SLC38A10"  "TALDO1"    "FN3K"      "POLR2F"    "SMCR7L"   
# [8] "ASCL2"     "C17orf101" "HMGA2

ul_10_least_t <- data.frame(t(ul_10_least))

png('Pairwise_10_Most_DE_in_UL.png', width=800, height=800)
splom(ul_10_most_t[1:10], main='Ten UL Genes Most DE in UL')
dev.off()

png('Pairwise_10_Least_DE_in_UL.png', width=800, height=800)
splom(ul_10_least_t[1:10], main='Ten UL Genes Least DE in UL')
dev.off()

#using row.names() for the genes, doesn't matter if ul, nonUL or all data frames
increased_DE <- ul_most[ul_most$Difference_UL_minus_non_means>0,] #70X71
increased_DE_ordered <- increased_DE[order(increased_DE$Difference_UL_minus_non_means,
                                           decreasing=TRUE),]
row.names(increased_DE_ordered) #ordered most to least increase in gene expression originally:
# [1] "CBX2"   "CANT1"    "IRF7"     "ARHGDIA"  "NOL12"    "SLC25A10" "ASPSCR1"     "SYNGR1"  
# [9] "LLPH"     "BET1L"    "POLR2L"   "LGALS1"   "HRAS"     "TNNI2"    "RPL3"     "DRD4"    
# [17] "CYTH1"    "CYTH4"    "ATF4"     "LGALS2"   "MRPL12"   "IGF2.AS"  "EIF4A3"   "BIRC5"   
# [25] "NPTX1"    "MFNG"     "ENGASE"   "C1QTNF1"  "RAC3"     "DMC1"     "GTPBP1"   "SGSM3"   
# [33] "CHMP6"    "PKP3"     "SIRT7"    "TOMM22"   "DUS1L"    "BAHCC1"   "GALR3"    "MAFG"    
# [41] "H1F0"     "ZNF750"   "CBX7"     "JOSD1"    "TBCD"     "DDX17"    "TRIOBP"   "SLC25A22"
# [49] "FN3KRP"   "TAB1"     "MICALL1"  "CD7"      "ASPSCR1"  "SOX10"    "APOBEC3G" "NPLOC4"  
# [57] "SCT"      "DEAF1"    "RFNG"     "DNAH17"   "TNRC6B"   "RAB40B"   "NPTXR"    "MRPL23"  
# [65] "DNAL4"    "SIRT3"    "RASSF7"   "GRAP2"    "AATK"     "KDELR3"

# BET1L, CYTH4, and TNRC6B are DE more in UL than nonUL tissue

# BUt this changed after info rownames change, to:
# [1] "CBX2"     "RAC3"     "KDELR3"   "SOCS3"    "PYCR1"    "TH"       "ASPSCR1"  "MICALL1" 
# [9] "EIF4A3"   "NPTX1"    "BET1L"    "ADSL"     "MGAT3"    "HRAS"     "APOBEC3C" "APOBEC3F"
# [17] "CEND1"    "POLR2L"   "ASCL2"    "MRPL12"   "FSCN2"    "RAB40B"   "DNAL4"    "TALDO1"  
# [25] "TOMM22"   "FN3KRP"   "NPTXR"    "PSMD13"   "TSSC4"    "SYNGR1"   "EIF3L"    "FASN"    
# [33] "IGF2.AS"  "LLPH"     "MAFG"     "CHMP6"    "LEMD3"    "NPLOC4"   "ATF4"     "CANT1"   
# [41] "CSNK1D"   "WDR45L"   "BIRC5"    "CCDC57"   "CTSD"     "CSNK1E"   "DUS1L"    "JOSD1"   
# [49] "CD81"     "LGALS1"   "NOL12"    "SGSM3"    "RFNG"     "TAB1"     "MRPL23"   "HMGA2"   
# [57] "SMCR7L"   "DEAF1"    "GNS"      "TNRC6B"   "ENGASE"   "GCGR"     "POLR2F"   "DDX17"   
# [65] "RPL3"     "IFITM3"   "PICK1"    "TBCD"     "SIRT3"    "PLA2G6" 
# 
# 
decreased_DE <- ul_most[ul_most$Difference_UL_minus_non_means<=0,]#60X71
decreased_DE_ordered <- decreased_DE[order(decreased_DE$Difference_UL_minus_non_means, 
                                           decreasing=FALSE),]

#the least expressed or most underexpressed of genes in UL tissue
row.names(decreased_DE_ordered)#original results:
# [1] "FSCN2"     "SOCS3"     "SLC38A10"  "TALDO1"    "FN3K"      "POLR2F"    "SMCR7L"   
# [8] "ASCL2"     "C17orf101" "HMGA2"     "PNPLA2"    "DCXR"      "PICK1"     "TSSC4"    
# [15] "IRAK3"     "MKL1"      "APOBEC3C"  "GNS"       "CCDC57"    "RPLP2"     "EPS8L2"   
# [22] "CEND1"     "CSNK1E"    "MGAT3"     "APOBEC3F"  "INS"       "ADSL"      "LEMD3"    
# [29] "SECTM1"    "EIF3L"     "CARD14"    "HGS"       "SLC16A8"   "C11orf21"  "GAA"      
# [36] "BAIAP2"    "TMEM80"    "CSNK1D"    "AZI1"      "PDE6G"     "CDC42EP1"  "TMC6"     
# [43] "PYCR1"     "RAC2"      "SOCS3"     "CBX2"      "IFITM3"    "CTSD"      "FASN"     
# [50] "CD81"      "PLA2G6"    "PSMD13"    "SIGIRR"    "FOXK2"     "ATHL1"     "CDHR5"    
# [57] "TMEM184B"  "WDR45L"    "TH"        "GCGR"  

# HMGA2, CCDC57, and FASN are under DE by means of UL compared to nonUL tissue

#results after correction to row names of MeanData from info error:
# [1] "ZNF750"    "SOCS3"     "CBX7"      "C1QTNF1"   "SLC38A10"  "CBX2"    "GRAP2"     "RAC2"     
# [9] "EPS8L2"    "TNNI2"     "IRAK3"     "CYTH1"     "SLC25A10"  "DNAH17"    "SCT"       "FN3K"     
# [17] "GALR3"     "INS"       "RPLP2"     "AATK"      "PKP3"      "PNPLA2"    "ATHL1"     "MFNG"     
# [25] "MKL1"      "SIGIRR"    "DRD4"      "CDHR5"     "DMC1"      "BAHCC1"    "TMC6"      "CARD14"   
# [33] "SECTM1"    "TMEM184B"  "C11orf21"  "PDE6G"     "CYTH4"     "C17orf101" "FOXK2"     "SLC25A22" 
# [41] "IRF7"      "ASPSCR1"      "GTPBP1"    "H1F0"      "CDC42EP1"  "TMEM80"    "APOBEC3G"  "HGS"      
# [49] "SOX10"     "ARHGDIA"   "LGALS2"    "DCXR"      "CD7"       "GAA"       "RASSF7"    "TRIOBP"   
# [57] "AZI1"      "BAIAP2"    "SLC16A8"   "SIRT7"

All_out <- read.csv('DE_data_unordered.csv', sep=',', header=TRUE)#130X125
colnames(All_out)[1] <- 'gene_symbol'
colnames(All_out)
# [1] "gene_symbol"                   "gsm1667144"                   
# [3] "gsm1667145"                    "gsm1667146"                   
# [5] "gsm336252"                     "gsm336253" ... 
# ... 
# [119] "gsm569430ul"                   "gsm569431ul"                  
# [121] "gsm569432ul"                   "gsm569433ul"                  
# [123] "nonUL_Mean"                    "UL_Mean"                      
# [125] "Difference_UL_minus_non_means" 

write.csv(All_out,'all_130_genes_symbol.csv', row.names=FALSE) #130X125
############################################################################################### 

############################################################################################### 

# Differential gene expression was calculated between ul and nonul samples
# the top ten under and over expressed genes in ul samples compared to nonul samples found
# Five data sets combined to give 130 genes in common along the same cytobands as top genes
# with 51 nonul and 70 ul samples, and
# 70 genes over expressed and 60 underexpressed in ul than nonul
# A chromosome map was built in middle section that should have gene names attached
# ggplots of scatters categorized as groups of ul/nonul or studies 1 through 5 of ul/nonul
# lattice plot was shown of few genes compared to others to find linear relationships
# boxplots by gene in each sample shown.
############################################################################################### 
# Next, attach the meta information with a merge of gene function and cytoband location
# to the latest data frame: 'all_130_genes_symbol.csv' #130 genes X 125 samples and mean/symbol
all <- read.csv('all_130_genes_symbol.csv', sep=',', 
                header=TRUE, na.strings=c('','NA'))  #130X125

All <- all[,c(1,123:125,2:122)] # puts the mean and symbol fields in front columns
All$gene_symbol <- gsub('IGF2.AS','IGF2-AS',All$gene_symbol) #change symbol name in META is IGF2-AS, else dropped

#the original meta fields were extracted from the GPL6480 txt file:
meta <- read.delim('GPL6480-9577.txt', sep='\t', quote="", header=TRUE,
                   comment.char='#', na.strings=c('','NA')) #41108X17
colnames(meta)
# [1] "ID"                   "SPOT_ID"              "CONTROL_TYPE"         "REFSEQ"              
# [5] "GB_ACC"               "GENE"                 "GENE_SYMBOL"          "GENE_NAME"           
# [9] "UNIGENE_ID"           "ENSEMBL_ID"           "TIGR_ID"              "ACCESSION_STRING"    
# [13] "CHROMOSOMAL_LOCATION" "CYTOBAND"             "DESCRIPTION"          "GO_ID"               
# [17] "SEQUENCE"      

# keep the GENE, GENE_SYMBOL, CYTOBAND, GENE_NAME, and DESCRIPTION FIELDS,
# many duplicated entries because of the chromosomal location field so exclude it
META <- meta[,c(6:8,14:15)]#41108X5

colnames(META)
# [1] "GENE"        "GENE_SYMBOL" "GENE_NAME"   "CYTOBAND"    "DESCRIPTION"

META <- META[complete.cases(META),]#30831X5, remove NAs in the meta fields selected

#combine the gene samples and meta fields by gene symbol
df <- merge(META, All, by.x='GENE_SYMBOL', by.y='gene_symbol') #226X129

#remove most of the duplicated entries
DF <- df[!duplicated(df$GENE_SYMBOL),]#130X129

library(dplyr)

data <- DF %>% group_by(GENE_SYMBOL) %>% summarise(n=n())#130X2, no duplicate genes in table
head(data,5)
# # A tibble: 5 x 2
# GENE_SYMBOL     n
# <fct>       <int>
#   1 AATK            1
# 2 ADSL            1
# 3 APOBEC3C        1
# 4 APOBEC3F        1
# 5 APOBEC3G        1

tail(data,5)
# # A tibble: 5 x 2
# GENE_SYMBOL     n
# <fct>       <int>
#   1 TOMM22          1
# 2 TRIOBP          1
# 3 TSSC4           1
# 4 WDR45L          1
# 5 ZNF750          1

# DF is the data of merged meta and gene samples and mean data information
colnames(DF)
# [1] "GENE_SYMBOL"                   "GENE"                         
# [3] "GENE_NAME"                     "CYTOBAND"                     
# [5] "DESCRIPTION"                   "nonUL_Mean"                   
# [7] "UL_Mean"                       "Difference_UL_minus_non_means"
# [9] "gsm1667144"                    "gsm1667145"                   
# [11] "gsm1667146"                    "gsm336252" 
# ... 
# [121] "gsm9094ul"                     "gsm9095ul"                    
# [123] "gsm9096ul"                     "gsm9097ul"                    
# [125] "gsm569429ul"                   "gsm569430ul"                  
# [127] "gsm569431ul"                   "gsm569432ul"                  
# [129] "gsm569433ul"

unique(DF$CYTOBAND) # chromosome bands are 11,12,17,22, when making the All data frame
# it was filtered for those chromosomes that the top six genes having SNPs associated with
# UL risk, 

# remove unused objects, other than the data set DF, meta, mean, and gene samples for genes
rm(all, All, data, df, meta, META)

# save DF file to use later in code when needed
write.csv(DF, 'Genes_Means_DE_Meta.csv', row.names=FALSE) #130 genes X 129 variables, 99.9KB
############################################################################################
# This could be a limitation of the study, as other chromosomes were not included          #
# to make the data smaller to work with, thus faster to extract, load, and transform (etl) #
############################################################################################

means <- read.csv('Genes_Means_DE_Meta.csv', sep=',',
                  header=TRUE, na.strings=c('','NA')) #130X129

Means <- means[order(means$Difference_UL_minus_non_means),] #130X129

Upregulated <- subset(Means, Means$Difference_UL_minus_non_means > 0) #70 genes upreg., 70X129
Downregulated <- subset(Means, Means$Difference_UL_minus_non_means <= 0) #60 genes downreg,x129

# save to csv up and down regulated genes
write.csv(Upregulated, 'upregulated.csv', row.names=FALSE) #55 kb, 70X129
write.csv(Downregulated, 'downregulated.csv', row.names=FALSE) #46.5 kb, 60X129
#################################################################################

# go back to the data of the ENSEMBL fields for Gviz chromosomal plotting

#################################################################################
# use this file ub_genes_ensembl_gviz.csv, this file excludes the exon field
Data <- read.csv('ub_genes_ensembl_gviz.csv', sep=',', header=TRUE, 
                 na.strings=c('','NA')) #149X129

library(dplyr)

df <- Data %>% group_by(symbol) %>% summarise(n=n()) #130X2 with some multiple counts

df_sort <- df[order(df$n, decreasing=TRUE),]
head(df_sort,20)
# # A tibble: 20 x 2
# symbol       n
# <fct>    <int>
#   1 AATK         3
# 2 SLC38A10     3
# 3 ATF4         2
# 4 BAIAP2       2
#... 
# 18 ADSL         1
# 19 APOBEC3C     1
# 20 APOBEC3F     1
rm(df_sort,df)

DF <- Data[!duplicated(Data$symbol),] #130X129, not same meta fields as previous section
DF$symbol

rm(Data)

# #if you don't have the Means object from 30 ealier lines of code, sectioned off ####:
means <- read.csv('Genes_Means_DE_Meta.csv', sep=',',
                  header=TRUE, na.strings=c('','NA')) #130X129

Means <- means[order(means$Difference_UL_minus_non_means),] #130X129

Means$GENE_SYMBOL

#compare genes in both data frames, both used genes in common with the 5 GSE series
#and filtered to only include the chromosomes of the top six genes ubiquitous to UL studies
#currently the latest research, chr:11,12,17, and 22
d <- data.frame(sort(DF$symbol)) #lists gene symbols in the data created for ensemble
m <- data.frame(sort(Means$GENE_SYMBOL)) # lists gene symbols in data created for ggplot2
md <- cbind(m,d)# the genes are the same 130 genes, combines both symbols, ordered incr. alpha.

head(md,5)
# sort.Means.GENE_SYMBOL. sort.DF.symbol.
# 1                    AATK            AATK
# 2                    ADSL            ADSL
# 3                APOBEC3C        APOBEC3C
# 4                APOBEC3F        APOBEC3F
# 5                APOBEC3G        APOBEC3G
tail(md)
# sort.Means.GENE_SYMBOL. sort.DF.symbol.
# 125                  TNRC6B          TNRC6B
# 126                  TOMM22          TOMM22
# 127                  TRIOBP          TRIOBP
# 128                   TSSC4           TSSC4
# 129                  WDR45L          WDR45L
# 130                  ZNF750          ZNF750

#rm(m,d,md,means,Data)#remove objects that won't be used again

# the symbols for genes indicates both data sets use only the same 130 genes,
# but the following shows the two data sets have different meta fields/variables
head(colnames(DF),12) #designed for Gviz with ENSEMBL fields retrieved/combined
# [1] "chromosome" "start"      "end"        "width"      "strand"     "gene"       "transcript"
# [8] "symbol"     "gsm1667144" "gsm1667145" "gsm1667146" "gsm336252" 

head(colnames(Means),12)#designed for plotting, DE of means, and informational purposes
# [1] "GENE_SYMBOL"                   "GENE"                         
# [3] "GENE_NAME"                     "CYTOBAND"                     
# [5] "DESCRIPTION"                   "nonUL_Mean"                   
# [7] "UL_Mean"                       "Difference_UL_minus_non_means"
# [9] "gsm1667144"                    "gsm1667145"                   
# [11] "gsm1667146"                    "gsm336252"  

Up <- subset(Means, Means$Difference_UL_minus_non_means>0)#70X129
Down <- subset(Means, Means$Difference_UL_minus_non_means <=0)#60X129

# Write data table DF that has been adjusted to have no duplicates and shown to 
# have the same genes as Means, but with ENSEMBL fields to use Gviz on
write.csv(DF, 'Data_130genes_Gviz.csv', row.names=FALSE)
write.csv(Up, 'upGviz.csv', row.names=FALSE)
write.csv(Down, 'downGviz.csv', row.names=FALSE)
#############################################################################################
#############################################################################################
#############################################################################################
# This section uses Gviz on the DF data frame and adds the symbol IDs to two plots
# read in the previous data frame, DF, of nonduplicated 130 genes X 121 samples and 
# ENSEMBL fields 
# use bioconductor install

# if (!requireNamespace("BiocManager", quietly = TRUE))

# install.packages("BiocManager")

# BiocManager::install("Gviz") #select 'a' to update all packages when asked

library(Gviz)
library(GenomicRanges) 

DF <- read.csv('Data_130genes_Gviz.csv', sep=',', header=TRUE, 
               na.strings=c('','NA')) #130X129

# read in the gviz up and down regulated data if this is a new session
Up <- read.csv('upGviz.csv', sep=',', header=TRUE, na.strings=c('','NA'))#70X129
Down <- read.csv('downGviz.csv', sep=',', header=TRUE, na.strings=c('','NA'))#60X129

#upregulated, genes expressed more in UL samples

data <- DF[,1:8] #130X8, only keep meta of gviz data frame to combine upregulated genes with
UP <- merge(data,Up, by.x='symbol', by.y='GENE_SYMBOL') #70X136
write.csv(UP, 'UPGviz.csv', row.names=FALSE)#59.9 kb to use later

#downregulated, genes expressed less in UL

Down <- merge(data, Down, by.x='symbol', by.y='GENE_SYMBOL') #60X136
write.csv(Down, 'DownGviz.csv', row.names=FALSE)#50.7 kb, to use in later sections

# get a count of the number of genes on each chromosome, the 5 GSE series have in common

length(grep('chr22', DF$chromosome)) #43 genes
length(grep('chr11', DF$chromosome)) #33 genes
length(grep('chr17', DF$chromosome)) #48 genes
length(grep('chr12', DF$chromosome)) #6 genes

# get a count of the number of genes on each chromosome, from the increased gene expression
# in UL compared to non-UL samples
length(grep('chr22', UP$chromosome)) # 25 genes
length(grep('chr11', UP$chromosome)) # 16 genes
length(grep('chr17', UP$chromosome)) # 24 genes
length(grep('chr12', UP$chromosome)) # 5 genes

# get a count of the number of genes on each chromosome, from the decreased gene expression
# in UL compared to non-UL samples
length(grep('chr22', Down$chromosome)) # 18 genes
length(grep('chr11', Down$chromosome)) # 17 genes
length(grep('chr17', Down$chromosome)) # 24 genes
length(grep('chr12', Down$chromosome)) # 1 genes

# you could also use dplyr
library(dplyr)

counts <- DF %>% group_by(chromosome) %>% summarise(n=n())
counts
# # A tibble: 4 x 2
# chromosome     n
# <fct>      <int>
# 1 chr11         33
# 2 chr12          6
# 3 chr17         48
# 4 chr22         43
colnames(counts) <- c('chromosome','all')

counts_up <- UP %>% group_by(chromosome) %>% summarise(n=n())
counts_up
# # A tibble: 4 x 2
# chromosome     n
# <fct>      <int>
# 1 chr11         16
# 2 chr12          5
# 3 chr17         24
# 4 chr22         25

colnames(counts_up) <- c('chromosome','up')

counts_up <- counts_up[,2] # keep only counts of up regulated

counts_down <- Down %>% group_by(chromosome) %>% summarise(n=n())
counts_down
# # A tibble: 4 x 2
# chromosome     n
# <fct>      <int>
# 1 chr11         17
# 2 chr12          1
# 3 chr17         24
# 4 chr22         18
colnames(counts_down) <- c('chromosome','down')

counts_down <- counts_down[,2] # keep only counts of down regulated

All_counts <- cbind(counts, counts_up, counts_down)
All_counts
# chromosome   all up down
# 1      chr11  33 16   17
# 2      chr12   6  5    1
# 3      chr17  48 24   24
# 4      chr22  43 25   18

# number of genes up or down regulated among the samples is almost split down the middle
# for chromosomes 11 and 17, but chromosome 12 has more down regulated, while 
# chromosome 22 has more up regulated

# write this table to csv
write.csv(All_counts, 'All_chr_counts.csv', row.names=FALSE)

up_genes <- data.frame(sort(unique(UP$symbol)), type='up') #70X2
colnames(up_genes) <- c('genes','type')

down_genes <- data.frame(sort(unique(Down$symbol)), type='down') #60X2
colnames(down_genes) <- c('genes','type')

gene_class <- rbind(up_genes,down_genes) #130X2, rowbinds the up and down regulated genes

genes_chr <- merge(gene_class, DF[,c(1,8)], by.x='genes',by.y='symbol') #130X3, chromosome included

genes_chr_counts <- merge(genes_chr, All_counts, by.x='chromosome',
                          by.y='chromosome', all.x=TRUE) #130X6, 
# The above table is of counts for all genes belonging to each of the four chromosomes,
# next to the chromosomes and genes in the data of all, increased expression is the number
# of up genes on that chromosome, and the decreased expression is the number of down genes
# on that chromosome,decreased expression is down, all is the total counts of genes 
# on that chromosome, This is interesting to view where the top six genes ubiquitous to UL
# studies are grouped in this data, to point to a relationship of each gene being over or
# under expressed in UL tissue when compared to normal myometrium tissue

# write the above table out to csv to use later in plotting with ggplot2 and analysis
write.csv(genes_chr_counts, 'genes_up_down_per_chr.csv', row.names=FALSE)

colnames(genes_chr_counts)
#[1] "chromosome" "genes"      "type"       "all"        "up"         "down" 
# the chromosome is where the gene is located
# the gene is the gene in common between the 5 GSE series
# the type is the group the gene was located in for either over or under expressed in UL
# the all is how many genes there are for that chromosome
# the up is how many up or over expressed genes for that chromosome
# the down is how many down or under expressed genes for that chromosome

# from the above table, list all the genes that are not among the group majority of genes
# that are up or down regulated for that chromosome

# using dplyr, group by chromome and filter
Maj_down <- filter(genes_chr_counts, down > up)#33X6
Maj_up <- filter(genes_chr_counts, up>down) #49X6
Equal <- filter(genes_chr_counts, up==down)#48X6

#adds a new Boolean field if gene is part of the majority for that chromosome, dplyr, mutate()
Maj_down <- mutate(Maj_down, majority=(type=='down'))#33X7
Maj_up <- mutate(Maj_up, majority=(type=='up'))#49X7
Equal <- mutate(Equal, majority='Equal')#48X7

#combine by obs, the two tables to show if each gene is part of the majority for that chromosome
member <- rbind(Maj_up, Maj_down, Equal) #130X7

#write this file out to csv
write.csv(member, 'member.csv', row.names=FALSE) #4.8kb


# get the majority genes that are part of the majority of either up or down regulated by chromosome
majority <- filter(member, majority=='TRUE') #47X7

# get the genes not part of the majority of up/down genes for that chromosome
minority <- filter(member, majority=='FALSE') #35X7

equal <- filter(member, majority=='Equal')#48X7

#write these last two table out to csv
write.csv(majority, 'majority.csv', row.names=FALSE)#1.8 kb
write.csv(minority, 'minority.csv', row.names=FALSE)#1.4 kb
write.csv(equal, 'equal.csv', row.names=FALSE)#1.9kb
############################
# create a minorityMerged and majorityMerged file to attach samples to from the down and
# up regulated or expressed UL tables

#clear your environment of clutter and objects then 
# Read in the files if you don't have them:
minority <- read.csv('minority.csv', sep=',', header=TRUE, na.strings=c('','NA'))#35X7
majority <- read.csv('majority.csv', sep=',', header=TRUE, na.strings=c('','NA'))#47X7
equal <- read.csv('equal.csv', sep=',', header=TRUE, na.strings=c('','NA')) #48X7
Up <- read.csv('UPGviz.csv', sep=',', header=TRUE, na.strings=c('','NA'))#70X136
Down <- read.csv('DownGviz.csv', sep=',', header=TRUE, na.strings=c('','NA'))#60X136

colnames(Up)
# [1] "symbol"                        "chromosome"                   
# [3] "start"                         "end"                          
# [5] "width"                         "strand"                       
# [7] "gene"                          "transcript"                   
# [9] "GENE"                          "GENE_NAME"                    
# [11] "CYTOBAND"                      "DESCRIPTION"                  
# [13] "nonUL_Mean"                    "UL_Mean"                      
# [15] "Difference_UL_minus_non_means" "gsm1667144"                   
# [17] "gsm1667145"                    "gsm1667146"                   
# [19] "gsm336252"                     "gsm336253"  
# ... 
# [131] "gsm9097ul"                     "gsm569429ul"                  
# [133] "gsm569430ul"                   "gsm569431ul"                  
# [135] "gsm569432ul"                   "gsm569433ul"  
# 

#remove the Up and Down 'chromosome' field so there aren't duplicates in the merge
Up <- Up[,-2]# 70X135
Down <- Down[,-2]#60X135


colnames(minority)
# [1] "chromosome" "genes"      "type"       "all"        "up"         "down"       "majority" 

# Make the minority of the up data table
minorityUp <- merge(minority, Up, by.x='genes', by.y='symbol')#16X141, 
#when the majority of genes on that chromosome are underexpressed

minorityDown <- merge(minority, Down, by.x='genes', by.y='symbol') #19X141
#chromosome are over expressed in UL samples

majorityUp <- merge(majority, Up, by.x='genes', by.y='symbol') #30X141
# are the majority of genes on that chromosome

majorityDown  <- merge(majority, Down, by.x='genes', by.y='symbol')#17X141
#of genes on that chromosome underexpressed

equalDown <- merge(equal, Down, by.x='genes', by.y='symbol')#24X141
16+19+30+17+24 #106

equalUp <- merge(equal, Up, by.x='genes', by.y='symbol')#24X141
106+24 #130

#merge all together

allMembers <- rbind(minorityUp,minorityDown,majorityUp,majorityDown,equalUp,equalDown)#130X141

#write to csv
write.csv(minorityUp, 'minorityUp.csv', row.names=FALSE)
write.csv(minorityDown, 'minorityDown.csv', row.names=FALSE)
write.csv(majorityUp, 'majorityUp.csv', row.names=FALSE)
write.csv(majorityDown, 'majorityDown.csv', row.names=FALSE)
write.csv(equalDown, 'equalDown.csv', row.names=FALSE)
write.csv(equalUp, 'equalUp.csv', row.names=FALSE)
write.csv(allMembers, 'allMembers.csv', row.names=FALSE)


minorityUp <- read.csv('minorityUp.csv', sep=',', header=TRUE)
minorityDown <- read.csv('minorityDown.csv', sep=',', header=TRUE)
majorityUp <- read.csv('majorityUp.csv', sep=',', header=TRUE)
majorityDown <- read.csv('majorityDown.csv', sep=',', header=TRUE)
equalUp <- read.csv('equalUp.csv', sep=',', header=TRUE)
equalDown <- read.csv('equalDown.csv', sep=',', header=TRUE)
allMembers <- read.csv('allMembers.csv', sep=',', header=TRUE)

# now plot the reverse and forward strands of cytoband locations of the UL risk genes using Gviz
library(Gviz)
library(GenomicRanges) 


# chromosome 11 reverse strand of all genes
rev11 <- subset(allMembers, allMembers$strand=='-')#66X141
chr <- as.character(unique(rev11$chromosome)) #[1] "chr11"  "chr22" "chr12" "chr17"
grtrack <- GeneRegionTrack(rev11, genome = "hg38", chromosome = chr[1], 
                           strand=rev11$strand,
                           gene=rev11$genes, symbol=rev11$genes, 
                           transcriptAnnotation = "genes", #adds the genes to bottom track
                           background.title = "brown", lwd=1, cex=1,
                           stackHeight=0.5,
                           name = "Reverse Strand Cytoband 11p15.5 of all Genes")

itrack <- IdeogramTrack(genome = "hg38", chromosome = chr[1])
gtrack <- GenomeAxisTrack()

png('allRevChr11.png',width=500,height=400)
plotTracks(list(itrack, gtrack, grtrack),#omitted atrack--messy
           featureAnnotation = "id",  lwd=1, cex=2,
           shape='arrow',
           background.panel = "beige", reverseStrand=FALSE, from = 1e+05, to = 3e+05,
           # add35=TRUE, add53=TRUE, #adds 3' to 5' and 5' to 3' on ideogram header
           showBandId = TRUE, cex.bands = 2,
           scale=.5,showOverplotting=TRUE,
           background.title = "#005544")# title color green-teal
dev.off()

# Forward strand genes on chromosome 12
fwd12 <- subset(allMembers, allMembers$strand=='+')
chr <- as.character(unique(fwd12$chromosome)) #[1] "chr11" "chr22"  "chr12" "chr17"
grtrack <- GeneRegionTrack(fwd12, genome = "hg38", chromosome = chr[3], 
                           strand=fwd12$strand,
                           gene=fwd12$genes, symbol=fwd12$genes, 
                           transcriptAnnotation = "genes", #adds the genes to bottom track
                           background.title = "brown", lwd=1, cex=1,
                           stackHeight=0.5,
                           name = "Forward Strand Cytoband 12q14.3 All Genes")

itrack <- IdeogramTrack(genome = "hg38", chromosome = chr[3])
gtrack <- GenomeAxisTrack()


png('allFwdChr12.png',width=500,height=400)
plotTracks(list(itrack, gtrack, grtrack),
           featureAnnotation = "id", shape='arrow', lwd=1, cex=1,
           background.panel = "beige",
           add35=TRUE, add53=TRUE, #adds 3' to 5' and 5' to 3' on ideogram header
           showBandId = TRUE, cex.bands = 0.5, from = 6.5e+07, to = 6.7e+07,
           scale=1,showOverplotting=TRUE,
           background.title = "#005544")# title color green-teal
dev.off()


# all genes on reverse strand chromosome 17
rev <- subset(allMembers, allMembers$strand=='-')
chr <- as.character(unique(rev$chromosome)) #[1] "chr11" "chr22" "chr12" "chr17"
grtrack <- GeneRegionTrack(rev, genome = "hg38", chromosome = chr[4], 
                           strand=rev$strand,
                           gene=rev$genes, symbol=rev$genes, 
                           transcriptAnnotation = "genes", #adds the genes to bottom track
                           background.title = "brown", lwd=1, cex=2,
                           stackHeight=0.5,
                           name = "Reverse Strand 17.25q.3 All Genes")

itrack <- IdeogramTrack(genome = "hg38", chromosome = chr[4])
gtrack <- GenomeAxisTrack()


png('allGenesRevChr17.png',width=500,height=400)
plotTracks(list(itrack, gtrack, grtrack),
           featureAnnotation = "id", shape='arrow', lwd=1, cex=1,
           background.panel = "beige", from = 8.15e+07, to = 8.3e+07,
           add35=TRUE, add53=TRUE, #adds 3' to 5' and 5' to 3' on ideogram header
           showBandId = TRUE, cex.bands = 0.5, reverseStrand=FALSE,
           scale=1,showOverplotting=TRUE,
           background.title = "#005544")# title color green-teal
dev.off()

# all Forward Strand genes on chromosome 22
fwd <- subset(allMembers, allMembers$strand=='+')
chr <- as.character(unique(fwd$chromosome)) #[1] "chr11" "chr22" "chr12" "chr17"
grtrack <- GeneRegionTrack(fwd, genome = "hg38", chromosome = chr[2], 
                           strand=fwd$strand,
                           gene=fwd$genes, symbol=fwd$genes, 
                           transcriptAnnotation = "genes", #adds the genes to bottom track
                           background.title = "brown", lwd=1, cex=2,
                           stackHeight=0.5,
                           name = "Forward Strand 22q13.1 All Genes")

itrack <- IdeogramTrack(genome = "hg38", chromosome = chr[2])
gtrack <- GenomeAxisTrack()


png('fwdAllChr22.png',width=500,height=400)
plotTracks(list(itrack, gtrack, grtrack),
           featureAnnotation = "id", shape='arrow', lwd=1, cex=1,
           background.panel = "beige", from = 3.7e+07, to = 4.1e+07,
           add35=TRUE, add53=TRUE, #adds 3' to 5' and 5' to 3' on ideogram header
           showBandId = TRUE, cex.bands = 0.5, reverseStrand=FALSE,
           scale=1,showOverplotting=TRUE,
           background.title = "#005544")# title color green-teal
dev.off()


# save to csv the allMembers data set of all meta with all genes 
# in same cytobands of 6 UL risk genes
write.csv(allMembers, 'MemberGviz_130_141.csv', row.names=FALSE) #111.5 kb

MemberGviz <- read.csv('MemberGviz_130_141.csv', sep=',', header=TRUE)#130X141

B <- grep('BET1L', MemberGviz$genes) #2
T <- grep('TNRC6B', MemberGviz$genes)#64
CY <- grep('CYTH4', MemberGviz$genes) #21
CC <- grep('CCDC57', MemberGviz$genes) #87
F <- grep('FASN', MemberGviz$genes) #93
H <- grep('HMGA2', MemberGviz$genes) #46

#create a data table of only these genes ubiquitous to the UL studies currently published
ubiq <- MemberGviz[c(B,T,CY,CC,F,H),] #6x141

ubiq[,c(1:7)]
#     genes chromosome type all up down majority
# 2   BET1L      chr11   up  33 16   17    FALSE
# 64 TNRC6B      chr22   up  43 25   18     TRUE
# 21  CYTH4      chr22 down  43 25   18    FALSE
# 87 CCDC57      chr17   up  48 24   24    Equal
# 93   FASN      chr17   up  48 24   24    Equal
# 46  HMGA2      chr12   up   6  5    1     TRUE

library(dplyr)

MemberMagnitude <- mutate(MemberGviz, 
                          Magnitude= abs(MemberGviz$Difference_UL_minus_non_means))
#130X142, has a magnitude field to see up/down change in UL

# place Magnitude field by DE field
Magnitude <- MemberMagnitude[,c(1:20,142,21:141)] #130X142, 

# order the table to list highest DE change
Magnitude_Ordered <- Magnitude[order(Magnitude$Magnitude, decreasing=TRUE),] #130X142

# table of top genes ordered by DE change in UL from nonUL, Gviz and Mean DE field and Samples
write.csv(Magnitude_Ordered, 'MemberMagnitude_130_142.csv', row.names=FALSE)

Magnitude_Ordered <- read.csv('MemberMagnitude_130_142.csv', sep=',', 
                              header=TRUE, na.strings=c('','NA')) #130X142

B <- grep('BET1L', Magnitude_Ordered$genes) #2
T <- grep('TNRC6B', Magnitude_Ordered$genes)#110
CY <- grep('CYTH4', Magnitude_Ordered$genes) #73
CC <- grep('CCDC57', Magnitude_Ordered$genes) #85
F <- grep('FASN', Magnitude_Ordered$genes) #62
H <- grep('HMGA2', Magnitude_Ordered$genes) #101

ubiqOrdered <- Magnitude_Ordered[c(B,H,CY,CC,F,T),] #6X142, order of highest up/down change
# for the top six genes ubiquitous to current UL research studies

ubiqOrdered$genes
# [1] BET1L  HMGA2  CYTH4  CCDC57 FASN   TNRC6B

Top10 <- Magnitude_Ordered[1:10,] #10X142
Top10$genes
### Not anymore: [1] FSCN2    CBX2   SOCS3    CANT1    IRF7     ARHGDIA  NOL12    SLC25A10 SLC38A10 ASPSCR1
##new values:
#[1] ZNF750  CBX2    SOCS3   RAC3    KDELR3  SOCS3   PYCR1   TH      CBX7    ASPSCR1


# combine the ubiquitous genes and the top 10 genes with 
# the most up/down gene expression in UL

B <- grep('BET1L', Magnitude_Ordered$genes) #22
T <- grep('TNRC6B', Magnitude_Ordered$genes)#110
CY <- grep('CYTH4', Magnitude_Ordered$genes) #73
CC <- grep('CCDC57', Magnitude_Ordered$genes) #85
F <- grep('FASN', Magnitude_Ordered$genes) #62
H <- grep('HMGA2', Magnitude_Ordered$genes) #101

ubiq_top_10 <- rbind(Top10,ubiqOrdered)#16X142

write.csv(ubiq_top_10, 'top_10_and_ubiq_6.csv', row.names=FALSE)#16X142

# It would be interesting to take these 16 genes and run a 10,000 sample bootstrap
# simulation with replacement on the 121 samples per gene to determine a simulated
# population mean for each gene, instead of the sample mean, then to test those same 
# 121 samples as partitions, so that 70% of the 121 samples are trained, while the 30%
# remaining are used to test a fitted model to see how well that model will be able
# to predict whether the sample of those 16 genes is UL or not UL, 
# 9 genes are part of the majority of UL DE change and the other 7 are the minority

# to get started with this remove the fields other than the gene name and the sample IDs

ubiq_top10_samples_only <- ubiq_top_10[,-c(2:21)] #16X122

write.csv(ubiq_top10_samples_only, 'ubiq_and_top10_samples_only.csv', row.names=FALSE)#10.9kb

#clear out the environment, and work with the last table saved to csv

samples_16 <- read.csv('ubiq_and_top10_samples_only.csv', sep=',', header=TRUE,
                       na.strings=c('','NA')) #16X122
samples_16$genes
# ZNF750  CBX2    SOCS3   RAC3    KDELR3  SOCS3   PYCR1   TH      CBX7    ASPSCR1 BET1L   HMGA2  
# CYTH4   CCDC57  FASN    TNRC6B 
top16_names <- data.frame(samples_16$genes)
colnames(top16_names) <- 'genes'
write.csv(top16_names)

# separate this data table into UL and nonUL samples
colnames(samples_16)#genes, then 51 nonUL samples(1:52), and 70 UL samples(1,53:122)

UL <- samples_16[,c(1,53:122)] #16X71, includes gene names
nonUL <- samples_16[,c(1:52)] #16X52, includes gene names

# create row names of the genes field and remove genes field
row.names(UL) <- UL[,1]
UL <- UL[,-1] #16X70, removes gene name field, only ul samples, genes are row names

row.names(nonUL) <- nonUL[,1]
nonUL <- nonUL[,-1] #16X51, 51 nonUL samples and gene names are the row names

# Transpose the data, so samples are observations for each gene as a field

UL_t <- data.frame(t(UL)) #genes are the header, and samples are row names, 70X16 matrix
nonUL_t <- data.frame(t(nonUL)) # genes are the header and samples are row names, 51X16 matrix
# retrieve the bootstrap code used in a previous script SaveAsTablesInRStudio.R on other
# data in a different folder to get stats on each of these 16 genes

# write these out to use in bootstrapping simulations for each of 16 genes

write.csv(UL_t, 'UL_t_16.csv', row.names=TRUE)
write.csv(nonUL_t, 'nonUL_t_16.csv', row.names=TRUE)

#############################################################################################
#############################################################################################
#############################################################################################
#BOOTSTRAPPING Top 10 + 6 ubiquitous genes
# keeps the row names as row names instead of adding a new 'X' field
UL_t <- read.csv('UL_t_16.csv', sep=',', header=TRUE, row.names=1)
nonUL_t <- read.csv('nonUL_t_16.csv', sep=',', header=TRUE, row.names=1)

# install.packages("UsingR") #for version 3.6 R
library(UsingR)

set.seed(45623489)# set this seed so the sampling is the same output
# Example:

ZNF750 <- as.numeric(nonUL_t[,1]) # this gene, FSCN2, in the 51 nonUL samples
n <- length(ZNF750)
B <- 10000

resamples <- matrix(sample(ZNF750,n*B,replace=TRUE),B,n)#large matrix (51000 elements, 3.9MB)
#10000 simulations X 51 samples, sampled and replaced in each matrix row and column

resampledRowMeans <- apply(resamples, 1, mean) #num 1:10000 values, applies over rows,
# each row simulates an independent, identical, replaced, sample 
#from the 51 samples of FSCN2

# get the column means for each sample spanning 10,000 simulated draws in sampling
resampledColMeans <- apply(resamples,2,mean)# sample means of 10k simulated nonUL FSCN2 gene
# numeric length 51

sd(resampledRowMeans)#[1] 0.23882 for the standard deviation of mean by row
quantile(resampledRowMeans, c(0.025,0.975)) #confidence interval for 5% two tail distribution
#   2.5%     97.5% 
# 0.8768431 1.8205882 

sd(resampledColMeans) # [1] 0.01619119, 
quantile(resampledColMeans, c(0.025,0.975))
#   2.5%    97.5% 
# 1.325118 1.378639

# plot histograms of the means by row, resampledRowMeans 
# and means by columns, resampledColMeans

library(lattice)
######################################################################### IGNORE all comments/results

png('ZNF750_hist_row_means.png', width=500, height=600)
histogram(resampledRowMeans, type='count', xlab='10,000 Simulated ZNF750 Row Means',
          ylab='Counts of ZNF750',  )
dev.off()

png('ZNF750_hist_col_means.png', width=800, height=1200)
histogram(resampledColMeans, type='count', xlab='10,000 Simulated ZNF750 Col Means',
          ylab='Counts of ZNF750',  )
dev.off()

# the row means of resampledRowMeans of nonUL_t is a better fit for a gaussian bell curve 
mean(resampledRowMeans)# 1.353323
sd(resampledRowMeans)#0.2419191
quantile(resampledRowMeans, c(0.025, 0.9725))
#     2.5%    97.25% 
#   0.8768431 1.8096348 

ZNF750_mean <- data.frame(mean(resampledRowMeans), row.names = 'ZNF750')
ZNF750_sd <- data.frame(sd(resampledRowMeans), row.names='ZNF750')
ZNF750_qu <- t(data.frame(quantile(resampledRowMeans, c(0.025, 0.9725))))
row.names(ZNF750_qu) <- 'ZNF750'
ZNF750_stats <- cbind(ZNF750_mean, ZNF750_sd, ZNF750_qu) #1X4

# find the stats for the UL samples of ZNF750
ZNF750_ul <- as.numeric(UL_t[,1]) # this gene, ZNF750, in the 51 nonUL samples
n <- length(ZNF750_ul)
B <- 10000

resamples <- matrix(sample(ZNF750_ul,n*B,replace=TRUE),B,n)#large matrix (70000 elements, 
# 5.3MB)
# 10000 simulations X 70 samples, sampled and replaced in each matrix row and column

resampledRowMeans <- apply(resamples, 1, mean) #num 1:10000 values, applies over rows,
# each row simulates an independent, identical, replaced, sample 
#from the 70 samples of ZNF750 ul

png('ZNF750_ul_hist_row_means.png', width=500, height=600)
histogram(resampledRowMeans, type='count', xlab='10,000 Simulated ZNF750 in UL Row Means',
          ylab='Counts of ZNF750 in UL',  )
dev.off()


sd(resampledRowMeans)#[1] 0.2157723 for the standard deviation of mean by row
quantile(resampledRowMeans, c(0.025,0.9725)) #confidence interval for 5% two tail distribution
#   2.5%     97.5% 
# 0.1152786 0.9594429 

ZNF750_ul_mean <- data.frame(mean(resampledRowMeans), row.names = 'ZNF750_ul')
ZNF750_ul_sd <- data.frame(sd(resampledRowMeans), row.names='ZNF750_ul')
ZNF750_ul_qu <- t(data.frame(quantile(resampledRowMeans, c(0.025,0.9725))))
row.names(ZNF750_ul_qu) <- 'ZNF750_ul'
ZNF750_ul_stats <- cbind(ZNF750_ul_mean, ZNF750_ul_sd, ZNF750_ul_qu) #1X2

ZNF750 <- rbind(ZNF750_stats, ZNF750_ul_stats) #2X4



# now for the next gene in the set of 16 genes, CBX2
CBX2 <- as.numeric(nonUL_t[,2]) 
n <- length(CBX2)
B <- 10000

resamples <- matrix(sample(CBX2,n*B,replace=TRUE),B,n)

resampledRowMeans <- apply(resamples, 1, mean) 

png('CBX2_hist_row_means.png', width=500, height=600)
histogram(resampledRowMeans, type='count', xlab='10,000 Simulated CBX2 Row Means',
          ylab='Counts of CBX2',  )
dev.off()



sd(resampledRowMeans)
quantile(resampledRowMeans, c(0.025,0.9725)) 

CBX2_mean <- data.frame(mean(resampledRowMeans), row.names = 'CBX2')
CBX2_sd <- data.frame(sd(resampledRowMeans), row.names='CBX2')
CBX2_qu<- t(data.frame(quantile(resampledRowMeans, c(0.025,0.9725))))
row.names(CBX2_qu) <- 'CBX2'
CBX2_stats <- cbind(CBX2_mean, CBX2_sd, CBX2_qu) 

CBX2_ul <- as.numeric(UL_t[,2]) 
n <- length(CBX2_ul)
B <- 10000

resamples <- matrix(sample(CBX2_ul,n*B,replace=TRUE),B,n)

resampledRowMeans <- apply(resamples, 1, mean) 

png('CBX2_ul_hist_row_means.png', width=500, height=600)
histogram(resampledRowMeans, type='count', xlab='10,000 Simulated CBX2 in UL Row Means',
          ylab='Counts of CBX2 in UL',  )
dev.off()


sd(resampledRowMeans)
quantile(resampledRowMeans, c(0.025,0.9725)) 

CBX2_ul_mean <- data.frame(mean(resampledRowMeans), row.names = 'CBX2_ul')
CBX2_ul_sd <- data.frame(sd(resampledRowMeans), row.names='CBX2_ul')
CBX2_ul_qu<- t(data.frame(quantile(resampledRowMeans, c(0.025,0.9725))))
row.names(CBX2_qu) <- 'CBX2_ul'
CBX2_ul_stats <- cbind(CBX2_ul_mean, CBX2_ul_sd, CBX2_ul_qu)

CBX2 <- rbind(CBX2_stats, CBX2_ul_stats) 

#bind the two gene stats to one table to bind other 14 genes to as each stats is built

geneStats16 <- rbind(ZNF750,CBX2)#4X4

# now get the third gene table of stats
SOCS3 <- as.numeric(nonUL_t[,3]) 
n <- length(SOCS3)
B <- 10000

resamples <- matrix(sample(SOCS3,n*B,replace=TRUE),B,n)

resampledRowMeans <- apply(resamples, 1, mean) 

png('SOCS3_hist_row_means.png', width=500, height=600)
histogram(resampledRowMeans, type='count', xlab='10,000 Simulated SOCS3 Row Means',
          ylab='Counts of SOCS3',  )
dev.off()


sd(resampledRowMeans)
quantile(resampledRowMeans, c(0.025,0.9725)) 

SOCS3_mean <- data.frame(mean(resampledRowMeans), row.names = 'SOCS3')
SOCS3_sd <- data.frame(sd(resampledRowMeans), row.names='SOCS3')
SOCS3_qu<- t(data.frame(quantile(resampledRowMeans, c(0.025,0.9725))))
row.names(SOCS3_qu) <- 'SOCS3'
SOCS3_stats <- cbind(SOCS3_mean, SOCS3_sd, SOCS3_qu) 

SOCS3_ul <- as.numeric(UL_t[,3]) 
n <- length(SOCS3)
B <- 10000

resamples <- matrix(sample(SOCS3_ul,n*B,replace=TRUE),B,n)

resampledRowMeans <- apply(resamples, 1, mean) 

png('SOCS3_ul_hist_row_means.png', width=500, height=600)
histogram(resampledRowMeans, type='count', xlab='10,000 Simulated SOCS3 in UL Row Means',
          ylab='Counts of SOCS3 in UL',  )
dev.off()


sd(resampledRowMeans)
quantile(resampledRowMeans, c(0.025,0.9725)) 

SOCS3_ul_mean <- data.frame(mean(resampledRowMeans), row.names = 'SOCS3_ul')
SOCS3_ul_sd <- data.frame(sd(resampledRowMeans), row.names='SOCS3_ul')
SOCS3_ul_qu<- t(data.frame(quantile(resampledRowMeans, c(0.025,0.9725))))
row.names(SOCS3_ul_qu) <- 'SOCS3_ul'
SOCS3_ul_stats <- cbind(SOCS3_ul_mean, SOCS3_ul_sd, SOCS3_ul_qu) 

SOCS3 <- rbind(SOCS3_stats, SOCS3_ul_stats) 

geneStats16 <- rbind(geneStats16,SOCS3) #6X4


# Now the 4th gene, RAC3:
RAC3 <- as.numeric(nonUL_t[,4]) 
n <- length(RAC3)
B <- 10000

resamples <- matrix(sample(RAC3,n*B,replace=TRUE),B,n)

resampledRowMeans <- apply(resamples, 1, mean) 

png('RAC3_hist_row_means.png', width=500, height=600)
histogram(resampledRowMeans, type='count', xlab='10,000 Simulated RAC3 Row Means',
          ylab='Counts of RAC3',  )
dev.off()


sd(resampledRowMeans)
quantile(resampledRowMeans, c(0.025,0.9725)) 

RAC3_mean <- data.frame(mean(resampledRowMeans), row.names = 'RAC3')
RAC3_sd <- data.frame(sd(resampledRowMeans), row.names='RAC3')
RAC3_qu<- t(data.frame(quantile(resampledRowMeans, c(0.025,0.9725))))
row.names(RAC3_qu) <- 'RAC3'
RAC3_stats <- cbind(RAC3_mean, RAC3_sd, RAC3_qu) 

RAC3_ul <- as.numeric(UL_t[,4]) 
n <- length(RAC3)
B <- 10000

resamples <- matrix(sample(RAC3_ul,n*B,replace=TRUE),B,n)

resampledRowMeans <- apply(resamples, 1, mean) 

png('RAC3_ul_hist_row_means.png', width=500, height=600)
histogram(resampledRowMeans, type='count', xlab='10,000 Simulated RAC3 in UL Row Means',
          ylab='Counts of RAC3 in UL',  )
dev.off()



sd(resampledRowMeans)
quantile(resampledRowMeans, c(0.025,0.9725)) 

RAC3_ul_mean <- data.frame(mean(resampledRowMeans), row.names = 'RAC3_ul')
RAC3_ul_sd <- data.frame(sd(resampledRowMeans), row.names='RAC3_ul')
RAC3_ul_qu<- t(data.frame(quantile(resampledRowMeans, c(0.025,0.9725))))
row.names(RAC3_ul_qu) <- 'RAC3_ul'
RAC3_ul_stats <- cbind(RAC3_ul_mean, RAC3_ul_sd, RAC3_ul_qu) 

RAC3 <- rbind(RAC3_stats, RAC3_ul_stats) 

geneStats16 <- rbind(geneStats16,RAC3) #8X4


# Now the fifth gene KDELR3:
KDELR3 <- as.numeric(nonUL_t[,5]) 
n <- length(KDELR3)
B <- 10000

resamples <- matrix(sample(KDELR3,n*B,replace=TRUE),B,n)

resampledRowMeans <- apply(resamples, 1, mean) 

png('KDELR3_hist_row_means.png', width=500, height=600)
histogram(resampledRowMeans, type='count', xlab='10,000 Simulated KDELR3 Row Means',
          ylab='Counts of CANT1',  )
dev.off()



sd(resampledRowMeans)
quantile(resampledRowMeans, c(0.025,0.9725)) 

KDELR3_mean <- data.frame(mean(resampledRowMeans), row.names = 'KDELR3')
KDELR3_sd <- data.frame(sd(resampledRowMeans), row.names='KDELR3')
KDELR3_qu<- t(data.frame(quantile(resampledRowMeans, c(0.025,0.9725))))
row.names(KDELR3_qu) <- 'KDELR3'
KDELR3_stats <- cbind(KDELR3_mean, KDELR3_sd, KDELR3_qu) 

KDELR3_ul <- as.numeric(UL_t[,5]) 
n <- length(KDELR3)
B <- 10000

resamples <- matrix(sample(KDELR3_ul,n*B,replace=TRUE),B,n)

resampledRowMeans <- apply(resamples, 1, mean) 

png('KDELR3_ul_hist_row_means.png', width=500, height=600)
histogram(resampledRowMeans, type='count', xlab='10,000 Simulated KDELR3 in UL Row Means',
          ylab='Counts of KDELR3 in UL',  )
dev.off()



sd(resampledRowMeans)
quantile(resampledRowMeans, c(0.025,0.9725)) 

KDELR3_ul_mean <- data.frame(mean(resampledRowMeans), row.names = 'KDELR3_ul')
KDELR3_ul_sd <- data.frame(sd(resampledRowMeans), row.names='KDELR3_ul')
KDELR3_ul_qu<- t(data.frame(quantile(resampledRowMeans, c(0.025,0.9725))))
row.names(KDELR3_ul_qu) <- 'KDELR3_ul'
KDELR3_ul_stats <- cbind(KDELR3_ul_mean, KDELR3_ul_sd, KDELR3_ul_qu) 

KDELR3 <- rbind(KDELR3_stats, KDELR3_ul_stats) 

geneStats16 <- rbind(geneStats16,KDELR3) #10X4

# Now the sixth gene GRIP1:
GRIP1 <- as.numeric(nonUL_t[,6]) 
n <- length(GRIP1)
B <- 10000

resamples <- matrix(sample(GRIP1,n*B,replace=TRUE),B,n)

resampledRowMeans <- apply(resamples, 1, mean) 

png('GRIP1_hist_row_means.png', width=500, height=600)
histogram(resampledRowMeans, type='count', xlab='10,000 Simulated GRIP1 Row Means',
          ylab='Counts of GRIP1',  )
dev.off()



sd(resampledRowMeans)
quantile(resampledRowMeans, c(0.025,0.9725)) 

GRIP1_mean <- data.frame(mean(resampledRowMeans), row.names = 'GRIP1')
GRIP1_sd <- data.frame(sd(resampledRowMeans), row.names='GRIP1')
GRIP1_qu<- t(data.frame(quantile(resampledRowMeans, c(0.025,0.9725))))
row.names(GRIP1_qu) <- 'GRIP1'
GRIP1_stats <- cbind(GRIP1_mean, GRIP1_sd, GRIP1_qu) 

GRIP1_ul <- as.numeric(UL_t[,6]) 
n <- length(GRIP1)
B <- 10000

resamples <- matrix(sample(GRIP1_ul,n*B,replace=TRUE),B,n)

resampledRowMeans <- apply(resamples, 1, mean) 

png('GRIP1_ul_hist_row_means.png', width=500, height=600)
histogram(resampledRowMeans, type='count', xlab='10,000 Simulated GRIP1 in UL Row Means',
          ylab='Counts of GRIP1 in UL',  )
dev.off()


sd(resampledRowMeans)
quantile(resampledRowMeans, c(0.025,0.9725)) 

GRIP1_ul_mean <- data.frame(mean(resampledRowMeans), row.names = 'GRIP1_ul')
GRIP1_ul_sd <- data.frame(sd(resampledRowMeans), row.names='GRIP1_ul')
GRIP1_ul_qu<- t(data.frame(quantile(resampledRowMeans, c(0.025,0.9725))))
row.names(GRIP1_ul_qu) <- 'IRF7_ul'
GRIP1_ul_stats <- cbind(GRIP1_ul_mean, GRIP1_ul_sd, GRIP1_ul_qu) 

GRIP1 <- rbind(GRIP1_stats, GRIP1_ul_stats) 

geneStats16 <- rbind(geneStats16,GRIP1) #12X4

# Now the seventh gene PYCR1:
PYCR1 <- as.numeric(nonUL_t[,7]) 
n <- length(PYCR1)
B <- 10000

resamples <- matrix(sample(PYCR1,n*B,replace=TRUE),B,n)

resampledRowMeans <- apply(resamples, 1, mean) 

png('PYCR1_hist_row_means.png', width=500, height=600)
histogram(resampledRowMeans, type='count', xlab='10,000 Simulated PYCR1 Row Means',
          ylab='Counts of PYCR1',  )
dev.off()


sd(resampledRowMeans)
quantile(resampledRowMeans, c(0.025,0.9725)) 

PYCR1_mean <- data.frame(mean(resampledRowMeans), row.names = 'PYCR1')
PYCR1_sd <- data.frame(sd(resampledRowMeans), row.names='PYCR1')
PYCR1_qu<- t(data.frame(quantile(resampledRowMeans, c(0.025,0.9725))))
row.names(PYCR1_qu) <- 'PYCR1'
PYCR1_stats <- cbind(PYCR1_mean, PYCR1_sd, PYCR1_qu) 

PYCR1_ul <- as.numeric(UL_t[,7]) 
n <- length(PYCR1)
B <- 10000

resamples <- matrix(sample(PYCR1_ul,n*B,replace=TRUE),B,n)

resampledRowMeans <- apply(resamples, 1, mean) 

png('PYCR1_ul_hist_row_means.png', width=500, height=600)
histogram(resampledRowMeans, type='count', xlab='10,000 Simulated PYCR1 in UL Row Means',
          ylab='Counts of PYCR1 in UL',  )
dev.off()


sd(resampledRowMeans)
quantile(resampledRowMeans, c(0.025,0.9725)) 

PYCR1_ul_mean <- data.frame(mean(resampledRowMeans), row.names = 'PYCR1_ul')
PYCR1_ul_sd <- data.frame(sd(resampledRowMeans), row.names='PYCR1_ul')
PYCR1_ul_qu<- t(data.frame(quantile(resampledRowMeans, c(0.025,0.9725))))
row.names(PYCR1_ul_qu) <- 'PYCR1_ul'
PYCR1_ul_stats <- cbind(PYCR1_ul_mean, PYCR1_ul_sd, PYCR1_ul_qu) 

PYCR1 <- rbind(PYCR1_stats, PYCR1_ul_stats) 

geneStats16 <- rbind(geneStats16,PYCR1) #14X4

# Now the eight gene TH:
TH <- as.numeric(nonUL_t[,8]) 
n <- length(TH)
B <- 10000

resamples <- matrix(sample(TH,n*B,replace=TRUE),B,n)

resampledRowMeans <- apply(resamples, 1, mean) 

png('TH_hist_row_means.png', width=500, height=600)
histogram(resampledRowMeans, type='count', xlab='10,000 Simulated TH Row Means',
          ylab='Counts of TH',  )
dev.off()


sd(resampledRowMeans)
quantile(resampledRowMeans, c(0.025,0.9725)) 

TH_mean <- data.frame(mean(resampledRowMeans), row.names = 'TH')
TH_sd <- data.frame(sd(resampledRowMeans), row.names='TH')
TH_qu<- t(data.frame(quantile(resampledRowMeans, c(0.025,0.9725))))
row.names(TH_qu) <- 'TH'
TH_stats <- cbind(TH_mean, TH_sd, TH_qu) 

TH_ul <- as.numeric(UL_t[,8]) 
n <- length(TH)
B <- 10000

resamples <- matrix(sample(TH_ul,n*B,replace=TRUE),B,n)

resampledRowMeans <- apply(resamples, 1, mean) 

png('TH_ul_hist_row_means.png', width=500, height=600)
histogram(resampledRowMeans, type='count', xlab='10,000 Simulated TH in UL Row Means',
          ylab='Counts of TH in UL',  )
dev.off()


sd(resampledRowMeans)
quantile(resampledRowMeans, c(0.025,0.9725)) 

TH_ul_mean <- data.frame(mean(resampledRowMeans), row.names = 'TH_ul')
TH_ul_sd <- data.frame(sd(resampledRowMeans), row.names='TH_ul')
TH_ul_qu<- t(data.frame(quantile(resampledRowMeans, c(0.025,0.9725))))
row.names(TH_ul_qu) <- 'TH_ul'
TH_ul_stats <- cbind(TH_ul_mean, TH_ul_sd, TH_ul_qu) 

TH <- rbind(TH_stats, TH_ul_stats) 

geneStats16 <- rbind(geneStats16,TH) #16X4



# Now the ninth gene CBX7
CBX7 <- as.numeric(nonUL_t[,9]) 
n <- length(CBX7)
B <- 10000

resamples <- matrix(sample(CBX7,n*B,replace=TRUE),B,n)

resampledRowMeans <- apply(resamples, 1, mean) 

png('CBX7_hist_row_means.png', width=500, height=600)
histogram(resampledRowMeans, type='count', xlab='10,000 Simulated CBX7 Row Means',
          ylab='Counts of CBX7',  )
dev.off()


sd(resampledRowMeans)
quantile(resampledRowMeans, c(0.025,0.9725)) 

CBX7_mean <- data.frame(mean(resampledRowMeans), row.names = 'CBX7')
CBX7_sd <- data.frame(sd(resampledRowMeans), row.names='CBX7')
CBX7_qu<- t(data.frame(quantile(resampledRowMeans, c(0.025,0.9725))))
row.names(CBX7_qu) <- 'CBX7'
CBX7_stats <- cbind(CBX7_mean, CBX7_sd, CBX7_qu) 

CBX7_ul <- as.numeric(UL_t[,9]) 
n <- length(CBX7)
B <- 10000

resamples <- matrix(sample(CBX7_ul,n*B,replace=TRUE),B,n)

resampledRowMeans <- apply(resamples, 1, mean) 

png('CBX7_ul_hist_row_means.png', width=500, height=600)
histogram(resampledRowMeans, type='count', xlab='10,000 Simulated CBX7 in UL Row Means',
          ylab='Counts of CBX7 in UL',  )
dev.off()


sd(resampledRowMeans)
quantile(resampledRowMeans, c(0.025,0.9725)) 

CBX7_ul_mean <- data.frame(mean(resampledRowMeans), row.names = 'CBX7_ul')
CBX7_ul_sd <- data.frame(sd(resampledRowMeans), row.names='CBX7_ul')
CBX7_ul_qu<- t(data.frame(quantile(resampledRowMeans, c(0.025,0.9725))))
row.names(CBX7_ul_qu) <- 'CBX7_ul'
CBX7_ul_stats <- cbind(CBX7_ul_mean, CBX7_ul_sd, CBX7_ul_qu) 

CBX7 <- rbind(CBX7_stats, CBX7_ul_stats) 

geneStats16 <- rbind(geneStats16,CBX7) #18X4

# Now the 10th gene ASPSCR1:
ASPSCR1 <- as.numeric(nonUL_t[,10]) 
n <- length(ASPSCR1)
B <- 10000

resamples <- matrix(sample(ASPSCR1,n*B,replace=TRUE),B,n)

resampledRowMeans <- apply(resamples, 1, mean) 

png('ASPSCR1_hist_row_means.png', width=500, height=600)
histogram(resampledRowMeans, type='count', xlab='10,000 Simulated ASPSCR1 Row Means',
          ylab='Counts of ASPSCR1',  )
dev.off()


sd(resampledRowMeans)
quantile(resampledRowMeans, c(0.025,0.9725)) 

ASPSCR1_mean <- data.frame(mean(resampledRowMeans), row.names = 'ASPSCR1')
ASPSCR1_sd <- data.frame(sd(resampledRowMeans), row.names='ASPSCR1')
ASPSCR1_qu<- t(data.frame(quantile(resampledRowMeans, c(0.025,0.9725))))
row.names(ASPSCR1_qu) <- 'ASPSCR1'
ASPSCR1_stats <- cbind(ASPSCR1_mean, ASPSCR1_sd, ASPSCR1_qu) 

ASPSCR1_ul <- as.numeric(UL_t[,10]) 
n <- length(ASPSCR1)
B <- 10000

resamples <- matrix(sample(ASPSCR1_ul,n*B,replace=TRUE),B,n)

resampledRowMeans <- apply(resamples, 1, mean) 

png('ASPSCR1_ul_hist_row_means.png', width=500, height=600)
histogram(resampledRowMeans, type='count', xlab='10,000 Simulated ASPSCR1 in UL Row Means',
          ylab='Counts of ASPSCR1 in UL',  )
dev.off()


sd(resampledRowMeans)
quantile(resampledRowMeans, c(0.025,0.9725)) 

ASPSCR1_ul_mean <- data.frame(mean(resampledRowMeans), row.names = 'ASPSCR1_ul')
ASPSCR1_ul_sd <- data.frame(sd(resampledRowMeans), row.names='ASPSCR1_ul')
ASPSCR1_ul_qu<- t(data.frame(quantile(resampledRowMeans, c(0.025,0.9725))))
row.names(ASPSCR1_ul_qu) <- 'ASPSCR1_ul'
ASPSCR1_ul_stats <- cbind(ASPSCR1_ul_mean, ASPSCR1_ul_sd, ASPSCR1_ul_qu) 

ASPSCR1 <- rbind(ASPSCR1_stats, ASPSCR1_ul_stats) 

geneStats16 <- rbind(geneStats16,ASPSCR1) #20X4

#Now the 11th gene BET1L:
BET1L <- as.numeric(nonUL_t[,11]) 
n <- length(BET1L)
B <- 10000

resamples <- matrix(sample(BET1L,n*B,replace=TRUE),B,n)

resampledRowMeans <- apply(resamples, 1, mean) 

png('BET1L_hist_row_means.png', width=500, height=600)
histogram(resampledRowMeans, type='count', xlab='10,000 Simulated BET1L Row Means',
          ylab='Counts of BET1L',  )
dev.off()


sd(resampledRowMeans)
quantile(resampledRowMeans, c(0.025,0.9725)) 

BET1L_mean <- data.frame(mean(resampledRowMeans), row.names = 'BET1L')
BET1L_sd <- data.frame(sd(resampledRowMeans), row.names='BET1L')
BET1L_qu<- t(data.frame(quantile(resampledRowMeans, c(0.025,0.9725))))
row.names(BET1L_qu) <- 'BET1L'
BET1L_stats <- cbind(BET1L_mean, BET1L_sd, BET1L_qu) 

BET1L_ul <- as.numeric(UL_t[,11]) 
n <- length(BET1L)
B <- 10000

resamples <- matrix(sample(BET1L_ul,n*B,replace=TRUE),B,n)

resampledRowMeans <- apply(resamples, 1, mean) 

png('BET1L_ul_hist_row_means.png', width=500, height=600)
histogram(resampledRowMeans, type='count', xlab='10,000 Simulated BET1L in UL Row Means',
          ylab='Counts of BET1L in UL',  )
dev.off()


sd(resampledRowMeans)
quantile(resampledRowMeans, c(0.025,0.9725)) 

BET1L_ul_mean <- data.frame(mean(resampledRowMeans), row.names = 'BET1L_ul')
BET1L_ul_sd <- data.frame(sd(resampledRowMeans), row.names='BET1L_ul')
BET1L_ul_qu<- t(data.frame(quantile(resampledRowMeans, c(0.025,0.9725))))
row.names(BET1L_ul_qu) <- 'BET1L_ul'
BET1L_ul_stats <- cbind(BET1L_ul_mean, BET1L_ul_sd, BET1L_ul_qu) 

BET1L <- rbind(BET1L_stats, BET1L_ul_stats) 

geneStats16 <- rbind(geneStats16,BET1L) #22X4


# Now the 12th gene HMGA2:
HMGA2 <- as.numeric(nonUL_t[,12]) 
n <- length(HMGA2)
B <- 10000

resamples <- matrix(sample(HMGA2,n*B,replace=TRUE),B,n)

resampledRowMeans <- apply(resamples, 1, mean) 

png('HMGA2_hist_row_means.png', width=500, height=600)
histogram(resampledRowMeans, type='count', xlab='10,000 Simulated HMGA2 Row Means',
          ylab='Counts of HMGA2',  )
dev.off()


sd(resampledRowMeans)
quantile(resampledRowMeans, c(0.025,0.9725)) 

HMGA2_mean <- data.frame(mean(resampledRowMeans), row.names = 'HMGA2')
HMGA2_sd <- data.frame(sd(resampledRowMeans), row.names='HMGA2')
HMGA2_qu<- t(data.frame(quantile(resampledRowMeans, c(0.025,0.9725))))
row.names(HMGA2_qu) <- 'HMGA2'
HMGA2_stats <- cbind(HMGA2_mean, HMGA2_sd, HMGA2_qu) 

HMGA2_ul <- as.numeric(UL_t[,12]) 
n <- length(HMGA2)
B <- 10000

resamples <- matrix(sample(HMGA2_ul,n*B,replace=TRUE),B,n)

resampledRowMeans <- apply(resamples, 1, mean) 

png('HMGA2_ul_hist_row_means.png', width=500, height=600)
histogram(resampledRowMeans, type='count', xlab='10,000 Simulated HMGA2 in UL Row Means',
          ylab='Counts of HMGA2 in UL',  )
dev.off()


sd(resampledRowMeans)
quantile(resampledRowMeans, c(0.025,0.9725)) 

HMGA2_ul_mean <- data.frame(mean(resampledRowMeans), row.names = 'HMGA2_ul')
HMGA2_ul_sd <- data.frame(sd(resampledRowMeans), row.names='HMGA2_ul')
HMGA2_ul_qu<- t(data.frame(quantile(resampledRowMeans, c(0.025,0.9725))))
row.names(HMGA2_ul_qu) <- 'HMGA2_ul'
HMGA2_ul_stats <- cbind(HMGA2_ul_mean, HMGA2_ul_sd, HMGA2_ul_qu) 

HMGA2 <- rbind(HMGA2_stats, HMGA2_ul_stats) 

geneStats16 <- rbind(geneStats16,HMGA2) #24X4

# Now the 13th gene CYTH4:
CYTH4 <- as.numeric(nonUL_t[,13]) 
n <- length(CYTH4)
B <- 10000

resamples <- matrix(sample(CYTH4,n*B,replace=TRUE),B,n)

resampledRowMeans <- apply(resamples, 1, mean) 

png('CYTH4_hist_row_means.png', width=500, height=600)
histogram(resampledRowMeans, type='count', xlab='10,000 Simulated CYTH4 Row Means',
          ylab='Counts of CYTH4',  )
dev.off()


sd(resampledRowMeans)
quantile(resampledRowMeans, c(0.025,0.9725)) 

CYTH4_mean <- data.frame(mean(resampledRowMeans), row.names = 'CYTH4')
CYTH4_sd <- data.frame(sd(resampledRowMeans), row.names='CYTH4')
CYTH4_qu<- t(data.frame(quantile(resampledRowMeans, c(0.025,0.9725))))
row.names(CYTH4_qu) <- 'CYTH4'
CYTH4_stats <- cbind(CYTH4_mean, CYTH4_sd, CYTH4_qu) 

CYTH4_ul <- as.numeric(UL_t[,13]) 
n <- length(CYTH4)
B <- 10000

resamples <- matrix(sample(CYTH4_ul,n*B,replace=TRUE),B,n)

resampledRowMeans <- apply(resamples, 1, mean) 

png('CYTH4_ul_hist_row_means.png', width=500, height=600)
histogram(resampledRowMeans, type='count', xlab='10,000 Simulated CYTH4 in UL Row Means',
          ylab='Counts of CYTH4 in UL',  )
dev.off()


sd(resampledRowMeans)
quantile(resampledRowMeans, c(0.025,0.9725)) 

CYTH4_ul_mean <- data.frame(mean(resampledRowMeans), row.names = 'CYTH4_ul')
CYTH4_ul_sd <- data.frame(sd(resampledRowMeans), row.names='CYTH4_ul')
CYTH4_ul_qu<- t(data.frame(quantile(resampledRowMeans, c(0.025,0.9725))))
row.names(CYTH4_ul_qu) <- 'CYTH4_ul'
CYTH4_ul_stats <- cbind(CYTH4_ul_mean, CYTH4_ul_sd, CYTH4_ul_qu) 

CYTH4 <- rbind(CYTH4_stats, CYTH4_ul_stats) 

geneStats16 <- rbind(geneStats16,CYTH4) #26X4

# Now the 14th gene CCDC57:
CCDC57 <- as.numeric(nonUL_t[,14]) 
n <- length(CCDC57)
B <- 10000

resamples <- matrix(sample(CCDC57,n*B,replace=TRUE),B,n)

resampledRowMeans <- apply(resamples, 1, mean) 

png('CcDC57_hist_row_means.png', width=500, height=600)
histogram(resampledRowMeans, type='count', xlab='10,000 Simulated CCDC57 Row Means',
          ylab='Counts of CCDC57',  )
dev.off()


sd(resampledRowMeans)
quantile(resampledRowMeans, c(0.025,0.9725)) 

CCDC57_mean <- data.frame(mean(resampledRowMeans), row.names = 'CCDC57')
CCDC57_sd <- data.frame(sd(resampledRowMeans), row.names='CCDC57')
CCDC57_qu<- t(data.frame(quantile(resampledRowMeans, c(0.025,0.9725))))
row.names(CCDC57_qu) <- 'CCDC57'
CCDC57_stats <- cbind(CCDC57_mean, CCDC57_sd, CCDC57_qu) 

CCDC57_ul <- as.numeric(UL_t[,14]) 
n <- length(CCDC57)
B <- 10000

resamples <- matrix(sample(CCDC57_ul,n*B,replace=TRUE),B,n)

resampledRowMeans <- apply(resamples, 1, mean) 

png('CcDC57_ul_hist_row_means.png', width=500, height=600)
histogram(resampledRowMeans, type='count', xlab='10,000 Simulated CCDC57 in UL Row Means',
          ylab='Counts of CCDC57 in UL',  )
dev.off()


sd(resampledRowMeans)
quantile(resampledRowMeans, c(0.025,0.9725)) 

CCDC57_ul_mean <- data.frame(mean(resampledRowMeans), row.names = 'CCDC57_ul')
CCDC57_ul_sd <- data.frame(sd(resampledRowMeans), row.names='CCDC57_ul')
CCDC57_ul_qu<- t(data.frame(quantile(resampledRowMeans, c(0.025,0.9725))))
row.names(CCDC57_ul_qu) <- 'CCDC57_ul'
CCDC57_ul_stats <- cbind(CCDC57_ul_mean, CCDC57_ul_sd, CCDC57_ul_qu) 

CCDC57 <- rbind(CCDC57_stats, CCDC57_ul_stats) 

geneStats16 <- rbind(geneStats16,CCDC57) #28X4

# Now the 15th gene FASN:
FASN <- as.numeric(nonUL_t[,15]) 
n <- length(FASN)
B <- 10000

resamples <- matrix(sample(FASN,n*B,replace=TRUE),B,n)

resampledRowMeans <- apply(resamples, 1, mean) 

png('FASN_hist_row_means.png', width=500, height=600)
histogram(resampledRowMeans, type='count', xlab='10,000 Simulated FASN Row Means',
          ylab='Counts of FASN',  )
dev.off()


sd(resampledRowMeans)
quantile(resampledRowMeans, c(0.025,0.9725)) 

FASN_mean <- data.frame(mean(resampledRowMeans), row.names = 'FASN')
FASN_sd <- data.frame(sd(resampledRowMeans), row.names='FASN')
FASN_qu<- t(data.frame(quantile(resampledRowMeans, c(0.025,0.9725))))
row.names(FASN_qu) <- 'FASN'
FASN_stats <- cbind(FASN_mean, FASN_sd, FASN_qu) 

FASN_ul <- as.numeric(UL_t[,15]) 
n <- length(FASN)
B <- 10000

resamples <- matrix(sample(FASN_ul,n*B,replace=TRUE),B,n)

resampledRowMeans <- apply(resamples, 1, mean) 

png('FASN_ul_hist_row_means.png', width=500, height=600)
histogram(resampledRowMeans, type='count', xlab='10,000 Simulated FASN in UL Row Means',
          ylab='Counts of FASN in UL',  )
dev.off()


sd(resampledRowMeans)
quantile(resampledRowMeans, c(0.025,0.9725)) 

FASN_ul_mean <- data.frame(mean(resampledRowMeans), row.names = 'FASN_ul')
FASN_ul_sd <- data.frame(sd(resampledRowMeans), row.names='FASN_ul')
FASN_ul_qu<- t(data.frame(quantile(resampledRowMeans, c(0.025,0.9725))))
row.names(FASN_ul_qu) <- 'FASN_ul'
FASN_ul_stats <- cbind(FASN_ul_mean, FASN_ul_sd, FASN_ul_qu) 

FASN <- rbind(FASN_stats, FASN_ul_stats) 

geneStats16 <- rbind(geneStats16,FASN) #30X4

# Now the 16th gene TNRC6B:
TNRC6B <- as.numeric(nonUL_t[,16]) 
n <- length(TNRC6B)
B <- 10000

resamples <- matrix(sample(TNRC6B,n*B,replace=TRUE),B,n)

resampledRowMeans <- apply(resamples, 1, mean) 

png('TNRC6B_hist_row_means.png', width=500, height=600)
histogram(resampledRowMeans, type='count', xlab='10,000 Simulated TNRC6B Row Means',
          ylab='Counts of TNRC6B',  )
dev.off()


sd(resampledRowMeans)
quantile(resampledRowMeans, c(0.025,0.9725)) 

TNRC6B_mean <- data.frame(mean(resampledRowMeans), row.names = 'TNRC6B')
TNRC6B_sd <- data.frame(sd(resampledRowMeans), row.names='TNRC6B')
TNRC6B_qu<- t(data.frame(quantile(resampledRowMeans, c(0.025,0.9725))))
row.names(TNRC6B_qu) <- 'TNRC6B'
TNRC6B_stats <- cbind(TNRC6B_mean, TNRC6B_sd, TNRC6B_qu) 

TNRC6B_ul <- as.numeric(UL_t[,16]) 
n <- length(TNRC6B)
B <- 10000

resamples <- matrix(sample(TNRC6B_ul,n*B,replace=TRUE),B,n)

resampledRowMeans <- apply(resamples, 1, mean) 

png('TNRC6B_ul_hist_row_means.png', width=500, height=600)
histogram(resampledRowMeans, type='count', xlab='10,000 Simulated TNRC6B in UL Row Means',
          ylab='Counts of TNRC6B in UL',  )
dev.off()


sd(resampledRowMeans)
quantile(resampledRowMeans, c(0.025,0.9725)) 

TNRC6B_ul_mean <- data.frame(mean(resampledRowMeans), row.names = 'TNRC6B_ul')
TNRC6B_ul_sd <- data.frame(sd(resampledRowMeans), row.names='TNRC6B_ul')
TNRC6B_ul_qu<- t(data.frame(quantile(resampledRowMeans, c(0.025,0.9725))))
row.names(TNRC6B_ul_qu) <- 'TNRC6B_ul'
TNRC6B_ul_stats <- cbind(TNRC6B_ul_mean, TNRC6B_ul_sd, TNRC6B_ul_qu) 

TNRC6B <- rbind(TNRC6B_stats, TNRC6B_ul_stats) 

geneStats16 <- rbind(geneStats16,TNRC6B) #32X4

# change the '2.5%' and '97,25%' column names
colnames(geneStats16)
# [1] "mean.resampledRowMeans." "sd.resampledRowMeans."   "2.5%"                   
# [4] "97.25%"

colnames(geneStats16) <- c('simulatedMean10k','simulatedSD10K', 'leftTail2.5',
                           'rightTail97.25')
ul <- grep('_ul', row.names(geneStats16)) #16 obs
ULStats <- geneStats16[ul,] #16X4
nonULStats <- geneStats16[-ul,]#16X4

non <- data.frame(c(rep('nonUL',16)))
u <- data.frame(c(rep('UL',16)))

UL_stats <- cbind(ULStats,u) #16X5, adds a field identifying gene as UL
colnames(UL_stats)[5] <- 'ulStatus'
row.names(UL_stats) <- gsub('_ul','', row.names(UL_stats))
UL_stats[6] <- row.names(UL_stats)
colnames(UL_stats)[6] <- 'gene'

nonUL_stats <- cbind(nonULStats,non) #16X5, adds a field identifying gene as nonUL
colnames(nonUL_stats)[5] <- 'ulStatus'
nonUL_stats[6] <- row.names(nonUL_stats)
colnames(nonUL_stats)[6] <- 'gene'

Stats16 <- rbind(UL_stats, nonUL_stats) #32X6

# save the stats on the top 10 most differentially expressed genes in UL compared to nonUL
# and the most differentially expressed genes ubiquitous to UL research studies
write.csv(Stats16,'Stats16.csv', row.names=FALSE)
#########################################################################################
#########################################################################################
#########################################################################################


# The most differentially expressed genes in UL, the top 10, and the most differentially
# expressed six genes in decreasing order of most expressed in UL in up/down change:
Stats16 <- read.csv('Stats16.csv', sep=',', header=TRUE, na.strings=c('','NA'))#32X6

# Read in the table of genes that keeps track of whether the gene is part of the majority
majority <- read.csv("MemberGviz_130_141.csv", sep=',', header=TRUE, 
                     na.strings=c('','NA'))

# Merge these two tables together to add Gviz data and samples to the stats of the 16 genes
StatsGviz16 <- merge(Stats16, majority, by.x='gene', by.y='genes') #32X146##CYTOBAND##included

write.csv(StatsGviz16, 'Gviz16.csv', row.names=FALSE)#32X146##CYTOBAND##included

colnames(StatsGviz16)
# [1] "gene"                          "simulatedMean10k"              "simulatedSD10K"               
# [4] "leftTail2.5"                   "rightTail97.25"                "ulStatus"                     
# [7] "chromosome"                    "type"                          "all"                          
# [10] "up"                            "down"                          "majority"                     
# [13] "start"                         "end"                           "width"                        
# [16] "strand"                        "gene.y"                        "transcript"                   
# [19] "GENE"                          "GENE_NAME"                     "CYTOBAND"                     
# [22] "DESCRIPTION"                   "nonUL_Mean"                    "UL_Mean"                      
# [25] "Difference_UL_minus_non_means" "gsm1667144"                    "gsm1667145"                   
# [28] "gsm1667146"                    "gsm336252"                     "gsm336253"                    
# [31] "gsm336254"                     "gsm336255"                     "gsm336256"                    
# ... 
# [136] "gsm38695ul"                    "gsm9093ul"                     "gsm9094ul"                    
# [139] "gsm9095ul"                     "gsm9096ul"                     "gsm9097ul"                    
# [142] "gsm569429ul"                   "gsm569430ul"                   "gsm569431ul"                  
# [145] "gsm569432ul"                   "gsm569433ul" 

#########################################################################################
#########################################################################################
# Using lattice, add a pairwise gene comparison of the 16 genes in different sets
# those up, down, majority, and minority


# now get those genes part of the majority of DE change as up or down in UL on each chromosome

StatsGviz16 <- read.csv('Gviz16.csv', sep=',', header=TRUE, 
                        na.strings=c('', 'NA')) #32X146

library(dplyr)

majorGroup <- filter(StatsGviz16, majority==TRUE) #18X146, 9 genes
unique(majorGroup$gene)
# [1] GRIP1  HMGA2  KDELR3 TNRC6B

up_maj <- filter(majorGroup, type=='up')
unique(up_maj$gene)#majority of genes up-regulated
# [1] GRIP1  HMGA2  KDELR3 TNRC6B

up_maj <- up_maj[,-c(2:25)]
up_maj <- up_maj[!duplicated(up_maj$gene),]
up_maj <- t(up_maj)
names <- as.character(up_maj[1,])
up_maj <- data.frame(up_maj)
up_maj <- up_maj[-1,]
colnames(up_maj) <- names #121X7 genes up regulated in majority of each chromosome


down_maj <- filter(majorGroup, type=='down')
unique(down_maj$gene)#majority of genes down-regulated
#SOCS3 HMGA2

down_maj <- down_maj[,-c(2:25)]
down_maj <- down_maj[!duplicated(down_maj$gene),]
down_maj <- t(down_maj)
names <- as.character(down_maj[1,])
down_maj <- data.frame(down_maj)
down_maj <- down_maj[-1,]
colnames(down_maj) <- names#121X2


maj_latt <- majorGroup[,-c(2:25)]#18X122
maj_latt <- maj_latt[!duplicated(majorGroup$gene),]#9X122, gets rid of duplicates
maj_latt <- t(maj_latt)#122X9
names <- as.character(maj_latt[1,])#9 values
maj_latt <- data.frame(maj_latt)#122X9
maj_latt <- maj_latt[-1,]#121X9
colnames(maj_latt) <- names

minorGroup <- filter(StatsGviz16, majority==FALSE) #14X146, 7 genes
unique(minorGroup$gene)
# BET1L    CCDC57   FASN     FSCN2    IRF7     ASPSCR1     SLC38A10

up_min <- filter(minorGroup, type=='up')
unique(up_min$gene)#minority of genes up-regulated
#BET1L IRF7  ASPSCR1

up_min <- up_min[,-c(2:25)]
up_min <- up_min[!duplicated(up_min$gene),]
up_min <- t(up_min)
names <- as.character(up_min[1,])
up_min <- data.frame(up_min)
up_min <- up_min[-1,]
colnames(up_min) <- names #121X3 genes up regulated in minority of each chromosome


down_min <- filter(minorGroup, type=='down')
unique(down_min$gene)#minority of genes down-regulated
#CCDC57   FASN     FSCN2    SLC38A10

down_min <- down_min[,-c(2:25)]
down_min <- down_min[!duplicated(down_min$gene),]
down_min <- t(down_min)
names <- as.character(down_min[1,])
down_min <- data.frame(down_min)
down_min <- down_min[-1,]
colnames(down_min) <- names #121X4

min_latt <- minorGroup[,-c(2:25)]#14X122
min_latt <- min_latt[!duplicated(minorGroup$gene),]#7X122, gets rid of duplicates
min_latt <- t(min_latt)#122X7
names <- as.character(min_latt[1,])#7 values
min_latt <- data.frame(min_latt)#122X7
min_latt <- min_latt[-1,]#121X7
colnames(min_latt) <- names

# Create a pairwise comparison plot with lattice
library(lattice)

png('up_maj_splom.png', width=1200, height=1200)
splom(up_maj,
      main="Genes in Majority Group for Up Gene Expression in Each Chromosome")
dev.off()

# png('up_min_splom.png', width=900, height=900)#terrible looks like swastikas
# splom(up_min,
#       main="Genes in Minority Group for Up Gene Expression in Each Chromosome")
# dev.off()

png('down_maj_splom.png', width=900, height=900)
splom(down_maj,
      main="Genes in Majority Group for Down Gene Expression in Each Chromosome")
dev.off()

# png('down_min_splom.png', width=900, height=900)#terrible looks like swastikas
# splom(down_min,
#       main="Genes in Minority Group for Down Gene Expression in Each Chromosome")
# dev.off()

# note for the most part there is a positive correlation, except for the randomness
# in the minority group of up-regulated genes, and somewhat for the down regulated 
# majority. For mid level values of the other two groups there is more positive correlation
# but for low and high values there is more variation

#########################################################################################
#########################################################################################
# Build a heatmap of the differential gene expression with heatmaply
library(dplyr)
library(heatmaply)

StatsGviz16 <- read.csv('Gviz16.csv', sep=',', header=TRUE, 
                        na.strings=c('', 'NA')) #32X146

UL <- filter(StatsGviz16, ulStatus=='UL')#16X146
nonUL <- filter(StatsGviz16, ulStatus=='nonUL') #16X146

# sort each by simulated standard deviation from lowest to highest
UL_sd_sort <- UL[order(UL$simulatedSD10K, decreasing=FALSE),]
UL_sd_sort$gene
# [1] TNRC6B   ASPSCR1     IRF7     FASN     CBX2   NOL12    SLC38A10 CYTH4    SOCS3    HMGA2   
# [11] ARHGDIA  FSCN2    BET1L    CANT1    SLC25A10 CCDC57

nonUL_sd_sort <- nonUL[order(nonUL$simulatedSD10K),] #16X146 also increasing
nonUL_sd_sort$gene
# [1] TNRC6B   IRF7     ASPSCR1     NOL12    FASN     CBX2   SLC38A10 CANT1    BET1L    CYTH4   
# [11] ARHGDIA  FSCN2    SLC25A10 SOCS3    HMGA2    CCDC57

#sort by simulated mean
UL_mean <- UL[order(UL$simulatedMean10k, decreasing=TRUE),]
UL_mean$gene
# [1] CYTH4    SLC38A10 IRF7     TNRC6B   NOL12    ASPSCR1     SOCS3    CCDC57   CBX2   BET1L   
# [11] FASN     SLC25A10 ARHGDIA  CANT1    HMGA2    FSCN2 

nonUL_mean <- nonUL[order(nonUL$simulatedMean10k, decreasing=TRUE),]
nonUL_mean$gene
# [1] SLC38A10 CYTH4    TNRC6B   IRF7     NOL12    SOCS3    ASPSCR1     CCDC57   CBX2   BET1L   
# [11] FASN     SLC25A10 ARHGDIA  CANT1    HMGA2    FSCN2

# from above the first 7 genes move around by neighboring genes in mean simulated expressions
# between the two groups of UL or nonUL, but the last 9 genes stay the same for ranked
# gene expression levels.

# create a field for differential gene expression between the two groups of UL or nonUL
# simulated means

meanSimUL <- data.frame(UL$simulatedMean10k) #16X1
meanSimNonUL <- data.frame(nonUL$simulatedMean10k)
meanSimAll <- cbind(meanSimUL,meanSimNonUL)
row.names(meanSimAll) <- as.character(UL$gene)
DE <- mutate(meanSimAll, 
             DE=meanSimAll$UL.simulatedMean10k-meanSimAll$nonUL.simulatedMean10k)#16X3
row.names(DE) <- row.names(meanSimAll)
SimDE <- data.frame(DE$DE)
row.names(SimDE) <- as.character(UL$gene)
colnames(SimDE) <- 'DE'
genes16_non <- nonUL[,2:3] #keep simulated mean and standard deviation per gene in nonUL
colnames(genes16_non) <- c('nonUL_means', 'nonUL_sd')
genes16_ul <- UL[,2:3] #keep simulated mean and sd in UL
colnames(genes16_ul) <- c('ul_means','ul_sd')
SimDE16 <- cbind(SimDE, genes16_non, genes16_ul)#DE, non means and sd, ul means and sd,16X5
SimDE16 <- data.frame(SimDE16)
write.csv(SimDE16, 'SimDE16.csv', row.names=TRUE)

png('heatmapTop16_DE_mean_sd_sims.png', width=500, height=700)
heatmap(as.matrix(SimDE16), margins = c(8,7),#x-axis, y-axis label distance from grid
        cexCol = 1.1, cexRow = 1.1, keep.dendro=FALSE,
        main=NULL,
        xlab='Statistical Simulations for UL and nonUL',
        ylab='Genes Most* DE for UL')

# the above heatmap doesn't have a title, nor a legend
#base because heatmaply() uses viewer not plots panel
dev.off()

# Load this package in Bioconductor for 'ComplexHeatmap'
# source('http://www.bioconductor.org/biocLite.R')
# BiocManager::install('ComplexHeatmap')
library(ComplexHeatmap)

#see the heatmap for the stats of 16 genes, DE, means, and sd using Heatmap()
png('heatmap_SimDE16_stats.png', width=500, height=500)
Heatmap(as.matrix(t(SimDE16)), cluster_columns = FALSE,
        row_names_side='left',
        row_dend_side = 'left',
        clustering_distance_rows='maximum',#'euclidean' (default), 'manhattan', 
        #'carberra', 'bing'
        #, 'minkowski', 'pearson', 'spearman', 'kendall'
        
        clustering_method_rows='single') #'single', 'complete', 'average' (UPGMA),
# 'mcquitty'(UPGMA), 'median'(UPGMA), 'central' (UPGMA)
# a legend included, 16 genes as rows and 5 columns of stats

round(range(SimDE16$DE),3) #[1] -0.821  0.795
dev.off()


Samples16genes <- UL[,-c(1:25)]
row.names(Samples16genes) <- as.character(UL$gene)# 16 genes X 121 samples
Samples16genes <- as.matrix(Samples16genes)
Samples16genes <- data.frame(Samples16genes)

write.csv(Samples16genes,'Samples16genes.csv', row.names=TRUE)

##REPORT##
png('heatmap_sample16genesLRG.png', width=500,height=800)
Heatmap(as.matrix(t(Samples16genes)),#shows side label of 121 sample IDs, and 16 genes
        #as columns, cluttered, red-blue color scheme, most to least values of gene expression
        cluster_columns = TRUE, #default
        row_names_side='left', 
        show_row_names=FALSE,
        row_dend_side = 'left',
        clustering_distance_rows='euclidean', #(default) 
        #'carberra', 'bing','maximum' , 'manhattan',
        #, 'minkowski', 'pearson', 'spearman', 'kendall'
        clustering_method_rows='single', #'single', 'complete', 'average' (UPGMA),
        # 'mcquitty'(UPGMA), 'median'(UPGMA), 'central' (UPGMA))
        row_title_side='right',
        column_title_side='bottom',
        row_title_gp=gpar(fontsize=6),
        #row_title = character(heatmap 121 samples),
        heatmap_height=1200)

dev.off()

#########################################################################################
#########################################################################################
# Use ggplot2 on the 32 UL and nonUL samples

StatsGviz16 <- read.csv('Gviz16.csv', sep=',', header=TRUE, na.strings=c('','NA'))
samples <- StatsGviz16[,c(1,26:146)]#32X122

library(ggplot2)


meta <- StatsGviz16[,c(1, 7:16)]#32X11
colnames(meta)
# [1] "gene"       "chromosome" "type"       "all"        "up"         "down"      
# [7] "majority"   "start"      "end"        "width"      "strand" 

Meta <- meta[!duplicated(meta$gene),]#16X11
gene <- data.frame(row.names(SimDE16))
colnames(gene) <- 'gene'
DE <- cbind(gene,SimDE16)
Meta_DE_Stats <- merge(DE, Meta, by.x='gene', by.y='gene')

#add the cytoband location
cyto <- read.csv('MemberGviz_130_141.csv', sep=',', header=TRUE)
cyto <- cyto[,c(1,16)]
Meta_DE_Stats <- merge(cyto, Meta_DE_Stats, by.x='genes', by.y='gene')

write.csv(Meta_DE_Stats, 'metaGviz_DE_Stats_16.csv', row.names=FALSE)#16X17

##REPORT##
png('Sim_UL_nonUL_means_chr.png', width=800, height=600)
g <- ggplot(Meta_DE_Stats, aes(x=Meta_DE_Stats$nonUL_means, y=Meta_DE_Stats$ul_means))
g= g+xlab('Non-UL Simulated Mean Top 10 Plus 6 DE Genes')
x=Meta_DE_Stats$ul_means
y=Meta_DE_Stats$ul_means
g= g+ xlim(min(x)-1,max(x)+1)+ ylim(min(y)-1,max(y)+1)
g= g+ ylab('UL Simulated Mean Top 10 Plus 6 DE Genes')
g= g+ geom_point(aes(colour=CYTOBAND),size=6, alpha=0.9)
g= g+ geom_text(aes(label=genes), size=3,hjust=1, vjust=.5)
g= g+ geom_smooth(method='lm', colour='red')
g
dev.off()


#########################################################################################
#########################################################################################
#########################################################################################
#clear out environment

##REPORT##

# Machine Learning on data

StatsGviz16 <- read.csv('Gviz16.csv', sep=',', header=TRUE, na.strings=c('','NA'))#32X146
Meta_DE_Stats <- read.csv('metaGviz_DE_Stats_16.csv', sep=',', 
                          na.strings=c('','NA'))#16X16

DE <- Meta_DE_Stats[,c(1,2)] #16X2
StatsDE <- merge(DE, StatsGviz16, by.x='gene', by.y='gene') #32X147

# install.packages('caret')
# install.packages('randomForest')
# install.packages('MASS')
# install.packages('gbm')

library(caret)
library(randomForest)
library(MASS)
library(gbm)

# There is a problem with dplyr after those installs above, fix with this:
# install.packages('devtools')
# devtools::install_github("tidyverse/dplyr")
# install.packages('pillar')
library(dplyr)

ul <- filter(StatsDE, StatsDE$ulStatus=='UL') #16X147
non <- filter(StatsDE, StatsDE$ulStatus=='nonUL')#16X147

write.csv(ul, 'ul_16_147.csv', row.names=FALSE)#same tables that include the UL and non
write.csv(non, 'nonUL_16_147.csv', row.names=FALSE)

UL_df <- ul[,c(1,27:147)] #16X122
names <- UL_df[,1]
UL_df <- UL_df[,-1] #16X121
row.names(UL_df) <- names
ul_t <- t(UL_df) #121X16

samples_t <- data.frame(ul_t) #121 samples of UL and nonUL X 16 genes

#add a field that is meta for if the sample is ul or not, by grep() row.names
#there are 51 nonUL and 70 UL samples listed in order
u <- t(c(rep('UL',70)))
n <- t(c(rep('nonUL',51)))

type <- t(cbind(n,u))#121X1
colnames(type) <- toupper('type')

samplesType <- data.frame(cbind(type,samples_t))#121X17

summary(samplesType$TYPE)
# nonUL    UL 
# 51    70 

write.csv(samplesType, 'samplesType.csv', row.names=TRUE)
samplesType <- read.csv('samplesType.csv', header=TRUE, row.names=1)#121X17
colnames(samplesType)
# [1] "TYPE"     "ARHGDIA"  "BET1L"    "CANT1"    "CBX2"   "CCDC57"   "CYTH4"    "FASN"    
# [9] "FSCN2"    "SOCS3"    "HMGA2"    "IRF7"     "NOL12"    "ASPSCR1"     "SLC25A10" "SLC38A10"
# [17] "TNRC6B"

set.seed(189678345) # this will reproduce the same numbers each time sampling is done
# run set.seed value to get exact results, otherwise results differ

# creates a partition of the data by indices
inTrain <- createDataPartition(y=samplesType$SOCS3, p=0.7, list=FALSE)
# lenght of 86 out of 121 samples sampled

trainingSet <- samplesType[inTrain,]#86X16
testingSet <- samplesType[-inTrain,]#35X16

#install.packages('e1071')
library(e1071)

# randomForest, cross-validation (cv) = 5
rfMod <- train(TYPE~., method='rf', data=(trainingSet), 
               trControl=trainControl(method='cv'), number=5)#list of 23
png('RandomForest.png')
plot(rfMod)
dev.off()



# generalizedBoostedModel
gbmMod <- train(TYPE~., method='gbm', data=trainingSet, verbose=FALSE )

png('globalBoostedModel.png')
plot(gbmMod)
dev.off()


# linkage dirichlet allocation model
ldaMod <- train(TYPE~., method='lda', data=trainingSet)

# the ldaMod cannot be plotted, 'no tuning parameters for this model'

# run predictions on the testing set
predRF <- predict(rfMod, testingSet)
predGbm <- predict(gbmMod, testingSet)
predlda <- predict(ldaMod, testingSet)

predDF <- data.frame(predRF, predGbm, predlda, type=testingSet$TYPE)
predDF
#predRF predGbm predlda  type
# 1      UL      UL   nonUL nonUL
# 2   nonUL   nonUL   nonUL nonUL
# 3   nonUL   nonUL   nonUL nonUL
# 4   nonUL   nonUL   nonUL nonUL
# 5   nonUL   nonUL   nonUL nonUL
# 6   nonUL   nonUL   nonUL nonUL
# 7   nonUL   nonUL   nonUL nonUL
# 8   nonUL      UL   nonUL nonUL
# 9   nonUL   nonUL   nonUL nonUL
# 10     UL   nonUL   nonUL nonUL
# 11     UL      UL   nonUL nonUL
# 12  nonUL   nonUL   nonUL nonUL
# 13     UL      UL   nonUL nonUL
# 14  nonUL   nonUL   nonUL nonUL
# 15  nonUL   nonUL   nonUL nonUL
# 16  nonUL   nonUL   nonUL    UL
# 17     UL      UL      UL    UL
# 18     UL      UL      UL    UL
# 19     UL      UL      UL    UL
# 20     UL      UL      UL    UL
# 21     UL      UL      UL    UL
# 22  nonUL      UL      UL    UL
# 23     UL      UL      UL    UL
# 24     UL      UL      UL    UL
# 25     UL      UL      UL    UL
# 26     UL      UL   nonUL    UL
# 27     UL      UL   nonUL    UL
# 28     UL      UL      UL    UL
# 29  nonUL   nonUL   nonUL    UL
# 30  nonUL   nonUL      UL    UL
# 31  nonUL   nonUL   nonUL    UL
# 32     UL      UL   nonUL    UL
# 33  nonUL   nonUL   nonUL    UL
# 34     UL   nonUL   nonUL    UL
# 35     UL      UL   nonUL    UL

CombinedModels <- train(type~., method='gam', data=predDF)
CombinedPredictions <- predict(CombinedModels, predDF)
CombinedPredictions
# [1] UL    nonUL nonUL nonUL nonUL nonUL nonUL nonUL nonUL UL    UL    nonUL UL    nonUL nonUL
# [16] nonUL UL    UL    UL    UL    UL    UL    UL    UL    UL    UL    UL    UL    nonUL UL   
# [31] nonUL UL    nonUL UL    UL 

predRF
# [1] UL    nonUL nonUL nonUL nonUL nonUL nonUL nonUL nonUL UL    UL    nonUL UL    nonUL nonUL
# [16] nonUL UL    UL    UL    UL    UL    nonUL UL    UL    UL    UL    UL    UL    nonUL nonUL
# [31] nonUL UL    nonUL UL    UL   

predGbm
# [1] UL    nonUL nonUL nonUL nonUL nonUL nonUL UL    nonUL nonUL UL    nonUL UL    nonUL nonUL
# [16] nonUL UL    UL    UL    UL    UL    UL    UL    UL    UL    UL    UL    UL    nonUL nonUL
# [31] nonUL UL    nonUL nonUL UL   

predlda
# [1] nonUL nonUL nonUL nonUL nonUL nonUL nonUL nonUL nonUL nonUL nonUL nonUL nonUL nonUL nonUL
# [16] nonUL UL    UL    UL    UL    UL    UL    UL    UL    UL    nonUL nonUL UL    nonUL UL   
# [31] nonUL nonUL nonUL nonUL nonUL

sum <- sum(CombinedPredictions==testingSet$TYPE)#28
length <- length(CombinedPredictions)#35

sum <- sum(predRF==testingSet$TYPE) #24
length <- length(predRF)#35
accuracy_rfMOd <- (sum/length) #[1] 0.6857143

sum <- sum(predGbm==testingSet$TYPE)#25
accuracy_gbmMod <- sum/length # [1] 0.6857143


sum <- sum(predlda==testingSet$TYPE) #26
accuracy_ldaMod <- sum/length #[1] 0.7428571

rf2 <- randomForest(TYPE~., data=trainingSet, method='class')

png('randomForest2.png')
plot(rf2)
dev.off()

predRF2 <- predict(rf2, testingSet, type='class')
predRF2
# gsm1667145   gsm336254   gsm336258   gsm336260   gsm336270   gsm336273   gsm336276    gsm52662 
# UL       nonUL       nonUL       nonUL          UL       nonUL       nonUL       nonUL 
# gsm52663    gsm52665    gsm52667    gsm52669     gsm9099   gsm569425   gsm569427 gsm336202ul 
# nonUL          UL          UL       nonUL          UL       nonUL       nonUL       nonUL 
# gsm336208ul gsm336209ul gsm336214ul gsm336215ul gsm336218ul gsm336220ul gsm336229ul gsm336232ul 
# UL          UL          UL          UL          UL       nonUL          UL          UL 
# gsm336234ul gsm336238ul gsm336239ul gsm336240ul gsm336241ul gsm336245ul gsm336248ul  gsm38689ul 
# UL          UL          UL          UL       nonUL       nonUL       nonUL       nonUL 
# gsm38692ul   gsm9094ul gsm569429ul 
# nonUL          UL          UL 

sum <- sum(predRF2==testingSet$TYPE)#23
accuracy_RF2 <- sum/length #[1] 0.6571429

confusionMatrix(predRF2, testingSet$TYPE)
# Confusion Matrix and Statistics
# 
# Reference
# Prediction nonUL UL
# nonUL    10  7
# UL        5 13
# 
# Accuracy : 0.6571          
# 95% CI : (0.4779, 0.8087)
# No Information Rate : 0.5714          
# P-Value [Acc > NIR] : 0.1974          
# 
# Kappa : 0.3115          
# 
# Mcnemar's Test P-Value : 0.7728          
#                                           
#             Sensitivity : 0.6667          
#             Specificity : 0.6500          
#          Pos Pred Value : 0.5882          
#          Neg Pred Value : 0.7222          
#              Prevalence : 0.4286          
#          Detection Rate : 0.2857          
#    Detection Prevalence : 0.4857          
#       Balanced Accuracy : 0.6583          
#                                           
#        'Positive' Class : nonUL     

predDF2 <- data.frame(predRF,predRF2,predlda, predGbm, testingSet$TYPE)
colnames(predDF2)[5] <- 'TYPE'

CombinedModels <- train(TYPE ~ ., method='gam', data=predDF2)
CombinedPredictions2 <- predict(CombinedModels, predDF)
sum <- sum(CombinedPredictions2==testingSet$TYPE)#29
length(CombinedPredictions2)#35
accuracy_CombinedPredictions2 <- sum/length #[1] 0.8285714

predDF3 <- data.frame(predRF,predRF2,predlda, predGbm, 
                      CombinedPredictions2, testingSet$TYPE)#35X6
colnames(predDF3)[6] <- 'TYPE'


results <- c(round(accuracy_rfMOd,2), round(accuracy_RF2,2), 
             round(accuracy_ldaMod, 2), round(accuracy_gbmMod,2), 
             round(accuracy_CombinedPredictions2,2), round(100,2))
results <- as.factor(results)
results <- t(data.frame(results))#1X6
colnames(results) <- colnames(predDF3)
Results <- rbind(predDF3, results) #36X6

write.csv(Results, 'Results_predictions_TOP16.csv', row.names=TRUE)

#try 'rpart', 'knn', 'glm', 
knnMod <- train(TYPE ~ .,
                method='knn', preProcess=c('center','scale'),
                tuneLength=10, trControl=trainControl(method='cv'), data=trainingSet)#list 23
plot(knnMod)

rpartMod <- train(TYPE ~ ., method='rpart', tuneLength=9, data=trainingSet) #lists 23
plot(rpart)

glmMod <- train(TYPE ~ ., 
                method='glm', data=trainingSet) #list of 23, 32 warnings
#plot(glmMod)#errors

predKNN <- predict(knnMod, testingSet)
predRPART <- predict(rpartMod, testingSet)
predGLM <- predict(glmMod, testingSet)



length=length(testingSet$TYPE)#36

sumKNN <- sum(predKNN==testingSet$TYPE)#25
sumRPart <- sum(predRPART==testingSet$TYPE)#19
sumGLM <- sum(predGLM==testingSet$TYPE)#26

accuracy_KNN <- sumKNN/length #0.7143
accuracy_RPART <- sumRPart/length #0.543
accuracy_GLM <- sumGLM/length #0.743

predDF3 <- data.frame(predDF2[,1:4],predKNN,predRPART,predGLM, 
                      TYPE=testingSet$TYPE)

CombinedModels <- train(TYPE ~ ., method='gam', data=predDF3)
CombinedPredictions2 <- predict(CombinedModels, predDF3)
accuracy_CP2 <- sum(CombinedPredictions2==testingSet$TYPE)/length #0.826

predDF4 <- data.frame(predDF3[,1:7], CombinedPredictions2, TYPE=testingSet$TYPE)#35X9
colnames(predDF4)
# [1] "predRF"               "predRF2"              "predlda"              "predGbm"             
# [5] "predKNN"              "predRPART"            "predGLM"              "CombinedPredictions2"
# [9] "TYPE"

results <- c(round(accuracy_rfMOd,2), round(accuracy_RF2,2), 
             round(accuracy_ldaMod, 2), round(accuracy_gbmMod,2), 
             round(accuracy_KNN,2), round(accuracy_RPART,2),
             round(accuracy_GLM,2), 
             round(accuracy_CP2,2), round(100,2))
results <- as.factor(results)
results <- t(data.frame(results))#1X6
colnames(results) <- colnames(predDF4)
Results <- rbind(predDF4, results) #36X9

write.csv(Results, 'Results_predictions_TOP16_8_algorithms.csv', row.names=TRUE)


####################################################################################

# Test out Principal Component Analysis (PCA) on this data table of obs. X genes

samplesType <- read.csv('samplesType.csv', sep=',', row.names=1, header=TRUE,
                        na.strings=c('','NA')) #121 samples X 17 
#(16 genes and a type field for 'UL' or 'nonUL')

# When a lot of feature/field/column variation in the data, this eliminates the 
# noise better using Singular Value Decomposition (SVD) matrices of the data
# scale variable, if not done so already (these are), PCA reduces rank or features 
# to explain data best

install.packages('kernlab')
library(kernlab)
library(caret)
library(randomForest)
library(e1071)
library(MASS)
library(gbm)

set.seed(2348765)
inTrain <- createDataPartition(y=samplesType$TYPE, p=0.70, list=FALSE)
trainingSet <- samplesType[inTrain,] #85X17
testingSet <- samplesType[-inTrain,] #36X17

# Leave out the 'TYPE' outcome
M <- abs(cor(training[,-1])) #16X16
diag(M) <- 0
which(M>0.75, arr.ind=T)#squeeze the limits until 2 left
#       row col
# ASPSCR1   13   6
# CYTH4   6  13

png('plotScatter45_PCA_CYTH4_ASPSCR1.png')
plot(trainingSet[,6], trainingSet[,13]) #scatter plot not showing on the 45 degree line,
# if the plot did show the perfect or almost perfect alignment on the 45, then leave
# these two genes out for PCA
dev.off()

cyth4ASPSCR1 <- trainingSet[, c(7,14)] #add a column for the one omitted in table of cor() above
prComp <- prcomp(cyth4ASPSCR1)

plot(prComp$x[,1], prComp$x[,2]) #plots 1st principal component against 2nd 

prComp$rotation # explains most variation by adding columns
# PC1        PC2
# CYTH4 -0.8063751 -0.5914044
# ASPSCR1  -0.5914044  0.8063751

typeColor <- unique(factor(trainingSet$TYPE)) #UL or nonUL, red or black
prComp <- prcomp(trainingSet[,-1])#remove outcome field 'TYPE', add 1 for symm.
plot(prComp$x[,1], prComp$x[,2], col=typeColor, xlab='PC1', ylab='PC2')


preProc <- preProcess(trainingSet, method='pca',pcaComp=2)#list of 21
trainPC <- predict(preProc, trainingSet)#85X3

#KNN model
modelFit <- train(TYPE ~ .,
                  method='knn', preProcess=c('center','scale'),
                  tuneLength=10, trControl=trainControl(method='cv'), data=trainPC)
testPC <- predict(preProc, testingSet[,-1])#36X2
confusionMatrix(testingSet$TYPE, predict(modelFit, testPC)) #Accuracy of 56%
# Confusion Matrix and Statistics
# 
# Reference
# Prediction nonUL UL
# nonUL     2 13
# UL        3 18
# 
# Accuracy : 0.5556         
# 95% CI : (0.381, 0.7206)
# No Information Rate : 0.8611         
# P-Value [Acc > NIR] : 1.00000        
# 
# Kappa : -0.0105        
# 
# Mcnemar's Test P-Value : 0.02445        
#                                          
#             Sensitivity : 0.40000        
#             Specificity : 0.58065        
#          Pos Pred Value : 0.13333        
#          Neg Pred Value : 0.85714        
#              Prevalence : 0.13889        
#          Detection Rate : 0.05556        
#    Detection Prevalence : 0.41667        
#       Balanced Accuracy : 0.49032        
#                                          
#        'Positive' Class : nonUL          
#                                        
modelFit <- train(TYPE ~ ., 
                  method='glm', data=trainPC) #list of 23

testPC <- predict(preProc, testing[,-1])
confusionMatrix(testing$TYPE, predict(modelFit, testPC)) #produces only 65% accuracy
# Confusion Matrix and Statistics
# 
# Reference
# Prediction nonUL UL
# nonUL     3  9
# UL        1 16
# 
# Accuracy : 0.6552          
# 95% CI : (0.4567, 0.8206)
# No Information Rate : 0.8621          
# P-Value [Acc > NIR] : 0.99893         
# 
# Kappa : 0.212           
# 
# Mcnemar's Test P-Value : 0.02686 
# 
# install.packages('rpart')
library(rpart)

# same as above, but using 'rpart' method instead of 'glm':
preProc <- preProcess(training, method='pca',pcaComp=2)#list of 21
trainPC <- predict(preProc, training)#92X2

modelFit <- train(TYPE ~ ., 
                  method='rpart', tuneLength=9, data=trainPC) 

testPC <- predict(preProc, testing[,-1])
testPC

confusionMatrix(testing$TYPE, predict(modelFit, testPC)) #accuracy of 48.28%
# Confusion Matrix and Statistics
# 
# Reference
# Prediction nonUL UL
# nonUL     4  8
# UL        7 10
# 
# Accuracy : 0.4828          
# 95% CI : (0.2945, 0.6747)
# No Information Rate : 0.6207          
# P-Value [Acc > NIR] : 0.9557          
# 
# Kappa : -0.0794         
# 
# Mcnemar's Test P-Value : 1.0000          
#                                           
#             Sensitivity : 0.3636          
#             Specificity : 0.5556          
#          Pos Pred Value : 0.3333          
#          Neg Pred Value : 0.5882          
#              Prevalence : 0.3793          
#          Detection Rate : 0.1379          
#    Detection Prevalence : 0.4138          
#       Balanced Accuracy : 0.4596          
#                                           
#        'Positive' Class : nonUL  
# 
# using lm, will run previous results, and 'wrong model type for classification':
preProc <- preProcess(training, method='pca',pcaComp=2)#list of 21
trainPC <- predict(preProc, training)#92X2

modelFit <- train(TYPE ~ ., 
                  method='lm',data=trainPC) 

testPC <- predict(preProc, testing[,-1])
testPC

confusionMatrix(testing$TYPE, predict(modelFit, testPC))




# results of prediction using PCA and 'glm' or 'rpart' on the sample prediction by
# gene values was 65% and 48% respectively. The previous predictions not using PCA
# had 71%, 71%, 74%, and 68.57% for rf, gbm, and lda from the methods of rpart pkg, 
# and the 68.57% was from the randomForest pkg using the randomForest()
####################################################################################
#########################################################################################
#########################################################################################
# What if the six genes were left out and only the top 16 genes with highest magnitude
# of chanage in general were used
# to build the lda model on? Would the accuracy be better?

# Let us find out if the top 16 magnitude or differentially expressed genes of the 130
# in common and on the chromosomes 11,12,17, and 22 that the six genes ubiquitous to 
# current literature are better predictors of UL

# go to the 'MemberMagnitude_130_142.csv' file to get the same data before filtering
# for the top 10 and adding the six genes
All <- read.csv('MemberMagnitude_130_142.csv', sep=',', header=TRUE,
                na.string=c('','NA')) # 130 X 142

# The magnitude of gene expression change up or down between UL and nonUL samples is 
# already sorted, choose the first 16 genes listed here as well as all the meta info
# per gene to work with, without using simulated means or changes from mean of UL
# to nonUL

magnitude <- All[1:16,c(1,21,22:142)] #16X123, samples, gene, and magnitude fields

row.names(magnitude) <- magnitude$genes

all_16 <- magnitude[,-c(1:2)] #16X121, samples only, genes as row names

all_16_t <- data.frame(t(all_16)) #121 samples by 16 genes having most change in UL

write.csv(all_16_t, 'top_16_overall.csv', row.names=TRUE)

# use those previous algorithms on this set of 16 most differentially expressed genes
# the caret lda and rf algorithms

# add in the TYPE field for UL and nonUL
row.names(all_16_t) # 1:51 are nonUL, 52:121 are UL

ul <- data.frame(rep('UL', 70)) #70
non <- data.frame(rep('nonUL', 51)) #51
colnames(ul) <- 'TYPE'
colnames(non) <- 'TYPE'

TYPE <- rbind(non,ul)#121X1

all_16 <- cbind(TYPE, all_16_t) #121X16
write.csv(all_16, 'most_DE_ml_ready_130.csv', row.names=TRUE)
colnames(all_16)
# [1] "TYPE"     "FSCN2"    "CBX2"   "SOCS3"    "CANT1"    "IRF7"     "ARHGDIA"  "NOL12"   
# [9] "SLC25A10" "SLC38A10" "ASPSCR1"     "SYNGR1"   "TALDO1"   "FN3K"     "POLR2F"   "LLPH"    
# [17] "SMCR7L"

library(caret)
set.seed(189678345)

inTrain <- createDataPartition(all_16$TYPE, p=0.7, list=FALSE)

trainingSet <- all_16[inTrain,]#85X17
testingSet <- all_16[-inTrain,]#36X17


# latent dirichlet allocation model, best model with same results of 74%
ldaMod <- train(TYPE~., method='lda', data=trainingSet)

# run predictions on the testing set for the lda model
predlda <- predict(ldaMod, testingSet)


sum <- sum(predlda==testingSet$TYPE) #28
accuracy_lda <- sum/length(testingSet$TYPE) # [1] 0.7777778


# next best model with rf for random forest method of caret:
rfMod <- train(TYPE ~ ., method='rf', data=trainingSet)

predRF <- predict(rfMod, testingSet)

sum <- sum(predRF==testingSet$TYPE) #25
accuracy_rfMod <- sum/length(testingSet$TYPE) # [1] 0.694

gbmMod <- train(TYPE~., method='gbm', data=trainingSet, verbose=FALSE )
predGbm <- predict(gbmMod, testingSet)
sum <- sum(predGbm==testingSet$TYPE)#24
accuracy_gbmMod <- sum/length(testingSet$TYPE)#0.6666

predDF <- data.frame(predRF, predGbm, predlda, type=testingSet$TYPE)#36X4

CombinedModels <- train(type~., method='gam', data=predDF)
CombinedPredictions <- predict(CombinedModels, predDF)
CombinedPredictions

accuracy_CombinedPredictions <- sum(CombinedPredictions==testingSet$TYPE)/
  length(testingSet$TYPE)#0.77777777

rf2 <- randomForest(TYPE~., data=trainingSet, method='class')

plot(rf2)

predRF2 <- predict(rf2, testingSet, type='class')

sum <- sum(predRF2==testingSet$TYPE)#24
accuracy_RF2 <- sum/length(testingSet$TYPE) # 0.66666667

confusionMatrix(predRF2, testingSet$TYPE)
predDF2 <- data.frame(predRF, predRF2, predlda, predGbm, TYPE=testingSet$TYPE)#36X5

CombinedModels2 <- train(TYPE~., method='gam', data=predDF2)
CombinedPredictions2 <- predict(CombinedModels2, predDF2)
CombinedPredictions2

sum <- sum(CombinedPredictions2==testingSet$TYPE)#29
accuracy_CombinedPredictions2 <- sum/length(testingSet$TYPE) #[1] 0.8055556

predDF3 <- data.frame(predRF, predRF2, predlda, predGbm, 
                      CombinedPredictions2, TYPE=testingSet$TYPE)

results <- c(round(accuracy_rfMod,2), round(accuracy_RF2,2), 
             round(accuracy_lda, 2), round(accuracy_gbmMod,2), 
             round(accuracy_CombinedPredictions2,2), round(100,2))

results <- as.factor(results)
results <- t(data.frame(results))
colnames(results) <- colnames(predDF3)
Results <- rbind(predDF3, results) #37X6


write.csv(Results,'Results_predictions_DE16.csv', row.names=TRUE)#37X6, 1.7kb

#try 'rpart', 'knn', 'glm', 
knnMod <- train(TYPE ~ .,
                method='knn', preProcess=c('center','scale'),
                tuneLength=10, trControl=trainControl(method='cv'), data=trainingSet)#list 23
plot(knnMod)

rpartMod <- train(TYPE ~ ., method='rpart', tuneLength=9, data=trainingSet) #lists 23
plot(rpart)

glmMod <- train(TYPE ~ ., 
                method='glm', data=trainingSet) #list of 23, 32 warnings
#plot(glmMod)#errors

predKNN <- predict(knnMod, testingSet)
predRPART <- predict(rpartMod, testingSet)
predGLM <- predict(glmMod, testingSet)



length=length(testingSet$TYPE)#36

sumKNN <- sum(predKNN==testingSet$TYPE)#29
sumRPart <- sum(predRPART==testingSet$TYPE)#22
sumGLM <- sum(predGLM==testingSet$TYPE)#26

accuracy_KNN <- sumKNN/length #0.9167
accuracy_RPART <- sumRPart/length #0.6111
accuracy_GLM <- sumGLM/length #0.7222

predDF3 <- data.frame(predDF2[,1:4],predKNN,predRPART,predGLM, 
                      TYPE=testingSet$TYPE)

CombinedModels <- train(TYPE ~ ., method='gam', data=predDF3)
CombinedPredictions2 <- predict(CombinedModels, predDF3)
accuracy_CP2 <- sum(CombinedPredictions2==testingSet$TYPE)/length #0.9167

predDF4 <- data.frame(predDF3[,1:7], CombinedPredictions2, TYPE=testingSet$TYPE)#35X9
colnames(predDF4)
# [1] "predRF"               "predRF2"              "predlda"              "predGbm"             
# [5] "predKNN"              "predRPART"            "predGLM"              "CombinedPredictions2"
# [9] "TYPE"

results <- c(round(accuracy_rfMod,2), round(accuracy_RF2,2), 
             round(accuracy_lda, 2), round(accuracy_gbmMod,2), 
             round(accuracy_KNN,2), round(accuracy_RPART,2),
             round(accuracy_GLM,2), 
             round(accuracy_CP2,2), round(100,2))
results <- as.factor(results)
results <- t(data.frame(results))#1X6
colnames(results) <- colnames(predDF4)
Results <- rbind(predDF4, results) #37X9

write.csv(Results,'Results_predictions_DE16_8_algorithms_used.csv', 
          row.names=TRUE)#37X9, 1.7kb
#####################################################################################
# what about on the 16 genes with the least amount of change to predict UL, is it the
# same?

All <- read.csv('MemberMagnitude_130_142.csv', sep=',', header=TRUE,
                na.string=c('','NA')) # 130 X 142

# The magnitude of gene expression change up or down between UL and nonUL samples is 
# already sorted, choose the first 16 genes listed here as well as all the meta info
# per gene to work with, without using simulated means or changes from mean of UL
# to nonUL

magnitude <- All[115:130,c(1,21,22:142)] #16X123, samples, gene, and magnitude fields

row.names(magnitude) <- magnitude$genes

all_16 <- magnitude[,-c(1:2)] #16X121, samples only, genes as row names

all_16_t <- data.frame(t(all_16)) #121 samples by 16 genes having most change in UL

row.names(all_16_t) # 1:51 are nonUL, 52:121 are UL

ul <- data.frame(rep('UL', 70)) #70
non <- data.frame(rep('nonUL', 51)) #51
colnames(ul) <- 'TYPE'
colnames(non) <- 'TYPE'

TYPE <- rbind(non,ul)#121X1

all_16 <- cbind(TYPE, all_16_t) #121X17
colnames(all_16)
# [1] "TYPE"     "PSMD13"   "DNAL4"    "ATHL1"    "CDHR5"    "SIGIRR"   "MRPL23"   "FOXK2"   
# [9] "RASSF7"   "SIRT3"    "GRAP2"    "AATK"     "TH"       "TMEM184B" "WDR45L"   "KDELR3"  
# [17] "GCGR"
write.csv(all_16, 'least_DE16_ml_ready_130.csv', row.names=TRUE)
library(caret)
set.seed(189678345)

inTrain <- createDataPartition(all_16$TYPE, p=0.7, list=FALSE)

trainingSet <- all_16[inTrain,]#85X17
testingSet <- all_16[-inTrain,]#36X17

# latent dirichlet allocation model, best model with same results of 74%
ldaMod <- train(TYPE~., method='lda', data=trainingSet)

# run predictions on the testing set for the lda model
predlda <- predict(ldaMod, testingSet)


sum <- sum(predlda==testingSet$TYPE) #15
accuracy_lda <- sum/length(testingSet$TYPE) # [1] 0.41667


# next best model with rf for random forest method of caret:
rfMod <- train(TYPE ~ ., method='rf', data=trainingSet)

predRF <- predict(rfMod, testingSet)

sum <- sum(predRF==testingSet$TYPE) #17
accuracy_rfMod <- sum/length(testingSet$TYPE) # [1] 0.4722

gbmMod <- train(TYPE~., method='gbm', data=trainingSet, verbose=FALSE )
predGbm <- predict(gbmMod, testingSet)
sum <- sum(predGbm==testingSet$TYPE)#17
accuracy_gbmMod <- sum/length(testingSet$TYPE)#0.4722

predDF <- data.frame(predRF, predGbm, predlda, TYPE=testingSet$TYPE)#36X4

CombinedModels <- train(TYPE~., method='gam', data=predDF)
CombinedPredictions <- predict(CombinedModels, predDF)
CombinedPredictions

accuracy_CombinedPredictions <- sum(CombinedPredictions==testingSet$TYPE)/
  length(testingSet$TYPE)#0.5833

rf2 <- randomForest(TYPE~., data=trainingSet, method='class')

plot(rf2)

predRF2 <- predict(rf2, testingSet, type='class')

sum <- sum(predRF2==testingSet$TYPE)#19
accuracy_RF2 <- sum/length(testingSet$TYPE) # 0.5277

confusionMatrix(predRF2, testingSet$TYPE)
predDF2 <- data.frame(predRF, predRF2, predlda, predGbm, TYPE=testingSet$TYPE)#36X5

CombinedModels2 <- train(TYPE~., method='gam', data=predDF2)
CombinedPredictions2 <- predict(CombinedModels2, predDF2)
CombinedPredictions2

sum <- sum(CombinedPredictions2==testingSet$TYPE)#23
accuracy_CombinedPredictions2 <- sum/length(testingSet$TYPE) #[1] 0.6388

predDF3 <- data.frame(predRF, predRF2, predlda, predGbm, 
                      CombinedPredictions2, TYPE=testingSet$TYPE)

results <- c(round(accuracy_rfMod,2), round(accuracy_RF2,2), 
             round(accuracy_lda, 2), round(accuracy_gbmMod,2), 
             round(accuracy_CombinedPredictions2,2), round(100,2))

results <- as.factor(results)
results <- t(data.frame(results))
colnames(results) <- colnames(predDF3)
Results <- rbind(predDF3, results) #37X6


write.csv(Results,'Results_predictions_Least_DE16.csv', row.names=TRUE)#37X6, 1.7kb

#try 'rpart', 'knn', 'glm', 
knnMod <- train(TYPE ~ .,
                method='knn', preProcess=c('center','scale'),
                tuneLength=10, trControl=trainControl(method='cv'), data=trainingSet)#list 23
plot(knnMod)

rpartMod <- train(TYPE ~ ., method='rpart', tuneLength=9, data=trainingSet) #lists 23
plot(rpart)

glmMod <- train(TYPE ~ ., 
                method='glm', data=trainingSet) #list of 23, 32 warnings
#plot(glmMod)#errors

predKNN <- predict(knnMod, testingSet)
predRPART <- predict(rpartMod, testingSet)
predGLM <- predict(glmMod, testingSet)



length=length(testingSet$TYPE)#36

sumKNN <- sum(predKNN==testingSet$TYPE)#18
sumRPart <- sum(predRPART==testingSet$TYPE)#21
sumGLM <- sum(predGLM==testingSet$TYPE)#12

accuracy_KNN <- sumKNN/length #0.5
accuracy_RPART <- sumRPart/length #0.5833
accuracy_GLM <- sumGLM/length #0.3333

predDF3 <- data.frame(predDF2[,1:4],predKNN,predRPART,predGLM, 
                      TYPE=testingSet$TYPE)

CombinedModels <- train(TYPE ~ ., method='gam', data=predDF3)
CombinedPredictions2 <- predict(CombinedModels, predDF3)
accuracy_CP2 <- sum(CombinedPredictions2==testingSet$TYPE)/length #0.7222

predDF4 <- data.frame(predDF3[,1:7], CombinedPredictions2, TYPE=testingSet$TYPE)#35X9
colnames(predDF4)
# [1] "predRF"               "predRF2"              "predlda"              "predGbm"             
# [5] "predKNN"              "predRPART"            "predGLM"              "CombinedPredictions2"
# [9] "TYPE"

results <- c(round(accuracy_rfMod,2), round(accuracy_RF2,2), 
             round(accuracy_lda, 2), round(accuracy_gbmMod,2), 
             round(accuracy_KNN,2), round(accuracy_RPART,2),
             round(accuracy_GLM,2), 
             round(accuracy_CP2,2), round(100,2))
results <- as.factor(results)
results <- t(data.frame(results))#1X6
colnames(results) <- colnames(predDF4)
Results <- rbind(predDF4, results) #37X9

write.csv(Results,'Results_predictions_Least_DE16_8_algorithms_used.csv', 
          row.names=TRUE)#37X9, 1.7kb

###########################################################################################
###########################################################################################

# Using a different data set of top fold change, most genes are the same genes, but the
# accuracy in prediction using the top 10 genes with the highest fold change and
# the 6 ubiquitous genes in the ratio
# of UL mean to non-UL mean is worse than the top 10 highest magnitude of change and 6 ubiq
# genes

# Use fold change on the 130 genes to see the top 10 and six ubiq genes have a better
# prediction accuracy

member <- read.csv('MemberMagnitude_130_142.csv', sep=',', header=TRUE, na.strings=c('','NA'))

library(dplyr)

fold <- mutate(member, foldChange= UL_Mean/nonUL_Mean)#130X143

colnames(fold)
# [1] "genes"                         "chromosome"                   
# [3] "type"                          "all"                          
# [5] "up"                            "down"                         
# [7] "majority"                      "start"                        
# [9] "end"                           "width"                        
# [11] "strand"                        "gene"                         
# [13] "transcript"                    "GENE"                         
# [15] "GENE_NAME"                     "CYTOBAND"                     
# [17] "DESCRIPTION"                   "nonUL_Mean"                   
# [19] "UL_Mean"                       "Difference_UL_minus_non_means"
# [21] "Magnitude"                     "gsm1667144"                   
# [23] "gsm1667145"                    "gsm1667146"                   
# [25] "gsm336252"                     "gsm336253"                    
# [27] "gsm336254"                     "gsm336255"                    
# [29] "gsm336256"                     "gsm336257"                    
# ...                 
# [127] "gsm38690ul"                    "gsm38691ul"                   
# [129] "gsm38692ul"                    "gsm38693ul"                   
# [131] "gsm38694ul"                    "gsm38695ul"                   
# [133] "gsm9093ul"                     "gsm9094ul"                    
# [135] "gsm9095ul"                     "gsm9096ul"                    
# [137] "gsm9097ul"                     "gsm569429ul"                  
# [139] "gsm569430ul"                   "gsm569431ul"                  
# [141] "gsm569432ul"                   "gsm569433ul"                  
# [143] "foldChange"  

# put the 'foldChange' field next to the magnitude field
fold <- fold[,c(1:21,143,22:142)]

head(colnames(fold),24)
# [1] "genes"                         "chromosome"                   
# [3] "type"                          "all"                          
# [5] "up"                            "down"                         
# [7] "majority"                      "start"                        
# [9] "end"                           "width"                        
# [11] "strand"                        "gene"                         
# [13] "transcript"                    "GENE"                         
# [15] "GENE_NAME"                     "CYTOBAND"                     
# [17] "DESCRIPTION"                   "nonUL_Mean"                   
# [19] "UL_Mean"                       "Difference_UL_minus_non_means"
# [21] "Magnitude"                     "foldChange"                   
# [23] "gsm1667144"                    "gsm1667145"


write.csv(fold, 'fold_magnitude_member_gviz_130_143.csv', row.names=FALSE) #114.3 kb

# order the foldChange most to least
FOLD <- fold[order(fold$foldChange, decreasing=TRUE),] 

# get top 10 from FOLD by highest fold change in UL/nonUL samples
FOLD10 <- FOLD[1:10,]

FOLD10$genes
#[1] CANT1    ARHGDIA  SLC25A10 ATF4     CBX2   MRPL12   BET1L    ASPSCR1     NOL12    IRF7    

#BET1l is in the top 10 genes having the highest fold change
# but none of the other six ubiquitous genes are in the top 10 of highest magnitude of change
# or highest fold change
# since, BET1L is already included, get the first 11, and grep() the other 5 ubuiquitous genes

t <- grep('TNRC6B', FOLD$genes)#62
cy <- grep('CYTH4', FOLD$genes)#29
cc <- grep('CCDC57', FOLD$genes)#106
F <- grep('FASN', FOLD$genes)#90
H <- grep('HMGA2', FOLD$genes)#128

FOLD5 <- FOLD[c(t,cy,cc,F,H),]#5X143

FOLD11 <- FOLD[1:11,] #11X143

TOP16_fold <- rbind(FOLD11, FOLD5) #16X143
TOP16_fold$genes
# [1] CANT1    ARHGDIA  SLC25A10 ATF4     CBX2   MRPL12   BET1L    ASPSCR1     NOL12    IRF7    
# [11] DRD4     TNRC6B   CYTH4    CCDC57   FASN     HMGA2  

#This is a table of ml ready TOP16 by DE, made in other script
samplesType <- read.csv('samplesType.csv', sep=',', header=TRUE, row.names=1) 

#make a ml ready table of the TOP16_fold data set
row.names(TOP16_fold) <- TOP16_fold$genes

TOP16_fold_ml <- TOP16_fold[,c(23:143)] #16X121

TOP16_fold_ml_t <- data.frame(t(TOP16_fold_ml)) #121X16

#add the 'TYPE' field of 'UL' or 'nonUL' to the foldchange top 16 genes
ul <- data.frame(rep('UL',70)) #70X1
non <- data.frame(rep('nonUL', 51))#51X1
colnames(ul) <- 'TYPE'
colnames(non) <- 'TYPE'

TYPE <- rbind(non,ul)

FOLD16_ml_ready <- cbind(TYPE,TOP16_fold_ml_t) #121X17

#save FOLD16_ml_ready.csv and TOP16_ml_ready.csv for the DE selections
write.csv(FOLD16_ml_ready, 'FOLD16_ml_ready.csv', row.names=TRUE)
write.csv(samplesType, 'TOP16_ml_ready.csv', row.names=TRUE)

#######################################################################

# ML on the top 10 genes of the 130 with highest fold change and the 6 ubiquitous genes

fold <- read.csv('FOLD16_ml_ready.csv', sep=',', header=TRUE, row.names=1)
colnames(fold)

library(caret)
library(randomForest)
library(MASS)
library(gbm)

# set the same seed used earlier in ml on TOP16 genes

set.seed(189678345) # this will reproduce the same numbers each time sampling is done
# run set.seed value to get exact results, otherwise results differ

# creates a partition of the data by indices
inTrain <- createDataPartition(y=fold$TYPE, p=0.7, list=FALSE)
# length of 85 out of 121 samples sampled

trainingSet <- fold[inTrain,]#85X17
testingSet <- fold[-inTrain,]#36X17

#install.packages('e1071')
library(e1071)

# randomForest, cross-validation (cv) = 5
rfMod <- train(TYPE~., method='rf', data=(trainingSet), 
               trControl=trainControl(method='cv'), number=5)
plot(rfMod)

# generalizedBoostedModel
gbmMod <- train(TYPE~., method='gbm', data=trainingSet, verbose=FALSE )

plot(gbmMod)


# linkage dirichlet allocation model
ldaMod <- train(TYPE~., method='lda', data=trainingSet)

# the ldaMod cannot be plotted, 'no tuning parameters for this model'

# run predictions on the testing set
predRF <- predict(rfMod, testingSet)
predGbm <- predict(gbmMod, testingSet)
predlda <- predict(ldaMod, testingSet)

predDF <- data.frame(predRF, predGbm, predlda, type=testingSet$TYPE)
predDF

CombinedModels <- train(type~., method='gam', data=predDF)
CombinedPredictions <- predict(CombinedModels, predDF)
CombinedPredictions

sum <- sum(CombinedPredictions==testingSet$TYPE)#27
length <- length(testingSet$TYPE)#36
accuracy_CP1 <- sum/length #0.75

sum <- sum(predRF==testingSet$TYPE) #23
length <- length(testingSet$TYPE)#36
accuracy_rfMod <- (sum/length) #0.6389

sum <- sum(predGbm==testingSet$TYPE)#26
accuracy_gbmMod <- sum/length # 0.7222

sum <- sum(predlda==testingSet$TYPE) #24
accuracy_ldaMod <- sum/length #0.6667

rf2 <- randomForest(TYPE~., data=trainingSet, method='class')

plot(rf2)

predRF2 <- predict(rf2, testingSet, type='class')

sum <- sum(predRF2==testingSet$TYPE)#25
accuracy_RF2 <- sum/length # 0.6944

confusionMatrix(predRF2, testingSet$TYPE)
predDF2 <- data.frame(predRF, predRF2, predlda, predGbm, TYPE=testingSet$TYPE)

CombinedModels2 <- train(TYPE~., method='gam', data=predDF2)
CombinedPredictions2 <- predict(CombinedModels2, predDF2)
CombinedPredictions2

sum <- sum(CombinedPredictions2==testingSet$TYPE)
accuracy_CombinedPredictions2 <- sum/length

predDF3 <- data.frame(predRF, predRF2, predlda, predGbm, 
                      CombinedPredictions2, TYPE=testingSet$TYPE)

results <- c(round(accuracy_rfMOd,2), round(accuracy_RF2,2), 
             round(accuracy_ldaMod, 2), round(accuracy_gbmMod,2), 
             round(accuracy_CombinedPredictions2,2), round(100,2))

results <- as.factor(results)
results <- t(data.frame(results))
colnames(results) <- colnames(predDF3)
Results <- rbind(predDF3, results) #37X6


write.csv(Results,'Results_predictions_FOLD16.csv', row.names=TRUE)#37X6, 1.7kb

#try 'rpart', 'knn', 'glm', 
knnMod <- train(TYPE ~ .,
                method='knn', preProcess=c('center','scale'),
                tuneLength=10, trControl=trainControl(method='cv'), data=trainingSet)#list 23
plot(knnMod)

rpartMod <- train(TYPE ~ ., method='rpart', tuneLength=9, data=trainingSet) #lists 23
plot(rpart)

glmMod <- train(TYPE ~ ., 
                method='glm', data=trainingSet) #list of 23, 32 warnings
#plot(glmMod)#errors

predKNN <- predict(knnMod, testingSet)
predRPART <- predict(rpartMod, testingSet)
predGLM <- predict(glmMod, testingSet)



length=length(testingSet$TYPE)#36

sumKNN <- sum(predKNN==testingSet$TYPE)#25
sumRPart <- sum(predRPART==testingSet$TYPE)#21
sumGLM <- sum(predGLM==testingSet$TYPE)#26

accuracy_KNN <- sumKNN/length #0.694
accuracy_RPART <- sumRPart/length #0.5833
accuracy_GLM <- sumGLM/length #0.7222

predDF3 <- data.frame(predDF2[,1:4],predKNN,predRPART,predGLM, 
                      TYPE=testingSet$TYPE)

CombinedModels <- train(TYPE ~ ., method='gam', data=predDF3)
CombinedPredictions2 <- predict(CombinedModels, predDF3)
accuracy_CP2 <- sum(CombinedPredictions2==testingSet$TYPE)/length #0.7778

predDF4 <- data.frame(predDF3[,1:7], CombinedPredictions2, TYPE=testingSet$TYPE)#36X9
colnames(predDF4)
# [1] "predRF"               "predRF2"              "predlda"              "predGbm"             
# [5] "predKNN"              "predRPART"            "predGLM"              "CombinedPredictions2"
# [9] "TYPE"

results <- c(round(accuracy_rfMod,2), round(accuracy_RF2,2), 
             round(accuracy_ldaMod, 2), round(accuracy_gbmMod,2), 
             round(accuracy_KNN,2), round(accuracy_RPART,2),
             round(accuracy_GLM,2), 
             round(accuracy_CP2,2), round(100,2))
results <- as.factor(results)
results <- t(data.frame(results))#1X6
colnames(results) <- colnames(predDF4)
Results <- rbind(predDF4, results) #37X9

write.csv(Results,'Results_predictions_FOLD16_8_algorithms_used.csv', 
          row.names=TRUE)#37X9, 2.7kb

########################################################################################

#### What about those genes part of the majority of genes that are up/down expressed? 
#### Are those good indicators of UL or nonUL?

member <- read.csv('MemberMagnitude_130_142.csv', sep=',', 
                   header=TRUE,na.strings=c('','NA'))
majority <- subset(member, member$majority=='TRUE')#77X142
down <- subset(majority, type=='down') #23
up <- subset(majority, type=='up') #54X142

#these are ordered by magnitude alread, select top 5 from both then combine

down <- down[1:5,]
up <- up[1:5,]

majority10 <- rbind(down,up)
row.names(majority10) <- majority10$genes

maj10 <- data.frame(t(majority10[,-c(1:21)]))#121X10, removes meta fields

UL <- t(data.frame(rep('UL', 70)))
nonUL <- t(data.frame(rep('nonUL', 51)))
TYPE <- cbind(nonUL,UL)
TYPE <- data.frame(t(TYPE))
colnames(TYPE) <- 'TYPE'
MAJ10 <- cbind(TYPE, maj10)#121X11
colnames(MAJ10)

write.csv(MAJ10, 'majority_ml_ready_10_total.csv')

colnames(MAJ10)
# "SOCS3"    "TALDO1"   "ASCL2"    "HMGA2"    "PNPLA2"   "CBX2"   "CANT1"   
# "ARHGDIA"  "NOL12"    "SLC25A10"


#Now the analytics can be done, using the same packages
set.seed(189678345) # this will reproduce the same numbers each time sampling is done
# run set.seed value to get exact results, otherwise results differ

# creates a partition of the data by indices
inTrain <- createDataPartition(y=MAJ10$TYPE, p=0.7, list=FALSE)#85

trainingSet <- MAJ10[inTrain,]#85X7
testingSet <- MAJ10[-inTrain,]#36X7

# randomForest, cross-validation (cv) = 5
rfMod <- train(TYPE~., method='rf', data=(trainingSet), 
               trControl=trainControl(method='cv'), number=5)
plot(rfMod)

# generalizedBoostedModel
gbmMod <- train(TYPE~., method='gbm', data=trainingSet, verbose=FALSE )

plot(gbmMod)


# linkage dirichlet allocation model
ldaMod <- train(TYPE~., method='lda', data=trainingSet)

# the ldaMod cannot be plotted, 'no tuning parameters for this model'

# run predictions on the testing set
predRF <- predict(rfMod, testingSet)
predGbm <- predict(gbmMod, testingSet)
predlda <- predict(ldaMod, testingSet)

predDF <- data.frame(predRF, predGbm, predlda, type=testingSet$TYPE)
predDF

CombinedModels <- train(type~., method='gam', data=predDF)
CombinedPredictions <- predict(CombinedModels, predDF)
CombinedPredictions

sum <- sum(CombinedPredictions==testingSet$TYPE)#30
length <- length(testingSet$TYPE)#36
accuracy_CP1 <- sum/length #0.833

sum <- sum(predRF==testingSet$TYPE) #26
length <- length(testingSet$TYPE)#36
accuracy_rfMod <- (sum/length) #0.722

sum <- sum(predGbm==testingSet$TYPE)#21
accuracy_gbmMod <- sum/length # 0.583

sum <- sum(predlda==testingSet$TYPE) #29
accuracy_ldaMod <- sum/length #0.805

rf2 <- randomForest(TYPE~., data=trainingSet, method='class')

plot(rf2)

predRF2 <- predict(rf2, testingSet, type='class')

sum <- sum(predRF2==testingSet$TYPE)#25
accuracy_RF2 <- sum/length # 0.6944

confusionMatrix(predRF2, testingSet$TYPE)
predDF2 <- data.frame(predRF, predRF2, predlda, predGbm, TYPE=testingSet$TYPE)

CombinedModels2 <- train(TYPE~., method='gam', data=predDF2)
CombinedPredictions2 <- predict(CombinedModels2, predDF2)
CombinedPredictions2

sum <- sum(CombinedPredictions2==testingSet$TYPE)#30
accuracy_CombinedPredictions2 <- sum/length #0.8333

predDF3 <- data.frame(predRF, predRF2, predlda, predGbm, 
                      CombinedPredictions2, TYPE=testingSet$TYPE)

results <- c(round(accuracy_rfMod,2), round(accuracy_RF2,2), 
             round(accuracy_ldaMod, 2), round(accuracy_gbmMod,2), 
             round(accuracy_CombinedPredictions2,2), round(100,2))

results <- as.factor(results)
results <- t(data.frame(results))
colnames(results) <- colnames(predDF3)
Results <- rbind(predDF3, results) #37X6


write.csv(Results,'Results_predictions_majority10.csv', row.names=TRUE)#37X6, 1.7kb

#try 'rpart', 'knn', 'glm', 
knnMod <- train(TYPE ~ .,
                method='knn', preProcess=c('center','scale'),
                tuneLength=10, trControl=trainControl(method='cv'), data=trainingSet)#list 23
plot(knnMod)

rpartMod <- train(TYPE ~ ., method='rpart', tuneLength=9, data=trainingSet) #lists 23
plot(rpart)

glmMod <- train(TYPE ~ ., 
                method='glm', data=trainingSet) #list of 23, 32 warnings
#plot(glmMod)#errors

predKNN <- predict(knnMod, testingSet)
predRPART <- predict(rpartMod, testingSet)
predGLM <- predict(glmMod, testingSet)



length=length(testingSet$TYPE)#36

sumKNN <- sum(predKNN==testingSet$TYPE)#29
sumRPart <- sum(predRPART==testingSet$TYPE)#24
sumGLM <- sum(predGLM==testingSet$TYPE)#29

accuracy_KNN <- sumKNN/length #0.8055
accuracy_RPART <- sumRPart/length #0.6667
accuracy_GLM <- sumGLM/length #0.8055

predDF3 <- data.frame(predDF2[,1:4],predKNN,predRPART,predGLM, 
                      TYPE=testingSet$TYPE)

CombinedModels <- train(TYPE ~ ., method='gam', data=predDF3)
CombinedPredictions2 <- predict(CombinedModels, predDF3)
accuracy_CP2 <- sum(CombinedPredictions2==testingSet$TYPE)/length #0.8333

predDF4 <- data.frame(predDF3[,1:7], CombinedPredictions2, TYPE=testingSet$TYPE)#36X9
colnames(predDF4)
# [1] "predRF"               "predRF2"              "predlda"              "predGbm"             
# [5] "predKNN"              "predRPART"            "predGLM"              "CombinedPredictions2"
# [9] "TYPE"

results <- c(round(accuracy_rfMod,2), round(accuracy_RF2,2), 
             round(accuracy_ldaMod, 2), round(accuracy_gbmMod,2), 
             round(accuracy_KNN,2), round(accuracy_RPART,2),
             round(accuracy_GLM,2), 
             round(accuracy_CP2,2), round(100,2))
results <- as.factor(results)
results <- t(data.frame(results))#1X6
colnames(results) <- colnames(predDF4)
Results <- rbind(predDF4, results) #37X9

write.csv(Results,'Results_predictions_majority10_8_algorithms_used.csv', 
          row.names=TRUE)#37X9, 2.7kb

########################################################################################
# What about the 12,173 genes all in common among the five series before filtering for 
# the four cytoband locations, what genes are the top genes in fold change being the 
# greatest in UL compared to non-UL?

universeFold <- read.csv('orderedCommonMagnitude12173X129.csv', sep=',', header=TRUE, 
                         na.strings=c('','NA'), row.names=1)#12173X129

library(dplyr)

universe <- mutate(universeFold, foldChange=UL_mean/nonUL_mean) #12173X130
universe <- universe[,c(1:8,130,9:129)]

write.csv(universe, 'all_common_12173_130_fold_magnitude.csv', row.names=TRUE)

universeOrdered <- universe[order(universe$foldChange, decreasing=TRUE),]
row.names(universeOrdered) <- universe$GENE_SYMBOL

# according to fold change rank following lineage from most fold change as 1, to least as 12,173:
t <- grep('TNRC6B', universeOrdered$GENE_SYMBOL)#3841
b <- grep('BET1L', universeOrdered$GENE_SYMBOL)#2100
cy <- grep('CYTH4', universeOrdered$GENE_SYMBOL)#8848
cc <- grep('CCDC57', universeOrdered$GENE_SYMBOL)#2585
h <- grep('HMGA2', universeOrdered$GENE_SYMBOL)#253
F <- grep('FASN', universeOrdered$GENE_SYMBOL)#2351

universe16 <- universeOrdered[1:16,-c(1:9)]

universe_t <- data.frame(t(universe16)) # 121 X 16 genes

row.names(universe_t) #they are ordered with the first 51 non-UL and next 70 UL samples

# create the 'TYPE' field to identify each sample as UL or nonUL
ul <- data.frame(TYPE=rep('UL', 70))
non <- data.frame(TYPE=rep('nonUL',51))
TYPE <- rbind(non,ul) #121X1

# add the TYPE field to the universal 16 genes having the most fold change
UNIVERSE <- cbind(TYPE, universe_t) #121X17

write.csv(UNIVERSE, 'most_universe_fold.csv', row.names=TRUE)

colnames(UNIVERSE)
# [1] "TYPE"    "HSPB1"   "DSTN"    "S100A6"  "CNN1"    "ACTG2"   "VIM"     "SPARCL1" "TPM2"    "ACTA2"  
#[11] "PCP4"    "TAGLN"   "DES"     "RAMP1"   "CYR61"   "UBC"     "ACTB"

#none of those genes are any of the six ubiquitous genes to UL risk studies nor are they
#in the TOP16

# Now see if these top fold change genes (16) are good predictors for UL or nonUL 
# using the same algorithms

set.seed(189678345) # this will reproduce the same numbers each time sampling is done
# run set.seed value to get exact results, otherwise results differ

# creates a partition of the data by indices
inTrain <- createDataPartition(y=UNIVERSE$TYPE, p=0.7, list=FALSE)#85

trainingSet <- UNIVERSE[inTrain,]#85X7
testingSet <- UNIVERSE[-inTrain,]#36X7

# randomForest, cross-validation (cv) = 5
rfMod <- train(TYPE~., method='rf', data=(trainingSet), 
               trControl=trainControl(method='cv'), number=5)
plot(rfMod)

# generalizedBoostedModel
gbmMod <- train(TYPE~., method='gbm', data=trainingSet, verbose=FALSE )

plot(gbmMod)


# latent dirichlet allocation model
ldaMod <- train(TYPE~., method='lda', data=trainingSet)

# the ldaMod cannot be plotted, 'no tuning parameters for this model'

# run predictions on the testing set
predRF <- predict(rfMod, testingSet)
predGbm <- predict(gbmMod, testingSet)
predlda <- predict(ldaMod, testingSet)

predDF <- data.frame(predRF, predGbm, predlda, type=testingSet$TYPE)
predDF

CombinedModels <- train(type~., method='gam', data=predDF)
CombinedPredictions <- predict(CombinedModels, predDF)
CombinedPredictions

sum <- sum(CombinedPredictions==testingSet$TYPE)#33
length <- length(testingSet$TYPE)#36
accuracy_CP1 <- sum/length #0.9167

sum <- sum(predRF==testingSet$TYPE) #31
length <- length(testingSet$TYPE)#36
accuracy_rfMod <- (sum/length) #0.8611

sum <- sum(predGbm==testingSet$TYPE)#31
accuracy_gbmMod <- sum/length # 0.8611

sum <- sum(predlda==testingSet$TYPE) #31
accuracy_ldaMod <- sum/length #0.8611

rf2 <- randomForest(TYPE~., data=trainingSet, method='class')

plot(rf2)

predRF2 <- predict(rf2, testingSet, type='class')

sum <- sum(predRF2==testingSet$TYPE)#31
accuracy_RF2 <- sum/length # 0.8611

confusionMatrix(predRF2, testingSet$TYPE)
predDF2 <- data.frame(predRF, predRF2, predlda, predGbm, TYPE=testingSet$TYPE)

CombinedModels2 <- train(TYPE~., method='gam', data=predDF2)
CombinedPredictions2 <- predict(CombinedModels2, predDF2)
CombinedPredictions2

sum <- sum(CombinedPredictions2==testingSet$TYPE)#33
accuracy_CombinedPredictions2 <- sum/length #0.9167

predDF3 <- data.frame(predRF, predRF2, predlda, predGbm, 
                      CombinedPredictions2, TYPE=testingSet$TYPE)

results <- c(round(accuracy_rfMod,2), round(accuracy_RF2,2), 
             round(accuracy_ldaMod, 2), round(accuracy_gbmMod,2), 
             round(accuracy_CombinedPredictions2,2), round(100,2))

results <- as.factor(results)
results <- t(data.frame(results))
colnames(results) <- colnames(predDF3)
Results <- rbind(predDF3, results) #37X6


write.csv(Results,'Results_predictions_universe16_fold.csv', row.names=TRUE)#37X6, 1.7kb

#try 'rpart', 'knn', 'glm', 
knnMod <- train(TYPE ~ .,
                method='knn', preProcess=c('center','scale'),
                tuneLength=10, trControl=trainControl(method='cv'), data=trainingSet)#list 23
plot(knnMod)

rpartMod <- train(TYPE ~ ., method='rpart', tuneLength=9, data=trainingSet) #lists 23
plot(rpart)

glmMod <- train(TYPE ~ ., 
                method='glm', data=trainingSet) #list of 23, 44 warnings
#plot(glmMod)#errors

predKNN <- predict(knnMod, testingSet)
predRPART <- predict(rpartMod, testingSet)
predGLM <- predict(glmMod, testingSet)



length=length(testingSet$TYPE)#36

sumKNN <- sum(predKNN==testingSet$TYPE)#31
sumRPart <- sum(predRPART==testingSet$TYPE)#26
sumGLM <- sum(predGLM==testingSet$TYPE)#31

accuracy_KNN <- sumKNN/length #0.8611
accuracy_RPART <- sumRPart/length #0.8611
accuracy_GLM <- sumGLM/length #0.7222

predDF3 <- data.frame(predDF2[,1:4],predKNN,predRPART,predGLM, 
                      TYPE=testingSet$TYPE)

CombinedModels <- train(TYPE ~ ., method='gam', data=predDF3)
CombinedPredictions2 <- predict(CombinedModels, predDF3)
accuracy_CP2 <- sum(CombinedPredictions2==testingSet$TYPE)/length #0.9167

predDF4 <- data.frame(predDF3[,1:7], CombinedPredictions2, TYPE=testingSet$TYPE)#36X9
colnames(predDF4)
# [1] "predRF"               "predRF2"              "predlda"              "predGbm"             
# [5] "predKNN"              "predRPART"            "predGLM"              "CombinedPredictions2"
# [9] "TYPE"

results <- c(round(accuracy_rfMod,2), round(accuracy_RF2,2), 
             round(accuracy_ldaMod, 2), round(accuracy_gbmMod,2), 
             round(accuracy_KNN,2), round(accuracy_RPART,2),
             round(accuracy_GLM,2), 
             round(accuracy_CP2,2), round(100,2))
results <- as.factor(results)
results <- t(data.frame(results))#1X6
colnames(results) <- colnames(predDF4)
Results <- rbind(predDF4, results) #37X9

write.csv(Results,'Results_predictions_universe16_fold_8_algorithms_used.csv', 
          row.names=TRUE)#37X9, 2.7kb

###############################################################################

# What about the 12,173 genes all in common among the five series before filtering for 
# the four cytoband locations, what genes are the top genes in DE being the 
# greatest in UL compared to non-UL?

universeDE <- read.csv('all_common_12173_130_fold_magnitude.csv', 
                       sep=',', header=TRUE, 
                       na.strings=c('','NA'), row.names=1)#12173X129


universeOrdered <- universeDE[order(universeDE$Magnitude, decreasing=TRUE),]
row.names(universeOrdered) <- universeOrdered$GENE_SYMBOL

universe16 <- universeOrdered[1:16,-c(1:9)]

universe_t <- data.frame(t(universe16)) # 121 X 16 genes

row.names(universe_t) #they are ordered with the first 51 non-UL and next 70 UL samples

# create the 'TYPE' field to identify each sample as UL or nonUL
ul <- data.frame(TYPE=rep('UL', 70))
non <- data.frame(TYPE=rep('nonUL',51))
TYPE <- rbind(non,ul) #121X1

# add the TYPE field to the universal 16 genes having the most DE
UNIVERSE <- cbind(TYPE, universe_t) #121X17

write.csv(UNIVERSE, 'most_universe_DE.csv', row.names=TRUE)

colnames(UNIVERSE)
# [1] "TYPE"    "HSPB1"   "DSTN"    "S100A6"  "CNN1"    "ACTG2"   "VIM"     "SPARCL1" "TPM2"    "ACTA2"  
#[11] "PCP4"    "TAGLN"   "DES"     "RAMP1"   "CYR61"   "UBC"     "ACTB"


set.seed(189678345) # this will reproduce the same numbers each time sampling is done
# run set.seed value to get exact results, otherwise results differ

# creates a partition of the data by indices
inTrain <- createDataPartition(y=UNIVERSE$TYPE, p=0.7, list=FALSE)#85

trainingSet <- UNIVERSE[inTrain,]#85X7
testingSet <- UNIVERSE[-inTrain,]#36X7

# randomForest, cross-validation (cv) = 5
rfMod <- train(TYPE~., method='rf', data=(trainingSet), 
               trControl=trainControl(method='cv'), number=5)
plot(rfMod)

# generalizedBoostedModel
gbmMod <- train(TYPE~., method='gbm', data=trainingSet, verbose=FALSE )

plot(gbmMod)


# linkage dirichlet allocation model
ldaMod <- train(TYPE~., method='lda', data=trainingSet)

# the ldaMod cannot be plotted, 'no tuning parameters for this model'

# run predictions on the testing set
predRF <- predict(rfMod, testingSet)
predGbm <- predict(gbmMod, testingSet)
predlda <- predict(ldaMod, testingSet)

predDF <- data.frame(predRF, predGbm, predlda, type=testingSet$TYPE)
predDF

CombinedModels <- train(type~., method='gam', data=predDF)
CombinedPredictions <- predict(CombinedModels, predDF)
CombinedPredictions

sum <- sum(CombinedPredictions==testingSet$TYPE)#29
length <- length(testingSet$TYPE)#36
accuracy_CP1 <- sum/length #0.8056

sum <- sum(predRF==testingSet$TYPE) #24
length <- length(testingSet$TYPE)#36
accuracy_rfMod <- (sum/length) #0.6667

sum <- sum(predGbm==testingSet$TYPE)#25
accuracy_gbmMod <- sum/length # 0.6944

sum <- sum(predlda==testingSet$TYPE) #29
accuracy_ldaMod <- sum/length #0.8055

rf2 <- randomForest(TYPE~., data=trainingSet, method='class')

plot(rf2)

predRF2 <- predict(rf2, testingSet, type='class')

sum <- sum(predRF2==testingSet$TYPE)#23
accuracy_RF2 <- sum/length # 0.6388

confusionMatrix(predRF2, testingSet$TYPE)
predDF2 <- data.frame(predRF, predRF2, predlda, predGbm, TYPE=testingSet$TYPE)

CombinedModels2 <- train(TYPE~., method='gam', data=predDF2)
CombinedPredictions2 <- predict(CombinedModels2, predDF2)
CombinedPredictions2

sum <- sum(CombinedPredictions2==testingSet$TYPE)#30
accuracy_CombinedPredictions2 <- sum/length #0.8333

predDF3 <- data.frame(predRF, predRF2, predlda, predGbm, 
                      CombinedPredictions2, TYPE=testingSet$TYPE)

results <- c(round(accuracy_rfMod,2), round(accuracy_RF2,2), 
             round(accuracy_ldaMod, 2), round(accuracy_gbmMod,2), 
             round(accuracy_CombinedPredictions2,2), round(100,2))

results <- as.factor(results)
results <- t(data.frame(results))
colnames(results) <- colnames(predDF3)
Results <- rbind(predDF3, results) #37X6


write.csv(Results,'Results_predictions_universe16_DE.csv', row.names=TRUE)#37X6, 1.7kb

#try 'rpart', 'knn', 'glm', 
knnMod <- train(TYPE ~ .,
                method='knn', preProcess=c('center','scale'),
                tuneLength=10, trControl=trainControl(method='cv'), data=trainingSet)#list 23
plot(knnMod)

rpartMod <- train(TYPE ~ ., method='rpart', tuneLength=9, data=trainingSet) #lists 23
plot(rpart)

glmMod <- train(TYPE ~ ., 
                method='glm', data=trainingSet) #list of 23, 44 warnings
#plot(glmMod)#errors

predKNN <- predict(knnMod, testingSet)
predRPART <- predict(rpartMod, testingSet)
predGLM <- predict(glmMod, testingSet)



length=length(testingSet$TYPE)#36

sumKNN <- sum(predKNN==testingSet$TYPE)#26
sumRPart <- sum(predRPART==testingSet$TYPE)#19
sumGLM <- sum(predGLM==testingSet$TYPE)#28

accuracy_KNN <- sumKNN/length #0.7222
accuracy_RPART <- sumRPart/length #0.5277
accuracy_GLM <- sumGLM/length #0.7777

predDF3 <- data.frame(predDF2[,1:4],predKNN,predRPART,predGLM, 
                      TYPE=testingSet$TYPE)

CombinedModels <- train(TYPE ~ ., method='gam', data=predDF3)
CombinedPredictions2 <- predict(CombinedModels, predDF3)
accuracy_CP2 <- sum(CombinedPredictions2==testingSet$TYPE)/length #0.8611

predDF4 <- data.frame(predDF3[,1:7], CombinedPredictions2, TYPE=testingSet$TYPE)#36X9
colnames(predDF4)
# [1] "predRF"               "predRF2"              "predlda"              "predGbm"             
# [5] "predKNN"              "predRPART"            "predGLM"              "CombinedPredictions2"
# [9] "TYPE"

results <- c(round(accuracy_rfMod,2), round(accuracy_RF2,2), 
             round(accuracy_ldaMod, 2), round(accuracy_gbmMod,2), 
             round(accuracy_KNN,2), round(accuracy_RPART,2),
             round(accuracy_GLM,2), 
             round(accuracy_CP2,2), round(100,2))
results <- as.factor(results)
results <- t(data.frame(results))#1X6
colnames(results) <- colnames(predDF4)
Results <- rbind(predDF4, results) #37X9

write.csv(Results,'Results_predictions_universe16_DE_8_algorithms_used.csv', 
          row.names=TRUE)#37X9, 2.7kb
################################################################################
# What about the 12,173 genes all in common among the five series before filtering for 
# the four cytoband locations, what genes are the lowest genes in differential
# expression being the 
# lowest in UL compared to non-UL?

universeDE <- read.csv('all_common_12173_130_fold_magnitude.csv', 
                       sep=',', header=TRUE, 
                       na.strings=c('','NA'), row.names=1)#12173X129


universeOrdered <- universeDE[order(universeDE$Magnitude, decreasing=FALSE),]
row.names(universeOrdered) <- universeOrdered$GENE_SYMBOL

universe16 <- universeOrdered[1:16,-c(1:9)]

universe_t <- data.frame(t(universe16)) # 121 X 16 genes

row.names(universe_t) #they are ordered with the first 51 non-UL and next 70 UL samples

# create the 'TYPE' field to identify each sample as UL or nonUL
ul <- data.frame(TYPE=rep('UL', 70))
non <- data.frame(TYPE=rep('nonUL',51))
TYPE <- rbind(non,ul) #121X1

# add the TYPE field to the universal 16 genes having the most fold change
UNIVERSE <- cbind(TYPE, universe_t) #121X17

write.csv(UNIVERSE, 'least_universe_DE.csv', row.names=TRUE)

colnames(UNIVERSE)
# [1] "TYPE"    "USP32P2" "RCVRN"   "SYNGR3"  "MORC1"   "KLK2"    "SUV39H1" "LIG4"   
# [9] "KLHDC4"  "GRIK4"   "FABP1"   "TLX3"    "LAMB4"   "DNTT"    "VN1R1"   "LEFTY1" 
# [17] "C7orf64"

set.seed(189678345) # this will reproduce the same numbers each time sampling is done
# run set.seed value to get exact results, otherwise results differ

# creates a partition of the data by indices
inTrain <- createDataPartition(y=UNIVERSE$TYPE, p=0.7, list=FALSE)#85

trainingSet <- UNIVERSE[inTrain,]#85X7
testingSet <- UNIVERSE[-inTrain,]#36X7

# randomForest, cross-validation (cv) = 5
rfMod <- train(TYPE~., method='rf', data=(trainingSet), 
               trControl=trainControl(method='cv'), number=5)
plot(rfMod)

# generalizedBoostedModel
gbmMod <- train(TYPE~., method='gbm', data=trainingSet, verbose=FALSE )

plot(gbmMod)


# linkage dirichlet allocation model
ldaMod <- train(TYPE~., method='lda', data=trainingSet)

# the ldaMod cannot be plotted, 'no tuning parameters for this model'

# run predictions on the testing set
predRF <- predict(rfMod, testingSet)
predGbm <- predict(gbmMod, testingSet)
predlda <- predict(ldaMod, testingSet)

predDF <- data.frame(predRF, predGbm, predlda, type=testingSet$TYPE)
predDF

CombinedModels <- train(type~., method='gam', data=predDF)
CombinedPredictions <- predict(CombinedModels, predDF)
CombinedPredictions

sum <- sum(CombinedPredictions==testingSet$TYPE)#27
length <- length(testingSet$TYPE)#36
accuracy_CP1 <- sum/length #0.75

sum <- sum(predRF==testingSet$TYPE) #11
length <- length(testingSet$TYPE)#36
accuracy_rfMod <- (sum/length) #0.3056

sum <- sum(predGbm==testingSet$TYPE)#17
accuracy_gbmMod <- sum/length # 0.4722

sum <- sum(predlda==testingSet$TYPE) #12
accuracy_ldaMod <- sum/length #0.3333

rf2 <- randomForest(TYPE~., data=trainingSet, method='class')

plot(rf2)

predRF2 <- predict(rf2, testingSet, type='class')

sum <- sum(predRF2==testingSet$TYPE)#12
accuracy_RF2 <- sum/length # 0.3333

confusionMatrix(predRF2, testingSet$TYPE)
predDF2 <- data.frame(predRF, predRF2, predlda, predGbm, TYPE=testingSet$TYPE)

CombinedModels2 <- train(TYPE~., method='gam', data=predDF2)
CombinedPredictions2 <- predict(CombinedModels2, predDF2)
CombinedPredictions2

sum <- sum(CombinedPredictions2==testingSet$TYPE)#27
accuracy_CombinedPredictions2 <- sum/length #0.75

predDF3 <- data.frame(predRF, predRF2, predlda, predGbm, 
                      CombinedPredictions2, TYPE=testingSet$TYPE)

results <- c(round(accuracy_rfMod,2), round(accuracy_RF2,2), 
             round(accuracy_ldaMod, 2), round(accuracy_gbmMod,2), 
             round(accuracy_CombinedPredictions2,2), round(100,2))

results <- as.factor(results)
results <- t(data.frame(results))
colnames(results) <- colnames(predDF3)
Results <- rbind(predDF3, results) #37X6


write.csv(Results,'Results_predictions_universe16_DE_least.csv', row.names=TRUE)#37X6, 1.7kb

#try 'rpart', 'knn', 'glm', 
knnMod <- train(TYPE ~ .,
                method='knn', preProcess=c('center','scale'),
                tuneLength=10, trControl=trainControl(method='cv'), data=trainingSet)#list 23
plot(knnMod)

rpartMod <- train(TYPE ~ ., method='rpart', tuneLength=9, data=trainingSet) #lists 23
plot(rpart)

glmMod <- train(TYPE ~ ., 
                method='glm', data=trainingSet) #list of 23, 44 warnings
#plot(glmMod)#errors

predKNN <- predict(knnMod, testingSet)
predRPART <- predict(rpartMod, testingSet)
predGLM <- predict(glmMod, testingSet)



length=length(testingSet$TYPE)#36

sumKNN <- sum(predKNN==testingSet$TYPE)#20
sumRPart <- sum(predRPART==testingSet$TYPE)#21
sumGLM <- sum(predGLM==testingSet$TYPE)#15

accuracy_KNN <- sumKNN/length #0.5556
accuracy_RPART <- sumRPart/length #0.5833
accuracy_GLM <- sumGLM/length #0.4167

predDF3 <- data.frame(predDF2[,1:4],predKNN,predRPART,predGLM, 
                      TYPE=testingSet$TYPE)

CombinedModels <- train(TYPE ~ ., method='gam', data=predDF3)
CombinedPredictions2 <- predict(CombinedModels, predDF3)
accuracy_CP2 <- sum(CombinedPredictions2==testingSet$TYPE)/length #0.8611

predDF4 <- data.frame(predDF3[,1:7], CombinedPredictions2, TYPE=testingSet$TYPE)#36X9
colnames(predDF4)
# [1] "predRF"               "predRF2"              "predlda"              "predGbm"             
# [5] "predKNN"              "predRPART"            "predGLM"              "CombinedPredictions2"
# [9] "TYPE"

results <- c(round(accuracy_rfMod,2), round(accuracy_RF2,2), 
             round(accuracy_ldaMod, 2), round(accuracy_gbmMod,2), 
             round(accuracy_KNN,2), round(accuracy_RPART,2),
             round(accuracy_GLM,2), 
             round(accuracy_CP2,2), round(100,2))
results <- as.factor(results)
results <- t(data.frame(results))#1X6
colnames(results) <- colnames(predDF4)
Results <- rbind(predDF4, results) #37X9

write.csv(Results,'Results_predictions_universe16_DE_least_8_algorithms_used.csv', 
          row.names=TRUE)#37X9, 2.7kb



################################################################################
# Table of all the results for the different variations, used as UL predictors
TOP16 <- read.csv('Results_predictions_TOP16_8_algorithms.csv',
                  sep=',', header=TRUE, row.names=1)
TOP16 <- TOP16[36,]
row.names(TOP16) <- 'TOP16_results'

# the 16 genes having the most magnitude of change in 130
DE16_130 <- read.csv('Results_predictions_DE16_8_algorithms_used.csv',
                     sep=',', header=TRUE, row.names=1)
DE16_130 <- DE16_130[37,]
row.names(DE16_130) <- 'DE16_most_130_results'

# the 16 genes having the least magnitude of change in 130
DE16_130_least <- read.csv('Results_predictions_Least_DE16_8_algorithms_used.csv',
                           sep=',', header=TRUE, row.names=1)
DE16_130_least <- DE16_130_least[37,]
row.names(DE16_130_least) <- 'DE16_least_130_results'


# the 10 genes with the highest fold change and the 6 ubiq genes out of 130
FOLD16_130 <- read.csv('Results_predictions_FOLD16_8_algorithms_used.csv',
                       sep=',', header=TRUE, row.names=1)
FOLD16_130 <- FOLD16_130[37,]
row.names(FOLD16_130) <- 'FOLD16_130_results'


# the majority of the 130 genes on 4 chromosomes for 6 ubiq genes, 5 up and 5 down
majority10_5_5_130 <- read.csv('Results_predictions_majority10_8_algorithms_used.csv',
                               sep=',', header=TRUE, row.names=1)
majority10_130 <- majority10_5_5_130[37,]
row.names(majority10_130) <- 'majority_10_results'

# all 12173 genes by 16 highest fold change
universe16_FOLD <- read.csv('Results_predictions_universe16_fold_8_algorithms_used.csv',
                            sep=',', header=TRUE, row.names=1)
universe16_fold <- universe16_FOLD[37,]
row.names(universe16_fold) <- 'universe16_fold_results'


# the 16 genes in the 12173 genes with the highest magnitude of change
universe16_DE <- read.csv('Results_predictions_universe16_DE_8_algorithms_used.csv',
                          sep=',', header=TRUE, row.names=1)
universe16_DE <- universe16_DE[37,]
row.names(universe16_DE) <- 'universe16_DE_most_results'

# the 16 genes in the 12173 genes with the least magnitude of change
universe16_DE_least <- read.csv('Results_predictions_universe16_DE_least_8_algorithms_used.csv',
                                sep=',', header=TRUE, row.names=1)
universe16_DE_least <- universe16_DE_least[37,]
row.names(universe16_DE_least) <- 'universe16_DE_least_results'


RESULTS <- rbind(TOP16, DE16_130, DE16_130_least, FOLD16_130, majority10_130,
                 universe16_fold, universe16_DE, universe16_DE_least)
write.csv(RESULTS, 'results_8_algorithms_8_data_sets.csv', row.names=TRUE)
########################################################################################

# This is the script for mapping out the populations in ggplot2 used in the presentation

library(ggplot2)
#install.packages('maps')
#library(maps)

# https://www.countries-ofthe-world.com
Europe <- c('Albania',
            'Andorra',
            'Armenia',
            'Austria',
            'Azerbaijan',
            'Belarus',
            'Belgium',
            'Bosnia',
            'Bulgaria',
            'Croatia',
            'Cyprus',
            'Czechia',
            'Denmark',
            'Estonia',
            'Finland',
            'France',
            'Georgia',
            'Germany',
            'Greece',
            'Hungary',
            'Iceland',
            'Ireland',
            'Italy',
            'Kazakhstan',
            'Kosovo',
            'Latvia',
            'Liechtenstein',
            'Lithuania',
            'Luxembourg',
            'Malta',
            'Moldova',
            'Monaco',
            'Montenegro',
            'Netherlands',
            'North Macedonia',
            'Norway',
            'Poland',
            'Portugal',
            'Romania',
            'Russia',
            'San Marino',
            'Serbia',
            'Slovakia',
            'Slovenia',
            'Spain',
            'Sweden',
            'Switzerland',
            'Turkey',
            'Ukraine',
            'UK',
            'Vatican City')
Africa <- c('Algeria',
            'Angola',
            'Benin',
            'Botswana',
            'Burkina Faso',
            'Burundi',
            'Cabo Verde',
            'Cameroon',
            'Central African Republic',
            'Chad',
            'Comoros',
            'Congo', 
            'Democratic Republic of the Congo', 
            'Republic of the Cote d Ivoire',
            'Djibouti',
            'Egypt',
            'Equatorial Guinea',
            'Eritrea',
            'Eswatini',
            'Ethiopia',
            'Gabon',
            'Gambia',
            'Ghana',
            'Guinea',
            'Guinea-Bissau',
            'Kenya',
            'Lesotho',
            'Liberia',
            'Libya',
            'Madagascar',
            'Malawi',
            'Mali',
            'Mauritania',
            'Mauritius',
            'Morocco',
            'Mozambique',
            'Namibia',
            'Niger',
            'Nigeria',
            'Rwanda',
            'Sao Tome and Principe',
            'Senegal',
            'Seychelles',
            'Sierra Leone',
            'Somalia',
            'South Africa',
            'South Sudan',
            'Sudan',
            'Tanzania',
            'Togo',
            'Tunisia',
            'Uganda',
            'Zambi',
            'Zimbabwe')
Asia <- c('Afghanistan',
          'Armenia',
          'Azerbaijan',
          'Bahrain',
          'Bangladesh',
          'Bhutan',
          'Brunei',
          'Cambodia',
          'China',
          'Cyprus',
          'Georgia',
          'India',
          'Indonesia',
          'Iran',
          'Iraq',
          'Israel',
          'Japan',
          'Jordan',
          'Kazakhstan',
          'Kuwait',
          'Kyrgyzstan',
          'Laos',
          'Lebanon',
          'Malaysia',
          'Maldives',
          'Mongolia',
          'Myanmar',
          'Nepal',
          'North Korea',
          'Oman',
          'Pakistan',
          'Palestine',
          'Philippines',
          'Qatar',
          'Russia',
          'Saudi Arabia',
          'Singapore',
          'South Korea',
          'Sri Lanka',
          'Syria',
          'Taiwan',
          'Tajikistan',
          'Thailand',
          'Timor-Leste',
          'Turkey',
          'Turkmenistan',
          'UAE',
          'Uzbekistan',
          'Vietnam',
          'Yemen'
)

TNRC6B= c(Europe,'Australia','Saudi Arabia','China',
          'Japan', 'USA')
BET1L = c(Europe, 'China','Japan','USA','Saudi Arabia')
FASN = c(Europe, 'Australia','USA')
HMGA2 = c('USA',Europe, Asia)# cell samples 7 EuroAmerican 2 AsianAmerican, 9 total
CYTH4 = c('USA',Africa)
CCDC57= c(Europe, 'USA', 'Australia')


T <- data.frame(TNRC6B,gene='TNRC6B')
colnames(T)[1] <- 'country'
B <- data.frame(BET1L, gene='BET1L')
colnames(B)[1] <- 'country'
Cy <- data.frame(CYTH4, gene='CYTH4')
colnames(Cy)[1] <- 'country'
Cc <- data.frame(CCDC57, gene='CCDC57')
colnames(Cc)[1] <- 'country'
H <- data.frame(HMGA2, gene='HMGA2')
colnames(H)[1] <- 'country'
F <- data.frame(FASN, gene='FASN')
colnames(F)[1] <- 'country'

GENES <- rbind(T,B,Cy,Cc,H,F)
genes <- GENES[complete.cases(GENES),]

world <- map_data("world") 

mapGenes <- merge(world, genes, by.x='region', by.y='country') 
dat <- mapGenes

world <- world[((world$region != "Antarctica") & (world$region != 'South America')),] 

png('UL_Risk_Genes_Population_Maps_Facet.png', width=800, height=600)
gg <- ggplot()
gg <- gg + geom_map(data=world, map=world,
                    aes(x=long, y=lat, map_id=region),
                    color="white", fill="#7f7f7f", size=0.05, alpha=1/4)
gg <- gg + geom_point(data=dat, 
                      aes(x=long, y=lat, color=gene), 
                      size=.5, alpha=.2)
gg <- gg + facet_wrap(~gene)
gg <- gg + theme(legend.position="none")
gg <- gg + xlab('Longitude') + ylab('Latitude') + ggtitle('UL Risk Genes Confirmed in Populations')
gg
dev.off()




