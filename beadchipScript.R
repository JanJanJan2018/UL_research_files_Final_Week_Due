# Examining the Beadchip data

GPL13376 <- read.delim('GPL13376-11269.txt', sep='\t', quote="", header=TRUE,
                      comment.char='#', na.strings=c('','NA'))#48701X30

colnames(GPL13376)
# [1] "ID"                    "Species"               "Source"               
# [4] "Search_Key"            "Transcript"            "ILMN_Gene"            
# [7] "Source_Reference_ID"   "RefSeq_ID"             "Unigene_ID"           
# [10] "Entrez_Gene_ID"        "GI"                    "Accession"            
# [13] "Symbol"                "Protein_Product"       "Probe_Id"             
# [16] "Array_Address_Id"      "Probe_Type"            "Probe_Start"          
# [19] "SEQUENCE"              "Chromosome"            "Probe_Chr_Orientation"
# [22] "Probe_Coordinates"     "Cytoband"              "Definition"           
# [25] "Ontology_Component"    "Ontology_Process"      "Ontology_Function"    
# [28] "Synonyms"              "Obsolete_Probe_Id"     "GB_ACC" 


GPL23976 <- read.delim('GPL23976-16627.txt', sep='\t', quote="", header=TRUE,
                       comment.char='#', na.strings=c('','NA')) #866836X2

colnames(GPL23976)
# [1] "ID"      "SPOT_ID"

#The following platform is to add meta to the GPL23976 platform that excludes it
#using GPL23976 'ID' and GPL13534 'IlmnID'


GPL13534 <- read.csv('GPL13534_HumanMethylation450_15017482_v.1.1.csv', 
                     sep=',', quote="", header=TRUE,skip=7,
                       comment.char='#', na.strings=c('','NA')) # 486428 X 33

colnames(GPL13534)
# "IlmnID"                      "Name"                       
# [3] "AddressA_ID"                 "AlleleA_ProbeSeq"           
# [5] "AddressB_ID"                 "AlleleB_ProbeSeq"           
# [7] "Infinium_Design_Type"        "Next_Base"                  
# [9] "Color_Channel"               "Forward_Sequence"           
# [11] "Genome_Build"                "CHR"                        
# [13] "MAPINFO"                     "SourceSeq"                  
# [15] "Chromosome_36"               "Coordinate_36"              
# [17] "Strand"                      "Probe_SNPs"                 
# [19] "Probe_SNPs_10"               "Random_Loci"                
# [21] "Methyl27_Loci"               "UCSC_RefGene_Name"          
# [23] "UCSC_RefGene_Accession"      "UCSC_RefGene_Group"         
# [25] "UCSC_CpG_Islands_Name"       "Relation_to_UCSC_CpG_Island"
# [27] "Phantom"                     "DMR"                        
# [29] "Enhancer"                    "HMM_Island"                 
# [31] "Regulatory_Feature_Name"     "Regulatory_Feature_Group"   
# [33] "DHS"

platformGPL13534 <- GPL13534[,c(1,10,12,17,22,23,32)]
p13534 <- merge(GPL23976, platformGPL13534, by.x='ID', by.y='IlmnID')
platform13534 <- p13534[,-c(2,7,8)]
pl13534 <- na.omit(platform13534)

gene <- strsplit(as.character(pl13534$UCSC_RefGene_Name), ';') #list of 41108
geneName <- lapply(gene, '[',1) #this takes the 1st item in each field row

colnames(pl13534)
# [1] "ID"                "Forward_Sequence"  "CHR"               "Strand"           
# [5] "UCSC_RefGene_Name"

pl13534$Gene <- geneName

colnames(pl13534)
# [1] "ID"                "Forward_Sequence"  "CHR"               "Strand"           
# [5] "UCSC_RefGene_Name" "Gene"

#all 6 UL risk genes are in the above data set

pl13534 <- pl13534[,-c(5)]

colnames(pl13534)
# [1] "ID"               "Forward_Sequence" "CHR"              "Strand"          
# [5] "Gene" 

pl13534$Gene <- as.factor(as.character(pl13534$Gene))

# Rename PSCD4 to CYTH4, it is the GPL13376 Symbolfield, 
#CYTH4's alternate/previous symbol was PSCD4, 
#it shows up index 37170 but not in pl13534 as PSCD4
GPL13376$Symbol <- gsub('PSCD4','CYTH4', GPL13376$Symbol)
grep('CYTH4', GPL13376$Symbol)#37170


#merge the two platforms by 'Symbol' and 'Gene'
both <-merge(pl13534, GPL13376, by.x='Gene', by.y='Symbol')#400712 X 34

write.csv(both, 'Meta_all_beadchip.csv', row.names=FALSE) # 1.1 GB

Both <- both[,-c(10)]

TNRC6B <- Both[grep('TNRC6B', Both$Gene),]
BET1L <- Both[grep('BET1L', Both$Gene),]
HMGA2 <- Both[grep('HMGA2', Both$Gene),]
FASN <- Both[grep('FASN', Both$Gene),]
CCDC57 <- Both[grep('CCDC57', Both$Gene),]
CYTH4 <- Both[grep('CYTH4', Both$Gene),] 

UL_risk <- rbind(BET1L, CCDC57, CYTH4, FASN, HMGA2, TNRC6B) #267X33
write.csv(UL_risk, 'UL_risk_genes_beadchip.csv', row.names=FALSE)


GSE120854 <- read.delim('GSE120854_series_matrix.txt', sep='\t', header=TRUE,
                          na.strings=c('','NA'), comment.char='!')
GSE95101 <- read.delim('GSE95101_series_matrix.txt', sep='\t', header=TRUE,
                        na.strings=c('','NA'), comment.char='!')


##Keep all the originals and just merge the GPL13534 'IlmnID' with GSE120854  'ID_REF'
##and the GPL13376 with GSE95101

table1 <- merge(GPL13534, GSE120854, by.x='IlmnID', by.y='ID_REF')#452453X67
#GSM3417135-GSM3417144 nonUL, GSM2417145-GSM2417168 UL

#colnames(table1)
nonUL <- table1[,c(10,12,14,17,22,23,34:43)]
UL <- table1[,c(10,12,14,17,22,23,44:67)]
colnames(UL)[7:30] <- paste(colnames(UL)[7:30], sep='', 'UL')

Access <- strsplit(as.character(UL$UCSC_RefGene_Accession),';',fixed=TRUE)
Accession <- lapply(Access, '[',1)

UL$Access <- as.factor(as.character(Accession))
                       
UL <- UL[,c(1:6,31,7:30)]

Gene <- strsplit(as.character(UL$UCSC_RefGene_Name),';', fixed=TRUE)
Symbol <- lapply(Gene, '[',1)
UL$GeneSymbol <-  as.factor(as.character(Symbol))
UL <- UL[,c(1:5,32,6:31)]

Table1 <- cbind(UL,nonUL[,7:16])

table2 <- merge(GPL13376, GSE95101, by.x='ID', by.y='ID_REF')#48701X68
#GSM2496185-GSM2496193 UL, GSM2496194-GSM2496202 nonUL, GSM2496203-GSM249609 UL,
#GSM2496210-GSM2496216 nonUL, GSM2496217-GSM2496219 UL, GSM2496220-GSM249622 nonUL

nonUL2 <- table2[,c(10:13,19:21,23:28,40:48,56:62,66:68)]
UL2 <- table2[,c(10:13,19:21,23:28,31:39,49:55,63:65)]
colnames(UL2)[14:32] <- paste(colnames(UL2)[14:32], sep='','UL')

Accession <- strsplit(as.character(UL2$Accession), '.', fixed=TRUE)
Extension <- lapply(Accession, '[',2)
Access <- lapply(Accession,'[',1)

UL2$Access <-  as.factor(as.character(Access))
UL2$Extension <-  as.factor(as.character(Extension))
UL2 <- UL2[,c(1:3,33,34,4:32)]

Table2 <- cbind(UL2,nonUL2[,14:32])

TABLE <- merge(Table1, Table2, by.x='Access', 
               by.y='Access')

table <- TABLE[,c(1:8,43:56,33:42,76:94,9:32,57:75)]#22 meta, 29 nonUL, 43 UL
#262,519 X 94

table_less <- table[,c(1:5,12:14,17:18,23:94)]#262,519 X 82


write.csv(table, 'table_more.csv', row.names=FALSE)
write.csv(table_less_complete, 'table_less.csv', row.names=FALSE)

table_less <- read.csv('table_less.csv', header=TRUE, na.strings=c('', 'NA','null'),)
table <- na.omit(table_less)

write.csv(table, 'table_less.csv', row.names=FALSE)

######################################################################################
table_more <- read.csv('table_more.csv', header=TRUE, sep=',',
                       na.strings=c('','null','NA'))

table_less <- read.csv('table_less.csv', sep=',', header=TRUE,
                      na.strings=c('','null','NA'))

library(dplyr)

table <- table_less %>% group_by(Symbol) %>% summarise(n=n()) #16,152 genes
accessTable <- table_less %>% group_by(Access) %>% summarise(n=n())#17496 accession IDs

ul <- grep('UL', colnames(table_less))#43
UL <- table_less[,c(7,ul)]#43 no meta
nonUL <- table_less[,-ul]#39 includes meta
meta <- nonUL[,1:10]
names <- meta$Symbol #factor 16,152
nonUL <- nonUL[,-c(1:6,8:10)]#29 no meta

#creates data frames of means and medians for UL and nonUL samples
samples <- colnames(UL[2:44]) #GSE sample IDs
samples <- as.vector(samples)
UL_means <- UL %>% group_by(Symbol) %>% summarise_at(samples, mean, na.rm=TRUE)
row.names(UL_means) <- UL_means$Symbol

samples2 <- colnames(nonUL[2:30]) #GSE sample IDs
samples2 <- as.vector(samples2)
nonUL_means <- nonUL %>% group_by(Symbol) %>% summarise_at(samples2, mean, na.rm=TRUE)
row.names(nonUL_means) <- nonUL_means$Symbol

samples <- colnames(UL[2:44]) #GSE sample IDs
samples <- as.vector(samples)
UL_medians <- UL %>% group_by(Symbol) %>% summarise_at(samples, median, na.rm=TRUE)
row.names(UL_medians) <- UL_medians$Symbol

samples2 <- colnames(nonUL[2:30]) #GSE sample IDs
samples2 <- as.vector(samples2)
nonUL_medians <- nonUL %>% group_by(Symbol) %>% summarise_at(samples2, median, na.rm=TRUE)
row.names(nonUL_medians) <- nonUL_medians$Symbol

UL_rowMeans <- data.frame(rowMeans(UL_means[2:44],1))
colnames(UL_rowMeans) <- 'UL_Means'
nonUL_rowMeans <- data.frame(rowMeans(nonUL_means[2:30],1))
colnames(nonUL_rowMeans) <- 'nonUL_Means'
UL_rowMedians <- data.frame(apply(UL_medians[2:44],1, median))
colnames(UL_rowMedians) <- 'UL_Medians'
nonUL_rowMedians <- data.frame(apply(nonUL_medians[2:30],1,median))
colnames(nonUL_rowMedians) <- 'nonUL_Medians'


Stats <- cbind(nonUL_rowMeans, UL_rowMeans, nonUL_rowMedians, UL_rowMedians)
row.names(Stats) <- UL_means$Symbol

StatsFC <- mutate(Stats, FoldChange=UL_Means/nonUL_Means)
StatsDE <- mutate(StatsFC, DifferentialExpression=UL_Means-nonUL_Means)
row.names(StatsDE) <- UL_means$Symbol

StatsMagnitude <- mutate(StatsDE, Magnitude=abs(DifferentialExpression))

StatsMedianFC <- mutate(StatsMagnitude, 
                        FoldChangeMedian=UL_Medians/nonUL_Medians)
StatsMedianDE <- mutate(StatsMedianFC, Median_DE=UL_Medians-nonUL_Medians)
StatsMedianMagnitude <- mutate(StatsMedianDE, MedianMagnitude=abs(Median_DE))
row.names(StatsMedianMagnitude) <- UL_means$Symbol

StatsAll <- cbind(StatsMedianMagnitude,nonUL_means, UL_means[2:44])

write.csv(StatsAll, 'StatsAllCompressed.csv', row.names=TRUE)

StatsAll <- read.csv('StatsAllCompressed.csv', sep=',', header=TRUE, row.names=1)

StatsOrdered <- StatsAll[order(StatsAll$Magnitude, decreasing=TRUE),]
Mag20 <- StatsOrdered[1:20,]
write.csv(Mag20, 'magnitude_means_20.csv', row.names=TRUE)

FC_ordered <- StatsAll[order(StatsAll$FoldChange, decreasing=TRUE),]
FC20 <- FC_ordered[1:20,]
write.csv(FC20, 'FoldChange_means_20.csv', row.names=TRUE)

FC_med_ordered <- StatsAll[order(StatsAll$FoldChangeMedian, decreasing=TRUE),]
FC_med_20 <- FC_med_ordered[1:20,]
write.csv(FC_med_20, 'FoldChange_median_20.csv', row.names=TRUE)

mag_med_ordered <- StatsAll[order(StatsAll$MedianMagnitude, decreasing=TRUE),]
mag_med_20 <- mag_med_ordered[1:20,]
write.csv(mag_med_20, 'magnitude_median_20.csv', row.names=TRUE)

names_mag_med_20 <- row.names(mag_med_20)
names_FC_med_20 <- row.names(FC_med_20)
names_FC20 <- row.names(FC20)
names_Mag20 <- row.names(Mag20)

m_mag_20 <- t(mag_med_20[,-c(1:11)])
m_fc_20 <- t(FC_med_20[,-c(1:11)])
FC_20 <- t(FC20[,-c(1:11)])
mag_20 <- t(Mag20[,-c(1:11)])

ul <- data.frame(as.character(rep('UL',43)))
colnames(ul) <- 'TYPE'
non <- data.frame(as.character(rep('nonUL', 29)))
colnames(non) <- 'TYPE'

TYPE <- rbind(non,ul)

med_mag_20 <- cbind(TYPE, m_mag_20)
med_FC_20 <- cbind(TYPE, m_fc_20)
mean_FC_20 <- cbind(TYPE, FC_20)
mean_mag_20 <- cbind(TYPE, mag_20)

write.csv(med_mag_20, 'ml_ready_med_mag_20.csv', row.names=TRUE)
write.csv(med_FC_20, 'ml_ready_med_FC_20.csv', row.names=TRUE)
write.csv(mean_FC_20, 'ml_ready_mean_FC_20.csv', row.names=TRUE)
write.csv(mean_mag_20, 'ml_ready_mean_mag_20.csv', row.names=TRUE)

library(caret)
library(randomForest)
library(MASS)
library(gbm)
library(e1071)

set.seed(189678345)

# creates a partition of the data by indices
inTrain <- createDataPartition(y=mean_FC_20$TYPE, p=0.7, list=FALSE)


trainingSet <- mean_FC_20[inTrain,]#52X21
testingSet <- mean_FC_20[-inTrain,]#20X21


rfMod <- train(TYPE~., method='rf', data=(trainingSet), 
               trControl=trainControl(method='cv'), number=5)

gbmMod <- train(TYPE~., method='gbm', data=trainingSet, verbose=FALSE )
ldaMod <- train(TYPE~., method='lda', data=trainingSet)

predRF <- predict(rfMod, testingSet)
predGbm <- predict(gbmMod, testingSet)
predlda <- predict(ldaMod, testingSet)

predDF <- data.frame(predRF, predGbm, predlda, type=testingSet$TYPE)
predDF
# predRF predGbm predlda  type
# 1   nonUL      UL      UL nonUL
# 2   nonUL      UL      UL nonUL
# 3   nonUL      UL      UL nonUL
# 4      UL      UL      UL nonUL
# 5   nonUL   nonUL      UL nonUL
# 6   nonUL   nonUL   nonUL nonUL
# 7   nonUL   nonUL   nonUL nonUL
# 8   nonUL   nonUL   nonUL nonUL
# 9   nonUL      UL      UL    UL
# 10     UL      UL      UL    UL
# 11     UL      UL      UL    UL
# 12     UL      UL      UL    UL
# 13     UL      UL      UL    UL
# 14     UL      UL      UL    UL
# 15     UL      UL      UL    UL
# 16     UL      UL      UL    UL
# 17  nonUL      UL      UL    UL
# 18     UL      UL      UL    UL
# 19  nonUL   nonUL   nonUL    UL
# 20     UL   nonUL   nonUL    UL

CombinedModels <- train(type~., method='gam', data=predDF)
CombinedPredictions <- predict(CombinedModels, predDF)

sum <- sum(CombinedPredictions==testingSet$TYPE)#16
length <- length(CombinedPredictions)#20

sum <- sum(predRF==testingSet$TYPE) #16
length <- length(predRF)#20
accuracy_rfMOd <- (sum/length) #[1] 0.8

sum <- sum(predGbm==testingSet$TYPE)#14
accuracy_gbmMod <- sum/length # [1] 0.7


sum <- sum(predlda==testingSet$TYPE) #13
accuracy_ldaMod <- sum/length #[1] 0.65

rf2 <- randomForest(TYPE~., data=trainingSet, method='class')
predRF2 <- predict(rf2, testingSet, type='class')

sum <- sum(predRF2==testingSet$TYPE)#17
accuracy_RF2 <- sum/length #[1] 0.85

confusionMatrix(predRF2, testingSet$TYPE)

predDF2 <- data.frame(predRF,predRF2,predlda, predGbm, testingSet$TYPE)
colnames(predDF2)[5] <- 'TYPE'

CombinedModels <- train(TYPE ~ ., method='gam', data=predDF2)
CombinedPredictions2 <- predict(CombinedModels, predDF)
sum <- sum(CombinedPredictions2==testingSet$TYPE)#17
length(CombinedPredictions2)#20
accuracy_CombinedPredictions2 <- sum/length #[1] 0.85

predDF3 <- data.frame(predRF,predRF2,predlda, predGbm, 
                      CombinedPredictions2, testingSet$TYPE)#20X6
colnames(predDF3)[6] <- 'TYPE'

knnMod <- train(TYPE ~ .,
                method='knn', preProcess=c('center','scale'),
                tuneLength=10, trControl=trainControl(method='cv'), 
                data=trainingSet)#list 23
plot(knnMod)

rpartMod <- train(TYPE ~ ., method='rpart', tuneLength=9, data=trainingSet) #lists 23
plot(rpart)

glmMod <- train(TYPE ~ ., 
                method='glm', data=trainingSet) #list of 23

predKNN <- predict(knnMod, testingSet)
predRPART <- predict(rpartMod, testingSet)
predGLM <- predict(glmMod, testingSet)



length=length(testingSet$TYPE)#20

sumKNN <- sum(predKNN==testingSet$TYPE)#17
sumRPart <- sum(predRPART==testingSet$TYPE)#15
sumGLM <- sum(predGLM==testingSet$TYPE)#11

accuracy_KNN <- sumKNN/length #0.85
accuracy_RPART <- sumRPart/length #0.75
accuracy_GLM <- sumGLM/length #0.55

predDF3 <- data.frame(predDF2[,1:4],predKNN,predRPART,predGLM, 
                      TYPE=testingSet$TYPE)

CombinedModels <- train(TYPE ~ ., method='gam', data=predDF3)
CombinedPredictions2 <- predict(CombinedModels, predDF3)
accuracy_CP2 <- sum(CombinedPredictions2==testingSet$TYPE)/length #1

predDF4 <- data.frame(predDF3[,1:7], CombinedPredictions2, TYPE=testingSet$TYPE)#20X9
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
results <- t(data.frame(results))#1X9
colnames(results) <- colnames(predDF4)
Results <- rbind(predDF4, results) #21X9


write.csv(Results, 'FC20means_ml_results_8_algorithms.csv', row.names=TRUE)

##############################################################################################
library(caret)
library(randomForest)
library(MASS)
library(gbm)
library(e1071)

set.seed(189678345)

# creates a partition of the data by indices
inTrain <- createDataPartition(y=med_FC_20$TYPE, p=0.7, list=FALSE)


trainingSet <- med_FC_20[inTrain,]#52X21
testingSet <- med_FC_20[-inTrain,]#20X21

library(e1071)

rfMod <- train(TYPE~., method='rf', data=(trainingSet), 
               trControl=trainControl(method='cv'), number=5)

gbmMod <- train(TYPE~., method='gbm', data=trainingSet, verbose=FALSE )
ldaMod <- train(TYPE~., method='lda', data=trainingSet)

predRF <- predict(rfMod, testingSet)
predGbm <- predict(gbmMod, testingSet)
predlda <- predict(ldaMod, testingSet)

predDF <- data.frame(predRF, predGbm, predlda, type=testingSet$TYPE)


CombinedModels <- train(type~., method='gam', data=predDF)
CombinedPredictions <- predict(CombinedModels, predDF)

sum <- sum(CombinedPredictions==testingSet$TYPE)
length <- length(CombinedPredictions)

sum <- sum(predRF==testingSet$TYPE) 
length <- length(predRF)
accuracy_rfMOd <- (sum/length) 

sum <- sum(predGbm==testingSet$TYPE)
accuracy_gbmMod <- sum/length 


sum <- sum(predlda==testingSet$TYPE) 
accuracy_ldaMod <- sum/length 

rf2 <- randomForest(TYPE~., data=trainingSet, method='class')
predRF2 <- predict(rf2, testingSet, type='class')

sum <- sum(predRF2==testingSet$TYPE)
accuracy_RF2 <- sum/length 

confusionMatrix(predRF2, testingSet$TYPE)

predDF2 <- data.frame(predRF,predRF2,predlda, predGbm, testingSet$TYPE)
colnames(predDF2)[5] <- 'TYPE'

CombinedModels <- train(TYPE ~ ., method='gam', data=predDF2)
CombinedPredictions2 <- predict(CombinedModels, predDF)
sum <- sum(CombinedPredictions2==testingSet$TYPE)
length(CombinedPredictions2)
accuracy_CombinedPredictions2 <- sum/length 

predDF3 <- data.frame(predRF,predRF2,predlda, predGbm, 
                      CombinedPredictions2, testingSet$TYPE)
colnames(predDF3)[6] <- 'TYPE'

knnMod <- train(TYPE ~ .,
                method='knn', preProcess=c('center','scale'),
                tuneLength=10, trControl=trainControl(method='cv'), 
                data=trainingSet)
plot(knnMod)

rpartMod <- train(TYPE ~ ., method='rpart', tuneLength=9, data=trainingSet) 
plot(rpart)

glmMod <- train(TYPE ~ ., 
                method='glm', data=trainingSet) 

predKNN <- predict(knnMod, testingSet)
predRPART <- predict(rpartMod, testingSet)
predGLM <- predict(glmMod, testingSet)



length=length(testingSet$TYPE)

sumKNN <- sum(predKNN==testingSet$TYPE)
sumRPart <- sum(predRPART==testingSet$TYPE)
sumGLM <- sum(predGLM==testingSet$TYPE)

accuracy_KNN <- sumKNN/length 
accuracy_RPART <- sumRPart/length 
accuracy_GLM <- sumGLM/length 

predDF3 <- data.frame(predDF2[,1:4],predKNN,predRPART,predGLM, 
                      TYPE=testingSet$TYPE)

CombinedModels <- train(TYPE ~ ., method='gam', data=predDF3)
CombinedPredictions2 <- predict(CombinedModels, predDF3)
accuracy_CP2 <- sum(CombinedPredictions2==testingSet$TYPE)/length 

predDF4 <- data.frame(predDF3[,1:7], CombinedPredictions2, TYPE=testingSet$TYPE)
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
results <- t(data.frame(results))
colnames(results) <- colnames(predDF4)
Results <- rbind(predDF4, results) 


write.csv(Results, 'FC20medians_ml_results_8_algorithms.csv', row.names=TRUE)

##############################################################################################
library(caret)
library(randomForest)
library(MASS)
library(gbm)
library(e1071)

set.seed(189678345)

# creates a partition of the data by indices
inTrain <- createDataPartition(y=mean_mag_20$TYPE, p=0.7, list=FALSE)


trainingSet <- mean_mag_20[inTrain,]#52X21
testingSet <- mean_mag_20[-inTrain,]#20X21

library(e1071)

rfMod <- train(TYPE~., method='rf', data=(trainingSet), 
               trControl=trainControl(method='cv'), number=5)

gbmMod <- train(TYPE~., method='gbm', data=trainingSet, verbose=FALSE )
ldaMod <- train(TYPE~., method='lda', data=trainingSet)

predRF <- predict(rfMod, testingSet)
predGbm <- predict(gbmMod, testingSet)
predlda <- predict(ldaMod, testingSet)

predDF <- data.frame(predRF, predGbm, predlda, type=testingSet$TYPE)


CombinedModels <- train(type~., method='gam', data=predDF)
CombinedPredictions <- predict(CombinedModels, predDF)

sum <- sum(CombinedPredictions==testingSet$TYPE)
length <- length(CombinedPredictions)

sum <- sum(predRF==testingSet$TYPE) 
length <- length(predRF)
accuracy_rfMOd <- (sum/length) 

sum <- sum(predGbm==testingSet$TYPE)
accuracy_gbmMod <- sum/length 


sum <- sum(predlda==testingSet$TYPE) 
accuracy_ldaMod <- sum/length 

rf2 <- randomForest(TYPE~., data=trainingSet, method='class')
predRF2 <- predict(rf2, testingSet, type='class')

sum <- sum(predRF2==testingSet$TYPE)
accuracy_RF2 <- sum/length 

confusionMatrix(predRF2, testingSet$TYPE)

predDF2 <- data.frame(predRF,predRF2,predlda, predGbm, testingSet$TYPE)
colnames(predDF2)[5] <- 'TYPE'

CombinedModels <- train(TYPE ~ ., method='gam', data=predDF2)
CombinedPredictions2 <- predict(CombinedModels, predDF)
sum <- sum(CombinedPredictions2==testingSet$TYPE)
length(CombinedPredictions2)
accuracy_CombinedPredictions2 <- sum/length 

predDF3 <- data.frame(predRF,predRF2,predlda, predGbm, 
                      CombinedPredictions2, testingSet$TYPE)
colnames(predDF3)[6] <- 'TYPE'

knnMod <- train(TYPE ~ .,
                method='knn', preProcess=c('center','scale'),
                tuneLength=10, trControl=trainControl(method='cv'), 
                data=trainingSet)
plot(knnMod)

rpartMod <- train(TYPE ~ ., method='rpart', tuneLength=9, data=trainingSet) 
plot(rpart)

glmMod <- train(TYPE ~ ., 
                method='glm', data=trainingSet) 

predKNN <- predict(knnMod, testingSet)
predRPART <- predict(rpartMod, testingSet)
predGLM <- predict(glmMod, testingSet)



length=length(testingSet$TYPE)

sumKNN <- sum(predKNN==testingSet$TYPE)
sumRPart <- sum(predRPART==testingSet$TYPE)
sumGLM <- sum(predGLM==testingSet$TYPE)

accuracy_KNN <- sumKNN/length 
accuracy_RPART <- sumRPart/length 
accuracy_GLM <- sumGLM/length 

predDF3 <- data.frame(predDF2[,1:4],predKNN,predRPART,predGLM, 
                      TYPE=testingSet$TYPE)

CombinedModels <- train(TYPE ~ ., method='gam', data=predDF3)
CombinedPredictions2 <- predict(CombinedModels, predDF3)
accuracy_CP2 <- sum(CombinedPredictions2==testingSet$TYPE)/length 

predDF4 <- data.frame(predDF3[,1:7], CombinedPredictions2, TYPE=testingSet$TYPE)
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
results <- t(data.frame(results))
colnames(results) <- colnames(predDF4)
Results <- rbind(predDF4, results) 


write.csv(Results, 'mag20means_ml_results_8_algorithms.csv', row.names=TRUE)

##############################################################################################
library(caret)
library(randomForest)
library(MASS)
library(gbm)
library(e1071)

set.seed(189678345)

# creates a partition of the data by indices
inTrain <- createDataPartition(y=med_mag_20$TYPE, p=0.7, list=FALSE)


trainingSet <- med_mag_20[inTrain,]#52X21
testingSet <- med_mag_20[-inTrain,]#20X21

library(e1071)

rfMod <- train(TYPE~., method='rf', data=(trainingSet), 
               trControl=trainControl(method='cv'), number=5)

gbmMod <- train(TYPE~., method='gbm', data=trainingSet, verbose=FALSE )
ldaMod <- train(TYPE~., method='lda', data=trainingSet)

predRF <- predict(rfMod, testingSet)
predGbm <- predict(gbmMod, testingSet)
predlda <- predict(ldaMod, testingSet)

predDF <- data.frame(predRF, predGbm, predlda, type=testingSet$TYPE)


CombinedModels <- train(type~., method='gam', data=predDF)
CombinedPredictions <- predict(CombinedModels, predDF)

sum <- sum(CombinedPredictions==testingSet$TYPE)
length <- length(CombinedPredictions)

sum <- sum(predRF==testingSet$TYPE) 
length <- length(predRF)
accuracy_rfMOd <- (sum/length) 

sum <- sum(predGbm==testingSet$TYPE)
accuracy_gbmMod <- sum/length 


sum <- sum(predlda==testingSet$TYPE) 
accuracy_ldaMod <- sum/length 

rf2 <- randomForest(TYPE~., data=trainingSet, method='class')
predRF2 <- predict(rf2, testingSet, type='class')

sum <- sum(predRF2==testingSet$TYPE)
accuracy_RF2 <- sum/length 

confusionMatrix(predRF2, testingSet$TYPE)

predDF2 <- data.frame(predRF,predRF2,predlda, predGbm, testingSet$TYPE)
colnames(predDF2)[5] <- 'TYPE'

CombinedModels <- train(TYPE ~ ., method='gam', data=predDF2)
CombinedPredictions2 <- predict(CombinedModels, predDF)
sum <- sum(CombinedPredictions2==testingSet$TYPE)
length(CombinedPredictions2)
accuracy_CombinedPredictions2 <- sum/length 

predDF3 <- data.frame(predRF,predRF2,predlda, predGbm, 
                      CombinedPredictions2, testingSet$TYPE)
colnames(predDF3)[6] <- 'TYPE'

knnMod <- train(TYPE ~ .,
                method='knn', preProcess=c('center','scale'),
                tuneLength=10, trControl=trainControl(method='cv'), 
                data=trainingSet)
plot(knnMod)

rpartMod <- train(TYPE ~ ., method='rpart', tuneLength=9, data=trainingSet) 
plot(rpart)

glmMod <- train(TYPE ~ ., 
                method='glm', data=trainingSet) 

predKNN <- predict(knnMod, testingSet)
predRPART <- predict(rpartMod, testingSet)
predGLM <- predict(glmMod, testingSet)



length=length(testingSet$TYPE)

sumKNN <- sum(predKNN==testingSet$TYPE)
sumRPart <- sum(predRPART==testingSet$TYPE)
sumGLM <- sum(predGLM==testingSet$TYPE)

accuracy_KNN <- sumKNN/length 
accuracy_RPART <- sumRPart/length 
accuracy_GLM <- sumGLM/length 

predDF3 <- data.frame(predDF2[,1:4],predKNN,predRPART,predGLM, 
                      TYPE=testingSet$TYPE)

CombinedModels <- train(TYPE ~ ., method='gam', data=predDF3)
CombinedPredictions2 <- predict(CombinedModels, predDF3)
accuracy_CP2 <- sum(CombinedPredictions2==testingSet$TYPE)/length 

predDF4 <- data.frame(predDF3[,1:7], CombinedPredictions2, TYPE=testingSet$TYPE)
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
results <- t(data.frame(results))
colnames(results) <- colnames(predDF4)
Results <- rbind(predDF4, results) 


write.csv(Results, 'mag20medians_ml_results_8_algorithms.csv', row.names=TRUE)

############################# READ IN RESULTS TO COMPARE ##################################################

FC20 <- read.csv("FC20means_ml_results_8_algorithms.csv", sep=',', header=TRUE,
                 row.names=1)
FC20_median <- read.csv("FC20medians_ml_results_8_algorithms.csv", sep=',',
                        header=TRUE, row.names=1)
mag20 <- read.csv("mag20means_ml_results_8_algorithms.csv", sep=',', header=TRUE,
                  row.names=1)
mag20_median <- read.csv("mag20medians_ml_results_8_algorithms.csv", sep=',',
                         header=TRUE, row.names=1)

##results##

FC20[21,]
# predRF predRF2 predlda predGbm predKNN predRPART predGLM CombinedPredictions2 TYPE
# results    0.8    0.85    0.65     0.7    0.85      0.75    0.55                    1  100

FC20_median[21,]
# predRF predRF2 predlda predGbm predKNN predRPART predGLM CombinedPredictions2 TYPE
# results   0.85     0.9    0.65     0.5    0.75       0.8    0.55                    1  100

mag20[21,]
# predRF predRF2 predlda predGbm predKNN predRPART predGLM CombinedPredictions2 TYPE
# results   0.85       1     0.7    0.85    0.85       0.7     0.7                    1  100

mag20_median[21,]
# predRF predRF2 predlda predGbm predKNN predRPART predGLM CombinedPredictions2 TYPE
# results   0.85    0.85     0.7     0.8    0.85       0.7     0.7                 0.95  100

################################################################################################
################################################################################################

StatsAll <- read.csv('StatsAllCompressed.csv', sep=',', header=TRUE, row.names=1)
table_less <- read.csv('table_less.csv', header=TRUE, na.strings=c('', 'NA','null'),)

mag20mean <- read.csv('magnitude_means_20.csv', sep=',', header=TRUE, row.names=1)
mag20median <- read.csv('magnitude_median_20.csv', sep=',', header=TRUE, row.names=1)
FC20mean <- read.csv('FoldChange_means_20.csv', sep=',', header=TRUE, row.names=1)
FC20median <- read.csv('FoldChange_median_20.csv', sep=',', header=TRUE, row.names=1)

FC20 <- FC20mean[,c(11,5)]
FC20_median <- FC20median[,c(11,8)]
mag20 <- mag20mean[,c(11,7)]
mag20_median <- mag20median[,c(11,10)]

fold <- merge(FC20, table_less, by.x='Symbol', by.y='Symbol')
fold_median <- merge(FC20_median, table_less, by.x='Symbol', by.y='Symbol')
magnitude <- merge(mag20, table_less, by.x='Symbol', by.y='Symbol')
magnitude_median <- merge(mag20_median, table_less, by.x='Symbol', by.y='Symbol')

Fold <- fold[,c(1:3,8,4:7,9:83)]
Fold_median <- fold_median[,c(1:3,8,4:7,9:83)]
Magnitude <- magnitude[,c(1:3,8,4:7,9:83)]
Magnitude_median <- magnitude_median[,c(1:3,8,4:7,9:83)]
  
write.csv(Fold, 'fold_sequences_20.csv', row.names=FALSE)
write.csv(Fold_median, 'fold_median_sequences_20.csv', row.names=FALSE)
write.csv(Magnitude, 'magnitude_sequences_20.csv', row.names=FALSE)
write.csv(Magnitude_median, 'magnitude_median_sequences_20.csv', row.names=FALSE)

#############################################################################################
#############################################################################################
#############################################################################################


fold <- read.csv('fold_sequences_20.csv', sep=',', header=TRUE)
fold_median <- read.csv('fold_median_sequences_20.csv', sep=',', header=TRUE)
magnitude <- read.csv('magnitude_sequences_20.csv', sep=',', header=TRUE)
magnitude_median <- read.csv('magnitude_median_sequences_20.csv', sep=',', header=TRUE)

table_more <- read.csv('table_more.csv', sep=',', header=TRUE)
StatsAll <- read.csv('StatsAllCompressed.csv', sep=',', header=TRUE, row.names=1)

mag20mean <- read.csv('magnitude_means_20.csv', sep=',', header=TRUE, row.names=1)
mag20median <- read.csv('magnitude_median_20.csv', sep=',', header=TRUE, row.names=1)
FC20mean <- read.csv('FoldChange_means_20.csv', sep=',', header=TRUE, row.names=1)
FC20median <- read.csv('FoldChange_median_20.csv', sep=',', header=TRUE, row.names=1)

fold <- merge(FC20mean, table_more[,1:22], by.x='Symbol', by.y='Symbol')
fold_median <- merge(FC20median, table_more[,1:22], by.x='Symbol', by.y='Symbol')
magnitude <- merge(mag20mean, table_more[,1:22], by.x='Symbol', by.y='Symbol')
magnitude_median <- merge(mag20median, table_more[,1:22], by.x='Symbol', by.y='Symbol')

#add more meta information to these analytics data sets
Fold <- fold[,c(1:11,84:104,12:83)]
Fold_median <- fold_median[,c(1:11,84:104,12:83)]
Magnitude <- magnitude[,c(1:11,84:104,12:83)]
Magnitude_median <- magnitude_median[,c(1:11,84:104,12:83)]

write.csv(Fold,'Fold_Allmeta.csv',row.names=FALSE)
write.csv(Fold_median,'Fold_median_Allmeta.csv', row.names=FALSE)
write.csv(Magnitude, 'Magnitude_Allmeta.csv', row.names=FALSE)
write.csv(Magnitude_median, 'Magnitude_Median_Allmeta.csv', row.names=FALSE)

# Extract which human tissue sees the most expression or inhibition of these genes from
# NCBIgene.org, or UCSC, or ENSEMBLE, or genenames.org
# go to each symbol gene ID in https://www.ncbi.nlm.nih.gov/gene search the symbol, and
# go to 'expression' and 'summary' sections to see a summary and bar chart of tissues the 
# gene is seen in

#you can download the expression values, by ENTREZ gene ID (numeric ID) for each

Fold <- unique(Fold$Entrez_Gene_ID) #look each up at https://www.ncbi.nlm.nih.gov/gene/
# for gene expression HPA RNA-Seq normal tissues
# 70  140458  169044   1515   8788   3321  27190  57214   3872  10319   4319   4320   4826   4885
# 5179   5354   2834  387914   6999  347733
# these files were individually downloaded from expression values in normal tissue of 95 humans
# in 27 tissues beginning with 'GeneID' no separator and ID value i.e. 'GeneID70'

FoldList <- as.character(Fold)
foldpaste <- paste('GeneID', FoldList[1:20], '.txt', sep='')

# Switch to the folder you stored these downloaded files

list.files()
# [1] "GeneID10319.txt"  "GeneID140458.txt" "GeneID1515.txt"   "GeneID169044.txt" "GeneID27190.txt" 
# [6] "GeneID2834.txt"   "GeneID3321.txt"   "GeneID347733.txt" "GeneID3872.txt"   "GeneID387914.txt"
# [11] "GeneID4319.txt"   "GeneID4320.txt"   "GeneID4826.txt"   "GeneID4885.txt"   "GeneID5179.txt"  
# [16] "GeneID5354.txt"   "GeneID57214.txt"  "GeneID6999.txt"   "GeneID70.txt"     "GeneID8788.txt"

d1 <- read.delim(foldpaste[1], sep='\t', quote="", header=TRUE,
                 skip=1, na.strings=c('','NA'))
d2 <- read.delim(foldpaste[2], sep='\t', quote="", header=TRUE,
                 skip=1, na.strings=c('','NA'))
d3 <- read.delim(foldpaste[3], sep='\t', quote="", header=TRUE,
                 skip=1, na.strings=c('','NA'))
d4 <- read.delim(foldpaste[4], sep='\t', quote="", header=TRUE,
                 skip=1, na.strings=c('','NA'))
d5 <- read.delim(foldpaste[5], sep='\t', quote="", header=TRUE,
                 skip=1, na.strings=c('','NA'))
d6 <- read.delim(foldpaste[6], sep='\t', quote="", header=TRUE,
                 skip=1, na.strings=c('','NA'))
d7 <- read.delim(foldpaste[7], sep='\t', quote="", header=TRUE,
                 skip=1, na.strings=c('','NA'))
d8 <- read.delim(foldpaste[8], sep='\t', quote="", header=TRUE,
                 skip=1, na.strings=c('','NA'))
d9 <- read.delim(foldpaste[9], sep='\t', quote="", header=TRUE,
                 skip=1, na.strings=c('','NA'))
d10 <- read.delim(foldpaste[10], sep='\t', quote="", header=TRUE,
                  skip=1, na.strings=c('','NA'))
d11 <- read.delim(foldpaste[11], sep='\t', quote="", header=TRUE,
                  skip=1, na.strings=c('','NA'))
d12 <- read.delim(foldpaste[12], sep='\t', quote="", header=TRUE,
                  skip=1, na.strings=c('','NA'))
d13 <- read.delim(foldpaste[13], sep='\t', quote="", header=TRUE,
                  skip=1, na.strings=c('','NA'))
d14 <- read.delim(foldpaste[14], sep='\t', quote="", header=TRUE,
                  skip=1, na.strings=c('','NA'))
d15 <- read.delim(foldpaste[15], sep='\t', quote="", header=TRUE,
                  skip=1, na.strings=c('','NA'))
d16 <- read.delim(foldpaste[16], sep='\t', quote="", header=TRUE,
                  skip=1, na.strings=c('','NA'))
d17 <- read.delim(foldpaste[17], sep='\t', quote="", header=TRUE,
                  skip=1, na.strings=c('','NA'))
d18 <- read.delim(foldpaste[18], sep='\t', quote="", header=TRUE,
                  skip=1, na.strings=c('','NA'))
d19 <- read.delim(foldpaste[19], sep='\t', quote="", header=TRUE,
                  skip=1, na.strings=c('','NA'))
d20 <- read.delim(foldpaste[20], sep='\t', quote="", header=TRUE,
                  skip=1, na.strings=c('','NA'))

# Switch back to parent directory of all other files

Fold_Entrez_df <- rbind(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,
                        d13,d14,d15,d16,d17,d18,d19,d20)

colnames(Fold_Entrez_df)[1] <- 'Entrez_Gene_ID'
Fold_Entrez_df <- Fold_Entrez_df[,-29]
write.csv(Fold_Entrez_df,'Entrez_Gene_Expressions_Fold20.csv',row.names=FALSE)

ExpressionsFold <- read.csv('Entrez_Gene_Expressions_Fold20.csv', sep=',', header=TRUE)

Fold_meta <- read.csv('Fold_Allmeta.csv', sep=',', header=TRUE)

row.names(ExpressionsFold) <- ExpressionsFold$Entrez_Gene_ID
Expressions <- ExpressionsFold[,-1]
names <- row.names(Expressions)

normalMax <- apply(Expressions, 1, max)

library(dplyr)

FoldExpressions <- mutate(Expressions, NormalMax=normalMax)
row.names(FoldExpressions) <- names
FoldExpressions <- round(FoldExpressions,2)

t1 <- order(FoldExpressions[1,], decreasing=TRUE)[1]
t2 <- order(FoldExpressions[2,], decreasing=TRUE)[1]
t3 <- order(FoldExpressions[3,], decreasing=TRUE)[1]
t4 <- order(FoldExpressions[4,], decreasing=TRUE)[1]
t5 <- order(FoldExpressions[5,], decreasing=TRUE)[1]
t6 <- order(FoldExpressions[6,], decreasing=TRUE)[1]
t7 <- order(FoldExpressions[7,], decreasing=TRUE)[1]
t8 <- order(FoldExpressions[8,], decreasing=TRUE)[1]
t9 <- order(FoldExpressions[9,], decreasing=TRUE)[1]
t10 <- order(FoldExpressions[10,], decreasing=TRUE)[1]
t11 <- order(FoldExpressions[11,], decreasing=TRUE)[1]
t12 <- order(FoldExpressions[12,], decreasing=TRUE)[1]
t13 <- order(FoldExpressions[13,], decreasing=TRUE)[1]
t14 <- order(FoldExpressions[14,], decreasing=TRUE)[1]
t15 <- order(FoldExpressions[15,], decreasing=TRUE)[1]
t16 <- order(FoldExpressions[16,], decreasing=TRUE)[1]
t17 <- order(FoldExpressions[17,], decreasing=TRUE)[1]
t18 <- order(FoldExpressions[18,], decreasing=TRUE)[1]
t19 <- order(FoldExpressions[19,], decreasing=TRUE)[1]
t20 <- order(FoldExpressions[20,], decreasing=TRUE)[1]

TissueMax <- as.vector(c(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20))
Tissue <- colnames(FoldExpressions)[c(TissueMax)]

FoldExpressionValues <- mutate(FoldExpressions, TissueMax=Tissue)
row.names(FoldExpressionValues) <- row.names(FoldExpressions)
FoldExpressionValues$Entrez_Gene_ID <- row.names(FoldExpressionValues)

write.csv(FoldExpressionValues, 'FoldValuesTissueMaxExpressions.csv', row.names=TRUE)

Fold_meta_tissues <- merge(FoldExpressionValues, Fold_meta, by.x='Entrez_Gene_ID',
                           by.y='Entrez_Gene_ID')
write.csv(Fold_meta_tissues, 'Fold_meta_tissues_mean.csv', row.names=TRUE)

################################################################################################

Fold_median <- read.csv('Fold_median_Allmeta.csv', sep=',', header=TRUE)

# Extract which human tissue sees the most expression or inhibition of these genes from
# NCBIgene.org, or UCSC, or ENSEMBLE, or genenames.org
# go to each symbol gene ID in https://www.ncbi.nlm.nih.gov/gene search the symbol, and
# go to 'expression' and 'summary' sections to see a summary and bar chart of tissues the 
# gene is seen in

#you can download the expression values, by ENTREZ gene ID (numeric ID) for each

Fold <- unique(Fold_median$Entrez_Gene_ID) #look each up at https://www.ncbi.nlm.nih.gov/gene/
# for gene expression HPA RNA-Seq normal tissues

# 8312 127343   2893   8369(NA-histone)   3605 386618   3886   3950 283551 407835(NA-provisional)
# 128360(NA-olfactory) 401667(NA-olfactory) 119678(NA-olfactory) 
# 171558 137868 346673 126637  84076(NA-validated)  7743 127665
Fold <- Fold[-c(4,10,11,12,13,18)]
# these files were individually downloaded from expression values in normal tissue of 95 humans
# in 27 tissues beginning with 'GeneID' no separator and ID value i.e. 'GeneID70'

FoldList <- as.character(Fold)
foldpaste <- paste('GeneID', FoldList[1:14], '.txt', sep='')

# Switch to subdirectory that you stored these downloaded files
list.files()
# [1] "GeneID126637.txt" "GeneID127343.txt" "GeneID127665.txt" "GeneID137868.txt" "GeneID171558.txt"
# [6] "GeneID283551.txt" "GeneID2893.txt"   "GeneID346673.txt" "GeneID3605.txt"   "GeneID386618.txt"
# [11] "GeneID3886.txt"   "GeneID3950.txt"   "GeneID7743.txt"   "GeneID8312.txt" 


d1 <- read.delim(foldpaste[1], sep='\t', quote="", header=TRUE,
                 skip=1, na.strings=c('','NA'))
d2 <- read.delim(foldpaste[2], sep='\t', quote="", header=TRUE,
                 skip=1, na.strings=c('','NA'))
d3 <- read.delim(foldpaste[3], sep='\t', quote="", header=TRUE,
                 skip=1, na.strings=c('','NA'))
d4 <- read.delim(foldpaste[4], sep='\t', quote="", header=TRUE,
                 skip=1, na.strings=c('','NA'))
d5 <- read.delim(foldpaste[5], sep='\t', quote="", header=TRUE,
                 skip=1, na.strings=c('','NA'))
d6 <- read.delim(foldpaste[6], sep='\t', quote="", header=TRUE,
                 skip=1, na.strings=c('','NA'))
d7 <- read.delim(foldpaste[7], sep='\t', quote="", header=TRUE,
                 skip=1, na.strings=c('','NA'))
d8 <- read.delim(foldpaste[8], sep='\t', quote="", header=TRUE,
                 skip=1, na.strings=c('','NA'))
d9 <- read.delim(foldpaste[9], sep='\t', quote="", header=TRUE,
                 skip=1, na.strings=c('','NA'))
d10 <- read.delim(foldpaste[10], sep='\t', quote="", header=TRUE,
                  skip=1, na.strings=c('','NA'))
d11 <- read.delim(foldpaste[11], sep='\t', quote="", header=TRUE,
                  skip=1, na.strings=c('','NA'))
d12 <- read.delim(foldpaste[12], sep='\t', quote="", header=TRUE,
                  skip=1, na.strings=c('','NA'))
d13 <- read.delim(foldpaste[13], sep='\t', quote="", header=TRUE,
                  skip=1, na.strings=c('','NA'))
d14 <- read.delim(foldpaste[14], sep='\t', quote="", header=TRUE,
                  skip=1, na.strings=c('','NA'))

# Switch back to parent directory of all other files

Fold_Entrez_df <- rbind(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,
                        d13,d14)

colnames(Fold_Entrez_df)[1] <- 'Entrez_Gene_ID'
Fold_Entrez_df <- Fold_Entrez_df[,-29]
write.csv(Fold_Entrez_df,'Entrez_Gene_Expressions_Fold_median_14.csv',row.names=FALSE)

ExpressionsFold <- read.csv('Entrez_Gene_Expressions_Fold_median_14.csv', sep=',', header=TRUE)

Fold_meta <- read.csv('Fold_median_Allmeta.csv', sep=',', header=TRUE)

row.names(ExpressionsFold) <- ExpressionsFold$Entrez_Gene_ID
Expressions <- ExpressionsFold[,-1]
names <- row.names(Expressions)

normalMax <- apply(Expressions, 1, max)

library(dplyr)

FoldExpressions <- mutate(Expressions, NormalMax=normalMax)
row.names(FoldExpressions) <- names
FoldExpressions <- round(FoldExpressions,2)

t1 <- order(FoldExpressions[1,], decreasing=TRUE)[1]
t2 <- order(FoldExpressions[2,], decreasing=TRUE)[1]
t3 <- order(FoldExpressions[3,], decreasing=TRUE)[1]
t4 <- order(FoldExpressions[4,], decreasing=TRUE)[1]
t5 <- order(FoldExpressions[5,], decreasing=TRUE)[1]
t6 <- order(FoldExpressions[6,], decreasing=TRUE)[1]
t7 <- order(FoldExpressions[7,], decreasing=TRUE)[1]
t8 <- order(FoldExpressions[8,], decreasing=TRUE)[1]
t9 <- order(FoldExpressions[9,], decreasing=TRUE)[1]
t10 <- order(FoldExpressions[10,], decreasing=TRUE)[1]
t11 <- order(FoldExpressions[11,], decreasing=TRUE)[1]
t12 <- order(FoldExpressions[12,], decreasing=TRUE)[1]
t13 <- order(FoldExpressions[13,], decreasing=TRUE)[1]
t14 <- order(FoldExpressions[14,], decreasing=TRUE)[1]
# t15 <- order(FoldExpressions[15,], decreasing=TRUE)[1]
# t16 <- order(FoldExpressions[16,], decreasing=TRUE)[1]
# t17 <- order(FoldExpressions[17,], decreasing=TRUE)[1]
# t18 <- order(FoldExpressions[18,], decreasing=TRUE)[1]
# t19 <- order(FoldExpressions[19,], decreasing=TRUE)[1]
# t20 <- order(FoldExpressions[20,], decreasing=TRUE)[1]

TissueMax <- as.vector(c(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14))
Tissue <- colnames(FoldExpressions)[c(TissueMax)]

FoldExpressionValues <- mutate(FoldExpressions, TissueMax=Tissue)
row.names(FoldExpressionValues) <- row.names(FoldExpressions)
FoldExpressionValues$Entrez_Gene_ID <- row.names(FoldExpressionValues)

write.csv(FoldExpressionValues, 'Fold_median_ValuesTissueMaxExpressions.csv', row.names=TRUE)

Fold_meta_tissues <- merge(FoldExpressionValues, Fold_meta, by.x='Entrez_Gene_ID',
                           by.y='Entrez_Gene_ID')
write.csv(Fold_meta_tissues, 'Fold_meta_tissues_median.csv', row.names=TRUE)

################################################################################################

Magnitude <- read.csv('Magnitude_Allmeta.csv', header=TRUE, sep=',')

# Extract which human tissue sees the most expression or inhibition of these genes from
# NCBIgene.org, or UCSC, or ENSEMBLE, or genenames.org
# go to each symbol gene ID in https://www.ncbi.nlm.nih.gov/gene search the symbol, and
# go to 'expression' and 'summary' sections to see a summary and bar chart of tissues the 
# gene is seen in

#you can download the expression values, by ENTREZ gene ID (numeric ID) for each

Fold <- unique(Magnitude$Entrez_Gene_ID) #look each up at https://www.ncbi.nlm.nih.gov/gene/
# for gene expression HPA RNA-Seq normal tissues

# 59     60     72    975   1490   1634   1843   1915   3315   4240   4256   
# 4629  284119   6135   6169   6210   6230   6194   7117(NA-pseudo)   7178

# these files were individually downloaded from expression values in normal tissue of 95 humans
# in 27 tissues beginning with 'GeneID' no separator and ID value i.e. 'GeneID70'
Fold <- Fold[-19]
FoldList <- as.character(Fold)
foldpaste <- paste('GeneID', FoldList[1:19], '.txt', sep='')

# Switch to subdirectory you stored the downloaded files or wherever you put them in own folder
list.files()
# [1] "GeneID1490.txt"   "GeneID1634.txt"   "GeneID1843.txt"   "GeneID1915.txt"   "GeneID284119.txt"
# [6] "GeneID3315.txt"   "GeneID4240.txt"   "GeneID4256.txt"   "GeneID4629.txt"   "GeneID59.txt"    
# [11] "GeneID60.txt"     "GeneID6135.txt"   "GeneID6169.txt"   "GeneID6194.txt"   "GeneID6210.txt"  
# [16] "GeneID6230.txt"   "GeneID7178.txt"   "GeneID72.txt"     "GeneID975.txt"

d1 <- read.delim(foldpaste[1], sep='\t', quote="", header=TRUE,
                 skip=1, na.strings=c('','NA'))
d2 <- read.delim(foldpaste[2], sep='\t', quote="", header=TRUE,
                 skip=1, na.strings=c('','NA'))
d3 <- read.delim(foldpaste[3], sep='\t', quote="", header=TRUE,
                 skip=1, na.strings=c('','NA'))
d4 <- read.delim(foldpaste[4], sep='\t', quote="", header=TRUE,
                 skip=1, na.strings=c('','NA'))
d5 <- read.delim(foldpaste[5], sep='\t', quote="", header=TRUE,
                 skip=1, na.strings=c('','NA'))
d6 <- read.delim(foldpaste[6], sep='\t', quote="", header=TRUE,
                 skip=1, na.strings=c('','NA'))
d7 <- read.delim(foldpaste[7], sep='\t', quote="", header=TRUE,
                 skip=1, na.strings=c('','NA'))
d8 <- read.delim(foldpaste[8], sep='\t', quote="", header=TRUE,
                 skip=1, na.strings=c('','NA'))
d9 <- read.delim(foldpaste[9], sep='\t', quote="", header=TRUE,
                 skip=1, na.strings=c('','NA'))
d10 <- read.delim(foldpaste[10], sep='\t', quote="", header=TRUE,
                  skip=1, na.strings=c('','NA'))
d11 <- read.delim(foldpaste[11], sep='\t', quote="", header=TRUE,
                  skip=1, na.strings=c('','NA'))
d12 <- read.delim(foldpaste[12], sep='\t', quote="", header=TRUE,
                  skip=1, na.strings=c('','NA'))
d13 <- read.delim(foldpaste[13], sep='\t', quote="", header=TRUE,
                  skip=1, na.strings=c('','NA'))
d14 <- read.delim(foldpaste[14], sep='\t', quote="", header=TRUE,
                  skip=1, na.strings=c('','NA'))
d15 <- read.delim(foldpaste[15], sep='\t', quote="", header=TRUE,
                  skip=1, na.strings=c('','NA'))
d16 <- read.delim(foldpaste[16], sep='\t', quote="", header=TRUE,
                  skip=1, na.strings=c('','NA'))
d17 <- read.delim(foldpaste[17], sep='\t', quote="", header=TRUE,
                  skip=1, na.strings=c('','NA'))
d18 <- read.delim(foldpaste[18], sep='\t', quote="", header=TRUE,
                  skip=1, na.strings=c('','NA'))
d19 <- read.delim(foldpaste[19], sep='\t', quote="", header=TRUE,
                  skip=1, na.strings=c('','NA'))
# d20 <- read.delim(foldpaste[20], sep='\t', quote="", header=TRUE,
#                   skip=1, na.strings=c('','NA'))
Fold_Entrez_df <- rbind(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,
                        d13,d14,d15,d16,d17,d18,d19)

# Switch back to parent folder of all other files

colnames(Fold_Entrez_df)[1] <- 'Entrez_Gene_ID'
Fold_Entrez_df <- Fold_Entrez_df[,-29]
write.csv(Fold_Entrez_df,'Entrez_Gene_Expressions_magnitude19.csv',row.names=FALSE)

ExpressionsFold <- read.csv('Entrez_Gene_Expressions_magnitude19.csv', sep=',', header=TRUE)

Fold_meta <- read.csv('Magnitude_Allmeta.csv', sep=',', header=TRUE)

row.names(ExpressionsFold) <- ExpressionsFold$Entrez_Gene_ID
Expressions <- ExpressionsFold[,-1]
names <- row.names(Expressions)

normalMax <- apply(Expressions, 1, max)

library(dplyr)

FoldExpressions <- mutate(Expressions, NormalMax=normalMax)
row.names(FoldExpressions) <- names
FoldExpressions <- round(FoldExpressions,2)

t1 <- order(FoldExpressions[1,], decreasing=TRUE)[1]
t2 <- order(FoldExpressions[2,], decreasing=TRUE)[1]
t3 <- order(FoldExpressions[3,], decreasing=TRUE)[1]
t4 <- order(FoldExpressions[4,], decreasing=TRUE)[1]
t5 <- order(FoldExpressions[5,], decreasing=TRUE)[1]
t6 <- order(FoldExpressions[6,], decreasing=TRUE)[1]
t7 <- order(FoldExpressions[7,], decreasing=TRUE)[1]
t8 <- order(FoldExpressions[8,], decreasing=TRUE)[1]
t9 <- order(FoldExpressions[9,], decreasing=TRUE)[1]
t10 <- order(FoldExpressions[10,], decreasing=TRUE)[1]
t11 <- order(FoldExpressions[11,], decreasing=TRUE)[1]
t12 <- order(FoldExpressions[12,], decreasing=TRUE)[1]
t13 <- order(FoldExpressions[13,], decreasing=TRUE)[1]
t14 <- order(FoldExpressions[14,], decreasing=TRUE)[1]
t15 <- order(FoldExpressions[15,], decreasing=TRUE)[1]
t16 <- order(FoldExpressions[16,], decreasing=TRUE)[1]
t17 <- order(FoldExpressions[17,], decreasing=TRUE)[1]
t18 <- order(FoldExpressions[18,], decreasing=TRUE)[1]
t19 <- order(FoldExpressions[19,], decreasing=TRUE)[1]
#t20 <- order(FoldExpressions[20,], decreasing=TRUE)[1]

TissueMax <- as.vector(c(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19))
Tissue <- colnames(FoldExpressions)[c(TissueMax)]

FoldExpressionValues <- mutate(FoldExpressions, TissueMax=Tissue)
row.names(FoldExpressionValues) <- row.names(FoldExpressions)
FoldExpressionValues$Entrez_Gene_ID <- row.names(FoldExpressionValues)

write.csv(FoldExpressionValues, 'MagnitudeMeanValuesTissueMaxExpressions.csv', row.names=TRUE)

Fold_meta_tissues <- merge(FoldExpressionValues, Fold_meta, by.x='Entrez_Gene_ID',
                           by.y='Entrez_Gene_ID')
write.csv(Fold_meta_tissues, 'MagnitudeMean_meta_tissues.csv', row.names=TRUE)

################################################################################################

Magnitude_median <- read.csv('Magnitude_Median_Allmeta.csv', sep=',',
                             header=TRUE)

# Extract which human tissue sees the most expression or inhibition of these genes from
# NCBIgene.org, or UCSC, or ENSEMBLE, or genenames.org
# go to each symbol gene ID in https://www.ncbi.nlm.nih.gov/gene search the symbol, and
# go to 'expression' and 'summary' sections to see a summary and bar chart of tissues the 
# gene is seen in

#you can download the expression values, by ENTREZ gene ID (numeric ID) for each

Fold <- unique(Magnitude_median$Entrez_Gene_ID) #look each up at 
# https://www.ncbi.nlm.nih.gov/gene/
# for gene expression HPA RNA-Seq normal tissues

# 59    60    72   975  1915 10399  4629  6169  6125  6204  6205  
# 6210  6217  6230  6235  6194 6876  7117(NA-pseudo)  7169  7178

Fold <- Fold[-18]

# these files were individually downloaded from expression values in normal tissue of 95 humans
# in 27 tissues beginning with 'GeneID' no separator and ID value i.e. 'GeneID70'

FoldList <- as.character(Fold)
foldpaste <- paste('GeneID', FoldList[1:20], '.txt', sep='')

# switch to folder you downloaded these files:
list.files()
# [1] "GeneID10399.txt" "GeneID1915.txt"  "GeneID4629.txt"  "GeneID59.txt"    "GeneID60.txt"   
# [6] "GeneID6125.txt"  "GeneID6169.txt"  "GeneID6194.txt"  "GeneID6204.txt"  "GeneID6205.txt" 
# [11] "GeneID6210.txt"  "GeneID6217.txt"  "GeneID6230.txt"  "GeneID6235.txt"  "GeneID6876.txt" 
# [16] "GeneID7169.txt"  "GeneID7178.txt"  "GeneID72.txt"    "GeneID975.txt" 

d1 <- read.delim(foldpaste[1], sep='\t', quote="", header=TRUE,
                 skip=1, na.strings=c('','NA'))
d2 <- read.delim(foldpaste[2], sep='\t', quote="", header=TRUE,
                 skip=1, na.strings=c('','NA'))
d3 <- read.delim(foldpaste[3], sep='\t', quote="", header=TRUE,
                 skip=1, na.strings=c('','NA'))
d4 <- read.delim(foldpaste[4], sep='\t', quote="", header=TRUE,
                 skip=1, na.strings=c('','NA'))
d5 <- read.delim(foldpaste[5], sep='\t', quote="", header=TRUE,
                 skip=1, na.strings=c('','NA'))
d6 <- read.delim(foldpaste[6], sep='\t', quote="", header=TRUE,
                 skip=1, na.strings=c('','NA'))
d7 <- read.delim(foldpaste[7], sep='\t', quote="", header=TRUE,
                 skip=1, na.strings=c('','NA'))
d8 <- read.delim(foldpaste[8], sep='\t', quote="", header=TRUE,
                 skip=1, na.strings=c('','NA'))
d9 <- read.delim(foldpaste[9], sep='\t', quote="", header=TRUE,
                 skip=1, na.strings=c('','NA'))
d10 <- read.delim(foldpaste[10], sep='\t', quote="", header=TRUE,
                  skip=1, na.strings=c('','NA'))
d11 <- read.delim(foldpaste[11], sep='\t', quote="", header=TRUE,
                  skip=1, na.strings=c('','NA'))
d12 <- read.delim(foldpaste[12], sep='\t', quote="", header=TRUE,
                  skip=1, na.strings=c('','NA'))
d13 <- read.delim(foldpaste[13], sep='\t', quote="", header=TRUE,
                  skip=1, na.strings=c('','NA'))
d14 <- read.delim(foldpaste[14], sep='\t', quote="", header=TRUE,
                  skip=1, na.strings=c('','NA'))
d15 <- read.delim(foldpaste[15], sep='\t', quote="", header=TRUE,
                  skip=1, na.strings=c('','NA'))
d16 <- read.delim(foldpaste[16], sep='\t', quote="", header=TRUE,
                  skip=1, na.strings=c('','NA'))
d17 <- read.delim(foldpaste[17], sep='\t', quote="", header=TRUE,
                  skip=1, na.strings=c('','NA'))
d18 <- read.delim(foldpaste[18], sep='\t', quote="", header=TRUE,
                  skip=1, na.strings=c('','NA'))
d19 <- read.delim(foldpaste[19], sep='\t', quote="", header=TRUE,
                  skip=1, na.strings=c('','NA'))
# d20 <- read.delim(foldpaste[20], sep='\t', quote="", header=TRUE,
#                   skip=1, na.strings=c('','NA'))

#switch to directory of all other files

Fold_Entrez_df <- rbind(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,
                        d13,d14,d15,d16,d17,d18,d19)

colnames(Fold_Entrez_df)[1] <- 'Entrez_Gene_ID'
Fold_Entrez_df <- Fold_Entrez_df[,-29]
write.csv(Fold_Entrez_df,'Entrez_Gene_Expressions_Magnitude_Median_19.csv',row.names=FALSE)

ExpressionsFold <- read.csv('Entrez_Gene_Expressions_Magnitude_Median_19.csv', sep=',', header=TRUE)

Fold_meta <- read.csv('Magnitude_Median_Allmeta.csv', sep=',', header=TRUE)

row.names(ExpressionsFold) <- ExpressionsFold$Entrez_Gene_ID
Expressions <- ExpressionsFold[,-1]
names <- row.names(Expressions)

normalMax <- apply(Expressions, 1, max)

library(dplyr)

FoldExpressions <- mutate(Expressions, NormalMax=normalMax)
row.names(FoldExpressions) <- names
FoldExpressions <- round(FoldExpressions,2)

t1 <- order(FoldExpressions[1,], decreasing=TRUE)[1]
t2 <- order(FoldExpressions[2,], decreasing=TRUE)[1]
t3 <- order(FoldExpressions[3,], decreasing=TRUE)[1]
t4 <- order(FoldExpressions[4,], decreasing=TRUE)[1]
t5 <- order(FoldExpressions[5,], decreasing=TRUE)[1]
t6 <- order(FoldExpressions[6,], decreasing=TRUE)[1]
t7 <- order(FoldExpressions[7,], decreasing=TRUE)[1]
t8 <- order(FoldExpressions[8,], decreasing=TRUE)[1]
t9 <- order(FoldExpressions[9,], decreasing=TRUE)[1]
t10 <- order(FoldExpressions[10,], decreasing=TRUE)[1]
t11 <- order(FoldExpressions[11,], decreasing=TRUE)[1]
t12 <- order(FoldExpressions[12,], decreasing=TRUE)[1]
t13 <- order(FoldExpressions[13,], decreasing=TRUE)[1]
t14 <- order(FoldExpressions[14,], decreasing=TRUE)[1]
t15 <- order(FoldExpressions[15,], decreasing=TRUE)[1]
t16 <- order(FoldExpressions[16,], decreasing=TRUE)[1]
t17 <- order(FoldExpressions[17,], decreasing=TRUE)[1]
t18 <- order(FoldExpressions[18,], decreasing=TRUE)[1]
t19 <- order(FoldExpressions[19,], decreasing=TRUE)[1]
# t20 <- order(FoldExpressions[20,], decreasing=TRUE)[1]

TissueMax <- as.vector(c(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19))
Tissue <- colnames(FoldExpressions)[c(TissueMax)]

FoldExpressionValues <- mutate(FoldExpressions, TissueMax=Tissue)
row.names(FoldExpressionValues) <- row.names(FoldExpressions)
FoldExpressionValues$Entrez_Gene_ID <- row.names(FoldExpressionValues)

write.csv(FoldExpressionValues, 'Magnitude_Median_ValuesTissueMaxExpressions.csv', row.names=TRUE)

Fold_meta_tissues <- merge(FoldExpressionValues, Fold_meta, by.x='Entrez_Gene_ID',
                           by.y='Entrez_Gene_ID')
write.csv(Fold_meta_tissues, 'Magnitude_meta_tissues_median.csv', row.names=TRUE)

################################################################################################
# there are no Entrez_Gene_ID names in common between the fold change for mean and median values
# there are some genes in common among the mean and median values for magnitude of DE:
# 59,60,72, 975,1915,4629,6169,6194,6210,6230,7178
DE_median <- read.csv('Magnitude_meta_tissues_median.csv', sep=',', header=TRUE, row.names=1)
DE_mean <- read.csv('MagnitudeMean_meta_tissues.csv', sep=',', header=TRUE, row.names=1)
FC_median <- read.csv('Fold_meta_tissues_median.csv', sep=',', header=TRUE, row.names=1)
FC_mean <- read.csv('Fold_meta_tissues_mean.csv', sep=',', header=TRUE, row.names=1)

DE_mean_2 <- DE_mean[,c(31,30,37,38,44,55,56,57:60,62:133)]
DE_MEAN <- DE_mean_2[!duplicated(DE_mean_2),]

DE_median_2 <- DE_median[,c(31,30,40,41,44,55,56,57:60,62:133)]
DE_MEDIAN <- DE_median_2[!duplicated(DE_median_2),]

FC_mean_2 <- FC_mean[,c(31,30,36,37,44,55,56,57:60,62:133)]
FC_MEAN <- FC_mean_2[!duplicated(FC_mean_2),]

FC_median_2 <- FC_median[,c(31,30,39,40,44,55,56,57:60,62:133)]
FC_MEDIAN <- FC_median_2[!duplicated(FC_median_2),]

# for researching the genes further and involvement in Fold Change and DE 
# in median and mean samples
write.csv(DE_MEAN, 'DE_MEAN_21_83.csv', row.names=FALSE)
write.csv(DE_MEDIAN, 'DE_MEDIAN_21_83.csv', row.names=FALSE)
write.csv(FC_MEAN, 'FC_MEAN_20_83.csv', row.names=FALSE)
write.csv(FC_MEDIAN, 'FC_MEDIAN_15_83.csv', row.names=FALSE)

####################################################################################################
####################################################################################################
####################################################################################################

# How about testing the genes that have a regulation in Vitamin D, iron, and alchol catabolization?
# The first pick of listed genes at ncbi.nlm.nih.gov
# Vitamin D genes: VDR (receptor) and GC (binding)
# Alcohol dehydrogenase genes: indicators for alcohol presence, ADH1B, ADH1C, ADH4
# Iron genes: TF, HFR, SLC40A1
# glucose genes:G6PC2, SLC2A1, TCF4, GCK
# glycogen genes: 

# alcohol protein genes:
ADH1B <- StatsAll[grep('ADH1B', StatsAll$Symbol),]
ADH1C <- StatsAll[grep('ADH1C', StatsAll$Symbol),]
ADH4 <- StatsAll[grep('ADH4', StatsAll$Symbol),]

# vitamin D genes:
GC <- StatsAll[grep('GC', StatsAll$Symbol),]
GC <- GC[8,]
VDR <- StatsAll[grep('VDR', StatsAll$Symbol),]

# iron genes:
TF <- StatsAll[grep('TF', StatsAll$Symbol),]
TF <- TF[68,]
# HFR <- StatsAll[grep('HFR', StatsAll$Symbol),] #none are HFR

# Glycogen genes:
PYGM <- StatsAll[grep('PYGM', StatsAll$Symbol),]
AGL <- StatsAll[grep('AGL', StatsAll$Symbol),]
AGL <- AGL[1,]
GSK3B <- StatsAll[grep('GSK3B', StatsAll$Symbol),]
GYS1 <- StatsAll[grep('GYS1', StatsAll$Symbol),]

# Glucose genes:
G6PC2 <- StatsAll[grep('G6PC2', StatsAll$Symbol),]
SLC2A1 <- StatsAll[grep('SLC2A1', StatsAll$Symbol),]
SLC2A1 <- SLC2A1[1,]
TCF4 <- StatsAll[grep('TCF4', StatsAll$Symbol),]
GCK <- StatsAll[grep('GCK', StatsAll$Symbol),]
GCK <- GCK[1,]

# data frame of these genes involved/mentioned in tumorigenesis and UL pathogenesis research
OtherGenes <- rbind(ADH4,ADH1B,ADH1C,AGL,PYGM,GSK3B,GYS1,G6PC2,SLC2A1,TCF4,GCK, 
                    GC, VDR, TF)
write.csv(OtherGenes, 'otherULgenes.csv', row.names=TRUE)

OtherGenes <- read.csv('otherULgenes.csv', sep=',', header=TRUE, row.names=1,
                       na.strings=c('','NA'))

StatsAll <- read.csv('StatsAllCompressed.csv', sep=',', header=TRUE, row.names=1)
DE_MEAN <- read.csv('DE_MEAN_21_83.csv', sep=',', header=TRUE)
DE_MEDIAN <- read.csv('DE_MEDIAN_21_83.csv', sep=',', header=TRUE)
FC_MEAN <- read.csv('FC_MEAN_20_83.csv', sep=',', header=TRUE)
FC_MEDIAN <-  read.csv('FC_MEDIAN_15_83.csv', sep=',', header=TRUE)

table_less <- read.csv('table_less.csv', header=TRUE, na.strings=c('', 'NA','null'),)
table_more <- read.csv('table_more.csv', header=TRUE, sep=',',
                       na.strings=c('','null','NA'))
meta <- table_more[,1:22]

Symbol <- row.names(OtherGenes)
OtherGenes$Symbol <- Symbol

OtherDF <- merge(OtherGenes, meta, by.x='Symbol', by.y='Symbol')

DF_others <- OtherDF[,c(1:11,84:104,12:83)]

DF_others_sm <- DF_others[,-c(12,13:19,21,23,24)]#keep Entrez_Gene_ID

DF_others_sm <- DF_others_sm[!duplicated(DF_others_sm),]

write.csv(DF_others_sm, 'Other_Genes_meta_stats.csv', row.names=FALSE)

# Extract which human tissue sees the most expression or inhibition of these genes from
# NCBIgene.org, or UCSC, or ENSEMBLE, or genenames.org
# go to each symbol gene ID in https://www.ncbi.nlm.nih.gov/gene search the symbol, and
# go to 'expression' and 'summary' sections to see a summary and bar chart of tissues the 
# gene is seen in

#you can download the expression values, by ENTREZ gene ID (numeric ID) for each

ID <- unique(DF_others_sm$Entrez_Gene_ID) #look each up at 
# https://www.ncbi.nlm.nih.gov/gene/
# for gene expression HPA RNA-Seq normal tissues

# 125   126   127   178 57818  2638  2645  2932  2997  5837  6513  6925  7018  7421

# these files were individually downloaded from expression values in normal tissue of 95 humans
# in 27 tissues beginning with 'GeneID' no separator and ID value i.e. 'GeneID70'

FoldList <- as.character(ID)
foldpaste <- paste('GeneID', FoldList[1:20], '.txt', sep='')
#switch to directory you downloaded these files

list.files()#the names as is from download directly to subfolder, there are 14 files
# [1] "GeneID125.txt"   "GeneID126.txt"   "GeneID127.txt"   "GeneID178.txt"   "GeneID2638.txt" 
# [6] "GeneID2645.txt"  "GeneID2932.txt"  "GeneID2997.txt"  "GeneID57818.txt" "GeneID5837.txt" 
# [11] "GeneID6513.txt"  "GeneID6925.txt"  "GeneID7018.txt"  "GeneID7421.txt"

d1 <- read.delim(foldpaste[1], sep='\t', quote="", header=TRUE,
                 skip=1, na.strings=c('','NA'))
d2 <- read.delim(foldpaste[2], sep='\t', quote="", header=TRUE,
                 skip=1, na.strings=c('','NA'))
d3 <- read.delim(foldpaste[3], sep='\t', quote="", header=TRUE,
                 skip=1, na.strings=c('','NA'))
d4 <- read.delim(foldpaste[4], sep='\t', quote="", header=TRUE,
                 skip=1, na.strings=c('','NA'))
d5 <- read.delim(foldpaste[5], sep='\t', quote="", header=TRUE,
                 skip=1, na.strings=c('','NA'))
d6 <- read.delim(foldpaste[6], sep='\t', quote="", header=TRUE,
                 skip=1, na.strings=c('','NA'))
d7 <- read.delim(foldpaste[7], sep='\t', quote="", header=TRUE,
                 skip=1, na.strings=c('','NA'))
d8 <- read.delim(foldpaste[8], sep='\t', quote="", header=TRUE,
                 skip=1, na.strings=c('','NA'))
d9 <- read.delim(foldpaste[9], sep='\t', quote="", header=TRUE,
                 skip=1, na.strings=c('','NA'))
d10 <- read.delim(foldpaste[10], sep='\t', quote="", header=TRUE,
                  skip=1, na.strings=c('','NA'))
d11 <- read.delim(foldpaste[11], sep='\t', quote="", header=TRUE,
                  skip=1, na.strings=c('','NA'))
d12 <- read.delim(foldpaste[12], sep='\t', quote="", header=TRUE,
                  skip=1, na.strings=c('','NA'))
d13 <- read.delim(foldpaste[13], sep='\t', quote="", header=TRUE,
                  skip=1, na.strings=c('','NA'))
d14 <- read.delim(foldpaste[14], sep='\t', quote="", header=TRUE,
                  skip=1, na.strings=c('','NA'))

# switch back to folder of all other files 

Fold_Entrez_df <- rbind(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,
                        d13,d14)

colnames(Fold_Entrez_df)[1] <- 'Entrez_Gene_ID'
Fold_Entrez_df <- Fold_Entrez_df[,-29]
write.csv(Fold_Entrez_df,'Entrez_Gene_Expressions_Other_Genes_14.csv',row.names=FALSE)

ExpressionsOthers <- read.csv('Entrez_Gene_Expressions_Other_Genes_14.csv', sep=',',
                              header=TRUE)

row.names(ExpressionsOthers) <- ExpressionsOthers$Entrez_Gene_ID
Expressions <- ExpressionsOthers[,-1]
names <- row.names(Expressions)

normalMax <- apply(Expressions, 1, max)

library(dplyr)

OthersExpressions <- mutate(Expressions, NormalMax=normalMax)
row.names(OthersExpressions) <- names
OthersExpressions <- round(OthersExpressions,2)

t1 <- order(OthersExpressions[1,], decreasing=TRUE)[1]
t2 <- order(OthersExpressions[2,], decreasing=TRUE)[1]
t3 <- order(OthersExpressions[3,], decreasing=TRUE)[1]
t4 <- order(OthersExpressions[4,], decreasing=TRUE)[1]
t5 <- order(OthersExpressions[5,], decreasing=TRUE)[1]
t6 <- order(OthersExpressions[6,], decreasing=TRUE)[1]
t7 <- order(OthersExpressions[7,], decreasing=TRUE)[1]
t8 <- order(OthersExpressions[8,], decreasing=TRUE)[1]
t9 <- order(OthersExpressions[9,], decreasing=TRUE)[1]
t10 <- order(OthersExpressions[10,], decreasing=TRUE)[1]
t11 <- order(OthersExpressions[11,], decreasing=TRUE)[1]
t12 <- order(OthersExpressions[12,], decreasing=TRUE)[1]
t13 <- order(OthersExpressions[13,], decreasing=TRUE)[1]
t14 <- order(OthersExpressions[14,], decreasing=TRUE)[1]

TissueMax <- as.vector(c(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14))
Tissue <- colnames(OthersExpressions)[c(TissueMax)]

OthersExpressionValues <- mutate(OthersExpressions, TissueMax=Tissue)
row.names(OthersExpressionValues) <- row.names(OthersExpressions)
OthersExpressionValues$Entrez_Gene_ID <- row.names(OthersExpressionValues)

write.csv(OthersExpressionValues, 'OtherGenesExpression_ValuesTissueMax.csv', row.names=TRUE)

Others_meta_tissues <- merge(OthersExpressionValues, DF_others_sm, by.x='Entrez_Gene_ID',
                           by.y='Entrez_Gene_ID')
write.csv(Others_meta_tissues, 'OthersGenesALLMetaTissueExpressions.csv', row.names=TRUE)

ml_others <- Others_meta_tissues[,c(31,51:122)]
ml_others <- ml_others[!duplicated(ml_others),]
row.names(ml_others) <- ml_others$Symbol

ml_others <- ml_others[,-1]
names <- row.names(ml_others)
names2 <- colnames(ml_others)

ml_others_t <- t(ml_others)

ul <- data.frame(rep('UL',43))
colnames(ul) <- 'TYPE'
non <- data.frame(rep('nonUL',29))
colnames(non) <- 'TYPE'
TYPE <- rbind(non,ul)
row.names(TYPE) <- names2

ml_others_ready <- cbind(TYPE, ml_others_t)
write.csv(ml_others_ready, 'otherGenesMLready.csv', row.names=TRUE)
##############################################################################################

Others <- read.csv('otherGenesMLready.csv', sep=',', header=TRUE, row.names=1)

library(caret)
library(randomForest)
library(MASS)
library(gbm)
library(e1071)

set.seed(189678345)

# creates a partition of the data by indices
inTrain <- createDataPartition(y=Others$TYPE, p=0.7, list=FALSE)


trainingSet <- Others[inTrain,]#52X15
testingSet <- Others[-inTrain,]#20X15

rfMod <- train(TYPE~., method='rf', data=(trainingSet), 
               trControl=trainControl(method='cv'), number=5)

gbmMod <- train(TYPE~., method='gbm', data=trainingSet, verbose=FALSE )
ldaMod <- train(TYPE~., method='lda', data=trainingSet)

predRF <- predict(rfMod, testingSet)
predGbm <- predict(gbmMod, testingSet)
predlda <- predict(ldaMod, testingSet)

predDF <- data.frame(predRF, predGbm, predlda, type=testingSet$TYPE)


CombinedModels <- train(type~., method='gam', data=predDF)
CombinedPredictions <- predict(CombinedModels, predDF)

sum <- sum(CombinedPredictions==testingSet$TYPE)
length <- length(CombinedPredictions)

sum <- sum(predRF==testingSet$TYPE) 
length <- length(predRF)
accuracy_rfMOd <- (sum/length) 

sum <- sum(predGbm==testingSet$TYPE)
accuracy_gbmMod <- sum/length 


sum <- sum(predlda==testingSet$TYPE) 
accuracy_ldaMod <- sum/length 

rf2 <- randomForest(TYPE~., data=trainingSet, method='class')
predRF2 <- predict(rf2, testingSet, type='class')

sum <- sum(predRF2==testingSet$TYPE)
accuracy_RF2 <- sum/length 

confusionMatrix(predRF2, testingSet$TYPE)

predDF2 <- data.frame(predRF,predRF2,predlda, predGbm, testingSet$TYPE)
colnames(predDF2)[5] <- 'TYPE'

CombinedModels <- train(TYPE ~ ., method='gam', data=predDF2)
CombinedPredictions2 <- predict(CombinedModels, predDF)
sum <- sum(CombinedPredictions2==testingSet$TYPE)
length(CombinedPredictions2)
accuracy_CombinedPredictions2 <- sum/length 

predDF3 <- data.frame(predRF,predRF2,predlda, predGbm, 
                      CombinedPredictions2, testingSet$TYPE)
colnames(predDF3)[6] <- 'TYPE'

knnMod <- train(TYPE ~ .,
                method='knn', preProcess=c('center','scale'),
                tuneLength=10, trControl=trainControl(method='cv'), 
                data=trainingSet)
plot(knnMod)

rpartMod <- train(TYPE ~ ., method='rpart', tuneLength=9, data=trainingSet) 
plot(rpart)

glmMod <- train(TYPE ~ ., 
                method='glm', data=trainingSet) 

predKNN <- predict(knnMod, testingSet)
predRPART <- predict(rpartMod, testingSet)
predGLM <- predict(glmMod, testingSet)



length=length(testingSet$TYPE)

sumKNN <- sum(predKNN==testingSet$TYPE)
sumRPart <- sum(predRPART==testingSet$TYPE)
sumGLM <- sum(predGLM==testingSet$TYPE)

accuracy_KNN <- sumKNN/length 
accuracy_RPART <- sumRPart/length 
accuracy_GLM <- sumGLM/length 

predDF3 <- data.frame(predDF2[,1:4],predKNN,predRPART,predGLM, 
                      TYPE=testingSet$TYPE)

CombinedModels <- train(TYPE ~ ., method='gam', data=predDF3)
CombinedPredictions2 <- predict(CombinedModels, predDF3)
accuracy_CP2 <- sum(CombinedPredictions2==testingSet$TYPE)/length 

predDF4 <- data.frame(predDF3[,1:7], CombinedPredictions2, TYPE=testingSet$TYPE)
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
results <- t(data.frame(results))
colnames(results) <- colnames(predDF4)
Results <- rbind(predDF4, results) 


write.csv(Results, 'OtherGenes_ml_results_8_algorithms.csv', row.names=TRUE)


Results[21,]
# predRF predRF2 predlda predGbm predKNN predRPART predGLM CombinedPredictions2 TYPE
# results      1       1     0.7    0.85    0.95      0.75     0.7                    1  100


#########################################################################################

# explore other genes that are involved in inflammation, interleukins, cytokins for pain,
# genes for radiation poisoning if they exist or are known, genes involved in the 
# bpa deliveries of cancers from bottled water if they exist, how? search NCBI gene and 
# keyword search

# and/or grep() the fields for ontology process, function for 'inflammation', 'pain'

table_more <- read.csv('table_more.csv', sep=',', header=TRUE, na.strings=c('','NA'))

Genefunction <- strsplit(as.character(table_more$Ontology_Function), ';')
functionGene <- lapply(Genefunction,'[',1)

functional <- data.frame(functionGene)
func_t <- t(functional)
row.names(func_t) <- NULL
functional <- data.frame(func_t)
colnames(functional) <- 'Gene_Function'
symbol <- data.frame(table_more$Symbol)
colnames(symbol) <- 'Symbol'

Genes <- cbind(symbol,functional)
genes <- Genes[!duplicated(Genes$Symbol),]

write.csv(genes, 'Gene_functions_16565.csv', row.names=FALSE)

Genes <- read.csv('Gene_functions_16565.csv', header=TRUE, sep=',')

alcohol <- Genes[grep('alcohol', Genes$Gene_Function),]
type0 <- data.frame(rep('alcohol',39))
colnames(type0) <- 'Type'
a <- cbind(alcohol, type0)

interleukin <- Genes[grep('interleukin', Genes$Gene_Function),]
type <- data.frame(rep('interleukin',17))
colnames(type) <- 'Type'
i <- cbind(interleukin,type)

cytokin <- Genes[grep('cytokin', Genes$Gene_Function),]
type1 <- data.frame(rep('cytokin', 30))
colnames(type1) <- 'Type'
c <- cbind(cytokin,type1)

collagen <- Genes[grep('collagen', Genes$Gene_Function),]
type2 <- data.frame(rep('collagen',4))
colnames(type2) <- 'Type'
cn <- cbind(collagen,type2)

pituitary <- Genes[grep('pituitary', Genes$Gene_Function),]
type3 <- data.frame(rep('pituitary',2))
colnames(type3) <- 'Type'
p <- cbind(pituitary,type3)

fat <- Genes[grep('fat', Genes$Gene_Function),]
type4 <- data.frame(rep('fat', 77))
colnames(type4) <- 'Type'
ft <- cbind(fat,type4)

iron <- Genes[grep('iron', Genes$Gene_Function),]
type5 <- data.frame(rep('iron',83))
colnames(type5) <- 'Type'
ir <- cbind(iron,type5)

sodium <- Genes[grep('sodium', Genes$Gene_Function),]
type6 <- data.frame(rep('sodium',5))
colnames(type6) <- 'Type'
s <- cbind(sodium,type6)

potassium <- Genes[grep('potassium', Genes$Gene_Function),]
type7 <- data.frame(rep('potassium', 28))
colnames(type7) <- 'Type'
po <- cbind(potassium, type7)

nitrogen <- Genes[grep('nitrogen', Genes$Gene_Function),]
type8 <- data.frame(rep('nitrogen',47))
colnames(type8) <- 'Type'
ni <- cbind(nitrogen, type8)

GenesInfo <- rbind(c,i,p,a,cn,ir,ft,s,po,ni)

StatsAll <- read.csv('StatsAllCompressed.csv', sep=',', header=TRUE)
colnames(StatsAll)[1] <- 'symbol'

Inflammation <- merge(GenesInfo, StatsAll, by.x="Symbol", by.y="symbol")

write.csv(Inflammation, 'minerals_inflammation_pain_fat_alcohol_Stats_UL_beadchip.csv', row.names=FALSE)

#both beadchip studies used UL and non-UL tissue from the uterus's in all samples
#These two studies focused on HOXA13 and SFRP4 separately for GSE120854 and GSE95101 respectively
#HOXA13 is associated with inducing UL formation when over-expressed in non-UL samples 
#according to GSE120854
#And SFRP4 is a gene product up-regulated functionally linked to 
#progesterone-regulated PR activation

# find other minerals essential to human health: calcium, magnesium, selenium,
# copper, chlorine, iodine, manganese, phosphorus, molybdenum, zinc, (iron, potassium,sodium)

calcium <- Genes[grep('calcium', Genes$Gene_Function),]
type9 <- data.frame(rep('calcium',327))
colnames(type9) <- 'Type'
ca <- cbind(calcium,type9)

magnesium <- Genes[grep('magnesium', Genes$Gene_Function),]
type10 <- data.frame(rep('magnesium', 168))
colnames(type10) <- 'Type'
mg <- cbind(magnesium,type10)

selenium <- Genes[grep('selenium', Genes$Gene_Function),]
type11 <- data.frame(rep('selenium', 9))
colnames(type11) <- 'Type'
se <- cbind(selenium,type11)

copper <- Genes[grep('copper', Genes$Gene_Function),]
type12 <- data.frame(rep('copper',29))
colnames(type12) <- 'Type'
cu <- cbind(copper, type12)

chlorine <- Genes[grep('chlorine', Genes$Gene_Function),]# 0 observations

iodine <- Genes[grep('iodine', Genes$Gene_Function),] # 0 observations

manganese <- Genes[grep('manganese', Genes$Gene_Function),] #0 observations

phosphorus <- Genes[grep('phosphorus', Genes$Gene_Function),]
type13 <- data.frame(rep('phosphorus', 26))
colnames(type13) <- 'Type'
ph <- cbind(phosphorus,type13)

molybdenum <- Genes[grep('molybdenum', Genes$Gene_Function),] #0 observations

zinc <- Genes[grep('zinc', Genes$Gene_Function),]
type14 <- data.frame(rep('zinc', 176))
colnames(type14) <- 'Type'
z <- cbind(zinc,type14)

genesMinerals <- rbind(GenesInfo,z,ph,cu,se,mg,ca)

Inflammation2 <- merge(genesMinerals, StatsAll, by.x="Symbol", by.y="symbol")
Inflammation3 <- Inflammation2[order(Inflammation2$FoldChange, decreasing=TRUE),]
write.csv(Inflammation3, 'minerals_inflammation_FC_ordered.csv', row.names=FALSE)

Inflammation4 <- Inflammation2[order(Inflammation2$FoldChangeMedian, decreasing=TRUE),]
write.csv(Inflammation4, 'minerals_inflammation_FC_med_ordered.csv', row.names=FALSE)

Inflammation5 <- Inflammation2[order(Inflammation2$DifferentialExpression, decreasing=TRUE),]
write.csv(Inflammation5, 'minerals_inflammation_DE_ordered.csv', row.names=FALSE)

Inflammation6 <- Inflammation2[order(Inflammation2$Median_DE, decreasing=TRUE),]
write.csv(Inflammation6, 'minerals_inflammation_DE_med_ordered.csv', row.names=FALSE)


#######################################################################################
StatsAll <- read.csv('StatsAllCompressed.csv', sep=',', header=TRUE)
Genes <- read.csv('Gene_functions_16565.csv', header=TRUE, sep=',')


grep('PPAR',Genes$Symbol)#receptor involved in CBD studies,Alzheimer & IBS studies
# PPAR gamma or PPARG is the target, it is less in UL by 1/3
grep('GW9662', Genes$Symbol)#0

grep('S100B', Genes$Symbol)#inflammatory protein, less in UL
grep('CB1', Genes$Symbol)#CBD receptor, not specific, 'CBD' is in >5 symbols
grep('NO', Genes$Symbol)#inflammatory response, >10 symbols with NO, no 'NO' exactly
grep('IL1B', Genes$Symbol)#inflammatory response, ~1/2 less in UL
grep('TNFA', Genes$Symbol)#inflammatory response, ~1/3 less in UL and >5 symbol variations

#######################################################################################
#I made and moved these 2 files, one is endocrine and the other various pain, cytokines,
#inflammation, some estrogen, and CBD receptor genes

GenesCBD_pain <- read.csv('GenesSumm1.csv', header=TRUE, sep=',')
GenesEndocrine <- read.csv('GenesSumm2.csv', header=TRUE, sep=',')
colnames(GenesCBD_pain)[1] <- 'Gene'
colnames(GenesEndocrine)[1] <- 'Gene'

Genes <- read.csv('Gene_functions_16565.csv', header=TRUE, sep=',')

CBD <- merge(GenesCBD_pain,Genes, by.x='Gene', by.y='Symbol')
endocrine <- merge(GenesEndocrine, Genes, by.x='Gene', by.y='Symbol')

DE_MEAN <- read.csv('DE_MEAN_21_83.csv', sep=',', header=TRUE)
DE_MEDIAN <- read.csv('DE_MEDIAN_21_83.csv', sep=',', header=TRUE)
FC_MEAN <- read.csv('FC_MEAN_20_83.csv', sep=',', header=TRUE)
FC_MEDIAN <-  read.csv('FC_MEDIAN_15_83.csv', sep=',', header=TRUE)

#table_more <- read.csv('table_more.csv', sep=',', header=TRUE)#262519X94
StatsAll <- read.csv('StatsAllCompressed.csv', sep=',', header=TRUE, row.names=1)#16152X83

CBD_all <- merge(CBD, StatsAll, by.x='Gene', by.y='Symbol')
endocrine_all <- merge(endocrine, StatsAll, by.x='Gene',  by.y='Symbol')

endoTopFC <- endocrine_all[order(endocrine_all$FoldChange, decreasing=TRUE),]
CBDTopFC <- CBD_all[order(CBD_all$FoldChange,decreasing=TRUE),]

write.csv(endoTopFC,'endocrine_bc_top_fc.csv', row.names=FALSE)
write.csv(CBDTopFC, 'CBD_pain_bc_top_fc.csv', row.names=FALSE)

library(dplyr)

expression <- mutate(CBDTopFC, expression = CBDTopFC$DifferentialExpression > 0)
expression$expression <- gsub('TRUE','up-regulated', expression$expression)
expression$expression <- gsub('FALSE','down-regulated', expression$expression)
expression <- expression[,c(1,87,2:86)]

png('CBD_pain_allGenes_UL.png', width=768, height=576)
g <- ggplot(expression, aes(x=FoldChange, y=as.factor(Gene)))
g= g+ xlab('Fold Change in UL to Non-UL')
g= g+ylab('CBD Inflammation Pain Genes')
g= g+ geom_point(aes(colour=expression),size=6, alpha=0.9)
g= g + geom_vline(xintercept=1)
g
dev.off()


expression2 <- mutate(endoTopFC, expression = endoTopFC$DifferentialExpression > 0)
expression2$expression <- gsub('TRUE','up-regulated', expression2$expression)
expression2$expression <- gsub('FALSE','down-regulated', expression2$expression)
expression2 <- expression2[,c(1,87,2:86)]

png('endocrine_bc_allGenes_UL.png', width=768, height=576)
g <- ggplot(expression2, aes(x=FoldChange, y=as.factor(Gene)))
g= g+ xlab('Fold Change in UL to Non-UL')
g= g+ylab('Endocrine Genes')
g= g+ geom_point(aes(colour=expression),size=6, alpha=0.9)
g= g + geom_vline(xintercept=1)
g
dev.off()
