
#read in the files for all genes in common and a new one on hormones from endocrine system with
#summaries

all <- read.csv('all_common_12173_130_fold_magnitude.csv', sep=',', header=TRUE, row.names=1)
endocrine <-read.csv('GenesSumm2.csv', sep=',', header=TRUE)
Inflammation <- read.csv('GenesSumm1.csv', sep=',', header=TRUE)
colnames(Inflammation)[1] <- 'Gene'
colnames(endocrine)[1] <- 'Gene'

ENDOCRINE <- merge(endocrine, all, by.x='Gene', by.y='GENE_SYMBOL')
INFLAMMATION <-merge(Inflammation, all, by.x='Gene', by.y='GENE_SYMBOL')

write.csv(ENDOCRINE, 'endocrine_UL_stats.csv', row.names=FALSE)
write.csv(INFLAMMATION,'pain_inflm_cbd_UL.csv', row.names=FALSE)

thresh50_endocrine <- subset(ENDOCRINE, foldChange>1.3 | foldChange<.7)
thresh50_pain <- subset(INFLAMMATION, foldChange>1.4 | foldChange<.6)

library(dplyr)

end50 <- mutate(thresh50_endocrine, expression=thresh50_endocrine$DE>0)
end50$expression <- gsub('TRUE','up-regulated',end50$expression)
end50$expression <- gsub('FALSE','down-regulated', end50$expression)
write.csv(end50,'endocrineMost.csv', row.names=FALSE)
library(ggplot2)

png('endocrine_UL.png', width=768, height=576)
g <- ggplot(end50, aes(x=as.factor(Gene), y=foldChange))
g= g+xlab('Endocrine Genes')
g= g+ ylab('Fold Change in UL to Non-UL')
g= g+ geom_point(aes(colour=expression),size=6, alpha=0.9)
g
dev.off()

expressionAll <- mutate(ENDOCRINE, expression=ENDOCRINE$DE>0)
expressionAll$expression <- gsub('TRUE','up-regulated',expressionAll$expression)
expressionAll$expression <- gsub('FALSE','down-regulated', expressionAll$expression)
expressionAll <- expressionAll[,c(1:9,133,10:132)]
write.csv(expressionAll, 'expression_All_endo_genes_UL.csv', row.names=FALSE)

png('endocrine_allGenes_UL.png', width=768, height=576)
g <- ggplot(expressionAll, aes(x=foldChange, y=as.factor(Gene)))
g= g+ xlab('Fold Change in UL to Non-UL')
g= g+ylab('Endocrine Genes')
g= g+ geom_point(aes(colour=expression),size=6, alpha=0.9)
g
dev.off()


expressionAll_inf <- mutate(INFLAMMATION, expression=INFLAMMATION$DE>0)
expressionAll_inf$expression <- gsub('TRUE','up-regulated',expressionAll_inf$expression)
expressionAll_inf$expression <- gsub('FALSE','down-regulated', expressionAll_inf$expression)
expressionAll_inf <- expressionAll_inf[,c(1:9,133,10:132)]
write.csv(expressionAll_inf, 'expression_All_endo_genes_UL.csv', row.names=FALSE)

png('inflammation_allGenes_UL.png', width=768, height=576)
g <- ggplot(expressionAll_inf, aes(x=foldChange, y=as.factor(Gene)))
g= g+ xlab('Fold Change in UL to Non-UL')
g= g+ylab('CBD Inflammation Pain Genes')
g= g+ geom_point(aes(colour=expression),size=6, alpha=0.9)
g
dev.off()

