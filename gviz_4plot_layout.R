

# all Forward Strand genes on chromosome 22
fwd <- subset(allMembers, allMembers$strand=='+')
chr <- as.character(unique(fwd$chromosome)) #[1] "chr11" "chr22" "chr12" "chr17"
grtrack <- GeneRegionTrack(fwd, genome = "hg38", chromosome = chr[2], 
                           strand=fwd$strand,
                           gene=fwd$genes, symbol=fwd$genes, 
                           transcriptAnnotation = "genes", #adds the genes to bottom track
                           background.title = "brown", lwd=1, cex=2,
                           stackHeight=0.5,
                           name = "Forward Strand 22q13.1")

itrack <- IdeogramTrack(genome = "hg38", chromosome = chr[2])
gtrack <- GenomeAxisTrack()
aTrack <- AnnotationTrack(start = c(37282027,40044817), 
                          width = c(33315,290992), 
                          chromosome = "chr22", strand = rep(c("+", "+"), 
                          c(1, 1)), group = rep(c("CYTH4", "TNRC6B"), 
                          c(1, 1)), genome = "hg38", + name = "Forward")
#plotTracks(aTrack, groupAnnotation = "group")


png('fwdAllChr22.png',width=500,height=400)
plotTracks(list(itrack, gtrack, grtrack, aTrack),
           featureAnnotation = "id",  lwd=1, cex=1, #shape='arrow',
           background.panel = "beige", from = 3.7e+07, to = 4.05e+07,
           add53=TRUE,  littleTicks = TRUE, labelPos='alternating',
           #add35=TRUE,adds 3' to 5' and 5' to 3' on ideogram header
           showBandId = TRUE, cex.bands = 0.5, reverseStrand=FALSE,
           scale=1,showOverplotting=TRUE,groupAnnotation = "group",
           background.title = "#005544")# title color green-teal
dev.off()

