#As so often before, we will rely heavily on Aaron Luns work here. See chapter
#20 of his OSCA book, currently available from 
#https://bioconductor.org/books/release/OSCA/integrating-with-protein-abundance.html
library(SingleCellExperiment)
library(scran)
library(scater)
library(Seurat)
library(DropletUtils)
library(DepecheR)
library(uwot)
#After having made a successful test with data from 10X, I found a dataset with
#50 000 cells from healthy individuals and with 82 CITE-seq antibodies, from
#this paper: 
#Broad immune activation underlies shared set point signatures for 
#vaccine responsiveness in healthy individuals and disease activity in patients
#with lupus. 
#The data is available from: 
#https://nih.figshare.com/articles/dataset/CITE-seq_protein-mRNA_single_cell_data_from_high_and_low_vaccine_responders_to_reproduce_Figs_4-6_and_associated_Extended_Data_Figs_/11349761?backTo=/collections/Data_and_software_code_repository_for_Broad_immune_activation_underlies_shared_set_point_signatures_for_vaccine_responsiveness_in_healthy_individuals_and_disease_activity_in_patients_with_lupus_Kotliarov_Y_Sparks_R_et_al_Nat_Med_DOI_https_d/4753772
raw_seurat_data <- readRDS("Data/H1_day0_demultilexed_singlets.RDS")

#It does for some reason not work to just convert this object into a singleCellExperiment. 
#Instead, we will have to take the raw data from the file, and put it together with
#some metadata. 
countsTranscript <- raw_seurat_data@data
metaDataObj <- raw_seurat_data@meta.data

#This metadata contains a lot of potentially confusing information that will 
#overlap with the analyses in the OSCA suite, and also a number of variables that
#do not vary, so only a small fraction is kept
slimMetaDataObj <- metaDataObj[,c("nUMI", "batch", "sampleid")]
#And then the protein assay: 
countsProt <- raw_seurat_data@assay$CITE@raw.data

#And based on this, we can now build our own SingleCellExperiment
sce <- SingleCellExperiment(assays = list(counts = countsTranscript))

altExp(sce) <- SingleCellExperiment(assays = list(counts = countsProt))

colData(sce) <- cbind(colData(sce), slimMetaDataObj)

sce

#This looks very neat. 

#As the ADT data is not sparse, we here simplify downstream analyses by saving
#it non-sparsely.
counts(altExp(sce)) <- as.matrix(counts(altExp(sce)))
counts(altExp(sce))[,1:10] # sneak peek

#This does indeed look good. 

#Now, we remove cells that do not show any background or true binding for at least
#50% of the antibodies. 
df.ab <- perCellQCMetrics(altExp(sce))
n.nonzero <- sum(!rowAlls(counts(altExp(sce)), value=0L))
ab.discard <- df.ab$detected <= n.nonzero/2
summary(ab.discard)
#   Mode   FALSE    TRUE 
#logical    58618      36 

#So we lose 36 cells by this procedure. 
dir.create("Diagnostics")
pdf("Diagnostics/Number_of_detected_ADTs.pdf")
hist(df.ab$detected, col='grey', main="", xlab="Number of detected ADTs", 
     breaks = 50)
abline(v=n.nonzero/2, col="red", lty=2)
dev.off()

#We also try a different strategy for exclusion, namely to use only the isotype
#controls, where we expect all cells to behave reasonably similarly: 
controls <- grep(".+sotype_PROT", rownames(altExp(sce)))
df.ab2 <- perCellQCMetrics(altExp(sce), subsets=list(controls=controls))
con.discard <- isOutlier(df.ab2$subsets_controls_sum, log=TRUE, type="lower", 
                         nmads = 4)
summary(con.discard)
pdf("Diagnostics/Number_of_detected_isotype_ADTs.pdf")
hist(log1p(df.ab2$subsets_controls_sum), col='grey', breaks=50,
     main="", xlab="Log-total count for isotype ADT controls per cell")
abline(v=log1p(attr(con.discard, "thresholds")["lower"]), col="red", lty=2)
dev.off()

#This really seems not to be a good strategy for this dataset, as as many as 
#5% of the cells are excluded here, not showing any expression at all. We
#instead hope that the more standardized transcriptome exclusions will be helpful

#And now, of course, we will also run exclusions on the transcriptome side of
#the data
mito <- grep("^MT-", rownames(sce))
df <- perCellQCMetrics(sce, subsets=list(Mito=mito))
mito.discard <- isOutlier(df$subsets_Mito_percent, type="higher")
detected.discard <- isOutlier(df$detected, type="lower")
summary(mito.discard)
#   Mode   FALSE    TRUE 
#logical    53739    4915 

#Here, we exclude quite a fraction of the cells. The detected.discard is not
#as powerful, and he overlap to the mito.discard is about 90%, so we will go
#for the simple setup of excluding the ones with extremely few markers detected
#and the ones with a high mito percentage. 
pdf("Diagnostics/Detected_prots_vs_MT_content.pdf")
plot(df$altexps_unnamed1_detected, df$subsets_Mito_percent)
dev.off()


discard <- ab.discard | mito.discard
summary(discard)
#Mode   FALSE    TRUE 
#logical    53707    4947 

sce <- sce[,!discard]

#

#Now, three different methods for normalization are presented by Aaron. The first
#normalization to library size, fails, as he correctly points out earlier in the
#tutorial, on the fact that we do expect very different protein abundance in 
#different cells. The last instead builds on the idea of isotype controls to be
#identically binding, which is of course also a fallacy, shown clearly in this 
#dataset. Therefore, we will go
#for the third approach that assumes that a majority of the markers are not 
#abundantly expressed in each individual cell. 

baseline <- inferAmbience(counts(altExp(sce)))
head(baseline)

#To be able to plot this, we need to subset the dataset
sceSubset <- sce[,sample(ncol(sce), 5000)]
plotExpression(altExp(sceSubset), features=rownames(altExp(sceSubset)), exprs_values="counts") +
    scale_y_log10() + 
    geom_point(data=data.frame(x=names(baseline), y=baseline), 
    mapping=aes(x=x, y=y), cex=3)
dev.copy(pdf, "Diagnostics/Protein_counts_all.pdf", width = 30, height = 10)
dev.off()

#Now the median expression is calculated for each cell
sf.amb <- medianSizeFactors(altExp(sce), reference=baseline)
summary(sf.amb)

#Now, we visually show that there are compositional divergences between the 
#protein and the transcriptome data, which makes sense. This works badly for a
#dataset of this size, and is not necessary as it is mainly included to prove 
#a point in Aarons discussion. 
#tagdata <- logNormCounts(altExp(sceSubset)) # library size factors by default.
#g <- buildSNNGraph(tagdata, k=20, d=NA) # no need for PCA, see below.
#clusters <- igraph::cluster_walktrap(g)$membership
#
#sf.lib <- librarySizeFactors(altExp(sce))
#summary(sf.lib)
#
#plot(sf.lib, sf.amb, log="xy", col=clusters, 
#     xlab="Library size factors (tag)",
#     ylab="DESeq-like size factors (tag)")
#abline(0, 1, col="grey", lty=2)
#dev.copy(pdf, "Diagnostics/Lib_size_factors_vs_protein_size_factors.pdf")
#dev.off()

#It is possible that there are small divergences there, but they are clearly not
#as visible as for the other dataset used in the tutorial. 

#Now over to normalization
sizeFactors(altExp(sce)) <- sf.amb
sce <- logNormCounts(sce)
altExpNames(sce) <- "Antibody_capture"

#There has been a change to the code since the release of the book, so
#this is what one is supposed to do to get the normalization done for the altExp
#sce1 <- logNormCounts(sce, use.altexps=TRUE)
sce <- applySCE(sce, logNormCounts, WHICH = "Antibody_capture")

# Checking that we have normalized values:
assayNames(sce)
assayNames(altExp(sce))

#It works. Now we have a normalized assay, that is ready to be exported
#for gating elsewhere. 
#Before doing so, though, we need to name the columns 
colData(sce)$Event <- seq(1, ncol(sce))

library(flowCore)
#Before exporting, we need to reorganize the columns as the current order makes
#it impossible to find anything. They will now be ordered alphabetically and 
#numerically
flowMat <- t(logcounts(altExp(sce)))
rawColNames <- colnames(flowMat)

#Now, we pick out the ones that need extra zeros to work fine. 
cDColNames <- rawColNames[grep("CD", rawColNames)]
nInts <- nchar(gsub("[^0-9]+", "", cDColNames))
cDColNames[which(nInts == 2)] <- gsub("CD", "CD0", cDColNames[which(nInts == 2)])
cDColNames[which(nInts == 1)] <- gsub("CD", "CD00", cDColNames[which(nInts == 1)])
updatedColNames <- rawColNames
updatedColNames[grep("CD", rawColNames)] <- cDColNames

flowMatOrdered <- flowMat[,order(updatedColNames)]

#And here, the flow matrix and the whole sce are exported
saveRDS(sce, "Data/sce_post_normalization.rds")
write.csv(flowMatOrdered, "Data/FlowMatOrdered.csv")

