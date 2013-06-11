### R code from vignette source 'apcluster.Rnw'

###################################################
### code chunk number 1: Init
###################################################
options(width=75)
set.seed(0)
library(apcluster)
apclusterVersion <- packageDescription("apcluster")$Version
apclusterDateRaw <- packageDescription("apcluster")$Date
apclusterDateYear <- as.numeric(substr(apclusterDateRaw, 1, 4))
apclusterDateMonth <- as.numeric(substr(apclusterDateRaw, 6, 7))
apclusterDateDay <- as.numeric(substr(apclusterDateRaw, 9, 10))
apclusterDate <- paste(month.name[apclusterDateMonth], " ",
                     apclusterDateDay, ", ",
                     apclusterDateYear, sep="")


###################################################
### code chunk number 2: InstallAPCluster (eval = FALSE)
###################################################
## install.packages("apcluster")


###################################################
### code chunk number 3: LoadAPCluster (eval = FALSE)
###################################################
## library(apcluster)


###################################################
### code chunk number 4: OpenVignette (eval = FALSE)
###################################################
## vignette("apcluster")


###################################################
### code chunk number 5: ShowHelp (eval = FALSE)
###################################################
## help(apcluster)


###################################################
### code chunk number 6: CreateDataSet1
###################################################
cl1 <- cbind(rnorm(30, 0.3, 0.05), rnorm(30, 0.7, 0.04))
cl2 <- cbind(rnorm(30, 0.7, 0.04), rnorm(30, 0.4, .05))
x1 <- rbind(cl1, cl2)
plot(x1, xlab="", ylab="", pch=19, cex=0.8)


###################################################
### code chunk number 7: APClusterDataSet1
###################################################
apres1a <- apcluster(negDistMat(r=2), x1)


###################################################
### code chunk number 8: APClusterDataSet1b
###################################################
s1 <- negDistMat(x1, r=2)
apres1b <- apcluster(s1)


###################################################
### code chunk number 9: ShowHelpAPResult (eval = FALSE)
###################################################
## help(APResult)


###################################################
### code chunk number 10: ShowResultAPClusterDataSet1
###################################################
apres1a


###################################################
### code chunk number 11: PlotResultAPClusterDataSet1
###################################################
plot(apres1a, x1)


###################################################
### code chunk number 12: HeatmapResultAPClusterDataSet1
###################################################
heatmap(apres1a)


###################################################
### code chunk number 13: HeatmapResultAPClusterDataSet1b
###################################################
heatmap(apres1b, s1)


###################################################
### code chunk number 14: APClusterDataSet1Details
###################################################
apres1c <- apcluster(s1, details=TRUE)


###################################################
### code chunk number 15: PlotAPClusterDataSet1Details
###################################################
plot(apres1c)


###################################################
### code chunk number 16: APClusterDataSet1convits15
###################################################
apres1c <- apcluster(s1, convits=15, details=TRUE)
apres1c


###################################################
### code chunk number 17: CreateDataSet2
###################################################
cl3 <- cbind(rnorm(20, 0.50, 0.03), rnorm(20, 0.72, 0.03))
cl4 <- cbind(rnorm(25, 0.50, 0.03), rnorm(25, 0.42, 0.04))
x2 <- rbind(x1, cl3, cl4)
plot(x2, xlab="", ylab="", pch=19, cex=0.8)


###################################################
### code chunk number 18: APClusterDataSet2
###################################################
apres2a <- apcluster(negDistMat(r=2), x2)
plot(apres2a, x2)


###################################################
### code chunk number 19: APClusterDataSet2q0
###################################################
apres2b <- apcluster(negDistMat(r=2), x2, q=0)
plot(apres2b, x2)


###################################################
### code chunk number 20: PlotAPClusterDataSet2q08
###################################################
apres2c <- apcluster(negDistMat(r=2), x2, q=0.8)
plot(apres2c, x2)


###################################################
### code chunk number 21: APClusterDataSet2q08showp
###################################################
apres2c@p


###################################################
### code chunk number 22: HeatmapResultAPClusterDataSet2q08
###################################################
heatmap(apres2c)


###################################################
### code chunk number 23: PreferenceRangeDataSet2
###################################################
preferenceRange(apres2b@sim)


###################################################
### code chunk number 24: APClusterKDataSet2
###################################################
apres2d <- apclusterK(negDistMat(r=2), x2, K=2, verbose=TRUE)
plot(apres2d, x2)


###################################################
### code chunk number 25: IrisData1
###################################################
data(iris)
apIris1 <- apcluster(negDistMat(r=2), iris)
apIris1


###################################################
### code chunk number 26: IrisDataPlot1
###################################################
plot(apIris1, iris)


###################################################
### code chunk number 27: IrisDataHeatmap1
###################################################
heatmap(apIris1)


###################################################
### code chunk number 28: IrisData2
###################################################
data(iris)
apIris2 <- apcluster(negDistMat(r=2), iris, q=0)
apIris2


###################################################
### code chunk number 29: IrisDataPlot
###################################################
plot(apIris2, iris)


###################################################
### code chunk number 30: IrisDataHeatmap2
###################################################
heatmap(apIris2)


###################################################
### code chunk number 31: AggExClusterDataSet1
###################################################
aggres1a <- aggExCluster(negDistMat(r=2), x1)
aggres1a


###################################################
### code chunk number 32: DendrogramAggExClusterDataSet1
###################################################
plot(aggres1a)


###################################################
### code chunk number 33: HeatmapAggExClusterDataSet1
###################################################
heatmap(aggres1a, s1)


###################################################
### code chunk number 34: ExtractAggExClustersDataSet1
###################################################
cl1a <- cutree(aggres1a, k=2)
cl1a
plot(cl1a, x1)


###################################################
### code chunk number 35: AggExClusterAPDataSet2q08
###################################################
aggres2a <- aggExCluster(x=apres2c)
aggres2a


###################################################
### code chunk number 36: DendrogramAggExAPDataSet2
###################################################
plot(aggres2a)


###################################################
### code chunk number 37: DendrogramAggExAPDataSet2b
###################################################
plot(aggres2a, showSamples=TRUE, nodePar=list(pch=NA, lab.cex=0.4))


###################################################
### code chunk number 38: HeatmapAggExAPDataSet2
###################################################
heatmap(aggres2a)


###################################################
### code chunk number 39: PlotAggExAPDataSet2k25
###################################################
par(mfrow=c(2,2))
for (k in 5:2)
    plot(aggres2a, x2, k=k, main=paste(k, "clusters"))


###################################################
### code chunk number 40: APClusterLevDataSet3
###################################################
cl5 <- cbind(rnorm(100, 0.3, 0.05), rnorm(100, 0.7, 0.04))
cl6 <- cbind(rnorm(100, 0.70, 0.04), rnorm(100, 0.4, 0.05))
x3 <- rbind(cl5, cl6)
apres3 <- apclusterL(s=negDistMat(r=2), x=x3, frac=0.1, sweeps=5, p=-0.2)
apres3
plot(apres3, x3)


###################################################
### code chunk number 41: APClusterLevResultDataSet3
###################################################
dim(apres3@sim)
apres3@sel
apres3@netsimLev


###################################################
### code chunk number 42: APClusterLevDataSet3Heat
###################################################
heatmap(apres3)


###################################################
### code chunk number 43: LoadCh22Promoters
###################################################
data(ch22Promoters)
names(ch22Promoters)[1:5]
substr(ch22Promoters[1:5], 951, 1000)


###################################################
### code chunk number 44: SimCh22Promoters
###################################################
library(kernlab)
promSim <- kernelMatrix(stringdot(length=6, type="spectrum"), ch22Promoters)
rownames(promSim) <- names(ch22Promoters)
colnames(promSim) <- names(ch22Promoters)


###################################################
### code chunk number 45: APCh22Promoters
###################################################
promAP <- apcluster(promSim, q=0)
promAP


###################################################
### code chunk number 46: HeatMapAPCh22Promoters
###################################################
heatmap(promAP, promSim, Rowv=FALSE, Colv=FALSE)


###################################################
### code chunk number 47: aggExCh22Promoters
###################################################
promAgg <- aggExCluster(promSim, promAP)


###################################################
### code chunk number 48: DendrogramAPCh22Promoters
###################################################
plot(promAgg)


###################################################
### code chunk number 49: ExtractAggCh22Promoters
###################################################
prom5 <- cutree(promAgg, k=5)
prom5


###################################################
### code chunk number 50: HeatMap5Ch22Promoters
###################################################
heatmap(prom5, promSim, Rowv=FALSE, Colv=FALSE)


###################################################
### code chunk number 51: NegDistMatDataSet2
###################################################
s <- negDistMat(x2)


###################################################
### code chunk number 52: CreateToyData
###################################################
ex <- matrix(c(0, 0.5, 0.8, 1, 0, 0.2, 0.5, 0.7,
               0.1, 0, 1, 0.3, 1, 0.8, 0.2), 5, 3,byrow=TRUE)
ex


###################################################
### code chunk number 53: NegEuclDistMatToyData
###################################################
negDistMat(ex)


###################################################
### code chunk number 54: NegSqEuclDistMatToyData
###################################################
negDistMat(ex, r=2)


###################################################
### code chunk number 55: NegMaxDistToyData
###################################################
negDistMat(ex, method="maximum")


###################################################
### code chunk number 56: NegManhattanDistToyData
###################################################
negDistMat(ex, method="manhattan")


###################################################
### code chunk number 57: NegCanberraDistToyData
###################################################
negDistMat(ex, method="canberra")


###################################################
### code chunk number 58: NegMinkowskiDistToyData
###################################################
negDistMat(ex, method="minkowski", p=3)


###################################################
### code chunk number 59: GetFunction
###################################################
sim <- negDistMat(r=2)
is.function(sim)
apcluster(sim, x1)


###################################################
### code chunk number 60: RBFKernelToyData
###################################################
expSimMat(ex)


###################################################
### code chunk number 61: LaplaceKernelToyData
###################################################
expSimMat(ex, r=1)


###################################################
### code chunk number 62: PearsonToyData
###################################################
corSimMat(ex, method="pearson")


###################################################
### code chunk number 63: SpearmanToyData
###################################################
corSimMat(ex, method="spearman")


###################################################
### code chunk number 64: TruncDistToyData
###################################################
linSimMat(ex, w=1.2)


###################################################
### code chunk number 65: LinKernelToyData
###################################################
linKernel(ex[2:5,])


###################################################
### code chunk number 66: NormLinKernelToyData
###################################################
linKernel(ex[2:5,], normalize=TRUE)


###################################################
### code chunk number 67: RectangularNegDistMatDataSet1
###################################################
sel <- sort(sample(1:nrow(x1), ceiling(0.08 * nrow(x1))))
sel
s1r <- negDistMat(x1, sel, r=2)
dim(s1r)
s1r[1:7,]


###################################################
### code chunk number 68: SimFunCh22Promoters
###################################################
spectrumK6 <- function(x, sel=NA)
{
    if (any(is.na(sel)))
    {
        s <- kernelMatrix(stringdot(length=6, type="spectrum"), x)
        rownames(s) <- names(x)
        colnames(s) <- names(x)
    }
    else
    {
        s <- kernelMatrix(stringdot(length=6, type="spectrum"), x, x[sel])
        rownames(s) <- names(x)
        colnames(s) <- names(x)[sel]
    }

    s
}


###################################################
### code chunk number 69: APCh22Promoters
###################################################
promAPL <- apclusterL(s=spectrumK6, ch22Promoters, frac=0.1, sweeps=10,
                        p=promAP@p)
promAPL


###################################################
### code chunk number 70: CreateLabeledToyData
###################################################
x3 <- c(1, 2, 3, 7, 8, 9)
names(x3) <- c("a", "b", "c", "d", "e", "f")
s3 <- negDistMat(x3, r=2)


###################################################
### code chunk number 71: ShowToyDataLabels
###################################################
s3
colnames(s3)


###################################################
### code chunk number 72: ClusterLabeledToyData
###################################################
apres3a <-apcluster(s3)
apres3a
apres3a@exemplars
apres3a@clusters


###################################################
### code chunk number 73: ExtractLabelsFromClusterToyData
###################################################
apres3a@exemplars
labels(apres3a, type="names")
labels(apres3a, type="exemplars")
labels(apres3a, type="enum")


###################################################
### code chunk number 74: HeatmapResultAPClusterDataSetq08b
###################################################
heatmap(apres2c, sideColors=c("darkgreen", "yellowgreen"),
        col=terrain.colors(12), Rowv=FALSE, dendScale=0.6)


###################################################
### code chunk number 75: HeatmapResultAPClusterDataSet2q08c
###################################################
heatmap(apres2c, sideColors=rainbow(length(apres2c)), Rowv=FALSE, Colv=FALSE,
        cexRow=(0.2 + 1 / log10(nrow(apres2c@sim))),
        cexCol=(0.2 + 1 / log10(nrow(apres2c@sim))))


###################################################
### code chunk number 76: PlotAddLegend
###################################################
plot(apres2a, x2)
legend("bottomleft", legend=paste("Cluster", 1:length(apres2a)),
       col=rainbow(length(apres2a)), pch=19)


###################################################
### code chunk number 77: PlotOnlyLegend
###################################################
plot.new()
par(oma=rep(0, 4), mar=rep(0, 4))
legend("center", legend=paste("Cluster", 1:length(apres2c)),
       col=rainbow(length(apres2c)), pch=19)


###################################################
### code chunk number 78: GetBibTeX (eval = FALSE)
###################################################
## toBibtex(citation("apcluster"))


