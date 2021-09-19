## ----echo=FALSE---------------------------------------------------------------
.standalone=FALSE

## ---- message=FALSE, warning=FALSE, include=!.standalone----------------------
library("AnnotationDbi")
library("abind")
library("beeswarm")
library("Biobase")
library("biomaRt")
library("broom")
library("colorspace")
library("cowplot")
library("dendsort")
library("DESeq2")
library("doParallel")
library("dplyr")
library("foreach")
library("forestplot")
library("genefilter")
library("ggbeeswarm")
library("ggdendro")
library("ggplot2")
#library("ggtern")
library("glmnet")
library("grid")
library("gridExtra")
library("gtable")
library("hexbin")
library("IHW")
library("ipflasso")
library("knitr")
library("limma")
library("magrittr")
library("maxstat")
library("nat")
library("org.Hs.eg.db")
library("BloodCancerMultiOmics2017")
library("pheatmap")
library("piano")
library("readxl")
library("RColorBrewer")
library("reshape2")
library("Rtsne")
library("scales")
library("SummarizedExperiment")
library("survival")
library("tibble")
library("tidyr")
library("xtable")

## ---- message=FALSE, include=!exists(".standalone"), eval=!exists(".standalone")----
#  library("BloodCancerMultiOmics2017")
#  library("Biobase")
#  library("ggplot2")
#  library("gtable")
#  library("grid")
#  library("dplyr")
#  library("gridExtra")

## ----echo=FALSE---------------------------------------------------------------
plotDir = ifelse(exists(".standalone"), "", "part01/")
if(plotDir!="") if(!file.exists(plotDir)) dir.create(plotDir)

## -----------------------------------------------------------------------------
options(stringsAsFactors=FALSE)

## -----------------------------------------------------------------------------
data("drpar", "drugs", "patmeta", "mutCOM")

## -----------------------------------------------------------------------------
# PATIENTS
patM = colnames(drpar)

# DRUGS
drM = rownames(drpar)
drM = drM[!drM %in% "D_CHK"] # remove combintation of 2 drugs: D_CHK

## -----------------------------------------------------------------------------
bwScale = c("0"="white","1"="black","N.A."="grey90")
lfsize = 16 # legend font size

## -----------------------------------------------------------------------------
drugs$target_category = as.character(drugs$target_category)
drugs$group = NA
drugs$group[which(drugs$approved_042016==1)] = "FDA approved"
drugs$group[which(drugs$devel_042016==1)] = "clinical development/\ntool compound"

## ----echo=FALSE---------------------------------------------------------------
res = table(drugs[,c("target_category","group")])
knitr::kable(res[order(res[,1], decreasing=TRUE),])

## ----echo=FALSE---------------------------------------------------------------
goM = BloodCancerMultiOmics2017:::plotPathways(dat=drugs[drM,])

## ----dr_desc_M, echo=FALSE, fig.path=plotDir, dev=c("png", "pdf"), fig.width=goM[["figure"]][["width"]], fig.height=goM[["figure"]][["height"]]----
#FIG# 1C
grid.draw(goM[["figure"]][["plot"]])

## ---- echo=FALSE--------------------------------------------------------------
knitr::knit_hooks$set(crop=NULL)

## ----dr_desc_M_leg, echo=FALSE, fig.path=plotDir, dev=c("png", "pdf"), fig.width=goM[["legend"]][["width"]], fig.height=goM[["legend"]][["height"]], out.height=120, out.width=240----
#FIG# 1C
grid.draw(goM[["legend"]][["plot"]])

## ---- echo=FALSE--------------------------------------------------------------
knitr::knit_hooks$set(crop=knitr:::hook_pdfcrop)

## ----echo=FALSE---------------------------------------------------------------
knitr::kable(data.frame(sort(table(patmeta[patM, "Diagnosis"]), decreasing=TRUE)))

## ----echo=FALSE---------------------------------------------------------------
goM = BloodCancerMultiOmics2017:::plotPatientStat(pats=patM, gap=c(30,160))

## ----pt_desc_M, echo=FALSE, fig.path=plotDir, dev=c("png", "pdf"), fig.width=goM[["figure"]][["width"]], fig.height=goM[["figure"]][["height"]]----
#FIG# 1B
grid.draw(goM[["figure"]][["plot"]])

## ---- echo=FALSE--------------------------------------------------------------
knitr::knit_hooks$set(crop=NULL)

## ----pt_desc_M_leg, echo=FALSE, fig.path=plotDir, dev=c("png", "pdf"), fig.width=goM[["legend"]][["width"]], fig.height=goM[["legend"]][["height"]], out.height=120, out.width=240----
#FIG# 1B
grid.draw(goM[["legend"]][["plot"]])

## ---- echo=FALSE--------------------------------------------------------------
knitr::knit_hooks$set(crop=knitr:::hook_pdfcrop)

## -----------------------------------------------------------------------------
# select CLL samples
patM = patM[patmeta[patM,"Diagnosis"]=="CLL"]

ighv = factor(setNames(patmeta[patM,"IGHV"], nm=patM), levels=c("U","M"))

mut1 = c("del17p13", "del11q22.3", "trisomy12", "del13q14_any")
mut2 = c("TP53", "ATM", "SF3B1", "NOTCH1", "MYD88")

mc = assayData(mutCOM)$binary[patM,]

## SELECTION OF MUTATIONS
# # include mutations with at least incidence of 4
mut2plot = names(which(sort(colSums(mc, na.rm=TRUE), decreasing=TRUE)>3))

# remove chromothrypsis
mut2plot = mut2plot[-grep("Chromothripsis", mut2plot)]
# divide mutations into gene mut and cnv
mut2plotSV = mut2plot[grep("[[:lower:]]", mut2plot)]
mut2plotSP = mut2plot[grep("[[:upper:]]", mut2plot)]

# remove some other things (it is quite manual thing, so be careful)
# IF YOU WANT TO REMOVE SOME MORE MUTATIONS JUST ADD THE LINES HERE!
mut2plotSV = mut2plotSV[-grep("del13q14_mono", mut2plotSV)]
mut2plotSV = mut2plotSV[-grep("del13q14_bi", mut2plotSV)]
mut2plotSV = mut2plotSV[-grep("del14q24.3", mut2plotSV)]

# rearrange the top ones to match the order in mut1 and mut2
mut2plotSV = c(mut1, mut2plotSV[!mut2plotSV %in% mut1])
mut2plotSP = c(mut2, mut2plotSP[!mut2plotSP %in% mut2])

factors = data.frame(assayData(mutCOM)$binary[patM, c(mut2plotSV, mut2plotSP)],
                     check.names=FALSE)
# change del13q14_any to del13q14
colnames(factors)[which(colnames(factors)=="del13q14_any")] = "del13q14"
mut2plotSV = gsub("del13q14_any", "del13q14", mut2plotSV)
# change it to factors
for(i in 1:ncol(factors)) {
  factors[,i] = factor(factors[,i], levels=c(1,0))
}

ord = order(factors[,1], factors[,2], factors[,3], factors[,4], factors[,5],
            factors[,6], factors[,7], factors[,8], factors[,9], factors[,10],
            factors[,11], factors[,12], factors[,13], factors[,14],
            factors[,15], factors[,16], factors[,17], factors[,18],
            factors[,19], factors[,20], factors[,21], factors[,22],
            factors[,23], factors[,24], factors[,25], factors[,26],
            factors[,27], factors[,28], factors[,29], factors[,30],
            factors[,31], factors[,32])

factorsord = factors[ord,]
patM = patM[ord]

(c(mut2plotSV, mut2plotSP))

## -----------------------------------------------------------------------------
plotDF = meltWholeDF(factorsord)
plotDF$Mut =
  ifelse(sapply(plotDF$X,
                function(x) grep(x, list(mut2plotSV, mut2plotSP)))==1,"SV","SP")
plotDF$Status = "N.A."
plotDF$Status[plotDF$Measure==1 & plotDF$Mut=="SV"] = "1a"
plotDF$Status[plotDF$Measure==1 & plotDF$Mut=="SP"] = "1b"
plotDF$Status[plotDF$Measure==0] = "0"
plotDF$Status = factor(plotDF$Status, levels=c("1a","1b","0","N.A."))

plotDF$Y = factor(plotDF$Y, levels=patM)
plotDF$X = factor(plotDF$X, levels=rev(colnames(factorsord)))

mutPL = ggplotGrob(
  ggplot(data=plotDF, aes(x=Y, y=X, fill=Status)) + geom_tile() +
    scale_fill_manual(
      values=c("0"="white","1a"="forestgreen","1b"="navy","N.A."="grey90"),
      name="Mutation", labels=c("CNV","Gene mutation","WT","NA")) +
    ylab("") + xlab("") +
    geom_vline(xintercept=seq(0.5,length(patM)+1,5), colour="grey60") +
    geom_hline(yintercept=seq(0.5,ncol(factorsord)+1,1), colour="grey60") +
    scale_y_discrete(expand=c(0,0)) + scale_x_discrete(expand=c(0,0)) +
    theme(axis.ticks=element_blank(), axis.text.x=element_blank(),
          axis.text.y=element_text(
            size=60, face=ifelse(levels(plotDF$X) %in% mut2plotSV,
                                 "plain","italic")),
          axis.text=element_text(margin=unit(0.5,"cm"), colour="black"),
          legend.key = element_rect(colour = "black"),
          legend.text=element_text(size=lfsize),
          legend.title=element_text(size=lfsize)))

res = table(plotDF[,c("X","Measure")])
knitr::kable(res[order(res[,2], decreasing=TRUE),])

## -----------------------------------------------------------------------------
ageDF = data.frame(Factor="Age",
                   PatientID=factor(patM, levels=patM),
                   Value=patmeta[patM,c("Age4Main")])

agePL = ggplotGrob(
  ggplot(ageDF, aes(x=PatientID, y=Factor, fill=Value)) + geom_tile() +
    scale_fill_gradient(low = "gold", high = "#3D1F00", na.value="grey92",
                        name="Age", breaks=c(40,60,80)) +
    theme(axis.ticks=element_blank(),
          axis.text=element_text(size=60, colour="black",
                                 margin=unit(0.5,"cm")),
          legend.text=element_text(size=lfsize),
          legend.title=element_text(size=lfsize)))

hist(ageDF$Value, col="slategrey", xlab="Age", main="")

## -----------------------------------------------------------------------------
sexDF = data.frame(Factor="Sex", PatientID=factor(patM, levels=patM),
                   Value=patmeta[patM, "Gender"])

sexPL = ggplotGrob(
  ggplot(sexDF, aes(x=PatientID, y=Factor, fill=Value)) + geom_tile() +
    scale_fill_manual(values=c("f"="maroon","m"="royalblue4","N.A."="grey90"),
                      name="Sex", labels=c("Female","Male","NA")) +
    theme(axis.ticks=element_blank(),
          axis.text=element_text(size=60, colour="black",
                                 margin=unit(0.5,"cm")),
          legend.key = element_rect(colour = "black"),
          legend.text=element_text(size=lfsize),
          legend.title=element_text(size=lfsize)))

table(sexDF$Value)

## -----------------------------------------------------------------------------
treatDF = data.frame(Factor="Treated", PatientID=factor(patM, levels=patM),
                     Value=ifelse(patmeta[patM, "IC50beforeTreatment"], 0, 1))
treatDF$Value[is.na(treatDF$Value)] = "N.A."
treatDF$Value = factor(treatDF$Value, levels=c("0","1","N.A."))

treatPL = ggplotGrob(
  ggplot(treatDF, aes(x=PatientID, y=Factor, fill=Value)) +geom_tile() +
    scale_fill_manual(values=bwScale, name="Treated",
                      labels=c("0"="No","1"="Yes","N.A."="NA")) +
    theme(axis.ticks=element_blank(),
          axis.text=element_text(size=60, colour="black",
                                 margin=unit(0.5,"cm")),
          legend.key = element_rect(colour = "black"),
          legend.text=element_text(size=lfsize),
          legend.title=element_text(size=lfsize)))

table(treatDF$Value)

## -----------------------------------------------------------------------------
ighvDF = data.frame(Factor="IGHV", PatientID=factor(patM, levels=patM),
                    Value=patmeta[patM, "IGHV"])
ighvDF$Value = ifelse(ighvDF$Value=="M", 1, 0)
ighvDF$Value[is.na(ighvDF$Value)] = "N.A."
ighvDF$Value = factor(ighvDF$Value, levels=c("0","1","N.A."))

ighvPL = ggplotGrob(
  ggplot(ighvDF, aes(x=PatientID, y=Factor, fill=Value)) + geom_tile() +
    scale_fill_manual(values=bwScale, name="IGHV",
                      labels=c("0"="Unmutated","1"="Mutated","N.A."="NA")) +
    theme(axis.ticks=element_blank(), 
          axis.text=element_text(size=60, colour="black", margin=unit(0.5,"cm")),
          legend.key=element_rect(colour = "black"),
          legend.text=element_text(size=lfsize),
          legend.title=element_text(size=lfsize)))

table(ighvDF$Value)

## ----echo=FALSE---------------------------------------------------------------
nX = length(patM)
nY = ncol(factorsord)
unY1 = 0.6*1.6
unY2 = 0.6*1.8
unX = 0.2
sp = 0.001
wdths = c(6, unX*nX, sp)
hghts = c(sp, unY1,unY1,unY1,unY1, 0.8, sp, sp ,unY2*nY, sp)
gt = gtable(widths=unit(wdths, "in"), heights=unit(hghts, "in"))

# add the plots
gt = gtable_add_grob(gt, sexPL$grobs[[whichInGrob(sexPL, "panel")]], 2, 2)
gt = gtable_add_grob(gt, treatPL$grobs[[whichInGrob(treatPL, "panel")]], 3, 2)
gt = gtable_add_grob(gt, agePL$grobs[[whichInGrob(agePL, "panel")]], 4, 2)
gt = gtable_add_grob(gt, ighvPL$grobs[[whichInGrob(ighvPL, "panel")]], 5, 2)
gt = gtable_add_grob(gt, mutPL$grobs[[whichInGrob(mutPL, "panel")]], 9, 2)

# add x axis
gt = gtable_add_grob(gt, mutPL$grobs[[whichInGrob(mutPL, "axis-b")]], 10, 2)

# add y axis
gt = gtable_add_grob(gt, sexPL$grobs[[whichInGrob(sexPL, "axis-l")]], 2, 1)
gt = gtable_add_grob(gt, treatPL$grobs[[whichInGrob(treatPL, "axis-l")]], 3, 1)
gt = gtable_add_grob(gt, agePL$grobs[[whichInGrob(agePL, "axis-l")]], 4, 1)
gt = gtable_add_grob(gt, ighvPL$grobs[[whichInGrob(ighvPL, "axis-l")]], 5, 1)
gt = gtable_add_grob(gt, mutPL$grobs[[whichInGrob(mutPL, "axis-l")]], 9, 1)

## ---- echo=FALSE--------------------------------------------------------------
knitr::knit_hooks$set(crop=NULL)

## ----part1, echo=FALSE, fig.path=plotDir, dev=c("png", "pdf"), fig.width=sum(wdths), fig.height=sum(hghts)----
#FIG# 1D
grid.draw(gt)

## ---- echo=FALSE--------------------------------------------------------------
knitr::knit_hooks$set(crop=knitr:::hook_pdfcrop)

## ----charLegend, echo=FALSE, fig.path=plotDir, dev=c("png", "pdf"), fig.width=10, fig.height=2----
#FIG# 1D
BloodCancerMultiOmics2017:::drawLegends(plobj=list(agePL,sexPL,treatPL,ighvPL,mutPL))

## -----------------------------------------------------------------------------
data("lpdAll")

## -----------------------------------------------------------------------------
plotTab = pData(lpdAll) %>%
  transmute(x=log10(ATPday0), y=log10(ATP48h), diff=ATP48h/ATPday0) %>%
  filter(!is.na(x))

## ---- fig.width=8, fig.height=7-----------------------------------------------
lm_eqn <- function(df){
  m <- lm(y ~ 1, df, offset = x)
  ypred <- predict(m, newdata = df)
  r2 = sum((ypred - df$y)^2)/sum((df$y - mean(df$y)) ^ 2)
  eq <- substitute(italic(y) == italic(x) + a*","~~italic(r)^2~"="~r2, 
                   list(a = format(coef(m)[1], digits = 2), 
                        r2 = format(r2, digits = 2)))
  as.character(as.expression(eq))                 
}

plotTab$ypred <- predict(lm(y~1,plotTab, offset = x), newdata = plotTab)
sca <- ggplot(plotTab, aes(x= x, y = y)) + geom_point(size=3) + 
  geom_smooth(data = plotTab, mapping = aes(x=x, y = ypred), method = lm, se = FALSE, formula = y ~ x) + 
  geom_text(x = 5.2, y = 6.2, label = lm_eqn(plotTab), parse = TRUE, size =8) + 
  xlab("log10(day0 ATP luminescence)") + ylab("log10(48h ATP luminescence)") +
  theme_bw() + theme(axis.title = element_text(size = 15, face = "bold"), 
                     axis.text = element_text(size=15), legend.position = "none") +
  coord_cartesian(xlim = c(4.6,6.3), ylim = c(4.6,6.3))

## -----------------------------------------------------------------------------
histo <- ggplot(plotTab, aes(x = diff)) + geom_histogram(col = "red", fill = "red", bins=30, alpha = 0.5) + theme_bw() +
  theme(axis.title = element_text(size = 15, face = "bold"), axis.text = element_text(size=15), legend.position = "none") +
  xlab("(48h ATP luminescence) / (day0 ATP luminescence)")

## ----ATPcount_combine, fig.path=plotDir, dev=c("png", "pdf"), fig.width=13, fig.height=6----
grid.arrange(sca, histo, ncol=2)

## -----------------------------------------------------------------------------
data("day23rep")

## -----------------------------------------------------------------------------
maxXY = 125

plottingDF = do.call(rbind, lapply(c("day2","day3"), function(day) {
  tmp = merge(
    meltWholeDF(assayData(day23rep)[[paste0(day,"rep1")]]),
    meltWholeDF(assayData(day23rep)[[paste0(day,"rep2")]]),
    by=c("X","Y"))
  colnames(tmp) = c("PatientID", "DrugID", "ViabX", "ViabY")
  tmp[,c("ViabX", "ViabY")] = tmp[,c("ViabX", "ViabY")] * 100
  tmp$Day = ifelse(day=="day2", "48h", "72h")
  tmp
}))

plottingDF$Shape =
  ifelse(plottingDF$ViabX > maxXY | plottingDF$ViabY > maxXY, "B", "A")

## -----------------------------------------------------------------------------
annotation = 
  do.call(rbind,
          tapply(1:nrow(plottingDF),
                 paste(plottingDF$PatientID,
                       plottingDF$Day, sep="_"),
                 function(idx) {
                   data.frame(X=110, Y=10,
                              Shape="A",
                              PatientID=plottingDF$PatientID[idx[1]],
                              Day=plottingDF$Day[idx[1]],
                              Cor=cor(plottingDF$ViabX[idx],
                                      plottingDF$ViabY[idx],
                                      method="pearson"))
                 }))

## ----pilotRep, fig.path=plotDir, dev=c("png", "pdf"), fig.width=7, fig.height=5----
#FIG# S31
ggplot(data=plottingDF, 
       aes(x=ifelse(ViabX>maxXY,maxXY,ViabX), y=ifelse(ViabY>maxXY,maxXY,ViabY),
           shape=Shape)) +
  facet_grid(Day ~ PatientID) + theme_bw() +
  geom_hline(yintercept=100, linetype="dashed",color="darkgrey") +
  geom_vline(xintercept=100, linetype="dashed",color="darkgrey") +
  geom_abline(intercept=0, slope=1, colour="grey") +
  geom_point(size=1.5, alpha=0.6) +
  scale_x_continuous(limits=c(0,maxXY), breaks=seq(0,maxXY,25)) +
  scale_y_continuous(limits=c(0,maxXY), breaks=seq(0,maxXY,25)) +
  xlab("% viability - replicate 1") + ylab("% viability - replicate 2") +
  coord_fixed() + expand_limits(x = 0, y = 0) +
  theme(axis.title.x=element_text(size = rel(1), vjust=-1),
        axis.title.y=element_text(size = rel(1), vjust=1),
        strip.background=element_rect(fill="gainsboro")) +
  guides(shape=FALSE, size=FALSE) +
  geom_text(data=annotation,
            aes(x=X, y=Y, label=format(Cor, digits=2), size=1.2),
            colour="maroon", hjust=0.2)

## ---- include=!exists(".standalone"), eval=!exists(".standalone")-------------
#  sessionInfo()

## ---- message=FALSE, warning=FALSE, include=FALSE-----------------------------
rm(list=ls())

## ---- message=FALSE, include=!exists(".standalone"), eval=!exists(".standalone")----
#  library("BloodCancerMultiOmics2017")
#  library("reshape2") # melt
#  library("Biobase")
#  library("dplyr")
#  library("RColorBrewer")
#  library("ggplot2")
#  library("ggdendro")
#  library("gtable")
#  library("grid")
#  library("Rtsne")
#  library("ggbeeswarm")

## ----echo=FALSE---------------------------------------------------------------
plotDir = ifelse(exists(".standalone"), "", "part02/")
if(plotDir!="") if(!file.exists(plotDir)) dir.create(plotDir)

## -----------------------------------------------------------------------------
options(stringsAsFactors=FALSE)

## -----------------------------------------------------------------------------
data("lpdAll")

## -----------------------------------------------------------------------------
#select drug screening data on patient samples
lpd <- lpdAll[fData(lpdAll)$type == "viab", pData(lpdAll)$Diagnosis != "hMNC"]
viabTab <- Biobase::exprs(lpd)
viabTab <- viabTab[,complete.cases(t(viabTab))]
viabTab <- reshape2::melt(viabTab)
viabTab$Concentration <- fData(lpd)[viabTab$Var1,"subtype"]
viabTab <- viabTab[viabTab$Concentration %in% c("1","2","3","4","5"),]
viabTab$drugName <- fData(lpd)[viabTab$Var1,"name"]
viabTab <- viabTab[order(viabTab$Concentration),]

#order drug by mean viablitity
drugOrder <- group_by(viabTab, drugName) %>%
  summarise(meanViab = mean(value)) %>%
  arrange(meanViab)
viabTab$drugName <- factor(viabTab$drugName, levels = drugOrder$drugName)

## ----ViabilityScatter_main, echo=FALSE, fig.path=plotDir, dev=c("png", "pdf"), fig.width=10, fig.height=5, warning=FALSE, out.width=560, out.height=280----
#FIG# S1

#Color for each concentration
colorCode <- rev(brewer.pal(7,"Blues")[3:7])
names(colorCode) <- unique(viabTab$Concentration)

ggplot(viabTab, aes(x=drugName,y=value, color=Concentration)) +
  geom_jitter(size=1, na.rm = TRUE, alpha=0.8, shape =16) +
  scale_color_manual(values = colorCode) +
  ylab("Viability") + ylim(c(0,1.2)) + xlab("") +
  guides(color = guide_legend(override.aes = list(size=3,alpha=1),
                              title = "concentration index")) +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
        legend.key = element_blank())

## -----------------------------------------------------------------------------
data("drpar", "patmeta", "drugs")

## -----------------------------------------------------------------------------
givePatientID4Diag = function(pts, diag=NA) {
  pts = if(is.na(diag)) {
    names(pts)
  } else {
    names(pts[patmeta[pts,"Diagnosis"]==diag])
  }
  pts
}

## -----------------------------------------------------------------------------
giveViabMatrix = function(diag, screen, chnnl) {
  data = if(screen=="main") drpar
  else print("Incorrect screen name.")
  pid = colnames(data)
  if(!is.na(diag))
    pid = pid[patmeta[pid,"Diagnosis"]==diag]
  
  return(assayData(data)[[chnnl]][,pid])
}

## -----------------------------------------------------------------------------
palette.cor1 = c(rev(brewer.pal(9, "Blues"))[1:8],
                 "white","white","white","white",brewer.pal(7, "Reds"))
palette.cor2 = c(rev(brewer.pal(9, "Blues"))[1:8],
                 "white","white","white","white",brewer.pal(7, "YlOrRd"))

## -----------------------------------------------------------------------------
main.cll.tpll = BloodCancerMultiOmics2017:::makeCorrHeatmap(
  mt=giveViabMatrix(diag="CLL", screen="main", chnnl="viaraw.4_5"),
  mt2=giveViabMatrix(diag="T-PLL", screen="main", chnnl="viaraw.4_5"),
  colsc=palette.cor2, concNo="one")

## ----main.CLL.T-PLL, echo=FALSE, fig.path=plotDir, dev=c("png", "pdf"), fig.width=main.cll.tpll[["figure"]][["width"]], fig.height=main.cll.tpll[["figure"]][["height"]]----
#FIG# 2A
grid.draw(main.cll.tpll[["figure"]][["plot"]])

## ----main.cll.tpll.legend, echo=FALSE, fig.path=plotDir, dev=c("png", "pdf"), fig.width=main.cll.tpll[["legend"]][["width"]], fig.height=main.cll.tpll[["legend"]][["height"]], out.width=300, out.height=150----
#FIG# 2A
grid.draw(main.cll.tpll[["legend"]][["plot"]])

## -----------------------------------------------------------------------------
# select the data
mtcll = as.data.frame(t(giveViabMatrix(diag="CLL",
                                       screen="main",
                                       chnnl="viaraw.4_5")))
colnames(mtcll) = drugs[colnames(mtcll),"name"]

# function which plots the scatter plot
scatdr = function(drug1, drug2, coldot, mtNEW, min){
  
  dataNEW = mtNEW[,c(drug1, drug2)]
  colnames(dataNEW) = c("A", "B")
  
  p = ggplot(data=dataNEW,  aes(A,  B)) + geom_point(size=3, col=coldot, alpha=0.8) +
    labs(x = drug1, y = drug2) + ylim(c(min, 1.35)) +  xlim(c(min, 1.35)) +
    theme(panel.background = element_blank(),
          axis.text = element_text(size = 15),
          axis.title = element_text(size = rel(1.5)),
          axis.line.x = element_line(colour = "black", size = 0.5),
          axis.line.y = element_line(colour = "black", size = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    geom_smooth(method=lm) +
    geom_text(x=1, y=min+0.1,
              label=paste0("Pearson-R = ",
                           round(cor(dataNEW$A, dataNEW$B ), 2)),
              size = 5)
  
  return(p)
}

## ----cor_scatter, echo=FALSE, fig.path=plotDir, dev=c("png", "pdf"), , fig.width=6, fig.height=4, warning=FALSE, out.width=420, out.height=280----
#FIG# 2A
scatdr("ibrutinib", "spebrutinib", coldot="deeppink1", mtNEW=mtcll, min=0.4)
scatdr("ibrutinib", "PRT062607 HCl", coldot="deeppink1", mtNEW=mtcll, min=0.4)
scatdr("ibrutinib", "idelalisib", coldot="deeppink1", mtNEW=mtcll, min=0.4)
scatdr("venetoclax", "navitoclax", coldot="goldenrod2", mtNEW=mtcll, min=0.2)
scatdr("SD51", "MIS-43", coldot="dodgerblue3", mtNEW=mtcll, min=0.2)

## -----------------------------------------------------------------------------
main.tpll = BloodCancerMultiOmics2017:::makeCorrHeatmap(
  mt=giveViabMatrix(diag="T-PLL", screen="main", chnnl="viaraw.4_5"),
  colsc=palette.cor1, concNo="one")

## ----main.T-PLL, echo=FALSE, fig.path=plotDir, dev=c("png", "pdf"), fig.width=main.tpll[["figure"]][["width"]], fig.height=main.tpll[["figure"]][["height"]]----
#FIG# S6 B
grid.draw(main.tpll[["figure"]][["plot"]])

## ----main.tpll.legend, echo=FALSE, fig.path=plotDir, dev=c("png", "pdf"), fig.width=main.tpll[["legend"]][["width"]], fig.height=main.tpll[["legend"]][["height"]], out.width=300, out.height=150----
#FIG# S6 B
grid.draw(main.tpll[["legend"]][["plot"]])

## -----------------------------------------------------------------------------
main.mcl = BloodCancerMultiOmics2017:::makeCorrHeatmap(
  mt=giveViabMatrix(diag="MCL", screen="main", chnnl="viaraw.4_5"),
  colsc=palette.cor1, concNo="one")

## ----main.MCL, echo=FALSE, fig.path=plotDir, dev=c("png", "pdf"), fig.width=main.mcl[["figure"]][["width"]], fig.height=main.mcl[["figure"]][["height"]]----
#FIG# S6 A
grid.draw(main.mcl[["figure"]][["plot"]])

## ----main.mcl.legend, echo=FALSE, fig.path=plotDir, dev=c("png", "pdf"), fig.width=main.mcl[["legend"]][["width"]], fig.height=main.mcl[["legend"]][["height"]], out.width=300, out.height=150----
#FIG# S6 A
grid.draw(main.mcl[["legend"]][["plot"]])

## -----------------------------------------------------------------------------
data(list=c("lpdAll", "conctab", "patmeta"))

## -----------------------------------------------------------------------------
#Select rows contain drug response data
lpdSub <- lpdAll[fData(lpdAll)$type == "viab",]

#Only use samples with complete values
lpdSub <- lpdSub[,complete.cases(t(Biobase::exprs(lpdSub)))]

#Transformation of the values
Biobase::exprs(lpdSub) <- log(Biobase::exprs(lpdSub))
Biobase::exprs(lpdSub) <- t(scale(t(Biobase::exprs(lpdSub))))

#annotation for drug ID
anno <- sprintf("%s(%s)",fData(lpdSub)$name,fData(lpdSub)$subtype)
names(anno) <- rownames(lpdSub)

## -----------------------------------------------------------------------------
tsneRun <- function(distMat,perplexity=10,theta=0,max_iter=5000, seed = 1000) {
  set.seed(seed)
  tsneRes <- Rtsne(distMat, perplexity = perplexity, theta = theta, 
                   max_iter = max_iter, is_distance = TRUE, dims =2)
  tsneRes <- tsneRes$Y
  rownames(tsneRes) <- labels(distMat)
  colnames(tsneRes) <- c("x","y")
  tsneRes
}

## -----------------------------------------------------------------------------
colDiagFill = c(`CLL` = "grey80",
                `U-CLL` = "grey80",
                `B-PLL`="grey80",
                `T-PLL`="#cc5352",
                `Sezary`="#cc5352",
                `PTCL-NOS`="#cc5352",
                `HCL`="#b29441",
                `HCL-V`="mediumaquamarine",
                `AML`="#addbaf",
                `MCL`="#8e65ca",
                `MZL`="#c95e9e",
                `FL`="darkorchid4",
                `LPL`="#6295cd",
                `hMNC`="pink")

colDiagBorder <- colDiagFill
colDiagBorder["U-CLL"] <- "black"

## -----------------------------------------------------------------------------
annoDiagNew <- function(patList, lpdObj = lpdSub) {
  Diagnosis <- pData(lpdObj)[patList,c("Diagnosis","IGHV Uppsala U/M")]
  DiagNew <- c()
  
  for (i in seq(1:nrow(Diagnosis))) {
    if (Diagnosis[i,1] == "CLL") {
      if (is.na(Diagnosis[i,2])) {
        DiagNew <- c(DiagNew,"CLL")
      } else if (Diagnosis[i,2] == "U") {
        DiagNew <- c(DiagNew,sprintf("%s-%s",Diagnosis[i,2],Diagnosis[i,1]))
      } else if (Diagnosis[i,2] == "M") {
        DiagNew <- c(DiagNew,"CLL")
      }
    } else DiagNew <- c(DiagNew,Diagnosis[i,1])
  }
  DiagNew
}

## -----------------------------------------------------------------------------
#prepare distance matrix
distLpd <- dist(t(Biobase::exprs(lpdSub)))

#run t-SNE
plotTab <- data.frame(tsneRun(distLpd,perplexity=25, max_iter=5000, seed=338))

#annotated patient sample
plotTab$Diagnosis <- pData(lpdSub[,rownames(plotTab)])$Diagnosis
plotTab$Diagnosis <- annoDiagNew(rownames(plotTab,lpdSub)) #consider IGHV status
plotTab$Diagnosis <- factor(plotTab$Diagnosis,levels = names(colDiagFill))

## ----tSNE, echo=FALSE, fig.path=plotDir, dev=c("png", "pdf"), fig.width=10, fig.height=8, warning=FALSE, out.width=600, out.height=480----
#FIG# 2 C
p <- (ggplot(plotTab, aes(x=x,y=y)) +
        geom_point(size=3, shape= 21, aes(col = Diagnosis, fill = Diagnosis)) +
        theme_classic() +
        theme(axis.ticks=element_line(color="black",size=0.5),
              text=element_text(size=20),
              axis.line.x = element_line(color="black",size=0.5),
              axis.line.y = element_line(color="black",size=0.5),
              legend.position="right") +
        scale_fill_manual(values = colDiagFill) +
        scale_color_manual(values = colDiagBorder) +
        xlab("Component 1") + ylab("Component 2")) +
  coord_cartesian(xlim = c(-20,20),ylim=c(-20,20))

print(p)

## -----------------------------------------------------------------------------
lpdPlot <- lpdAll[fData(lpdAll)$type == "viab",]
concList <- c()
for (drugID in rownames(fData(lpdPlot))) {
  concIndex <- as.character(fData(lpdPlot)[drugID,"subtype"])
  concSplit <- unlist(strsplit(as.character(concIndex),":"))
  ID <- substr(drugID,1,5)
  if (length(concSplit) == 1) {
    realConc <- conctab[ID,as.integer(concSplit)]
    concList <- c(concList,realConc)
  } else {
    realConc <- sprintf("%s:%s",
                        conctab[ID,as.integer(concSplit[1])],
                        conctab[ID,as.integer(concSplit[2])])
    concList <- c(concList,realConc)
  }
}

fData(lpdPlot)$concValue <- concList
lpdPlot <- lpdPlot[,complete.cases(t(Biobase::exprs(lpdPlot)))]

## -----------------------------------------------------------------------------
patDiag <- c("CLL","T-PLL","HCL","MCL")
drugID <- c("D_012_5","D_017_4","D_039_3","D_040_5","D_081_4","D_083_5")

lpdBee <- lpdPlot[drugID,pData(lpdPlot)$Diagnosis %in% patDiag]

## -----------------------------------------------------------------------------
lpdCurve <-
  lpdPlot[fData(lpdPlot)$name %in% fData(lpdBee)$name,
          pData(lpdPlot)$Diagnosis %in% patDiag]
lpdCurve <- lpdCurve[fData(lpdCurve)$subtype %in% seq(1,5),]
dataCurve <- data.frame(Biobase::exprs(lpdCurve))
dataCurve <- cbind(dataCurve,fData(lpdCurve)[,c("name","concValue")])
tabCurve <- melt(dataCurve,
                 id.vars = c("name","concValue"), variable.name = "patID")
tabCurve$Diagnosis <- factor(pData(lpdCurve[,tabCurve$patID])$Diagnosis,
                             levels = patDiag)
tabCurve$value <- tabCurve$value
tabCurve$concValue <- as.numeric(tabCurve$concValue)

# set order
tabCurve$name <- factor(tabCurve$name, levels = fData(lpdBee)$name)

#calculate the mean and mse for each drug+cencentration in different disease
tabGroup <- group_by(tabCurve,name,concValue,Diagnosis)
tabSum <- summarise(tabGroup,meanViab = mean(value))

## ----viabilityCurve, fig.path=plotDir, dev=c("png", "pdf"), fig.width=4, fig.height=3----
#FIG# 2 C
tconc = expression("Concentration [" * mu * "M]")
fmt_dcimals <- function(decimals=0){
  # return a function responpsible for formatting the 
  # axis labels with a given number of decimals 
  function(x) as.character(round(x,decimals))
}

for (drugName in unique(tabSum$name)) {
  tabDrug <- filter(tabSum, name == drugName)
  p <- (ggplot(data=tabDrug, aes(x=concValue,y=meanViab, col=Diagnosis)) +
          geom_line() + geom_point(pch=16, size=4) +
          scale_color_manual(values = colDiagFill[patDiag])
        + theme_classic() +
          theme(panel.border=element_blank(),
                axis.line.x=element_line(size=0.5,
                                         linetype="solid", colour="black"),
                axis.line.y = element_line(size = 0.5,
                                           linetype="solid", colour="black"),
                legend.position="none",
                plot.title = element_text(hjust = 0.5, size=20),
                axis.text = element_text(size=15),
                axis.title = element_text(size=20)) +
          ylab("Viability") + xlab(tconc) + ggtitle(drugName) +
          scale_x_log10(labels=fmt_dcimals(2)) +
          scale_y_continuous(limits = c(0,1.3), breaks = seq(0,1.3,0.2)))
  plot(p)
}

## ----viabilityBee, fig.path=plotDir, dev=c("png", "pdf"), fig.width=5, fig.height=10, warning=FALSE----
#FIG# 2 D
lpdDiag <- lpdAll[,pData(lpdAll)$Diagnosis %in% c("CLL", "MCL", "HCL", "T-PLL")]
dr <- c("D_012_5", "D_083_5", "D_081_3", "D_040_4", "D_039_3")

m <- data.frame(t(Biobase::exprs(lpdDiag)[dr, ]), diag=pData(lpdDiag)$Diagnosis)
m <- melt(m)
m$lable <- 1
for (i in 1:nrow(m )) {
  m[i, "lable"] <- giveDrugLabel(as.character(m[i, "variable"]), conctab, drugs)
} 


ggplot( m, aes(diag, value, color=factor(diag) ) ) +
  ylim(0,1.3) + ylab("Viability") + 
  xlab("") +
  geom_boxplot(outlier.shape = NA) +
  geom_beeswarm(cex=1.4, size=1.4,alpha=0.5, color="grey80") +
  scale_color_manual("diagnosis", values=c(colDiagFill["CLL"], colDiagFill["MCL"], 
                                           colDiagFill["HCL"], colDiagFill["T-PLL"])) +
  theme_bw() +
  theme(legend.position="right") +
  theme(
    panel.background =  element_blank(), 
    panel.grid.minor.x =  element_blank(),
    axis.text = element_text(size=15),
    axis.title = element_text(size=15),
    strip.text = element_text(size=15)
  ) +
  facet_wrap(~ lable, ncol=1) 

## ---- include=!exists(".standalone"), eval=!exists(".standalone")-------------
#  sessionInfo()

## ---- message=FALSE, warning=FALSE, include=FALSE-----------------------------
rm(list=ls())

## ---- message=FALSE, include=!exists(".standalone"), eval=!exists(".standalone")----
#  library("BloodCancerMultiOmics2017")
#  library("readxl")
#  library("dplyr")
#  library("ggplot2")
#  library("reshape2")
#  library("xtable")

## ----echo=FALSE---------------------------------------------------------------
plotDir = ifelse(exists(".standalone"), "", "part14/")
if(plotDir!="") if(!file.exists(plotDir)) dir.create(plotDir)

## ---- echo=FALSE--------------------------------------------------------------

# make trasperent 
makeTransparent = function(..., alpha=0.18) {
  
  if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")
  
  alpha = floor(255*alpha)  
  newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
  
  .makeTransparent = function(col, alpha) {
    rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
  }
  
  newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)
  
  return(newColor)
}

## ----load_affinity_data-------------------------------------------------------
# AZD7762 binding affinity constants
azd = read_excel(system.file("extdata","TargetProfiling.xlsx",
                             package="BloodCancerMultiOmics2017"), sheet = 1)

# PF477736 binding affinity constants
pf = read_excel(system.file("extdata","TargetProfiling.xlsx",
                            package="BloodCancerMultiOmics2017"), sheet = 2)

# BCR tagets Proc Natl Acad Sci U S A. 2016 May 17;113(20):5688-93
pProt = read_excel(system.file("extdata","TargetProfiling.xlsx",
                               package="BloodCancerMultiOmics2017"),sheet = 3)

## -----------------------------------------------------------------------------
p <- full_join(azd, pf )
p <- full_join(p, pProt )

pp <- p[p$BCR_effect=="Yes",]
pp <- data.frame(pp[-which(is.na(pp$BCR_effect)),])

## ----azd-pf, fig.path=plotDir, dev=c('png','pdf'), fig.width=7, fig.height=4----
#FIG# 2B
rownames(pp) <- 1:nrow(pp)
pp <- as.data.frame(pp)
pp <- melt(pp)
colnames(pp)[3] <- "Drugs"
colnames(pp)[4] <- "Score"


ggplot(pp, aes(x= reorder(gene, Score), Score, colour=Drugs ) )+ geom_point(size=3) +
  
  scale_colour_manual(values = c(makeTransparent("royalblue1", alpha = 0.75),
                                 makeTransparent("royalblue4", alpha = 0.75), 
                                 makeTransparent("brown1", alpha = 0.55),
                                 makeTransparent("brown3", alpha = 0.35)),
                      
                      breaks = c("az10", "az2", "pf10", "pf2"),
                      labels = c("AZD7762 10 µM","AZD7762 2 µM","PF477736 10 µM","PF477736 2 µM") ) +
  
  ylab("Binding affinity") +
  
  theme_bw() + geom_hline(yintercept = 0.5) +
  
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x=element_blank() )  

## ----kinobead-----------------------------------------------------------------
j <- apply(p[,c("az10", "az2", "pf10", "pf2")], 1, function (x) { min(x, na.rm=FALSE) } )

p <-  p[which(j<0.5), ]
p <-  unique(p, by = p$gene)

knitr::kable(p)

## ----kinobead_write_file, echo=FALSE, results='hide', eval=FALSE--------------
#  write(print(p), file=paste0(plotDir,"kinobead.tex"))

## ---- include=!exists(".standalone"), eval=!exists(".standalone")-------------
#  sessionInfo()

## ---- message=FALSE, warning=FALSE, include=FALSE-----------------------------
rm(list=ls())

## ----setup, message=FALSE, include=!exists(".standalone"), eval=!exists(".standalone")----
#  options(error = recover)
#  knitr::opts_chunk$set(tidy = FALSE, # cache = TRUE, autodep = TRUE,
#                        message = FALSE, error = FALSE, warning = TRUE)
#  
#  library("BloodCancerMultiOmics2017")
#  library("Biobase")
#  library("dplyr")
#  library("RColorBrewer")
#  library("dendsort")
#  library("nat")  # provides nlapply
#  library("grid")
#  library("magrittr")

## ----echo=FALSE---------------------------------------------------------------
plotDir = ifelse(exists(".standalone"), "", "part03/")
if(plotDir!="") if(!file.exists(plotDir)) dir.create(plotDir)

## ----load---------------------------------------------------------------------
data("conctab", "lpdAll", "drugs", "patmeta")
lpdCLL <- lpdAll[, lpdAll$Diagnosis=="CLL"]
lpdAll = lpdAll[, lpdAll$Diagnosis!="hMNC"]

## -----------------------------------------------------------------------------
someMatch <- function(...) {
  rv <- match(...)
  if (all(is.na(rv)))
    stop(sprintf("`match` failed to match any of the following: %s",
                 paste(x[is.na(rv)], collapse=", ")))
  rv  
}

## ----geneusage----------------------------------------------------------------
colnames(pData(lpdAll))
gu    <- pData(lpdAll)$`IGHV Uppsala gene usage`
tabgu <- sort(table(gu), decreasing = TRUE)
biggestGroups <- names(tabgu)[1:5]
gu[ is.na(match(gu, biggestGroups)) & !is.na(gu) ] <- "other"
pData(lpdAll)$`IGHV gene usage` <- factor(gu, levels = c(biggestGroups, "other"))

## ----targetedDrugNames--------------------------------------------------------
stopifnot(is.null(drugs$id))
drugs$id <- rownames(drugs)

targetedDrugNames <- c("ibrutinib", "idelalisib",  "PRT062607 HCl", 
                       "duvelisib", "spebrutinib", "selumetinib", "MK-2206",  
                       "everolimus", "encorafenib")
id1 <- safeMatch(targetedDrugNames, drugs$name)
targetedDrugs <- paste( rep(drugs[id1, "id"], each = 2), 4:5, sep="_" )

chemoDrugNames <- c("fludarabine", "doxorubicine",  "nutlin-3")
id2 <- safeMatch(chemoDrugNames, drugs$name)
chemoDrugs <- paste( rep(drugs[id2, "id"], each = 5), 3:5, sep="_" )

tzselDrugNames <- c("ibrutinib", "idelalisib", "duvelisib", "selumetinib", 
                    "AZD7762", "MK-2206", "everolimus", "venetoclax", "thapsigargin", 
                    "AT13387", "YM155", "encorafenib", "tamatinib", "ruxolitinib", 
                    "PF 477736", "fludarabine", "nutlin-3")
id3 <- safeMatch(tzselDrugNames, drugs$name)
tzselDrugs <- unlist(lapply(seq(along = tzselDrugNames), function(i)
  paste(drugs[id3[i], "id"], 
        if (tzselDrugNames[i] %in% c("fludarabine", "nutlin-3")) 2:3 else 4:5, 
        sep = "_" )))

## ----addBCR-------------------------------------------------------------------
bcrDrugs <- c("ibrutinib", "idelalisib", "PRT062607 HCl", "spebrutinib")
everolID <- drugs$id[ safeMatch("everolimus",  drugs$name)]
bcrID    <- drugs$id[ safeMatch(bcrDrugs, drugs$name)]

is_BCR  <- 
  (fData(lpdAll)$id %in% bcrID)    & (fData(lpdAll)$subtype %in% paste(4:5))
is_mTOR <- 
  (fData(lpdAll)$id %in% everolID) & (fData(lpdAll)$subtype %in% paste(4:5))

myin <- function(x, y) as.numeric( (x %in% y) & !is.na(x) )

weights1 <- data.frame(
  hclust = rep(1, nrow(lpdAll)) + 1.75 * is_mTOR, 
  score  = as.numeric( is_BCR ), 
  row.names = rownames(lpdAll))

weights2 <- data.frame(
  row.names = tzselDrugs,
  hclust = myin(drugs$target_category[id3], "B-cell receptor") * 0.3 + 0.7,
  score  = rep(1, length(tzselDrugs))) 

## ----badDrugs-----------------------------------------------------------------
badDrugs <- c(bortezomib = "D_008", `NSC 74859` = "D_025")
stopifnot(identical(drugs[ badDrugs, "name"], names(badDrugs)))

candDrugs <- rownames(lpdAll)[  
  fData(lpdAll)$type=="viab" & !(fData(lpdAll)$id %in% badDrugs) &
    fData(lpdAll)$subtype %in% paste(2:5)
]

## ----threshs------------------------------------------------------------------
thresh <- list(effectVal = 0.7, effectNum = 4, viab = 0.6, maxval = 1.1)
overallMean  <- rowMeans(Biobase::exprs(lpdAll)[candDrugs, ])
nthStrongest <- apply(Biobase::exprs(lpdAll)[candDrugs, ], 1,
                      function(x) sort(x)[thresh$effectNum])

## ----hists, fig.width = 8, fig.height = 2.6-----------------------------------
par(mfrow = c(1, 3))
hist(overallMean,  breaks = 30, col = "pink")
abline(v = thresh$viab,      col="blue")
hist(nthStrongest, breaks = 30, col = "pink")
abline(v = thresh$effectVal, col="blue")
plot(overallMean, nthStrongest)
abline(h = thresh$effectVal, v = thresh$viab, col = "blue")

## ----seldrugs1----------------------------------------------------------------
seldrugs1 <- candDrugs[ overallMean >= thresh$viab &
                          nthStrongest <= thresh$effectVal ] %>%
  union(targetedDrugs) %>% 
  union(chemoDrugs)

d1 <- Biobase::exprs(lpdAll[seldrugs1,, drop = FALSE ]) %>%
  deckel(lower = 0, upper = thresh$maxval)

d2 <- Biobase::exprs(lpdAll[tzselDrugs,, drop = FALSE ]) %>%
  deckel(lower = 0, upper = thresh$maxval)

## ----differentspread, fig.width = 4, fig.height = 4---------------------------
spreads <- sapply(list(mad = mad, `Q95-Q05` = function(x)
  diff(quantile(x, probs = c(0.05, 0.95)))), function(s) apply(d1, 1, s))
plot( spreads )

jj <- names(which( spreads[, "mad"] < 0.15 & spreads[, "Q95-Q05"] > 0.7))
jj
drugs[ stripConc(jj), "name" ]

## ----transf-------------------------------------------------------------------
medianCenter_MadScale <- function(x) {
  s <- median(x)
  (x - s) / deckel(mad(x, center = s), lower = 0.05, upper = 0.2)
}

scaleDrugResp  <- function(x) t(apply(x, 1, medianCenter_MadScale)) 

scd1 <- scaleDrugResp(d1)
scd2 <- scaleDrugResp(d2)

## ----defdisgrp----------------------------------------------------------------
sort(table(lpdAll$Diagnosis), decreasing = TRUE)

diseaseGroups <- list(
  `CLL` = c("CLL"), 
  `MCL` = c("MCL"),
  `HCL` = c("HCL", "HCL-V"),
  `other B-cell` = c("B-PLL", "MZL", "LPL", "FL"), 
  `T-cell` = c("T-PLL", "Sezary", "PTCL-NOS"),
  `myeloid` = c("AML"))
stopifnot(setequal(unlist(diseaseGroups), unique(lpdAll$Diagnosis)))

fdg <- factor(rep(NA, ncol(lpdAll)), levels = names(diseaseGroups))
for (i in names(diseaseGroups))
  fdg[ lpdAll$Diagnosis %in% diseaseGroups[[i]] ] <- i
lpdAll$`Disease Group` <- fdg

## ----matclust-----------------------------------------------------------------
matClust <- function(x, 
                     rowweights,
                     colgroups   = factor(rep("all", ncol(x))),
                     reorderrows = FALSE) {
  
  stopifnot(is.data.frame(rowweights),
            c("hclust", "score") %in% colnames(rowweights),
            !is.null(rownames(rowweights)),
            !is.null(rownames(x)), 
            all(rownames(x) %in% rownames(rowweights)),
            is.factor(colgroups), 
            !any(is.na(colgroups)), 
            length(colgroups) == ncol(x))
  
  wgt <- rowweights[ rownames(x), ]
  
  columnsClust <- function(xk) {
    score   <- -svd(xk * wgt[, "score"])$v[, 1]
    cmns    <- colSums(xk * wgt[, "score"])
    ## make sure that high score = high sensitivity  
    if (cor(score, cmns) > 0) score <- (-score)
    
    ddraw  <- as.dendrogram(hclust(dist(t(xk * wgt[, "hclust"]), 
                                        method = "euclidean"), 
                                   method = "complete"))
    dd  <- reorder(ddraw, wts = -score, agglo.FUN = mean)
    ord <- order.dendrogram(dd)
    list(dd = dd, ord = ord, score = score)
  }
  
  sp <- split(seq(along = colgroups), colgroups)
  cc <- lapply(sp, function(k) columnsClust(x[, k, drop=FALSE]))
  cidx <- unlist(lapply(seq(along = cc), function (i)
    sp[[i]][ cc[[i]]$ord ]))
  csc  <- unlist(lapply(seq(along = cc), function (i)
    cc[[i]]$score[ cc[[i]]$ord ]))
  
  rddraw <- as.dendrogram(hclust(dist(x, method = "euclidean"),
                                 method = "complete"))
  ridx  <- if (reorderrows) {
    ww <- (colgroups == "CLL")
    stopifnot(!any(is.na(ww)), any(ww))
    rowscore <- svd(t(x) * ww)$v[, 1]
    dd <- reorder(rddraw, wts = rowscore, agglo.FUN = mean)
    order.dendrogram(dd)
  } else {
    rev(order.dendrogram(dendsort(rddraw)))
  }
  
  res <- x[ridx, cidx]
  
  stopifnot(identical(dim(res), dim(x)))
  attr(res, "colgap") <- cumsum(sapply(cc, function(x) length(x$score)))
  res
}

## -----------------------------------------------------------------------------
translation = list(IGHV=c(U=0, M=1),
                   Methylation_Cluster=c(`LP-CLL`=0, `IP-CLL`=1, `HP-CLL`=2))

## ----annosamp-----------------------------------------------------------------
make_pd <- function(cn, ...) {
  df <- function(...) data.frame(..., check.names = FALSE)
  
  x <- lpdAll[, cn]
  pd <- df(
    t(Biobase::exprs(x)[    c("del17p13", "TP53", "trisomy12"), , drop = FALSE]) %>% 
      `colnames<-`(c("del 17p13", "TP53", "trisomy 12")))
  
  # pd <- df(pd,
  #    t(Biobase::exprs(x)[  c("SF3B1", "del11q22.3", "del13q14_any"),, drop = FALSE]) %>% 
  #    `colnames<-`(c("SF3B1", "del11q22.3", "del13q14"))) 
  
  pd <- df(pd,
           cbind(as.integer(Biobase::exprs(x)["KRAS",] | Biobase::exprs(x)["NRAS",])) %>% 
             `colnames<-`("KRAS | NRAS")) 
  
  pd <- df(pd,
           # IGHV = Biobase::exprs(x)["IGHV Uppsala U/M", ],
           `IGHV (%)` = cut(x[["IGHV Uppsala % SHM"]],
                            breaks = c(0, seq(92, 100, by = 2), Inf), right = FALSE),
           `Meth. Cluster` = names(translation$Methylation_Cluster)[
             someMatch(paste(Biobase::exprs(x)["Methylation_Cluster", ]),
                       translation$Methylation_Cluster)],
           `Gene usage` = x$`IGHV gene usage`)
  
  if(length(unique(x$Diagnosis)) > 1)  
    pd <- df(pd, Diagnosis = x$Diagnosis)
  
  pd <- df(pd, 
           pretreated  = ifelse(patmeta[colnames(x),"IC50beforeTreatment"],"no","yes"), 
           alive       = ifelse(patmeta[colnames(x),"died"]>0, "no", "yes"),
           sex         = factor(x$Gender))
  
  rownames(pd) <- colnames(Biobase::exprs(x))
  
  for (i in setdiff(colnames(pd), "BCR score")) {
    if (!is.factor(pd[[i]]))
      pd[[i]] <- factor(pd[[i]])
    if (any(is.na(pd[[i]]))) {
      levels(pd[[i]]) <- c(levels(pd[[i]]), "n.d.")
      pd[[i]][ is.na(pd[[i]]) ] <- "n.d."   
    }
  }
  
  pd
}

## ----annokey------------------------------------------------------------------
gucol <- rev(brewer.pal(nlevels(lpdAll$`IGHV gene usage`), "Set3")) %>% 
  `names<-`(sort(levels(lpdAll$`IGHV gene usage`)))
gucol["IGHV3-21"] <- "#E41A1C"

make_ann_colors <- function(pd) {
  bw <- c(`TRUE` = "darkblue", `FALSE` = "#ffffff")
  res <- list(
    Btk = bw, Syk = bw, PI3K = bw, MEK = bw)
  
  if ("exptbatch" %in% colnames(pd))
    res$exptbatch <- brewer.pal(nlevels(pd$exptbatch), "Set2") %>%
    `names<-`(levels(pd$exptbatch))
  
  if ("IGHV (%)" %in% colnames(pd))
    res$`IGHV (%)` <- 
    c(rev(colorRampPalette(
      brewer.pal(9, "Blues"))(nlevels(pd$`IGHV (%)`)-1)), "white") %>% 
    `names<-`(levels(pd$`IGHV (%)`))
  
  if ("CD38" %in% colnames(pd))
    res$CD38 <- colorRampPalette(
      c("blue", "yellow"))(nlevels(pd$CD38)) %>% `names<-`(levels(pd$CD38))
  
  if("Gene usage" %in% colnames(pd)) 
    res$`Gene usage` <- gucol 
  
  if("Meth. Cluster" %in% colnames(pd)) 
    res$`Meth. Cluster` <- brewer.pal(9, "Blues")[c(1, 5, 9)] %>% 
    `names<-`(names(translation$Methylation_Cluster))
  
  res <- c(res, BloodCancerMultiOmics2017:::sampleColors)   # from addons.R  
  
  if("Diagnosis" %in% colnames(pd)) 
    res$Diagnosis <- BloodCancerMultiOmics2017:::colDiagS[
      names(BloodCancerMultiOmics2017:::colDiagS) %in% levels(pd$Diagnosis) ] %>% 
    (function(x) x[order(names(x))])
  
  for(i in names(res)) {
    whnd <- which(names(res[[i]]) == "n.d.")
    if(length(whnd)==1)
      res[[i]][whnd] <- "#e0e0e0" else 
        res[[i]] <- c(res[[i]], `n.d.` = "#e0e0e0")
      stopifnot(all(pd[[i]] %in% names(res[[i]])))
  }
  
  res
}

## ----dm-----------------------------------------------------------------------
theatmap <- function(x, cellwidth = 7, cellheight = 11) {
  
  stopifnot(is.matrix(x))
  patDat <- make_pd(colnames(x))
  
  bpp <- brewer.pal(9, "Set1")
  breaks <- 2.3 * c(seq(-1, 1, length.out = 101)) %>% `names<-`(
    colorRampPalette(c(rev(brewer.pal(7, "YlOrRd")), 
                       "white", "white", "white", 
                       brewer.pal(7, "Blues")))(101))
  
  if (!is.null(attr(x, "colgap")))
    stopifnot(last(attr(x, "colgap")) == ncol(x))
  
  pheatmapwh(deckel(x, lower = first(breaks), upper = last(breaks)), 
             cluster_rows = FALSE,  
             cluster_cols = FALSE,
             gaps_col = attr(x, "colgap"),
             gaps_row = attr(x, "rowgap"),
             scale = "none",
             annotation_col = patDat,
             annotation_colors = make_ann_colors(patDat),
             color = names(breaks),
             breaks = breaks, 
             show_rownames = TRUE, 
             show_colnames = !TRUE, 
             cellwidth = cellwidth, cellheight = cellheight, 
             fontsize = 10, fontsize_row = 11, fontsize_col = 8,
             annotation_legend = TRUE, drop_levels = TRUE)
}

## ----doheatmaps---------------------------------------------------------------
clscd1 <- matClust(scd1, rowweights = weights1, 
                   colgroups = lpdAll$`Disease Group`)
clscd2 <- matClust(scd2, rowweights = weights2, 
                   colgroups = lpdAll$`Disease Group`, reorderrows = TRUE)

## ----gapsrow------------------------------------------------------------------
setGapPositions <- function(x, gapAt) {
  rg <- if (missing(gapAt)) c(0) else {
    s <- strsplit(gapAt, split = "--")
    stopifnot(all(listLen(s) == 2L))
    s <- strsplit(unlist(s), split = ":")
    spname <- drugs[safeMatch(sapply(s, `[`, 1), drugs$name), "id"]
    spconc <- as.numeric(sapply(s, `[`, 2))
    spi <- mapply(function(d, cc) {
      i <- which(cc == conctab[d, ])
      stopifnot(length(i) == 1)
      i
    }, spname, spconc)
    spdrug <- paste(spname, spi, sep = "_")
    mt <- safeMatch(spdrug, rownames(x))
    igp <- seq(1, length(mt), by = 2L)
    stopifnot(all( mt[igp] - mt[igp + 1] == 1))
    #stopifnot(all( mt[igp] - mt[igp + 1] == 1))
    c(mt[igp + 1], 0)
  }
  attr(x, "rowgap") <- rg
  x
}

clscd1 %<>% setGapPositions(gapAt = c(
  "PF 477736:10--idelalisib:10",
  "spebrutinib:2.5--PF 477736:2.5",
  "PRT062607 HCl:10--selumetinib:2.5",
  "selumetinib:10--MK-2206:2.5",
  "MK-2206:0.156--tipifarnib:10",
  "AT13387:0.039--encorafenib:10",
  "encorafenib:2.5--SD07:1.111",
  "doxorubicine:0.016--encorafenib:0.625",
  "encorafenib:0.156--rotenone:2.5",
  "SCH 900776:0.625--everolimus:0.625",
  "everolimus:10--afatinib:1.667",
  "arsenic trioxide:1--thapsigargin:5",
  "thapsigargin:0.313--fludarabine:0.156"
))

clscd2 %<>% setGapPositions(gapAt = c(
  "AT13387:0.039--everolimus:0.156",
  "everolimus:0.625--nutlin-3:10", 
  "fludarabine:10--thapsigargin:0.078",
  "thapsigargin:0.313--encorafenib:0.625", 
  "encorafenib:0.156--ruxolitinib:0.156"
))

## ----Fig3AheatmapV1, dev = c("png", "pdf"), fig.path = plotDir, fig.height = 8.0 + nrow(clscd1) * 0.12, fig.width = 31, fig.wide = TRUE----
#FIG# S8
rownames(clscd1) <- with(fData(lpdAll)[ rownames(clscd1),, drop = FALSE], 
                         paste0(drugs[id, "name"], " ", conctab[cbind(id, paste0("c", subtype))], "uM")) 
rownames(clscd1)
theatmap(clscd1)

## ----Fig3AheatmapV2, dev = c("png", "pdf"), fig.path = plotDir, fig.height = 8.0 + nrow(clscd2) * 0.12, fig.width = 31, fig.wide = TRUE----
#FIG# 3A
rownames(clscd2) <- with(fData(lpdAll)[ rownames(clscd2),, drop = FALSE], 
                         paste0(drugs[id, "name"], " ", conctab[cbind(id, paste0("c", subtype))], "uM")) 

rownames(clscd2)
theatmap(clscd2)

## ---- include=!exists(".standalone"), eval=!exists(".standalone")-------------
#  devtools::session_info()

## ---- message=FALSE, warning=FALSE, include=FALSE, eval=exists(".standalone")----
rm(list=ls())

## ---- message=FALSE, include=!exists(".standalone"), eval=!exists(".standalone")----
#  library("BloodCancerMultiOmics2017")
#  library("Biobase")
#  library("RColorBrewer")
#  library("colorspace") # hex2RGB
#  #library("ggtern") # ggtern
#  library("reshape2") # melt
#  library("limma")
#  library("pheatmap")
#  library("beeswarm")
#  library("dplyr")
#  library("ggplot2")
#  library("tibble") # as_tibble
#  library("survival")
#  library("SummarizedExperiment")
#  library("DESeq2")

## ----echo=FALSE---------------------------------------------------------------
plotDir = ifelse(exists(".standalone"), "", "part10/")
if(plotDir!="") if(!file.exists(plotDir)) dir.create(plotDir)

## -----------------------------------------------------------------------------
options(stringsAsFactors=TRUE)

## -----------------------------------------------------------------------------
data("conctab", "lpdAll", "drugs", "patmeta")

## -----------------------------------------------------------------------------
lpdCLL <- lpdAll[, lpdAll$Diagnosis=="CLL"]

## -----------------------------------------------------------------------------
makeTransparent = function(..., alpha=0.18) {
  
  if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")
  
  alpha = floor(255*alpha)  
  newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
  
  .makeTransparent = function(col, alpha) {
    rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
  }
  newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)
  
  return(newColor)
}

giveColors = function(idx, alpha=1) {
  bp = brewer.pal(12, "Paired")
  makeTransparent(
    sequential_hcl(12, h = coords(as(hex2RGB(bp[idx]), "polarLUV"))[1, "H"])[1],
    alpha=alpha)
}

## ----calculate-input----------------------------------------------------------
# calculate  (x+c)/(s+3c), (y+c)/(s+3c), (z+c)/(s+3c)
prepareTernaryData = function(lpd, targets, invDrugs) {
  
  # calculate values for ternary
  df = sapply(targets, function(tg) {
    dr = paste(invDrugs[tg], c(4,5), sep="_")
    tmp = 1-Biobase::exprs(lpd)[dr,]
    tmp = colMeans(tmp)
    pmax(tmp, 0)
  })
  df = data.frame(df, sum=rowSums(df), max=rowMax(df))
  
  tern = apply(df[,targets], 2, function(x) {
    (x+0.005) / (df$sum+3*0.005)
  })
  colnames(tern) = paste0("tern", 1:3)
  # add IGHV status
  cbind(df, tern, IGHV=patmeta[rownames(df),"IGHV"],
        treatNaive=patmeta[rownames(df),"IC50beforeTreatment"])
}

## ----plottingFunction, eval=FALSE---------------------------------------------
#  makeTernaryPlot = function(td=ternData, targets, invDrugs) {
#  
#    drn = setNames(drugs[invDrugs[targets],"name"], nm=targets)
#  
#    plot = ggtern(data=td, aes(x=tern1, y=tern2, z=tern3)) +
#          #countours
#          stat_density_tern(geom='polygon', aes(fill=..level..),
#                            position = "identity", contour=TRUE, n=400,
#                            weight = 1, base = 'identity', expand = c(1.5, 1.5)) +
#          scale_fill_gradient(low='lightblue', high='red', guide = FALSE) +
#  
#          #points
#          geom_mask() +
#          geom_point(size=35*td[,"max"],
#                     fill=ifelse(td[,"treatNaive"],"green","yellow"),
#                     color="black", shape=21) +
#  
#          #themes
#          theme_rgbw( ) +
#          theme_custom(
#            col.T=giveColors(2),
#            col.L=giveColors(10),
#            col.R=giveColors(4),
#            tern.plot.background="white", base_size = 18 ) +
#  
#      labs( x = targets[1], xarrow = drn[targets[1]],
#                y = targets[2], yarrow = drn[targets[2]],
#                z = targets[3], zarrow = drn[targets[3]] ) +
#                theme_showarrows() + theme_clockwise() +
#  
#          # lines
#          geom_Tline(Tintercept=.5, colour=giveColors(2)) +
#          geom_Lline(Lintercept=.5, colour=giveColors(10)) +
#          geom_Rline(Rintercept=.5, colour=giveColors(4))
#  
#     plot
#  }

## ---- eval=FALSE--------------------------------------------------------------
#  # RUN TERNARY
#  makeTernary = function(lpd, targets, ighv=NA) {
#  
#    # list of investigated drugs and their targets
#    invDrugs = c("PI3K"="D_003", "BTK"="D_002", "SYK"="D_166",
#                 "MTOR"="D_063", "MEK"="D_012")
#  
#    ternData = prepareTernaryData(lpd, targets, invDrugs)
#    if(!is.na(ighv)) ternData = ternData[which(ternData$IGHV==ighv),]
#  
#    print(table(ternData$treatNaive))
#    ternPlot = makeTernaryPlot(ternData, targets, invDrugs)
#  
#    ternPlot
#  }

## ----BCR-tern-CLL, fig.path=plotDir, dev=c('png','pdf'), fig.width = 7, fig.height = 7, eval=FALSE----
#  #FIG# 3B
#  makeTernary(lpdCLL, c("PI3K", "BTK", "SYK"), ighv=NA)

## ----BCR-tern-MCLL, fig.path=plotDir, dev=c('png','pdf'), fig.width = 7, fig.height = 7, eval=FALSE----
#  #FIG# 3B
#  makeTernary(lpdCLL, c("PI3K", "BTK", "SYK"), ighv="M")

## ----BCR-tern-UCLL, fig.path=plotDir, dev=c('png','pdf'), fig.width = 7, fig.height = 7, eval=FALSE----
#  #FIG# 3B
#  makeTernary(lpdCLL, c("PI3K", "BTK", "SYK"), ighv="U")

## ----MEK-mTOR-tern, fig.path=plotDir, dev=c('png','pdf'), fig.width = 7, fig.height = 7, eval=FALSE----
#  #FIG# 3BC
#  makeTernary(lpdCLL, c("MTOR", "BTK", "MEK"), ighv=NA)

## ----MEK-mTOR-tern-MCLL, fig.path=plotDir, dev=c('png','pdf'), fig.width = 7, fig.height = 7, eval=FALSE----
#  #FIG# 3BC
#  makeTernary(lpdCLL, c("MTOR", "BTK", "MEK"), ighv="M")

## ----MEK-mTOR-tern-UCLL, fig.path=plotDir, dev=c('png','pdf'), fig.width = 7, fig.height = 7, eval=FALSE----
#  #FIG# 3BC
#  makeTernary(lpdCLL, c("MTOR", "BTK", "MEK"), ighv="U")

## ----PI3K-MEK-mTOR-tern-CLL, fig.path=plotDir, dev=c('png','pdf'), fig.width = 7, fig.height = 7, eval=FALSE----
#  #FIG# S9 left
#  makeTernary(lpdCLL, c("MTOR", "PI3K", "MEK"), ighv=NA)

## ----SYK-MEK-mTOR-tern-CLL, fig.path=plotDir, dev=c('png','pdf'), fig.width = 7, fig.height = 7, eval=FALSE----
#  #FIG# S9 right
#  makeTernary(lpdCLL, c("MTOR", "SYK", "MEK"), ighv=NA)

## -----------------------------------------------------------------------------
data("exprTreat", "drugs")

## -----------------------------------------------------------------------------
e <- exprTreat
colnames( pData(e) ) <- sub( "PatientID", "Patient", colnames( pData(e) ) )
colnames( pData(e) ) <- sub( "DrugID", "Drug", colnames( pData(e) ) )
pData(e)$Drug[ is.na(pData(e)$Drug) ] <- "none"
pData(e)$Drug <- relevel( factor( pData(e)$Drug ), "none" )
pData(e)$SampleID <- colnames(e)
colnames(e) <- paste( pData(e)$Patient, pData(e)$Drug, sep=":" )

head( pData(e) )

## -----------------------------------------------------------------------------
fData(e) <- fData(e)[ , c( "ProbeID", "Entrez_Gene_ID", "Symbol", 
                           "Cytoband", "Definition" ) ]

## -----------------------------------------------------------------------------
pheatmap( cor(Biobase::exprs(e)), symm=TRUE, cluster_rows = FALSE, cluster_cols = FALSE, 
          color = colorRampPalette(c("gray10","lightpink"))(100) )

## -----------------------------------------------------------------------------
mm <- model.matrix( ~ 0 + Patient + Drug, pData(e) )
colnames(mm) <- sub( "Patient", "", colnames(mm) )
colnames(mm) <- sub( "Drug", "", colnames(mm) )
head(mm)

## -----------------------------------------------------------------------------
fit <- lmFit( e, mm )
fit <- eBayes( fit )

## -----------------------------------------------------------------------------
a <- decideTests( fit, p.value = 0.1 )
t( apply( a[ , grepl( "D_...", colnames(a) ) ], 2, 
          function(x) table( factor(x,levels=c(-1,0,1)) ) ) )

## ----overlap------------------------------------------------------------------
a <-
  sapply( levels(pData(e)$Drug)[-1], function(dr1) 
    sapply( levels(pData(e)$Drug)[-1], function(dr2) 
      
      100*( length( intersect(
        unique( topTable( fit, coef=dr1, p.value=0.1,
                          number=Inf )$`Entrez_Gene_ID` ),
        unique( topTable( fit, coef=dr2, p.value=0.1,
                          number=Inf )$`Entrez_Gene_ID` ) ) ) / 
          length(unique( topTable( fit, coef=dr1, p.value=0.1,
                                   number=Inf )$`Entrez_Gene_ID`)))
    ) 
  )

rownames(a) <-drugs[ rownames(a), "name" ]
colnames(a) <-rownames(a)
a <- a[-4, -4]
a

## ----corrGenes----------------------------------------------------------------
extractGenes = function(fit, coef) {
  tmp = topTable(fit, coef=coef, number=Inf )[, c("ProbeID", "logFC")]
  rownames(tmp) = tmp$ProbeID
  colnames(tmp)[2] = drugs[coef,"name"]
  tmp[order(rownames(tmp)),2, drop=FALSE]
}

runExtractGenes = function(fit, drs) {
  tmp = do.call(cbind, lapply(drs, function(dr) {
    extractGenes(fit, dr)
  }))
  as.matrix(tmp)
}

mt = runExtractGenes(fit, drs=c("D_002","D_003","D_012","D_063"))

## -----------------------------------------------------------------------------
mt <- cbind( mt, median=rowMedians(mt))
mt <- mt[order(mt[,"median"]), ]
mt <- mt[1:2000, ]
mt <- mt[,-5]

## -----------------------------------------------------------------------------
(mtcr = cor(mt))

## ----plot-corrGenes, fig.path=plotDir, fig.width = 4, fig.height = 4, dev = c("png", "pdf")----
#FIG# 3D
pheatmap(mtcr, cluster_rows = FALSE, cluster_cols = FALSE,
         col=colorRampPalette(c("white", "lightblue","darkblue") )(100))

## -----------------------------------------------------------------------------
data("lpdAll", "drugs")

## -----------------------------------------------------------------------------
lpdCLL <- lpdAll[ , lpdAll$Diagnosis== "CLL"]

## ----z_factor-----------------------------------------------------------------
z_factor <- qnorm(0.05, lower.tail = FALSE)
z_factor

## ----seldrugs-----------------------------------------------------------------
drugnames <- c( "ibrutinib", "everolimus", "selumetinib")
ib  <- "D_002_4:5" 
ev  <- "D_063_4:5" 
se  <- "D_012_4:5"
stopifnot(identical(fData(lpdAll)[c(ib, ev, se), "name"], drugnames))

## -----------------------------------------------------------------------------
df  <- Biobase::exprs(lpdAll)[c(ib, ev, se), lpdAll$Diagnosis=="CLL"] %>%
  t %>% as_tibble %>% `colnames<-`(drugnames)

## -----------------------------------------------------------------------------
mdf <- melt(data.frame(df))

## ----scatter2, fig.width = 7, fig.height = 3.5, fig.path=plotDir--------------
grid.arrange(ncol = 2,
             ggplot(df, aes(x = 1-ibrutinib,  y = 1-everolimus )) + geom_point(),
             ggplot(df, aes(x = 1-everolimus, y = 1-selumetinib)) + geom_point()  
)

## ----mirror-------------------------------------------------------------------
pmdf <- filter(mdf, value >= 1)
ssd  <- mean( (pmdf$value - 1) ^ 2 ) ^ 0.5
ssd

## ----dnorm--------------------------------------------------------------------
dn <- tibble(
  x = seq(min(mdf$value), max(mdf$value), length.out = 100),
  y = dnorm(x, mean = 1, sd = ssd) * 2 * nrow(pmdf) / nrow(mdf) 
)

## ----hist1, fig.width = 10, fig.height = 5, dev = c("png", "pdf"), fig.path=plotDir----
#FIG# S10 A
thresh   <- 1 - z_factor * ssd
thresh
hh <- ggplot() +
  geom_histogram(aes(x = value, y = ..density..),
                 binwidth = 0.025, data = mdf) +
  theme_minimal() +
  geom_line(aes(x = x, y = y), col = "darkblue", data = dn) +
  geom_vline(col = "red", xintercept = thresh)
hh

## ----hist2, fig.width = 10, fig.height = 5, dev = c("png", "pdf"), fig.path=plotDir----
hh + facet_grid( ~ variable)

## ----dec----------------------------------------------------------------------
thresh 

df <- mutate(df,
             group = ifelse(ibrutinib < thresh, "BTK",
                            ifelse(everolimus < thresh, "mTOR",
                                   ifelse(selumetinib < thresh, "MEK", "Non-responder")))
)

## ----scatter3, fig.width = 10, fig.height = 5, dev = c("png", "pdf"), fig.path=plotDir, warning=FALSE----
#FIG# S10 B
mycol <- c(`BTK` = "Royalblue4",
           `mTOR` = "chartreuse4",
           `MEK` = "mediumorchid4",
           `Non-responder` = "grey61")

plots <- list(
  ggplot(df, aes(x = 1-ibrutinib,  y = 1-everolimus)),
  ggplot(filter(df, group != "BTK"), aes(x = 1-selumetinib, y = 1-everolimus))
)

plots <- lapply(plots, function(x)
  x + geom_point(aes(col = group), size = 1.5) + theme_minimal() +
    coord_fixed() + 
    scale_color_manual(values = mycol) +
    geom_hline(yintercept = 1 - thresh) + 
    geom_vline(xintercept = 1- thresh) +
    ylim(-0.15, 0.32) + xlim(-0.15, 0.32) 
  
)

grid.arrange(ncol = 2, grobs  = plots)

## -----------------------------------------------------------------------------
sel = defineResponseGroups(lpd=lpdAll)

## ----co-sens-all--------------------------------------------------------------
# colors
c1 = giveColors(2, 0.5)
c2 = giveColors(4, 0.85)
c3 = giveColors(10, 0.75)

# vectors
p <- vector(); d <- vector();
pMTOR <- vector(); pBTK <- vector(); pMEK <- vector(); pNONE <- vector()
dMTOR <- vector(); dBTK <- vector(); dMEK <- vector(); dNONE <- vector()
dMTOR_NONE <- vector(); pMTOR_NONE <- vector()

# groups                               
sel$grMTOR_NONE <- ifelse(sel$group=="mTOR", "mTOR", NA)
sel$grMTOR_NONE <- ifelse(sel$group=="none", "none", sel$grMTOR_NONE)

sel$grMTOR <- ifelse(sel$group=="mTOR", "mTOR", "rest")
sel$col <- ifelse(sel$group=="mTOR", c3, "grey")
sel$grBTK  <- ifelse(sel$group=="BTK",  "BTK", "rest")
sel$col <- ifelse(sel$group=="BTK",  c1, sel$col)
sel$grMEK  <- ifelse(sel$group=="MEK",  "MEK", "rest")
sel$col <- ifelse(sel$group=="MEK",  c2, sel$col)
sel$grNONE <- ifelse(sel$group=="none", "none", "rest")



for (i in 1: max(which(fData(lpdCLL)$type=="viab"))) {
  
  fit <- aov(Biobase::exprs(lpdCLL)[i, rownames(sel)] ~ sel$group)
  p[i] <- summary(fit)[[1]][["Pr(>F)"]][1]
  
  
  calc_p = function(clmn) {
    p.adjust(
      t.test(Biobase::exprs(lpdCLL)[i, rownames(sel)] ~ sel[,clmn],
             alternative = c("two.sided") )$p.value, "BH" )
  }
  
  calc_d = function(clmn) {
    diff(
      t.test(Biobase::exprs(lpdCLL)[i, rownames(sel)] ~ sel[,clmn])$estimate,
      alternative = c("two.sided") )
  }
  
  pMTOR_NONE[i] <- calc_p("grMTOR_NONE")
  dMTOR_NONE[i] <- calc_d("grMTOR_NONE")
  
  pMTOR[i] <- calc_p("grMTOR")
  dMTOR[i] <- calc_d("grMTOR")
  
  pBTK[i] <- calc_p("grBTK")
  dBTK[i] <- calc_d("grBTK")
  
  pMEK[i] <- calc_p("grMEK")
  dMEK[i] <- calc_d("grMEK")
  
  pNONE[i] <- calc_p("grNONE")
  dNONE[i] <- calc_d("grNONE")
  
  # drugnames
  d[i] <- rownames(lpdCLL)[i]
}

## ----heatmap_mean_4_5, fig.path=plotDir, fig.width = 6, fig.height = 8, dev = c("png", "pdf")----
#FIG# 3F
#construct data frame
ps <- data.frame(drug=d, pMTOR, pBTK, pMEK, pNONE, p )
ds <- data.frame(dMTOR, dBTK, dMEK, dNONE)
rownames(ps) <- ps[,1]; rownames(ds) <- ps[,1]

# selcet only rows for singel concentrations, set non-sig to zero 
ps45 <- ps[rownames(ps)[grep(rownames(ps), pattern="_4:5")],2:5 ]
for (i in 1:4) { ps45[,i] <- ifelse(ps45[,i]<0.05, ps45[,i], 0) }

ds45 <- ds[rownames(ds)[grep(rownames(ds), pattern="_4:5")],1:4 ]
for (i in 1:4) { ds45[,i] <- ifelse(ps45[,i]<0.05, ds45[,i], 0) }

# exclude non-significant rows
selDS <- rownames(ds45)[rowSums(ps45)>0]
selPS <- rownames(ps45)[rowSums(ps45)>0]
ps45 <- ps45[selPS, ]
ds45 <- ds45[selDS, ]

groupMean = function(gr) {
  rowMeans(Biobase::exprs(lpdCLL)[rownames(ps45), rownames(sel)[sel$group==gr]])
}

MBTK <- groupMean("BTK")
MMEK <- groupMean("MEK")
MmTOR <- groupMean("mTOR")
MNONE <- groupMean("none")

# create data frame, new colnames
ms <- data.frame(BTK=MBTK, MEK=MMEK, mTOR=MmTOR, NONE=MNONE)
colnames(ms) <- c("BTK", "MEK", "mTOR", "WEAK")
rownames(ms) <- drugs[substr(selPS, 1,5), "name"]

# select rows with effect sizes group vs. rest >0.05
ms <- ms[ which(rowMax(as.matrix(ds45)) > 0.05 ) ,  ]

# exclude some drugs
ms <- ms[-c(
  which(rownames(ms) %in%
          c("everolimus", "ibrutinib", "selumetinib", "bortezomib"))),]

pheatmap(ms[, c("MEK", "BTK","mTOR", "WEAK")], cluster_cols = FALSE,
         cluster_rows =TRUE, clustering_method = "centroid",
         scale = "row",
         color=colorRampPalette(
           c(rev(brewer.pal(7, "YlOrRd")), "white", "white", "white",
             brewer.pal(7, "Blues")))(101)
)

## ----sel-beeswarm, fig.path=plotDir, fig.width = 7, fig.height = 5.5, dev = c("png", "pdf")----
#FIG# 3G
# drug label
giveDrugLabel3 <- function(drid) {
  vapply(strsplit(drid, "_"), function(x) {
    k <- paste(x[1:2], collapse="_")
    paste0(drugs[k, "name"])
  }, character(1))
}

groups = sel[,"group", drop=FALSE]
groups[which(groups=="none"), "group"] = "WEAK"

# beeswarm function
beeDrug <- function(xDrug) {
  
  par(bty="l", cex.axis=1.5)
  beeswarm(
    Biobase::exprs(lpdCLL)[xDrug, rownames(sel)] ~ groups$group,
    axes=FALSE, cex.lab=1.5, ylab="Viability", xlab="", pch = 16,
    pwcol=sel$col, cex=1,
    ylim=c(min( Biobase::exprs(lpdCLL)[xDrug, rownames(sel)] ) - 0.05, 1.2) )
  
  boxplot(Biobase::exprs(lpdCLL)[xDrug, rownames(sel)] ~ groups$group,  add = T,
          col="#0000ff22", cex.lab=2, outline = FALSE) 
  
  mtext(side=3, outer=F, line=0, 
        paste0(giveDrugLabel3(xDrug) ), cex=2) 
}


beeDrug("D_001_4:5")
beeDrug("D_081_4:5")
beeDrug("D_013_4:5")
beeDrug("D_003_4:5")
beeDrug("D_020_4:5")
beeDrug("D_165_3")

## ----surv, fig.path=plotDir, fig.width = 8, fig.height = 8, dev = c("png", "pdf")----
#FIG# S11 A
patmeta[, "group"] <- sel[rownames(patmeta), "group"]

c1n <- giveColors(2)
c2n <- giveColors(4)
c3n <- giveColors(10)
c4n <- "lightgrey"

survplot(Surv(patmeta[ , "T5"],
              patmeta[ , "treatedAfter"] == TRUE)  ~ patmeta$group ,
         lwd=3, cex.axis = 1.2, cex.lab=1.5, col=c(c1n, c2n, c3n, c4n), 
         data = patmeta,  
         legend.pos = 'bottomleft', stitle = 'Drug response', 
         xlab = 'Time (years)', ylab = 'Patients without treatment (%)',
)

## -----------------------------------------------------------------------------
data(lpdAll)

## ----selectCLL----------------------------------------------------------------
lpdCLL <- lpdAll[ , lpdAll$Diagnosis== "CLL"]

## -----------------------------------------------------------------------------
sel = defineResponseGroups(lpd=lpdAll)

## -----------------------------------------------------------------------------
genes <- data.frame(
  t(Biobase::exprs(lpdCLL)[fData(lpdCLL)$type %in% c("gen", "IGHV"), rownames(sel)]),
  group = factor(sel$group)
)

genes <- genes[!(is.na(rownames(genes))), ]
colnames(genes) %<>%
  sub("del13q14_any", "del13q14", .) %>%
  sub("IGHV.Uppsala.U.M", "IGHV", .)

Nmut = rowSums(genes[, colnames(genes) != "group"], na.rm = TRUE)

mf <- sapply(c("BTK", "MEK", "mTOR", "none"), function(i) 
  mean(Nmut[genes$group==i]) 
)

## ----bp, fig.width = 4, fig.height = 3.5--------------------------------------
barplot(mf, ylab="Total number of mutations/CNVs per patient", col="darkgreen")

## ----mutationTests------------------------------------------------------------
mutsUse <- setdiff( colnames(genes), "group" )
mutsUse <- mutsUse[ colSums(genes[, mutsUse], na.rm = TRUE) >= 4 ]
mutationTests <- lapply(mutsUse, function(m) {
  tibble(
    mutation = m, 
    p = fisher.test(genes[, m], genes$group)$p.value)
}) %>% bind_rows  %>% mutate(pBH = p.adjust(p, "BH")) %>% arrange(p)

mutationTests <- mutationTests %>% filter(pBH < 0.1)

## ----numMutationTests, fig.width = 3, fig.height = 2--------------------------
nrow(mutationTests)

## ----groupTests---------------------------------------------------------------
groupTests <- lapply(unique(genes$group), function(g) {
  tibble(
    group = g, 
    p = fisher.test(
      colSums(genes[genes$group == g, mutsUse], na.rm=TRUE),
      colSums(genes[genes$group != g, mutsUse], na.rm=TRUE),
      simulate.p.value = TRUE, B = 10000)$p.value)
}) %>% bind_rows %>% arrange(p)

groupTests 

## ----ft-----------------------------------------------------------------------
fisher.genes <- function(g, ref) {
  stopifnot(length(ref) == 1)
  ggg <- ifelse(g$group == ref, ref, "other")
  idx <- which(colnames(g) != "group") 
  
  lapply(idx, function(i)
    if (sum(g[, i], na.rm = TRUE) > 2) {
      ft  <- fisher.test(ggg, g[, i])
      tibble(
        p    = ft$p.value,
        es   = ft$estimate, 
        prop = sum((ggg == ref) & !is.na(g[, i]), na.rm = TRUE),
        mut1 = sum((ggg == ref) & (g[, i] != 0),  na.rm = TRUE),
        gene = colnames(g)[i])
    }  else {
      tibble(p = 1, es = 1, prop = 0, mut1 = 0, gene = colnames(g)[i])
    }
  ) %>% bind_rows
}

## -----------------------------------------------------------------------------
pMTOR <- fisher.genes(genes, ref="mTOR") 
pBTK  <- fisher.genes(genes, ref="BTK") 
pMEK  <- fisher.genes(genes, ref="MEK") 
pNONE <- fisher.genes(genes, ref="none")

p    <- cbind(pBTK$p,    pMEK$p,    pMTOR$p,    pNONE$p)
es   <- cbind(pBTK$es,   pMEK$es,   pMTOR$es,   pNONE$es)
prop <- cbind(pBTK$prop, pMEK$prop, pMTOR$prop, pNONE$prop)
mut1 <- cbind(pBTK$mut1, pMEK$mut1, pMTOR$mut1, pNONE$mut1)

## ----prepmat------------------------------------------------------------------
p <- ifelse(p < 0.05, 1, 0) 
p <- ifelse(es <= 1, p, -p) 
rownames(p) <- pMTOR$gene
colnames(p) <- c("BTK", "MEK", "mTOR", "NONE")

pM <- p[rowSums(abs(p))!=0, ]
propM <- prop[rowSums(abs(p))!=0, ]

## -----------------------------------------------------------------------------
N <- cbind( paste0(mut1[,1],"/",prop[,1] ),
            paste0(mut1[,2],"/",prop[,2] ), 
            paste0(mut1[,3],"/",prop[,3] ),
            paste0(mut1[,4],"/",prop[,4] )
)

rownames(N) <- rownames(p)

## ----heatmap2_1, fig.path=plotDir, fig.width = 7, fig.height = 7, dev = c("png", "pdf")----
#FIG# S11 B
mutationTests <-
  mutationTests[which(!(mutationTests$mutation %in%
                          c("del13q14_bi", "del13q14_mono"))), ]
pMA <- p[mutationTests$mutation, ]
pMA

pheatmap(pMA, cluster_cols = FALSE, 
         cluster_rows = FALSE, legend=TRUE, annotation_legend = FALSE, 
         color = c("red", "white", "lightblue"),
         display_numbers = N[ rownames(pMA), ]
)

## -----------------------------------------------------------------------------
data("dds")

sel = defineResponseGroups(lpd=lpdAll)
sel$group = gsub("none","weak", sel$group)

## -----------------------------------------------------------------------------
# select patients with CLL in the main screen data 
colnames(dds) <- colData(dds)$PatID
pat <- intersect(colnames(lpdCLL), colnames(dds))
dds_CLL <- dds[,which(colData(dds)$PatID %in% pat)]

# add group label 
colData(dds_CLL)$group <- factor(sel[colnames(dds_CLL), "group"])
colData(dds_CLL)$IGHV = factor(patmeta[colnames(dds_CLL),"IGHV"])

## -----------------------------------------------------------------------------
vsd <- varianceStabilizingTransformation( assay(dds_CLL) )
colnames(vsd) = colData(dds_CLL)$PatID
rowVariance <- setNames(rowVars(vsd), nm=rownames(vsd))
sortedVar <- sort(rowVariance, decreasing=TRUE)
mostVariedGenes <- sortedVar[1:10000]
dds_CLL <- dds_CLL[names(mostVariedGenes), ]

## -----------------------------------------------------------------------------
cb <- combn(unique(colData(dds_CLL)$group), 2)

gg <- list(); ggM <- list(); ggU <- list()
DE <- function(ighv) {
  for (i in 1:ncol(cb)) {
    dds_sel <- dds_CLL[,which(colData(dds_CLL)$IGHV %in% ighv)]
    dds_sel <- dds_sel[,which(colData(dds_sel)$group %in% cb[,i])]
    design(dds_sel) = ~group
    dds_sel$group <- droplevels(dds_sel$group)
    dds_sel$group <- relevel(dds_sel$group, ref=as.character(cb[2,i]) )
    dds_sel <- DESeq(dds_sel)
    res <- results(dds_sel)
    gg[[i]] <- res[order(res$padj, decreasing = FALSE), ]
    names(gg)[i] <- paste0(cb[1,i], "_", cb[2,i])
  }
  return(gg)
}

## ----eval=FALSE---------------------------------------------------------------
#  ggM <- DE(ighv="M")
#  ggU <- DE(ighv="U")
#  gg <- DE(ighv=c("M", "U"))

## ----echo=FALSE, eval=FALSE---------------------------------------------------
#  save(ggM, ggU, gg, file=paste0(plotDir,"gexGroups.RData"))

## ----echo=TRUE, eval=TRUE-----------------------------------------------------
load(system.file("extdata","gexGroups.RData", package="BloodCancerMultiOmics2017"))

## ---- eval=FALSE--------------------------------------------------------------
#  library("biomaRt")
#  # extract all ensembl ids
#  allGenes = unique(unlist(lapply(gg, function(x) rownames(x))))
#  # get gene ids for ensembl ids
#  genSymbols = getBM(filters="ensembl_gene_id",
#                     attributes=c("ensembl_gene_id", "hgnc_symbol"),
#                     values=allGenes, mart=mart)
#  # select first id if more than one is present
#  genSymbols = genSymbols[!duplicated(genSymbols[,"ensembl_gene_id"]),]
#  # set rownames to ens id
#  rownames(genSymbols) = genSymbols[,"ensembl_gene_id"]

## ----echo=TRUE, eval=TRUE-----------------------------------------------------
load(system.file("extdata","genSymbols.RData", package="BloodCancerMultiOmics2017"))

## ----IL10, fig.width=6.5, fig.height=6.5, dev = c("png", "pdf"), fig.path=plotDir----
#FIG# S14
gen="ENSG00000136634" #IL10
drug   <- "D_063_4:5" 

patsel <- intersect( rownames(sel)[sel$group %in% c("mTOR")], colnames(vsd) )

c <- cor.test( Biobase::exprs(lpdCLL)[drug, patsel], vsd[gen, patsel] )

# get hgnc_symbol for gen
# mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
# genSym = getBM(filters="ensembl_gene_id", attributes="hgnc_symbol",
# values=gen, mart=mart)
genSym = genSymbols[gen, "hgnc_symbol"]

plot(vsd[gen, patsel], Biobase::exprs(lpdCLL)[drug, patsel],
     xlab=paste0(genSym, " expression"),
     ylab="Viability (everolimus)", pch=19, ylim=c(0.70, 0.92), col="purple",
     main = paste0("mTOR-group", "\n cor = ", round(c$estimate, 3),
                   ", p = ", signif(c$p.value,2 )),
     cex=1.2)
abline(lm(Biobase::exprs(lpdCLL)[drug, patsel] ~ vsd[gen, patsel]))

## -----------------------------------------------------------------------------
c1 = giveColors(2, 0.4)
c2 = giveColors(4, 0.7)
c3 = giveColors(10, 0.6)

## -----------------------------------------------------------------------------
sigEx <- function(real) {
  
  ggsig = lapply(real, function(x) {
    x = data.frame(x)
    x = x[which(!(is.na(x$padj))),]
    x = x[x$padj<0.1,]
    x = x[order(x$padj, decreasing = TRUE),]
    x = data.frame(x[ ,c("padj","log2FoldChange")], ensg=rownames(x) )
    x$hgnc1 = genSymbols[rownames(x), "hgnc_symbol"]
    x$hgnc2 = ifelse(x$hgnc1=="" | x$hgnc1=="T" | is.na(x$hgnc1),
                     as.character(x$ensg), x$hgnc1)
    x[-(grep(x[,"hgnc2"], pattern="IG")),]
  })
  return(ggsig)
}

## -----------------------------------------------------------------------------
barplot1 <- function(real, tit) { 
  
  # process real diff genes
  sigExPlus = sigEx(real)
  ng <- lapply(sigExPlus, function(x){ cbind(
    up=nrow(x[x$log2FoldChange>0, ]),
    dn=nrow(x[x$log2FoldChange<0, ]) ) } ) 
  ng = melt(ng)
  
  p <- ggplot(ng, aes(reorder(L1, -value)), ylim(-500:500)) + 
    geom_bar(data = ng, aes(y =  value, fill=Var2), stat="identity",
             position=position_dodge() ) +
    scale_fill_brewer(palette="Paired", direction = -1,
                      labels = c("up", "down")) +
    ggtitle(label=tit) +
    
    geom_hline(yintercept = 0,colour = "grey90") + 
    theme(
      panel.grid.minor = element_blank(),
      
      axis.title.x = element_blank(),
      axis.text.x = element_text(size=14, angle = 60, hjust = 1),
      axis.ticks.x = element_blank(),
      
      axis.title.y  = element_blank(),
      axis.text.y = element_text(size=14, colour="black"),
      axis.line.y = element_line(colour = "black",
                                 size = 0.5, linetype = "solid"),
      
      legend.key = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white"),
      legend.title = element_blank(),
      legend.text = element_text(size=14, colour="black"),
      
      panel.background = element_rect(fill = "white", color="white")
    ) 
  
  plot(p) 
}

## ----deGroups, fig.width = 10, fig.height = 7, dev = c("png", "pdf"), fig.path=plotDir----
#FIG# S11 C
barplot1(real=gg, tit="")

## -----------------------------------------------------------------------------
# beeswarm funtion
beefun <- function(df, sym) {
  par(bty="l", cex.axis=1.5)
  beeswarm(df$x ~ df$y, axes=FALSE, cex.lab=1.5, col="grey", ylab=sym, xlab="",
           pch = 16, pwcol=sel[colnames(vsd),"col"], cex=1.3)
  
  boxplot(df$x ~ df$y, add = T, col="#0000ff22", cex.lab=1.5) 
}

## ----groups-------------------------------------------------------------------
sel = defineResponseGroups(lpdCLL)
sel[,1:3] = 1-sel[,1:3]
sel$IGHV = pData(lpdCLL)[rownames(sel), "IGHV Uppsala U/M"]

c1 = giveColors(2, 0.5)
c2 = giveColors(4, 0.85)
c3 = giveColors(10, 0.75)

# add colors
sel$col <- ifelse(sel$group=="mTOR", c3, "grey")
sel$col <- ifelse(sel$group=="BTK",  c1, sel$col)
sel$col <- ifelse(sel$group=="MEK",  c2, sel$col)

## -----------------------------------------------------------------------------
cytokines <- c("CXCL2","TGFB1","CCL2","IL2","IL12B","IL4","IL6","IL10","CXCL8",
               "TNF")
cyEN =  sapply(cytokines, function(i) {
  genSymbols[which(genSymbols$hgnc_symbol==i)[1],"ensembl_gene_id"]
})

makeEmpty = function() {
  data.frame(matrix(ncol=3, nrow=length(cyEN),
                    dimnames=list(names(cyEN), c("BTK", "MEK", "mTOR"))) )
}
p = makeEmpty()
ef = makeEmpty()

## ----CYTOKINES, fig.path=plotDir, fig.width=7, fig.height=5.5, dev=c("png", "pdf")----
for (i in 1:length(cyEN) ) {
  
  geneID <- cyEN[i]
  df <- data.frame(x=vsd[geneID,  ], y=sel[colnames(vsd) ,"group"])
  df$y <- as.factor(df$y)
  
  beefun(df, sym=names(geneID))
  
  df <- within(df, y <- relevel(y, ref = "none"))
  fit <- lm(x ~y, data=df)
  
  p[i,] <- summary(fit)$coefficients[ 2:4, "Pr(>|t|)"]
  
  abtk = mean(df[df$y=="BTK", "x"], na.rm=TRUE) - mean(df[df$y=="none", "x"],
                                                       na.rm=TRUE)
  amek = mean(df[df$y=="MEK", "x"], na.rm=TRUE) - mean(df[df$y=="none", "x"],
                                                       na.rm=TRUE)
  amtor= mean(df[df$y=="mTOR", "x"], na.rm=TRUE) - mean(df[df$y=="none", "x"],
                                                        na.rm=TRUE)
  
  ef[i,] <-  c(as.numeric(abtk), as.numeric(amek), as.numeric(amtor)) 
  
  mtext( paste(    "pBTK=", summary(fit)$coefficients[ 2, "Pr(>|t|)"], 
                   "\npMEK=", summary(fit)$coefficients[ 3, "Pr(>|t|)"], 
                   "\npMTOR=", summary(fit)$coefficients[ 4, "Pr(>|t|)"], 
                   side=3 ))
}

## ----CYTOKINES-summary, fig.path=plotDir, fig.width=4.5, fig.height=3.5, dev = c("png", "pdf")----
#FIG# S11 D
# log p-values
plog <- apply(p, 2, function(x){-log(x)} )
plog_m <- melt(as.matrix(plog))
ef_m <- melt(as.matrix(ef))

# introduce effect direction
plog_m$value <- ifelse(ef_m$value>0, plog_m$value, -plog_m$value)

rownames(plog_m) <- 1:nrow(plog_m)

# fdr
fdrmin = min( p.adjust(p$mTOR, "fdr") )

### plot ####
colnames(plog_m) <- c("cytokine", "group", "p")

lev = names(sort(tapply(plog_m$p, plog_m$cytokine, function(p) min(p))))

plog_m$cytokine <- factor(plog_m$cytokine, levels=lev)


ggplot(data=plog_m, mapping=aes(x=cytokine, y=p, color=group)) + 
  scale_colour_manual(values = c(c1, c2, c3)) +
  geom_point( size=3 ) + 
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.9),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        axis.line.x=element_line(),
        axis.line.y=element_line(),
        legend.position="none"
  ) +
  scale_y_continuous(name="-log(p-value)", breaks=seq(-3,7.5,2),
                     limits=c(-3,7.5)) +
  xlab("") +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = -log(0.004588897), color="purple", linetype="dashed") +
  geom_hline(yintercept = (-log(0.05)), color="grey", linetype="dashed") 

## ---- include=!exists(".standalone"), eval=!exists(".standalone")-------------
#  sessionInfo()

## ---- message=FALSE, warning=FALSE, include=FALSE-----------------------------
rm(list=ls())

## ---- message=FALSE, include=!exists(".standalone"), eval=!exists(".standalone")----
#  library("BloodCancerMultiOmics2017")
#  library("ggplot2")
#  library("dplyr")
#  library("gridExtra")

## ----echo=FALSE---------------------------------------------------------------
plotDir = ifelse(exists(".standalone"), "", "part16/")
if(plotDir!="") if(!file.exists(plotDir)) dir.create(plotDir)

## -----------------------------------------------------------------------------
data("cytokineViab")

## ----cytokinesC0, fig.width=8, fig.height=6, warning=FALSE, fig.path=plotDir, dev=c("png", "pdf")----
cond <- c("IL-2", "IL-4", "IL-10", "IL-21", "LPS", "IgM")

for (i in cond){  
  plot = ggplot(
    filter(cytokineViab, Duplicate%in%c("1"), Stimulation==i, Timepoint=="48h"),
    aes(x=as.factor(Cytokine_Concentration2), y=Normalized_DMSO, colour=mtor,
        group=interaction(Patient))) +
    ylab("viability") + xlab("c(stimulation)") + ylim(c(0.8, 1.4)) +
    geom_line() + geom_point() + ggtitle(i) + theme_bw() + guides(color="none")
  
  assign(paste0("p",i), plot)
}

grid.arrange(`pIL-2`,`pIL-10`,`pIL-4`,`pIL-21`,pLPS, pIgM, nrow=2)

## ---- include=!exists(".standalone"), eval=!exists(".standalone")-------------
#  sessionInfo()

## ---- message=FALSE, warning=FALSE, include=FALSE-----------------------------
rm(list=ls())

## ---- message=FALSE, include=!exists(".standalone"), eval=!exists(".standalone")----
#  library("BloodCancerMultiOmics2017")
#  library("DESeq2")
#  library("reshape2")
#  library("dplyr")
#  library("tibble")
#  library("Biobase")
#  library("SummarizedExperiment")
#  library("genefilter")
#  library("piano") # loadGSC
#  library("ggplot2")
#  library("gtable")
#  library("grid")

## ----echo=FALSE---------------------------------------------------------------
plotDir = ifelse(exists(".standalone"), "", "part07/")
if(plotDir!="") if(!file.exists(plotDir)) dir.create(plotDir)

## -----------------------------------------------------------------------------
options(stringsAsFactors=FALSE)

## -----------------------------------------------------------------------------
data(list=c("dds", "lpdAll"))

gmts = list(H=system.file("extdata","h.all.v5.1.symbols.gmt",
                          package="BloodCancerMultiOmics2017"),
            C6=system.file("extdata","c6.all.v5.1.symbols.gmt",
                           package="BloodCancerMultiOmics2017"))

## -----------------------------------------------------------------------------
patGroup = defineResponseGroups(lpd=lpdAll)

## ----subset-------------------------------------------------------------------
lpdCLL <- lpdAll[fData(lpdAll)$type=="viab",
                 pData(lpdAll)$Diagnosis %in% c("CLL")]

ddsCLL <- dds[,colData(dds)$PatID %in% colnames(lpdCLL)]

## ----pat_group----------------------------------------------------------------
ddsCLL <- ddsCLL[,colData(ddsCLL)$PatID %in% rownames(patGroup)]

#remove genes without gene symbol annotations
ddsCLL <- ddsCLL[!is.na(rowData(ddsCLL)$symbol),]
ddsCLL <- ddsCLL[rowData(ddsCLL)$symbol != "",]

#add drug sensitivity annotations to coldata
colData(ddsCLL)$group <- factor(patGroup[colData(ddsCLL)$PatID, "group"])

## ----remove_low_counts--------------------------------------------------------
#only keep genes that have counts higher than 10 in any sample
keep <- apply(counts(ddsCLL), 1, function(x) any(x >= 10)) 
ddsCLL <- ddsCLL[keep,]
dim(ddsCLL)

## ----filtering----------------------------------------------------------------
ddsCLL <- estimateSizeFactors(ddsCLL)
sds <- rowSds(counts(ddsCLL, normalized = TRUE))
sh <- shorth(sds)
ddsCLL <- ddsCLL[sds >= sh,]

## ---- cache=TRUE--------------------------------------------------------------
ddsCLL.norm <- varianceStabilizingTransformation(ddsCLL)

## ----DE, cache=TRUE-----------------------------------------------------------
DEres <- list()
design(ddsCLL) <- ~ group

rnaRaw <- DESeq(ddsCLL, betaPrior = FALSE)

#extract results for different comparisons
# responders versus weak-responders
DEres[["BTKnone"]] <- results(rnaRaw, contrast = c("group","BTK","none"))
DEres[["MEKnone"]] <- results(rnaRaw, contrast = c("group","MEK","none"))
DEres[["mTORnone"]] <- results(rnaRaw, contrast = c("group","mTOR","none"))

## -----------------------------------------------------------------------------
pCut = 0.05

## ----gsea_function------------------------------------------------------------
runGSEA <- function(inputTab, gmtFile, GSAmethod="gsea", nPerm=1000){
  inGMT <- loadGSC(gmtFile,type="gmt")
  #re-rank by score
  rankTab <- inputTab[order(inputTab[,1],decreasing = TRUE),,drop=FALSE] 
  if (GSAmethod == "gsea"){
    #readin geneset database
    #GSEA analysis
    res <- runGSA(geneLevelStats = rankTab,
                  geneSetStat = GSAmethod,
                  adjMethod = "fdr", gsc=inGMT,
                  signifMethod = 'geneSampling', nPerm = nPerm)
    GSAsummaryTable(res)
  } else if (GSAmethod == "page"){
    res <- runGSA(geneLevelStats = rankTab,
                  geneSetStat = GSAmethod,
                  adjMethod = "fdr", gsc=inGMT,
                  signifMethod = 'nullDist')
    GSAsummaryTable(res)
  }
}

## -----------------------------------------------------------------------------
runGSE = function(gmt) {
  
  Res <- list()
  for (i in names(DEres)) {
    dataTab <- data.frame(DEres[[i]])
    dataTab$ID <- rownames(dataTab)
    #filter using pvalues
    dataTab <- filter(dataTab, pvalue <= pCut) %>%
      arrange(pvalue) %>% 
      mutate(Symbol = rowData(ddsCLL[ID,])$symbol)
    dataTab <- dataTab[!duplicated(dataTab$Symbol),]
    statTab <- data.frame(row.names = dataTab$Symbol, stat = dataTab$stat)
    resTab <- runGSEA(inputTab=statTab, gmtFile=gmt, GSAmethod="page")
    Res[[i]] <- arrange(resTab,desc(`Stat (dist.dir)`))
  }
  Res
}

## ----get_genes_function-------------------------------------------------------
getGenes <- function(inputTab, gmtFile){
  geneList <- loadGSC(gmtFile,type="gmt")$gsc
  enrichedUp <- lapply(geneList, function(x) 
    intersect(rownames(inputTab[inputTab[,1] >0,,drop=FALSE]),x))
  enrichedDown <- lapply(geneList, function(x)
    intersect(rownames(inputTab[inputTab[,1] <0,,drop=FALSE]),x))
  return(list(up=enrichedUp, down=enrichedDown))
}

## ----intersection_heatmap-----------------------------------------------------
plotSetHeatmap <- 
  function(geneTab, enrichTab, topN, gmtFile, tittle="",
           asterixList = NULL, anno=FALSE) {
    
    if (nrow(enrichTab) < topN) topN <- nrow(enrichTab)
    enrichTab <- enrichTab[seq(1,topN),]
    
    geneList <- getGenes(geneTab,gmtFile)
    
    geneList$up <- geneList$up[enrichTab[,1]]
    geneList$down <- geneList$down[enrichTab[,1]]
    
    #form a table 
    allGenes <- unique(c(unlist(geneList$up),unlist(geneList$down)))
    allSets <- unique(c(names(geneList$up),names(geneList$down)))
    plotTable <- matrix(data=NA,ncol = length(allSets),
                        nrow = length(allGenes),
                        dimnames = list(allGenes,allSets))
    for (setName in names(geneList$up)) {
      plotTable[geneList$up[[setName]],setName] <- 1
    }
    for (setName in names(geneList$down)) {
      plotTable[geneList$down[[setName]],setName] <- -1
    }
    
    if(is.null(asterixList)) {
      #if no correlation table specified, order by the number of
      # significant gene
      geneOrder <- rev(
        rownames(plotTable[order(rowSums(plotTable, na.rm = TRUE),
                                 decreasing =FALSE),]))
    } else {
      #otherwise, order by the p value of correlation
      asterixList <- arrange(asterixList, p)
      geneOrder <- filter(
        asterixList, symbol %in% rownames(plotTable))$symbol
      geneOrder <- c(
        geneOrder, rownames(plotTable)[! rownames(plotTable) %in% geneOrder])
    }
    
    plotTable <- melt(plotTable)
    colnames(plotTable) <- c("gene","set","value")
    plotTable$gene <- as.character(plotTable$gene)
    
    if(!is.null(asterixList)) {
      #add + if gene is positivily correlated with sensitivity, else add "-"
      plotTable$ifSig <- asterixList[
        match(plotTable$gene, asterixList$symbol),]$coef
      plotTable <- mutate(plotTable, ifSig =
                            ifelse(is.na(ifSig) | is.na(value), "",
                                   ifelse(ifSig > 0, "-", "+")))
    }
    plotTable$value <- replace(plotTable$value,
                               plotTable$value %in% c(1), "Up")
    plotTable$value <- replace(plotTable$value,
                               plotTable$value %in%  c(-1), "Down")
    
    allSymbols <- plotTable$gene
    
    geneSymbol <- geneOrder
    
    if (anno) { #if add functional annotations in addition to gene names
      annoTab <- tibble(symbol = rowData(ddsCLL)$symbol, 
                        anno = sapply(rowData(ddsCLL)$description,
                                      function(x) unlist(strsplit(x,"[[]"))[1]))
      annoTab <- annoTab[!duplicated(annoTab$symbol),]
      annoTab$combine <- sprintf("%s (%s)",annoTab$symbol, annoTab$anno)
      plotTable$gene <- annoTab[match(plotTable$gene,annoTab$symbol),]$combine
      geneOrder <- annoTab[match(geneOrder,annoTab$symbol),]$combine
      geneOrder <- rev(geneOrder)
    }
    
    plotTable$gene <- factor(plotTable$gene, levels =geneOrder)
    plotTable$set <- factor(plotTable$set, levels = enrichTab[,1])
    
    
    g <- ggplot(plotTable, aes(x=set, y = gene)) +
      geom_tile(aes(fill=value), color = "black") +
      scale_fill_manual(values = c("Up"="red","Down"="blue")) +
      xlab("") + ylab("") + theme_classic() +
      theme(axis.text.x=element_text(size=7, angle = 60, hjust = 0),
            axis.text.y=element_text(size=7),
            axis.ticks = element_line(color="white"),
            axis.line = element_line(color="white"),
            legend.position = "none") +
      scale_x_discrete(position = "top") +
      scale_y_discrete(position = "right")
    
    if(!is.null(asterixList)) {
      g <- g + geom_text(aes(label = ifSig), vjust =0.40)
    }
    
    # construct the gtable
    wdths = c(0.05, 0.25*length(levels(plotTable$set)), 5)
    hghts = c(2.8, 0.1*length(levels(plotTable$gene)), 0.05)
    gt = gtable(widths=unit(wdths, "in"), heights=unit(hghts, "in"))
    ## make grobs
    ggr = ggplotGrob(g)
    ## fill in the gtable
    gt = gtable_add_grob(gt, gtable_filter(ggr, "panel"), 2, 2)
    gt = gtable_add_grob(gt, ggr$grobs[[5]], 1, 2) # top axis
    gt = gtable_add_grob(gt, ggr$grobs[[9]], 2, 3) # right axis
    
    return(list(list(plot=gt,
                     width=sum(wdths),
                     height=sum(hghts),
                     genes=geneSymbol)))
  }    

## -----------------------------------------------------------------------------
statTab = setNames(lapply(c("mTORnone","BTKnone","MEKnone"), function(gr) {
  dataTab <- data.frame(DEres[[gr]])
  dataTab$ID <- rownames(dataTab)
  #filter using pvalues
  dataTab <- filter(dataTab, pvalue <= pCut) %>%
    arrange(pvalue) %>%
    mutate(Symbol = rowData(ddsCLL[ID,])$symbol) %>%
    filter(log2FoldChange > 0)
  dataTab <- dataTab[!duplicated(dataTab$Symbol),]
  data.frame(row.names = dataTab$Symbol, stat = dataTab$stat)
}), nm=c("mTORnone","BTKnone","MEKnone"))

## ----run_GSE_hallmark, warning=FALSE, message= FALSE--------------------------
hallmarkRes = runGSE(gmt=gmts[["H"]])

## ----run_GSE_C6, warning=FALSE, message= FALSE--------------------------------
c6Res = runGSE(gmt=gmts[["C6"]])

## -----------------------------------------------------------------------------
ddsCLL.mTOR <- ddsCLL.norm[,ddsCLL.norm$group %in% "mTOR"]
viabMTOR <- Biobase::exprs(lpdCLL["D_063_4:5", ddsCLL.mTOR$PatID])[1,]
stopifnot(all(ddsCLL.mTOR$PatID == colnames(viabMTOR)))  

## -----------------------------------------------------------------------------
#only keep genes that have counts higher than 10 in any sample
keep <- apply(assay(ddsCLL.mTOR), 1, function(x) any(x >= 10)) 
ddsCLL.mTOR <- ddsCLL.mTOR[keep,]
dim(ddsCLL.mTOR)

## -----------------------------------------------------------------------------
tmp = do.call(rbind, lapply(1:nrow(ddsCLL.mTOR), function(i) {
  res = cor.test(viabMTOR, assay(ddsCLL.mTOR[i,])[1,], method = "pearson")
  data.frame(coef=unname(res$estimate), p=res$p.value)
}))

corResult <- tibble(ID = rownames(ddsCLL.mTOR), 
                    symbol = rowData(ddsCLL.mTOR)$symbol,
                    coef = tmp$coef,
                    p = tmp$p)

corResult <- arrange(corResult, p) %>% mutate(p.adj = p.adjust(p, method="BH"))

## -----------------------------------------------------------------------------
pCut = 0.05
corResult.sig <- filter(corResult, p <= pCut)
c6Plot <- plotSetHeatmap(geneTab=statTab[["mTORnone"]],
                         enrichTab=c6Res[["mTORnone"]],
                         topN=5, gmtFile=gmts[["C6"]],
                         #add asterix in front of the overlapped genes
                         asterixList = corResult.sig, 
                         anno=TRUE, i)

## ----fig_mTOR_C6_asterix, echo=FALSE, fig.path=plotDir, dev=c("png", "pdf"), fig.height=c6Plot[[1]][["height"]], fig.width=c6Plot[[1]][["width"]]----
#FIG# S13 B left
grid.draw(c6Plot[[1]]$plot)

## -----------------------------------------------------------------------------
hallmarkPlot <- plotSetHeatmap(geneTab=statTab[["mTORnone"]],
                               enrichTab=hallmarkRes[["mTORnone"]],
                               topN=5, gmtFile=gmts[["H"]],
                               asterixList = corResult.sig,
                               anno=TRUE, i)

## ----fig_mTOR_H_asterix, echo=FALSE, fig.path=plotDir, dev=c("png", "pdf"), fig.height=hallmarkPlot[[1]][["height"]], fig.width=hallmarkPlot[[1]][["width"]]----
#FIG# S13 B right
grid.draw(hallmarkPlot[[1]]$plot)

## -----------------------------------------------------------------------------
ddsCLL.BTK <- ddsCLL.norm[,ddsCLL.norm$group %in% "BTK"]
viabBTK <- Biobase::exprs(lpdCLL["D_002_4:5", ddsCLL.BTK$PatID])[1,]
stopifnot(all(ddsCLL.BTK$PatID == colnames(viabBTK)))  

## -----------------------------------------------------------------------------
#only keep genes that have counts higher than 10 in any sample
keep <- apply(assay(ddsCLL.BTK), 1, function(x) any(x >= 10)) 
ddsCLL.BTK <- ddsCLL.BTK[keep,]
dim(ddsCLL.BTK)

## -----------------------------------------------------------------------------
tmp = do.call(rbind, lapply(1:nrow(ddsCLL.BTK), function(i) {
  res = cor.test(viabBTK, assay(ddsCLL.BTK[i,])[1,], method = "pearson")
  data.frame(coef=unname(res$estimate), p=res$p.value)
}))

corResult <- tibble(ID = rownames(ddsCLL.BTK), 
                    symbol = rowData(ddsCLL.BTK)$symbol,
                    coef = tmp$coef,
                    p = tmp$p)

corResult <- arrange(corResult, p) %>% mutate(p.adj = p.adjust(p, method="BH"))

## -----------------------------------------------------------------------------
pCut = 0.05
corResult.sig <- filter(corResult, p <= pCut)
c6Plot <- plotSetHeatmap(geneTab=statTab[["BTKnone"]],
                         enrichTab=c6Res[["BTKnone"]],
                         topN=5, gmtFile=gmts[["C6"]],
                         #add asterix in front of the overlapped genes
                         asterixList = corResult.sig, 
                         anno=TRUE, i)

## ----fig_BTK_C6_asterix, echo=FALSE, fig.path=plotDir, dev=c("png", "pdf"), fig.height=c6Plot[[1]][["height"]], fig.width=c6Plot[[1]][["width"]]----
#FIG# S12 left
grid.draw(c6Plot[[1]]$plot)

## -----------------------------------------------------------------------------
hallmarkPlot <- plotSetHeatmap(geneTab=statTab[["BTKnone"]],
                               enrichTab=hallmarkRes[["BTKnone"]],
                               topN=5, gmtFile=gmts[["H"]],
                               asterixList = corResult.sig,
                               anno=TRUE, i)

## ----fig_BTK_H_asterix, echo=FALSE, fig.path=plotDir, dev=c("png", "pdf"), fig.height=hallmarkPlot[[1]][["height"]], fig.width=hallmarkPlot[[1]][["width"]]----
#FIG# S12 right
grid.draw(hallmarkPlot[[1]]$plot)

## -----------------------------------------------------------------------------
ddsCLL.MEK <- ddsCLL.norm[,ddsCLL.norm$group %in% "MEK"]
viabMEK <- Biobase::exprs(lpdCLL["D_012_4:5", ddsCLL.MEK$PatID])[1,]
stopifnot(all(ddsCLL.MEK$PatID == colnames(viabMEK)))  

## -----------------------------------------------------------------------------
#only keep genes that have counts higher than 10 in any sample
keep <- apply(assay(ddsCLL.MEK), 1, function(x) any(x >= 10)) 
ddsCLL.MEK <- ddsCLL.MEK[keep,]
dim(ddsCLL.MEK)

## -----------------------------------------------------------------------------
tmp = do.call(rbind, lapply(1:nrow(ddsCLL.MEK), function(i) {
  res = cor.test(viabMEK, assay(ddsCLL.MEK[i,])[1,], method = "pearson")
  data.frame(coef=unname(res$estimate), p=res$p.value)
}))

corResult <- tibble(ID = rownames(ddsCLL.MEK), 
                    symbol = rowData(ddsCLL.MEK)$symbol,
                    coef = tmp$coef,
                    p = tmp$p)

corResult <- arrange(corResult, p) %>% mutate(p.adj = p.adjust(p, method="BH"))

## -----------------------------------------------------------------------------
pCut = 0.05
corResult.sig <- filter(corResult, p <= pCut)
c6Plot <- plotSetHeatmap(geneTab=statTab[["MEKnone"]],
                         enrichTab=c6Res[["MEKnone"]],
                         topN=5, gmtFile=gmts[["C6"]],
                         anno=TRUE, i)

## ----fig_MEK_C6_asterix, echo=FALSE, fig.path=plotDir, dev=c("png", "pdf"), fig.height=c6Plot[[1]][["height"]], fig.width=c6Plot[[1]][["width"]]----
#FIG# S13 A left
grid.draw(c6Plot[[1]]$plot)

## -----------------------------------------------------------------------------
hallmarkPlot <- plotSetHeatmap(geneTab=statTab[["MEKnone"]],
                               enrichTab=hallmarkRes[["MEKnone"]],
                               topN=5, gmtFile=gmts[["H"]],
                               asterixList = corResult.sig,
                               anno=TRUE, i)

## ----fig_MEK_H_asterix, echo=FALSE, fig.path=plotDir, dev=c("png", "pdf"), fig.height=hallmarkPlot[[1]][["height"]], fig.width=hallmarkPlot[[1]][["width"]]----
#FIG# S13 A right
grid.draw(hallmarkPlot[[1]]$plot)

## ---- include=!exists(".standalone"), eval=!exists(".standalone")-------------
#  sessionInfo()

## ---- message=FALSE, warning=FALSE, include=FALSE-----------------------------
rm(list=ls())

## ---- message=FALSE, include=!exists(".standalone"), eval=!exists(".standalone")----
#  library("BloodCancerMultiOmics2017")
#  library("Biobase")
#  library("ggplot2")
#  library("grid")

## ----echo=FALSE---------------------------------------------------------------
plotDir = ifelse(exists(".standalone"), "", "part04/")
if(plotDir!="") if(!file.exists(plotDir)) dir.create(plotDir)

## -----------------------------------------------------------------------------
options(stringsAsFactors=FALSE)

## -----------------------------------------------------------------------------
data(list=c("drpar", "patmeta", "drugs", "mutCOM", "conctab"))

## -----------------------------------------------------------------------------
testFactors = function(msrmnt, factors, test="student", batch=NA) {
  
  # cut out the data
  tmp = colnames(factors)
  factors = data.frame(factors[rownames(msrmnt),], check.names=FALSE)
  colnames(factors) = tmp
  for(cidx in 1:ncol(factors))
    factors[,cidx] = factor(factors[,cidx], levels=c(0,1))
  
  # calculate the group size
  groupSizes = do.call(rbind, lapply(factors, function(tf) {
    tmp = table(tf)
    data.frame(n.0=tmp["0"], n.1=tmp["1"])
  }))
  
  # remove the factors with less then 2 cases per group
  factors = factors[,names(which(apply(groupSizes, 1,
                                       function(i) all(i>2)))), drop=FALSE]
  
  # calculate the effect
  effect = do.call(rbind, lapply(colnames(factors), function(tf) {
    tmp = aggregate(msrmnt ~ fac, data=data.frame(fac=factors[,tf]), mean)
    rownames(tmp) = paste("mean", tmp$fac, sep=".")
    tmp = t(tmp[2:ncol(tmp)])
    data.frame(TestFac=tf,
               DrugID=rownames(tmp),
               FacDr=paste(tf, rownames(tmp), sep="."),
               n.0=groupSizes[tf,"n.0"], n.1=groupSizes[tf,"n.1"],
               tmp, WM=tmp[,"mean.0"]-tmp[,"mean.1"])
  }))
  
  # do the test
  T = if(test=="student") {
    do.call(rbind, lapply(colnames(factors), function(tf) {
      tmp = do.call(rbind, lapply(colnames(msrmnt), function(dr) {
        res = t.test(msrmnt[,dr] ~ factors[,tf], var.equal=TRUE)
        data.frame(DrugID=dr, TestFac=tf,
                   pval=res$p.value, t=res$statistic,
                   conf1=res$conf.int[1], conf2=res$conf.int[2])
      }))
      tmp
    }))
  } else if(test=="anova") {
    do.call(rbind, lapply(colnames(factors), function(tf) {
      tmp = do.call(rbind, lapply(colnames(msrmnt), function(dr) {
        # make sure that the order in batch is the same as in msrmnt
        stopifnot(identical(rownames(msrmnt), names(batch)))
        res = anova(lm(msrmnt[,dr] ~ factors[,tf]+batch))
        data.frame(DrugID=dr, TestFac=tf, pval=res$`Pr(>F)`[1],
                   f=res$`F value`[1], meanSq1=res$`Mean Sq`[1],
                   meanSq2=res$`Mean Sq`[2])
      }))
      tmp
    }))
  } else {
    NA
  }
  
  enhanceObject = function(obj) {
    # give nice drug names
    obj$Drug = giveDrugLabel(obj$DrugID, conctab, drugs)
    # combine the testfac and drug id
    obj$FacDr = paste(obj$TestFac, obj$DrugID, sep=".")
    # select just the drug name
    obj$DrugID2 = substring(obj$DrugID, 1, 5)
    obj
  }
  
  list(effect=effect, test=enhanceObject(T))
}

## -----------------------------------------------------------------------------
## VIABILITIES
## list of matrices; one matrix per screen/day
## each matrix contains all CLL patients
measurements=list()

### Main Screen
patM = colnames(drpar)[which(patmeta[colnames(drpar),"Diagnosis"]=="CLL")]
measurements[["main"]] =
  do.call(cbind,
          lapply(list("viaraw.1","viaraw.2","viaraw.3","viaraw.4","viaraw.5"),
                 function(viac) {
                   tmp = t(assayData(drpar)[[viac]][,patM])
                   colnames(tmp) = paste(colnames(tmp), conctab[colnames(tmp),
                                                                paste0("c",substring(viac,8,8))], sep="-")
                   tmp
                 }))

pats = sort(unique(patM))

## TESTING FACTORS
testingFactors = list()
# ighv
ighv = setNames(patmeta[pats, "IGHV"], nm=pats)
# mutations
tmp = cbind(IGHV=ifelse(ighv=="U",1,0), assayData(mutCOM)$binary[pats,])
testingFactors[["mutation"]] = tmp[,-grep("Chromothripsis", colnames(tmp))]

# BATCHES
j = which(pData(drpar)[patM, "ExpDate"] < as.Date("2014-01-01"))
k = which(pData(drpar)[patM, "ExpDate"] < as.Date("2014-08-01") &
            pData(drpar)[patM, "ExpDate"] > as.Date("2014-01-01"))
l = which(pData(drpar)[patM, "ExpDate"] > as.Date("2014-08-01"))

measurements[["main"]] = measurements[["main"]][c(patM[j], patM[k], patM[l]),]
batchvec = factor(
  setNames(c(rep(1, length(j)), rep(2, length(k)), rep(3, length(l))),
           nm=c(patM[j], patM[k], patM[l])))

# LABELS FOR GROUPING
beelabs = t(sapply(colnames(testingFactors[["mutation"]]), function(fac) {
  if(fac=="IGHV")
    c(`0`="IGHV mut", `1`="IGHV unmut")
  else if(grepl("[[:upper:]]",fac)) # if all letters are uppercase
    c(`0`=paste(fac, "wt"),`1`=paste(fac, "mt"))
  else
    c(`0`="wt",`1`=fac)
}))

## -----------------------------------------------------------------------------
allresultsT = testFactors(msrmnt=measurements[["main"]],
                          factors=testingFactors[["mutation"]],
                          test="student", batch=NA)

resultsT = allresultsT$test
resultsT$adj.pval = p.adjust(resultsT$pval, method="BH")

## -----------------------------------------------------------------------------
allresultsA = testFactors(msrmnt=measurements[["main"]],
                          factors=testingFactors[["mutation"]],
                          test="anova", batch=batchvec)

resultsA = allresultsA$test
resultsA$adj.pval = p.adjust(resultsA$pval, method="BH")

## ----batchEffect, results='asis', echo=FALSE, fig.path=plotDir, dev=c('png','pdf'), fig.height=20, fig.width=14----
#FIG# S30
xylim = 1e-8
tmp = merge(resultsT[,c("FacDr","Drug","DrugID","DrugID2","FacDr","pval")],
            resultsA[,c("FacDr","pval")], by.x="FacDr", by.y="FacDr")
tmp$DrugName = toCaps(drugs[tmp$DrugID2, "name"])
tmp$Shape = ifelse(tmp$pval.x < xylim | tmp$pval.y < xylim, "tri", "dot")
tmp$pval.cens.x = ifelse(tmp$pval.x < xylim, xylim, tmp$pval.x)
tmp$pval.cens.y = ifelse(tmp$pval.y < xylim, xylim, tmp$pval.y)

ggplot(tmp) + geom_abline(intercept=0, slope=1, colour="hotpink") +
  geom_point(aes(x=pval.cens.x, y=pval.cens.y, shape=Shape), alpha=0.6) +
  facet_wrap(~DrugName, ncol=7) +
  scale_x_log10(labels=scientific_10, limits=c(xylim,1)) +
  scale_y_log10(labels=scientific_10, limits=c(xylim,1)) +
  theme_bw() + coord_fixed() + xlab("P-value from Student t-test") +
  ylab("P-value from ANOVA (including batch group factor)") +
  theme(axis.text.x=element_text(angle=0, hjust=1, vjust=0.5, size=rel(1.5)),
        axis.title=element_text(size=18)) + guides(shape=FALSE)

## -----------------------------------------------------------------------------
badrugs = c("D_008", "D_025")

measurements = lapply(measurements, function(drres) {
  drres[,-grep(paste(badrugs, collapse="|"), colnames(drres))]
})

## -----------------------------------------------------------------------------
allresults1 = lapply(measurements, function(measurement) {
  testFactors(msrmnt=measurement, factors=testingFactors[["mutation"]],
              test="student", batch=NA)
})

effects1 = lapply(allresults1, function(res) res[["effect"]])
results1 = lapply(allresults1, function(res) res[["test"]])

results1 = lapply(results1, function(res) {
  res$adj.pval = p.adjust(res$pval, method="BH")
  res
})

## -----------------------------------------------------------------------------
measurements1 = measurements
testingFactors1 = testingFactors
beelabs1 = beelabs

## ----echo=FALSE---------------------------------------------------------------
## CREATE THE PLOTS
plotmp01 = BloodCancerMultiOmics2017:::run.ggvolcGr2(results=results1, effects=effects1,
                                                     screen="main", mts="IGHV", fdr=0.1, maxX=0.75,
                                                     maxY=NA, expY=0.05, hghBox=0.15, axisMarkY=4,
                                                     breaksX=c(-0.75,-0.5,-0.25,0,0.25,0.5,0.75),
                                                     arrowLength=0.5, Xhang=0.3, minConc=3)

plotmp02 = BloodCancerMultiOmics2017:::run.ggvolcGr2(results=results1, effects=effects1,
                                                     screen="main", mts="trisomy12", fdr=0.1,
                                                     maxX=0.75, maxY=NA, expY=0.05, hghBox=0.15,
                                                     axisMarkY=4,
                                                     breaksX=c(-0.75,-0.5,-0.25,0,0.25,0.5,0.75),
                                                     arrowLength=0.5, Xhang=0.3, minConc=1)

## ----volc_IGHV, echo=FALSE, fig.path=plotDir, dev=c("png", "pdf"), fig.width=plotmp01$IGHV$figure$width, fig.height=plotmp01$IGHV$figure$height----
#FIG# 4B
grid.draw(plotmp01$IGHV$figure$plot)

## ----volc_IGHV_legend, echo=FALSE, fig.path=plotDir, dev=c("png", "pdf"), fig.width=plotmp01$IGHV$legend$width, fig.height=plotmp01$IGHV$legend$height, out.height=300, out.width=150----
#FIG# 4B legend
grid.draw(plotmp01$IGHV$legend$plot)

## ----volc_trisomy12, echo=FALSE, fig.path=plotDir, dev=c("png", "pdf"), fig.width=plotmp02$trisomy12$figure$width, fig.height=plotmp02$trisomy12$figure$height----
#FIG# 4C
grid.draw(plotmp02$trisomy12$figure$plot)

## ----volc_trisomy12_legend, echo=FALSE, fig.path=plotDir, dev=c("png", "pdf"), fig.width=plotmp02$trisomy12$legend$width, fig.height=plotmp02$trisomy12$legend$height, out.height=300, out.width=150----
#FIG# 4C legend
grid.draw(plotmp02$trisomy12$legend$plot)

## -----------------------------------------------------------------------------
fac2test = lapply(measurements, function(mea) {
  tf = testingFactors[["mutation"]][rownames(mea),]
  names(which(apply(tf,2,function(i) {
    if(length(table(i,tf[,1]))!=4)
      FALSE
    else
      all(table(i,tf[,1])>2)
  })))
})

## -----------------------------------------------------------------------------
measurements2 = setNames(lapply(names(measurements), function(mea) {
  ig = testingFactors[["mutation"]][rownames(measurements[[mea]]),"IGHV"]
  patU = names(which(ig==1))
  patM = names(which(ig==0))
  list(U=measurements[[mea]][patU,], M=measurements[[mea]][patM,])
}), nm=names(measurements))

## -----------------------------------------------------------------------------
allresults2 = setNames(lapply(names(measurements2), function(mea) {
  list(U = testFactors(msrmnt=measurements2[[mea]]$U,
                       factors=testingFactors[["mutation"]][
                         rownames(measurements2[[mea]]$U),fac2test[[mea]]]),
       M = testFactors(msrmnt=measurements2[[mea]]$M,
                       factors=testingFactors[["mutation"]][
                         rownames(measurements2[[mea]]$M),fac2test[[mea]]]))
}), nm=names(measurements2))

## -----------------------------------------------------------------------------
results2 = lapply(allresults2, function(allres) {
  list(U=allres[["U"]][["test"]], M=allres[["M"]][["test"]])
})

effects2 = lapply(allresults2, function(allres) {
  list(U=allres[["U"]][["effect"]], M=allres[["M"]][["effect"]])
})

## -----------------------------------------------------------------------------
results2 = lapply(results2, function(res) {
  tmp = p.adjust(c(res$U$pval,res$M$pval), method="BH")
  l = length(tmp)
  res$U$adj.pval = tmp[1:(l/2)]
  res$M$adj.pval = tmp[(l/2+1):l]
  res
})

## -----------------------------------------------------------------------------
testingFactors2 = testingFactors
beelabs2 = beelabs

## ----echo=FALSE---------------------------------------------------------------
## CREATE THE PLOTS
plotmp03 = BloodCancerMultiOmics2017:::run.ggvolcGr2(results=results2$main, effects=effects2$main,
                                                     screen="U", mts="trisomy12", fdr=0.1,
                                                     maxX=0.75, maxY=7, expY=0.05, hghBox=NA,
                                                     axisMarkY=4,
                                                     breaksX=c(-0.75,-0.5,-0.25,0,0.25,0.5,0.75),
                                                     arrowLength=0.5, Xhang=0.3, minConc=1,
                                                     fixedHght=6)

plotmp04 = BloodCancerMultiOmics2017:::run.ggvolcGr2(results=results2$main, effects=effects2$main,
                                                     screen="M", mts="trisomy12", fdr=0.1,
                                                     maxX=0.75, maxY=7, expY=0.05, hghBox=NA,
                                                     axisMarkY=4,
                                                     breaksX=c(-0.75,-0.5,-0.25,0,0.25,0.5,0.75),
                                                     arrowLength=0.5, Xhang=0.3, minConc=1,
                                                     fixedHght=6)

## ----volc_trisomy12_U, echo=FALSE, fig.path=plotDir, dev=c("png", "pdf"), fig.width=plotmp03$trisomy12$figure$width, fig.height=plotmp03$trisomy12$figure$height----
#FIG# S21 right
grid.draw(plotmp03$trisomy12$figure$plot)

## ----volc_trisomy12_U_legend, echo=FALSE, fig.path=plotDir, dev=c("png", "pdf"), fig.width=plotmp03$trisomy12$legend$width, fig.height=plotmp03$trisomy12$legend$height, out.height=300, out.width=150----
#FIG# S21 right legend
grid.draw(plotmp03$trisomy12$legend$plot)

## ----volc_trisomy12_M, echo=FALSE, fig.path=plotDir, dev=c("png", "pdf"), fig.width=plotmp04$trisomy12$figure$width, fig.height=plotmp04$trisomy12$figure$height----
#FIG# S21 left
grid.draw(plotmp04$trisomy12$figure$plot)

## ----volc_trisomy12_M_legend, echo=FALSE, fig.path=plotDir, dev=c("png", "pdf"), fig.width=plotmp04$trisomy12$legend$width, fig.height=plotmp04$trisomy12$legend$height, out.height=300, out.width=150----
#FIG# S21 left legend
grid.draw(plotmp04$trisomy12$legend$plot)

## -----------------------------------------------------------------------------
## VIABILITIES
### main
pats = colnames(drpar)
# make the big matrix with viabilities
measureTMP = do.call(cbind,
                     lapply(list("viaraw.1","viaraw.2","viaraw.3",
                                 "viaraw.4","viaraw.5"), function(viac) {
                                   tmp = t(assayData(drpar)[[viac]][,pats])
                                   colnames(tmp) = paste(colnames(tmp),
                                                         conctab[colnames(tmp),
                                                                 paste0("c",substring(viac,8,8))], sep="-")
                                   tmp
                                 }))

# select diagnosis to work on
pats4diag = tapply(pats, patmeta[pats,"Diagnosis"], function(i) i)

diags = names(which(table(patmeta[pats,"Diagnosis"])>2))
diags = diags[-which(diags=="CLL")]
# there will be two lists: one with CLL and the second with other diagnosis
# (first one is passed as argument to the createObjects function)
pats4diag2 = pats4diag[diags]

# function that creates testingFactors, measurements and beelabs
createObjects = function(pats4diag1, beefix="") {
  
  measurements=list()
  testingFactors=list()
  # make the list for testing
  for(m in names(pats4diag1)) {
    for(n in names(pats4diag2)) {
      p1 = pats4diag1[[m]]
      p2 = pats4diag2[[n]]
      measurements[[paste(m,n,sep=".")]] = measureTMP[c(p1, p2),]
      testingFactors[[paste(m,n,sep=".")]] = setNames(c(rep(0,length(p1)),
                                                        rep(1,length(p2))),
                                                      nm=c(p1,p2))
    }
  }
  
  # reformat testingFactors to the df
  pats=sort(unique(c(unlist(pats4diag1),unlist(pats4diag2))))
  testingFactors = as.data.frame(
    do.call(cbind, lapply(testingFactors, function(tf) {
      setNames(tf[pats], nm=pats)
    })))
  
  # Labels for beeswarms
  beelabs = t(sapply(colnames(testingFactors), function(fac) {
    tmp = unlist(strsplit(fac, ".", fixed=TRUE))
    c(`0`=paste0(tmp[1], beefix),`1`=tmp[2])
  }))
  
  return(list(msrmts=measurements, tf=testingFactors, bl=beelabs))
}

# all CLL together
res = createObjects(pats4diag1=pats4diag["CLL"])
measurements3 = res$msrmts
testingFactors3 = res$tf
beelabs3 = res$bl

## ---- results='hide'----------------------------------------------------------
allresults3 = setNames(lapply(names(measurements3), function(mea) {
  tmp = data.frame(testingFactors3[,mea])
  colnames(tmp) = mea
  rownames(tmp) = rownames(testingFactors3)
  testFactors(msrmnt=measurements3[[mea]], factors=tmp)
}), nm=names(measurements3))

## -----------------------------------------------------------------------------
results3 = lapply(allresults3, function(res) res[["test"]])
effects3 = lapply(allresults3, function(res) res[["effect"]])

## -----------------------------------------------------------------------------
results3 = lapply(results3, function(res) {
  res$adj.pval = p.adjust(res$pval, method="BH")
  res
})

## ----echo=FALSE---------------------------------------------------------------
tmpheat = BloodCancerMultiOmics2017:::ggheat(results3, effects3)

## ----cll.diag, echo=FALSE, fig.path=plotDir, dev=c("png", "pdf"), fig.width=tmpheat$figure$width, fig.height=tmpheat$figure$height----
#FIG# S7 plot
grid.draw(tmpheat$figure$plot)

## ----cll.diag.legend, echo=FALSE, fig.path=plotDir, dev=c("png", "pdf"), fig.width=tmpheat$legend$width, fig.height=tmpheat$legend$height----
#FIG# S7 legend
grid.draw(tmpheat$legend$plot)

## -----------------------------------------------------------------------------
data(drugs, lpdAll, mutCOM, conctab)

## -----------------------------------------------------------------------------
lpdCLL = lpdAll[ , lpdAll$Diagnosis %in% "CLL"]

## ----beesMutMain, fig.path=plotDir, dev=c("png", "pdf"), fig.width=8, fig.height=10, out.width=280, out.height=350----
#FIG# 4D
BloodCancerMultiOmics2017:::beeF(diag="CLL", "D_010_2", "TP53", cs=T, y1=0, y2=1.2, custc=T)
BloodCancerMultiOmics2017:::beeF(diag="CLL", "D_006_3", "TP53", cs=T,y1=0, y2=1.2, custc=T)
BloodCancerMultiOmics2017:::beeF(diag="CLL", "D_063_5", "CREBBP", cs=T, y1=0, y2=1.2, custc=T)
BloodCancerMultiOmics2017:::beeF(diag="CLL", "D_056_5", "PRPF8", cs=T, y1=0, y2=1.2, custc=T)
BloodCancerMultiOmics2017:::beeF(diag="CLL", "D_012_5", "trisomy12", cs=F,y1=0.6, y2=1.2, custc=T)

## ----beesMutSupp, fig.path=plotDir, fig.width = 18, fig.height = 20, dev = c("png", "pdf")----
#FIG# S17
par(mfrow = c(3,4), mar=c(5,4.5,5,2))

BloodCancerMultiOmics2017:::beeF(diag="CLL", drug="D_159_3", mut="TP53", cs=T, y1=0, y2=1.2, custc=T)
BloodCancerMultiOmics2017:::beeF(diag="CLL", "D_006_2", "del17p13", cs=T, y1=0, y2=1.2, custc=T)
BloodCancerMultiOmics2017:::beeF(diag="CLL", "D_159_3", "del17p13", cs=T, y1=0, y2=1.2, custc=T)
BloodCancerMultiOmics2017:::beeF(diag="CLL", "D_010_2", "del17p13", cs=T, y1=0, y2=1.2, custc=T)
BloodCancerMultiOmics2017:::beeF(diag="MCL", "D_006_2", "TP53", cs=T,  y1=0, y2=1.2, custc=T)
BloodCancerMultiOmics2017:::beeF(diag="MCL", "D_010_2", "TP53", cs=T,  y1=0, y2=1.2, custc=T)
BloodCancerMultiOmics2017:::beeF(diag=c("HCL", "HCL-V"), "D_012_3", "BRAF", cs=T, y1=0, y2=1.2,
                                 custc=F)
BloodCancerMultiOmics2017:::beeF(diag=c("HCL", "HCL-V"), "D_083_4", "BRAF", cs=T, y1=0, y2=1.2,
                                 custc=F)
BloodCancerMultiOmics2017:::beeF(diag="CLL", "D_012_5", "KRAS", cs=T, y1=0.6, y2=1.2, custc=T)
BloodCancerMultiOmics2017:::beeF(diag="CLL", "D_083_5", "KRAS", cs=T, y1=0.6, y2=1.45, custc=T)
BloodCancerMultiOmics2017:::beeF(diag="CLL", "D_081_4", "UMODL1", cs=T,  y1=0, y2=1.2, custc=T)
BloodCancerMultiOmics2017:::beeF(diag="CLL", "D_001_4", "UMODL1", cs=T,  y1=0, y2=1.2, custc=T)

## ----colorbar, echo=FALSE, fig.path=plotDir, dev=c("png", "pdf"), fig.width=2, fig.height=4, out.height=200, out.width=100----

# the below function comes from the CRAN package named monogeneaGM
colorBar <-function(colpalette=NULL, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), 
                    tit="") {
  
  if(is.null(colpalette)) {
    scale <- (length(colpalette)-1)/(max-min)
    rwb <- c("#99000D","#FB6A4A","white","#6BAED6","#084594")
    colpalette<- colorRampPalette(rwb, space="Lab")(101)
  }
  
  scale <- (length(colpalette)-1)/(max-min)
  plot(c(0,10), c(min,max), type="n", bty="n", xaxt="n", xlab="", yaxt="n", ylab="", main=tit)
  axis(2, ticks, las=1)
  for (i in 1:(length(colpalette)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=colpalette[i], border=NA)
  }
}

colorBar(colorRampPalette(c('coral1','blue4'))(100), min=0, max = 1,
         ticks=c(0,0.5,1))

## ----bee-pretreatment, fig.path=plotDir, fig.width=15, fig.height=10, dev = c("png", "pdf")----
#FIG# S18
par(mfrow = c(2,3), mar=c(5,4.5,2,2)) 

BloodCancerMultiOmics2017:::beePretreatment(lpdCLL, "D_006_1:5", y1=0.2, y2=1.3, fac="TP53",
                                            val=c(0,1), name="Fludarabine")
BloodCancerMultiOmics2017:::beePretreatment(lpdCLL, "D_006_1:5", y1=0.2, y2=1.3, fac="TP53",
                                            val=c(0),   name="p53 wt:  Fludarabine")
BloodCancerMultiOmics2017:::beePretreatment(lpdCLL, "D_006_1:5", y1=0.2, y2=1.3, fac="TP53",
                                            val=c(1),   name="p53 mut: Fludarabine")

BloodCancerMultiOmics2017:::beePretreatment(lpdCLL, "D_002_4:5", y1=0.4, y2=1.3,
                                            fac="IGHV Uppsala U/M", val=c(0,1), name="Ibrutinib")
BloodCancerMultiOmics2017:::beePretreatment(lpdCLL, "D_002_4:5", y1=0.4, y2=1.3,
                                            fac="IGHV Uppsala U/M", val=c(0), name="U-CLL: Ibrutinib")
BloodCancerMultiOmics2017:::beePretreatment(lpdCLL, "D_002_4:5", y1=0.4, y2=1.3,
                                            fac="IGHV Uppsala U/M", val=c(1), name="M-CLL: Ibrutinib")

## ---- include=!exists(".standalone"), eval=!exists(".standalone")-------------
#  sessionInfo()

## ---- message=FALSE, warning=FALSE, include=FALSE-----------------------------
rm(list=ls())

## ---- message=FALSE, warning=FALSE, include=!exists(".standalone"), eval=!exists(".standalone")----
#  library("BloodCancerMultiOmics2017")
#  library("Biobase")
#  library("dplyr")
#  library("tidyr")
#  library("broom")
#  library("ggplot2")
#  library("grid")
#  library("gridExtra")
#  library("reshape2")
#  library("foreach")
#  library("doParallel")
#  library("scales")
#  library("knitr")
#  registerDoParallel()

## ----echo=FALSE---------------------------------------------------------------
plotDir = ifelse(exists(".standalone"), "", "part05/")
if(plotDir!="") if(!file.exists(plotDir)) dir.create(plotDir)

## -----------------------------------------------------------------------------
options(stringsAsFactors=FALSE)

## -----------------------------------------------------------------------------
# get drug responsee data
get.drugresp <- function(lpd) {
  drugresp = t(Biobase::exprs(lpd[fData(lpd)$type == 'viab'])) %>%
    dplyr::tbl_df() %>% dplyr::select(-ends_with(":5")) %>%
    dplyr::mutate(ID = colnames(lpd)) %>%
    tidyr::gather(drugconc, viab, -ID) %>%
    dplyr::mutate(drug = drugs[substring(drugconc, 1, 5), "name"],
                  conc = sub("^D_([0-9]+_)", "", drugconc)) %>%
    dplyr::mutate(conc = as.integer(gsub("D_CHK_", "", conc)))
  
  drugresp
}

# extract mutations and IGHV status
get.somatic <- function(lpd) {
  somatic = t(Biobase::exprs(lpd[Biobase::fData(lpd)$type == 'gen' | 
                                   Biobase::fData(lpd)$type == 'IGHV']))
  ## rename IGHV Uppsala to 'IGHV' (simply)
  colnames(somatic)[grep("IGHV", colnames(somatic))] = "IGHV"
  
  ## at least 3 patients should have this mutation
  min.samples = which(Matrix::colSums(somatic, na.rm = T) > 2)
  somatic = dplyr::tbl_df(somatic[, min.samples]) %>%
    dplyr::select(-one_of("del13q14_bi", "del13q14_mono", 
                          "Chromothripsis", "RP11-766F14.2")) %>%
    dplyr::rename(del13q14 = del13q14_any) %>% 
    dplyr::mutate(ID = colnames(lpd)) %>%
    tidyr::gather(mutation, mut.value, -ID)
  somatic
}

## ----ggTheme------------------------------------------------------------------
t1<-theme(                              
  plot.background = element_blank(), 
  panel.grid.major = element_line(),
  panel.grid.major.x = element_line(linetype = "dotted", colour = "grey"),
  panel.grid.minor = element_blank(), 
  panel.border = element_blank(), 
  panel.background = element_blank(),
  axis.line = element_line(size=.4),
  axis.line.x = element_line(),
  axis.line.y = element_line(),
  axis.text.x  = element_text(angle=90, size=12, 
                              face="bold", hjust = 1, vjust = 0.4),
  axis.text.y = element_text(size = 14),
  axis.ticks.x = element_line(linetype = "dotted"),
  axis.ticks.length = unit(0.3,"cm"),
  axis.title.x = element_text(face="bold", size=16), 
  axis.title.y = element_text(face="bold", size=20),
  plot.title = element_text(face="bold", size=16, hjust = 0.5)
)

## theme for the legend
t.leg <-  theme(legend.title = element_text(face='bold', 
                                            hjust = 1, size=11),
                legend.position = c(0, 0.76),
                legend.key = element_blank(),
                legend.text = element_text(size=12),
                legend.background = element_rect(color = "black"))

## ----colorPalette-------------------------------------------------------------
colors= c("#015872","#3A9C94","#99977D","#ffbf00","#5991C7","#99cc00",
          "#D5A370","#801416","#B2221C","#ff5050","#33bbff","#5c5cd6",
          "#E394BB","#0066ff","#C0C0C0")

## -----------------------------------------------------------------------------
get.pretreat <- function(patmeta, lpd) {
  patmeta = patmeta[rownames(patmeta) %in% colnames(lpd),]
  data.frame(ID=rownames(patmeta), pretreat=!patmeta$IC50beforeTreatment) %>% 
    mutate(pretreat = as.factor(pretreat))
  
}

## -----------------------------------------------------------------------------
make.dr <- function(resp, features, patmeta, lpd) {
  treat = get.pretreat(patmeta, lpd)
  dr = full_join(resp, features) %>% 
    inner_join(treat) 
}

## -----------------------------------------------------------------------------
get.medp <- function(drugresp) {
  tab = drugresp %>% group_by(drug, conc) %>% 
    do(data.frame(v = .$viab, ID = .$ID)) %>% spread(ID, v)
  
  med.p = foreach(n=unique(tab$drug), .combine = cbind) %dopar% {
    tb = filter(tab, drug == n) %>% ungroup() %>% dplyr::select(-(drug:conc)) %>% 
      as.matrix %>% `rownames<-`(1:5)
    mdp = stats::medpolish(tb)
    df = as.data.frame(mdp$col) + mdp$overall
    colnames(df) <- n
    df
  }
  
  medp.viab = dplyr::tbl_df(med.p) %>% dplyr::mutate(ID = rownames(med.p)) %>%
    tidyr::gather(drug, viab, -ID) 
  medp.viab
}

## -----------------------------------------------------------------------------
get.labels <- function(pvals) {
  lev = levels(factor(pvals$mutation))
  lev = gsub("^(gain)([0-9]+)([a-z][0-9]+)$", "\\1(\\2)(\\3)", lev)
  lev =  gsub("^(del)([0-9]+)([a-z].+)$", "\\1(\\2)(\\3)", lev)
  lev = gsub("trisomy12", "trisomy 12", lev)
  lev
}

## -----------------------------------------------------------------------------
get.mutation.order <- function(lev) {
  ord = c("trisomy 12", "TP53",
          "del(11)(q22.3)", "del(13)(q14)",
          "del(17)(p13)",
          "gain(8)(q24)",
          "BRAF", "CREBBP", "PRPF8",
          "KLHL6", "NRAS", "ABI3BP", "UMODL1")
  mut.order = c(match(ord, lev),
                grep("Other", lev), grep("Below", lev))
  
  mut.order
}

## -----------------------------------------------------------------------------
get.drug.order <- function(pvals, drugs) {
  ## determine drug order by column sums of log-p values
  dr.order = pvals %>% 
    mutate(logp = -log10(p.value)) %>% 
    group_by(drug) %>% summarise(logsum = sum(logp)) 
  
  dr.order = inner_join(dr.order, pvals %>%
                          group_by(drug) %>% 
                          summarise(n = length(unique(mutation)))) %>% 
    arrange(desc(n), desc(logsum))
  
  dr.order = inner_join(dr.order, drugs %>% rename(drug = name))
  
  dr.order = left_join(dr.order, dr.order %>% 
                         group_by(`target_category`) ) %>%
    arrange(`target_category`, drug) %>%
    filter(! `target_category` %in% c("ALK", "Angiogenesis", "Other")) %>%
    filter(!is.na(`target_category`))
  
  dr.order
}

## -----------------------------------------------------------------------------
make.annot <- function(g, dr.order) {
  # make a color palette for drug pathways
  drug.class = c("#273649", "#647184", "#B1B2C8",
                 "#A7755D", "#5D2E1C", "#38201C")
  pathways = c("BH3","B-cell receptor","DNA damage",
               "MAPK", "PI3K", "Reactive oxygen species")
  names(pathways) = c("BH3", "BCR inhibitors", "DNA damage",
                      "MAPK", "PI3K", "ROS")
  
  for (i in 1:6) {
    prange = grep(pathways[i], dr.order$`target_category`)
    path.grob <- grobTree(rectGrob(gp=gpar(fill=drug.class[i])),
                          textGrob(names(pathways)[i], 
                                   gp = gpar(cex =0.8, col = "white")))
    g = g + 
      annotation_custom(path.grob, 
                        xmin = min(prange) -0.25 - 0.1 * ifelse(i == 2, 1, 0), 
                        xmax = max(prange) + 0.25 + 0.1 * ifelse(i == 2, 1, 0), 
                        ymin = -0.52, ymax = -0.2)
  }
  g
}

## -----------------------------------------------------------------------------
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
} ## end define

## -----------------------------------------------------------------------------
data(list=c("conctab", "drugs", "lpdAll", "patmeta"))

## ----preprocesslpd------------------------------------------------------------
lpdCLL <- lpdAll[ , lpdAll$Diagnosis=="CLL"]
## extract viability data for 5 different concentrations
drugresp = get.drugresp(lpdCLL)

## ----preprocessMuts-----------------------------------------------------------
## extract somatic variants
somatic = get.somatic(lpdCLL) %>% 
  mutate(mut.value = as.factor(mut.value))

## ----medPolish, warning=FALSE, results='hide'---------------------------------
## compute median polish patient effects and recalculate p-values
medp.viab = get.medp(drugresp)
dr = make.dr(medp.viab, somatic, patmeta, lpdCLL)

## ----pvalsAndFDR--------------------------------------------------------------
pvals = dr %>% group_by(drug, mutation) %>%
  do(tidy(t.test(viab ~ mut.value, data = ., var.equal = T))) %>%
  dplyr::select(drug, mutation, p.value)

# compute the FDR threshold
fd.thresh = 10
padj = p.adjust(pvals$p.value, method = "BH")
fdr = max(pvals$p.value[which(padj <= fd.thresh/100)])

## ----subsetPvals--------------------------------------------------------------
# selected mutations
select.mutations = c("trisomy12", "TP53",
                     "del11q22.3", "del13q14",
                     "del17p13",
                     "gain8q24",
                     "BRAF", "CREBBP", "PRPF8",
                     "KLHL6", "NRAS", "ABI3BP", "UMODL1")

pvals = filter(pvals, mutation != 'IGHV')
pvals = pvals %>% ungroup() %>%
  mutate(mutation = ifelse(p.value > fdr, 
                           paste0("Below ", fd.thresh,"% FDR"), mutation)) %>%
  mutate(mutation = ifelse(!(mutation %in% select.mutations) & 
                             !(mutation == paste0("Below ", fd.thresh,"% FDR")), 
                           "Other", mutation)) %>%
  filter(drug != "bortezomib" & drug != "NSC 74859")

## ----renameMuts---------------------------------------------------------------
## order of mutations
lev = get.labels(pvals)
folge = get.mutation.order(lev)

## ----pvalsForGGplot-----------------------------------------------------------
drugs = drugs[,c("name", "target_category")]
# get the drug order
dr.order = get.drug.order(pvals, drugs)

## -----------------------------------------------------------------------------
plot.pvalues <- function(pvals, dr.order, folge, colors, shapes) {
  g = ggplot(data = filter(pvals, drug %in% dr.order$drug)) +
    geom_point(aes(x = factor(drug, levels = dr.order$drug), y = -log10(p.value), 
                   colour = factor(mutation, levels(factor(mutation))[folge]),
                   shape  = factor(mutation, levels(factor(mutation))[folge])),
               size=5, show.legend = T)  + 
    scale_color_manual(name = "Mutations",
                       values = colors,
                       labels = lev[folge]) + 
    scale_shape_manual(name = "Mutations",
                       values = shapes,
                       labels = lev[folge]) + t1 + 
    labs(x = "", y = expression(paste(-log[10], "p")), title = "") +
    scale_y_continuous(expression(italic(p)*"-value"),
                       breaks=seq(0,10,5),
                       labels=math_format(expr=10^.x)(-seq(0,10,5))) 
  g
}

## ----pvalsMain, fig.path=plotDir, dev=c("png", "pdf"), fig.width=14, fig.height=10----
#FIG# 4A
## plot the p-values 
g = plot.pvalues(pvals, dr.order, folge, 
                 colors, shapes = c(rep(16,13), c(1,1)))

## add FDR threshold
g = g + geom_hline(yintercept = -log10(fdr),
                   linetype="dashed", size=0.3)

g = g + 
  annotation_custom(grob = textGrob(label = paste0("FDR", fd.thresh, "%"), 
                                    hjust = 1, vjust = 1, 
                                    gp = gpar(cex = 0.5,
                                              fontface = "bold",
                                              fontsize = 25)),
                    ymin = -log10(fdr) - 0.2, 
                    ymax = -log10(fdr) + 0.5, 
                    xmin = -1.3, xmax = 1.5) + 
  theme(legend.position = "none")

# generate pathway/target annotations for certain drug classes
#g = make.annot(g, dr.order)

# legend guide
leg.guides <- guides(colour = guide_legend(ncol = 1, 
                                           byrow = TRUE,
                                           override.aes = list(size = 3),
                                           title = "Mutations",
                                           label.hjust = 0,
                                           keywidth = 0.4,
                                           keyheight = 0.8), 
                     shape = guide_legend(ncol = 1, 
                                          byrow = TRUE,
                                          title = "Mutations",
                                          label.hjust = 0,
                                          keywidth = 0.4,
                                          keyheight = 0.8))

# create a legend grob
legend = g_legend(g + t.leg +  leg.guides)

## arranget the main plot and the legend
# using grid graphics
gt <- ggplot_gtable(ggplot_build(g + theme(legend.position = 'none')))
gt$layout$clip[gt$layout$name == "panel"] <- "off"

grid.arrange(gt, legend,
             ncol=2, nrow=1, widths=c(0.92,0.08))

## -----------------------------------------------------------------------------
## lm(viab ~ mutation + pretreatment.status)
pvals = dr %>% group_by(drug, mutation) %>%
  do(tidy(lm(viab ~ mut.value + pretreat, data = .))) %>%
  filter(term == 'mut.value1') %>%
  dplyr::select(drug, mutation, p.value)

# compute the FDR threshold
fd.thresh = 10
padj = p.adjust(pvals$p.value, method = "BH")
fdr = max(pvals$p.value[which(padj <= fd.thresh/100)])


pvals = filter(pvals, mutation != 'IGHV')
pvals = pvals %>% ungroup() %>%
  mutate(mutation = ifelse(p.value > fdr,
                           paste0("Below ", fd.thresh,"% FDR"),
                           mutation)) %>%
  mutate(mutation = ifelse(!(mutation %in% select.mutations) &
                             !(mutation == paste0("Below ",
                                                  fd.thresh,"% FDR")), 
                           "Other", mutation)) %>%
  filter(drug != "bortezomib" & drug != "NSC 74859")


lev = get.labels(pvals)
folge = get.mutation.order(lev)

# get the drug order
dr.order = get.drug.order(pvals, drugs)

mut.order = folge[!is.na(folge)]

## ----pvalsSupp, fig.path=plotDir, dev=c("png", "pdf"), fig.width=14, fig.height=10----
#FIG# S19
## plot the p-values 
g = plot.pvalues(pvals, dr.order, mut.order, 
                 colors[which(!is.na(folge))], shapes = c(rep(16,9), c(1,1)))

## add FDR threshold
g = g + geom_hline(yintercept = -log10(fdr),
                   linetype="dashed", size=0.3)

g = g + 
  annotation_custom(grob = textGrob(label = paste0("FDR", fd.thresh, "%"), 
                                    hjust = 1, vjust = 1, 
                                    gp = gpar(cex = 0.5,
                                              fontface = "bold",
                                              fontsize = 25)),
                    ymin = -log10(fdr) - 0.2, 
                    ymax = -log10(fdr) + 0.5, 
                    xmin = -1.3, xmax = 1.5) + 
  theme(legend.position = "none")

# generate pathway/target annotations for certain drug classes
#g = make.annot(g, dr.order)

# legend guide
leg.guides <- guides(colour = guide_legend(ncol = 1, 
                                           byrow = TRUE,
                                           override.aes = list(size = 3),
                                           title = "Mutations",
                                           label.hjust = 0,
                                           keywidth = 0.4,
                                           keyheight = 0.8), 
                     shape = guide_legend(ncol = 1, 
                                          byrow = TRUE,
                                          title = "Mutations",
                                          label.hjust = 0,
                                          keywidth = 0.4,
                                          keyheight = 0.8))

# create a legend grob
legend = g_legend(g + t.leg +  leg.guides)

## arranget the main plot and the legend
# using grid graphics
gt <- ggplot_gtable(ggplot_build(g + theme(legend.position = 'none')))
gt$layout$clip[gt$layout$name == "panel"] <- "off"

grid.arrange(gt, legend,
             ncol=2, nrow=1, widths=c(0.92,0.08))

## -----------------------------------------------------------------------------
pvals.main = dr %>% group_by(drug, mutation) %>%
  do(tidy(t.test(viab ~ mut.value, data = ., var.equal = T))) %>%
  dplyr::select(drug, mutation, p.value)

p.main.adj = p.adjust(pvals.main$p.value, method = "BH")
fdr.main = max(pvals.main$p.value[which(p.main.adj <= fd.thresh/100)])

pvals.main = filter(pvals.main, mutation != "IGHV") %>%
  rename(p.main = p.value)


## lm(viab ~ mutation + pretreatment.status)
pvals.sup = dr %>% group_by(drug, mutation) %>%
  do(tidy(lm(viab ~ mut.value + pretreat, data = .))) %>%
  filter(term == 'mut.value1') %>%
  dplyr::select(drug, mutation, p.value)

p.sup.adj = p.adjust(pvals.sup$p.value, method = "BH")
fdr.sup = max(pvals.sup$p.value[which(p.sup.adj <= fd.thresh/100)])

pvals.sup = filter(pvals.sup, mutation != "IGHV") %>%
  rename(p.sup = p.value)


pvals = inner_join(pvals.main, pvals.sup)
pvals = mutate(pvals, signif = ifelse(p.main > fdr.main, 
                                      ifelse(p.sup > fdr.sup, 
                                             "Below 10% FDR in both models", 
                                             "Significant with pretreatment accounted"), 
                                      ifelse(p.sup > fdr.sup, 
                                             "Significant without pretreatment in the model", 
                                             "Significant in both models")))

t2<-theme(                              
  plot.background = element_blank(), 
  panel.grid.major = element_line(),
  panel.grid.major.x = element_line(),
  panel.grid.minor = element_blank(), 
  panel.border = element_blank(), 
  panel.background = element_blank(),
  axis.line = element_line(size=.4),
  axis.line.x = element_line(),
  axis.line.y = element_line(),
  axis.text.x  = element_text(size=12),
  axis.text.y = element_text(size = 12),
  axis.title.x = element_text(face="bold", size=12), 
  axis.title.y = element_text(face="bold", size=12),
  legend.title = element_text(face='bold', 
                              hjust = 1, size=10),
  legend.position = c(0.78, 0.11),
  legend.key = element_blank(),
  legend.text = element_text(size=10),
  legend.background = element_rect(color = "black")
)

## ----pvalComparisonScatterplot, fig.path=plotDir, dev=c("png", "pdf"), fig.width=10, fig.height=7----
#FIG# S19
ggplot(pvals, aes(-log10(p.main), -log10(p.sup), colour = factor(signif))) + 
  geom_point() + t2 + labs(x = expression(paste(-log[10], "p, pretreatment not considered", sep = "")),
                           y = expression(paste(-log[10], "p, accounting for pretreatment", sep = ""))) +
  coord_fixed() +
  scale_x_continuous(breaks = seq(0,9,by = 3)) + 
  scale_y_continuous(breaks = seq(0,9,by = 3)) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
  scale_color_manual(name = "Statistical Significance",
                     values = c("#F1BB7B","#669999", "#FD6467", "#5B1A18"))

## -----------------------------------------------------------------------------
signif.in.one = filter(pvals, 
                       signif %in% c("Significant with pretreatment accounted",
                                     "Significant without pretreatment in the model")) %>%
  arrange(signif)
kable(signif.in.one, digits = 4,
      align = c("l", "l", "c", "c", "c"),
      col.names = c("Drug", "Mutation", "P-value (Main)",
                    "P-value (Supplement)", "Statistical significance"),
      format.args = list(width = 14))

## ---- comment=NA, eval=FALSE--------------------------------------------------
#  print(kable(signif.in.one, format = "latex", digits = 4,
#         align = c("l", "l", "c", "c", "c"),
#        col.names = c("Drug", "Mutation", "P-value (Main)",
#                      "P-value (Supplement)", "Statistical significance")))

## ---- include=!exists(".standalone"), eval=!exists(".standalone")-------------
#  sessionInfo()

## ---- message=FALSE, warning=FALSE, include=FALSE-----------------------------
rm(list=ls())

## ---- message=FALSE, include=!exists(".standalone"), eval=!exists(".standalone")----
#  library("BloodCancerMultiOmics2017")
#  library("Biobase")
#  library("ggbeeswarm")
#  library("ggplot2")
#  library("gridExtra")
#  library("dplyr")

## ----echo=FALSE---------------------------------------------------------------
plotDir = ifelse(exists(".standalone"), "", "part15/")
if(plotDir!="") if(!file.exists(plotDir)) dir.create(plotDir)

## -----------------------------------------------------------------------------
data(list= c("validateExp","lpdAll"))

## -----------------------------------------------------------------------------
plotTab <- filter(validateExp, Drug %in% c("Ganetespib", "Onalespib")) %>%
  mutate(IGHV = Biobase::exprs(lpdAll)["IGHV Uppsala U/M", patientID]) %>%
  filter(!is.na(IGHV)) %>%
  mutate(IGHV = as.factor(ifelse(IGHV == 1, "M","U")),
         Concentration = as.factor(Concentration))

## -----------------------------------------------------------------------------
pTab <- group_by(plotTab, Drug, Concentration) %>%
  do(data.frame(p = t.test(viab ~ IGHV, .)$p.value)) %>%
  mutate(p = format(p, digits =2, scientific = TRUE))

## ----HSP90confirm, fig.width=14, fig.height=5, warning=FALSE, fig.path=plotDir, dev=c("png", "pdf")----
pList <- group_by(plotTab, Drug) %>% 
  do(plots = ggplot(., aes(x=Concentration, y = viab)) + 
       stat_boxplot(geom = "errorbar", width = 0.3,
                    position = position_dodge(width=0.6), 
                    aes(dodge = IGHV)) +
       geom_boxplot(outlier.shape = NA, position = position_dodge(width=0.6), 
                    col="black", width=0.5, aes(dodge = IGHV)) + 
       geom_beeswarm(size=1,dodge.width=0.6, aes(col=IGHV)) +
       theme_classic() +
       scale_y_continuous(expand = c(0, 0),breaks=seq(0,1.2,0.20)) +
       coord_cartesian(ylim = c(0,1.30)) +
       xlab("Concentration (µM)") + ylab("Viability") + 
       ggtitle(unique(.$Drug)) +
       geom_text(data=filter(pTab, Drug == unique(.$Drug)), y = 1.25, 
                 aes(x=Concentration, label=sprintf("p=%s",p)),
                 size = 4.5) + 
       theme(axis.line.x = element_blank(),
             axis.ticks.x = element_blank(),
             axis.text  = element_text(size=15),
             axis.title = element_text(size =15),
             legend.text = element_text(size=13),
             legend.title = element_text(size=15),
             plot.title = element_text(face="bold", hjust=0.5, size=17),
             plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))) 
grid.arrange(grobs = pList$plots, ncol =2)

## -----------------------------------------------------------------------------
plotTab <- filter(validateExp, Drug %in%
                    c("Cobimetinib","SCH772984","Trametinib")) %>%
  mutate(Trisomy12 = Biobase::exprs(lpdAll)["trisomy12", patientID]) %>%
  filter(!is.na(Trisomy12)) %>%
  mutate(Trisomy12 = as.factor(ifelse(Trisomy12 == 1, "present","absent")),
         Concentration = as.factor(Concentration))

## -----------------------------------------------------------------------------
pTab <- group_by(plotTab, Drug, Concentration) %>% 
  do(data.frame(p = t.test(viab ~ Trisomy12, .)$p.value)) %>%
  mutate(p = format(p, digits =2, scientific = FALSE))

## ----tris12confirm, fig.width=8, fig.height=12, warning=FALSE, fig.path=plotDir, dev=c("png", "pdf")----
pList <- group_by(plotTab, Drug) %>% 
  do(plots = ggplot(., aes(x=Concentration, y = viab)) + 
       stat_boxplot(geom = "errorbar", width = 0.3,
                    position = position_dodge(width=0.6), 
                    aes(dodge = Trisomy12)) +
       geom_boxplot(outlier.shape = NA, position = position_dodge(width=0.6), 
                    col="black", width=0.5, aes(dodge = Trisomy12)) + 
       geom_beeswarm(size=1,dodge.width=0.6, aes(col=Trisomy12)) +
       theme_classic() +
       scale_y_continuous(expand = c(0, 0),breaks=seq(0,1.2,0.2)) +
       coord_cartesian(ylim = c(0,1.3)) +
       xlab("Concentration (µM)") + ylab("Viability") + 
       ggtitle(unique(.$Drug)) +
       geom_text(data=filter(pTab, Drug == unique(.$Drug)), y = 1.25, 
                 aes(x=Concentration, label=sprintf("p=%s",p)), size = 5) + 
       theme(axis.line.x = element_blank(),
             axis.ticks.x = element_blank(),
             axis.text  = element_text(size=15),
             axis.title = element_text(size =15),
             legend.text = element_text(size=13),
             legend.title = element_text(size=15),
             plot.title = element_text(face="bold", hjust=0.5, size=17),
             plot.margin = unit(c(0.5,0,0.5,0), "cm"))) 

grid.arrange(grobs = pList$plots, ncol =1)

## ---- include=!exists(".standalone"), eval=!exists(".standalone")-------------
#  sessionInfo()

## ---- message=FALSE, warning=FALSE, include=FALSE-----------------------------
rm(list=ls())

## ---- echo=FALSE, include=!exists(".standalone")------------------------------
knitr::opts_chunk$set(cache = TRUE)

## ---- message=FALSE, include=!exists(".standalone"), eval=!exists(".standalone")----
#  library("BloodCancerMultiOmics2017")
#  library("DESeq2")
#  library("piano")
#  library("pheatmap")
#  library("genefilter")
#  library("grid")
#  library("gridExtra")
#  library("RColorBrewer")
#  library("cowplot")
#  library("dplyr")
#  library("ggplot2")
#  library("tibble")

## ----echo=FALSE---------------------------------------------------------------
plotDir = ifelse(exists(".standalone"), "", "part13/")
if(plotDir!="") if(!file.exists(plotDir)) dir.create(plotDir)

## -----------------------------------------------------------------------------
data(list=c("dds", "patmeta", "mutCOM"))

#load genesets
gmts = list(
  H=system.file("extdata","h.all.v5.1.symbols.gmt",
                package="BloodCancerMultiOmics2017"),
  C6=system.file("extdata","c6.all.v5.1.symbols.gmt",
                 package="BloodCancerMultiOmics2017"),
  KEGG=system.file("extdata","c2.cp.kegg.v5.1.symbols.gmt",
                   package="BloodCancerMultiOmics2017"))

## -----------------------------------------------------------------------------
#only choose CLL samples
colData(dds)$Diagnosis <- patmeta[match(dds$PatID,rownames(patmeta)),]$Diagnosis
ddsCLL <- dds[,dds$Diagnosis %in% "CLL"]

#add trisomy 12 and IGHV information
colData(ddsCLL)$trisomy12 <-
  factor(assayData(mutCOM[ddsCLL$PatID,])$binary[,"trisomy12"])
colData(ddsCLL)$IGHV <- factor(patmeta[ddsCLL$PatID,]$IGHV)

#remove samples that do not have trisomy 12 information
ddsCLL <- ddsCLL[,!is.na(ddsCLL$trisomy12)]

#how many genes and samples we have?
dim(ddsCLL)

## ---- cache=TRUE--------------------------------------------------------------
#remove genes without gene symbol annotations
ddsCLL <- ddsCLL[!is.na(rowData(ddsCLL)$symbol),]
ddsCLL <- ddsCLL[rowData(ddsCLL)$symbol != "",]

#only keep genes that have counts higher than 10 in any sample
keep <- apply(counts(ddsCLL), 1, function(x) any(x >= 10)) 
ddsCLL <- ddsCLL[keep,]

#Remove transcripts do not show variance across samples
ddsCLL <- estimateSizeFactors(ddsCLL)
sds <- rowSds(counts(ddsCLL, normalized = TRUE))
sh <- shorth(sds)
ddsCLL <- ddsCLL[sds >= sh,]

#variance stabilization
ddsCLL.norm <- varianceStabilizingTransformation(ddsCLL, blind=TRUE)

#how many genes left
dim(ddsCLL)

## ---- cache=TRUE--------------------------------------------------------------
design(ddsCLL) <- ~ trisomy12
ddsCLL <- DESeq(ddsCLL, betaPrior = FALSE)
DEres <- results(ddsCLL)
DEres.shr <- lfcShrink(ddsCLL, type="normal", contrast = c("trisomy12","1","0"),
                       res = DEres)

## ---- warning=FALSE-----------------------------------------------------------
#FIG# S23 A
plotTab <- as.data.frame(DEres)
plotTab$onChr12 <- rowData(ddsCLL)$chromosome == 12
dosePlot <- ggplot(plotTab) +
  geom_density(aes(x=log2FoldChange, col=onChr12, fill=onChr12), alpha=0.4) +
  xlim( -3, 3 )
dosePlot

## -----------------------------------------------------------------------------
#filter genes
fdrCut <- 0.1
fcCut <- 1.5

allDE <- data.frame(DEres.shr) %>%
  rownames_to_column(var = "ID") %>% 
  mutate(Symbol = rowData(ddsCLL[ID,])$symbol,
         Chr = rowData(ddsCLL[ID,])$chromosome) %>% 
  filter(padj <= fdrCut & abs(log2FoldChange) > fcCut) %>% 
  arrange(pvalue) %>% filter(!duplicated(Symbol)) %>%
  mutate(Chr12 = ifelse(Chr == 12, "yes", "no"))

#get the expression matrix
plotMat <- assay(ddsCLL.norm[allDE$ID,])
colnames(plotMat) <- ddsCLL.norm$PatID
rownames(plotMat) <- allDE$Symbol

#sort columns of plot matrix based on trisomy 12 status
plotMat <- plotMat[,order(ddsCLL.norm$trisomy12)]

#calculate z-score and scale
plotMat <- t(scale(t(plotMat)))
plotMat[plotMat >= 4] <- 4
plotMat[plotMat <= -4] <- -4

## ----  trisomy12_heatmap, dev = c("png", "pdf"), fig.path=plotDir, fig.width = 8, fig.height = 10----
#FIG# S23 B
#prepare colums and row annotations
annoCol <- data.frame(row.names=ddsCLL.norm$PatID, Tris12=ddsCLL.norm$trisomy12)
levels(annoCol$Tris12) <- list(wt = 0, mut =1)
annoRow <- data.frame(row.names = allDE$Symbol, Chr12 = allDE$Chr12)
annoColor <- list(Tris12 = c(wt = "grey80", mut = "black"),
                  Chr12 = c(yes="red", no = "grey80"))


pheatmap(plotMat,
         color=colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100),
         cluster_cols = FALSE,
         annotation_row = annoRow, annotation_col = annoCol,
         show_colnames = FALSE, fontsize_row = 3,
         breaks = seq(-4,4, length.out = 101),
         annotation_colors = annoColor, border_color = NA)


## -----------------------------------------------------------------------------
runGSEA <- function(inputTab, gmtFile, GSAmethod="gsea", nPerm=1000){
  inGMT <- loadGSC(gmtFile,type="gmt")
  #re-rank by score
  rankTab <- inputTab[order(inputTab[,1],decreasing = TRUE),,drop=FALSE] 
  if (GSAmethod == "gsea"){
    #readin geneset database
    #GSEA analysis
    res <- runGSA(geneLevelStats = rankTab,
                  geneSetStat = GSAmethod,
                  adjMethod = "fdr", gsc=inGMT,
                  signifMethod = 'geneSampling', nPerm = nPerm)
    GSAsummaryTable(res)
  } else if (GSAmethod == "page"){
    res <- runGSA(geneLevelStats = rankTab,
                  geneSetStat = GSAmethod,
                  adjMethod = "fdr", gsc=inGMT,
                  signifMethod = 'nullDist')
    GSAsummaryTable(res)
  }
}

## -----------------------------------------------------------------------------
plotEnrichmentBar <- function(resTab, pCut=0.05, ifFDR=FALSE,
                              setName="Signatures") {
  pList <- list()
  rowNum <- c()
  for (i in names(resTab)) {
    plotTab <- resTab[[i]]
    if (ifFDR) {
      plotTab <- dplyr::filter(
        plotTab, `p adj (dist.dir.up)` <= pCut | `p adj (dist.dir.dn)` <= pCut)
    } else {
      plotTab <- dplyr::filter(
        plotTab, `p (dist.dir.up)` <= pCut | `p (dist.dir.dn)` <= pCut)
    }
    if (nrow(plotTab) == 0) {
      print("No sets passed the criteria")
      next
    } else {
      #firstly, process the result table
      plotTab <- apply(plotTab, 1, function(x) {
        statSign <- as.numeric(x[3])
        data.frame(Name = x[1],
                   p = as.numeric(ifelse(statSign >= 0, x[4], x[6])),
                   geneNum = ifelse(statSign >= 0, x[8], x[9]),
                   Direction = ifelse(statSign > 0, "Up", "Down"),
                   stringsAsFactors = FALSE)
      }) %>% do.call(rbind,.)
      
      plotTab$Name <- sprintf("%s (%s)",plotTab$Name,plotTab$geneNum)
      plotTab <- plotTab[with(plotTab,order(Direction, p, decreasing=TRUE)),]
      plotTab$Direction <- factor(plotTab$Direction, levels = c("Down","Up"))
      plotTab$Name <- factor(plotTab$Name, levels = plotTab$Name)
      #plot the barplot
      pList[[i]] <- ggplot(data=plotTab, aes(x=Name, y= -log10(p),
                                             fill=Direction)) +
        geom_bar(position="dodge",stat="identity", width = 0.5) +
        scale_fill_manual(values=c(Up = "blue", Down = "red")) +
        coord_fixed(ratio = 0.5) + coord_flip() + xlab(setName) +
        ggtitle(i) + theme_bw() + theme(
          plot.title = element_text(face = "bold", hjust =0.5),
          axis.title = element_text(size=15))
      rowNum <-c(rowNum,nrow(plotTab))
    }
  }
  
  if (length(pList) == 0) {
    print("Nothing to plot")
  } else {
    rowNum <- rowNum
    grobList <- lapply(pList, ggplotGrob)
    grobList <- do.call(rbind,c(grobList,size="max"))
    panels <- grobList$layout$t[grep("panel", grobList$layout$name)]
    grobList$heights[panels] <- unit(rowNum, "null")
  }
  return(grobList)
}

## ----message=FALSE------------------------------------------------------------
pCut <- 0.05

dataTab <- data.frame(DEres)
dataTab$ID <- rownames(dataTab)

#filter using raw pvalues
dataTab <- filter(dataTab, pvalue <= pCut) %>%
  arrange(pvalue) %>%
  mutate(Symbol = rowData(ddsCLL[ID,])$symbol)
dataTab <- dataTab[!duplicated(dataTab$Symbol),]
statTab <- data.frame(row.names = dataTab$Symbol, stat = dataTab$stat)

## ----fig.width=8, fig.height=6, message=FALSE---------------------------------
hallmarkRes <- list()

#run PAGE
resTab <- runGSEA(statTab, gmts$H ,GSAmethod = "page")

#remove the HALLMARK_
resTab$Name <- gsub("HALLMARK_","",resTab$Name)

hallmarkRes[["Gene set enrichment analysis"]] <- 
  arrange(resTab,desc(`Stat (dist.dir)`))

hallBar <- plotEnrichmentBar(hallmarkRes, pCut = 0.01, ifFDR = TRUE,
                             setName = "Hallmark gene sets")

## ---- message=FALSE-----------------------------------------------------------
keggRes <- list()

resTab <- runGSEA(statTab,gmts$KEGG,GSAmethod = "page")

#remove the KEGG_
resTab$Name <- gsub("KEGG_","",resTab$Name)

keggRes[["Gene set enrichment analysis"]] <- resTab

keggBar <- plotEnrichmentBar(keggRes, pCut = 0.01, ifFDR = TRUE,
                             setName = "KEGG gene sets")

## -----------------------------------------------------------------------------
#select differentially expressed genes
fdrCut <- 0.05
cytoDE <- data.frame(DEres) %>% rownames_to_column(var = "ID") %>% 
  mutate(Symbol = rowData(ddsCLL[ID,])$symbol,
         Chr=rowData(ddsCLL[ID,])$chromosome) %>% 
  filter(padj <= fdrCut, log2FoldChange > 0) %>% 
  arrange(pvalue) %>% filter(!duplicated(Symbol)) %>%
  mutate(Chr12 = ifelse(Chr == 12, "yes", "no"))

#get the expression matrix
plotMat <- assay(ddsCLL.norm[cytoDE$ID,])
colnames(plotMat) <- ddsCLL.norm$PatID
rownames(plotMat) <- cytoDE$Symbol

#sort columns of plot matrix based on trisomy 12 status
plotMat <- plotMat[,order(ddsCLL.norm$trisomy12)]

#calculate z-score and sensor
plotMat <- t(scale(t(plotMat)))
plotMat[plotMat >= 4] <- 4
plotMat[plotMat <= -4] <- -4

annoCol <- data.frame(row.names = ddsCLL.norm$PatID,
                      Tris12 = ddsCLL.norm$trisomy12)
levels(annoCol$Tris12) <- list(wt = 0, mut =1)
annoRow <- data.frame(row.names = cytoDE$Symbol, Chr12 = cytoDE$Chr12)


## -----------------------------------------------------------------------------
gsc <- loadGSC(gmts$KEGG)
geneList <- gsc$gsc$KEGG_CHEMOKINE_SIGNALING_PATHWAY
plotMat.chemo <- plotMat[rownames(plotMat) %in% geneList,]
keggHeatmap <- pheatmap(plotMat.chemo,
                        color = colorRampPalette(
                          rev(brewer.pal(n=7, name="RdBu")))(100),
                        cluster_cols = FALSE, clustering_method = "ward.D2",
                        annotation_row = annoRow, annotation_col = annoCol,
                        show_colnames = FALSE, fontsize_row = 8,
                        breaks = seq(-4,4, length.out = 101), treeheight_row = 0,
                        annotation_colors = annoColor, border_color = NA,
                        main = "CHEMOKINE_SIGNALING_PATHWAY",silent = TRUE)$gtable

## ----geneEnrichment_result, dev = c("png", "pdf"), fig.path=plotDir, fig.width = 14, fig.height = 13----
#FIG# S24 ABC
ggdraw() + 
  draw_plot(hallBar, 0, 0.7, 0.5, 0.3) + 
  draw_plot(keggBar, 0.5, 0.7, 0.5, 0.3) +
  draw_plot(keggHeatmap, 0.1, 0, 0.8, 0.65) +
  draw_plot_label(c("A","B","C"), c(0, 0.5, 0.1), c(1, 1, 0.7), 
                  fontface = "plain", size=20)


## ---- include=!exists(".standalone"), eval=!exists(".standalone")-------------
#  sessionInfo()

## ---- message=FALSE, warning=FALSE, include=FALSE-----------------------------
rm(list=ls())

## ---- message=FALSE, include=!exists(".standalone"), eval=!exists(".standalone")----
#  library("BloodCancerMultiOmics2017")
#  library("Biobase")
#  library("SummarizedExperiment")
#  library("AnnotationDbi")
#  library("org.Hs.eg.db")
#  library("dplyr")
#  library("abind")
#  library("reshape2")
#  library("RColorBrewer")
#  library("glmnet")
#  library("ipflasso")
#  library("ggplot2")
#  library("grid")
#  library("DESeq2")

## ----echo=FALSE---------------------------------------------------------------
plotDir = ifelse(exists(".standalone"), "", "part11/")
if(plotDir!="") if(!file.exists(plotDir)) dir.create(plotDir)

## -----------------------------------------------------------------------------
options(stringsAsFactors=FALSE)

## -----------------------------------------------------------------------------
data(list=c("conctab", "drpar", "drugs", "patmeta", "lpdAll", "dds", "mutCOM",
            "methData"))

## -----------------------------------------------------------------------------
e<-dds
colnames(e)<-colData(e)$PatID


#only consider CLL patients
CLLPatients<-rownames(patmeta)[which(patmeta$Diagnosis=="CLL")]

#Methylation Data
methData = t(assay(methData)) 


#RNA Data
eCLL<-e[,colnames(e) %in% CLLPatients]
###
#filter out genes without gene namce
AnnotationDF<-data.frame(EnsembleId=rownames(eCLL),stringsAsFactors=FALSE)
AnnotationDF$symbol <- mapIds(org.Hs.eg.db,
                              keys=rownames(eCLL),
                              column="SYMBOL",
                              keytype="ENSEMBL",
                              multiVals="first")
eCLL<-eCLL[AnnotationDF$EnsembleId[!is.na(AnnotationDF$symbol)],]

#filter out low count genes
###
minrs <- 100
rs  <- rowSums(assay(eCLL))
eCLL<-eCLL[ rs >= minrs, ]
#variance stabilize the data
#(includes normalizing for library size and dispsersion estimation) 
vstCounts<-varianceStabilizingTransformation(eCLL)
vstCounts<-assay(vstCounts)
#no NAs in data
any(is.na(vstCounts))

#filter out low variable genes
ntop<-5000
vstCountsFiltered<-vstCounts[order(apply(vstCounts, 1, var, na.rm=T),
                                   decreasing = T)[1:ntop],]
eData<-t(vstCountsFiltered)
#no NAs
any(is.na(eData))

#genetics
#remove features with less than 5 occurences
mutCOMbinary<-channel(mutCOM, "binary")
mutCOMbinary<-mutCOMbinary[featureNames(mutCOMbinary) %in% CLLPatients,]
genData<-Biobase::exprs(mutCOMbinary)
idx <- which(colnames(genData) %in% c("del13q14_bi", "del13q14_mono"))
genData <- genData[,-idx]
colnames(genData)[which(colnames(genData)=="del13q14_any")] = "del13q14"
minObs <- 5
#remove feutes with less than 5 occurecnes
genData<-genData[,colSums(genData, na.rm=T)>=minObs]

#IGHV
translation <- c(`U` = 0, `M` = 1)
stopifnot(all(patmeta$IGHV %in% c("U","M", NA)))
IGHVData <- matrix(translation[patmeta$IGHV], 
                   dimnames = list(rownames(patmeta), "IGHV"), ncol = 1)
IGHVData<-IGHVData[rownames(IGHVData) %in% CLLPatients,,drop=F]
#remove patiente with NA IGHV status
IGHVData<-IGHVData[!is.na(IGHVData), ,drop=F]
any(is.na(IGHVData))

#demographics (age and sex)
patmeta<-subset(patmeta, Diagnosis=="CLL")
gender <- ifelse(patmeta[,"Gender"]=="m",0,1)


# impute missing values in age by mean
ImputeByMean <- function(x) {x[is.na(x)] <- mean(x, na.rm=TRUE); return(x)}
age<-ImputeByMean(patmeta[,"Age4Main"])


demogrData <- cbind(age=age,gender=gender)
rownames(demogrData) <- rownames(patmeta)

#Pretreatment
pretreated<-patmeta[,"IC50beforeTreatment", drop=FALSE]

##### drug viabilites
summaries <- c(paste("viaraw", 1:5, sep=".") %>% `names<-`(paste(1:5)), 
               `4:5` = "viaraw.4_5", `1:5` = "viaraw.1_5")
a <- do.call( abind, c( along=3, lapply( summaries, 
                                         function(x) assayData(drpar)[[x]]))) 
dimnames(a)[[3]] <- names(summaries)
names(dimnames(a)) <- c( "drug", "patient", "summary" )
viabData <- acast( melt(a), patient ~ drug + summary )
rownames(viabData)<-c(substr(rownames(viabData),1,4)[1:3],
                      substr(rownames(viabData),1,5)[4:nrow(viabData)])

## -----------------------------------------------------------------------------
# common patients 
Xlist<-list(RNA=eData, meth=methData, gen=genData, IGHV=IGHVData,
            demographics=demogrData, drugs=viabData, pretreated=pretreated)
PatientsPerOmic<-lapply(Xlist, rownames)
sapply(PatientsPerOmic, length)

allPatients<-Reduce(union, PatientsPerOmic)
PatientOverview<-sapply(Xlist, function(M) allPatients %in% rownames(M))
Patients <- (1:nrow(PatientOverview))
Omics <- (1:ncol(PatientOverview))
image(Patients,Omics, PatientOverview*1, axes=F, col=c("white", "black"),
      main="Sample overview across omics")
axis(2, at = 1:ncol(PatientOverview), labels=colnames(PatientOverview), tick=F)

commonPatients<-Reduce(intersect, PatientsPerOmic)
length(commonPatients)
XlistCommon<-lapply(Xlist, function(data) data[commonPatients,, drop=F])

#Take care of missing values (present in  genetic data)
ImputeByMean <- function(x) {x[is.na(x)] <- mean(x, na.rm=TRUE); return(x)}

#NAs in genetic
#remove feauters with less 90% completeness
RarlyMeasuredFeautres<-
  which(colSums(is.na(XlistCommon$gen))>0.1*nrow(XlistCommon$gen))
XlistCommon$gen<-XlistCommon$gen[,-RarlyMeasuredFeautres]
#remove patients with less than 90% of genetic feautres measured
IncompletePatients<-
  rownames(XlistCommon$gen)[
    (rowSums(is.na(XlistCommon$gen))>0.1*ncol(XlistCommon$gen))]
commonPatients<-commonPatients[!commonPatients %in% IncompletePatients]
XlistCommon<-lapply(XlistCommon, function(data) data[commonPatients,, drop=F])
#replace remaining NA by mean and round to 0 or 1
XlistCommon$gen<-round(apply(XlistCommon$gen, 2, ImputeByMean))

#NAs in methylation
#remove feauters with less 90% completeness
XlistCommon$meth<-
  XlistCommon$meth[,colSums(is.na(XlistCommon$meth))<0.1*nrow(methData)]
#impute remainin missing values by mean for each feautre across patients
XlistCommon$meth<-(apply(XlistCommon$meth, 2, ImputeByMean))

#final dimensions of the data
sapply(XlistCommon, dim)

## ---- fig.width=12, fig.height=10---------------------------------------------
pcaMeth<-prcomp(XlistCommon$meth, center=T, scale. = F)
XlistCommon$MethPCs<-pcaMeth$x[,1:20]
colnames(XlistCommon$MethPCs)<-
  paste("meth",colnames(XlistCommon$MethPCs), sep="")

pcaExpr<-prcomp(XlistCommon$RNA, center=T, scale. = F)
XlistCommon$RNAPCs<-pcaExpr$x[,1:20]
colnames(XlistCommon$RNAPCs)<-paste("RNA",colnames(XlistCommon$RNAPCs), sep="")

## -----------------------------------------------------------------------------
DOI <- c("D_006_1:5", "D_010_1:5", "D_159_1:5","D_002_4:5", "D_003_4:5",
         "D_012_4:5", "D_063_4:5", "D_166_4:5")
drugviab<-XlistCommon$drugs
drugviab<-drugviab[,DOI, drop=F]
colnames(drugviab) <- drugs[substr(colnames(drugviab),1,5),"name"]

## -----------------------------------------------------------------------------
ZPCs<-list(expression=XlistCommon$RNAPCs,
           genetic=XlistCommon$gen, 
           methylation= XlistCommon$MethPCs,
           demographics=XlistCommon$demographics, 
           IGHV=XlistCommon$IGHV,
           pretreated = XlistCommon$pretreated)
ZPCs$all<-do.call(cbind, ZPCs)
ZPCsunscaled<-ZPCs
ZPCsscaled<-lapply(ZPCs, scale)
lapply(ZPCsscaled, colnames)

## -----------------------------------------------------------------------------
set1 <- brewer.pal(9,"Set1")
colMod<-c(paste(set1[c(4,1,5,3,2,7)],"88",sep=""), "grey")
names(colMod) <-
  c("demographics", "genetic", "IGHV","expression", "methylation", "pretreated",
    "all")

## ---- echo=F------------------------------------------------------------------
#Function to calculate Var Explained for Penalized Regression
R2ForPenRegPlusViz<-function(Z, drugviab, nfolds=10, alpha=1, nrep=100,
                             Parmfrow=c(2,4), ylimMax=0.4, standardize=TRUE){
  Zlocal<-Z
  
  set.seed(1030)
  seeds<-sample(1:10000000, nrep)
  
  RepeatedR2list<-lapply(1:nrep, function(outer){
    #Use same folds for all omics to make comparable
    set.seed(seeds[outer])
    foldsAssignment<-sample(rep(seq(nfolds), length=nrow(drugviab)))
    
    R2echOmicadj<-sapply(colnames(drugviab), function(dname) {
      d<-drugviab[,dname]
      sapply(names(Zlocal), function(nameZ) {
        pred<-Zlocal[[nameZ]]
        #fit a lasso model for omics with more than one features
        if(ncol(pred)>1){
          fitcv<-cv.glmnet(pred,d,family="gaussian", standardize=standardize,
                           alpha=alpha, foldid=foldsAssignment)
          R2 <- 1-min(fitcv$cvm)/fitcv$cvm[1]
          
        }else {#fit a liner model for single feautres (IGHV)
          fitlm<-lm(d~., data.frame(pred))
          R2<- summary(fitlm)$r.squared
        }
        R2
      })
    }
    )
  })
  RepeatedR2<-RepeatedR2list
  #calculate mean and sd across repitions wiht different folds
  meanR2<-apply(simplify2array(RepeatedR2), 1:2, mean)
  sdR2<-apply(simplify2array(RepeatedR2), 1:2, sd)
  
  par(mfrow=Parmfrow, mar=c(7,5,3,3))
  for (i in 1: ncol(meanR2)) {
    barc<-barplot(meanR2[,i], main= colnames(meanR2)[i], ylim=c(0,ylimMax),
                  las=2, col=colMod[rownames(meanR2)], ylab="R2")
    segments(barc, meanR2[,i]-sdR2[,i],barc, meanR2[,i]+sdR2[,i])
  }
  RepeatedR2
}

## ----lasso_main, echo=FALSE, fig.path=plotDir, dev=c("png", "pdf"), fig.width=9, fig.height=8----
#FIG# 5A
resultLassoPCA<-R2ForPenRegPlusViz(ZPCsscaled, drugviab, nfold=10, alpha=1,
                                   nrep=100, ylimMax=0.4)

df_resultLassoPCA <- melt(resultLassoPCA)
colnames(df_resultLassoPCA) <- c("omic", "drug", "R2", "run")
summaryR2 <- df_resultLassoPCA %>% group_by(omic, drug) %>% 
  dplyr::summarise(meanR2=mean(R2),sdR2 = sd(R2), nR2 = length(R2)) %>%
  mutate(seR2 = sdR2/sqrt(nR2))

ggplot(summaryR2, aes(x=omic, y=meanR2, fill=omic, group=omic))+
  geom_bar(stat = "identity") +  scale_fill_manual(values=colMod) + 
  geom_errorbar(aes(ymax = meanR2 + sdR2,ymin = meanR2 - sdR2), position = "dodge", width = 0.25) +facet_wrap(~drug, ncol=4) +theme_bw(base_size = 18) +ylab(bquote(R^2)) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(colour = "black")) +xlab("")+guides(fill=guide_legend(title="Data type")) 

## ---- eval=FALSE, fig.path='supp_'--------------------------------------------
#  nfolds<-10
#  nrep<-100
#  
#  DOI <-
#    grepl("1:5",colnames(XlistCommon$drugs)) |
#    grepl("4:5",colnames(XlistCommon$drugs))
#  drugviabAll<-XlistCommon$drugs
#  drugviabAll<-drugviabAll[,DOI]
#  colnames(drugviabAll) <-
#    paste0(drugs[substr(colnames(drugviabAll),1,5),"name"],
#           substr(colnames(drugviabAll),6,9))
#  
#  R2ForPenReg(Zscaled, drugviabAll, nfold=nfolds, alpha=1, nrep=nrep,
#              Parmfrow=c(4,4), ylimMax=0.6)

## ----echo=FALSE---------------------------------------------------------------
dataType = c(M="Methylation_Cluster", V="viab", G="gen", I="IGHV", P="pretreat")

## ----echo=FALSE---------------------------------------------------------------
coldef<-list()
coldef["I"]<-brewer.pal(9, "Blues")[7]
coldef["M"]<-list(brewer.pal(9, "Blues")[c(1, 5, 9)])
coldef["G"]<-brewer.pal(8, "YlOrRd")[8]
coldef["P"]<-"chocolate4"

## -----------------------------------------------------------------------------
lpdCLL = lpdAll[ , lpdAll$Diagnosis=="CLL"]

## -----------------------------------------------------------------------------
lpdCLL = lpdAll[ , lpdAll$Diagnosis=="CLL"]
lpdCLL = BloodCancerMultiOmics2017:::prepareLPD(lpd=lpdCLL, minNumSamplesPerGroup=5)
(predictorList = BloodCancerMultiOmics2017:::makeListOfPredictors(lpdCLL))

## -----------------------------------------------------------------------------
drs = list("1:5"=c("D_006", "D_010", "D_159"),
           "4:5"=c("D_002", "D_003", "D_012", "D_063", "D_166"))

## -----------------------------------------------------------------------------
predvar = unlist(BloodCancerMultiOmics2017:::makePredictions(drs=drs,
                                                             lpd=lpdCLL,
                                                             predictorList=predictorList,
                                                             lim=0.15, std=FALSE, adaLasso=TRUE,
                                                             colors=coldef),
                 recursive=FALSE)

## ----echo=FALSE---------------------------------------------------------------
details = function(dr, what) {
  predvar[[dr]][["plot"]][[what]]
}

## ----prediction-D_006, echo=FALSE, fig.path=plotDir, dev=c("png", "pdf"), fig.width=details("D_006","width"), fig.height=details("D_006","height"), eval=TRUE----
#FIG# 5B
grid.draw(details("D_006","plot"))

## ----prediction-D_010, echo=FALSE, fig.path=plotDir, dev=c("png", "pdf"), fig.width=details("D_010","width"), fig.height=details("D_010","height"), eval=TRUE----
#FIG# 5B
grid.draw(details("D_010","plot"))

## ----prediction-D_159, echo=FALSE, fig.path=plotDir, dev=c("png", "pdf"), fig.width=details("D_159","width"), fig.height=details("D_159","height"), eval=TRUE----
#FIG# 5B
grid.draw(details("D_159","plot"))

## ----prediction-D_002, echo=FALSE, fig.path=plotDir, dev=c("png", "pdf"), fig.width=details("D_002","width"), fig.height=details("D_002","height"), eval=TRUE----
#FIG# 5B
grid.draw(details("D_002","plot"))

## ----prediction-D_003, echo=FALSE, fig.path=plotDir, dev=c("png", "pdf"), fig.width=details("D_003","width"), fig.height=details("D_003","height"), eval=TRUE----
#FIG# 5B
grid.draw(details("D_003","plot"))

## ----prediction-D_012, echo=FALSE, fig.path=plotDir, dev=c("png", "pdf"), fig.width=details("D_012","width"), fig.height=details("D_012","height"), eval=TRUE----
#FIG# 5B
grid.draw(details("D_012","plot"))

## ----prediction-D_063, echo=FALSE, fig.path=plotDir, dev=c("png", "pdf"), fig.width=details("D_063","width"), fig.height=details("D_063","height"), eval=TRUE----
#FIG# 5B
grid.draw(details("D_063","plot"))

## ----prediction-D_166, echo=FALSE, fig.path=plotDir, dev=c("png", "pdf"), fig.width=details("D_166","width"), fig.height=details("D_166","height"), eval=TRUE----
#FIG# 5B
grid.draw(details("D_166","plot"))

## -----------------------------------------------------------------------------
legends = BloodCancerMultiOmics2017:::makeLegends(legendFor=c("G","I","M", "P"),
                                                  coldef)

## ----legend, echo=FALSE, fig.path=plotDir, dev=c("png", "pdf"), fig.width=legends[["width"]], fig.height=legends[["height"]]----
#FIG# 5B legend
grid.draw(legends[["plot"]])

## -----------------------------------------------------------------------------
drs_rot = list("4:5"=c("D_067"))
predvar_rot = unlist(BloodCancerMultiOmics2017:::makePredictions(drs=drs_rot,
                                                                 lpd=lpdCLL,
                                                                 predictorList=predictorList,
                                                                 lim=0.23, std=FALSE, adaLasso=TRUE,
                                                                 colors=coldef),
                     recursive=FALSE)

## ----prediction-D_067, echo=FALSE, fig.path=plotDir, dev=c("png", "pdf"), fig.width=predvar_rot[["D_067"]][["plot"]][["width"]], fig.height=predvar_rot[["D_067"]][["plot"]][["height"]], eval=TRUE----
#FIG# S26
grid.draw(predvar_rot[["D_067"]][["plot"]][["plot"]])

## -----------------------------------------------------------------------------
alldrs = unique(fData(lpdCLL)[fData(lpdCLL)$type=="viab","id"])
drs = list("1:5"=alldrs, "4:5"=alldrs)

## ---- eval=TRUE---------------------------------------------------------------
predvar2 = BloodCancerMultiOmics2017:::makePredictions(drs=drs,
                                                       lpd=lpdCLL,
                                                       predictorList=predictorList,
                                                       lim=0.23,
                                                       colors=coldef)

## -----------------------------------------------------------------------------
givePreatreatSum = function(predNum) {
  
  idx = sapply(predvar2[[predNum]], function(x) length(x)==1)
  predvar2[[predNum]] = predvar2[[predNum]][!idx]
  # get model coefficients and reshape
  coeffs <- do.call(cbind,lapply(predvar2[[predNum]], "[[", 'coeffs'))
  coeffs <- coeffs[-1,]
  coeffs <- as.matrix(coeffs)
  # colnames(coeffs) <- unlist(drs["1:5"])
  colnames(coeffs) = names(predvar2[[predNum]])
  colnames(coeffs) <- drugs[colnames(coeffs),"name"]
  coeffDF <- melt(as.matrix(coeffs))
  colnames(coeffDF) <- c("predictor", "drug", "beta")
  coeffDF$selected <- coeffDF$beta !=0
  
  #sort by times being selected
  coeffDF$predictor <- factor(coeffDF$predictor, level=)
  
  # number of drugs a predictor is chosen for
  gg1 <- coeffDF %>% group_by(predictor) %>% 
    dplyr::summarize(selectedSum = sum(selected)) %>%
    mutate(predictor = factor(predictor,
                              levels=predictor[order(selectedSum)])) %>%
    ggplot(aes(x=predictor, y=selectedSum)) + 
    geom_bar(stat="identity")+ylab("# drugs selected for") +
    coord_flip()
  
  # boxplots of non-zero coeffients
  orderMedian <- filter(coeffDF, selected) %>% group_by(predictor) %>% 
    dplyr::summarize(medianBeta = median(abs(beta)))
  coeffDF$predictor <- factor(
    coeffDF$predictor,
    levels=orderMedian$predictor[order(orderMedian$medianBeta)] )
  gg2 <- ggplot(filter(coeffDF, selected), aes(x=predictor, y=abs(beta))) +
    geom_boxplot() +
    coord_flip() + ggtitle("Distribution of non-zero coefficients")
  gridExtra::grid.arrange(gg1,gg2, ncol=1)
  
  # coefficeints per drug
  ggplot(filter(coeffDF, selected), 
         aes(x= drug, y=abs(beta), col= predictor=="Pretreatment")) +
    geom_point() +
    coord_flip()
  
  #drugs pretreatment is selected for
  as.character(filter(coeffDF, predictor=="Pretreatment" & beta!=0)$drug)
  PselDrugs <- as.character(
    filter(coeffDF, predictor=="Pretreatment" & beta!=0)$drug)
  length(PselDrugs)
  # length(drs[[1]])
}

## ----pretreatment_c1-5, echo=FALSE, fig.path=plotDir, dev=c("png", "pdf"), fig.width=7, fig.height=10, eval=TRUE----
givePreatreatSum(predNum=1)

## ----pretreatment_c4-5, echo=FALSE, fig.path=plotDir, dev=c("png", "pdf"), fig.width=7, fig.height=10, eval=TRUE----
givePreatreatSum(predNum=2)

## ---- include=!exists(".standalone"), eval=!exists(".standalone")-------------
#  sessionInfo()

## ---- message=FALSE, warning=FALSE, include=FALSE-----------------------------
rm(list=ls())

## ---- message=FALSE, include=!exists(".standalone"), eval=!exists(".standalone")----
#  library("BloodCancerMultiOmics2017")
#  library("Biobase")
#  library("RColorBrewer")
#  library("grid")
#  library("ggplot2")
#  library("survival")
#  library("gtable")
#  library("forestplot")
#  library("xtable")
#  library("maxstat")

## ----echo=FALSE---------------------------------------------------------------
plotDir = ifelse(exists(".standalone"), "", "part08/")
if(plotDir!="") if(!file.exists(plotDir)) dir.create(plotDir)

## -----------------------------------------------------------------------------
options(stringsAsFactors=FALSE)

## -----------------------------------------------------------------------------
data(lpdAll, patmeta, drugs)

## -----------------------------------------------------------------------------
lpdCLL <- lpdAll[ , lpdAll$Diagnosis=="CLL"   ]

# data rearrangements
survT = patmeta[colnames(lpdCLL),]
survT[which(survT[,"IGHV"]=="U") ,"IGHV"] = 0
survT[which(survT[,"IGHV"]=="M") ,"IGHV"] = 1
survT$IGHV = as.numeric(survT$IGHV)

colnames(survT) = gsub("Age4Main", "age", colnames(survT))

survT$ibr45  <- 1-Biobase::exprs(lpdCLL)[ "D_002_4:5", rownames(survT)  ]  
survT$ide45  <- 1-Biobase::exprs(lpdCLL)[ "D_003_4:5", rownames(survT)  ] 
survT$prt45  <- 1-Biobase::exprs(lpdCLL)[ "D_166_4:5", rownames(survT)  ] 
survT$selu45 <- 1-Biobase::exprs(lpdCLL)[ "D_012_4:5", rownames(survT)  ]
survT$ever45 <- 1-Biobase::exprs(lpdCLL)[ "D_063_4:5", rownames(survT)  ]

survT$nut15 <- 1-Biobase::exprs(lpdCLL)[ "D_010_1:5", rownames(survT)  ] 
survT$dox15 <- 1-Biobase::exprs(lpdCLL)[ "D_159_1:5", rownames(survT)  ] 
survT$flu15 <- 1-Biobase::exprs(lpdCLL)[ "D_006_1:5", rownames(survT)  ] 

survT$SF3B1      <- Biobase::exprs(lpdCLL)[ "SF3B1",      rownames(survT)  ]
survT$NOTCH1     <- Biobase::exprs(lpdCLL)[ "NOTCH1",     rownames(survT)  ]
survT$BRAF       <- Biobase::exprs(lpdCLL)[ "BRAF",       rownames(survT)  ]
survT$TP53       <- Biobase::exprs(lpdCLL)[ "TP53",       rownames(survT)  ]
survT$del17p13   <- Biobase::exprs(lpdCLL)[ "del17p13",   rownames(survT)  ]
survT$del11q22.3 <- Biobase::exprs(lpdCLL)[ "del11q22.3", rownames(survT)  ]
survT$trisomy12 <-  Biobase::exprs(lpdCLL)[ "trisomy12", rownames(survT)  ]
survT$IGHV_cont <- patmeta[ rownames(survT) ,"IGHV Uppsala % SHM"]


# competinting risk endpoint fpr 
survT$compE <- ifelse(survT$treatedAfter == TRUE, 1, 0)
survT$compE <- ifelse(survT$treatedAfter == FALSE & survT$died==TRUE,
                      2, survT$compE )
survT$T7  <- ifelse(survT$compE == 1, survT$T5, survT$T6 )

## ----forest-------------------------------------------------------------------
forest <- function(Time, endpoint, title, sdrugs, split, sub) {  
  stopifnot(is.character(Time), is.character(title), is.character(split),
            is.character(endpoint), 
            all(c(Time, split, endpoint) %in% colnames(survT)),
            is.logical(survT[[endpoint]]),
            is.character(sdrugs), !is.null(names(sdrugs)))
  
  clrs <- fpColors(box="royalblue",line="darkblue", summary="royalblue")
  
  res <- lapply(sdrugs, function(g) { 
    drug <- survT[, g] * 10
    
    suse <- if (identical(sub, "none")) 
      rep(TRUE, nrow(survT)) 
    else 
      (survT[[split]] == sub)
    stopifnot(sum(suse, na.rm = TRUE) > 1)
    
    surv <- coxph(Surv(survT[,Time], survT[,endpoint]) ~ drug, subset=suse)  
    sumsu <- summary(surv) 
    c(p      = sumsu[["coefficients"]][, "Pr(>|z|)"], 
      coef   = sumsu[["coefficients"]][, "exp(coef)"], 
      lower  = sumsu[["conf.int"]][, "lower .95"], 
      higher = sumsu[["conf.int"]][, "upper .95"])
  })
  
  s <- do.call(rbind, res)
  rownames(s) <- names(sdrugs)
  
  tabletext <- list(c(NA, rownames(s)), append(list("p-value"),
                                               sprintf("%.4f", s[,"p"])))
  
  forestplot(tabletext, 
             rbind(
               rep(NA, 3), 
               s[, 2:4]), 
             page = new,
             clip = c(0.8,20), 
             xlog = TRUE, xticks = c(0.5,1, 1.5), title = title,
             col = clrs, 
             txt_gp = fpTxtGp(ticks = gpar(cex=1) ),
             new_page = TRUE)
}

## ----forest-together----------------------------------------------------------

com <- function( Time, endpoint, scaleX, sub, d, split, drug_names) {  
  
  res <- lapply(d, function(g)  { 
    
    drug <- survT[,g] * scaleX
    ## all=99, M-CLL=1, U-CLL=0
    if(sub==99) { surv <- coxph(Surv(survT[,paste0(Time)],
                                     survT[,paste0(endpoint)] == TRUE) ~ drug)} 
    if(sub<99)  { surv <- coxph(Surv(survT[,paste0(Time)],
                                     survT[,paste0(endpoint)] == TRUE) ~ drug,
                                subset=survT[,paste0(split)]==sub)}    
    
    c(summary(surv)[[7]][,5], summary(surv)[[7]][,2], 
      summary(surv)[[8]][,3], 
      summary(surv)[[8]][,4])
  })
  s <- do.call(rbind, res)
  colnames(s) <- c("p", "HR", "lower", "higher")
  rownames(s) <- drug_names
  s
}


fp <- function( sub, title, d, split, drug_names, a, b, scaleX) {  
  ttt <- com(Time="T5", endpoint="treatedAfter", sub=sub, d=d,
             split=split, drug_names=drug_names, scaleX=scaleX)
  rownames(ttt) <- paste0(rownames(ttt), "_TTT")
  
  os <-  com(Time="T6", endpoint="died", sub=sub, d=d, split=split,
             drug_names=drug_names, scaleX=scaleX)
  rownames(os) <- paste0(rownames(os), "_OS")
  
  n <- c( p=NA, HR=NA, lower=NA, higher=NA )
  nn <- t( data.frame( n ) )
  for (i in 1:(nrow(ttt)-1) ) { nn <-rbind(nn, n )  }
  rownames(nn) <- drug_names
  
  od <- order( c(seq(nrow(nn)), seq(nrow(ttt)), seq(nrow(os)) ))
  
  s <- data.frame( rbind(nn, ttt, os)[od,  ] )
  s$Name <- rownames(s)
  s$x <- 1:nrow(s)
  s$col <- rep(c("white", "black", "darkgreen"), nrow(ttt) )
  s$Endpoint <- factor( c(rep("nn", nrow(nn) ), rep("TTT", nrow(ttt) ),
                          rep("OS", nrow(os) ) )[od] )
  s$features <- "";  s[ which(s$Endpoint=="OS"),"features"] <- drug_names
  s[which(s$Endpoint=="nn"), "Endpoint"] <- "" 
  s <- rbind(s, rep(NA, 8))
  
  p <- ggplot(data=s ,aes(x=x, y=HR, ymin=lower, ymax=higher,
                          colour=Endpoint)) +  geom_pointrange() + 
    theme(legend.position="top", legend.text = element_text(size = 20) ) +
    scale_x_discrete(limits=s$x, labels=s$features ) +
    expand_limits(y=c(a,b)) +
    scale_y_log10(breaks=c(0.01,0.1,0.5,1,2,5,10),
                  labels=c(0.01,0.1,0.5,1,2,5,10)) +
    theme(
      panel.grid.minor = element_blank(),
      axis.title.x = element_text(size=16),
      axis.text.x = element_text(size=16, colour="black"),
      axis.title.y  = element_blank(),
      axis.text.y = element_text(size=12, colour="black"),
      legend.key = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white"),
      legend.title = element_blank(),
      panel.background = element_rect(fill = "white", color="black"),
      panel.grid.major = element_blank(),
      axis.ticks.y = element_blank() 
    ) +
    coord_flip() + 
    scale_color_manual(values=c("OS"="darkgreen", "TTT"="black"),
                       labels=c("OS", "TTT", "")) +
    geom_hline(aes(yintercept=1), colour="black", size=1.5,
               linetype="dashed", alpha=0.3)  +
    annotate("text", x = 1:nrow(s)+0.5, y = s$HR+0.003,
             label = ifelse( s$p<0.001, paste0("p<","0.001"), 
                             paste0("p=", round(s$p,3) ) ), colour=s$col)
  plot(p)
}

## ----Fig5A, fig.path=plotDir, fig.width=5.55, fig.height=( (1+8*1.2) ), dev = c("png", "pdf"), warning=FALSE----
#FIG# S27
d <- c("SF3B1", "NOTCH1", "BRAF", "TP53", "del17p13", "del11q22.3",
       "trisomy12", "IGHV")
drug_names <- c("SF3B1", "NOTCH1", "BRAF", "TP53", "del17p13", "del11q22.3",
                "Trisomy12" ,"IGHV")

fp(sub=99, d=d, drug_names=drug_names, split="IGHV", title="", a=0, b=10,
   scaleX=1)

## ----Fig5B, fig.path=plotDir, fig.width=6.0, fig.height=( (1+8*1.2) ), dev = c("png", "pdf"), warning=FALSE----
#FIG# 6A
d <- c("flu15", "nut15", "dox15", "ibr45", "ide45", "prt45", "selu45",
       "ever45")
drug_names <- c("Fludarabine",  "Nutlin-3", "Doxorubicine", "Ibrutinib",
                "Idelalisib", "PRT062607 HCl", "Selumetinib" ,"Everolimus")

fp(sub=99, d=d, drug_names=drug_names, split="TP53", title="", a=0, b=5,
   scaleX=10)

## ----SFig_genetics1, fig.path=plotDir, fig.width = 8.46, fig.height = 4.5, dev = c("png", "pdf")----
#FIG# S27 left (top+bottom)
par(mfcol=c(1,2))

for (fac in paste(c("IGHV", "TP53"))) {
  survplot( Surv(survT$T5, survT$treatedAfter == TRUE)  ~ as.factor(survT[,fac]), 
            snames=c("wt", "mut"),
            lwd=1.5, cex.axis = 1, cex.lab=1, col= c("darkmagenta", "dodgerblue4"),
            show.nrisk = FALSE,
            legend.pos = FALSE, stitle = "", hr.pos= "topright",
            main = paste(fac), 
            xlab = 'Time (Years)', ylab = 'Time to treatment')
}

## ----SFig_genetics2, fig.path=plotDir, fig.width = 8.46, fig.height = 4.5, dev = c("png", "pdf")----
#FIG# 6B left
#FIG# S27 right (top+bottom)
par(mfcol=c(1,2))

for (fac in paste(c("IGHV", "TP53"))) {
  survplot( Surv(survT$T6, survT$died == TRUE)  ~ as.factor(survT[,fac]), 
            snames=c("wt", "mut"),
            lwd=1.5, cex.axis = 1.0, cex.lab=1.0, col= c("darkmagenta", "dodgerblue4"),
            show.nrisk = FALSE,
            legend.pos = FALSE, stitle = "", hr.pos= "bottomleft",
            main = paste(fac), 
            xlab = 'Time (Years)', ylab = 'Overall survival')
}  

## ----KM, echo=FALSE-----------------------------------------------------------
km <- function(drug, split, title, t, hr, c) { 
  stopifnot(is.character(drug), length(drug)==1, is.character(title),
            length(title)==3)
  
  surv <- survT[!(is.na(survT[,split])), ]
  k <- Biobase::exprs(lpdCLL)[ drug, rownames(surv) ]
  
  ms5 <- maxstat.test(Surv(T5, treatedAfter)  ~ k, 
                      data = surv,
                      smethod = "LogRank",
                      minprop = 0.2, 
                      maxprop = 0.8, 
                      alpha = NULL)
  ms6 <- maxstat.test(Surv(T6, died) ~ k, 
                      data = surv,
                      smethod = "LogRank",
                      minprop = 0.2, 
                      maxprop = 0.8, 
                      alpha = NULL)
  
  
  # median & TTT
  if (c=="med" & t=="TTT") {    
    surv$cutA <- ifelse(
      k >= median( k[which(!(is.na(surv$T5)) ) ] ), "weak", "good")
    surv$cutM <- ifelse(
      k >= median( k[ which( surv[,paste0(split)]==1 & !(is.na(surv$T5)) ) ],
                   na.rm=TRUE ), "weak", "good") 
    surv$cutU <- ifelse(
      k >= median( k[ which( surv[,paste0(split)]==0 & !(is.na(surv$T5)) ) ],
                   na.rm=TRUE ), "weak", "good")
    
  }
  
  # median & OS
  if (c=="med" & t=="OS") {    
    surv$cutA <- ifelse(k >= median(k), "weak", "good")
    surv$cutM <- ifelse(k >= median( k[ which( surv[,paste0(split)]==1 ) ] ),
                        "weak", "good") 
    surv$cutU <- ifelse(k >= median( k[ which( surv[,paste0(split)]==0 ) ] ),
                        "weak", "good")
  }
  
  #TTT & maxstat
  if (c=="maxstat" & t=="TTT") {    
    surv$cutA <- surv$cut5 <- ifelse(k >= ms5$estimate, "weak", "good")
    surv$cutM <- surv$cut5 <- ifelse(k >= ms5$estimate, "weak", "good") 
    surv$cutU <- surv$cut5 <- ifelse(k >= ms5$estimate, "weak", "good")
  }
  
  #OS & maxstat
  if (c=="maxstat" & t=="OS") {    
    surv$cutA <- surv$cut5 <- ifelse(k >= ms6$estimate, "weak", "good")
    surv$cutM <- surv$cut5 <- ifelse(k >= ms6$estimate, "weak", "good") 
    surv$cutU <- surv$cut5 <- ifelse(k >= ms6$estimate, "weak", "good")
  }
  
  drName <- toCaps(drugs[stripConc(drug), "name"])
  
  sp <- function(...)
    survplot(..., 
             lwd = 3, cex.axis = 1.2, cex.lab = 1.5,
             col= c("royalblue", "darkred"), show.nrisk = FALSE,
             legend.pos = FALSE, stitle = "",
             hr.pos=ifelse(hr=="bl", "bottomleft", "topright" ),
             xlab = 'Time (Years)') 
  
  if (t=="TTT") {
    yl <- "Fraction w/o treatment"
    if (c=="med"){
      cat(sprintf("%s median-cutpoint for TTT: %5.2g\n",
                  drName, median(k) ) ) } else 
                  { cat(sprintf("%s cutpoint for TTT: %5.2g\n", drName,
                                ms5$estimate )) }
    
    sp(Surv(surv$T5, surv$treatedAfter) ~ surv$cutA,
       subset = rep(TRUE, nrow(surv)), ylab = yl, main = drName)
    sp(Surv(surv$T5, surv$treatedAfter) ~ surv$cutM,
       subset = surv[, split]==1, ylab = yl,
       main = paste(drName, title[1], title[3])) 
    sp(Surv(surv$T5, surv$treatedAfter) ~ surv$cutU,
       subset = surv[ ,split]==0, ylab = yl,
       main = paste(drName, title[1], title[2])) }
  # OS  
  else {
    yl <- "Fraction overall survival"
    if (c=="med"){
      cat(sprintf("%s median-cutpoint for OS: %5.2g\n",
                  drName, median(k) ) ) } else {
                    cat(sprintf("%s cutpoint for OS: %5.2g\n",
                                drName, ms6$estimate ))}
    sp(Surv(surv$T6, surv$died) ~ surv$cutA,
       subset = rep(TRUE, nrow(surv)), ylab = yl, main = drName)
    sp(Surv(surv$T6, surv$died) ~ surv$cutM,
       subset = surv[, split]==1, ylab = yl,
       main = paste(drName, title[1], title[3])) 
    sp(Surv(surv$T6, surv$died) ~ surv$cutU,
       subset = surv[ ,split]==0, ylab = yl,
       main = paste(drName, title[1], title[2]))
  }
}

## ----KM-TTT-maxstat, fig.path=plotDir, fig.width = 10, fig.height = 3.3, dev = c("png", "pdf")----
par(mfrow=c(1,3), mar=c(5,5,2,0.9))
km(drug = "D_006_1:5", split = "TP53", t="TTT",
   title=c("(TP53", "wt)", "mut)"),  hr="tr", c="maxstat")
km(drug = "D_159_1:5", split = "TP53", t="TTT",
   title=c("(TP53", "wt)", "mut)"), hr="tr",  c="maxstat")
km(drug = "D_010_1:5", split = "TP53", t="TTT",
   title=c("(TP53", "wt)", "mut)"), hr="tr",  c="maxstat")  

km(drug = "D_002_4:5", split = "IGHV", t="TTT",
   title=c("(IGHV",  "wt)" , "mut)"), hr="tr", c="maxstat" )
km(drug = "D_003_4:5", split = "IGHV", t="TTT",
   title=c("(IGHV",  "wt)" , "mut)"), hr="tr", c="maxstat" )
km(drug = "D_166_4:5", split = "IGHV", t="TTT",
   title=c("(IGHV",  "wt)" , "mut)"), hr="tr", c="maxstat" ) 

## ----KM-OS-maxstat, fig.path=plotDir, fig.width = 10, fig.height = 3.3, dev = c("png", "pdf")----
par(mfrow=c(1,3), mar=c(5,5,2,0.9))
km(drug = "D_006_1:5", split = "TP53", t="OS",
   title=c("(TP53", "wt)", "mut)"), hr="bl", c="maxstat")

#FIG# 6B right
#FIG# 6C
km(drug = "D_159_1:5", split = "TP53", t="OS", # doxorubicine
   title=c("(TP53", "wt)", "mut)"), hr="bl", c="maxstat" )

#FIG# 6B middle
km(drug = "D_010_1:5", split = "TP53", t="OS", # nutlin-3
   title=c("(TP53", "wt)", "mut)"), hr="bl", c="maxstat" )

km(drug = "D_002_4:5", split = "IGHV", t="OS",
   title=c("(IGHV",  "wt)" , "mut)"), hr="bl", c="maxstat" )
km(drug = "D_003_4:5", split = "IGHV", t="OS",
   title=c("(IGHV",  "wt)" , "mut)"), hr="bl", c="maxstat" )
km(drug = "D_166_4:5", split = "IGHV", t="OS",
   title=c("(IGHV",  "wt)" , "mut)"), hr="bl", c="maxstat" ) 

## ----extract------------------------------------------------------------------
extractSome <- function(x) {
  sumsu <- summary(x)
  data.frame(
    `p-value`      = 
      sprintf("%6.3g", sumsu[["coefficients"]][, "Pr(>|z|)"]),
    `HR`           = 
      sprintf("%6.3g", signif( sumsu[["coefficients"]][, "exp(coef)"], 2) ), 
    `lower 95% CI` = 
      sprintf("%6.3g", signif( sumsu[["conf.int"]][, "lower .95"], 2) ),
    `upper 95% CI` = 
      sprintf("%6.3g", signif( sumsu[["conf.int"]][, "upper .95"], 2),
              check.names = FALSE) )
}

## ----covariates, echo=FALSE---------------------------------------------------
survT$age <- survT$age/10
survT$IC50beforeTreatment <- ifelse(survT$IC50beforeTreatment==TRUE, 1, 0)
survT$IGHVwt <- ifelse(survT$IGHV==1, 0, 1)

survT$flu15 <- survT$flu15*10
survT$dox15 <- survT$dox15*10

survT$ibr45 <- survT$ibr45*10
survT$ide45 <- survT$ide45*10
survT$prt45 <- survT$prt45*10

## -----------------------------------------------------------------------------
surv1 <- coxph(
  Surv(T6, died) ~  
    age +
    as.factor(IC50beforeTreatment) +
    as.factor(trisomy12) +
    as.factor(del11q22.3) +
    as.factor(del17p13) +
    as.factor(TP53) +
    IGHVwt +
    flu15,       # continuous
  #dox15 +     # continuous
  #flu15:TP53,
  #TP53:dox15,
  data = survT )
extractSome(surv1)

cat(sprintf("%s patients considerd in the model; number of events %1g\n", 
            summary(surv1)$n, summary(surv1)[6] ) )

## ----echo=FALSE, results='hide', eval=FALSE-----------------------------------
#  write(print(xtable(extractSome(surv1))), file=paste0(plotDir,"flu_MD.tex"))

## -----------------------------------------------------------------------------
surv2 <- coxph(
  Surv(T6, died) ~   #as.factor(survT$TP53) , data=survT )
    age +
    as.factor(IC50beforeTreatment) +
    as.factor(trisomy12) +
    as.factor(del11q22.3) +
    as.factor(del17p13) +
    as.factor(TP53) +
    IGHVwt +
    #flu15 +    # continuous
    dox15 ,     # continuous
  #flu15:TP53 ,
  #TP53:dox15,
  data = survT )
extractSome(surv2)

cat(sprintf("%s patients considerd in the model; number of events %1g\n", 
            summary(surv2)$n, summary(surv2)[6] ) )

## ----echo=FALSE, results='hide', eval=FALSE-----------------------------------
#  write(print(xtable(extractSome(surv2))), file=paste0(plotDir,"dox_MD.tex"))

## -----------------------------------------------------------------------------
surv4 <- coxph(
  Surv(T5, treatedAfter) ~ 
    age +
    as.factor(IC50beforeTreatment) +
    as.factor(trisomy12) +
    as.factor(del11q22.3) +
    as.factor(del17p13) +
    IGHVwt +
    ibr45 +
    IGHVwt:ibr45,
  data = survT )

extractSome(surv4)

cat(sprintf("%s patients considerd in the model; number of events %1g\n", 
            summary(surv4)$n, summary(surv4)[6] ) )

## ----echo=FALSE, results='hide', eval=FALSE-----------------------------------
#  write(print(xtable(extractSome(surv4))), file=paste0(plotDir,"ibr_TTT.tex"))

## -----------------------------------------------------------------------------
surv6 <- coxph(
  Surv(T5, treatedAfter) ~ 
    age +
    as.factor(IC50beforeTreatment) +
    as.factor(trisomy12) +
    as.factor(del11q22.3) +
    as.factor(del17p13) +
    IGHVwt +
    ide45 +
    IGHVwt:ide45,
  data = survT )

extractSome(surv6)

cat(sprintf("%s patients considerd in the model; number of events %1g\n",
            summary(surv6)$n, summary(surv6)[6] ) )

## ----echo=FALSE, results='hide', eval=FALSE-----------------------------------
#  write(print(xtable(extractSome(surv6))), file=paste0(plotDir,"ide_TTT.tex"))

## -----------------------------------------------------------------------------
surv8 <- coxph(
  Surv(T5, treatedAfter) ~ 
    age +
    as.factor(IC50beforeTreatment) +
    as.factor(trisomy12) +
    as.factor(del11q22.3) +
    as.factor(del17p13) +
    IGHVwt +
    prt45 +
    IGHVwt:prt45,
  data = survT )

extractSome(surv8)

cat(sprintf("%s patients considerd in the model; number of events %1g\n", 
            summary(surv8)$n, summary(surv8)[6] ) )

## ----echo=FALSE, results='hide', eval=FALSE-----------------------------------
#  write(print(xtable(extractSome(surv8))), file=paste0(plotDir,"prt_TTT.tex"))

## ---- include=!exists(".standalone"), eval=!exists(".standalone")-------------
#  sessionInfo()

## ---- message=FALSE, warning=FALSE, include=FALSE-----------------------------
rm(list=ls())

## -----------------------------------------------------------------------------
sessionInfo()
