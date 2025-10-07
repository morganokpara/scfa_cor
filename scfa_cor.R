#install.packages("readxl")
#install.packages("vegan")

rm(list=ls())
genus <- readxl::read_excel("Morgan_Results.xlsx")
data <- data.frame(genus[,-1],row.names = genus$GENUS)
scfa <- readxl::read_excel("Morgan SCFA results_Excel.xlsx")
scfa <- data.frame(scfa[,-1],row.names = scfa$SAMPLE)

data.scfa <- data[,match(rownames(scfa),colnames(data))]

colnames(scfa)

acetate.cor.r <- c()
acetate.cor.p <- c()
propionate.cor.r <- c()
propionate.cor.p <- c()
butyrate.cor.r <- c()
butyrate.cor.p <- c()
valerate.cor.r <- c()
valerate.cor.p <- c()
caproate.cor.r <- c()
caproate.cor.p <- c()
isobutyrate.cor.r <- c()
isobutyrate.cor.p <- c()
isovalerate.cor.r <- c()
isovalerate.cor.p <- c()

for(i in 1:nrow(data.scfa)){
  scfa.1 <- cor.test(t(data.scfa[i,]),scfa$ACETIC.ACID..mM.100g.)
  scfa.2 <- cor.test(t(data.scfa[i,]),scfa$PROPIONIC.ACID.mM.100g.)
  scfa.3 <- cor.test(t(data.scfa[i,]),scfa$BUTYRIC.ACID.mM.100g.)
  scfa.4 <- cor.test(t(data.scfa[i,]),scfa$VALERIC.ACID.mM.100g.)
  scfa.5 <- cor.test(t(data.scfa[i,]),scfa$CAPROIC.ACID.mM.100g.)
  scfa.6 <- cor.test(t(data.scfa[i,]),scfa$ISOBUTYRIC.ACID.mM.100g.)
  scfa.7 <- cor.test(t(data.scfa[i,]),scfa$ISOVALERIC.ACID.mM.100g.)
  
  acetate.cor.r[i] <- scfa.1$estimate
  acetate.cor.p[i] <- scfa.1$p.value
  propionate.cor.r[i] <- scfa.2$estimate
  propionate.cor.p[i] <- scfa.2$p.value
  butyrate.cor.r[i] <- scfa.3$estimate
  butyrate.cor.p[i] <- scfa.3$p.value
  valerate.cor.r[i] <- scfa.4$estimate
  valerate.cor.p[i] <- scfa.4$p.value
  caproate.cor.r[i] <- scfa.5$estimate
  caproate.cor.p[i] <- scfa.5$p.value
  isobutyrate.cor.r[i] <- scfa.6$estimate
  isobutyrate.cor.p[i] <- scfa.6$p.value
  isovalerate.cor.r[i] <- scfa.7$estimate
  isovalerate.cor.p[i] <- scfa.7$p.value
}


acetate.cor.a <- p.adjust(acetate.cor.p,"fdr")
propionate.cor.a <- p.adjust(propionate.cor.p,"fdr")
butyrate.cor.a <- p.adjust(butyrate.cor.p,"fdr")
valerate.cor.a <- p.adjust(valerate.cor.p,"fdr")
caproate.cor.a <- p.adjust(caproate.cor.p,"fdr")
isobutyrate.cor.a <- p.adjust(isobutyrate.cor.p,"fdr")
isovalerate.cor.a <- p.adjust(isovalerate.cor.p,"fdr")


sig.idx <- which(acetate.cor.p < 0.05 | propionate.cor.p < 0.05 | butyrate.cor.p <0.05 | 
                   valerate.cor.p < 0.05 | caproate.cor.p < 0.05 | isobutyrate.cor.p < 0.05 | 
                   isovalerate.cor.p < 0.05)
scfa.cor <- data.frame(cbind(
  acetate.cor.r,
  propionate.cor.r,
  butyrate.cor.r,
  valerate.cor.r,
  caproate.cor.r,
  isobutyrate.cor.r,
  isovalerate.cor.r),row.names = rownames(data.scfa))
scfa.pval <- data.frame(cbind(
  acetate.cor.p,
  propionate.cor.p,
  butyrate.cor.p,
  valerate.cor.p,
  caproate.cor.p,
  isobutyrate.cor.p,
  isovalerate.cor.p),row.names = rownames(data.scfa))
scfa.pval[scfa.pval < 0.001] <- "***"
scfa.pval[scfa.pval < 0.01 & scfa.pval >= 0.001] <- "**"
scfa.pval[scfa.pval < 0.05 & scfa.pval >= 0.01] <- "*"
scfa.pval[scfa.pval >= 0.05] <- ""

pheatmap::pheatmap(scfa.cor[sig.idx,],cellwidth = 16,cellheight = 10,
                   cluster_rows = F,cluster_cols = F,
                   display_numbers = as.matrix(scfa.pval[sig.idx,]),number_color = "#FFFFFF",
                   border_color = NA)

library(vegan)
bc.dist <- vegdist(t(data.scfa))
pcoa <- cmdscale(bc.dist,eig = T)
groups <- factor(c("P","P","P","D","D","D","L","L","L"))
permanova <- adonis2(bc.dist ~ groups)
mod <- envfit(pcoa,scfa)

pct <- round(prop.table(pcoa$eig[pcoa$eig > 0]) * 100, 2)

plot(scores(pcoa), type = "n", las = 1,
     xlab = paste0("PCoA1 (",round(pct[1], 2), "%)"),
     ylab = paste0("PCoA2 (",round(pct[2], 2), "%)"))
abline(h=0,v=0,lty=3)
points(scores(pcoa), pch = 19, 
       cex = 1, col = as.numeric(groups))
plot(mod,p.max = 0.1)
legend("topright",legend = levels(groups),col = 1:3,pch=19,bty="n")
