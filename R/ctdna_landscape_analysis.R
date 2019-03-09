#!/usr/bin/env Rscript

## run analysis and generate all plots for Zill et al Clinical Cancer Research 2018

require(reshape2)
require(plotrix)
require(grid)
require(gridExtra)
require(GenomeGraphs)
source("R/scatterBarOZ.R")

# Figure 1A
# calculate somatic alteration detection rates for cancer types with >100 patients tested
# clean

detxn <- read.delim("data_tables_for_figures/detxn.csv", header=T, as.is=T, sep=',')
pdf(file="output/GH_detxn_rates_092216.pdf")
par(mar=c(7, 3, 2, 1))
barplot(detxn$pt_detxn, names=detxn$cancer, las=2, ylim=c(0,100), col="cornflowerblue")
dev.off()

# Figure 1C
# boxplots of ctDNA levels per indication
# clean
CTsum = readRDS('data_tables_for_figures/ct_sum.RDS')
df_fig1C = readRDS('data_tables_for_figures/df_fig1C.RDS')
pdf(file="output/sub_ctDNA_level_per_indication_adjVAF.pdf", useDingbats=F)
boxplot(df_fig1C, ylab="ctDNA level (max cfDNA%)", names=rownames(CTsum), pch=19, cex=0.5, las=2, outline=F, log="y")
dev.off()

# # calculate wilcoxon p-value matrix for all indication pairs
# 
# CT_wilcoxon <- matrix(nrow=13,ncol=13,data=NA)
# rownames(CT_wilcoxon) <- rownames(CTsum)
# #colnames(CT_wilcoxon) <- rownames(CTsum)
# rownames(CT_wilcoxon)[1] <- "Colo"
# rownames(CT_wilcoxon)[2] <- "Small Cell Lung"
# rownames(CT_wilcoxon)[4] <- "Hepat"
# rownames(CT_wilcoxon)[13] <- "Glio"
# colnames(CT_wilcoxon) <- rownames(CT_wilcoxon)
# rownames(CT_wilcoxon)[11] <- "Pancreatic Ductal"
# colnames(CT_wilcoxon)[11] <- "Pancreatic Ductal"
# for(j in 1:13){
# 	for(i in 1:length(rownames(CT_wilcoxon))){
# 		CT_wilcoxon[j,i] <- wilcox.test(SCT$max_pct[ grepl(rownames(CT_wilcoxon)[j], SCT$ct) ], SCT$max_pct[ grepl(rownames(CT_wilcoxon)[i], SCT$ct) ])$p.value
# 	}
# }
# 
# options(scipen=5)
# write.csv(CT_wilcoxon, file="output/ctDNA_levels_wilcoxon_pvalues.csv", quote=F)

# Figure 1D, Figure S1E and F
# boxplots of number of somatic SNVs detected vs ctDNA level
# clean

df_fig1D = readRDS('data_tables_for_figures/df_fig1D.RDS')
df_figS1E = readRDS('data_tables_for_figures/df_figS1E.RDS')
df_figS1F = readRDS('data_tables_for_figures/df_figS1F.RDS')

pdf(file="output/sub_ctDNA_level_vs_Num_SSNV_NSCLC.pdf", useDingbats=F)
boxplot(df_fig1D, ylim=c(0,15), pch=19, cex=0.5, names=c("0.1-0.5","0.5-1","1-5","5-10","10-20","20-50"), xlab="ctDNA level per sample", ylab="Somatic SNVs per sample", outline=F)
dev.off()

pdf(file="output/sub_ctDNA_level_vs_Num_SSNV_Breast.pdf", useDingbats=F)
boxplot(df_figS1E, ylim=c(0,15), pch=19, cex=0.5, names=c("0.1-0.5","0.5-1","1-5","5-10","10-20","20-50"), xlab="ctDNA level per sample", ylab="Somatic SNVs per sample", outline=F)
dev.off()

pdf(file="output/sub_ctDNA_level_vs_Num_SSNV_CRC.pdf", useDingbats=F)
boxplot(df_figS1F, ylim=c(0,15), pch=19, cex=0.5, names=c("0.1-0.5","0.5-1","1-5","5-10","10-20","20-50"), xlab="ctDNA level per sample", ylab="Somatic SNVs per sample", outline=F)
dev.off()


# Figure 2A and B

# plot TP53 mirrored frequency plot
load("data_tables_for_figures/tp53_counts.Rdata")
pdf(file="output/TP53_GH_vs_TCGA_freq.pdf")
plot.new()
plot.window(c(1,393), c(-20,20))
rect(0, -1, 393, 1, col="grey", border=NA)
rect(95, -1.5, 289, 1.5, col="firebrick1", border=NA)
rect(318, -1.5, 359, 1.5, col="deepskyblue", border=NA)
rect(5, -1.5, 29, 1.5, col="chartreuse", border=NA)
barplot(200*tp53_counts$cb_freq, width=0.833, col="black", border=NA, yaxt='n', add=T, offset=2)
barplot(-(200*tp53_counts$GH_freq), width=0.833, col="blue", border=NA, yaxt='n', add=T, offset=-2)
axis(2, pos=-2, at=c(2,6,10,14,18), labels=c(0,0.02,0.04,0.06,0.08), col="darkgrey", padj=0.8, cex.axis=0.8)
axis(2, pos=-2, at=-c(2,6,10,14,18), labels=c(0,0.02,0.04,0.06,0.08), col="darkgrey", padj=0.8, cex.axis=0.8)
dev.off()

# calculate correlation of TP53 SNVs per codon in cfDNA vs TCGA
tp53_cor <- cor(tp53_counts$GH_count, tp53_counts$cb_count)
tp53_spmn <- cor(tp53_counts$GH_count, tp53_counts$cb_count, method="spearman")
sprintf("Pearson correlation for TP53 (GH vs TCGA): %.2f", tp53_cor)
sprintf("Spearman correlation for TP53 (GH vs TCGA): %.2f", tp53_spmn)

# plot EGFR mirrored frequency plot
load("data_tables_for_figures/egfr_counts.Rdata")
pdf(file="output/EGFR_GH_vs_TCGA_freq.pdf")
plot.new()
plot.window(c(1,1210), c(-20,20))
rect(0, -1, 1210, 1, col="grey", border=NA)
rect(177, -1.5, 338, 1.5, col="firebrick1", border=NA)
rect(505, -1.5, 637, 1.5, col="deepskyblue", border=NA)
rect(57, -1.5, 168, 1.5, col="chartreuse", border=NA)
rect(361, -1.5, 481, 1.5, col="chartreuse", border=NA)
rect(712, -1.5, 968, 1.5, col="gold", border=NA)
barplot(60*egfr_counts$cb_freq, width=0.833, col="black", border=NA, yaxt='n', add=T, offset=2)
barplot(-(80*egfr_counts$GH_freq), width=0.833, col="blue", border=NA, yaxt='n', add=T, offset=-2)
axis(2, pos=-2, at=c(2,8,14,20), labels=c(0,0.10,0.20,0.30), col="darkgrey", padj=0.8, cex.axis=0.8)
axis(2, pos=-2, at=-c(2,6,10,14,18), labels=c(0,0.05,0.10,0.15,0.20), col="darkgrey", padj=0.8, cex.axis=0.8)
dev.off()

# calculate correlation of EGFR variants per codon in cfDNA vs TCGA
egfr_cor <- cor(egfr_counts$GH_count, egfr_counts$cb_count)
egfr_spmn <- cor(egfr_counts$GH_count, egfr_counts$cb_count, method="spearman")
sprintf("Pearson correlation for EGFR (GH vs TCGA): %.2f", egfr_cor)
sprintf("Spearman correlation for EGFR (GH vs TCGA): %.2f", egfr_spmn)
egfr_kin_cor <- cor(egfr_counts$GH_count[688:875], egfr_counts$cb_count[688:875])
sprintf("Pearson correlation for EGFR kinase domain, exons 18-21 (GH vs TCGA): %.2f", egfr_kin_cor)

# calculate correlation of EGFR variants per codon in cfDNA vs TCGA, with T790M removed
#n790cor1 <- cor(egfr_counts_no_t790m$GH_count, egfr_counts_no_t790m$cb_count)
#n790cor2 <- cor(egfr_counts_no_t790m$GH_count[688:875], egfr_counts_no_t790m$cb_count[688:875])
#sprintf("Pearson correlation for EGFR, no T790M (GH vs TCGA): %.2f", n790cor1)
#sprintf("Pearson correlation for EGFR kinase domain, exons 18-21, no T790M (GH vs TCGA): %.2f", n790cor2)


# Figure 2C - cna_ranks.r

cnv_ranks <- read.csv("data_tables_for_figures/cnv_ranks.csv", header=T, as.is=T)
B_cor <- round(cor(cnv_ranks$Grank, cnv_ranks$Trank, method="spearman"), digits=2)

pdf(file="output/Breast_CNA_ranks_TCGA_vs_GH.pdf")
plot(cnv_ranks$Trank, cnv_ranks$Grank, pch=19, cex=0.5, col="grey", xlim=c(0,20), ylim=c(0,20), xlab="CNA frequency rank (TCGA)", ylab="CNA frequency rank (cfDNA)")
abline(0,1)
points(cnv_ranks$Trank, cnv_ranks$Grank, pch=19, cex=0.5, col="grey", xlim=c(0,20), ylim=c(0,20), xlab="CNA frequency rank (TCGA)", ylab="CNA frequency rank (cfDNA)")
#text(1,19,labels="r=0.87")
text(1,19,labels=paste("r", B_cor, sep="="))
text(cnv_ranks$Trank, cnv_ranks$Grank, labels=cnv_ranks$BG, cex=0.8, pos=3)
dev.off()


# Figure 2D
# clean
fusion <- read.csv(file="data_tables_for_figures/fusion.csv", header=T, as.is=T)
# ALK-EML4 grid plot (table)
AE <- data.frame(EML4_intron=c("13","6","20","other"), ctDNA_bkpts=c("41%","35%","10%","14%"), COSMIC_bkpts=c("47%","35%","14%","4%"), stringsAsFactors=F)
pdf(file="output/ALK_EML4_grid_plot.pdf", height=8)
grid.newpage()
grid.table(AE, rows=NULL)
dev.off()

# ALK-EML4 GenomeGraphs gdPlot
listMarts(host="feb2014.archive.ensembl.org")
mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host="feb2014.archive.ensembl.org")

key <- paste(fusion$gene_a, fusion$gene_b, sep = "|")
i <- which(key == "ALK|EML4")
# plot ALK side
gene <- makeGene(id = "ALK", type="hgnc_symbol", biomart = mart) 
r <- range(fusion$position_A[i])
gr <- new("GeneRegion", start = r[1]-10000, end = r[2]+10000, biomart = mart, strand = "-", chromosome = "2", dp = DisplayPars(color = "black", size = 1))
genomeAxis <- makeGenomeAxis(add53=TRUE, add35=TRUE)
vaf <- makeBaseTrack(fusion$position_A[i], value = fusion$percent_fusion_ab[i], dp = DisplayPars(color = "purple", pch = 19, pointsize = 0.3 ))
pdf(file="output/ALK_fusion_side.pdf")
gdPlot(list(genomeAxis, "VAF" = vaf, "ALK" = gr))
dev.off()

# plot EML4 side
gene <- makeGene(id = "EML4", type="hgnc_symbol", biomart = mart) 
r <- range(fusion$position_B[i])
gr <- new("GeneRegion", start = 42394490, end = 42561688, biomart = mart, strand = "+", chromosome = "2", dp = DisplayPars(color = "black", size = 1))
genomeAxis <- makeGenomeAxis(add53=TRUE, add35=TRUE)
vaf <- makeBaseTrack(fusion$position_B[i], value = fusion$percent_fusion_ab[i], dp = DisplayPars(color = "purple", pch = 19, pointsize = 0.3 ))
pdf(file="output/EML4_fusion_side.pdf")
gdPlot(list(genomeAxis, "VAF" = vaf, "EML4" = gr))
dev.off()


# Figure 3A, Figure S3, Figure S6-9

# cn_vaf_list <- readRDS("data_tables_for_figures/cn_vaf_list.RDS")
# no_cnv <- cn_vaf_list$no_cnv
# p_cnv <- cn_vaf_list$p_cnv
# cnv_genes <- cn_vaf_list$cnv_genes
# 
# no_cnv$copynumber <- 2^(no_cnv$cn)
# 
# no_cnv$color <- "grey"
# no_cnv$color[ no_cnv$alteration == "L858R" ] <- "green"
# no_cnv$color[ no_cnv$gene == "EGFR" & no_cnv$type == "INDEL" ] <- "forestgreen"
# no_cnv$color[ no_cnv$alteration == "T790M" ] <- "firebrick"
# no_cnv$color[ no_cnv$alteration == "C797S" ] <- "orange"
# no_cnv$color[ no_cnv$copynumber == 2 ] <- "white"

# plot copy number vs VAF for EGFR in lung cancer - Figure 3A
# clean
df_fig3A = readRDS('data_tables_for_figures/df_fig3A.RDS')
pdf(file="output/sub_EGFR_CN_VAF_Figure3A.pdf", useDingbats=F)
plot(df_fig3A$percentage, df_fig3A$copynumber, col=df_fig3A$color, pch=19, cex=0.3, ylab="Copy number", xlab="VAF (cfDNA %)", mgp=c(1.5,0.5,0), cex.lab=0.7, cex.axis=0.7, xlim=c(0,100), cex.lab=1.3)
legend("bottomright", legend=c("L858R","e19 del","T790M","C797S"), fill=c("green", "forestgreen", "firebrick", "orange"))
dev.off()

# plot copy number vs VAF for all CNA genes in lung cancer - Figure S6
# clean
df_figS6to9 = readRDS('data_tables_for_figures/df_figS6to9.RDS')
cnv_genes = readRDS('data_tables_for_figures/df_cnv_genes.RDS')
pdf(file="output/sub_CN_vs_VAF_lung.pdf", width=15, height=8, useDingbats=F)
par(mfrow=c(3,6))
par(mar=c(2.5,2.5,2,2))
for(i in 1:length(cnv_genes)){
  plot(df_figS6to9$percentage[ grepl("Lung", df_figS6to9$cancertype) & df_figS6to9$gene == cnv_genes[i] ], 2^(df_figS6to9$cn[ grepl("Lung", df_figS6to9$cancertype) & df_figS6to9$gene == cnv_genes[i] ]), pch=19, cex=0.5, ylab="Copy number", xlab="VAF (cfDNA %)", mgp=c(1.5,0.5,0), cex.lab=0.7, cex.axis=0.7, xlim=c(0,100))
  title(main=cnv_genes[i], line=0.5)
}
dev.off()

# plot copy number vs VAF for all CNA genes in CRC - Figure S7
# clean
df_figS6to9 = readRDS('data_tables_for_figures/df_figS6to9.RDS')
pdf(file="output/sub_CN_vs_VAF_crc.pdf", width=15, height=8, useDingbats=F)
par(mfrow=c(3,6))
par(mar=c(2.5,2.5,2,2))
for(i in 1:length(cnv_genes)){
  plot(df_figS6to9$percentage[ grepl("Colo", df_figS6to9$cancertype) & df_figS6to9$gene == cnv_genes[i] ], 2^(df_figS6to9$cn[ grepl("Colo", df_figS6to9$cancertype) & df_figS6to9$gene == cnv_genes[i] ]), pch=19, cex=0.5, ylab="Copy number", xlab="VAF (cfDNA %)", mgp=c(1.5,0.5,0), cex.lab=0.7, cex.axis=0.7, xlim=c(0,100))
  title(main=cnv_genes[i], line=0.5)
}
dev.off()

# plot copy number vs VAF for all CNA genes in breast cancer - Figure S8
# clean
df_figS6to9 = readRDS('data_tables_for_figures/df_figS6to9.RDS')
pdf(file="output/sub_CN_vs_VAF_breast.pdf", width=15, height=8, useDingbats=F)
par(mfrow=c(3,6))
par(mar=c(2.5,2.5,2,2))
for(i in 1:length(cnv_genes)){
  plot(df_figS6to9$percentage[ grepl("Breast", df_figS6to9$cancertype) & df_figS6to9$gene == cnv_genes[i] ], 2^(df_figS6to9$cn[ grepl("Breast", df_figS6to9$cancertype) & df_figS6to9$gene == cnv_genes[i] ]), pch=19, cex=0.5, ylab="Copy number", xlab="VAF (cfDNA %)", mgp=c(1.5,0.5,0), cex.lab=0.7, cex.axis=0.7, xlim=c(0,100))
  title(main=cnv_genes[i], line=0.5)
}
dev.off()

# plot copy number vs VAF for all CNA genes in prostate cancer - Figure S9
# clean
df_figS6to9 = readRDS('data_tables_for_figures/df_figS6to9.RDS')
pdf(file="output/sub_CN_vs_VAF_prostate.pdf", width=15, height=8, useDingbats=F)
par(mfrow=c(3,6))
par(mar=c(2.5,2.5,2,2))
for(i in 1:length(cnv_genes)){
  plot(df_figS6to9$percentage[ grepl("Prost", df_figS6to9$cancertype) & df_figS6to9$gene == cnv_genes[i] ], 2^(df_figS6to9$cn[ grepl("Prost", df_figS6to9$cancertype) & df_figS6to9$gene == cnv_genes[i] ]), pch=19, cex=0.5, ylab="Copy number", xlab="VAF (cfDNA %)", mgp=c(1.5,0.5,0), cex.lab=0.7, cex.axis=0.7, xlim=c(0,100))
  title(main=cnv_genes[i], line=0.5)
}
dev.off()

# boxplots of copy number for all CNA genes in lung, breast, CRC, prostate - Figure S3
# clean
df_figS3_lung = readRDS('data_tables_for_figures/df_figS3_lung.RDS')
df_figS3_colo = readRDS('data_tables_for_figures/df_figS3_colo.RDS')
df_figS3_breast = readRDS('data_tables_for_figures/df_figS3_breast.RDS')
df_figS3_prost = readRDS('data_tables_for_figures/df_figS3_prost.RDS')

pdf(file="output/sub_CN_per_gene_lung.pdf", useDingbats=F)
bxp(df_figS3_lung,log='y',pch=19,cex=0.5,las=2)
dev.off()
pdf(file="output/sub_CN_per_gene_crc.pdf", useDingbats=F)
bxp(df_figS3_colo,log='y',pch=19,cex=0.5,las=2)
dev.off()
pdf(file="output/sub_CN_per_gene_breast.pdf", useDingbats=F)
bxp(df_figS3_breast,log='y',pch=19,cex=0.5,las=2)
dev.off()
pdf(file="output/sub_CN_per_gene_prostate.pdf", useDingbats=F)
bxp(df_figS3_prost,log='y',pch=19,cex=0.5,las=2)
dev.off()



# Figure 3B & C - pancan_clonal_frxn.r


# histogram for Figure 3B
df_fig3B = readRDS('data_tables_for_figures/df_fig3B.RDS')
pdf(file="output/EGFR_clonality_LUAD_CRC.pdf", useDingbats=F)
hist(df_fig3B$egfr_luad_df_clon2, breaks=100, col=rgb(0,0,1,alpha=0.6), ylim=c(0,250), main=NA, xlab="cfDNA clonality")
hist(df_fig3B$egfr_crc_df_clon2, breaks=100, col=rgb(1,0,0,alpha=0.6), add=T)
abline(v=0.1, lty=2, col="grey")
abline(v=0.9, lty=2, col="grey")
legend(x=0, y=250, legend=c("LUAD","CRC"), fill=c(rgb(0,0,1,alpha=0.6), rgb(1,0,0,alpha=0.6)))
dev.off()


df_fig3C = readRDS('data_tables_for_figures/df_fig3C.RDS')
# APC scatterBar plot for Figure 3C ("rainy day plot")
pdf(file="output/CRC_APC_muts_clonal_fraction.pdf", useDingbats=F)
scatterBarOZ(df_fig3C, xlab="CRC sample number", ylab="cfDNA clonality", pch=19, cex=0.5)
dev.off()


# Figure 3D
# clean
df_fig3D = readRDS('data_tables_for_figures/df_fig3D.RDS')
pdf(file="output/sub_NSCLC_EGFRmut_clonality_NEW.pdf")
boxplot(df_fig3D$clon2[df_fig3D$alteration == "L858R" ], df_fig3D$clon2[df_fig3D$alteration == "T790M"] , df_fig3D$clon2[df_fig3D$alteration == "C797S"], ylab="cfDNA clonality", names=c("L858R","T790M","C797S"), col=c(3,8,4), outline=F)
dev.off()


# Figure 4, Figure S11 & S12
# tables for use with the cBio oncoprinter web app
# just the files for the oncoprinter
# Figure 4 - lung adenocarcinoma oncoprint tables
# 
# read.csv(luad_out, file="output/luad_oncoprint_table_unfilt.csv", quote=F, row.names=F)
# # set the threshold for cfDNA clonality filtering:
# 
# read.csv(luad_filt_out, file="output/luad_oncoprint_table_filtered.csv", quote=F, row.names=F)
# 
# # Figure S11 - breast cancer oncoprint tables
# 
# read.csv(brca_out, file="output/brca_oncoprint_table_unfilt.csv", quote=F, row.names=F)
# # set the threshold for cfDNA clonality filtering:
# 
# read.csv(brca_filt_out, file="output/brca_oncoprint_table_filtered.csv", quote=F, row.names=F)
# 
# # Figure S12 - colorectal cancer oncoprint tables
# 
# read.csv(crc_out, file="output/crc_oncoprint_table_unfilt.csv", quote=F, row.names=F)
# # set the threshold for cfDNA clonality filtering:
# 
# read.csv(crc_filt_out, file="output/crc_oncoprint_table_filtered.csv", quote=F, row.names=F)
# 
# read.csv(crc_long_out, file="output/crc_oncoprint_table_unfilt_LONG.csv", quote=F, row.names=F)
# 
# read.csv(crc_long_filt_out, file="output/crc_oncoprint_table_filtered_LONG.csv", quote=F, row.names=F)


# Figure 5A & B - resistance landscape

resistance_freq <- readRDS(file="data_tables_for_figures/resistance_freq.RDS")

# NSCLC resistance barplot
nsclc_res <- resistance_freq$nsclc_res
nsclc_mat <- matrix(nrow=2,ncol=2,data=c(sum(nsclc_res[1:6,1]), nsclc_res["none",1], sum(nsclc_res[1:6,2]), nsclc_res["none",2]))
sprintf("NSCLC resistance enrichment, p=%s", chisq.test(nsclc_mat)$p.value)
x <- nsclc_res
x[,1] <- 100*(x[,1]/sum(x[,1]))
x[,2] <- 100*(x[,2]/sum(x[,2]))
pdf(file="output/nsclc_resistance_barplot.pdf")
barplot(as.matrix(x), names=colnames(x), col=c(6:1,8), width=0.1, space=2, xlim=c(0,1), ylab="Percentage of samples")
legend("topright", fill=c(6:1,8), legend=rownames(x))
dev.off()

# CRC resistance barplot
crc_res <- resistance_freq$crc_res
crc_mat <- matrix(nrow=2,ncol=2,data=c(sum(crc_res[1:4,1]), crc_res["none",1], sum(crc_res[1:4,2]), crc_res["none",2]))
sprintf("CRC resistance enrichment, p=%s", chisq.test(crc_mat)$p.value)
x <- crc_res
x[,1] <- 100*(x[,1]/sum(x[,1]))
x[,2] <- 100*(x[,2]/sum(x[,2]))
pdf(file="output/crc_res_barplot.pdf")
barplot(as.matrix(x), names=colnames(x), col=c(6:3,8), width=0.1, space=2, xlim=c(0,1), ylab="Percentage of samples")
legend("topright", fill=c(6:3,8), legend=rownames(x))
dev.off()

# Breast cancer resistance barplot
brca_res <- resistance_freq$brca_res
brca_mat <- matrix(nrow=2,ncol=2,data=c(sum(brca_res[1:2,1]), brca_res["none",1], sum(brca_res[1:2,2]), brca_res["none",2]))
sprintf("Breast cancer resistance enrichment, p=%s", chisq.test(brca_mat)$p.value)
x <- brca_res
x[,1] <- 100*(x[,1]/sum(x[,1]))
x[,2] <- 100*(x[,2]/sum(x[,2]))
pdf(file="output/brca_res_barplot.pdf")
barplot(as.matrix(x), names=colnames(x), col=c(6,5,8), width=0.1, space=2, xlim=c(0,1), ylab="Percentage of samples")
legend("topright", fill=c(6,5,8), legend=rownames(x))
dev.off()

# Prostate cancer resistance barplot
prad_res <- resistance_freq$prad_res
prad_mat <- as.matrix(prad_res)
sprintf("Prostate cancer resistance enrichment, p=%s", chisq.test(prad_mat)$p.value)
x <- prad_res
x[,1] <- 100*(x[,1]/sum(x[,1]))
x[,2] <- 100*(x[,2]/sum(x[,2]))
pdf(file="output/prad_res_barplot.pdf")
barplot(as.matrix(x), names=colnames(x), col=c(4,8), width=0.1, space=2, xlim=c(0,1), ylab="Percentage of samples")
legend("topright", fill=c(4,8), legend=rownames(x))
dev.off()

# Figure 5B, Table S13 - table and grid plot of on-label resistance alterations
comb_grid_table <- read.csv("data_tables_for_figures/combined_resistance_table.csv", header=T, as.is=T)
pdf(file="output/combined_grid_plot.pdf", height=36, width=8)
grid.newpage()
grid.table(as.data.frame(comb_grid_table), rows=NULL)
dev.off()


# Figure 5C,D,E - monitoring analysis

mon <- read.csv("data_tables_for_figures/monitoring_data.csv", header=T, as.is=T)

#Figure 5C plot
EL <- mon[ mon$patient == "p1", ]
EL_levels <- levels(as.factor(EL$var_key))
EL_df <- data.frame(vars=EL_levels, colors=1:length(EL_levels), stringsAsFactors = F)
EL_df$colors[7] <- 4
EL_df$colors[5] <- 2
EL_df$colors[2] <- 5
EL_df$colors[4] <- "orange"
pdf(file="output/EGFR_patient_monitoring_Figure5C.pdf")
plot(EL$time_days, EL$percentage, ylim=c(0.1,100), log="y", xlab="Time (days)", ylab="cfDNA VAF (%)", pch=19, cex=0.5, xlim=c(-20,420), col=EL_df$colors[match(EL$var_key, EL_df$vars)])
for(i in 1:nrow(EL_df)){
	lines(x=EL$time_days[ EL$var_key == EL_df$vars[i] ], y=EL$percentage[ EL$var_key == EL_df$vars[i] ], col=EL_df$colors[i])
}
text(x=EL$time_days, y=EL$percentage+0.1, labels=EL$var_key, cex=0.7, col=EL_df$colors[match(EL$var_key, EL_df$vars)])
dev.off()

#Figure 5D plot
PL <- mon[ mon$patient == "p2", ]
PL_levels <- levels(as.factor(PL$var_key))
PL_df <- data.frame(vars=PL_levels, colors=1:length(PL_levels), stringsAsFactors = F)
PL_df$colors[12] <- 3
PL_df$colors[13] <- 2
PL_df$colors[10] <- "orange"
PL_df$colors[11] <- 4
pdf(file="output/PIK3CA_patient_monitoring_Figure5D.pdf")
plot(PL$time_days, PL$percentage, ylim=c(0.1,100), log="y", xlab="Time (days)", ylab="cfDNA VAF (%)", pch=19, cex=0.5, xlim=c(-10,80), col=PL_df$colors[match(PL$var_key, PL_df$vars)])
for(i in 1:nrow(PL_df)){
	lines(x=PL$time_days[ PL$var_key == PL_df$vars[i] ], y=PL$percentage[ PL$var_key == PL_df$vars[i] ], col=PL_df$colors[i])
}
text(x=PL$time_days, y=PL$percentage+0.1, labels=PL$var_key, cex=0.7, col=PL_df$colors[match(PL$var_key, PL_df$vars)])
dev.off()

#Figure 5E plot
TL <- mon[ mon$patient == "p3", ]
TL_levels <- levels(as.factor(TL$var_key))
TL_df <- data.frame(vars=TL_levels, colors=1:length(TL_levels), stringsAsFactors = F)
pdf(file="output/TP53_patient_monitoring_Figure5E.pdf")
plot(TL$time_days, TL$percentage, ylim=c(0.1,100), log="y", xlab="Time (days)", ylab="cfDNA VAF (%)", pch=19, cex=0.5, xlim=c(-10,50), col=TL_df$colors[match(TL$var_key, TL_df$vars)])
for(i in 1:nrow(TL_df)){
	lines(x=TL$time_days[ TL$var_key == TL_df$vars[i] ], y=TL$percentage[ TL$var_key == TL_df$vars[i] ], col=TL_df$colors[i])
}
text(x=TL$time_days, y=TL$percentage+0.1, labels=TL$var_key, cex=0.7, col=TL_df$colors[match(TL$var_key, TL_df$vars)])
dev.off()


# Figure S1A,B,C,D

#Figure S1A
pts_df <- read.csv("data_tables_for_figures/samples_per_patient.csv", header=T, as.is=T)
pdf(file="output/cfDNA_tests_per_patient.pdf")
plot(table(pts_df$number_samples), ylim=c(0,1000), xlab="Number of tests", ylab="Number of patients")
dev.off()
pdf(file="output/cfDNA_tests_per_patient_FigureS1A.pdf")
plot(table(pts_df$number_samples[ pts_df$number_samples > 1 ]), xlab="Number of tests", ylab="Number of patients")
dev.off()

#Figure S1B - cfDNA TMB plot
# clean
df_figS1B = readRDS('data_tables_for_figures/df_figS1B.RDS')
df_figS1B$ct[ grepl("Glio", df_figS1B$ct) ] <- "Glioma/GBM"
canDF <- as.data.frame(table(df_figS1B$ct))
colnames(canDF) <- c("cancer","samples")
canDF <- canDF[ canDF$samples >= 100, ]
rownames(canDF) <- 1:19
canDF <- canDF[ canDF$cancer != "Other", ]
df_figS1B <- df_figS1B[ df_figS1B$ct %in% canDF$cancer, ]
df_figS1B <- df_figS1B[ df_figS1B$num_snv > 0, ]
TMB <- do.call(rbind, lapply(split(df_figS1B$num_snv, df_figS1B$ct), summary))
TMB <- as.data.frame(TMB)
TMB <- TMB[order(TMB$Mean),]
tmb_list <- split(df_figS1B$num_snv, df_figS1B$ct)
tmb_list <- tmb_list[match(rownames(TMB), names(tmb_list))]
tmb_means <- unlist(lapply(tmb_list, function(g) mean(g/0.107)))
pdf(file="output/sub_TMB_cfDNA_common_cancers.pdf", useDingbats=F)
par(mar=c(8, 4, 2, 2))
par(mgp=c(2,0.7,0))
bp <- boxplot(lapply(tmb_list, function(g) g/0.107), cex=0.3, pch=19, las=2, cex.axis=0.5, log="y", ylim=c(8,1000), outline=F, border="grey", ylab="Mutations/Mb")
stripchart(lapply(tmb_list, function(g) g/0.107), method = "jitter", vertical=T, cex=0.3, pch=19, las=2, cex.axis=0.5, log="y", add=T)
points(x=1:18,y=tmb_means,col="red",pch=17)
dev.off()

#Figure S1C
samples_df <- read.csv("data_tables_for_figures/alt_per_sample_hist_S1C.csv", header=T, as.is=T)
pdf(file="output/alterations_per_sample.pdf")
hist(as.numeric(samples_df$num_alts), breaks=0:167, xlim=c(1,50), col="grey", xlab="Number of alterations per sample", main=NA)
dev.off()

#Figure S1D
high_ct_df <- read.csv("data_tables_for_figures/alts_in_high_ctdna_samples_S1D.csv", header=T, as.is=T)
pdf(file="output/alterations_per_sample_high_ctDNA_overlay.pdf")
hist(as.numeric(samples_df$num_snv), breaks=0:167, xlim=c(1,30), col="grey", xlab="Number of SNVs per sample", main=NA)
hist(as.numeric(high_ct_df$clon_snv), breaks=0:41, xlim=c(1,50), col=rgb(1,0,0,alpha=0.5), add=T)
hist(as.numeric(high_ct_df$num_snv), breaks=0:131, xlim=c(1,50), col=rgb(0,0,1,alpha=0.5), add=T)
legend("topright", legend=c("All samples","High-ctDNA samples","High-ctDNA samples, Clonality filtered"), fill=c("grey", rgb(0,0,1,alpha=0.5), rgb(1,0,0,alpha=0.5)))
dev.off()


# Figure S2

# calculate proportions of clonal variants in samples grouped by ctDNA level (max_VAF):
# ( Number of NSCLC variants in samples with max_VAF <= 0.5% & with variant clonality > 0.9 ) / ( Number of NSCLC variants in samples with max_VAF <= 0.5% ) = 0.5236098
# ( Number of NSCLC variants in samples with 0.5% < max_VAF <= 1% & with variant clonality > 0.9 ) / ( Number of NSCLC variants in samples with 0.5% < max_VAF <= 1% ) = 0.3418711
# ( Number of NSCLC variants in samples with 1% < max_VAF <= 5% & with variant clonality > 0.9 ) / ( Number of NSCLC variants in samples with 1% < max_VAF <= 5% ) = 0.292123
# ( Number of NSCLC variants in samples with 5% < max_VAF <= 10% & with variant clonality > 0.9 ) / ( Number of NSCLC variants in samples with 5% < max_VAF <= 10% ) = 0.26391
# ( Number of NSCLC variants in samples with 10% < max_VAF <= 20% & with variant clonality > 0.9 ) / ( Number of NSCLC variants in samples with 10% < max_VAF <= 20% ) = 0.2561092
# ( Number of NSCLC variants in samples with 20% < max_VAF <= 50% & with variant clonality > 0.9 ) / ( Number of NSCLC variants in samples with 20% < max_VAF <= 50% ) = 0.2344256

# 95% CI:
# 0.5 - 0.5103, 0.5369
# 1 - 0.3267, 0.3575
# 5 - 0.2825, 0.3019
# 10 - 0.249, 0.2794
# 20 - 0.241, 0.2718
# 50 - 0.2188, 0.2507

x=c(0.5,1,5,10,20,50)  # upper ends of max_VAF bins
y=c(0.524,0.342,0.292,0.264,0.256,0.234)  # proportions of clonal variants
m <- nls(y~a*x/(b+x))
#cor(y,predict(m))

pdf(file="output/Max_VAF_vs_clonal_proportion.pdf")
plotCI(x=c(0.5,1,5,10,20,50), y=c(0.524,0.342,0.292,0.264,0.256,0.234), ui=c(0.537,0.358,0.302,0.279,0.272,0.251), li=c(0.510,0.327,0.283,0.249,0.241,0.219), cex=0.5, pch=19, xlab="Max VAF per sample", ylab="Proportion of clonal variants")
lines(x,predict(m),lty=2,col="grey")
dev.off()


# Figure S5, Table S6
comb_prev <- read.delim("data_tables_for_figures/driver_alt_prevalence_TableS6.tsv", header=T, as.is=T)
prev_cor <- cor.test(comb_prev$TCGA_freq, comb_prev$cfDNA_freq)
pdf(file="output/driver_prevalence_cfDNA_vs_TCGA.pdf", useDingbats=F)
plot(comb_prev$TCGA_freq, comb_prev$cfDNA_freq, cex=0.5, pch=19, col=comb_prev$color, xlim=c(0,0.6), ylim=c(0,0.6), xlab="TCGA prevalence", ylab="cfDNA prevalence")
abline(0, 1, col="grey")
text(0.58,0.01,labels=paste("r=", round(prev_cor$estimate, digits=2), sep=""))
legend("topleft", legend=c("LUAD","BRCA","CRC"), fill=c("green","red","blue"))
text(comb_prev$TCGA_freq[c(1:3,15,16,19,26:28)], comb_prev$cfDNA_freq[c(1:3,15,16,19,26:28)], labels=comb_prev$gene[c(1:3,15,16,19,26:28)], col=comb_prev$color[c(1:3,15,16,19,26:28)], cex=0.8, pos=3)
dev.off()


# Figure S10
# CODE DOESN'T WORK FOR THESE FIGURES DUE TO MISSING COLUMN p_cnv$max_pct 
# #fit linear model: log2(CN) ~ max_pct
# #entire cohort
# fit <- lm( log2(p_cnv$copynumber) ~ p_cnv$max_pct )
# pdf(file="output/CN_vs_maxVAF_lmfit_all.pdf", useDingbats=F)
# plot(p_cnv$max_pct, log2(p_cnv$copynumber), pch=19, cex=0.5, xlab="Max VAF (cfDNA %)", ylab="log2(copy number)", main="All CNV-positive samples")
# abline(fit, col="red")
# dev.off()
# 
# #NSCLC cases only
# nsclc <- unique(p_cnv$sampleid[ p_cnv$cancertype == "Non-small Cell Lung Carcinoma (NSCLC)" | p_cnv$cancertype == "Lung Adenocarcinoma" | p_cnv$cancertype == "Lung Squamous Cell Carcinoma" ])
# fit_nsclc <- lm( log2(p_cnv$copynumber[ p_cnv$sampleid %in% nsclc ]) ~ p_cnv$max_pct[ p_cnv$sampleid %in% nsclc ] )
# pdf(file="output/CN_vs_maxVAF_lmfit_NSCLC.pdf", useDingbats=F)
# plot(p_cnv$max_pct[ p_cnv$sampleid %in% nsclc ], log2(p_cnv$copynumber[ p_cnv$sampleid %in% nsclc ]), pch=19, cex=0.5, xlab="Max VAF (cfDNA %)", ylab="log2(copy number)", main="CNV-positive NSCLC samples")
# abline(fit_nsclc, col="red")
# dev.off()
# 
# #before and after CNA driver determination
# pdf(file="output/CN_vs_maxVAF_CNV_drivers_nondrivers.pdf", useDingbats=F)
# par(mfrow=c(1,3))
# plot(p_cnv$max_vaf2, p_cnv$copynumber, pch=19, cex=0.5, main="All CNVs", xlab="Max VAF per sample", ylab="Copy number", ylim=c(0,50))
# abline(h=2.2, col="red")
# plot(p_cnv$max_vaf2[ grepl("driver", p_cnv$candidate) ], p_cnv$copynumber[ grepl("driver", p_cnv$candidate) ], pch=19, cex=0.5, main="Driver CNV", xlab="Max VAF per sample", ylab="Copy number", ylim=c(0,50))
# abline(h=2.2, col="red")
# plot(p_cnv$max_vaf2[ p_cnv$candidate == "not_driv" ], p_cnv$copynumber[ p_cnv$candidate == "not_driv" ], pch=19, cex=0.5, main="Non-driver CNV", xlab="Max VAF per sample", ylab="Copy number", ylim=c(0,50))
# abline(h=2.2, col="red")
# dev.off()


# Figure S13

kras_plot <- read.csv("data_tables_for_figures/kras_plot.csv", header=T, as.is=T)
par(mgp=c(2,0.7,0))
pdf(file="output/CRC_KRASmut_clonality.pdf", useDingbats=F)
plot(kras_plot$max_pct[ grepl("Colo", kras_plot$cancertype) & kras_plot$type != "CNV" ], kras_plot$percentage[ grepl("Colo", kras_plot$cancertype) & kras_plot$type != "CNV" ], cex=0.5, pch=19, xlim=c(0,100), ylim=c(0,100), xlab="Maximum VAF in sample", ylab="KRASmut VAF", col=kras_plot$color[ grepl("Colo", kras_plot$cancertype) & kras_plot$type != "CNV" ], main="Colorectal")
abline(a=0,b=1,lty=2)
dev.off()

pdf(file="output/NSCLC_KRASmut_clonality.pdf", useDingbats=F)
plot(kras_plot$max_pct[ kras_plot$cancertype == "NSCLC" & kras_plot$type != "CNV" ], kras_plot$percentage[ kras_plot$cancertype == "NSCLC" & kras_plot$type != "CNV" ], cex=0.5, pch=19, xlim=c(0,100), ylim=c(0,100), xlab="Maximum VAF in sample", ylab="KRASmut VAF", col=kras_plot$color[ kras_plot$cancertype == "NSCLC" & kras_plot$type != "CNV" ], main="NSCLC")
abline(a=0,b=1,lty=2)
dev.off()

pdf(file="output/breast_KRASmut_clonality.pdf", useDingbats=F)
plot(kras_plot$max_pct[ grepl("Breast", kras_plot$cancertype) & kras_plot$type != "CNV" ], kras_plot$percentage[ grepl("Breast", kras_plot$cancertype) & kras_plot$type != "CNV" ], cex=0.5, pch=19, xlim=c(0,100), ylim=c(0,100), xlab="Maximum VAF in sample", ylab="KRASmut VAF", col=kras_plot$color[ grepl("Breast", kras_plot$cancertype) & kras_plot$type != "CNV" ], main="Breast")
abline(a=0,b=1,lty=2)
dev.off()

pdf(file="output/prostate_KRASmut_clonality.pdf", useDingbats=F)
plot(kras_plot$max_pct[ grepl("Prostate", kras_plot$cancertype) & kras_plot$type != "CNV" ], kras_plot$percentage[ grepl("Prostate", kras_plot$cancertype) & kras_plot$type != "CNV" ], cex=0.5, pch=19, xlim=c(0,100), ylim=c(0,100), xlab="Maximum VAF in sample", ylab="KRASmut VAF", col=kras_plot$color[ grepl("Prostate", kras_plot$cancertype) & kras_plot$type != "CNV" ], main="Prostate")
abline(a=0,b=1,lty=2)
dev.off()

pdf(file="output/pancreas_KRASmut_clonality.pdf", useDingbats=F)
plot(kras_plot$max_pct[ kras_plot$cancertype == "Pancreatic Ductal Adenocarcinoma" & kras_plot$type != "CNV" ], kras_plot$percentage[ kras_plot$cancertype == "Pancreatic Ductal Adenocarcinoma" & kras_plot$type != "CNV" ], cex=0.5, pch=19, xlim=c(0,100), ylim=c(0,100), xlab="Maximum VAF in sample", ylab="KRASmut VAF", col=kras_plot$color[ kras_plot$cancertype == "Pancreatic Ductal Adenocarcinoma" & kras_plot$type != "CNV" ], main="Pancreas")
abline(a=0,b=1,lty=2)
dev.off()

pdf(file="output/ovarian_KRASmut_clonality.pdf", useDingbats=F)
plot(kras_plot$max_pct[ kras_plot$cancertype == "Ovarian Carcinoma" & kras_plot$type != "CNV" ], kras_plot$percentage[ kras_plot$cancertype == "Ovarian Carcinoma" & kras_plot$type != "CNV" ], cex=0.5, pch=19, xlim=c(0,100), ylim=c(0,100), xlab="Maximum VAF in sample", ylab="KRASmut VAF", col=kras_plot$color[ kras_plot$cancertype == "Ovarian Carcinoma" & kras_plot$type != "CNV" ], main="Ovarian")
abline(a=0,b=1,lty=2)
dev.off()


# Tables S7-12 - mutual exclusivity stats

genes_muts <- readRDS("data_tables_for_figures/genes_muts.RDS")
cancer_types <- c("brca","brca_filt","crc","crc_filt","luad","luad_filt","crc_long","crc_long_filt")
pt_total <- c(genes_muts$brca_pt_total, genes_muts$brca_pt_total_filt, genes_muts$crc_pt_total, genes_muts$crc_pt_total_filt, genes_muts$luad_pt_total, genes_muts$luad_pt_total_filt, genes_muts$crc_pt_total, genes_muts$crc_pt_total_filt)

for(j in cancer_types){
  # prepare data for each cancer type
  p <- genes_muts[[j]]
  colnames(p) <- c("sampleid","gene","alteration","type")
  jj <- gsub("_filt","",j)
  genes <- paste(jj, "genes", sep="_")
  genes_list <- genes_muts[[genes]]

  # check that all genes in alterations table

  stopifnot(genes_list %in% p$gene)

  # build data frame of all possible pairs of genes from the above list

  temp <- NULL
  for(i in 1:length(genes_list)){
    temp <- c(temp, rep(genes_list[i], (length(genes_list) - i)))
  }

  temp2 <- NULL
  for(i in 2:length(genes_list)){
    temp2 <- c(temp2, genes_list[i:length(genes_list)])
  }

  d <- data.frame(gene1=temp, gene2=temp2, stringsAsFactors=F)
  d$pvalue <- NA
  d$logOR <- NA
  d$`++` <- NA
  d$`+-` <- NA
  d$`-+` <- NA
  d$`--` <- NA

  options(scipen=3)
  options(digits=3)

  # build 2x2 tables for each row in genes data frame, then run fisher.test

  for(i in 1:nrow(d)){
    g1 <- unique(p$sampleid[ p$gene == d$gene1[i] ])
    g2 <- unique(p$sampleid[ p$gene == d$gene2[i] ])
    A <- length(g1[ g1 %in% g2 ])
    B <- length(g1[ !(g1 %in% g2) ])
    C <- length(g2[ !(g2 %in% g1) ])
    D <- pt_total[ which(cancer_types == j) ] - (A+B+C)
    mat <- matrix(c(A,B,C,D), 2, 2)
    f <- fisher.test(mat)
    d$pvalue[i] <- f$p.value
    d$logOR[i] <- as.numeric(log(f$estimate))
    d$`++`[i] <- A
    d$`+-`[i] <- B
    d$`-+`[i] <- C
    d$`--`[i] <- D
  }

  # adjust p-value for multiple comparisons per cancer type (Benjamini-Hochberg)
  d$padj <- p.adjust(d$pvalue, method="BH")

  # write out each results table to CSV
  outfile <- paste(j, "mutex", "csv", sep=".")
  write.csv(d, file=file.path("output", outfile), quote=F, row.names=F)

}


# Table S14
# clean
df_jak2 = readRDS('data_tables_for_figures/df_jak2')
jak2_table <- as.data.frame(sort(table(df_jak2$alteration), decreasing = T))
colnames(jak2_table) <- c("SNV","Frequency")
write.table(jak2_table, file="output/sub_JAK2_SNVs_table.txt", sep="\t", quote=F, row.names=F)
