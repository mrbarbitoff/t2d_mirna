setwd("/media/barbitoff/DATA/Working issues/WES/microRNA")
library(caret)
library(ggplot2)
library(reshape2)
library(effsize)
library(DESeq2)
library(rgl)
options(rgl.printRglwidget = TRUE)
library(pROC)
library(VennDiagram)
library(RColorBrewer)
library(ggVennDiagram)
library(ggrepel)
library(cowplot)
library(clusterProfiler)
library(msigdbr)
library(limma)
library(org.Hs.eg.db)
myCol <- brewer.pal(3, "Pastel2")

### Basic data loading and preparation
all_mir = read.table('miRNA_all.tsv', sep='\t', header=T)
all_mir = all_mir[all_mir$miRNA != 'piRNA', ]

ggplot(na.omit(all_mir), aes(x=Z__K3_S55.UMIs, y=Z__K3_S55.READs)) + geom_point()
all_mir[all_mir$Z__K3_S55.UMIs > 500 & all_mir$Z__K3_S55.READs < 100, ]
colnames(all_mir)
all_mir = all_mir[, 1:67]
#all_mir = all_mir[, c(1, 68:133)]
write.table(colnames(all_mir), file='sample_ids.txt', quote=F, row.names=F)

design = read.table('design_umi.tsv', header=T, sep='\t')
head(design)
design$Diabetes = as.factor(ifelse(design$Diabetes == 0, 'control', 't2d'))
design$Obesity = as.factor(ifelse(design$Obesity == 0, 'normal', 'obese'))
design$FamHistory = as.factor(ifelse(design$FamHistory == 0, 'no', 'yes'))
rownames(design) = design$Sample
head(design)

#### DE analtsis with different contrasts

mir_data <- round((as.matrix(all_mir[, 2:67])), digits=0)
rownames(mir_data) <- all_mir$miRNA

dds <- DESeqDataSetFromMatrix(countData = na.omit(mir_data),
                              colData = design,
                              design = ~ Diabetes)
dds <- DESeq(dds)
res_full  <- results(dds, contrast=c("Diabetes","control","t2d"))
res_df <- na.omit(as.data.frame(res_full))
head(res_df)
table(res_df$padj < 0.05)
res_df$regulation = ifelse(res_df$padj < 0.05 & res_df$log2FoldChange >= 1, 'Up', 
          ifelse(res_df$padj < 0.05 & res_df$log2FoldChange <= -1, 'Down', 'none'))
res_df$trait = 'T22D'
table(res_df$regulation)

# Checking known marker miRNA
known_mirna <- c('hsa-miR-126', 'hsa-miR-375', 'hsa-miR-34a', 
                 'hsa-miR-29a', 'hsa-miR-144')
res_df[as.logical(rowSums(sapply(known_mirna, function(x) 
                  grepl(x, rownames(res_df))))), ]

# By obesity
dds_ob <- DESeqDataSetFromMatrix(countData = na.omit(mir_data),
                              colData = design,
                              design = ~ Obesity)
dds_ob <- DESeq(dds_ob)
res_full_ob  <- results(dds_ob, contrast=c("Obesity","normal","obese"))
res_df_ob <- na.omit(as.data.frame(res_full_ob))
head(res_df_ob)
table(res_df_ob$padj < 0.05)
res_df_ob$regulation = ifelse(res_df_ob$padj < 0.05 & res_df_ob$log2FoldChange >= 1, 
  'Up', ifelse(res_df_ob$padj < 0.05 & res_df_ob$log2FoldChange <= -1, 'Down', 'none'))
res_df_ob$trait = 'Obesity'


# By famhistory
dds_fh <- DESeqDataSetFromMatrix(countData = na.omit(mir_data),
                                 colData = design,
                                 design = ~ FamHistory)
dds_fh <- DESeq(dds_fh)
res_full_fh  <- results(dds_fh, contrast=c("FamHistory","no","yes"))
res_df_fh <- na.omit(as.data.frame(res_full_fh))
head(res_df_fh)
table(res_df_fh$padj < 0.05)
res_df_fh[res_df_fh$padj < 0.05, ]

res_df_fh$regulation = ifelse(res_df_fh$padj < 0.05 & res_df_fh$log2FoldChange >= 1, 
    'Up', ifelse(res_df_fh$padj < 0.05 & res_df_fh$log2FoldChange <= -1, 'Down', 'none'))
res_df_fh$trait = 'Family history'

all_res = rbind(res_df, res_df_ob, res_df_fh)
all_volc <- ggplot(all_res, aes(x=log2FoldChange, 
                                y=-log10(pvalue), 
                                col=regulation)) + 
  geom_point() + theme_bw() + 
  scale_color_manual(values=c('#1c61b6AA', 'gray', '#ff766aff')) +
  facet_wrap(~trait, nrow=1)
all_volc

label_df <- res_df[c('hsa-miR-5588-5p', 'hsa-miR-125b-2-3p', 
                     'hsa-miR-496', 'hsa-miR-1284'), ]
label_df$miRNA = rownames(label_df)
all_volc_wlab = all_volc + geom_text_repel(data=label_df, 
                           aes(label=miRNA))
all_volc_wlab

# PCA
pcaout <- prcomp(t(assay(rlog(dds))))
summary(pcaout)

pcax = as.data.frame(pcaout$x)
head(pcax)
pcax$T2D = sapply(rownames(pcax), function(x) design[x, 'Diabetes'])
pca_t2d <- ggplot(pcax, aes(x=PC1, y=PC2, col=T2D)) + geom_point() + 
  theme_bw() + scale_color_brewer(palette="Set1") + 
  xlab('PC1 - 26.4% of variance') + ylab('PC2 - 9.9% of variance')
pca_t2d
pca_t2d + stat_ellipse()

pcax$coloring = ifelse(as.character(pcax$T2D) == 't2d', 
                          'red', 'blue')
plot3d(x=pcax$PC1, y=pcax$PC2, z=pcax$PC3, 
       col=pcax$coloring, type='s', xlab='PC1', ylab='PC2', zlab='PC3')


### Venn diagram of significant miRNAs
signif_t2d <- rownames(res_df[res_df$padj < 0.05, ])
signif_ob <- rownames(res_df_ob[res_df_ob$padj < 0.05, ])
signif_fh <- rownames(res_df_fh[res_df_fh$padj < 0.05, ])
significant <- list(signif_t2d, signif_ob, signif_fh)
names(significant) = c('T2D', 'Obesity', 'Family history')

venn.diagram(
  x = list(signif_t2d, signif_ob, signif_fh),
  category.names = c("T2D" , "Obesity" , "Family History"),
  lwd = 2,  fill = myCol,
  filename = 'deg_venn_diagramm.pdf',
  output=TRUE
)

vd <- ggVennDiagram(significant, edge_color='black', lwd = 0.8, lty = 1) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")
vd

f1 <- plot_grid(plot_grid(pca_t2d, vd, labels=c('a', 'c')), all_volc, 
                labels=c('', 'b'), nrow=2)
f1

f1_x <- plot_grid(plot_grid(pca_t2d, vd, labels=c('a', 'c')), all_volc_wlab, 
          labels=c('', 'b'), nrow=2)
cairo_pdf(filename='Figure_1.pdf', width=10, height=6) 
print(f1_x)
dev.off()


# Candidates
selected_mirna = c('hsa-miR-5588-5p', 'hsa-miR-125b-2-3p', 
                   'hsa-miR-496', 'hsa-miR-1284')
sapply(selected_mirna, function(x) x %in% signif_t2d)
sapply(selected_mirna, function(x) x %in% signif_ob)
sapply(selected_mirna, function(x) x %in% signif_fh)

res_df[selected_mirna, ]


#### Playing around with a PCA plot of selected miRs
select_counts <- all_mir[all_mir$miRNA %in% selected_mirna, ]
pcaout <- prcomp(t(as.matrix(select_counts[, 2:67])))
pca_res = as.data.frame(pcaout$x)
pca_res$diabetes = design$Diabetes
ggplot(pca_res, aes(x=PC1, y=PC2, col=diabetes)) + geom_point(size=2) +
  theme_bw()
pca_res$coloring = ifelse(as.character(pca_res$diabetes) == 't2d', 
                          'red', 'blue')
plot3d(x=pca_res$PC1, y=pca_res$PC2, z=pca_res$PC3, 
       col=pca_res$coloring, type='s', xlab='PC1', ylab='PC2', zlab='PC3')



### Validating the relevance of the four selected miRs with simple ML
dumb_model <- data.frame(sample=design$Sample, diabetes=design$Diabetes,
            var1 = colSums(select_counts[c(1, 2, 4), 2:67]),
            var2 = as.numeric(select_counts[3, 2:67]))
head(dumb_model)
ggplot(dumb_model, aes(x=var1, y=var2, col=diabetes)) + geom_point()

#### Examining prediction of T2D status on the selected 4 miRs
for_ml <- as.data.frame(t(all_mir))
colnames(for_ml) = for_ml[1, ]
for_ml <- for_ml[-1, ]
for_ml <- apply(for_ml, 2, as.numeric)
for_ml <- as.data.frame(cbind(design, for_ml))
head(for_ml[, 1:10])

set.seed(32)
ml_select <- for_ml[c(1:20, sample(21:66, 20)), c('Diabetes', 'Obesity', 'FamHistory', 
           'hsa-miR-5588-5p', 'hsa-miR-125b-2-3p',
           'hsa-miR-496', 'hsa-miR-1284'), ]
str(ml_select)
colnames(ml_select) = c('Diabetes', 'Obesity', 'FamHistory',
                        'miR_5588_5p', 'miR_125b_2_3p', 'miR_496', 'miR_1284')
head(ml_select)
train_control <- trainControl(method="cv", number=3, classProbs=T,
          summaryFunction=twoClassSummary, savePredictions = T)

model_dumb = train(Diabetes ~ miR_5588_5p + miR_125b_2_3p + miR_496 + miR_1284,
                  data=ml_select, 
                  trControl=train_control, 
                  method="svmLinear", na.action = 'na.omit')
print(model_dumb)
print(model_dumb$results)

saveRDS(model_dumb, file="microRNA_SVM.RData")

rocobj <- plot.roc(model_dumb$pred$obs, print.auc = TRUE,
         model_dumb$pred$t2d, percent=TRUE, ci=TRUE)
ciobj <- ci.se(rocobj)
plot(ciobj, type = "shape", col = "#1c61b6AA") 
plot(ci(rocobj, of = "thresholds", thresholds = "best")) 


### Gene set enrichment analysis for selected and validated miRs
m_t2g <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  dplyr::select(gs_name, entrez_gene)
#m_t2g <- m_t2g[sapply(m_t2g$gs_name, function(x) 
#            grepl('KEGG|REACTOME|WP', x, perl=T)), ]
head(m_t2g)

# Reading targets and converting to Entrez
target_list <- read.table('./miRTar_all_targets.txt', header=F)$V1
ofs_targets = alias2Symbol(target_list)
target_eids = mapIds(org.Hs.eg.db, ofs_targets, 
                     column = "ENTREZID", keytype = "SYMBOL")
write.table(target_eids, 'entrez_targets.txt', row.names=F, col.names=F, quote=F)

# Running enricher and selecting pathways of interest
em <- enricher(target_eids, TERM2GENE=m_t2g, pvalueCutoff = 1)
canonical_res <- em@result[sapply(em@result$ID, function(x) 
                    grepl('KEGG|REACTOME|WP', x, perl=T)), ]
canonical_res$p.adjust = p.adjust(canonical_res$pvalue, method='fdr')

corrected_em = em
corrected_em@result = canonical_res
cp_dp <- dotplot(corrected_em)
cp_dp

go_enrich = enrichGO(target_eids, org.Hs.eg.db, ont="BP", 
                     pvalueCutoff=1, qvalueCutoff = 0.2)
bp_dp <- dotplot(go_enrich)
bp_dp

f2 <- plot_grid(bp_dp, cp_dp, nrow=1, labels=c('a', 'b'))

cairo_pdf('Figure_2.pdf', width=16, height=5)
print(f2)
dev.off()
