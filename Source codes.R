# survival analysis -- KM plot ---------------------------------------------

library(tidyverse)
library(survminer)
library(survival)

exp <- read.table("exp.mtx.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
range(exp)
surv <- read.table("survival.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T) #surv.time and surv.stat

a <- intersect(rownames(surv),rownames(exp))
exp <- exp[a,]
surv <- surv[a,]
surv <- cbind(surv,exp)

# optimal cutoff value

cutoff<-surv_cutpoint(surv,
                      time="surv.time",
                      event="surv.stat",
                      variables=c("gene")
);
summary(cutoff)
plot(cutoff, 
     "gene", 
     palette = "lancet")
groups<-surv_categorize(cutoff)
str(groups)
head(groups)
surv$group=surv$gene
table(surv$group)

fitd <- survdiff(Surv(surv.time, surv.stat) ~ group,
                 data      = surv,
                 na.action = na.exclude)
pValue <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit <- survfit(Surv(OS.time, OS)~ group, data = surv)
summary(fit)
p.lab <- paste0("P", ifelse(pValue < 0.001, " < 0.001", paste0(" = ",round(pValue, 3))))

ggsurvplot(fit,
           data = surv,
           pval = p.lab,
           conf.int = TRUE, 
           risk.table = TRUE, 
           risk.table.col = "strata",
           palette = "jco", 
           legend.labs = c("Low", "High"),
           size = 1,
           xlim = c(0,15), 
           break.time.by = 5, 
           legend.title = "gene",
           surv.median.line = "hv", 
           ylab = "Survival probability (%)",
           xlab = "Time (Years)", 
           ncensor.plot = TRUE, 
           ncensor.plot.height = 0.25,
           risk.table.y.text = FALSE)
dev.off()

# correlation analysis --------------------------------------------

df <- read.csv("exp.mtx.csv",row.names = 1)
colnames(df)
df$var1 <- as.numeric(df$var1)
df$var2 <- as.numeric(df$var2)

library(ggplot2)
library(ggpubr)
library(ggpointdensity)

colnames(df)

plot_with_statistics <- function(df, 
                                 x_col = "var1", 
                                 y_col = "var2", 
                                 color_palette = c("#dddddd", "#b9c3ce", "#95a9be", "#6486a9", "#4c749f"),
                                 show_fit_line = FALSE, 
                                 text_size = 14, 
                                 line_size = 0.5) {
  
  p <- ggplot(data = data, mapping = aes_string(x = x_col, y = y_col)) +
    geom_pointdensity(size = 2) +
    scale_color_gradientn(colors = color_palette) +
    theme_classic(base_size = text_size) +
    labs(x = x_col, y = y_col) +
    theme(
      legend.position = "",
      legend.title = element_text(size = text_size - 2),
      legend.text = element_text(size = text_size - 4),
      axis.title = element_text(size = text_size, face = "bold"),
      axis.text = element_text(size = text_size - 2, face = "bold"),
      plot.margin = unit(c(1, 1, 1, 1), "lines")  
    )
  
  if (show_fit_line) {
    p <- p +
      geom_smooth(method = "lm", formula = y ~ x, color = "black", size = line_size) +
      stat_cor(
        method = "spearman",
        aes(label = paste(..r.label.., ..p.label.., sep = "~`, `~")),
        label.x = 1.55,  
        size = text_size / 2.8,
        fontface = "bold"
      )
  }
  
  return(p)
}


plot <- plot_with_statistics(df, show_fit_line = TRUE)
print(plot)

# immune cell inference ------------------------------------

library(estimate)
library(MCPcounter)
library(dplyr)
library(tibble)
library(ggplot2)
library(tidyr)
library(mMCPcounter)
library(xCell)
library(readr)
library(limma)
library(tidyverse)
library(CIBERSORT)

## load data ------------------------------------------------------

expr_matrix <- read.table("exp.mtx.txt", header = T, sep = "\t")
expr_matrix <- as.matrix(expr_matrix)
rownames(expr_matrix) <- expr_matrix[,1]
expr_matrix <- expr_matrix[,2:ncol(expr_matrix)]
dimnames <- list(rownames(expr_matrix),colnames(expr_matrix))
expr_matrix <- matrix(as.numeric(as.matrix(expr_matrix)),nrow=nrow(expr_matrix),dimnames=dimnames)
expr_matrix <- avereps(expr_matrix)

## mouse to human gene ------------------------------------------------------

# 1. load Homologous gene mapping file

orthologs <- read_tsv("Human_Mouse_geneTrans.txt",  #downloaded from NCBI
                      na = c("", "NA", "N/A", "null", "NULL"),  # 将常见空值视为 NA
                      col_types = cols()  # 避免类型猜测警告
) 

colnames(orthologs) <- c("human_stable_id", "human_gene_name",
                         "mouse_stable_id", "mouse_gene_name")

# delete duplicated mouse_gene_name
orthologs_unique <- orthologs[!duplicated(orthologs$mouse_gene_name), ]

# mouse_gene_name -> human_gene_name
mouse_to_human_map <- setNames(
  orthologs_unique$human_gene_name,
  orthologs_unique$mouse_gene_name
)

# mouse gene names

current_mouse_genes <- rownames(expr_matrix)

# mapping

mapped_human_genes <- mouse_to_human_map[current_mouse_genes]

# exclude NA and duplicate

valid_idx <- !is.na(mapped_human_genes)
mouse_expr_mapped <- expr_matrix[valid_idx, , drop = FALSE]
new_rownames <- mapped_human_genes[valid_idx]

if (any(duplicated(new_rownames))) {
  warning("存在多个 mouse genes 映射到同一个人类基因！仅保留第一个出现的。")
  keep_idx <- !duplicated(new_rownames)
  mouse_expr_mapped <- mouse_expr_mapped[keep_idx, , drop = FALSE]
  new_rownames <- new_rownames[keep_idx]
}

rownames(mouse_expr_mapped) <- new_rownames

human_expr <- mouse_expr_mapped


## MCPcounter ----------------------------------------------------

mcp_results <- mMCPcounter.estimate(exp = expr_matrix, features = "Gene.Symbol")

# 转置结果以便于分析
mcp_scores <- as.data.frame(mcp_results)

## xCell ----------------------------------------------

xcell_results <- xCellAnalysis(human_expr)

xcell_scores <- as.data.frame(xcell_results)

## ESTIMATE ------------------------------------------------

# prepare input
write.table(human_expr, file = "expr_input.txt", sep = "\t", quote = FALSE, 
            row.names = TRUE, col.names = TRUE)

# run ESTIMATE
filterCommonGenes(input.f = "expr_input.txt", 
                  output.f = "expr_genes.gct", 
                  id = "GeneSymbol")

estimateScore(input.ds = "expr_genes.gct",
              output.ds = "ESTIMATE_scores.gct",
              platform = "affymetrix")

# load ESTIMATE results
estimate_results <- read.table("ESTIMATE_scores.gct", skip = 2, header = TRUE, 
                               row.names = 1, sep = "\t")
estimate_scores <- estimate_results[,-1]

# Gene-set enrichment analysis (GSEA) -------------------------------------

library(clusterProfiler)
library(edgeR)
library(limma)

## load data ------------------------------------------------------

exp_matrix <- read.table("exp.mtx.txt", header = T, sep = "\t")
exp_matrix <- as.matrix(exp_matrix)
rownames(exp_matrix) <- exp_matrix[,1]
exp_matrix <- exp_matrix[,2:ncol(exp_matrix)]
dimnames <- list(rownames(exp_matrix),colnames(exp_matrix))
exp_matrix <- matrix(as.numeric(as.matrix(exp_matrix)),nrow=nrow(exp_matrix),dimnames=dimnames)
exp_matrix <- avereps(exp_matrix)

## DEG ----------------------------------------------------------

group <- factor(c(rep("Group1", 3), rep("Group2", 3)))
dge <- DGEList(counts = exp_matrix, group = group)

# filter low-expressed genes
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]

dge <- calcNormFactors(dge)
dge <- estimateDisp(dge)
et <- exactTest(dge)

# DEG
results <- topTags(et, n = nrow(dge$counts))
diff_results <- results$table
diff_results$gene_id <- rownames(diff_results)

#ordering
diff_results <- diff_results[, c("gene_id", "logFC", "logCPM", "PValue", "FDR")]
diff_results <- diff_results[order(diff_results$logFC, diff_results$FDR, decreasing = c(TRUE, FALSE)), ]

## GSEA -------------------------------------------------------

library(msigdbr)

hallmark <- msigdbr(species = "Mus musculus", category = "H") %>% 
  dplyr::select(gs_name, gene_symbol) %>% as.data.frame()

PID <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "PID") %>% 
  dplyr::select(gs_name, gene_symbol) %>% as.data.frame()

REACTOME <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "REACTOME") %>% 
  dplyr::select(gs_name, gene_symbol) %>% as.data.frame()

BIOCARTA <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "BIOCARTA") %>% 
  dplyr::select(gs_name, gene_symbol) %>% as.data.frame()

geneList <- diff_results[order(diff_results$logFC, decreasing = T),]$logFC
names(geneList) <- rownames(diff_results[order(diff_results$logFC, decreasing = T),])

hallmark_GSEA <- GSEA(geneList = geneList, minGSSize = 5, pvalueCutoff = 1, pAdjustMethod = "fdr", TERM2GENE = hallmark)
hallmark_GSEA <- hallmark_GSEA@result

PID_GSEA <- GSEA(geneList = geneList, minGSSize = 5, pvalueCutoff = 1, pAdjustMethod = "fdr", TERM2GENE = PID)
PID_GSEA <- PID_GSEA@result

REACTOME_GSEA <- GSEA(geneList = geneList, minGSSize = 5, pvalueCutoff = 1, pAdjustMethod = "fdr", TERM2GENE = REACTOME)
REACTOME_GSEA <- REACTOME_GSEA@result

BIOCARTA_GSEA <- GSEA(geneList = geneList, minGSSize = 5, pvalueCutoff = 1, pAdjustMethod = "fdr", TERM2GENE = BIOCARTA)
BIOCARTA_GSEA <- BIOCARTA_GSEA@result
