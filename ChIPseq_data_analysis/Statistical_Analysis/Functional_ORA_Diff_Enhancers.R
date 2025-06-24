################################################################################
##        Title : ChIP-seq data - Functional analysis
##  Description : This script is for conducting Overrepresentation Analysis over
##                a list of genes annotated to diff enhancers (i.e. KO vs WT) or
##                diff binded regions
##   Researcher : Luis Galán
##         Date : 22nd Mar 2023
################################################################################

# Functional analysis will be jointly conducted for the six type of classified
# enhancers

# REMARK: Althoug this functional analysis is also valid for differential binded
# regions, code refers to diff enhancers

################################################################################
## 1. Load packages
################################################################################

require(clusterProfiler)
require(org.Mm.eg.db)
require(openxlsx)
require(readxl)
library(tidyverse)
require(enrichplot)
require(ggpubr)
library(UpSetR)

################################################################################
## 3. Import data: annotated genes
################################################################################

# Diff enhancers (check 'ChIPseq_Histones_Differential_Enhancers.R' script for 
# obtaining this excel file. Or consider the last 4 sheets from Suppl Table 2

path <- "results/Differential_Enhancers_in_IkBa_KO_samples.xlsx"

diff.enh <- path %>% 
  excel_sheets() %>% 
  set_names() %>% 
  map(read_excel, path = path)

# Each element: one type of diff enhancer class


################################################################################
## 4. WikiPathways overrepresentation analysis
################################################################################

ora.wiki <- lapply(diff.enh, function(genes) {
   res <- enrichWP(gene = genes$ENTREZID,
                   organism   = 'Mus musculus', 
                  pvalueCutoff = 1,
                  qvalueCutoff = 1,
                  pAdjustMethod = "BH")
   res <- setReadable(res, 'org.Mm.eg.db', 'ENTREZID')
})   


################################################################################
## FINAL. Store all results in an excel file
################################################################################

hs <- createStyle(textDecoration = "BOLD", fontColour = "#FFFFFF", fontSize = 12, fontName = "Arial Narrow", fgFill = "#4F80BD")

marks <- names(diff.enh)

sapply(marks, function(mark){
  
  results <- list("WikiPathways" = as.data.frame(ora.wiki[[mark]]))
  
  write.xlsx(x = results, 
             file= paste("results/ORA/ORA_results_DiffENHANCERS_KO_vs_WT_",mark,".xlsx", sep=""), 
             headerStyle=hs)
})


################################################################################
## Plot some graphs based on WikiPathways 
################################################################################

results = ora.wiki


# Dotplot (from clusterProfiler) of the subset scenario
plots <- lapply(names(results), function(scenario){
 # Results are already sorted per adjusted pvalue
  long = length(which(as.data.frame(results[[scenario]])$p.adjust < 0.15))
  dotplot(results[[scenario]], showCategory = ifelse(long >10, 10, long),
          title=scenario) +
    scale_colour_gradient(limits=c(0, 0.15), low="red") +
    theme_bw(base_size = 10) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.text.y = element_text(size=9, face="bold"),
          axis.text.x = element_text(size=10, face="bold")) +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 24)) +
    scale_x_continuous(expand = expansion(mult=0, add=0.01)) 
})

p <- ggarrange(plotlist = plots,
               ncol = 2, nrow = 3,labels = LETTERS[1:length(plots)])

pdf(file = "results/ORA/Dotplot_ORA_DiffBind_Wiki_All_Diff_Enhancers.pdf",
    width=9, height=10)
print(p)
dev.off()
