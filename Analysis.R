## libraries ====
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(ggpubr)
##===== Cure fraction ==========================================================
## relative bias in cure fraction violin plots
for (i in 1:4) {
  for (j in 1:4) {
    res = subset(results_cf_0.9, subset = sce == i & model == j)
    res$df = as.factor(res$df)
    assign(
      paste0("rbcf_0.9_sce_", sce, "_model_", model),
      ggplot(res,
             aes( y =
                 bias_relative_cure_frac, x = df
             )) + geom_violin(
               trim = TRUE,
               fill = 'seagreen',
               color = 'seagreen'
             )
    )
    
  }
}
