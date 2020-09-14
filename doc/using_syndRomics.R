## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(syndRomics)
library(ggplot2)
library(dplyr)
library(stringr)

## ---- fig.width=6, fig.height=5-----------------------------------------------
odc.df<-read.csv('odc-sci_26.csv', na.strings = '')
odc.df.filter<-odc.df[-1,str_detect(colnames(odc.df),"_wk6")]
odc.df.filter<-odc.df.filter[complete.cases(odc.df.filter),] #delete missing
colnames(odc.df.filter)<-str_remove(colnames(odc.df.filter), "_wk6")
colnames(odc.df.filter)<-str_remove(colnames(odc.df.filter), "AbsDev")#shorten var names

odc.df.filter[1:10,1:6]

## -----------------------------------------------------------------------------
pca_data<-as.matrix(apply(odc.df.filter, 2, as.numeric))
pca<-prcomp(pca_data, center = T, scale. = T)

## -----------------------------------------------------------------------------
original_loadings<-syndRomics:::stand_loadings(pca, pca_data)
original_loadings[1:10,1:5]

## -----------------------------------------------------------------------------
set.seed(500)
per<-permut_pc_test(pca, pca_data, P=100, ndim=5)
per$results

## ---- fig.width=4, fig.height=4-----------------------------------------------
per$per_plot

## -----------------------------------------------------------------------------
s_per<-permut_pc_test(pca, pca_data, P=100, ndim=3, statistic = 's.loadings', perm.method = 'permV')
s_per$results$PC1
s_per$results$PC2
s_per$results$PC3

## ---- fig.width=8-------------------------------------------------------------
s_per$per_plot

## -----------------------------------------------------------------------------
booted_PCA<-pc_stability(pca, pca_data,ndim = 3,B =200, test_similarity = T, similarity_metric = "all")
booted_PCA$results

## -----------------------------------------------------------------------------
booted_PCA$PC_similarity

## ----fig.width=7, fig.height=5------------------------------------------------
booted_PCA$boot_barmap_loadings

## ----fig.width=7, fig.height=5------------------------------------------------
booted_PCA$boot_barmap_communalities

## ----fig.height=9, fig.width=9------------------------------------------------
splot<-syndRomics::syndromic_plot(pca, pca_data, cutoff = 0.45, text_size = 7)
splot$PC1
splot$PC2
splot$PC3

## ----out.width=500, out.height=500--------------------------------------------
ggsave(plot = splot$PC1, filename = 'PC1.pdf', width = 9, height = 9)
knitr::include_graphics("PC1.pdf")

## ----fig.width=5, fig.height=4------------------------------------------------
h_plot<-heatmap_loading(pca, pca_data, ndim=3, cutoff = 0.45, star_values = T, text_values = F)
h_plot
#The plot can be saved as any other ggplot2 object
#ggsave(plot=h_plot, file="h_plot.png") 

## ----fig.width=8, fig.height=4------------------------------------------------
b_plot<-barmap_loading(pca, pca_data, ndim=3, cutoff = 0.45)
b_plot

#The plot can be saved as any other ggplot2 object
#ggsave(plot=b_plot, file="b_plot.png") 

## ----fig.width=8, fig.height=4------------------------------------------------
b_plot+theme(legend.position = 'bottom', legend.key.width = unit(15,'mm'),
           legend.key.height = unit(2,'mm'),legend.justification = 'left',
           legend.title = element_text(vjust = 1.3),panel.spacing = unit(5,'mm'))

## ----fig.width=6--------------------------------------------------------------
barmap_loading(pca, pca_data, ndim=3, cutoff = 0.45, plot_title = 'Loadings', legend_title = 'Loadings')

## ----fig.width=5, fig.height=4------------------------------------------------
heatmap_loading(pca, pca_data, ndim=3, cutoff = 0.45, star_values = F,text_values = T)

## ----fig.height=7, fig.width=7------------------------------------------------
syndromic_plot(pca, pca_data, ndim=1, cutoff = 0.45, text_size = 5, var_order = 'abs increasing')
syndromic_plot(pca, pca_data, ndim=1, cutoff = 0.45, text_size = 5, var_order = 'decreasing')

## -----------------------------------------------------------------------------
sessionInfo()

