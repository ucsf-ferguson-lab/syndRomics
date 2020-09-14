#install.packages(c("devtools", "roxygen2", "testthat", "knitr"))
#devtools::install_github("r-lib/devtools")

#libraries================================================
# library(tidyverse)
library(tidyverse)
library(magrittr)
library(devtools)
library(reshape2)
library(png)
library(magick)
library(ggrepel)
library(pdftools)
library(grid)
library(bootSVD)
library(ggnewscale)
library(Gifi)
library(progress)
library(profvis)
library(parallel)
library(pbapply)

pca<-prcomp(mtcars, center = T, scale. = T)
pca_data<-mtcars

a<-permut_pc_test(pca, pca_data,ndim=6, statistic = 'commun')


permut_pca_D<-function(i=1,pca_data, pb=NULL,...){

  if(!is.null(pb)){
    pb$tick()
  }
  x<-pca_data
  x<-apply(x, 2,function(x){base::sample(x, length(x))})

  if(class(pca)[1]=='prcomp'){
    if (is.numeric(pca$center)){
      center=T
    }else{
      center=F
    }

    if (is.numeric(pca$scale)){
      .scale=T
    }else{
      .scale=F
    }
    pca_per<-prcomp(x,scale. = .scale,center = center)
  }else if (class(pca)[1]=='princals'){
    pca$call$data<-quote(x)
    pca_per<-eval(pca$call)
  }

  if (output=='s.loadings'|| output=='commun'){
    return(stand_loadings(pca_per, x))
  }else if (output=='VAF'){
    if (class(pca)[1]=='prcomp'){
      VAF<-pca_per$sdev^2/sum(pca_per$sdev^2)
    }else if (class(pca)[1]%in%c('princals')){
      VAF<-pca_per$evals/sum(pca_per$evals)
    }
    return(VAF)
  }
}

per_samples<-replicate(P, permut_pca_D(pca_data = pca_data), simplify = F)
per_samples2<-pbsapply(1:P, FUN = permut_pca_D, pca_data,
                       simplify = F)


cl <- makeCluster(15)
clusterSetRNGStream(cl)

pb <- progress_bar$new(
  format = "[:bar] :percent ~remaining :eta",
  total = P+1, clear = FALSE, width= 80)

pca_data<-cbind(mtcars,mtcars,mtcars,mtcars,mtcars,
                mtcars,mtcars,mtcars,mtcars,mtcars,
                mtcars,mtcars,mtcars,mtcars,mtcars,
                mtcars,mtcars,mtcars,mtcars,mtcars,
                mtcars,mtcars,mtcars,mtcars,mtcars,
                mtcars,mtcars,mtcars,mtcars,mtcars,
                mtcars,mtcars,mtcars,mtcars,mtcars,
                mtcars,mtcars,mtcars,mtcars,mtcars,
                mtcars,mtcars,mtcars,mtcars,mtcars)
pca_data<-cbind(pca_data,pca_data,pca_data,pca_data,
                pca_data,pca_data,pca_data,pca_data,
                pca_data,pca_data,pca_data,pca_data,
                pca_data,pca_data,pca_data,pca_data)

microbenchmark::microbenchmark(replicate(P, permut_pca_D(pca_data = pca_data,
                                                         pb=pb), simplify = F),
                               {
                                 clusterExport(cl,"permut_pca_D")
                                 clusterExport(cl,"pca_data")
                                 clusterExport(cl,"pca")
                                 clusterExport(cl,"output")
                                 pbsapply(1:P, FUN = permut_pca_D, pca_data=pca_data,
                                      simplify = F,cl = cl)},times = 1)
P=100
stopCluster(cl)


load_df<-extract_loadings(pca, mtcars)

syndromic_plot(pca, mtcars, cutoff = 0.8,
               ndim = 1, colors=c("#FF0000",'purple'))

heatmap_loading(pca,mtcars)
barmap_loading(pca,mtcars, ndim = 3,colors=c('green','purple'))
barmap_commun(pca,mtcars, ndim = 3,color='green')

a<-pc_stability(pca, mtcars, colors=c('green','purple'))


## Progress bar
total <- 20
# create progress bar
pb <- txtProgressBar(min = 0, max = total, style = 3, width = 80)
for(i in 1:total){
  Sys.sleep(0.1)
  # update progress bar
  setTxtProgressBar(pb, i, label = 'hi')
  cat('hi')
}
close(pb)

##syndromics color:
extract_syndromic_plot<-function(load_df, pc,cutoff=0.5, VAF, arbitrary_var=NULL,arrow_size_multi=10,
                                 repel=F,plot_legend=T,text_size=6,
                                 var_order='abs decreasing', colors=c("steelblue1","firebrick1")){

  old_scipen<-getOption('scipen')
  on.exit(options(scipen=old_scipen))

  options(scipen=999)

  # img<-png::readPNG(system.file("PCVenn.png", package = "syndRomics"))
  # g<-grid::rasterGrob(img,interpolate=TRUE, just = 'center')

  # var_order=ran2
  p<-load_df

  if (length(p$loading)<1){
    stop(paste(pc, ' has no loadings above cutoff'))
  }

  if (is.character(var_order)&length(var_order)==1){
    if (var_order=='abs decreasing'){
      p<-p%>%dplyr::arrange(desc(abs(.data$loading)))
    }else if(var_order=='abs increasing'){
      p<-p%>%dplyr::arrange(abs(.data$loading))
    }else if(var_order=='decreasing'){
      p<-p%>%dplyr::arrange(desc(.data$loading))
    }else if(var_order=='increasing'){
      p<-p%>%dplyr::arrange(.data$loading)
    }
  }else if(is.numeric(var_order)){
    if (length(var_order)!=length(p$Variables)){
      stop('var_order is of different lenght than the number of variables')
    }
    p<-p[var_order,]
  }else if (is.character(var_order) & length(var_order)>1){
    order_var<-match(var_order, p$Variables)
    if (sum(is.na(order_var))>0){
      stop(paste('variable/s not found: ',
                 var_order[which(is.na(order_var))], sep = ''))
    }
    p<-p[order_var,]
  }

  p<-p%>%dplyr::filter(.data$component==pc, abs(.data$loading)>=cutoff)

  if (dim(p)[1]<1){
    stop(paste('There is no loading above cutoff for ',pc,sep = ''))
  }

  p<-p%>%
    dplyr::mutate(
      div=0:(n()-1),angle=.data$div*2*pi/n()+pi/2,
      angletext=div*2*pi/n()+pi/2-pi/22,
      xend=3.5*cos(.data$angle), yend=3.5*sin(.data$angle),
      x=7*cos(.data$angle), y=7*sin(.data$angle),
      xtext=9*cos(.data$angle), ytext=9*sin(.data$angle),
      xload=6*cos(.data$angletext), yload=6*sin(.data$angletext),
      angletext=(.data$angle*180/pi)+ifelse(.data$x<0, 180,0),#
      loading_txt=as.character(round(.data$loading,3)),
      arrow_weight=(round(.data$loading,3)),
      hadjust=ifelse(.data$x<0, 'right','left'),
      hadjust=ifelse(.data$y<(-3) | .data$y>3, 'center',hadjust)
    )

  p[p$Variables%in%arbitrary_var,'loading_txt']<-
    paste("|",
          abs(round(p[p$Variables%in%arbitrary_var,'loading'],3)),
          "|",sep = "")
  arbitrary_df<-p%>%
    mutate(arrow_weight=ifelse(.data$Variables%in%arbitrary_var,abs(.data$loading),NA))

  load_df<-p%>%
    mutate(loading=ifelse(.data$Variables%in%arbitrary_var,NA,.data$loading))

  angle1<-seq(0,1.2,length.out = 50)
  angle2<-seq(1.95,3.18,length.out = 50)
  angle3<-seq(4.02,5.4,length.out = 50)

  pol<-data.frame(x=c(2.8*cos(angle1)-1,2.8*cos(angle2)+1,2.8*cos(angle3)),
                  y=c(2.8*sin(angle1)-1,2.8*sin(angle2)-1,2.8*sin(angle3)+1))

  s_plot<-ggplot2::ggplot(p,aes(color=.data$loading, label=.data$Variables, x=.data$x, y=.data$y, xend=.data$xend, yend=.data$yend))+
    # ggplot2::annotation_custom(g, xmin = -2, xmax = 2, ymin=-2, ymax=2.4)+
    geom_polygon(data=pol,aes(x,y),inherit.aes = F, fill='grey')+
    ggplot2::scale_color_gradient2(name = "Loading",
                                   high = colors[2], mid = "white", low = colors[1],
                                   midpoint=0,na.value = 'transparent') +
    ggplot2::geom_segment(
      arrow = arrow(type='closed', length = unit(0.3, 'cm'), angle = 25),
      size=abs(p$arrow_weight)*arrow_size_multi, show.legend = F,
      linejoin = 'mitre')+
    ggplot2::ylab(NULL)+
    ggplot2::xlab(NULL)+
    ggplot2::theme(axis.text = element_blank(),panel.background = element_blank(),
                   axis.ticks= element_blank(),axis.line = element_blank(),
                   legend.text = element_blank(),legend.direction = 'horizontal',
                   legend.position = 'bottom',text = element_text(size=text_size*2))+
    ggplot2::xlim(-12,12)+
    ggplot2::ylim(-12,12)+
    coord_equal()


  legend_res<-0.005
  legend_df<-data.frame(x=seq(-3,3,legend_res),
                        z=seq(-1,1,legend_res/3))%>%
    dplyr::mutate(xend=ifelse(.data$x<=0, .data$x-legend_res, .data$x+legend_res), y=rep(-11, length(.data$x)),
                  yend=rep(-11, length(.data$x)))

  if (plot_legend){
    s_plot<-s_plot+ggplot2::geom_segment(data=legend_df, aes(x=.data$x, y=.data$y, xend=.data$xend,
                                                             yend=.data$yend,color=.data$z, size=abs(.data$z)*20),
                                         inherit.aes = F, show.legend = F)+
      ggplot2::geom_segment(aes(x=max(legend_df$x), y=-11, xend=max(legend_df$x)+0.1,
                                yend=-11),color=colors[2],arrow = arrow(type='closed',
                                                                           length = unit(0.1, 'cm'), angle = 25),
                            size=6, show.legend = F,linejoin = 'mitre')+
      ggplot2::geom_segment(aes(x=min(legend_df$x), y=-11, xend=min(legend_df$x)-0.1,
                                yend=-11),color=colors[1],arrow = arrow(type='closed',
                                                                           length = unit(0.1, 'cm'), angle = 25),
                            size=6, show.legend = F,linejoin = 'mitre')+
      ggplot2::annotate(geom='text', x=-4.1, y=-11, label="-1", size=text_size)+
      ggplot2::annotate(geom='text', x=4.1, y=-11, label="1", size=text_size)
  }

  if (!is.null(arbitrary_var)){
    s_plot<-s_plot+
      ggnewscale::new_scale_color() +
      ggplot2::geom_segment(data=arbitrary_df,aes(color=.data$arrow_weight),
                            arrow = arrow(type='closed', length = unit(0.3, 'cm'), angle = 25),
                            size=abs(p$arrow_weight)*arrow_size_multi, show.legend = F,
                            linejoin = 'mitre')+
      scale_color_gradient(low='white', high='grey30',limits=c(0,1), na.value =  "transparent")
  }

  s_plot<-s_plot+
    ggplot2::annotate(geom='text',x=0, y=0.25, vjust = 'center', hjust = 'center', color='black', label=pc, size=text_size,)+
    ggplot2::annotate(geom='text',x=0, y=-0.75, vjust = 'center', hjust = 'center',color='black', label=VAF,size=text_size)+
    ggplot2::geom_text(aes(label=.data$loading_txt, x=.data$xload, y=.data$yload, angle=.data$angletext),
                       color='black',size=text_size*0.7)

  if (repel){
    s_plot<-s_plot+ggrepel::geom_text_repel(aes(x=.data$xtext, y=.data$ytext,label=.data$Variables, hjust=hadjust),
                                            color='black',min.segment.length = 0.5,size=text_size)
  }else{
    s_plot<-s_plot+ggplot2::geom_text(aes(x=.data$xtext, y=.data$ytext,label=.data$Variables,hjust=hadjust),
                                      color='black',size=text_size)
  }
  return(s_plot)
}


install_github(repo = "ucsf-ferguson-lab/SyndRomics",
               auth_token=token_ucsf_ferguson_lab, build_vignettes = T)

library(syndRomics)
library(parallel)
cl <- makeCluster(detectCores()-1)
clusterEvalQ(cl,library(syndRomics))
clusterExport(cl,'stand_loadings')


syndromic_plot(pca, mtcars, cutoff = 0.4,
               ndim = 1, var_order = 'increasing')

barmap_loading(pca,mtcars, ndim = 3, vars=ran[1:5])
heatmap_loading(pca,mtcars, ndim = 3,vars = ran[1:5])

b<-pc_stability(pca, mtcars,B = 100,barmap_plot = T)
b$PC_similarity
b$boot_barmap_loadings

barmap_loading(pca, mtcars, load_list =  b$boot_samples, plot_list_original = T, plot_list_center = T)

system.time(b<-pc_stability2(pca, mtcars,B = 10000, inParallel = F))
system.time(b<-pc_stability2(pca, mtcars,B = 10000, inParallel = T))

b$boot_barmap
p2<-permut_pc_test(pca,mtcars, statistic = 'VAF',P=300, perm.method = 'permV')
p<-permut_pc_test(pca,mtcars, statistic = 's.loadings',P=1000,
                  perm.method = 'permV')

p$per_plot

p2$per_plot
p$results

pc<-permut_pc_test(pca,mtcars, statistic = 's.loadings',P=30,ndim = 2, perm.method = 'permV')
pc<-permut_pc_test(pca,mtcars, statistic = 'commun',P=30,ndim=2, perm.method = 'permV')

pc$per_plot

barmap_commun(pca, mtcars,load_list = pc$per_list, ndim = 2, plot_list_center = T, plot_original = F)

stand_loadings(pca, mtcars)

b_load2%>%
  ggplot(aes(.data$Variables, .data$value, color=.data$component, group=.data$component))+
  geom_line()+
  # scale_fill_gradient2(low='blue3', high='red3',limits=c(-1,1), na.value =  "transparent")+
  labs(title=plot_title,fill=legend_title)+
  ylab(NULL)+xlab(NULL)+
  theme_minimal()+
  # coord_flip()+
  # facet_grid(~component)+
  geom_hline(yintercept = 0, color='black')+
  geom_errorbar(data=b_load_ci, aes(ymin=.data$ymin, ymax=.data$ymax),width=0.5)

b_load2%>%
  ggplot(aes(.data$Variables, .data$component, fill=.data$value))+
  geom_tile()+
  # scale_fill_gradient2(low='blue3', high='red3',limits=c(-1,1), na.value =  "transparent")+
  labs(title=plot_title,fill=legend_title)+
  ylab(NULL)+xlab(NULL)+
  theme_minimal()+
  # coord_flip()+
  # facet_grid(~component)+
  geom_hline(yintercept = 0, color='black')+
  geom_errorbar(data=b_load_ci, aes(ymin=.data$ymin, ymax=.data$ymax),width=0.5)

heatmap(as.matrix(original_loadings),Colv = F)

