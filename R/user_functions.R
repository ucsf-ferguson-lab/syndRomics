guide_colourbar <- function(...) ggplot2::guide_colourbar(...)
guide_legend <- function(...) ggplot2::guide_legend(...)
guide_colorbar <- function(...) ggplot2::guide_colorbar(...)

#'@title Syndromic plot
#'
#'@description Extract the syndromic plots from a pca solution or from a table of loadings for all the
#'specified PCs.
#'
#'@param pca Object of class \emph{prcomp}, \emph{princals}, or \emph{data.frame}. If object is a \emph{prcomp} or \emph{princals} object, \emph{pca_data} argument is required, and
#'the loadings will be extracted. If object is a data.frame object, the dataframe needs to be formatted as: first column named \emph{Variables} and all other columns corresponding to a PC.
#' One row per variable. The values are the loadings.
#'
#'@param pca_data Data passed to the \emph{prcomp} or \emph{princals} function.
#'@param ndim Numeric. Number of PCs to plot.
#'@param cutoff Numeric or numeric vector of length \emph{ndim}. Value of the loadings threshold to plot. If only one value is passed,
#'the same cutoff will be used for all the PCs. If a vector passed, each value will be used for the corresponding PC.
#'@param VAF If \emph{pca} is from prcomp, VAF argument is not needed. Otherwise, VAF is a String vector with text for the centers of the syndromic plots. The text generally corresponds to the variance accounted for (VAF) for the respective PCs.
#'@param arbitrary_var Character or character vector with the names of the variables where loadings should be plot as absolute values.
#'This is the case for categorical variables in categorical PCA where variables do not have direction.
#'@param arrow_size_multi Numeric. Controls the size of the arrows proportional to the loading. Default=10
#'@param repel Boolean. Whether to repel the text for preventing text overlap. Default=TRUE
#'@param plot_legend Boolean. Whether to plot the legend or not. Default=TRUE
#'@param text_size Numeric. Controls for the size of the text. Default=9
#'@param var_order Character. Specify the order of the variables in the plot by the loading values, starting at 12 o’clock and moving counterclockwise. Possible values: 'abs decreasing': plot by decreasing absolute value;
#''abs increasing': plot by increasing absolute value; 'decreasing'; or 'increasing’.
#'@param ... Other arguments passed to \emph{extract_syndromic_plot()}
#'
#'@return Returns a list of \emph{ggplot2} objects with one element for each PC plot. It also renders and saves the plots as *.pdf
#'in the working directory or in the specified \emph{path}.
#'
#'@examples
#'pca<-prcomp(mtcars, center = TRUE, scale. = TRUE)
#'syndromic_plot(pca, mtcars, cutoff = 0.65, ndim=1)
#'
#'
#'@export
#'
#'@import dplyr tidyr ggnewscale ggplot2
#'@importFrom rlang .data
#'
syndromic_plot<-function(pca, pca_data=NULL, ndim=3, cutoff, VAF,arbitrary_var=NULL,arrow_size_multi=10,
                         repel=T, plot_legend=T, text_size=9,
                         var_order='abs decreasing',...){

  load_df<-extract_loadings(pca, pca_data)

  if (class(pca)[1]=='prcomp'){
    VAF<-paste(round(pca$sdev^2/sum(pca$sdev^2)*100,1), '%', sep = '')
  }else if (class(pca)[1]%in%c('princals')){
    VAF<-paste(round(pca$evals/sum(pca$evals)*100,1), '%', sep='')
  }

  ndim<-min(dim(load_df)[2]+1, ndim)

  load_df<-load_df%>%pivot_longer(-.data$Variables,
                                  names_to = "component", values_to = "loading")%>%
    filter(.data$component%in%unique(.data$component)[1:ndim])

  s_list<-list()
  PC<-unique(load_df$component)
  for (i in 1:length(PC)){
    pc<-PC[i]
    c<-ifelse(length(cutoff)>1,cutoff[i], cutoff)
    v<-VAF[i]
    s_plot<-NULL

    try(s_plot<-extract_syndromic_plot(load_df = load_df, pc=pc, cutoff=c, VAF=v, text_size=text_size,
                                       arrow_size_multi = arrow_size_multi, repel = repel,arbitrary_var=arbitrary_var,
                                       var_order=var_order, plot_legend = plot_legend,...))

    s_list[[pc]]<-s_plot
  }
  return(s_list)
}

#'@title Heatmap of standardized loadings
#'
#'@description Plot a heatmap of the standardized loadings from a PCA solution.
#'
#'@author Abel Torres Espin
#'
#'@param pca Object of class \emph{prcomp}, \emph{princals}, or \emph{data.frame}. If object is a \emph{prcomp} or \emph{princals} object, \emph{pca_data} is required, and
#'the loadings will be extracted. If object is a data.frame object, the dataframe needs to be formatted as: first column named \emph{Variables} and all other columns corresponding to a PC.
#' One row per variable. The values are the loadings.
#'
#'@param pca_data Data passed to the \emph{prcomp} or \emph{princals} function.
#'@param ndim Numeric. Number of PCs to plot
#'@param cutoff Numeric or numeric vector of length \emph{ndim}.
#'Value of the loadings threshold (i.e. |loadings| >= cutoff) to plot with stars. Default = 0.5
#'@param arbitrary_var Character or character vector with the names of the variables where the loadings should be plot as absolute values.
#'This is the case for categorical variables in categorical PCA where variables do not have direction.
#'@param plot_title String. Title of the plot. 'Standardized loadings' by default.
#'@param legend_title String. Title of the legend. 's. loading' by default.
#'@param text_values Boolean. Whether to plot the values of the loadings or not. Default=TRUE
#'@param star_values Boolean. Whether to plot a star in |loadings|>=cutoff. Only relevant if
#'\emph{text_values}=FALSE. Default=FALSE
#'@param text_size Numeric. Size of the text_values.
#'@param vars Character vector. Variables will be ordered as the provided variable names. Non-specified
#'variables will be excluded from the plot. By default variables are ordered in alphabetically by ggplot.
#'
#'@return Returns a \emph{ggplot2} object.
#'
#'
#'@export
#'
#'@import ggplot2 dplyr tidyr ggnewscale
#'@importFrom rlang .data
#'
heatmap_loading<-function(pca, pca_data, ndim=10, cutoff=0.5,arbitrary_var=NULL,plot_title='Standardized loadings',
                          legend_title='s. loading', text_values=T, star_values=F,
                          text_size=2, vars=NULL){

  old_scipen<-getOption('scipen')
  on.exit(options(scipen=old_scipen))

  options(scipen=999)

  load_df<-extract_loadings(pca, pca_data)
  ndim<-min(dim(load_df)[2]+1, ndim)

  if (!is.null(vars)){
    if(!is.character(vars)){
      stop('vars must be a character vector')
    }

    order_var<-match(vars, load_df$Variables)
    if (sum(is.na(order_var))>0){
      stop(paste('variable/s not found: ',
                 vars[which(is.na(order_var))], sep = ''))
    }

    load_df<-load_df%>%filter(.data$Variables%in%vars)%>%
      mutate(Variables=factor(.data$Variables,vars))
  }

  load_df<-load_df%>%pivot_longer(-(.data$Variables),
                                  names_to = "component", values_to = "loading")%>%
    mutate(component=factor(.data$component, levels=names(load_df)),
           star=ifelse(abs(.data$loading)>=cutoff,'*',''),
           loading_txt=as.character(round(.data$loading,2)),
           weight=(round(.data$loading,2)))

  load_df[load_df$Variables%in%arbitrary_var,'loading_txt']<-sapply(unique(load_df$component),function(x){
    paste("|",
          abs(round(load_df[load_df$Variables%in%arbitrary_var & load_df$component==x,'loading'],3)),
          "|",sep = "")
  })

  arbitrary_df<-load_df%>%
    filter(.data$component%in%paste('PC',1:ndim, sep = ''))%>%
    mutate(weight=ifelse(.data$Variables%in%arbitrary_var,abs(.data$loading),NA))

  load_df<-load_df%>%
    mutate(loading=ifelse(.data$Variables%in%arbitrary_var,NA,.data$loading))

  h_plot<-load_df%>%filter(.data$component%in%paste('PC',1:ndim, sep = ''))%>%
    ggplot(aes(.data$component, .data$Variables, fill=.data$loading))+
    geom_raster()+
    scale_fill_gradient2(low='blue3', high='red3',limits=c(-1,1), na.value =  "transparent")+
    labs(title=plot_title,fill=legend_title)+
    ylab(NULL)+xlab(NULL)+
    theme_minimal()

  if (!is.null(arbitrary_var)){

    h_plot<-h_plot+
      ggnewscale::new_scale_fill() +
      geom_raster(data=arbitrary_df,aes(fill=.data$weight))+
      scale_fill_gradient(low='white', high='grey30',limits=c(0,1), na.value =  "transparent")+
      labs(fill= paste("|",legend_title,"|", sep = ""))
  }

  if (text_values){
    h_plot<-h_plot+geom_text(aes(label=.data$loading_txt), size=2)
  }else if (star_values){
    h_plot<-h_plot+geom_text(aes(label=.data$star), color='black')
  }

  return(h_plot)
}

#'@title Barmap of standardized loadings
#'
#'@description Plot a Barmap of the standardized loadings from a PCA solution.
#'
#'@author Abel Torres Espin
#'
#'@param pca Object of class \emph{prcomp}, \emph{princals}, or \emph{data.frame}. If object is a \emph{prcomp} or \emph{princals} object, \emph{pca_data} is required, and
#'the loadings will be extracted. If object is a data.frame object, the dataframe needs to be formatted as: first column named \emph{Variables} and all other columns corresponding to a PC.
#' One row per variable. The values are the loadings.
#'
#'@param pca_data Data passed to the \emph{prcomp} or \emph{princals} function.
#'@param ndim Numeric. Number of PCs to plot
#'@param cutoff Numeric or numeric vector of length \emph{ndim}.
#'Value of the loadings threshold (i.e. |loadings| >= cutoff) to plot with stars. Default = 0.5
#'@param arbitrary_var Character or character vector with the names of the variables where the loadings should be plot as absolute values.
#'This is the case for categorical variables in categorical PCA where variables do not have direction.
#'@param load_list List. List of loading matrices used to plot percentile confidence intervals around the average. If NULL, no error bars are ploted. Default=NULL
#'@param conf Numeric. Confidence level used when \emph{load_list} is provided.
#'@param plot_list_original Boolean. Whether to plot the loadings obtained from the values passed to \emph{pca} and \emph{pca_data} when \emph{load_list} is not NULL.
#'@param plot_list_center Boolean. Whether to plot the average loadings obtained from \emph{load_list}.
#'@param plot_title String. Title of the plot. 'Standardized loadings' by default.
#'@param legend_title String. Title of the legend. 's. loading' by default.
#'@param text_values Boolean. Whether to plot the values of the loadings or not. Default=TRUE
#'@param star_values Boolean. Whether to plot a star in |loadings|>=cutoff. Only relevant if
#'\emph{text_values}=FALSE. Default=FALSE
#'@param text_size Numeric. Size of the text_values.
#'@param plot_cutoff Boolean. Whether to plot the cutoff lines or not.
#'@param vars Character vector. Variables will be ordered as the provided variable names. Non-specified
#'variables will be excluded from the plot. By default variables are ordered in alphabetically by ggplot.
#'
#'@return Returns a \emph{ggplot2} object.
#'
#'@export
#'
#'@import ggplot2 dplyr tidyr ggnewscale
#'@importFrom rlang .data
#'
barmap_loading<-function(pca, pca_data, ndim=10, cutoff=0.5,arbitrary_var=NULL,
                         load_list=NULL,conf=0.95, plot_list_original=F,plot_list_center=F,
                         plot_title='Standardized loadings', legend_title='s. loading',text_values=F, star_values=T,
                         text_size=2, plot_cutoff=T, vars=NULL){

  # load_list<-b$boot_samples
  # pca_data<-mtcars
  # ndim=3
  # cutoff=0.5
  # conf=0.95
  # plot_title='Standardized loadings'
  # legend_title='s. loading'
  # text_values=F
  # star_values=T
  # text_size=2
  # plot_cutoff=T
  # vars=NULL

  old_scipen<-getOption('scipen')
  on.exit(options(scipen=old_scipen))

  options(scipen=999)

  load_df<-extract_loadings(pca, pca_data)
  ndim<-min(dim(load_df)[2]+1, ndim)

  #No loadings list
  if (is.null(load_list)){
    if (!is.null(vars)){
      if(!is.character(vars)){
        stop('vars must be a character vector')
      }

      order_var<-match(vars, load_df$Variables)
      if (sum(is.na(order_var))>0){
        stop(paste('variable/s not found: ',
                   vars[which(is.na(order_var))], sep = ''))
      }

      load_df<-load_df%>%filter(.data$Variables%in%vars)%>%
        mutate(Variables=factor(.data$Variables,vars))
    }

    if (length(cutoff)!=1 & length(cutoff)!=ndim){
      stop('cutoff length must be one or equal to ndim')
    }
    cutoff_df<-load_df%>%pivot_longer(-(.data$Variables),
                                      names_to = "component", values_to = "loading")%>%
      filter(.data$component%in%paste('PC',1:ndim, sep = ''))%>%
      group_by(.data$component)%>%
      summarise(count=n())%>%mutate(cutoff=cutoff)

    load_df<-load_df%>%pivot_longer(-(.data$Variables),
                                    names_to = "component", values_to = "loading")%>%
      filter(.data$component%in%paste('PC',1:ndim, sep = ''))%>%
      left_join(cutoff_df, by='component')%>%
      group_by(.data$component)%>%
      mutate(cutoff=cutoff,
             star=ifelse(abs(.data$loading)>=cutoff,'*',''),
             loading_txt=as.character(round(.data$loading,2)),
             weight=(round(.data$loading,2)))%>%
      ungroup()%>%
      mutate(component=factor(.data$component, levels=names(load_df)))

    load_df[load_df$Variables%in%arbitrary_var,'loading_txt']<-sapply(unique(load_df$component),function(x){
      paste("|",
            abs(round(load_df[load_df$Variables%in%arbitrary_var & load_df$component==x,'loading'],3)),
            "|",sep = "")
    })

    arbitrary_df<-load_df%>%
      filter(.data$Variables%in%arbitrary_var)%>%
      mutate(weight=ifelse(.data$Variables%in%arbitrary_var,abs(.data$loading),NA))

    load_df<-load_df%>%
      mutate(loading=ifelse(.data$Variables%in%arbitrary_var,NA,.data$loading))

    b_plot<-load_df%>%
      ggplot(aes(.data$Variables, .data$loading, fill=.data$loading))+
      geom_col()+
      scale_fill_gradient2(low='blue3', high='red3',limits=c(-1,1), na.value =  "transparent")+
      labs(title=plot_title,fill=legend_title)+
      ylab(NULL)+xlab(NULL)+
      theme_minimal()+
      coord_flip()+
      facet_grid(~component)

    if (!is.null(arbitrary_var)){

      b_plot<-b_plot+
        ggnewscale::new_scale_fill() +
        geom_col(data=arbitrary_df,aes(fill=.data$weight))+
        geom_col(data=arbitrary_df,aes(y=.data$loading*-1,fill=.data$weight))+
        scale_fill_gradient(low='white', high='grey30',limits=c(0,1), na.value =  "transparent")+
        labs(fill= paste("|",legend_title,"|", sep = ""))
    }

    b_plot<-b_plot+geom_hline(yintercept = 0, color='black')

    if(plot_cutoff){
      b_plot<-b_plot+
        geom_hline(data=cutoff_df,aes(yintercept = cutoff), color='red3', alpha=0.6)+
        geom_hline(data=cutoff_df,aes(yintercept = -cutoff), color='blue3', alpha=0.6)
    }

    if (text_values){
      b_plot<-b_plot+geom_text(aes(label=.data$loading_txt,y=.data$loading*1.1),size=2)
      b_plot<-b_plot+geom_text(data=arbitrary_df, aes(label=.data$loading_txt,y=.data$loading*1.1),size=2)
      b_plot<-b_plot+geom_text(data=arbitrary_df, aes(label=.data$loading_txt,y=.data$loading*-1.1),size=2)
    }else if (star_values){
      b_plot<-b_plot+geom_text(aes(label=.data$star, y=.data$loading*1.1), color='black')
      b_plot<-b_plot+geom_text(data=arbitrary_df,aes(label=.data$star, y=.data$loading*1.1), color='black')
      b_plot<-b_plot+geom_text(data=arbitrary_df,aes(label=.data$star, y=.data$loading*-1.1), color='black')
    }
  }

  #PCA list
  if (!is.null(load_list)){

    ndim<-min(dim(load_list[[1]])[2], ndim)
    original_loadings<-stand_loadings(pca, pca_data)[,1:ndim]
    load_list<-lapply(load_list, function(x){
      x[,1:ndim]
    })

    array_per<-array(unlist(load_list), dim=c(dim(original_loadings)[1],
                                             dim(original_loadings)[2],
                                             length(load_list)))

    boot_mean<-apply(array_per, 1:2, mean)
    ci_low<-apply(array_per,1:2,function(x) quantile(x, (1-conf)/2, na.rm = T))
    ci_high<-apply(array_per,1:2,function(x) quantile(x, 1-(1-conf)/2, na.rm = T))

    if (plot_list_original==F){
      b_load<-as.data.frame(boot_mean[,1:ndim])
    }else{
      b_load<-original_loadings
    }

    colnames(b_load)<-colnames(original_loadings)
    b_load$Variables<-rownames(original_loadings)
    b_load2<-b_load%>%pivot_longer(-(.data$Variables))
    colnames(b_load2)[2]<-'component'

    b_load_low<-as.data.frame(ci_low[,1:ndim])
    colnames(b_load_low)<-colnames(original_loadings)
    b_load_low$Variables<-rownames(original_loadings)
    b_load_low<-b_load_low%>%pivot_longer(-(.data$Variables))
    colnames(b_load_low)[3]<-'ymin'

    b_load_high<-as.data.frame(ci_high[,1:ndim])
    colnames(b_load_high)<-colnames(original_loadings)
    b_load_high$Variables<-rownames(original_loadings)
    b_load_high<-b_load_high%>%pivot_longer(-(.data$Variables))
    colnames(b_load_high)[3]<-'ymax'

    b_load_ci<-b_load_low
    b_load_ci$ymax<-b_load_high$ymax
    b_load_ci$value<-b_load2$value
    colnames(b_load_ci)[2]<-'component'

    b_plot<-b_load2%>%
      ggplot(aes(.data$Variables, .data$value, fill=.data$value))+
      geom_col()+
      scale_fill_gradient2(low='blue3', high='red3',limits=c(-1,1), na.value =  "transparent")+
      labs(title=plot_title,fill=legend_title)+
      ylab(NULL)+xlab(NULL)+
      theme_minimal()+
      coord_flip()+
      facet_grid(~component)+
      geom_hline(yintercept = 0, color='black')+
      geom_errorbar(data=b_load_ci, aes(ymin=.data$ymin, ymax=.data$ymax),width=0.5)

    if (plot_list_center){
      b_load<-as.data.frame(boot_mean[,1:ndim])
      colnames(b_load)<-colnames(original_loadings)
      b_load$Variables<-rownames(original_loadings)
      b_load2<-b_load%>%pivot_longer(-(.data$Variables))
      colnames(b_load2)[2]<-'component'

      b_plot<-b_plot+
        geom_point(data=b_load2, aes(y=.data$value, x=.data$Variables),inherit.aes = F)
    }
  }
  return(b_plot)
}

#'@title Barmap of communalities
#'
#'@description Plot a Barmap of the communalities from a PCA solution given the first ndim PCs.
#'
#'@author Abel Torres Espin
#'
#'@param pca Object of class \emph{prcomp}, \emph{princals}, or \emph{data.frame}. If object is a \emph{prcomp} or \emph{princals} object, \emph{pca_data} is required, and
#'the loadings will be extracted. If object is a data.frame object, the dataframe needs to be formatted as: first column named \emph{Variables} and all other columns corresponding to a PC.
#' One row per variable. The values are the loadings.
#'
#'@param pca_data Data passed to the \emph{prcomp} or \emph{princals} function.
#'@param ndim Numeric. Number of PCs to plot
#'@param load_list List. List of loading matrices used to plot percentile confidence intervals around the average. If NULL, no error bars are ploted. Default=NULL
#'@param conf Numeric. Confidence level used when \emph{load_list} is provided.
#'@param plot_original Boolean. Whether to plot the communalities obtained from the values passed to \emph{pca} and \emph{pca_data}.
#'@param plot_list_center Boolean. Whether to plot the average communalities obtained from \emph{load_list}.
#'@param plot_title String. Title of the plot.
#'@param legend_title String. Title of the legend.
#'@param text_values Boolean. Whether to plot the values of the communalities or not. Default=FALSE
#'@param text_size Numeric. Size of the text_values.
#'@param vars Character vector. Variables will be ordered as the provided variable names. Non-specified
#'variables will be excluded from the plot. By default variables are ordered in alphabetically by ggplot.
#'@param var_order Character. Specify the order of the variables in the plot by the communalitys values, starting at 12 o’clock and moving counterclockwise. Possible values: 'abs decreasing': plot by decreasing absolute value;
#''abs increasing': plot by increasing absolute value; 'decreasing'; or 'increasing’.
#'
#'@return Returns a \emph{ggplot2} object.
#'
#'@export
#'
#'@import ggplot2 dplyr tidyr ggnewscale
#'@importFrom rlang .data
#'
barmap_commun<-function(pca, pca_data, ndim=10,load_list=NULL,conf=0.95, plot_original=T,plot_list_center=F,
                        plot_title='Communalities', legend_title='communa.',text_values=F,
                        text_size=2, var_order='increasing',vars=NULL){

  old_scipen<-getOption('scipen')
  on.exit(options(scipen=old_scipen))

  options(scipen=999)

  load_df<-extract_loadings(pca, pca_data)
  ndim<-min(dim(load_df)[2]-1, ndim)

  if (is.null(load_list)){
    plot_original==T
  }

  if (ndim>1){
    commun_df<-data.frame(commun=rowSums(load_df[,1:ndim]^2),
                          Variables=row.names(load_df))
  }else{
    commun_df<-data.frame(commun=load_df[,1:ndim]^2,
                          Variables=row.names(load_df))
  }

  if (!is.null(vars)){
    if(!is.character(vars)){
      stop('vars must be a character vector')
    }

    order_var<-match(vars, commun_df$Variables)
    if (sum(is.na(order_var))>0){
      stop(paste('variable/s not found: ',
                 vars[which(is.na(order_var))], sep = ''))
    }

    commun_df<-commun_df%>%filter(.data$Variables%in%vars)%>%
      mutate(Variables=factor(.data$Variables,vars))
  }else if (var_order=='increasing'){
    commun_df<-commun_df%>%arrange(.data$commun)%>%
      mutate(Variables=factor(.data$Variables, levels = unique(.data$Variables)))
  }else if (var_order=='decreasing'){
    commun_df<-commun_df%>%arrange(desc(.data$commun))%>%
      mutate(Variables=factor(.data$Variables, levels = unique(.data$Variables)))
  }

  c_plot<-commun_df%>%
    ggplot(aes(.data$Variables, .data$commun, fill=.data$commun))+
    scale_fill_gradient2(low='white', high='red3',limits=c(0,1), na.value =  "transparent")+
    labs(title=paste(plot_title, ' for total of ', ndim, ' PCs', sep = ''),
         fill=legend_title)+
    ylab(NULL)+xlab(NULL)+
    theme_minimal()+
    coord_flip()

  if (plot_original){
    c_plot<-c_plot+
      geom_col()
  }

  if (!is.null(load_list)){

    if (ndim>1){
      b_communalities<-lapply(load_list, function(x){
        rowSums(x[,1:ndim]^2)
      })
    }else{
      b_communalities<-lapply(load_list, function(x){
        x[,1:ndim]^2
      })
    }

    array_b<-do.call(rbind, b_communalities)

    mean_communalities<-apply(array_b, 2, mean)
    ci_low_c<-apply(array_b,2,function(x) quantile(x, (1-conf)/2, na.rm = T))
    ci_high_c<-apply(array_b,2,function(x) quantile(x, 1-(1-conf)/2, na.rm = T))

    results_c<-cbind(commun_df,mean_communalities, ci_low_c,ci_high_c)
    com_df<-as.data.frame(results_c)
    com_df$Variables<-row.names(results_c)

    c_plot<-c_plot+
        geom_errorbar(data=com_df, aes(ymin=.data$ci_low_c, ymax=.data$ci_high_c,
                                       x=.data$Variables), inherit.aes = F, width=0.5)
    if (plot_list_center){
      c_plot<-c_plot+
        geom_point(data=com_df, aes(y=.data$mean_communalities, x=.data$Variables), inherit.aes = F)
    }
  }

  if (text_values){
    c_plot<-c_plot+geom_text(aes(label=.data$commun,y=.data$commun*1.1),size=2)
  }
  return(c_plot)
}

#'@title PC stability by nonparametric bootstrapping
#'
#'@description Extract nonparametric estimate of the standardized loadings and the confident
#'region by means of bootstrapping. This function uses the \emph{boot} function from
#'the \strong{boot} package.
#'
#'@author Abel Torres Espin
#'
#'@param pca Object of class \emph{prcomp}, \emph{princals}, or \emph{data.frame}. If object is a \emph{prcomp} or \emph{princals} object, \emph{pca_data} is required, and
#'the loadings will be extracted. If object is a data.frame object, the dataframe needs to be formatted as: first column named \emph{Variables} and all other columns corresponding to a PC.
#' One row per variable. The values are the loadings.
#'
#'@param pca_data Data passed to the \emph{prcomp} or \emph{princals} function.
#'@param ndim Numeric. Number of PCs (1 to \emph{ndim}) to run the analysis on. Default = 3.
#'@param B Numeric. Number of bootstrapped samples passed to \emph{boot}.
#'@param sim Character. Determines the bootstrapping method for \emph{boot}.
#'From \emph{boot}: A character string indicating the type of simulation required.
#'Possible values are "ordinary" (the default), "parametric", "balanced",
#'"permutation", or "antithetic". Default="ordinary".
#'@param communalities Boolean. Whether to compute and return communalities.
#'@param test_similarity Boolean. Whether to compute component similarity by the
#'specified \emph{similarity_metric} and the nonparametric percentile confident interval.
#'See \emph{?component_similarity} for more information. Default=TRUE
#'@param similarity_metric character or character vector. Specify the similarity metric
#'to use when \emph{test_similarity}=TRUE. Possible values are "cc_index" (congruence coefficient),
#'"r_correlation" (Pearson's r), "rmse" (root mean squared error), "s_index' (Cattell's s metric), or "all".
#'See \emph{?component_similarity} for more details on the available metrics. Default="all".
#'@param s_cut_off Numeric. This is the loading cut off used to determine if a
#'variable is silent or not in Cattell's terms. See \emph{?extract_s} for more
#'information. Default=0.1.
#'@param ci_type Character. Type of confidence interval to compute. This argument is passed to
#'the \emph{boot.ci} function from the \emph{boot} package. See \emph{?boot.ci} for options. Given
#'that the BCA method has demonstrated good performance for bootstrapping PCAs,
#'we have set 'bca' as the default. See ref for more details.
#'@param conf Numeric. Level of confidence region for the confidence interval. E.g. 0.95 generates 95CI. Default=0.95
#'@param barmap_plot Boolean. Whether to generate a barmap plot of the bootstrapped loadings or not.
#'See \emph{?barmap_loadings} for details. Default=TRUE
#'@param ... Other arguments passed to the \emph{barmap_loadings} or \emph{barmap_commun} function when \emph{barmap_plot}=TRUE
#'
#'@return Returns a list object.
#'\describe{
#'   \item{\strong{boot_sample}}{A list of length B containing the resampled loadings.}
#'   \item{\strong{results}}{
#'     \describe{
#'        \item{\strong{boot_mean}}{A numeric matrix containing the mean loadings of the bootstrapped sample.}
#'        \item{\strong{ci_low and ci_high}}{A numeric matrix with the lower and upper CI respectively.}
#'     }}
#'   \item{\strong{pc_similarity}}{A list of results when \emph{test_similarity}=TRUE.}
#'   \describe{
#'    \item{\strong{similarity_mean}}{A numerix matrix with the mean of the chosen similarity metric.}
#'    \item{\strong{similarity_ci_low and similarity_ci_high}}{A numeric matrix with the lower and upper CI respectively.}
#'    }
#'  \item{\strong{boot_barmap_loadings}}{A ggplot2 object with a barmap plot of loadings of the bootstrapped solution.}
#'  \item{\strong{boot_barmap_loadings}}{A ggplot2 object with a barmap plot of communalities of the bootstrapped solution.}
#'}
#'
#'@details The number of bootstrap samples is set to 1000 by default, as it has been shown to be a robust number in most conditions of
#'data complexity and sample size. The user must be careful on setting such number too low which would reduce the performance
#'of the approximation. However, values that are too high might unnecessarily increase computing time with little gain (REFs).
#'
#'@references Efron B. Better Bootstrap Confidence Intervals. J Am Stat Assoc. 1987 Mar 1;82(397):171–85.
#'
#'@export
#'
#'@import ggplot2 dplyr boot stringr
#'@importFrom rlang .data
#'@importFrom stats na.omit
#'
pc_stability<-function(pca, pca_data, ndim=3, B=1000, sim='ordinary',communalities=T,
                       test_similarity=T, similarity_metric='all', s_cut_off=0.1,
                       ci_type='bca', conf=0.95,barmap_plot=T,...){

  # pca<-prcomp(mtcars[1:10,], center = T, scale. = T)
  # pca_data<-mtcars[1:10,]
  # ndim=3
  # B=200
  # sim='ordinary'
  # ci_type='bca'
  # conf=0.95
  # test_similarity=T
  # similarity_metric='all'
  # s_cut_off=0.4

  nvars<-dim(pca_data)[2]
  var_names<-colnames(pca_data)
  original_loadings<-stand_loadings(pca, pca_data)
  ndim<-min(dim(original_loadings)[2], ndim)

  cat(paste("Bootstrapping ", B, " times...", sep=''), '\n')

  pb <- dplyr::progress_estimated(B+1)
  # pb <- startpb(0, B)
  # pb_count<<-0
  b_pca<-boot::boot(data = pca_data, pca=pca, statistic = boot_pca_sample, R=B,
              ndim=ndim,original_loadings=original_loadings,
              sim = sim, pb=pb, ...)
  b_pca$t<-na.omit(b_pca$t)
  b_pca$R<-dim(b_pca$t)[1]

  pb$i<-pb$n
  pb$print()

  cat('\n',"Calculating confident intervals...",'\n')

  pb <- dplyr::progress_estimated(dim(b_pca$t)[2])
  ci_results<-list()

  try_ci<-try(for(i in 1:dim(b_pca$t)[2]){
    pb$tick()$print()
    ci<-boot::boot.ci(b_pca, index = i, type = ci_type,conf = conf,...)
    ci_results[[i]]<-ci[[names(ci)[str_detect(names(ci),ci_type)]]][,4:5]
  })

  if (class(try_ci)=='try-error'){
    stop('Computation of confident intervals has failed, probably because B
         is too small. Please increase B and try again.')
  }

  pb$i<-pb$n
  pb$print()

  ci_results<-do.call(rbind,ci_results)


  ci_low<-matrix(c(ci_results[,1]),ncol = ndim)
  row.names(ci_low)<-var_names
  colnames(ci_low)<-paste('PC',1:ndim, sep = '')
  ci_high<-matrix(c(ci_results[,2]),ncol = ndim)
  row.names(ci_high)<-var_names
  colnames(ci_high)<-paste('PC',1:ndim, sep = '')

  array_b<-array(t(b_pca$t), dim=c(dim(original_loadings)[1],
                                   ndim,b_pca$R))

  b_list<-lapply(1:dim(array_b)[3], function(x){
    temp<-array_b[ , 1:ndim, x]
    row.names(temp)<-var_names
    temp
  })

  boot_mean<-apply(array_b, 1:2, mean)
  row.names(boot_mean)<-var_names
  colnames(boot_mean)<-paste('PC',1:ndim, sep = '')

  results_list<-list()
  results_list[['results']]<-list('boot_mean'=boot_mean,
                    'ci_low'=ci_low,'ci_high'=ci_high)
  results_list[['boot_samples']]<-b_list

  if (test_similarity){
    cat('\n',"Calculating PC similarities...",'\n')
    pb <- dplyr::progress_estimated(length(b_list))
    similarity_res<-list()
    for(i in 1:length(b_list)){
      pb$tick()$print()
      load.list<-list(original_loadings[,1:ndim],b_list[[i]][,1:ndim])
      similarity_res[[i]]<-component_similarity(load.list, ndim=ndim, s_cut_off = s_cut_off,
                                                similarity_metric = similarity_metric)$index_mean
    }

    similarity_res<-do.call(rbind, similarity_res)

    ci_sim_low<-similarity_res%>%group_by(.data$PC)%>%
      summarise_all(.funs=function(x) quantile(x, (1-conf)/2, na.rm = T))
    ci_sim_high<-similarity_res%>%group_by(.data$PC)%>%
      summarise_all(.funs=function(x) quantile(x, 1-(1-conf)/2, na.rm = T))
    sim_mean<-similarity_res%>%group_by(.data$PC)%>%
      summarise_all(mean)

    results_list[['PC_similarity']]<-list('similarity_mean'=sim_mean,
                                     'similarity_ci_low'=ci_sim_low,
                                     'similarity_ci_high'=ci_sim_high)
    pb$i<-pb$n
    pb$print()
  }

  if (communalities){
    original_loadings<-stand_loadings(pca, pca_data)[,1:ndim]

    if (ndim>1){
      original_communalities<-rowSums(original_loadings^2)
      b_communalities<-lapply(b_list, function(x){
        rowSums(x[,1:ndim]^2)
      })
    }else{
      original_communalities<-original_loadings^2
      b_communalities<-lapply(b_list, function(x){
        x[,1:ndim]^2
      })
    }

    array_b<-do.call(rbind, b_communalities)

    mean_communalities<-apply(array_b, 2, mean)
    ci_low_c<-apply(array_b,2,function(x) quantile(x, (1-conf)/2, na.rm = T))
    ci_high_c<-apply(array_b,2,function(x) quantile(x, 1-(1-conf)/2, na.rm = T))

    results_c<-cbind(original_communalities,mean_communalities, ci_low_c,ci_high_c)
    colnames(results_c)<-c('Original_communalities', 'Mean_B','CI_low_B','CI_high_B')
    rownames(results_c)<-names(original_communalities)

    results_list[['communalities']]<-results_c
  }

  if (barmap_plot){
    results_list[['boot_barmap_loadings']]<-
      barmap_loading(pca, pca_data, load_list =  b_list,
                     star_values = F,
                     text_values = F,ndim = ndim, plot_cutoff = F,
                     plot_title = paste('Bootstrap s. loadings ','(',b_pca$R,' samples)',sep = ''),...)

    if (communalities){
      com_df<-as.data.frame(results_c)
      com_df$Variables<-row.names(results_c)

      com_plot<-barmap_commun(pca, pca_data, ndim=ndim,
                              plot_title = paste('Bootstrap communalities ','(',b_pca$R,' samples, ', ndim, ' PCs)',sep = ''),...)+
        geom_point(data=com_df, aes(y=.data$Mean_B, x=.data$Variables), inherit.aes = F)+
        geom_errorbar(data=com_df, aes(ymin=.data$CI_low_B, ymax=.data$CI_high_B,
                                       x=.data$Variables), inherit.aes = F, width=0.5)

      results_list[['boot_barmap_communalities']]<-com_plot
    }
  }
  cat('\n','=== DONE ===', '\n')
  cat('Final B iterations:', dim(b_pca$t)[1],'\n')
  if(B-dim(b_pca$t)[1]!=0){
    cat('=============================================','\n')
    cat(B-dim(b_pca$t)[1],'iterations could not be computed. Low number of rows might cause some bootstrap samples not suitable for PCA','\n')
  }
  return(results_list)
}

#'@title permutation test of PCA
#'
#'@description Compute a nonparametric permutation test for a PCA solution. Two options are
#'possible: hypothesis testing of the total variance accounted for (VAF) for each component as in XXX ref, and
#'hypothesis testing of the standardized loadings as in XXX ref.
#'
#'@author Abel Torres Espin
#'
#'@param pca Object of class \emph{prcomp}, \emph{princals}, or \emph{data.frame}. If object is a \emph{prcomp} or \emph{princals} object, \emph{pca_data} is required, and
#'the loadings will be extracted. If object is a data.frame object, the dataframe needs to be formatted as: first column named \emph{Variables} and all other columns corresponding to a PC.
#' One row per variable. The values are the loadings.
#'
#'@param pca_data Data passed to the \emph{prcomp} or \emph{princals} function.
#'@param ndim Numeric. Number of PCs (1 to \emph{ndim}) to run the analysis on. D
#'@param P Numeric. Number of permutations to run calling the \emph{permuted_pca} function. Default=1000
#'@param statistic Character. Determines the statistic to compute. Possible values are
#'Variance accounted for ("VAF") or the standardized loadings ("s.loadings"). Default="VAF"
#'@param conf Numeric. Level of confidence region for the confidence interval. E.g. 0.95 generates 95CI. Default=0.95
#'@param adj.method Character passed to the \emph{stats::p.adjust()} to adjust the p value for multiple comparisons.
#'See \emph{?p.adjust.methods}. Default="BH"
#'@param perm.method Character determining the permutation method to use as in Linting et al., 2011.
#'"permD" (Buja & Eyuboglu, 1992; Linting et al., 2011) where variables as permuted independently and
#'concomitantly (Fig. 2A) as opposite of the "permV" permutation strategy (Linting et al., 2011) where
#'variables are permuted one at the time. If statistic is set to "VAF", \emph{perm.method} "permD" will be
#'used. Default="permV".
#'@param plot Boolean. Whether to retrun a ggplot2 object containing a plot of "VAF", "loadings" or "communalities" as specified in \emph{statistic}.
#'@param ... Other arguments passed to the \emph{barmap_loadings} or \emph{barmap_commun} function when \emph{barmap_plot}=TRUE

#'@details Nonparametric permutation for hypothesis testing of the VAF of component, the loadings or communalities
#'have been studied (see refs). The hypothesis test is defined as:
#'H(null): PC metric (either VAF or loading) is indistinguishable from a random generation
#'H(alternative): PC metric (either VAF or loading) is different from random
#'The null distribution is generated by permuting the values of each variable several times (P)
#'and re-running the PCA on each permuted sample. Confidence intervals of the permuted distribution (null distribution) are
#'calculated using the percentile method. The p values are calculated as p = ((q+1))⁄((P+1)), where q is the number of times
#'the chosen metric is higher in the permuted distribution than in the original PCA solution and P is the number of permutations.
#'The user should note that the lowest p value that can be calculated is dependent on P. As an example, if P is set to a value of 10 (a relatively low value),
#'the smallest p value that can be detected is 0.09, considering q=0. Accordingly, P should be set high enough to reach the desired floor p value.
#'By default, we have set the number of permutations to 1000 (smallest p value approximately equal to 0.001 as a result) as this has been shown to be
#'high enough for approximating the null distribution in most cases.
#'
#'Permutation test of the loadings as in (Buja & Eyuboglu, 1992; Peres-Neto et al., 2003) that can serve to determine
#'the loading threshold, where the variables are permuted simultaneously and concomitantly. Linting et al., designed and tested an strategy
#'where only one variable is permuted at the time, showing great results in determining the contribution of variables using communalities (Linting et al., 2011).
#'This method has resulted in better determination of the significant contribution of variables on the PCA solution with higher statistical power and proper
#'type I error, and therefore has been incorporated in the package as the suggested method for loadings and communalities.
#'Following Linting et al., terminology, user can specify the permutation strategy for the loadings as one variable at the time (\emph{permV}, as in Linting et al., 2011)
#'or as all the variable together (\emph{permD}, as in Buja & Eyuboglu, 1992; Peres-Neto et al., 2003).
#'
#'@return List containing the following objects
#'\describe{
#'   \item{perm.methods}{Character for the used permutation method}
#'   \item{per_list}{List of P loading matrices for the P permuted samples}
#'   \item{statistic}{Character for the used statistic}
#'   \item{adj.method}{Character for the used method for adjusting p values}
#'   \item{per_plot}{ggplot2 object with the plot of VAF, loadings or communalities depending on the specified statistic}
#'   \item{results}{List of data.frame(s) with the results of the permutation test}
#' }
#'
#'@references
#'\enumerate{
#'  \item Buja A, Eyuboglu N. Remarks on Parallel Analysis. Multivar Behav Res. 1992 Oct 1;27(4):509–40
#'  \item Linting M, van Os BJ, Meulman JJ. Statistical Significance of the Contribution of Variables to the PCA solution: An Alternative Permutation Strategy. Psychometrika. 2011 Jul 1;76(3):440–60
#'}
#'
#'@export
#'
#'@import dplyr
#'@importFrom stats quantile
#'@importFrom stats p.adjust
#'

permut_pc_test<-function(pca, pca_data, P=1000, ndim=3, statistic='VAF', conf=0.95,
                         adj.method='BH', perm.method='permV', plot=T,...){

  # P=100
  # ndim=5
  # statistic='s.loadings'
  # conf=0.95
  # adj.method='BH'
  # perm.method='permD'

  results_list<-list()

  if (statistic=='VAF' | perm.method=='permD'){
    pb <- dplyr::progress_estimated(P+1)
    cat(paste("Permuting ", P, " times... for ", statistic, " using permD method", sep=''), '\n')
    results_list[['perm.method']]<-'permD'
  }else if ((statistic=='s.loadings'| statistic=='commun') & perm.method=='permV'){
    pb <- dplyr::progress_estimated((P+1)*dim(pca_data)[2])
    cat(paste("Permuting ", P, " times... for ", statistic, " using permV method", sep=''), '\n')
    results_list[['perm.method']]<-'permV'
  }

  per_list<-permut_pca(pca,pca_data, P,ndim = ndim, output = statistic, pb = pb, perm.method = perm.method)

  pb$i<-pb$n
  print(pb)

  results_list[['per_list']]<-per_list
  results_list[['statistic']]<-statistic
  results_list[['adj.method']]<-adj.method


  if (statistic=='VAF'){
    cat('\n','Calculating VAF...','\n')

    df_per<-do.call(rbind, per_list)

    if (class(pca)[1]=='prcomp'){
      original_VAF<-pca$sdev^2/sum(pca$sdev^2)
      ndim<-min(dim(pca$rotation)[2], ndim)
    }else if (class(pca)[1]%in%c('princals')){
      original_VAF<-pca$evals/sum(pca$evals)
      ndim<-min(pca$ndim, ndim)
    }

    mean_VAF<-colMeans(df_per)

    ci_low<-apply(df_per,2,function(x) quantile(x, (1-conf)/2, na.rm = T))
    ci_high<-apply(df_per,2,function(x) quantile(x, 1-(1-conf)/2, na.rm = T))

    pvalue<-sapply(1:ndim, function(x){
      (sum(df_per[,x]>original_VAF[x])+1)/(length(per_list)+1)#as in Buja & Eyuboglu, 1992;
    })

    if (adj.method=='none'){
      results<-cbind(original_VAF[1:ndim],mean_VAF[1:ndim],ci_low[1:ndim],ci_high[1:ndim],pvalue)
      colnames(results)<-c('Original_VAF', 'Mean_P', 'CI_low_P', 'CI_high_P','p_value')
      rownames(results)<-colnames(pca$rotation)[1:ndim]
    }else{
      adj.p.value<-p.adjust(pvalue, method=adj.method)
      results<-cbind(original_VAF[1:ndim],mean_VAF[1:ndim],ci_low[1:ndim],ci_high[1:ndim],pvalue,adj.p.value)
      colnames(results)<-c('Original_VAF', 'Mean_P', 'CI_low_P', 'CI_high_P','p_value','adj_p_value')
      rownames(results)<-colnames(pca$rotation)[1:ndim]
    }

    if (plot){
      plot_df<-as.data.frame(results)
      plot_df$pc<-row.names(plot_df)
      plot_df<-plot_df%>%select(.data$pc, .data$Original_VAF, .data$Mean_P, .data$adj_p_value)%>%pivot_longer(-c(.data$pc, .data$adj_p_value))%>%
        mutate(name=factor(.data$name, levels=c('Original_VAF', 'Mean_P')))
      per_CI<-as.data.frame(results)
      per_CI$pc<-colnames(pca$rotation)[1:ndim]
      per_CI<-per_CI%>%
        mutate(ymin=.data$CI_low_P,ymax=.data$CI_high_P)%>%
        select(.data$pc, .data$ymin, .data$ymax)

      per_plot<-ggplot(plot_df, aes(.data$pc,.data$value*100, color=.data$name, group=.data$name))+
        geom_point()+
        geom_line()+
        geom_errorbar(data=per_CI, aes(x=.data$pc, ymin=.data$ymin*100, ymax=.data$ymax*100), inherit.aes = F,
                      position = position_dodge(), color='orange', width=0.4)+
        theme_minimal()+
        scale_color_manual(values=c('steelblue','orange'),
                           name = "VAF", labels = c("Original", "Permuted"))+
        xlab(NULL)+ylab('Variance acounted for (VAF %)')+
        ggtitle(paste('VAF permutation test ','(',P,' samples)',sep = ''))+
        theme(legend.position = c(0.7,0.8), legend.background = element_rect(color='grey'))

      results_list[['per_plot']]<-per_plot
    }

  }

  if (statistic=='s.loadings'){
    cat('\n','Calculating loadings...','\n')

    original_loadings<-stand_loadings(pca, pca_data)
    nvars<-dim(pca_data)[2]
    var_names<-rownames(pca_data)

    array_per<-array(unlist(per_list), dim=c(dim(original_loadings)[1],
                                               dim(original_loadings)[2],P))

    mean_loadings<-apply(array_per, 1:2, mean)
    ci_low<-apply(array_per,1:2,function(x) quantile(x, (1-conf)/2, na.rm = T))
    ci_high<-apply(array_per,1:2,function(x) quantile(x, 1-(1-conf)/2, na.rm = T))

    pvalue<-sapply(1:ndim, function(x){
      sapply(1:nvars, function(z){
        (sum(abs(array_per[z,x,])>abs(original_loadings[z,x]))+1)/(length(per_list)+1)#as in Buja & Eyuboglu, 1992;
      })
    })

    results<-list()
    adj.p.value<-matrix(p.adjust(c(pvalue), method=adj.method), ncol = ndim)
    for (i in 1:ndim){
      if (adj.method=='none'){
        results_matrix<-cbind(original_loadings[,i],mean_loadings[,i], ci_low[,i],ci_high[,i],pvalue[,i])
        colnames(results_matrix)<-c('Original_loading', 'Mean_P','CI_low_P','CI_high_P','p_value')
        rownames(results_matrix)<-row.names(original_loadings)
        results[[colnames(original_loadings)[i]]]<-results_matrix
      }else{
        results_matrix<-cbind(original_loadings[,i],mean_loadings[,i], ci_low[,i],ci_high[,i],pvalue[,i],adj.p.value[,i])
        colnames(results_matrix)<-c('Original_loading', 'Mean_P','CI_low_P','CI_high_P','p_value', 'adj_p_value')
        rownames(results_matrix)<-row.names(original_loadings)
        results[[colnames(original_loadings)[i]]]<-results_matrix
      }
    }

    if (plot){
        per_plot<-
          barmap_loading(pca, pca_data, load_list = per_list,cutoff =0,
                         star_values = F,plot_list_original = T,plot_list_center = T,
                         text_values = F,ndim = ndim, plot_cutoff = F,
                         plot_title = paste('S. loadings permutation test ',
                                            '(',P,' samples)',sep = ''),...)

        results_list[['per_plot']]<-per_plot
      }
    }

    if (statistic=='commun'){
      cat('\n','Calculating communalities...','\n')

      original_loadings<-stand_loadings(pca, pca_data)[,1:ndim]
      original_communalities<-rowSums(original_loadings^2)
      nvars<-dim(pca_data)[2]
      var_names<-rownames(pca_data)

      if (ndim>1){
        per_communalities<-lapply(per_list, function(x){
          rowSums(x[,1:ndim]^2)
        })
      }else{
        per_communalities<-lapply(per_list, function(x){
          x[,1:ndim]^2
        })
      }

      array_per<-do.call(rbind, per_communalities)

      mean_communalities<-apply(array_per, 2, mean)
      ci_low<-apply(array_per,2,function(x) quantile(x, (1-conf)/2, na.rm = T))
      ci_high<-apply(array_per,2,function(x) quantile(x, 1-(1-conf)/2, na.rm = T))

      pvalue<-sapply(1:nvars, function(x){
          (sum(abs(array_per[,x])>abs(original_communalities[x]))+1)/(length(per_list)+1)#as in Buja & Eyuboglu, 1992;
      })

      adj.p.value<-p.adjust(c(pvalue), method=adj.method)

      if (adj.method=='none'){
        results<-cbind(original_communalities,mean_communalities, ci_low,ci_high,pvalue)
        colnames(results)<-c('Original_communalities', 'Mean_P','CI_low_P','CI_high_P','p_value')
        rownames(results)<-names(original_communalities)
      }else{
        results<-cbind(original_communalities,mean_communalities, ci_low,ci_high,pvalue,adj.p.value)
        colnames(results)<-c('Original_communalities', 'Mean_P','CI_low_P','CI_high_P','p_value','adj_p_value')
        rownames(results)<-names(original_communalities)
      }

      if (plot){
        per_df<-as.data.frame(results)
        per_df$Variables<-row.names(results)

        per_plot<-barmap_commun(pca, pca_data, ndim=ndim,
                                plot_title = paste('Communality permutation test ','(',P,' samples, ', ndim, ' PCs)',sep = ''),...)+
          geom_point(data=per_df, aes(y=.data$Mean_P, x=.data$Variables), inherit.aes = F)+
          geom_errorbar(data=per_df, aes(ymin=.data$CI_low_P, ymax=.data$CI_high_P,
                                         x=.data$Variables), inherit.aes = F, width=0.5)

        results_list[['per_plot']]<-per_plot
      }
    }

  results_list[['results']]<-results

  cat('\n','DONE', '\n')
  return(results_list)
}

#'
# pc_stability2<-function(pca, pca_data, ndim=3, B=1000, sim='ordinary',
#                         test_similarity=T, similarity_metric='all', s_cut_off=0.1,
#                         ci_type='bca', conf=0.95,barmap_plot=T,inParallel=F,...){
#
#   # pca<-prcomp(mtcars[1:10,], center = T, scale. = T)
#   # pca_data<-mtcars[1:10,]
#   # ndim=3
#   # B=200
#   # sim='ordinary'
#   # ci_type='bca'
#   # conf=0.95
#   # test_similarity=T
#   # similarity_metric='all'
#   # s_cut_off=0.4
#
#   cat('\n',paste("Bootstrapping ", B, " times...", sep=''), '\n')
#
#   if (inParallel){
#     m<-detectCores()-1
#     cl <- makeCluster(m)
#     clusterEvalQ(cl,library(syndRomics))
#
#   }else{
#     pb <- startpb(0, B)
#     pb_count<<-0
#   }
#   nvars<-dim(pca_data)[2]
#   var_names<-colnames(pca_data)
#   original_loadings<-stand_loadings(pca, pca_data)
#   ndim<-min(dim(original_loadings)[2], ndim)
#
#   if (inParallel){
#     b_pca<-boot::boot(data = pca_data, pca=pca, statistic = boot_pca_sample, R=B,
#                       ndim=ndim,original_loadings=original_loadings,
#                       sim = sim,parallel = 'snow', ncpus=m, cl=cl)
#   }else{
#     b_pca<-boot::boot(data = pca_data, pca=pca, statistic = boot_pca_sample, R=B,
#                       ndim=ndim,original_loadings=original_loadings,
#                       sim = sim, pb=pb)
#   }
#
#   b_pca$t<-na.omit(b_pca$t)
#   b_pca$R<-dim(b_pca$t)[1]
#   rm("pb_count")
#
#   cat('\n',"Calculating confident intervals...",'\n')
#
#   if (inParallel){
#     try_ci<-try(
#       ci_results<-pblapply(cl=cl,X=1:dim(b_pca$t)[2], FUN=function(i){
#         ci<-boot::boot.ci(b_pca, index = i, type = ci_type,conf = conf,...)
#         ci[[names(ci)[str_detect(names(ci),ci_type)]]][,4:5]
#       })
#     )
#   }else{
#     pb <- startpb(0, dim(b_pca$t)[2])
#     try_ci<-try(
#       ci_results<-lapply(1:dim(b_pca$t)[2],function(i){
#         setpb(pb,i)
#         ci<-boot::boot.ci(b_pca, index = i, type = ci_type,conf = conf,...)
#         ci[[names(ci)[str_detect(names(ci),ci_type)]]][,4:5]
#       })
#     )
#   }
#
#   if (class(try_ci)=='try-error'){
#     stop('Computation of confident intervals has failed, probably because B
#          is too small. Please increase B and try again.')
#   }
#
#   ci_results<-do.call(rbind,ci_results)
#
#
#   ci_low<-matrix(c(ci_results[,1]),ncol = ndim)
#   row.names(ci_low)<-var_names
#   colnames(ci_low)<-paste('PC',1:ndim, sep = '')
#   ci_high<-matrix(c(ci_results[,2]),ncol = ndim)
#   row.names(ci_high)<-var_names
#   colnames(ci_high)<-paste('PC',1:ndim, sep = '')
#
#   array_b<-array(t(b_pca$t), dim=c(dim(original_loadings)[1],
#                                    ndim,b_pca$R))
#
#   b_list<-lapply(1:dim(array_b)[3], function(x){
#     temp<-array_b[ , 1:ndim, x]
#     row.names(temp)<-var_names
#     temp
#   })
#
#   boot_mean<-apply(array_b, 1:2, mean)
#   row.names(boot_mean)<-var_names
#   colnames(boot_mean)<-paste('PC',1:ndim, sep = '')
#
#   if (test_similarity){
#     cat('\n',"Calculating PC similarities...",'\n')
#
#     if(inParallel){
#       similarity_res<-pblapply(cl=cl, X=b_list, FUN=function(i){
#         load.list<-list(original_loadings[,1:ndim],i[,1:ndim])
#         component_similarity(load.list, ndim=ndim, s_cut_off = s_cut_off,
#                              similarity_metric = similarity_metric)$index_mean
#       })
#
#       stopCluster(cl)
#     }else{
#       pb <- startpb(0, length(b_list))
#       similarity_res<-lapply(1:length(b_list),function(i){
#         setpb(pb,i)
#         load.list<-list(original_loadings[,1:ndim],b_list[[i]][,1:ndim])
#         component_similarity(load.list, ndim=ndim, s_cut_off = s_cut_off,
#                              similarity_metric = similarity_metric)$index_mean
#       })
#     }
#
#     similarity_res<-do.call(rbind, similarity_res)
#
#     ci_sim_low<-similarity_res%>%group_by(.data$PC)%>%
#       summarise_all(.funs=function(x) quantile(x, (1-conf)/2, na.rm = T))
#     ci_sim_high<-similarity_res%>%group_by(.data$PC)%>%
#       summarise_all(.funs=function(x) quantile(x, 1-(1-conf)/2, na.rm = T))
#     sim_mean<-similarity_res%>%group_by(.data$PC)%>%
#       summarise_all(mean)
#
#     return_list<-list('boot_samples'=b_list, 'boot_mean'=boot_mean,
#                       'ci_low'=ci_low,'ci_high'=ci_high,
#                       'PC_similarity'=list('similarity_mean'=sim_mean,
#                                            'similarity_ci_low'=ci_sim_low,
#                                            'similarity_ci_high'=ci_sim_high))
#   }else{
#     return_list<-list('boot_samples'=b_list, 'boot_mean'=boot_mean,
#                       'ci_low'=ci_low,'ci_high'=ci_high)
#   }
#
#   if (barmap_plot){
#     b_load<-as.data.frame(boot_mean)
#     b_load$Variables<-row.names(b_load)
#     b_load2<-b_load%>%pivot_longer(-(.data$Variables))
#
#     b_load_low<-as.data.frame(ci_low)
#     b_load_low$Variables<-row.names(b_load_low)
#     b_load_low<-b_load_low%>%pivot_longer(-(.data$Variables))
#     colnames(b_load_low)[3]<-'ymin'
#
#     b_load_high<-as.data.frame(ci_high)
#     b_load_high$Variables<-row.names(b_load_high)
#     b_load_high<-b_load_high%>%pivot_longer(-(.data$Variables))
#     colnames(b_load_high)[3]<-'ymax'
#
#     b_load_ci<-b_load_low
#     b_load_ci$ymax<-b_load_high$ymax
#     b_load_ci$loading<-b_load2$value
#     colnames(b_load_ci)[2]<-'component'
#
#     return_list[['boot_barmap']]<-
#       barmap_loading(b_load, cutoff =0, star_values = F,
#                      text_values = F,ndim = ndim, plot_cutoff = F,
#                      plot_title = paste('Bootstrap s. loadings ','(',b_pca$R,' samples)',sep = ''),...)+
#       geom_errorbar(data=b_load_ci, aes(ymin=.data$ymin, ymax=.data$ymax), width=0.5)+
#       labs(color=NULL)
#   }
#   cat('\n','=== DONE ===', '\n')
#   cat('Final B iterations:', dim(b_pca$t)[1],'\n')
#   if(B-dim(b_pca$t)[1]!=0){
#     cat('=============================================','\n')
#     cat(B-dim(b_pca$t)[1],'iterations could not be computed. Low number of rows might cause some bootstrap samples not suitable for PCA','\n')
#   }
#   return(return_list)
# }
#
