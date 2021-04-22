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
#'the same cutoff will be used for all the PCs. If a vector is passed, each value will be used for the corresponding PC.
#'@param VAF If \emph{pca} is from \emph{prcomp} or \emph{princals}, VAF argument is not needed. Otherwise, VAF is a String vector with text for the centers of the syndromic plots. The text generally corresponds to the variance accounted for (VAF) for the respective PCs.
#'@param arrow_size_multi Numeric. Controls the size of the arrows proportional to the loading. Default=10
#'@param repel Logical. Whether to repel the text for preventing text overlap. Default= TRUE
#'@param plot_legend Logical. Whether to plot the legend or not. Default= TRUE
#'@param text_size Numeric. Controls for the size of the text. Default=9
#'@param var_order Character. Specify the order of the variables in the plot by the loading values, starting at 12 o’clock and moving counterclockwise. Possible values: 'abs decreasing': plot by decreasing absolute value;
#''abs increasing': plot by increasing absolute value; 'decreasing'; or 'increasing’.
#'@param colors Character vector of length 3. Vector with the character name or
#'hexadecimal number (e.g. "#FF0000") of three colors, the lower color, the middle color and the higher color
#'for the gradient used in the plot. Hexadecimal number can be obtained using \emph{rgb} for example.
#'
#'@param ... Other arguments passed to \emph{extract_syndromic_plot()}
#'
#'@return Returns a list of \emph{ggplot2} objects with one element for each PC plot. Plots can be saved usign ggsave().
#'
#'@examples
#'data(mtcars)
#'pca_mtcars<-prcomp(mtcars, center = TRUE, scale. = TRUE)
#'
#'syndromic_plot(pca = pca_mtcars, pca_data = mtcars, ndim = 2, cutoff = 0.5)
#'
#'@export
#'
#'@import dplyr tidyr ggnewscale ggplot2
#'@importFrom rlang .data
#'
syndromic_plot<-function(pca, pca_data=NULL, ndim=3, cutoff, VAF,arrow_size_multi=10,
                         repel= TRUE, plot_legend= TRUE, text_size=9,
                         var_order='abs decreasing',colors=c("steelblue1","white","firebrick1"),
                         ...){

  load_df<-extract_loadings(pca, pca_data)

  if (class(pca)[1]=='prcomp'){
    VAF<-paste(round(pca$sdev^2/sum(pca$sdev^2)*100,1), '%', sep = '')
  }else if (class(pca)[1]%in%c('princals')){
    VAF<-paste(round(pca$evals/sum(pca$evals)*100,1), '%', sep='')
  }

  if (min(dim(load_df)[2]-1<ndim)){
    ndim<-min(dim(load_df)[2]-1, ndim)
    warning(sprintf("The specified ndim is bigger than number of dimensions, using ndim = %s", ndim))
  }

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

    try(s_plot<-extract_syndromic_plot(load_df = load_df, pc=pc, cutoff=c, VAF=v, text_size= text_size,
                                       arrow_size_multi = arrow_size_multi, repel = repel,
                                       var_order=var_order, plot_legend = plot_legend,
                                       colors=colors))

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
#'Value of the loadings threshold (i.e. |loadings| >= cutoff) to plot with stars. Default = 0
#'@param legend_title String. Title of the legend. 's. loading' by default.
#'@param text_values Logical. Whether to plot the values of the loadings or not. Default= TRUE
#'@param star_values Logical. Whether to plot a star in |loadings|>=cutoff. Only relevant if
#'\emph{text_values}=FALSE. Default=FALSE
#'@param text_size Numeric. Size of the text_values.
#'@param vars Character vector. Variables will be ordered as the provided variable names. Non-specified
#'variables will be excluded from the plot. By default variables are ordered in alphabetically by ggplot.
#'@param colors Character vector of length 3. Vector with the character name or
#'hexadecimal number (e.g. "#FF0000") of three colors, the lower color, the middle color and the higher color
#'for the gradient used in the plot. Hexadecimal number can be obtained using \emph{rgb} for example.
#'
#'@return Returns a \emph{ggplot2} object.
#'
#'@examples
#'data(mtcars)
#'pca_mtcars<-prcomp(mtcars, center = TRUE, scale. = TRUE)
#'
#'heatmap_loading(pca = pca_mtcars, pca_data = mtcars, ndim = 1:4)
#'
#'@export
#'
#'@import ggplot2 dplyr tidyr ggnewscale
#'@importFrom rlang .data
#'
heatmap_loading<-function(pca, pca_data, ndim=1:10, cutoff=0,
                          legend_title='s. loading', text_values= TRUE, star_values=F,
                          text_size=2, vars=NULL, colors=c("steelblue1","white","firebrick1")){

  old_scipen<-getOption('scipen')
  on.exit(options(scipen=old_scipen))

  options(scipen=999)

  load_df<-extract_loadings(pca, pca_data)

  if (max(ndim)>ncol(load_df[,-1])){
    stop('higher value for ndim can not be higher than the number of PCs')
  }


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

  # load_df[load_df$Variables%in%dif_var,'loading_txt']<-sapply(unique(load_df$component),function(x){
  #   paste("|",
  #         abs(round(load_df[load_df$Variables%in%dif_var & load_df$component==x,'loading'],3)),
  #         "|",sep = "")
  # })
#
#   arbitrary_df<-load_df%>%
#     filter(.data$component%in%paste('PC',1:ndim, sep = ''))%>%
#     mutate(weight=ifelse(.data$Variables%in%dif_var,abs(.data$loading),NA))

  # load_df<-load_df%>%
  #   mutate(loading=ifelse(.data$Variables%in%dif_var,NA,.data$loading))

  h_plot<-load_df%>%filter(.data$component%in%paste('PC',ndim, sep = ''))%>%
    ggplot(aes(.data$component, .data$Variables, fill=.data$loading))+
    geom_raster()+
    scale_fill_gradient2(low=colors[1], high=colors[3],mid =colors[2],
                         limits=c(-1,1), na.value =  "transparent")+
    labs(fill=legend_title)+
    ylab(NULL)+xlab(NULL)+
    theme_minimal()

  # if (!is.null(dif_var)){
  #
  #   h_plot<-h_plot+
  #     ggnewscale::new_scale_fill() +
  #     geom_raster(data=arbitrary_df,aes(fill=.data$weight))+
  #     scale_fill_gradient(low='white', high='grey30',limits=c(0,1), na.value =  "transparent")+
  #     labs(fill= paste("|",legend_title,"|", sep = ""))
  # }

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
#'@param resample_ci dataframe. Dataframe containing the columns "Variables", "component","original", "mean","ci_low" and "ci_high" containing the center, the lower bound and the upper bound of the confident intervals to plot.
#'These can be obtained by the pc_stability or permut_pc_test "results" element or computed independently.
#'@param conf Numeric. Confidence level used when \emph{load_list} is provided.
#'@param plot_list_center Logical. Whether to plot the average loadings obtained from \emph{load_list}.
#'@param plot_legend Logical. Whether legend should be plotted.
#'@param text_values Logical. Whether to plot the values of the loadings or not. Default= TRUE
#'@param star_values Logical. Whether to plot a star in |loadings|>=cutoff. Only relevant if
#'\emph{text_values}=FALSE. Default=FALSE
#'@param text_size Numeric. Size of the text_values.
#'@param plot_cutoff Logical. Whether to plot the cutoff lines or not.
#'@param vars Character vector. Variables will be ordered as the provided variable names. Non-specified
#'variables will be excluded from the plot. By default variables are ordered in alphabetically by ggplot.
#'@param colors Character vector of length 3. Vector with the character name or
#'hexadecimal number (e.g. "#FF0000") of three colors, the lower color, the middle color and the higher color
#'for the gradient used in the plot. Hexadecimal number can be obtained using \emph{rgb} for example.
#'@param gradient_color Logical. Whether colors should be plotted as a gradient proportional to the loadings.
#'
#'@return Returns a \emph{ggplot2} object.
#'
#'@examples
#'data(mtcars)
#'pca_mtcars<-prcomp(mtcars, center = TRUE, scale. = TRUE)
#'
#'barmap_loading(pca = pca_mtcars, pca_data = mtcars, ndim = 1:4)
#'
#'@export
#'
#'@import ggplot2 dplyr tidyr ggnewscale
#'@importFrom rlang .data
#'
barmap_loading<-function(pca, pca_data, ndim=1:5, cutoff=0,resample_ci=NULL,conf=0.95, plot_list_center=F,
                         plot_legend= TRUE,text_values=F, star_values=F,
                         text_size=2, plot_cutoff= TRUE, vars=NULL, colors=c("steelblue1","white","firebrick1"),
                         gradient_color= TRUE){

  # resample_ci<-b$boot_samples
  # pca_data<-mtcars
  # ndim=3
  # cutoff=0.5
  # conf=0.95
  # plot_title='Standardized loadings'
  # legend_title='s. loading'
  # text_values=F
  # star_values= TRUE
  # text_size=2
  # plot_cutoff= TRUE
  # vars=NULL
  # plot_legend= TRUE

  old_scipen<-getOption('scipen')
  on.exit(options(scipen=old_scipen))

  options(scipen=999)

  load_df<-extract_loadings(pca, pca_data)

  if (max(ndim)>ncol(load_df[,-1])){
    stop('higher value for ndim can not be higher than the number of PCs')
  }

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

  if (length(cutoff)!=1 & length(cutoff)!=length(ndim)){
    stop('cutoff length must be one or equal to ndim length')
  }

  cutoff_df<-load_df%>%pivot_longer(-(.data$Variables),
                                    names_to = "component", values_to = "loading")%>%
    filter(.data$component%in%paste('PC',ndim, sep = ''))%>%
    group_by(.data$component)%>%
    summarise(count=n())%>%mutate(cutoff=cutoff)%>%
    mutate(component=factor(.data$component, levels=names(load_df)))

    load_df<-load_df%>%pivot_longer(-(.data$Variables),
                                    names_to = "component", values_to = "loading")%>%
      filter(.data$component%in%paste('PC',ndim, sep = ''))%>%
      left_join(cutoff_df, by='component')%>%
      group_by(.data$component)%>%
      mutate(cutoff=cutoff,
             star=ifelse(abs(.data$loading)>=cutoff,'*',''),
             loading_txt=as.character(round(.data$loading,2)),
             weight=(round(.data$loading,2)))%>%
      ungroup()%>%
      mutate(component=factor(.data$component, levels=names(load_df)))

    b_plot<-load_df%>%
      ggplot(aes(.data$Variables, .data$loading))+
      geom_col(show.legend = plot_legend)+
      theme_minimal()+
      coord_flip()+
      facet_grid(~.data$component)

  #PCA list
  if (!is.null(resample_ci)){

      if (max(ndim)>length(unique(resample_ci$component))){
        stop(paste0('higher value for ndim can not be higher than the
                    number of PCs in resample_ci which is ',length(unique(resample_ci$component))))
      }

    resample_ci<-as.data.frame(resample_ci%>%
      filter(.data$component%in%unique(load_df$component)))

      b_plot<-b_plot+
        geom_errorbar(data=resample_ci, aes(x=.data$Variables,ymin=.data$ci_low, ymax=.data$ci_high),width=0.5, inherit.aes = F)

      if (plot_list_center){
        b_plot<-b_plot+
          geom_point(data=resample_ci, aes(y=.data$mean, x=.data$Variables),inherit.aes = F)
      }
    }

    if (text_values){
      b_plot<-b_plot+geom_text(aes(label=.data$loading_txt,y=.data$loading*1.1),size=2)
    }else if (star_values){
      b_plot<-b_plot+geom_text(aes(label=.data$star, y=.data$loading*1.1), color='black')
    }

    if (gradient_color){
      b_plot<-b_plot+aes(.data$Variables, .data$loading, fill=.data$loading)+
        scale_fill_gradient2(low=colors[1], high=colors[3],mid = colors[2],
                              limits=c(-1,1), na.value =  "transparent")
    }else{
      b_plot<-b_plot+aes(.data$Variables, .data$loading, fill=(.data$loading)>0)+
        scale_fill_manual(values=c(colors[1],colors[3]), na.value =  "transparent")+
        theme(legend.position = "none")
    }

    if(plot_cutoff){
      b_plot<-b_plot+
        geom_hline(data=cutoff_df,aes(yintercept = cutoff), color=colors[3], alpha=0.6)+
        geom_hline(data=cutoff_df,aes(yintercept = -cutoff), color=colors[1], alpha=0.6)
    }

  b_plot<-b_plot+ylim(-1.1,1.1)+
    labs(fill='s.loading')+
    geom_hline(yintercept = 0, color='black')+
    ylab(NULL)+xlab(NULL)

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
#'@param load_list List. List of loading matrices used to plot percentile confidence intervals around
#'the average. If NULL, no error bars are plotted. Default=NULL
#'@param conf Numeric. Confidence level used when \emph{load_list} is provided.
#'@param plot_original Logical. Whether to plot the communalities obtained from the values passed to \emph{pca} and \emph{pca_data}.
#'@param plot_list_center Logical. Whether to plot the average communalities obtained from \emph{load_list}.
#'@param plot_legend Logical. Whether legend should be plotted.
#'@param text_values Logical. Whether to plot the values of the communalities or not. Default=FALSE
#'@param text_size Numeric. Size of the text_values.
#'@param vars Character vector. Variables will be ordered as the provided variable names. Non-specified
#'variables will be excluded from the plot. By default variables are ordered in alphabetically by ggplot.
#'@param var_order Character. Specify the order of the variables in the plot by the communality values, starting at 12 o’clock and moving counterclockwise. Possible values: 'abs decreasing': plot by decreasing absolute value;
#''abs increasing': plot by increasing absolute value; 'decreasing'; or 'increasing’.
#'@param colors Character vector of length 2. Vector with the character name or
#'hexadecimal number (e.g. "#FF0000") of two colors, the lower color and the higher color
#'for the gradient used in the plot. Hexadecimal number can be obtained using \emph{rgb} for example.
#'@param gradient_color Logical. Whether colors should be plotted as a gradient proportional to the communalities.
#'
#'@return Returns a \emph{ggplot2} object.
#'
#'@examples
#'data(mtcars)
#'pca_mtcars<-prcomp(mtcars, center = TRUE, scale. = TRUE)
#'
#'barmap_commun(pca = pca_mtcars, pca_data = mtcars, ndim = 1:4)
#'
#'@export
#'
#'@import ggplot2 dplyr tidyr ggnewscale
#'@importFrom rlang .data
#'
barmap_commun<-function(pca, pca_data, ndim=1:5,load_list=NULL,conf=0.95, plot_original= TRUE,plot_list_center=F,
                        plot_legend= TRUE,text_values=F,
                        text_size=2, var_order='increasing',vars=NULL,
                        colors=c("white","firebrick1"), gradient_color= TRUE){

  # pca=prcomp(mtcars, scale. = T, center = T)
  # pca_data=mtcars
  # ndim = 1:3
  # load_list = ibooted_PCA$boot_samples
  # conf = 0.95
  # plot_legend = TRUE
  # text_values = T
  # text_size = 2
  # var_order = "increasing"
  # vars = NULL
  # colors = c("white", "firebrick1")
  # gradient_color = TRUE
  # plot_original = TRUE
  # plot_list_center = TRUE


  old_scipen<-getOption('scipen')
  on.exit(options(scipen=old_scipen))

  options(scipen=999)

  load_df<-stand_loadings(pca, pca_data)
  if (max(ndim)>ncol(load_df[,-1])){
    stop('higher value for ndim can not be higher than the number of PCs')
  }

  if (is.null(load_list)){
    plot_original== TRUE
  }

  if (length(ndim)>1){
    commun_df<-data.frame(commun=rowSums(load_df[,ndim]^2),
                          Variables=row.names(load_df))
  }else{
    commun_df<-data.frame(commun=load_df[,ndim]^2,
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
    ggplot(aes(.data$Variables, .data$commun))
    # scale_fill_gradient2(low=colors[1], high=colors[2],limits=c(0,1), na.value =  "transparent")+

  if (plot_original){
    if (gradient_color){
      c_plot<-c_plot+aes(.data$Variables, .data$commun, fill=.data$commun)+
        geom_col()+
        scale_fill_gradient2(low=colors[1], high=colors[2],limits=c(0,1), na.value =  "transparent")
    }else{
      c_plot<-c_plot+
        geom_col(fill=colors[1])
    }
  }


  if (!is.null(load_list)){

    if (length(ndim)>1){
      b_communalities<-lapply(load_list, function(x){
        rowSums(x[,ndim]^2)
      })
    }else{
      b_communalities<-lapply(load_list, function(x){
        x[,ndim]^2
      })
    }

    array_b<-do.call(rbind, b_communalities)

    mean_communalities<-apply(array_b, 2, mean)
    ci_low_c<-apply(array_b,2,function(x) quantile(x, (1-conf)/2, na.rm = TRUE))
    ci_high_c<-apply(array_b,2,function(x) quantile(x, 1-(1-conf)/2, na.rm = TRUE))
    summary_communalities<-data.frame(mean_communalities = mean_communalities,
                                      ci_low_c = ci_low_c,
                                      ci_high_c = ci_high_c,
                                      Variables = names (mean_communalities))

    results_c<-cbind(commun_df,mean_communalities, ci_low_c,ci_high_c)

    results_c<-commun_df%>%
      full_join(summary_communalities, by="Variables")

    c_plot<-c_plot+
        geom_errorbar(data=results_c, aes(ymin=.data$ci_low_c, ymax=.data$ci_high_c,
                                       x=.data$Variables), inherit.aes = F, width=0.5)
    if (plot_list_center){
      c_plot<-c_plot+
        geom_point(data=results_c, aes(y=.data$mean_communalities, x=.data$Variables), inherit.aes = F)
    }
  }

  c_plot<-c_plot+
    labs(fill="communa.")+
    ylab(NULL)+xlab(NULL)+
    theme_minimal()+
    coord_flip()

  if (!plot_legend){
   c_plot<-c_plot+theme(legend.position = "none")
  }

  if (text_values){
    c_plot<-c_plot+geom_text(aes(label=.data$commun,y=.data$commun*1.1),size=2)
  }

  return(c_plot)
}

#'@title Variance accounted for plot
#'@description Plot of VAF from a PCA solution given the first ndim PCs.
#'
#'@author Abel Torres Espin
#'
#'@param pca Object of class \emph{prcomp}, \emph{princals}, or \emph{data.frame}. If object is a \emph{prcomp} or \emph{princals} object, \emph{pca_data} is required, and
#'the loadings will be extracted. If object is a data.frame object, the dataframe needs to be formatted as: first column named \emph{Variables} and all other columns corresponding to a PC.
#' One row per variable. The values are the loadings.
#'
#'@param pca_data Data passed to the \emph{prcomp} or \emph{princals} function.
#'@param ndim Numeric. Number of PCs to plot
#'@param resample_ci dataframe. Dataframe containing the columns "original", "mean","ci_low" and "ci_high" containing the center, the lower bound and the upper bound of the confident intervals to plot. Each row contain the values of 1 PC in order (PC1 first row, PC2 second row, etc).
#'This can be obtained by the permut_pc_test "results" element with statistic="VAF" or computed independently.
#'@param style Character. There are two styles of VAF plots "line" by default or "reduced".
#'@param colors Character vector of length 2. Vector with the character name or
#'hexadecimal number (e.g. "#FF0000") of two colors, the lower color and the higher color
#'for the gradient used in the plot. Hexadecimal number can be obtained using \emph{rgb} for example.
#'
#'@return Returns a \emph{ggplot2} object.
#'
#'@examples
#'data(mtcars)
#'pca_mtcars<-prcomp(mtcars, center = TRUE, scale. = TRUE)
#'
#'VAF_plot(pca = pca_mtcars, pca_data = mtcars, ndim = 1:7)
#'
#'@export
#'
#'@import ggplot2 ggnewscale
#'@importFrom rlang .data
#'
VAF_plot<-function(pca, pca_data, ndim=1:5, resample_ci=NULL, style="line", colors=c("steelblue","orange")){

  if (class(pca)[1]=='prcomp'){
    original_VAF<-pca$sdev^2/sum(pca$sdev^2)
  }else if (class(pca)[1]%in%c('princals')){
    original_VAF<-pca$evals/sum(pca$evals)
  }

  plot_df<-data.frame("value"=original_VAF[ndim],
                      "component"=factor(paste0('PC',ndim), levels = paste0('PC',ndim)),
                      "name"="Original")

  if (style=="line"){

    vaf_plot<-ggplot(plot_df, aes(.data$component,.data$value*100))+
      geom_point(aes(color="Original"), show.legend = F)+
      geom_line(aes(group=.data$name, color="Original"), show.legend = F)+
      theme_minimal()+
      xlab(NULL)+ylab('Variance accounted for (VAF %)')+
      theme(legend.position = c(0.7,0.8), legend.background = element_rect(color='grey'))+
      scale_colour_manual(values=colors[1])

    if (!is.null(resample_ci)){
      resample_ci<-resample_ci%>%
        mutate(component=factor(rownames(resample_ci), levels = rownames(resample_ci)),
               name="Permuted")
      resample_ci<-resample_ci[ndim,]
      suppressMessages(vaf_plot<-vaf_plot+
        geom_point(data=resample_ci, aes(x=.data$component, y=.data$mean*100, color="Permuted"),
                   inherit.aes = F)+
        geom_line(data=resample_ci, aes(x=.data$component, y=.data$mean*100,group=.data$name,color="Permuted"),
                  inherit.aes = F, show.legend = TRUE)+
        geom_errorbar(data=resample_ci,
                      aes(x=.data$component, y=.data$mean*100, ymin=.data$ci_low*100,ymax=.data$ci_high*100,color="Permuted"),
                      inherit.aes = F,width=0.5)+
        scale_color_manual(values=c(colors[1],colors[2]))+
        theme(legend.title = element_blank()))

    }
  }
  if(style=="reduced"){
    plot_df$VAF_value<-
      factor(paste0(round(plot_df$value*100,1),"%"),
             levels=rev(paste0(round(plot_df$value*100,1),"%")))

    vaf_plot<-ggplot(plot_df, aes(.data$value, .data$VAF_value))+
      geom_segment(aes(xend=0, yend=.data$VAF_value), size=1.5, color=colors[1])+
      xlab(NULL)+ylab('VAF (%)')+
      scale_x_continuous(breaks = NULL)+
      theme_minimal()+
      theme(panel.grid = element_blank())

    if(!is.null(resample_ci)){
      resample_ci$component<-rownames(resample_ci)
      resample_ci$name<-"Permuted"
      resample_ci$VAF_value<-
        factor(paste0(round(plot_df$value*100,1),"%"),
               levels=rev(paste0(round(plot_df$value*100,1),"%")))
      vaf_plot<-vaf_plot+
        geom_errorbar(data = resample_ci, aes(xmin=.data$ci_low, xmax=.data$ci_high,
                                              x=.data$mean, y=.data$VAF_value),
                      width=0.4,inherit.aes = F, color=colors[2])+
        geom_point(data = resample_ci, aes(x=.data$mean, y=.data$VAF_value),
                   inherit.aes = F, color=colors[2])
    }
  }
  return(vaf_plot)
}


#'@title PC stability by nonparametric bootstrapping
#'
#'@description Extract nonparametric estimate of the standardized loadings and the confident
#'region by means of bootstrapping. This function uses the \emph{boot} function from
#'the \strong{boot} package.
#'
#'@author Abel Torres Espin
#'
#'@param pca Object of class \emph{prcomp}, \emph{princals}.
#'@param pca_data Data passed to the \emph{prcomp} or \emph{princals} function.
#'@param ndim Numeric. Number of PCs (1 to \emph{ndim}) to run the analysis on. Default = 3.
#'@param B Numeric. Number of bootstrapped samples passed to \emph{boot}.
#'@param sim Character. Determines the bootstrapping method for \emph{boot}.
#'From \emph{boot}: A character string indicating the type of simulation required.
#'Possible values are "ordinary" (the default), "parametric", "balanced",
#'"permutation", or "antithetic". Default="ordinary".
#'@param communalities Logical. Whether to compute and return communalities.
#'@param similarity_metric character or character vector. Specify the similarity metric
#'to use. Set "none" to not return similarity metrics. Possible values are "cc_index" (congruence coefficient),
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
#'@param ... Other arguments passed to the \emph{boot} function.
#'
#'@return Returns a list object of class "syndromics".
#'\describe{
#'   \item{\strong{method}}{Specify the method used to obtain the syndromics list object}
#'   \item{\strong{pca}}{Contains the object passed to the pca argument}
#'   \item{\strong{pca_data}}{Contains the object passed to the pca_data argument}
#'   \item{\strong{ndim}}{Value specified in the ndim argument}
#'   \item{\strong{ci_method}}{Method used to compute CIs}
#'   \item{\strong{conf}}{Confidence level used to compute CIs}
#'   \item{\strong{results}}{Object containing the results of the analysis}
#'   \item{\strong{boot_sample}}{A list of length B containing the resampled loadings.}
#'   \item{\strong{pc_similarity}}{A list of results when \emph{similarity_metric} is not "none".}
#'   \describe{
#'    \item{\strong{similarity_mean}}{A numerix matrix with the mean of the chosen similarity metric.}
#'    \item{\strong{similarity_ci_low and similarity_ci_high}}{A numeric matrix with the lower and upper CI respectively.}
#'    }
#'   \item{\strong{B}}{Number of resamples computed.}
#'   \item{\strong{communalities}}{If communalities= TRUE, the results are returned here.}
#'}
#'
#'@details The number of bootstrap samples is set to 1000 by default, as it has been shown to be a robust number in most conditions of
#'data complexity and sample size. The user must be careful on setting such number too low which would reduce the performance
#'of the approximation. However, values that are too high might unnecessarily increase computing time with little gain (REFs).
#'
#'@references Efron B. Better Bootstrap Confidence Intervals. J Am Stat Assoc. 1987 Mar 1;82(397):171–85.
#'
#'@examples
#'data(mtcars)
#'pca_mtcars<-prcomp(mtcars, center = TRUE, scale. = TRUE)
#'
#'pca_mtcars_stab<-pc_stability(pca = pca_mtcars, pca_data = mtcars, ndim = 3, B = 500)
#'plot(pca_mtcars_stab, plot_resample= TRUE)
#'
#'@export
#'
#'@import ggplot2 dplyr boot stringr progress
#'@importFrom rlang .data
#'@importFrom stats na.omit
#'
pc_stability<-function(pca, pca_data, ndim=3, B=1000, sim='ordinary',communalities= TRUE,
                       similarity_metric='all', s_cut_off=0.1,
                       ci_type='bca', conf=0.95,...){

  # pca<-pca
  # pca_data<-data
  # ndim=3
  # B=10
  # sim='ordinary'
  # ci_type='bca'
  # conf=0.95
  # test_similarity= TRUE
  # similarity_metric='all'
  # s_cut_off=0.1

  results_list<-list()
  results_list[["method"]]<-'bootstrap'
  results_list[["pca"]]<-pca
  results_list[["pca_data"]]<-pca_data
  results_list[["ndim"]]<-ndim
  results_list[["ci_method"]]<-ci_type
  results_list[["conf"]]<-conf

  similarity_metric<-match.arg(similarity_metric,
                               c("all","cc_index", "r_correlation","rmse",
                                 "s_index","none"), several.ok = TRUE)

  nvars<-dim(pca_data)[2]
  var_names<-colnames(pca_data)
  original_loadings<-stand_loadings(pca, pca_data)
  ndim<-min(dim(original_loadings)[2], ndim)

  message("Bootstrapping ", B, " times")
  pb <- progress_bar$new(
    format = "[:bar] :percent ~remaining :eta",
    total = B+1, clear = FALSE, width= 80)

  if(class(pca)[1]=='prcomp'){
    center=F
    .scale=F
    if (is.numeric(pca$center)) center= TRUE
    if (is.numeric(pca$scale)) .scale= TRUE
  }

  b_pca<-boot::boot(data = pca_data, pca=pca, statistic = boot_pca_sample, R=B,
              ndim=ndim,original_loadings=original_loadings,
              sim = sim, pb=NULL, center=center, .scale=.scale,...)
  b_pca$t<-na.omit(b_pca$t)
  b_pca$R<-dim(b_pca$t)[1]

  message("Calculating confident intervals by ", ci_type," method")
  pb <- progress_bar$new(
    format = "[:bar] :percent ~remaining :eta",
    total = dim(b_pca$t)[2], clear = FALSE, width= 80)
  ci_results<-list()

  try_ci<-try(for(i in 1:dim(b_pca$t)[2]){
    pb$tick()
    ci<-boot::boot.ci(b_pca, index = i, type = ci_type,conf = conf,...)
    ci_results[[i]]<-ci[[names(ci)[str_detect(names(ci),ci_type)]]][,4:5]
  })

  if (class(try_ci)=='try-error'){
    message("\n",'Computation of confident intervals using ',ci_type,' has failed, probably because B
         is too small. If ', ci_type, ' is the confident interval of interest, increase B and try again.')
    message("\n","Calculating confident intervals by percentage method")

    pb <- progress_bar$new(
      format = "[:bar] :percent ~remaining :eta",
      total = dim(b_pca$t)[2], clear = FALSE, width= 80)
    ci_results<-list()

    for(i in 1:dim(b_pca$t)[2]){
      pb$tick()
      ci<-boot::boot.ci(b_pca, index = i, type = "perc",conf = conf,...)
      ci_results[[i]]<-ci[[names(ci)[str_detect(names(ci),"perc")]]][,4:5]
    }

    results_list[["ci_method"]]<-"perc"
  }

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

  results<-list()
  for (i in 1:ndim){
      results_matrix<-cbind(row.names(original_loadings),original_loadings[,i],boot_mean[,i], ci_low[,i],ci_high[,i])
      colnames(results_matrix)<-c("Variables",'Original_loading', 'mean','ci_low','ci_high')
      rownames(results_matrix)<-row.names(original_loadings)
      results[[colnames(original_loadings)[i]]]<-as_tibble(results_matrix)
  }

  results<-bind_rows(results, .id="component")%>%
    relocate(.data$Variables)%>%
    mutate_at(.vars = c('Original_loading','mean','ci_low','ci_high'),
              .funs = as.numeric)
  rownames(results)<-NULL

  results_list[['results']]<-results
  results_list[['boot_samples']]<-b_list
  results_list[['B']]<-B

  if (!"none"%in%similarity_metric){
    message("Calculating PC similarities")
    pb <- progress_bar$new(
      format = "[:bar] :percent ~remaining :eta",
      total = length(b_list), clear = FALSE, width= 80)
    similarity_res<-list()

    for(i in 1:length(b_list)){
      pb$tick()
      load.list<-list(original_loadings[,1:ndim],b_list[[i]][,1:ndim])
      similarity_res[[i]]<-component_similarity(load.list, ndim=ndim, s_cut_off = s_cut_off,
                                                similarity_metric = similarity_metric)$index_mean
    }

    similarity_res<-do.call(rbind, similarity_res)

    ci_sim_low<-similarity_res%>%group_by(.data$PC)%>%
      summarise_all(.funs=function(x) quantile(x, (1-conf)/2, na.rm = TRUE))
    ci_sim_high<-similarity_res%>%group_by(.data$PC)%>%
      summarise_all(.funs=function(x) quantile(x, 1-(1-conf)/2, na.rm = TRUE))
    sim_mean<-similarity_res%>%group_by(.data$PC)%>%
      summarise_all(mean)

    results_list[['PC_similarity']]<-list('similarity_mean'=sim_mean,
                                     'similarity_ci_low'=ci_sim_low,
                                     'similarity_ci_high'=ci_sim_high)
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
    ci_low_c<-apply(array_b,2,function(x) quantile(x, (1-conf)/2, na.rm = TRUE))
    ci_high_c<-apply(array_b,2,function(x) quantile(x, 1-(1-conf)/2, na.rm = TRUE))

    results_c<-cbind(original_communalities,mean_communalities, ci_low_c,ci_high_c)
    colnames(results_c)<-c('Original_communalities', 'mean','ci_low','ci_high')
    rownames(results_c)<-names(original_communalities)

    results_list[['communalities']]<-results_c
  }

  message("\n",'\nFinal B iterations: ', dim(b_pca$t)[1])

  if(B-dim(b_pca$t)[1]!=0){
    message(paste0(rep("=",60),collapse = ""))
    message(B-dim(b_pca$t)[1],' iterations could not be computed. Low number of rows might cause some bootstrap samples not suitable for PCA')
  }

  return(new_syndromics(results_list))
}

#'@title permutation test of PCA
#'
#'@description Compute a nonparametric permutation test for a PCA solution. Two options are
#'possible: hypothesis testing of the total variance accounted for (VAF) for each component, and
#'hypothesis testing of the standardized loadings and communalities.
#'
#'@author Abel Torres Espin
#'
#'@param pca Object of class \emph{prcomp}, \emph{princals}.
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
#'
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
#'@return List of class "syndromics" containing the following objects
#'\describe{
#'   \item{\strong{methods}}{Character for the used method to obtain the object. "permutation"}
#'   \item{\strong{statistic}}{Statistic use in the permutation analysis}
#'   \item{\strong{perm.methods}}{Character for the used permutation method}
#'   \item{\strong{ndim}}{Value specified in the ndim argument}
#'   \item{\strong{conf}}{Confidence level used to compute CIs}
#'   \item{\strong{per_sample}}{List of P loading matrices for the P permuted samples}
#'   \item{\strong{adj.method}}{Character for the used method for adjusting p values}
#'   \item{\strong{pca}}{Contains the object passed to the pca argument}
#'   \item{\strong{pca_data}}{Contains the object passed to the pca_data argument}
#'   \item{\strong{results}}{Object containing the results of the analysis}
#' }
#'
#'@references
#'\enumerate{
#'  \item Buja A, Eyuboglu N. Remarks on Parallel Analysis. Multivar Behav Res. 1992 Oct 1;27(4):509–40
#'  \item Linting M, van Os BJ, Meulman JJ. Statistical Significance of the Contribution of Variables to the PCA solution: An Alternative Permutation Strategy. Psychometrika. 2011 Jul 1;76(3):440–60
#'}
#'
#'@examples
#'data(mtcars)
#'pca_mtcars<-prcomp(mtcars, center = TRUE, scale. = TRUE)
#'
#'pca_mtcars_perm<-permut_pc_test(pca = pca_mtcars, pca_data = mtcars, ndim = 3, P = 500)
#'plot(pca_mtcars_perm, plot_resample= TRUE)
#'
#'@export
#'
#'@import dplyr progress
#'@importFrom stats quantile
#'@importFrom stats p.adjust
#'
permut_pc_test<-function(pca, pca_data, P=1000, ndim=3, statistic='VAF', conf=0.95,
                         adj.method='BH', perm.method='permV'){

  # P=10
  # ndim=5
  # statistic='commun'
  # conf=0.95
  # adj.method='BH'
  # perm.method='permV'
  # pca<-nlpca
  # pca_data<-nlpca_data

  statistic<-match.arg(statistic, c('VAF','s.loadings', 'commun'))
  perm.method<-match.arg(perm.method, c('permV','permD'))

  original_loadings<-extract_loadings(pca,pca_data)
  original_loadings$Variables<-NULL

  results_list<-list()
  results_list[["method"]]<-"permutation"
  results_list[['statistic']]<-statistic
  results_list[['adj.method']]<-adj.method
  results_list[['pca']]<-pca
  results_list[['pca_data']]<-pca_data
  results_list[['ndim']]<-ndim
  results_list[['conf']]<-ndim

  if (statistic=='VAF' || perm.method=='permD'){
    message("Permuting ", P, " times for ", statistic,
            " using permD method")
    pb <- progress_bar$new(
      format = "[:bar] :percent ~remaining :eta",
      total = P, clear = FALSE, width= 80)

    perm.method<-'permD'
    results_list[['perm.method']]<-'permD'
  }else if ((statistic=='s.loadings'|| statistic=='commun') & perm.method=='permV'){
    message("Permuting ", P," x ",ncol(pca_data), " times for ", statistic,
            " using permV method")
    pb <- progress_bar$new(
      format = "[:bar] :percent ~remaining :eta",
      total = (P)*ncol(pca_data), clear = FALSE, width= 80)

    perm.method<-'permV'
    results_list[['perm.method']]<-'permV'
  }

  per_list<-list()

  if(class(pca)[1]=='prcomp'){
    center=F
    .scale=F
    if (is.numeric(pca$center)) center= TRUE
    if (is.numeric(pca$scale)) .scale= TRUE
  }

  if (perm.method=='permD'){
    per_list<-replicate(P,permut_pca_D(pca, x=pca_data, output = statistic,
                                       center=center, .scale=.scale, pb=pb),
                        simplify = F)
  }else if (perm.method=='permV'){

    per_list<-replicate(P,permut_pca_V(pca, x=pca_data, output = statistic,
                                       center=center, .scale=.scale, ndim = ndim,
                                       original_loadings = original_loadings,pb=pb),
                        simplify = F)
  }

  results_list[['per_samples']]<-per_list


  if (statistic=='VAF'){
    message("\nCalculating VAF...")
    df_per<-do.call(rbind, per_list)

    if (class(pca)[1]=='prcomp'){
      original_VAF<-pca$sdev^2/sum(pca$sdev^2)
      ndim<-min(dim(pca$rotation)[2], ndim)
    }else if (class(pca)[1]%in%c('princals')){
      original_VAF<-pca$evals/sum(pca$evals)
      ndim<-min(pca$ndim, ndim)
    }

    mean_VAF<-colMeans(df_per)

    ci_low<-apply(df_per,2,function(x) quantile(x, (1-conf)/2, na.rm = TRUE))
    ci_high<-apply(df_per,2,function(x) quantile(x, 1-(1-conf)/2, na.rm = TRUE))

    pvalue<-sapply(1:ndim, function(x){
      (sum(df_per[,x]>original_VAF[x])+1)/(length(per_list)+1)#as in Buja & Eyuboglu, 1992;
    })

    results<-as.data.frame(cbind(original_VAF[1:ndim],mean_VAF[1:ndim],
                                 ci_low[1:ndim],ci_high[1:ndim],pvalue))
    colnames(results)<-c('original', 'mean', 'ci_low', 'ci_high','pvalue')

    # if (adj.method=='none'){
    #   results$adj_p_value<-adj.p.value
    # }
    if(adj.method!='none'){
      results<-results%>%
        mutate(adj.p.value=p.adjust(.data$pvalue, method=adj.method))
    }
    rownames(results)<-colnames(original_loadings)[1:ndim]
  }

  if (statistic=='s.loadings'){
    message('\nCalculating loadings...')

    original_loadings<-extract_loadings(pca, pca_data)
    nvars<-dim(pca_data)[2]

    per_df<-as.data.frame(do.call(rbind, per_list))
    colnames(per_df)<-paste0("PC", 1:ncol(per_df))
    per_df$Variables<-rep(original_loadings$Variables,P)
    per_df$id<-rep(1:P, each=nrow(original_loadings))

    per_df<-per_df%>%pivot_longer(-c(.data$id,.data$Variables),
                                  names_to = "component", values_to = "loading")
    original_loadings_long<-original_loadings%>%
      pivot_longer(-.data$Variables, names_to = "component", values_to = "original_loading")

    results<-suppressMessages(per_df%>%left_join(original_loadings_long, by = c("Variables","component"))%>%
      group_by(.data$Variables, .data$component)%>%
      dplyr::summarise(original=mean(.data$original_loading),
                       mean=mean(.data$loading), ci_low=quantile(.data$loading, (1-conf)/2, na.rm = TRUE),
                ci_high=quantile(.data$loading, 1-(1-conf)/2, na.rm = TRUE),
                pvalue=(sum(abs(.data$loading)>abs(.data$original_loading))+1)/(length(per_list)+1))
      )
    if(adj.method!='none'){
      results<-results%>%group_by(.data$component)%>%
        mutate(adj.p.value=p.adjust(.data$pvalue, method=adj.method))
    }
  }

  if (statistic=='commun'){
    message('\nCalculating communalities...')

    original_loadings<-extract_loadings(pca, pca_data)
    original_communalities<-data.frame("original_communality"=rowSums(original_loadings[,1:ndim]^2),
                                       "Variables"=rownames(original_loadings[,1:ndim]))
    nvars<-dim(pca_data)[2]
    per_df<-as.data.frame(do.call(rbind, per_list))
    colnames(per_df)<-paste0("PC", 1:ncol(per_df))
    per_df$Variables<-rep(original_loadings$Variables,P)
    per_df$id<-rep(1:P, each=nrow(original_loadings))

    per_df<-per_df%>%pivot_longer(-c(.data$id,.data$Variables),
                                  names_to = "component", values_to = "loading")

    results<-suppressMessages(per_df%>%
      filter(.data$component%in%colnames(original_loadings))%>%
      left_join(original_communalities,by = "Variables")%>%
      group_by(.data$Variables, .data$id)%>%
      mutate(communality=sum(.data$loading^2))%>%
      group_by(.data$Variables)%>%
      dplyr::summarise(original=mean(.data$original_communality),mean=mean(.data$communality),
                       ci_low=quantile(.data$communality, (1-conf)/2, na.rm = TRUE),
                       ci_high=quantile(.data$communality, 1-(1-conf)/2, na.rm = TRUE),
                       pvalue=(sum(abs(.data$communality)>abs(.data$original_communality))+1)/(length(per_list)+1))
    )

    if(!adj.method%in%'none'){
      results<-results%>%
          mutate(adj.p.value=p.adjust(.data$pvalue, method=adj.method))
      }
  }

  results_list[['results']]<-results
  results_list[['P']]<-P

  message('\n','DONE')
  return(new_syndromics(results_list))
}


#'
# pc_stability2<-function(pca, pca_data, ndim=3, B=1000, sim='ordinary',
#                         test_similarity= TRUE, similarity_metric='all', s_cut_off=0.1,
#                         ci_type='bca', conf=0.95,barmap_plot= TRUE,inParallel=F,...){
#
#   # pca<-prcomp(mtcars[1:10,], center = TRUE, scale. = TRUE)
#   # pca_data<-mtcars[1:10,]
#   # ndim=3
#   # B=200
#   # sim='ordinary'
#   # ci_type='bca'
#   # conf=0.95
#   # test_similarity= TRUE
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
#       summarise_all(.funs=function(x) quantile(x, (1-conf)/2, na.rm = TRUE))
#     ci_sim_high<-similarity_res%>%group_by(.data$PC)%>%
#       summarise_all(.funs=function(x) quantile(x, 1-(1-conf)/2, na.rm = TRUE))
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


