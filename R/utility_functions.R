#'@title syndromic class constructor
#'
new_syndromics<-function(x){
  stopifnot(is.list(x))
  structure(x, class="syndromics")
}


#'@title standardized loadings
#'
#'@description Extract the standardized loadings from a \emph{prcomp} object by correlating the PC scores
#' and the original data.
#'@author Abel Torres Espin
#'
#'@param pca Object of class \emph{prcomp}, \emph{princals}, or \emph{data.frame}. If object is a \emph{prcomp} or \emph{princals} object, \emph{pca_data} is required, and
#'the loadings will be extracted. If object is a data.frame object, the dataframe needs to be formatted as: first column named \emph{Variables} and all other columns corresponding to a PC.
#' One row per variable. The values are the loadings.
#'
#'@param pca_data Data passed to the \emph{prcomp} or \emph{princals} function.
#'
#'@details The standardized loadings are calculated as the eigenvectors times the square roots of the respective eigenvalues and divided by
#' the variable standard deviation (which is 1 in case of standardized PCA (from correlation matrix)). These are equivalent to
#' the Pearson's correlation between the pca scores and the original dataset. This is the correlations of the PC
#' with the variables and the same as the correlation of vector coefficients suggested by Jackson and Hearne in 1987.
#'
#' This function extracts the standardized loadings from the output of the \emph{prcomp()} or the \emph{princals()} functions. In the case of the
#' \emph{prcomp()} solution, the standardized loadings are calculated as: s.loadings = eigenvectors x sqrt(eigenvalues) if the PCA was performed on the
#' standardized (scaled to unit variance) data or s.loadings=(eigenvector x sqrt(eigenvalues))⁄S where S is the vector of the variables standard deviation.
#' In the case of \emph{princals()}, standardized loadings are returned directly in its output and therefore \emph{stand_loadings()} returns those.
#' @return A data.frame with the standardized loadings in the form of variables as rows and components as columns.
#'
#'@export
#'
#'@references Jackson JE, Hearne FT. Relationships Among Coefficients of Vectors Used In Principal Components. Technometrics. 1973 Aug 1;15(3):601–10.
#'
#'@importFrom stats sd
#'
stand_loadings<-function(pca,pca_data){
  if (class(pca)[1]=='prcomp'){
    if (is.numeric(pca$scale)){
      loadings<-as.data.frame((pca$rotation%*%diag(pca$sdev)))
    }else{
      loadings<-as.data.frame((pca$rotation%*%diag(pca$sdev))/apply(pca_data,2,sd))
    }
    colnames(loadings)<-paste('PC', 1:dim(pca$x)[2],sep = '')
  }else if(class(pca)[1]%in%"princals"){
    loadings<-as.data.frame(pca$loadings)
    colnames(loadings)<-paste('PC', 1:pca$ndim,sep = '')
  }else{
    stop('pca must be a prcomp object')
  }
  return(loadings)
}
#'@title extract loadings
#'
#'@description This is a wrapper function for \emph{stand_loadings()}
#'with added functionalities such as error breakers that is used by
#'most functions in the package.
#'
#'@author Abel Torres Espin
#'
#'@param pca Object of class \emph{prcomp}, \emph{princals} or \emph{data.frame}. If the object is a \emph{prcomp} or \emph{princals} object, \emph{pca_data} is required, and
#'the loadings will be extracted. If the object is a data.frame object, the dataframe needs to be formatted as: first column named \emph{Variables} and all other columns corresponding to a PC.
#' One row per variable. The values are the loadings.
#'
#'@param pca_data Data passed to the \emph{prcomp} or \emph{princals} function.
#'
#'@return A data.frame with the standardized loadings in the form of variables as rows and components as columns.
#'
extract_loadings<-function(pca, pca_data){

  if (class(pca)[1]%in%c('prcomp','princals')){
    if (is.null(pca_data)){
      stop('pca_data argument needed when pca is the result of the prcomp function')
    }
    load_df<-stand_loadings(pca, pca_data)
    load_df$Variables<-colnames(pca_data)
  }else if(class(pca)=='data.frame'){
    load_df=pca
    col_names<-colnames(load_df)
    if (!('Variables'%in%col_names)){
      stop('Variables column in the pca data.frame argument do not exist')
    }
    load_df<-load_df[, c('Variables', col_names[!str_detect(col_names,'Variables')])]
    colnames(load_df)[2:length(load_df)]<-paste('PC', 1:(length(load_df)-1), sep ='')
  }else{
    stop('pca must be a prcomp object resulting from prcomp function or a data.frame object with
         the loadings. See documentation for details')
  }

  return(load_df)
}

#'@title Extract the syndromic plot
#'@description Extract the syndromic plot of a given list of loadings and specified pc.
#'
#'@author Abel Torres Espin
#'
#'@param load_df data.frame with the loadings. Format: one column named \emph{Variable} and a column
#'for each PC. It can be extracted from a \emph{prcomp} or \emph{princals} object using \emph{stand_loadings()}.
#'@param pc String. Name of the PC in \emph{load_df} to plot.
#'@param cutoff Numeric. Value of the loadings threshold (i.e. |loadings| > cutoff) to plot. Default = 0.5
#'@param VAF String. Text for the center of the syndromic plot. The text generally corresponds to the variance accounted
#'for (VAF) of the PC. If pca is the results of \emph{prcomp} or \emph{princals}, VAF is internally calculated.
#'@param arrow_size_multi Numeric. Controls the size of the arrows proportional to the loading. Default=10
#'@param repel Logical. Whether to repel the text for preventing text overlap. Default=FALSE
#'@param plot_legend Logical. Whether to plot the legend or not. Default=TRUE
#'@param text_size Numeric. Controls for the size of the text. Default=9
#'@param var_order Character, character vector or numeric vector. Specify the order of the variables in the plot by the loading values, starting at 12 o’clock and moving counterclockwise.
#'Possible values: 'abs decreasing': plot by decreasing absolute value;
#''abs increasing': plot by increasing absolute value; 'decreasing'; or 'increasing’. A vector can be specified with a custom order of the variables
#'to plot. If a numeric vector is provided, the variables will be ordered as specified by the numbers. If a character
#'vector is provided, the variables will be ordered in the same order as the variable names provided. In case of vector, non-specified
#'variables will be excluded from the plot.
#'@param colors Character vector of length 3. Vector with the character name or
#'hexadecimal number (e.g. "#FF0000") of three colors, the lower color, the middle color and the higher color
#'for the gradient used in the plot. Hexadecimal number can be obtained using \emph{rgb} for example.
#'
#'@return A \emph{ggplot2} object of the syndromic plot.
#'
#'
#'@import ggplot2 dplyr ggrepel
#'@importFrom rlang .data
#'
extract_syndromic_plot<-function(load_df, pc,cutoff=0.5, VAF,arrow_size_multi=10,
                                 repel=F,plot_legend=T,text_size=6,
                                 var_order='abs decreasing', colors=c("steelblue1","white","firebrick1")){

  old_scipen<-getOption('scipen')
  on.exit(options(scipen=old_scipen))

  options(scipen=999)

  p<-load_df

  if (length(p$loading)<1){
    stop(paste(pc, ' has no loadings above cutoff'))
  }

  if (is.character(var_order)&length(var_order)==1){
    var_order<-match.arg(var_order, c("abs decreasing","abs increasing","decreasing",
                                      "increasing"))

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
      angletext=.data$div*2*pi/n()+pi/2-pi/22,
      xend=3.5*cos(.data$angle), yend=3.5*sin(.data$angle),
      x=7*cos(.data$angle), y=7*sin(.data$angle),
      xtext=9*cos(.data$angle), ytext=9*sin(.data$angle),
      xload=6*cos(.data$angletext), yload=6*sin(.data$angletext),
      angletext=(.data$angle*180/pi)+ifelse(.data$x<0, 180,0),#
      loading_txt=as.character(round(.data$loading,3)),
      arrow_weight=(round(.data$loading,3)),
      hadjust=ifelse(.data$x<0, 'right','left'),
      hadjust=ifelse(.data$y<(-3) | .data$y>3, 'center',.data$hadjust)
    )

  # p[p$Variables%in%dif_var,'loading_txt']<-
  #   paste("|",
  #         abs(round(p[p$Variables%in%dif_var,'loading'],3)),
  #         "|",sep = "")
  # arbitrary_df<-p%>%
  #   mutate(arrow_weight=ifelse(.data$Variables%in%dif_var,abs(.data$loading),NA))
  #
  # load_df<-p%>%
  #   mutate(loading=ifelse(.data$Variables%in%dif_var,NA,.data$loading))

  angle1<-seq(0,1.2,length.out = 50)
  angle2<-seq(1.95,3.18,length.out = 50)
  angle3<-seq(4.02,5.4,length.out = 50)

  pol<-data.frame(x=c(2.8*cos(angle1)-1,2.8*cos(angle2)+1,2.8*cos(angle3)),
                  y=c(2.8*sin(angle1)-1,2.8*sin(angle2)-1,2.8*sin(angle3)+1))

  s_plot<-ggplot2::ggplot(p,aes(color=.data$loading, label=.data$Variables, x=.data$x, y=.data$y, xend=.data$xend, yend=.data$yend))+
    geom_polygon(data=pol,aes(.data$x,.data$y),inherit.aes = F, fill='grey')+
    ggplot2::scale_color_gradient2(name = "Loading",
                                   high = colors[3], mid = colors[2], low = colors[1],
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
                                yend=-11),color=colors[3],arrow = arrow(type='closed',
                                                                        length = unit(0.1, 'cm'), angle = 25),
                            size=6, show.legend = F,linejoin = 'mitre')+
      ggplot2::geom_segment(aes(x=min(legend_df$x), y=-11, xend=min(legend_df$x)-0.1,
                                yend=-11),color=colors[1],arrow = arrow(type='closed',
                                                                        length = unit(0.1, 'cm'), angle = 25),
                            size=6, show.legend = F,linejoin = 'mitre')+
      ggplot2::annotate(geom='text', x=-4.1, y=-11, label="-1", size=text_size)+
      ggplot2::annotate(geom='text', x=4.1, y=-11, label="1", size=text_size)
  }

  # if (!is.null(dif_var)){
  #   s_plot<-s_plot+
  #     ggnewscale::new_scale_color() +
  #     ggplot2::geom_segment(data=arbitrary_df,aes(color=.data$arrow_weight),
  #                           arrow = arrow(type='closed', length = unit(0.3, 'cm'), angle = 25),
  #                           size=abs(p$arrow_weight)*arrow_size_multi, show.legend = F,
  #                           linejoin = 'mitre')
  #     # scale_color_gradient(low='white', high=dif_var_colors,limits=c(0,1), na.value =  "transparent")
  #     # scale_color_manual(values = dif_var_colors)
  # }

  s_plot<-s_plot+
    ggplot2::annotate(geom='text',x=0, y=0.25, vjust = 'center', hjust = 'center', color='black', label=pc, size=text_size,)+
    ggplot2::annotate(geom='text',x=0, y=-0.75, vjust = 'center', hjust = 'center',color='black', label=VAF,size=text_size)+
    ggplot2::geom_text(aes(label=.data$loading_txt, x=.data$xload, y=.data$yload, angle=.data$angletext),
                       color='black',size=text_size*0.7)

  if (repel){
    s_plot<-s_plot+ggrepel::geom_text_repel(aes(x=.data$xtext, y=.data$ytext,label=.data$Variables, hjust=.data$hadjust),
                                            color='black',min.segment.length = 0.5,size=text_size)
  }else{
    s_plot<-s_plot+ggplot2::geom_text(aes(x=.data$xtext, y=.data$ytext,label=.data$Variables,hjust=.data$hadjust),
                                      color='black',size=text_size)
  }
  return(s_plot)
}


#'@title Extracts different component/factor similarity (matching) indices
#'
#'@description Given a list of loadings for a set of factors or components,
#'computes the Pearson's coefficient of determination (r), the coefficient of congruence (CC),
#'Cattell's S-statistic, and the root mean square error (RMSE) between them.
#'
#'@author Abel Torres Espin
#'
#'@param load.list List of factors to match. Each element of the list is a matrix \emph{p x m} where p are the
#'variables and m the factors or components. All matrices must have variables and components in the same order.
#'@param s_cut_off Numerical value for the loading cut off used to determine if a
#'variable is silent or not in Cattell's terms.
#'@param ndim Numeric. Number of PCs to compute the similarity from. Default=5
#'@param similarity_metric Character or character vector. Possible values are "cc_index" (congruence coefficient),
#'"r_correlation" (Pearson's r), "rmse" (root mean squared error), "s_index' (Cattell's s metric), or "all".
#'Default="all". See below for details on calculations.
#'
#'@details This function is internally called by \emph{pc_stability()}.
#'Each metric is computed using an external function:
#'\describe{
#'  \item{\strong{"cc_index"}\emph{(extract_cc()) function}}{The congurence
#'   coefficient is calculated as: \deqn{CC_{x,y} = sum(x_{i} X y_{i}) / sqrt(sum(x_{i}^2) X sum(y{i})^2)}
#'   Where \eqn{x_{i}} and \eqn{y_{i}} are the loadings of the variable \emph{i} on the
#'   component or factor \emph{x} and \emph{y} respectively. CC is equivalent to the cosine
#'   of the angle between two vectors (the cosine similarity metric) and has a numerical range
#'   from -1 to 1. The sign of a component is arbitrary and can be flipped
#'   without affecting its interpretation. Here we consider the absolute value of CC (0 to 1).
#'   The closer the CC is to 1, the more similar the two components are. (see refs 1,2)}
#'  \item{\strong{"r_correlation"}\emph{(cor()) function}}{The Pearson's r between
#'   two vectors of component loadings has also been used as a similarity metric for
#'   component/factor matching(ref 3). We calculate it here using the \emph{cor()} function.}
#'  \item{\strong{"rmse"}\emph{(extract_rmse()) function}}{RMSE has been also used as a metric
#'   for factor matching (see ref 3). It is calculated as:
#'   \deqn{RMSE_{x,y} = sqrt( sum((x_{i}-y_{i})^2) / n)} Where \emph{n} is the number of variables
#'   in both components \emph{x} and \emph{y}. A RMSE of 0 corresponds to a perfect match.
#'   The smaller the RMSE is, the more equivalent two components are.}
#'  \item{\strong{"s_index"}\emph{(extract_s()) function}}{The s index was first suggested by Cattell et al.
#'   It is based on the factor mandate matrix (ref 4) where loadings are either 1
#'   if a component is considered to act on a variable, called a salient variable, or 0 if not
#'   (forming the hyperplane space). Cattell’s suggested an arbitrary ±0.1 cut-off to be considered as salient variables. In practice,
#'   one might want to alter the threshold depending on the experimental conditions.}
#'
#'  }
#'
#'@return Returns a list of three objects. \strong{Index_all} contains all the comparisons between
#'all the elements of the \emph{load.list}. In general, similarity is calculated between two matrices of
#'loadings, but the user can extract the all the comparisons in case length (load.list) is > 2. \strong{index_mean}
#'is the average of the similarity metrics between all the comparisons. It will be the same as
#'the individual metric (index_all) when length(load.list)==2, because there is only a single comparison made in that scenario.
#'\strong{index_sd} is the standard deviation of the index in case length (load.list) is > 2.
#'
#'@references
#'\enumerate{
#'   \item Burt C. The Factorial Study of Temperamental Traits. Br J Stat Psychol. 1948;1(3):178–203.
#'   \item Tucker, L. R. A method for synthesis of factor analysis studies. Personnel Research Section Report No.984. Washington D.C.: Department of the Army.; 1951.
#'   \item Guadagnoli E, Velicer W. A Comparison of Pattern Matching Indices. Multivar Behav Res. 1991 Apr;26(2):323–43
#'   \item Cattell RB, Balcar KR, Horn JL, Nesselroade JR. Factor Matching Procedures: an Improvement of the s Index; with Tables. Educ Psychol Meas. 1969 Dec;29(4):781–92
#' }
#'
#'@export
#'
#'@importFrom utils combn
#'@importFrom stats cor sd
#'@importFrom rlang .data
#'
component_similarity<-function(load.list, s_cut_off=0.4, ndim=5, similarity_metric='all'){

  similarity_metric<-match.arg(similarity_metric,
                               c("all","cc_index", "r_correlation","rmse",
                                 "s_index"), several.ok = T)

  compar<-t(combn(1:length(load.list),m = 2))
  results<-list()

  if('all'%in%similarity_metric){
    similarity_metric<-c("cc_index", "r_correlation", "rmse","s_index")
  }
  if ('s_index'%in%similarity_metric){
    similarity_metric<-c(similarity_metric,"s_HP")
  }

  results<-lapply(1:dim(compar)[1], function(l){
    lname<-compar[l,]
    temp1_name<-lname[1]
    temp1<-load.list[[temp1_name]][,1:ndim]
    temp1<-as.matrix(apply(temp1, 2, as.numeric))

    temp2_name<-lname[2]
    temp2<-load.list[[temp2_name]][,1:ndim]
    temp2<-as.matrix(apply(temp2, 2, as.numeric))

    if (!identical(dim(temp1), dim(temp2))){
      stop(sprintf("Error: loading matrix %s and %s have different shape",
                   temp1_name,temp2_name))
    }

    rmse<-vector()
    cc_index<-vector()
    s_index<-vector()
    s_HP<-vector()
    for (pc in 1:ndim){
      if('rmse'%in%similarity_metric){
        rmse[pc]<-extract_rmse(temp1[,pc], temp2[,pc])
      }
      if ('cc_index'%in%similarity_metric){
        cc_index[pc]<-extract_cc(temp1[,pc], temp2[,pc])
      }
      if ('s_index'%in%similarity_metric){
        temp_s<-extract_s(temp1[,pc],temp2[,pc],s_cut_off)
        s_index[pc]<-temp_s[[1]]
        s_HP[pc]<-temp_s[[2]]
      }

    }
    if ("r_correlation"%in%similarity_metric){
      r_correlation<-diag(cor(temp1, temp2))
    }

    results_comparison<-list()
    for(s in similarity_metric){
      results_comparison[[s]]<-round(abs(get(s)),4)
    }
    results_comparison<-do.call(cbind, results_comparison)
    comparison<-paste(temp1_name,'vs.',temp2_name, sep = ' ')
    factor<-1:dim(temp1)[2]
    cbind(comparison, factor,results_comparison)
  })

  results.all<-as.data.frame(do.call(rbind,results))
  colnames(results.all)<-c('Comparison','PC',similarity_metric)
  results.all[,3:length(results.all)]<-sapply(results.all[,3:length(results.all)], function(x){as.numeric(as.character(x))})
  results.mean<-results.all[,2:length(results.all)]%>%group_by(.data$PC)%>%
    summarise_all(.funs = function(x){mean(x, na.rm=T)})
  results.sd<-results.all[,2:length(results.all)]%>%group_by(.data$PC)%>%
    summarise_all(.funs = function(x){sd(x, na.rm=T)})

  return(list('index_all'=results.all,'index_mean'=as.data.frame(results.mean),
              'index_sd'=as.data.frame(results.sd)))
}

#'@title Extracts Coefficient of congruence
#'
#'@description Given two vectors, generates the coefficient of congruence between them. This is equivalent to
#'the cosine of the angle between both vectors.
#'
#'@author Abel Torres Espin
#'
#'@param vector1 First numerical vector for the calculation.
#'@param vector2 Second numerical vector for the calculation.
#'
#'@details
#'The congurence
#'   coefficient is calculated as: \deqn{CC_{x,y} = sum(x_{i} X y_{i}) / sqrt(sum(x_{i}^2) X sum(y{i})^2)}
#'   Where \eqn{x_{i}} and \eqn{y_{i}} are the loadings of the variable \emph{i} on the
#'   component or factor \emph{x} and \emph{y} respectively. CC is equivalent to the cosine
#'   of the angle between two vectors (the cosine similarity metric) and has a numerical range
#'   from -1 to 1. The sign of a component is arbitrary and can be flipped
#'   without affecting its interpretation. Here we consider the absolute value of CC (0 to 1).
#'   The closer the CC is to 1, the more similar the two components are.
#'
#'@return Returns the coefficent of congruence (CC) between vector1 and vector2.
#'
#'@references
#'\enumerate{
#'   \item Burt C. The Factorial Study of Temperamental Traits. Br J Stat Psychol. 1948;1(3):178–203.
#'   \item Tucker, L. R. A method for synthesis of factor analysis studies. Personnel Research Section Report No.984. Washington D.C.: Department of the Army.; 1951.
#' }
#'
#'@export
#'
extract_cc<-function(vector1, vector2){
  x<-vector1
  y<-vector2
  x%*%y/sqrt(sum(x^2)*sum(y^2))
}

##'@title Extracts Cattell's S-statistic
#'
#'@description Given two vectors of loadings, computes the Cattell's S-statistic between them with a specified cut off.
#'
#'@author Abel Torres Espin
#'
#'@param vector1 First numerical vector of loadings for the calculation.
#'@param vector2 Second numerical vector of loadings for the calculation.
#'@param cut_off Numerical value for the loading cut off to determine if a
#'variable is silent or not in Cattell's terms. Default = 0.1
#'
#'@details The s index was first suggested by Cattell et al.
#'   It is based on the factor mandate matrix (see ref) where loadings are either 1
#'   if a component is considered to act on a variable, called a salient variable, or 0 if not
#'   (forming the hyperplane space). Cattell’s suggested an arbitrary ±0.1 cut-off to be considered as salient variables. In practice,
#'   one might want to alter the threshold depending on the experimental conditions.
#'
#'@return Returns the Cattell's S-statistic between vector1 and vector2 at cut_off.
#'
#'@references Cattell RB, Balcar KR, Horn JL, Nesselroade JR. Factor Matching Procedures: an Improvement of the s Index; with Tables. Educ Psychol Meas. 1969 Dec;29(4):781–92
#'
#'@export
#'
extract_s<-function(vector1,vector2, cut_off=0.1){

  Positive_salient1<-ifelse(vector1>=cut_off,1,0)
  Negative_salient1<-ifelse(vector1<=-cut_off,1,0)
  Hyperplane1<-ifelse((abs(vector1)<cut_off),1,0)

  Positive_salient2<-ifelse(vector2>=cut_off,1,0)
  Negative_salient2<-ifelse(vector2<=-cut_off,1,0)
  Hyperplane2<-ifelse((abs(vector2)<cut_off),1,0)

  # m.H1<-mean(Hyperplane1)
  # m.H2<-mean(Hyperplane2)

  h.p<-mean(c(Hyperplane1,Hyperplane2))

  c11<-sum(Positive_salient1*Positive_salient2)
  c12<-sum(Positive_salient1*Hyperplane2)
  c13<-sum(Positive_salient1*Negative_salient2)
  c21<-sum(Hyperplane1*Positive_salient2)
  c22<-sum(Hyperplane1*Hyperplane2)
  c23<-sum(Hyperplane1*Negative_salient2)
  c31<-sum(Negative_salient1*Positive_salient2)
  c32<-sum(Negative_salient1*Hyperplane2)
  c33<-sum(Negative_salient1*Negative_salient2)

  s<-(c11+c33-c13-c31)/(c11+c33+c13+c31+0.5*(c12+c21+c23+c32))

  return(list(s, h.p))
}

#'@title Extracts root mean square error (RMSE)
#'
#'@description Given two vectors of loadings, computes the root mean square between them.
#'
#'@author Abel Torres Espin
#'
#'@param vector1 First numerical vector of loadings for the calculation.
#'@param vector2 Second numerical vector of loadings for the calculation.
#'
#'@details RMSE has been also used as a metric
#'   for factor matching (see ref). It is calculated as:
#'   \deqn{RMSE_{x,y} = sqrt( sum((x_{i}-y_{i})^2) / n)} Where \emph{n} is the number of variables
#'   in both components \emph{x} and \emph{y}. A RMSE of 0 corresponds to a perfect match.
#'   The smaller the RMSE is, the more equivalent two components are.
#'
#'@return Returns the root mean square (RMS) between vector1 and vector2.
#'
#'@references Guadagnoli E, Velicer W. A Comparison of Pattern Matching Indices. Multivar Behav Res. 1991 Apr;26(2):323–43
#'
#'@export
#'
extract_rmse<-function(vector1, vector2){
  sqrt(mean((vector1-vector2)^2))
}

#'@title Generic Bootstrapping PCA sample
#'@export
boot_pca_sample<-function(data, indices, pca,...){
  UseMethod("boot_pca_sample", pca)
}

#'@title Bootstrapping PCA sample prcomp method
#'
#'@description This function is passed to the \emph{statistic} argument of the \emph{boot}
#'function from the \strong{boot} package. It generates a bootstrapped PCA sample, extracting the
#'standardized loadings using linear PCA (\emph{prcomp()}).
#'
#'@author Abel Torres Espin
#'
#'@param data The data argument passed from the \emph{boot} function
#'@param indices The indices of the rows for the bootstrapped sample passed from the \emph{boot} function
#'@param pca Object of class \emph{prcomp}.
#'@param original_loadings The data.frame containing the standardized loadings of the original sample.
#'@param ndim Numeric. Number of PCs to save (1 to ndim).
#'@param pb Object of class "Progress_bar" "R6" generated by \emph{progess::progress_bar$new()}. Not required.
#'@param center Logical. Whether pca is conducted on the centered data.
#'@param .scale Logical. Whether pca is conducted on the scaled data.
#'
#'@details A major problem of performing bootstrapping procedures in PCA is what is known as sign reflection:
#' the change of the sign (positive/negative) on the component loadings in a PC given slight variations in the data. In addition,
#'component/factor translocation can occur, meaning that the position of one component can change in a particular
#'PCA solution, especially when two components have similar VAF. Another problem on performing PCAs with variations
#'in the data is the possibility of rotation indeterminacy when the PCA solution of a resampled data presents with a
#'different rotation of the original PCA solution. These issues generate artificially biased bootstrapped distributions,
#'reducing the performance of the procedure. We have implemented a step of procrustes rotation between the original
#'loadings (target) and the bootstrapped sample, which has previously been demonstrated to be a reasonable method to deal with such
#'issues. The procrustes rotation is obtained by the \emph{pracma::procrustes()} function.
#'
#'@return A matrix with one sample bootstrapped standardized loadings. This will be returned to the
#'\emph{boot} function when \emph{pc_stability} function is called.
#'
#'@references
#'\enumerate{
#'  \item Linting M, Meulman JJ, Groenen PJF, van der Kooij AJ. Stability of nonlinear principal components analysis: An empirical study using the balanced bootstrap. Psychol Methods. 2007;12(3):359–79.
#'  \item	Babamoradi H, van den Berg F, Rinnan Å. Bootstrap based confidence limits in principal component analysis — A case study. Chemom Intell Lab Syst. 2013 Jan 15;120:97–105.
#'  \item	Timmerman ME, Kiers HAL, Smilde AK. Estimating confidence intervals for principal component loadings: a comparison between the bootstrap and asymptotic results. Br J Math Stat Psychol. 2007 Nov;60(Pt 2):295–314
#'}
#'
#'
#'@importFrom pracma procrustes
#'@importFrom stats prcomp
#'
#'@export
#'
boot_pca_sample.prcomp<-function(data, indices, pca, original_loadings,ndim,pb=NULL, center, .scale,...){
  d<-data[indices,]

  if(!is.null(pb)){
    pb$tick()
  }

  error<-try(pca_per<-prcomp(d,scale. = .scale,center = center),silent = T)

  if (class(error)!='try-error'){
    load<-stand_loadings(pca_per, d)
    dimension<-min(dim(original_loadings)[2],dim(load)[2])
    original_loadings<-original_loadings[,1:dimension]
    load<-load[,1:dimension]
    r<-pracma::procrustes(as.matrix(original_loadings),as.matrix(load))
    result<-r$P[,1:ndim]
    colnames(result)<-paste('PC',1:ndim, sep = '')

    return(result)
  }else{
    return(NA)
  }
}

#'@title Bootstrapping PCA sample princals method
#'
#'@description This function is passed to the \emph{statistic} argument of the \emph{boot}
#'function from the \strong{boot} package. It generates a bootstrapped PCA sample, extracting the
#'standardized loadings using non-linear PCA (\emph{Gifi::princals()}).
#'
#'@author Abel Torres Espin
#'
#'@param data The data argument passed from the \emph{boot} function
#'@param indices The indices of the rows for the bootstrapped sample passed from the \emph{boot} function
#'@param pca Object of class \emph{princals}
#'@param original_loadings The data.frame containing the standardized loadings of the original sample.
#'@param ndim Numeric. Number of PCs to save (1 to ndim).
#'@param pb Object of class "Progress_bar" "R6" generated by \emph{progess::progress_bar$new()}. Not required.
#'
#'@details A major problem of performing bootstrapping procedures in PCA is what is known as sign reflection:
#' the change of the sign (positive/negative) on the component loadings in a PC given slight variations in the data. In addition,
#'component/factor translocation can occur, meaning that the position of one component can change in a particular
#'PCA solution, especially when two components have similar VAF. Another problem on performing PCAs with variations
#'in the data is the possibility of rotation indeterminacy when the PCA solution of a resampled data presents with a
#'different rotation of the original PCA solution. These issues generate artificially biased bootstrapped distributions,
#'reducing the performance of the procedure. We have implemented a step of procrustes rotation between the original
#'loadings (target) and the bootstrapped sample, which has previously been demonstrated to be a reasonable method to deal with such
#'issues. The procrustes rotation is obtained by the \emph{pracma::procrustes()} function.
#'
#'@return A matrix with one sample bootstrapped standardized loadings. This will be returned to the
#'\emph{boot} function when \emph{pc_stability} function is called.
#'
#'@references
#'\enumerate{
#'  \item Linting M, Meulman JJ, Groenen PJF, van der Kooij AJ. Stability of nonlinear principal components analysis: An empirical study using the balanced bootstrap. Psychol Methods. 2007;12(3):359–79.
#'  \item	Babamoradi H, van den Berg F, Rinnan Å. Bootstrap based confidence limits in principal component analysis — A case study. Chemom Intell Lab Syst. 2013 Jan 15;120:97–105.
#'  \item	Timmerman ME, Kiers HAL, Smilde AK. Estimating confidence intervals for principal component loadings: a comparison between the bootstrap and asymptotic results. Br J Math Stat Psychol. 2007 Nov;60(Pt 2):295–314
#'}
#'
#'
#'@importFrom pracma procrustes
#'@importFrom stats prcomp
#'
#'@export
#'
boot_pca_sample.princals<-function(data, indices, pca, original_loadings, ndim,pb=NULL,...){
  d<-data[indices,]

  if(!is.null(pb)){
    pb$tick()
  }

  pca$call$data<-quote(d)
  error<-try(pca_per<-eval(pca$call), silent = T)

  if (class(error)!='try-error'){
    load<-stand_loadings(pca_per, d)
    dimension<-min(dim(original_loadings)[2],dim(load)[2])
    original_loadings<-original_loadings[,1:dimension]
    load<-load[,1:dimension]
    r<-pracma::procrustes(as.matrix(original_loadings),as.matrix(load))
    result<-r$P[,1:ndim]
    colnames(result)<-paste('PC',1:ndim, sep = '')

    return(result)
  }else{
    return(NA)
  }
}

#'@title Generic for permutating PCA using permD method
permut_pca_D<-function(pca, ...){
  UseMethod("permut_pca_D")
}

#'@title Generic for permuting PCA using permV method
permut_pca_V<-function(pca, ...){
  UseMethod("permut_pca_V")
}

#'@title Permutation PCA sample permD prcomp
#'
#'@description Generates a permuted sample from a PCA that can be used for
#'nonparametric permutation test (see \emph{permut_pc_test}).
#'
#'@author Abel Torres Espin
#'
#'@param pca Object of class \emph{prcomp}
#'@param x Data passed to the \emph{prcomp}.
#'@param center Logical. Whether pca is conducted on the centered data.
#'@param .scale Logical. Whether pca is conducted on the scaled data.
#'@param output Character. Determines the output to compute. Possible values are
#'Variance accounted for ("VAF"), the standardized loadings ("s.loadings") or communalities.
#'@param pb Object of class "Progress_bar" "R6" generated by \emph{progress::progress_bar$new()}. Not required.
#'
#'@details This is a helper function internally called by
#'\emph{permut_pca_test()} to produce a permutation sample of the given output of
#'the \emph{prcomp()}
#'
#'@references
#'\enumerate{
#'  \item Buja A, Eyuboglu N. Remarks on Parallel Analysis. Multivar Behav Res. 1992 Oct 1;27(4):509–40
#'  \item Linting M, van Os BJ, Meulman JJ. Statistical Significance of the Contribution of Variables to the PCA solution: An Alternative Permutation Strategy. Psychometrika. 2011 Jul 1;76(3):440–60
#'}
#'
#'@export
permut_pca_D.prcomp<-function(pca, x, center, .scale,output, pb=NULL,...){

  if(!is.null(pb)){
    pb$tick()
  }
  x<-apply(x, 2,function(x){base::sample(x, length(x))})

  error<-try(pca_per<-prcomp(x,scale. = .scale,center = center),silent = T)

  while (class(error)[1]=='try-error'){
    x<-apply(x, 2,function(x){base::sample(x, length(x))})
    error<-try(pca_per<-prcomp(x,scale. = .scale,center = center),silent = T)
  }

  if (output=='s.loadings'|| output=='commun'){
    return(stand_loadings(pca_per, x))
  }else if (output=='VAF'){
    VAF<-pca_per$sdev^2/sum(pca_per$sdev^2)
    return(VAF)
  }
}

#'@title Permutation PCA sample permD princals
#'
#'@description Generates a permuted sample from a PCA that can be used for
#'nonparametric permutation test (see \emph{permut_pc_test}).
#'
#'@author Abel Torres Espin
#'
#'@param pca Object of class \emph{princals}.
#'@param x Data passed to the \emph{princals}.
#'@param output Character. Determines the output to compute. Possible values are
#'Variance accounted for ("VAF"), the standardized loadings ("s.loadings") or communalities.
#'@param pb Object of class "Progress_bar" "R6" generated by \emph{progress::progress_bar$new()}. Not required.
#'
#'@details This is a helper function internally called by
#'\emph{permut_pca_test()} to produce a permutation sample of the given output of
#'the \emph{prcomp()}
#'
#'@references
#'\enumerate{
#'  \item Buja A, Eyuboglu N. Remarks on Parallel Analysis. Multivar Behav Res. 1992 Oct 1;27(4):509–40
#'  \item Linting M, van Os BJ, Meulman JJ. Statistical Significance of the Contribution of Variables to the PCA solution: An Alternative Permutation Strategy. Psychometrika. 2011 Jul 1;76(3):440–60
#'}
#'
#'@export
#'
permut_pca_D.princals<-function(pca, x, output, pb=NULL,...){

  if(!is.null(pb)){
    pb$tick()
  }
  x<-apply(x, 2,function(x){base::sample(x, length(x))})

  pca$call$data<-quote(x)

  if (is.null(pca$call$knots)){
    pca$call$knots<-quote(knotsGifi(x, "D"))
  }else{
    type<-pca$call$knots$type
    n<-ifelse(is.null(pca$call$knots$n),3,pca$call$knots$n)
    pca$call$knots<-quote(knotsGifi(x, type, n))
  }
  error<-try(pca_per<-eval(pca$call),silent = T)

  while (class(error)[1]=='try-error'){
    x<-apply(x, 2,function(x){base::sample(x, length(x))})
    pca$call$data<-quote(x)

    error<-try(pca_per<-eval(pca$call),silent = T)
  }

  if (output=='s.loadings'|| output=='commun'){
    return(stand_loadings(pca_per, x))
  }else if (output=='VAF'){
    VAF<-pca_per$evals/sum(pca_per$evals)
    return(VAF)
  }
}

#'@title Permutation PCA sample permV prcomp
#'
#'@description Generates a permuted sample from a PCA that can be used for
#'nonparametric permutation test (see \emph{permut_pc_test}).
#'
#'@author Abel Torres Espin
#'
#'@param pca Object of class \emph{prcomp}
#'@param x Data passed to the \emph{prcomp}.
#'@param center Logical. Whether pca is conducted on the centered data.
#'@param .scale Logical. Whether pca is conducted on the scaled data.
#'@param output Character. Determines the output to compute. Possible values are
#'Variance accounted for ("VAF"), the standardized loadings ("s.loadings") or communalities.
#'@param pb Object of class "Progress_bar" "R6" generated by \emph{progress::progress_bar$new()}. Not required.
#'
#'@details This is a helper function internally called by
#'\emph{permut_pca_test()} to produce a permutation sample of the given output of
#'the \emph{prcomp()}
#'
#'@references
#'\enumerate{
#'  \item Buja A, Eyuboglu N. Remarks on Parallel Analysis. Multivar Behav Res. 1992 Oct 1;27(4):509–40
#'  \item Linting M, van Os BJ, Meulman JJ. Statistical Significance of the Contribution of Variables to the PCA solution: An Alternative Permutation Strategy. Psychometrika. 2011 Jul 1;76(3):440–60
#'}
#'
#'@export
#'
permut_pca_V.prcomp<-function(pca, x, output, center,.scale,
                              original_loadings,ndim, pb=NULL,...){

  list_perm_loadings<-list()
  i<-1
  while(i<=ncol(x)){

    perm_x<-x
    perm_x[,i]<-sample(x[,i], length(x[,i]))

    error<-try(pca_per<-prcomp(perm_x,scale. = .scale,center = center),silent = T)

    if (class(error)[1]!='try-error'){

      if(!is.null(pb)){
        pb$tick()
      }

      temp_loading<-stand_loadings(pca_per, perm_x)[,1:ndim]
      if (output=='s.loadings'){
        r<-pracma::procrustes(as.matrix(original_loadings[,1:ndim]),as.matrix(temp_loading))
        list_perm_loadings[[i]]<-r$P[i,]
      }else{
        list_perm_loadings[[i]]<-temp_loading[i,]
      }

      i<-i+1
    }
  }
  list_perm_loadings<-do.call(rbind, list_perm_loadings)
  rownames(list_perm_loadings)<-colnames(x)
  colnames(list_perm_loadings)<-paste0("PC", 1:ndim)

  return(as.data.frame(list_perm_loadings))
}

#'@title Permutation PCA sample permV princals
#'
#'@description Generates a permuted sample from a PCA that can be used for
#'nonparametric permutation test (see \emph{permut_pc_test}).
#'
#'@author Abel Torres Espin
#'
#'@param pca Object of class \emph{princals}.
#'@param x Data passed to the \emph{princals}.
#'@param output Character. Determines the output to compute. Possible values are
#'Variance accounted for ("VAF"), the standardized loadings ("s.loadings") or communalities.
#'@param pb Object of class "Progress_bar" "R6" generated by \emph{progress::progress_bar$new()}. Not required.
#'
#'@details This is a helper function internally called by
#'\emph{permut_pca_test()} to produce a permutation sample of the given output of
#'the \emph{prcomp()}
#'
#'@references
#'\enumerate{
#'  \item Buja A, Eyuboglu N. Remarks on Parallel Analysis. Multivar Behav Res. 1992 Oct 1;27(4):509–40
#'  \item Linting M, van Os BJ, Meulman JJ. Statistical Significance of the Contribution of Variables to the PCA solution: An Alternative Permutation Strategy. Psychometrika. 2011 Jul 1;76(3):440–60
#'}
#'
#'@export
#'
permut_pca_V.princals<-function(pca, x, output,original_loadings,ndim, pb=NULL,...){

  list_perm_loadings<-list()
  i<-1

  while(i<=ncol(x)){

    perm_x<-x
    perm_x[,i]<-sample(x[,i], length(x[,i]))
    pca$call$data<-quote(perm_x)

    if (is.null(pca$call$knots)){
      pca$call$knots<-quote(knotsGifi(perm_x, "D"))
    }else{
      type<-pca$call$knots$type
      n<-ifelse(is.null(pca$call$knots$n),3,pca$call$knots$n)
      pca$call$knots<-quote(knotsGifi(perm_x, type, n))
    }

    error<-try(pca_per<-eval(pca$call),silent = T)


    if (class(error)[1]!='try-error'){

      if(!is.null(pb)){
        pb$tick()
      }

      temp_loading<-stand_loadings(pca_per, perm_x)[,1:ndim]
      if (output=='s.loadings'){
        r<-pracma::procrustes(as.matrix(original_loadings[,1:ndim]),as.matrix(temp_loading))
        list_perm_loadings[[i]]<-r$P[i,]
      }else{
        list_perm_loadings[[i]]<-temp_loading[i,]
      }

      i<-i+1
    }
  }
  list_perm_loadings<-do.call(rbind, list_perm_loadings)
  return(list_perm_loadings)
}

#'@title Method to plot results from pc_stability or permut_pca_test
#'@export
#'
plot.syndromics<-function(x,plot_type="barmap",ndim=NULL,plot_resample=F,
                          communalities=F,...){

  if (is.null(ndim)){
    ndim<-1:x$ndim
  }

  if (!plot_resample){
    x$results<-NULL
    x$communalities<-NULL
  }

  if (x$method=="permutation"){
    if(x$statistic=="VAF"){
      VAF_plot(x$pca, x$pca_data, ndim=ndim,resample_ci=x$results,...)
    }else if(x$statistic=="s.loadings"){
      if (plot_type=="barmap"){
        barmap_loading(x$pca, x$pca_data,resample_ci = x$results,ndim=ndim,...)
      }else if (plot_type=="heatmap"){
        if (plot_resample){
          message("Heatmaps do not plot resample metrics. Try plot_type=\"barmap\"")
        }
        heatmap_loading(x$pca,x$pca_data,ndim=ndim,...)
      }else if (plot_type=="VAF"){
        if (plot_resample){
          message("No permuted VAF available for s.loadings permutation")
        }
        VAF_plot(x$pca, x$pca_data, ndim=ndim,...)
      }else if (plot_type=="syndromics"){
        syndromic_plot(x$pca, x$pca_data,ndim=ndim,...)
      }
    }else if(x$statistic=="commun"){
      if (plot_type!="barmap"){
        message("Barmap are the only available plot type for permuted communalities")
      }
      barmap_commun(x$pca, x$pca_data, load_list = x$per_samples,ndim=ndim,...)
    }
  }else if (x$method=="bootstrap" & !communalities){
    if (plot_type=='barmap'){
      barmap_loading(x$pca, x$pca_data,resample_ci = x$results,ndim=ndim,...)
    }else if(plot_type=="heatmap"){
      if (plot_resample){
        message("Heatmaps do not plot resample metrics. Try plot_type=\"barmap\"")
      }
      heatmap_loading(x$pca,x$pca_data,ndim=ndim,...)
    }else if (plot_type=="VAF"){
      if (plot_resample){
        message("No permuted VAF available for s.loadings permutation")
      }
      VAF_plot(x$pca, x$pca_data, ndim=ndim,...)
    }
  }else if (x$method=="bootstrap" & communalities){
    barmap_commun(x$pca, x$pca_data, load_list = x$boot_samples,ndim=ndim,...)
  }
}

