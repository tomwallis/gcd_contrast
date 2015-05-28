## script to generate correlation figure (correlations between model params)

# load stuff we need:
rm(list=ls())
library(rstan)
library(ggplot2)
library(psybayes)
library(grid)
library(gridExtra)
library(reshape2)
library(scales)

load(paste0(getwd(),'/output/csf_data_subset.RData'))

model_files <- c(paste0(getwd(),'/output/single_level_glm.RData'),
                 paste0(getwd(),'/output/multilevel_glm.RData'),
                 paste0(getwd(),'/output/single_level_transducer.RData'),
                 paste0(getwd(),'/output/single_level_transducer_B.RData'))
model_names <- c("GLM", "Multi GLM", "Transducer A", "Transducer B")


# function to do scatterplot matrix for a given model ----------

scatterplot_matrix <- function(model_num, params, subject = 1, base_size = 8){
  n_samples <- dim(params[[1]])[1]
  subsample_rows <- sample.int(n_samples, size = 2000)
  
  if (model_num == 1 | model_num == 2){
    n_params <- 3
  } else {
    n_params <- 4
  }

  for (i in 1 : (n_params-1)){
    for (j in (i + 1) : n_params){
      this_plot <- scatterplot_object(model_num = model_num,
                                      params = params, 
                                      to_plot_num = c(j, i),
                                      subject = subject, 
                                      subsamp = subsample_rows)
      
      this_plot <- this_plot + theme_minimal(base_size=base_size) +
        theme(plot.margin = unit(rep(0.04, times = 4), "inches") ) + 
        theme(strip.background = element_rect(fill = "grey80", colour = "White"))
      
      text <- paste0("fig_",i,j," <- this_plot")
      eval(parse(text = text))
    }
  }
  
  # make a blank rect:
  blank <- grid.rect(gp=gpar(col="white"))
  
  # arrange:
  if (model_num == 1 | model_num == 2){
    fig <- arrangeGrob(fig_12, fig_13,
                   blank, fig_23,
                   ncol = 2,
                   main = model_names[model_num])
  } else {
    fig <- arrangeGrob(fig_12, fig_13, fig_14,
                     blank, fig_23, fig_24,
                     blank, blank, fig_34,
                     ncol = 3,
                     main = model_names[model_num])
  }

  return(fig)
  
}


scatterplot_object <- function(model_num, params, to_plot_num, subject = 1, subsamp = NULL){
  if(model_num == 1 | model_num == 2){
    # change param names into index:
    x <- params$beta[, subject, to_plot_num[1]]
    y <- params$beta[, subject, to_plot_num[2]]
    
    par_names <- c("Intercept", "Pedestal", "Increment")
  }
  
  if(model_num == 3 | model_num == 4){
    par_names <- c("p", "q", "z", "rmax")
    text <- paste0("x <- params$",par_names[to_plot_num[1]],"[,", subject, "]")
    eval(parse(text=text))
    text <- paste0("y <- params$",par_names[to_plot_num[2]],"[,", subject, "]")
    eval(parse(text=text))
    
    
  }
  
  df <- data.frame(x = x, y = y)
  # to reduce file size, subset a random points from df:
  if (!is.null(subsamp)){
    df_sub <- df[subsamp, ]
  } else {
    df_sub <- df
  }    
  
  fig <- ggplot(df, aes(x = x, y = y)) + 
#     geom_point(data = df_sub, alpha = 0.2, size = 0.8) + 
    stat_binhex(bins = 20) +
    scale_fill_continuous(guide = FALSE, low = "#636363", high = "#f0f0f0") +
    xlab(par_names[to_plot_num[1]]) + ylab(par_names[to_plot_num[2]])

    return(fig)
}

# function to generate plot correlation matrix of all parameters ----------
correlation_object <- function(model_num, params, subject = 1){
  if(model_num == 1 | model_num == 2){
    par_names <- c("Intcpt", "Ped", "Inc")
    n_samples <- dim(params[[1]])[1]
    # extract params to a matrix:
    df <- as.data.frame(matrix(rep(NA, times = n_samples * length(par_names) ), ncol = length(par_names) ))
    names(df) <- par_names
    for (p in 1 : length(par_names)){
      this_par <- params$beta[, subject, p]
      df[, p] <- this_par
    }
  }
  
  if(model_num == 3 | model_num == 4){
    par_names <- c("p", "q", "z", "rmax")
    n_samples <- dim(params[[1]])[1]
    # extract params to a matrix:
    df <- as.data.frame(matrix(rep(NA, times = n_samples * length(par_names) ), ncol = length(par_names) ))
    names(df) <- par_names
    for (p in 1 : length(par_names)){
      text <- paste0("this_par <- params$", par_names[p],"[, ", subject, "]")
      eval(parse(text = text))
      df[, p] <- this_par
    }    
  }
  
  cor_frame <- melt(cor(df, method = "spearman"))
  # reorder factor levels for plotting:
  cor_frame$Var1 <- factor(cor_frame$Var1, levels = par_names)
  cor_frame$Var2 <- factor(cor_frame$Var2, levels = par_names)
  
  fig <- ggplot(data = cor_frame, aes(x=Var1, y=Var2, fill=value, label=round(value, 2))) + 
    geom_tile() + geom_text(size=2) +
    xlab('') + ylab('') + 
    scale_fill_gradient2(name='Correlation',limits=c(-1,1), guide=FALSE,
                         low = muted("red"),
                         mid = "white", high = muted("blue"))
  
  return(fig)
}



# Do correlation coefficient plot ----------

for (m in 1 : 4){
  load(model_files[m])
  pars <- extract(fit)
  this_fig <- correlation_object(model_num = m, params = pars)
  
  this_fig <- this_fig + 
    theme_minimal(base_size=8) + 
    theme(plot.margin =unit(rep(0.02, times = 4), "inches") ) + 
    theme(strip.background = element_rect(fill = "grey80", colour = "White")) +
    ggtitle(model_names[m])
  
  text <- paste0("fig_", m, " <- this_fig")
  eval(parse(text = text))
}

pdf(file=paste0(getwd(),'/figs/correlation_coeffs.pdf'),width=3.5,height=3.5)

grid.arrange(fig_1, fig_2, fig_3, fig_4,
             ncol = 2, 
             widths = c(0.5, 0.5), heights = c(0.5, 0.5))
dev.off()


# all model scatters in one file? --------

for (m in 1 : 4){
  load(model_files[m])
  pars <- extract(fit)
  this_fig <- scatterplot_matrix(model_num = m, params = pars, base_size = 6)
  
  text <- paste0("fig_", m, " <- this_fig")
  eval(parse(text = text))
}

pdf(file=paste0(getwd(),'/figs/posterior_glm.pdf'),width=3.4,height=5.5)

grid.arrange(fig_1, fig_2,
             ncol = 1, 
             widths = c(0.5, 0.5), heights = c(0.5, 0.5))
dev.off()

pdf(file=paste0(getwd(),'/figs/posterior_transducer.pdf'),width=3.4,height=6)

grid.arrange(fig_3, fig_4,
             ncol = 1, 
             widths = c(0.5, 0.5), heights = c(0.5, 0.5))
dev.off()




# # Do scatterplots for select models----------
# model_num <- 3
# 
# load(model_files[model_num])
# pars <- extract(fit)
# 
# fig <- scatterplot_matrix(model_num = model_num, params = pars)
# 
# pdf(file=paste0(getwd(),'/figs/posterior_',model_names[model_num],".pdf"),width=6.5,height=6.5)
# fig
# dev.off()
# 

