# function to put new ggplot layer below the existing layers
`-.gg` <- function(plot, layer) {
  if (missing(layer)) {
    stop("Cannot use `-.gg()` with a single argument. Did you accidentally put - on a new line?")
  }
  if (!is.ggplot(plot)) {
    stop('Need a plot on the left side')
  }
  plot$layers = c(layer, plot$layers)
  plot
}

# 
# # trim function, for trimming off outliers for plots that give a better picture of the bulk of the data
# # trimming at 1.5x the interquartile range above and below the first and third quartile
# trim <- function(df, trimcol){
#   filter(df, {{ trimcol }} > (quantile({{ trimcol }}, 0.25) - (1.5*IQR({{ trimcol }}))) &
#            {{ trimcol }} < (quantile({{ trimcol }}, 0.75) + (1.5*IQR({{ trimcol }}))))
# }
# 
# # boxplot function for macrofauna data (requiring the existence of a "session" column)
# myboxplot <- function(df, xcol, ycol, ylab, title) {
#   plot <- ggplot(data = df,
#                  aes(x = {{xcol}},
#                      y = {{ycol}})) +
#     geom_boxplot() +
#     theme_minimal() +
#     ggtitle(title) +
#     labs(y = ylab)
#     print(plot)
# }

# bar graph function
# mybargraph <- function(df, xcol, ycol, ylab, title) {
#   plot <- ggplot(data = df,
#                   aes(x = {{xcol}},
#                   y = {{ycol}})) +
#           geom_bar() +
#           theme_minimal() +
#           ggtitle(title) +
#           labs(y = ylab)
#   if(!(all(is.na(df$session)))) {
#     print(plot + facet_wrap(~df$session))
#   } else {
#     print(plot)
#   }
# }


# # 'histogram' function (bar graph)
# myhistogram <- function(df, xcol, xlab) {
#   plot <- ggplot(df, aes(x = {{xcol}})) +
#     geom_histogram() +
#     labs(title = df.name,
#          x = xlab)
#     print(plot)
# }
# 
# # flextable function
# myflextable <- function(df, cap){
#   df %>% 
#     flextable() %>%
#     set_caption(caption = cap) %>%
#     bold(bold = TRUE, part = "header") %>%
#     colformat_double(digits = 1)
#   
# }

# # to show number of observations on boxplots
# give.n <- function(x){
#   return(c(y = median(x)*1.2, label = length(x))) 
#   # experiment with the multiplier to find the perfect position
# }
# 
# # theme for ggplots
# model_theme <- function() {
#   theme(
#     # add border 1)
#     # panel.border = element_rect(colour = "blue", fill = NA, linetype = 2),
#     # color background 2)
#     # panel.background = element_rect(fill = "aliceblue"),
#     # modify grid 3)
#     # panel.grid.major.x = element_line(colour = "steelblue", linetype = 3, size = 0.5),
#     panel.grid = element_blank(),
#     # axis.text = element_text(colour = "steelblue", face = "italic", family = "Times New Roman"),
#     axis.title = element_text(colour = "black", family = "Calibri"),
#     # axis.ticks = element_line(colour = "steelblue"),
#     legend.position = "none"
#   )
# }

# ggplot the model
# wrangling
library(data.table)
library(stringr)
library(forcats)

# analysis
library(nlme)
library(lmerTest)
library(emmeans)
# library(MASS)

# plot function
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(ggforce)
library(insight)
library(cowplot)
library(lazyWeave) #pretty p-values
# library(rstatix) # even prettier p-values -- rather using base R format.pval

pal_okabe_ito <- c(
  "#E69F00",
  "#56B4E9",
  "#009E73",
  "#F0E442",
  "#0072B2",
  "#D55E00",
  "#CC79A7"
)

pal_okabe_ito_FF <- pal_okabe_ito[c(2,6,3,4)]

# ----twoway-fmt_p_value_2, echo=FALSE, warning=FALSE, message=FALSE-------------
fmt_p_value_2 <- function(p, digits = 0.0001){
  p_char <- ifelse(p < 0.06 | is.na(p),
                   scales::pvalue(p,
                                  accuracy = ifelse(
                                    (10^-(ceiling(-log10(p)) + 1)) < digits,
                                    digits,
                                    (10^-(ceiling(-log10(p)) + 1)))),
                   scales::pvalue(p, accuracy = 0.01))
}
fmt_p_value_rmd <- function(x, digits = 0.0001){
  return(sapply(x, fmt_p_value_2, digits))
}


## -------------------------------------------------------------------------------
fmt_table <- function(df, digits = 0.0001){
  
}


## ----odd-even-------------------------------------------------------------------
odd <- function(x) x%%2 != 0
even <- function(x) x%%2 == 0


## ----not_in---------------------------------------------------------------------
'%not_in%' <- Negate('%in%')


## -------------------------------------------------------------------------------
factor_wrap <- function(gg, sep = " "){
  x_labels <- layer_scales(gg)$x$get_limits()
  label_x <- str_replace(x_labels, sep, "\n")
  gg <- gg +
    scale_x_discrete(labels = label_x)
  return(gg)
}



## -------------------------------------------------------------------------------
factor_wrap_1 <- function(gg, sep = " "){
  x_labels <- layer_scales(gg)$x$get_limits()
  label_x <- str_replace(x_labels, sep, "\n")
  return(label_x)
}



## ----remove-parentheses---------------------------------------------------------
remove_parentheses <- function(x){
  if(substr(x, 1, 1) == "("){
    x <- substr(x, 2, nchar(x))
  }
  if(substr(x, nchar(x), nchar(x)) == ")"){
    x <- substr(x, 1, nchar(x)-1)
  }
  return(x)
}


## ----ggcheck_the_qq, warning = FALSE--------------------------------------------
ggcheck_the_qq = function(m1,
                          line = "robust",
                          n_boot = 200){
  n <- nobs(m1)
  m1_res <- residuals(m1)
  #sigma_m1_res <- sigma(m1)
  
  normal_qq <- ppoints(n) %>%
    qnorm()
  sample_qq <- m1_res[order(m1_res)]
  
  # mean + sd
  parametric_slope <- sd(sample_qq)
  parametric_intercept <- mean(sample_qq)
  
  # quartiles
  m1_quartiles <- quantile(m1_res, c(.25, .75))
  
qnorm_quartiles <- qnorm( c(.25, .75))
m1_diff <- m1_quartiles[2] - m1_quartiles[1]
qnorm_diff <- qnorm_quartiles[2] - qnorm_quartiles[1] # = 1.349
quartile_slope <- m1_diff/qnorm_diff
quartile_intercept <- median(m1_quartiles) # median of quartiles not quantiles

# robust uses MASS:rlm (default arguments?)
qq_observed <- data.table(normal_qq = normal_qq,
                          sample_qq = sample_qq)
m2 <- rlm(sample_qq ~ normal_qq, data = qq_observed)
robust_intercept <- coef(m2)[1]
robust_slope <- coef(m2)[2]

# re-sample ribbon
set.seed(1)
resample_qq_model <- numeric(n_boot*n)
Y <- simulate(m1, n_boot)
fd <- model.frame(m1) %>%
  data.table
inc <- 1:n
for(sim_i in 1:n_boot){
  # parametric bound
  fd[, (1) := Y[,sim_i]]
  m1_class <- class(m1)[1]
  if(m1_class == "lm"){
    ff <- lm(formula(m1), data = fd) 
  }
  if(m1_class == "lmerModLmerTest" | m1_class == "lmerMod"){
    ff <- lmer(formula(m1), data = fd)
  }
  y_res <- residuals(ff)
  resample_qq <- y_res[order(y_res)]
  resample_qq_model[inc] <- resample_qq
  inc <- inc + n
  
  # robust bound
  qq_resampled <- data.table(normal_qq = normal_qq,
                             resample_qq = resample_qq)
  m2_resample <- rlm(resample_qq ~ normal_qq, data = qq_resampled)
  
}

qq_sim <- data.table(normal_qq = normal_qq,
                     resample_qq_model = resample_qq_model)

qq_ci_model <- qq_sim[, .(median = median(resample_qq_model),
                          lower = quantile(resample_qq_model, 0.025),
                          upper = quantile(resample_qq_model, 0.975)),
                      by = normal_qq]
m2_boot <- rlm(median ~ normal_qq, data = qq_ci_model)
robust_intercept_boot <- coef(m2_boot)[1]
robust_slope_boot <- coef(m2_boot)[2]

ggplot(data = qq_observed,
       aes(x = normal_qq, y = sample_qq)) +
  
  # ribbon
  geom_ribbon(data = qq_ci_model,
              aes(ymin = lower,
                  ymax = upper,
                  y = median,
                  fill = "band"),
              fill = "gray",
              alpha = 0.6) +
  # draw points
  geom_point() +
  
  # robust
  geom_abline(aes(intercept = robust_intercept,
                  slope = robust_slope,
                  color = "robust"),
              show.legend = FALSE,
              size = 0.75) +
  # robust_boot
  # geom_abline(aes(intercept = robust_intercept_boot,
  #                 slope = robust_slope_boot,
  #                 color = "robust boot"),
  #             show.legend = TRUE,
  #             size = 0.75) +
  xlab("Normal Quantiles") +
  ylab("Sample Quantiles") +
  
  scale_color_manual(values = pal_okabe_ito[c(1:2,5:6)]) +
  theme_minimal_grid() +
  NULL

}


## ----ggcheck_the_spreadlevel----------------------------------------------------
ggcheck_the_spreadlevel <- function(m1,
                                    n_boot = 200){
  n <- nobs(m1)
  m1_res <- residuals(m1)
  m1_scaled <- m1_res/sd(m1_res)
  m1_root <- sqrt(abs(m1_scaled))
  m1_fitted <- fitted(m1)
  
  m2 <- lm(m1_root ~ m1_fitted)
  m2_intercept <- coef(m2)[1]
  m2_slope <- coef(m2)[2]
  
  plot_data <- data.table(
    m1_res = sqrt(abs(m1_scaled)),
    fitted = m1_fitted
  )
  
  ggplot(data = plot_data,
         aes(x = fitted, y = m1_res)) +
    
    # ribbon
    # geom_ribbon(data = qq_ci_model,
    #             aes(ymin = lower,
    #                 ymax = upper,
    #                 y = median,
    #                 fill = "band"),
    #             fill = "gray",
    #             alpha = 0.6) +
    # # draw points
    geom_point() +
    
    geom_smooth(method = lm) +
    # robust 
    # geom_abline(aes(intercept = robust_intercept,
    #                 slope = robust_slope,
    #                 color = "robust"),
    #             show.legend = TRUE,
    #             size = 0.75) +
    # robust_boot
    # geom_abline(aes(intercept = robust_intercept_boot,
    #                 slope = robust_slope_boot,
    #                 color = "robust boot"),
    #             show.legend = TRUE,
  #             size = 0.75) +
  xlab("Fitted") +
    ylab("root abs-scaled-residual") +
    
    
    scale_color_manual(values = pal_okabe_ito[c(1:2,5:6)]) +
    theme_minimal_grid() +
    NULL
}


## ----ggcheck_the_model----------------------------------------------------------
ggcheck_the_model <- function(m1){
  gg1 <- ggcheck_the_qq(m1)
  gg2 <- ggcheck_the_spreadlevel(m1)
  cowplot::plot_grid(gg1, gg2, nrow = 1)
}


## ----emm_table------------------------------------------------------------------
emm_table <- function(fit_emm){
  if(is.data.frame(fit_emm) == TRUE){
    table_out <- data.table(fit_emm)
  }else{
    table_out <- summary(fit_emm) %>%
      data.table()
  }

  if("response" %in% colnames(table_out)){
    # if the lm with log(y) or glm with log link, use " / "
    setnames(table_out,
             old = c("response"),
             new = c("emmean"))
  }
  if("asymp.LCL" %in% colnames(table_out)){
    # if the lm with log(y) or glm with log link, use " / "
    setnames(table_out,
             old = c("asymp.LCL", "asymp.UCL"),
             new = c("lower.CL", "upper.CL"))
  }
  return(table_out)
}


## ----effects_table--------------------------------------------------------------
effects_table <- function(fit_pairs,
                          digits = 3,
                          accuracy = 1e-04,
                          trailing.zero = TRUE,
                          add.p = FALSE){
  if(is.data.frame(fit_pairs) == TRUE){
    pairs_dt <- data.table(fit_pairs)
  }else{
    pairs_dt <- summary(fit_pairs) %>%
      data.table()
  }
  
  if("asymp.LCL" %in% colnames(pairs_dt)){
    # if the lm with log(y) or glm with log link, use " / "
    setnames(pairs_dt,
             old = c("asymp.LCL", "asymp.UCL"),
             new = c("lower.CL", "upper.CL"))
  }
  
  # for simple = each,
  if(which(colnames(pairs_dt)=="contrast") == 3){
    # replace "." with blank
    g_col <- colnames(fit_pairs)[1]
    x_col <- colnames(fit_pairs)[2]
    pairs_dt[get(g_col) == ".", 
             (1) := ""]
    pairs_dt[get(x_col) == ".", 
             (2) := ""]
    
    # beautify contrast column and group1, group2
    pairs_dt[, contrast := paste0(get(names(pairs_dt)[1]),
                                  get(names(pairs_dt)[2]),
                                  ": ",
                                  contrast)]
  }
  
  # create a column of nicely formatted p-values for display.
  # pairs_dt[, p := pvalString(p.value)]
  pairs_dt[, p := p_format(p.value,
                           digits = digits,
                           accuracy = accuracy,
                           trailing.zero = trailing.zero,
                           add.p = add.p)]
  pairs_dt[, contrast := factor(contrast, contrast)]
  return(pairs_dt)
}




## ----pvalue_table---------------------------------------------------------------
pvalue_table <- function(fit_pairs,
                         digits = 1,
                         eps = 0.001,
                         nsmall = 3){
  if(is.data.frame(fit_pairs) == TRUE){
    pairs_dt <- data.table(fit_pairs)
  }else{
    pairs_dt <- summary(fit_pairs) %>%
      data.table()
  }
  
  if("ratio" %in% colnames(pairs_dt)){
    # if the lm with log(y) or glm with log link, use " / "
    groups <- unlist(str_split(pairs_dt$contrast, " / "))
    setnames(pairs_dt, old = "ratio", new = "estimate")
  }else{
    if("odds.ratio" %in% colnames(pairs_dt)){
      # but also if logit(y) as for beta regression (odds.ratio) 
      groups <- unlist(str_split(pairs_dt$contrast, " / "))
      setnames(pairs_dt, old = "odds.ratio", new = "estimate")
    }else{
    # if lm use " - "
    groups <- unlist(str_split(pairs_dt$contrast, " - "))
  }
  }
  
  # remove parentheses if they exist from group cols
  groups <- lapply(groups, remove_parentheses) %>%
    unlist
  
  # add the group1 and group 2 columns
  pairs_dt[, group1 := groups[odd(1:length(groups))]]
  pairs_dt[, group2 := groups[even(1:length(groups))]]
  
  # match for simple and remove parentheses if they exist
  if(names(pairs_dt)[1] != "contrast"){
    x_col <- names(pairs_dt)[2]
    g_col <- names(pairs_dt)[1]
    x_col_levels <- unique(pairs_dt[, get(x_col)])
    x_col_levels <- x_col_levels[x_col_levels != "."]
    g_col_levels <- unique(pairs_dt[, get(g_col)])
    g_col_levels <- g_col_levels[g_col_levels != "."]
    repair <- pairs_dt[, group1] %in% c(x_col_levels, g_col_levels)
    pairs_dt[, group1 := ifelse(repair,
                                group1,
                                gsub(")", "", group1, fixed = TRUE))]
    pairs_dt[, group1 := ifelse(repair,
                                group1,
                                gsub("(", "", group1, fixed = TRUE))]
    repair <- pairs_dt[, group2] %in% c(x_col_levels, g_col_levels)
    pairs_dt[, group2 := ifelse(repair,
                                group2,
                                gsub(")", "", group2, fixed = TRUE))]
    pairs_dt[, group2 := ifelse(repair,
                                group2,
                                gsub("(", "", group2, fixed = TRUE))]
  }
  
  
  if("asymp.LCL" %in% colnames(pairs_dt)){
    # if the lm with log(y) or glm with log link, use " / "
    setnames(pairs_dt,
             old = c("asymp.LCL", "asymp.UCL"),
             new = c("lower.CL", "upper.CL"))
  }
  
  # for simple = each,
  if(which(colnames(pairs_dt)=="contrast") == 3){
    # replace "." with blank
    pairs_dt[get(g_col) == ".", 
             (1) := ""]
    pairs_dt[get(x_col) == ".", 
             (2) := ""]
    
    # beautify contrast column and group1, group2
    pairs_dt[, contrast := paste0(get(names(pairs_dt)[1]),
                                  get(names(pairs_dt)[2]),
                                  ": ",
                                  contrast)]
    pairs_dt[get(g_col) == "",
             group1 := paste0(get(names(pairs_dt)[1]),
                              get(names(pairs_dt)[2]),
                              ",",
                              group1)]
    pairs_dt[get(g_col) == "",
             group2 := paste0(get(names(pairs_dt)[1]),
                              get(names(pairs_dt)[2]),
                              ",",
                              group2)]
    pairs_dt[get(x_col) == "",
             group1 := paste0(group1,
                              ",",
                              get(names(pairs_dt)[1]),
                              get(names(pairs_dt)[2]))]
    pairs_dt[get(x_col) == "",
             group2 := paste0(group2,
                              ",",
                              get(names(pairs_dt)[1]),
                              get(names(pairs_dt)[2]))]
  }
  
  # create a column of nicely formatted p-values for display.
  # pairs_dt[, p := pvalString(p.value)]
  pairs_dt[, p := format.pval(p.value,
                           digits = digits,
                           eps = eps,
                           nsmall = nsmall
                           )]
  pairs_dt[, contrast := factor(contrast, contrast)]
  return(pairs_dt)
}




## ----response-plot--------------------------------------------------------------
response_plot <- function(
    fit, # model fit from lm, lmer, nlme, glmmTMB
    fit_emm, # data frame with means, error
    x_label = "none",
    y_label = "Response (units)",
    g_label = NULL,
    dots = "sina1", # default, in case of many dots
    dodge_width = 0.8,
    adjust = 0.5, # used to be 0.5, when you change it 
    palette = pal_okabe_ito,
    legend_position = "top",
    flip_horizontal = FALSE,
    group_lines = FALSE
){
  dt <- model.frame(fit) %>% # correctly subsets
    data.table()
  y_col <- insight::find_response(fit)
  
  fit_emm_dt <- emm_table(fit_emm)
  
  # get number of pre emmean/ratio/rate columns
  if("emmean" %in% names(fit_emm_dt)){response <- "emmean"}
  if("rate" %in% names(fit_emm_dt)){response <- "rate"}
  if("ratio" %in% names(fit_emm_dt)){response <- "ratio"}
  
  if(which(names(fit_emm_dt) == response) == 2){
    x_col <- names(fit_emm_dt)[1]
    g_col <- NA
  }else{
    x_col <- names(fit_emm_dt)[1]
    g_col <- names(fit_emm_dt)[2]
  }
  
  rand_int <- find_random(fit)
  
  if(is.na(g_col) & group_lines==TRUE){
    group_lines <- FALSE
  }
  if(is.factor(dt[, get(x_col)]) == FALSE){
    dt[, (x_col) := factor(get(x_col))]
  }
  if(!is.na(g_col)){
    if(is.factor(dt[, get(g_col)]) == FALSE){
      dt[, (g_col) := factor(get(g_col))]
    }
  } 
  
  if(is.na(g_col)){
    pd <- position_dodge(width = 0)
  }else{
    pd <- position_dodge(width = dodge_width)
  }
  
  if(is.na(g_col)){
    color_col <- x_col
  }else{
    color_col <- g_col
  }
  
  gg <- ggplot(data=dt, aes(x = get(x_col),
                            y = get(y_col),
                            color = get(color_col))) 
  
  # plot points
  if(dots == "sina1"){
    gg <- gg + geom_sina(alpha = 0.6,
                         position = pd,
                         size = 2.5,
                         jitter_y = FALSE, # to avoid points with value zero to be displayed below zero line
                         adjust = adjust)
  }
  if(dots == "sina2"){ # for bigger dots, when there are less data points
    gg <- gg + geom_sina(alpha = 0.6,
                         position = pd,
                         size = 3.5,
                         adjust = adjust)
  }
  if(dots == "jitter"){
    gg <- gg + geom_point(alpha = 0.4,
                          position = position_jitter(seed = 1))
  }
  if(dots == "dotplot"){
    gg <- gg + geom_dotplot(binaxis='y',
                            stackdir='center',
                            alpha = 0.4,
                            position = pd)
  }
  
  # plot means and CI
  gg <- gg +
    geom_errorbar(data = fit_emm_dt, aes(y = get(response),
                                         ymin = lower.CL,
                                         ymax = upper.CL,
                                         color = get(color_col)),
                  width = 0.1,
                  position = pd,
                  colour = "#333333" #still too dark? Still seems kind of black..
    ) +
    geom_point(data = fit_emm_dt, aes(y = get(response),
                                      color = get(color_col)),
               size = 2.8,
               position = pd,
               colour = "#333333"
    ) +
    
    # aesthetics
    ylab(y_label) +
    scale_color_manual(values = palette,
                       name = color_col) +
    theme_pubr() +
    theme(legend.position = legend_position,
          axis.title = element_text(size = 20),
          axist.text = element_text(size = 18)) + # see next line
    # for some reason, changing size of axis.text here does not change axis label size
    
      NULL
  
  # if(group_lines == TRUE | group_lines == g_col){
  #   gg <- gg +
  #     geom_line(data = fit_emm_dt,
  #               aes(y = emmean,
  #                   group = get(g_col)),
  #               position = pd)
  # }
  # if(group_lines == x_col){
  #   gg <- gg +
  #     geom_line(data = fit_emm_dt,
  #               aes(y = emmean,
  #                   group = get(x_col)))
  # }
  
  if(is.na(g_col)){
    gg <- gg + theme(legend.position="none")
  }
  
  if(is.null(g_label)){
    gg <- gg + guides(color = guide_legend(title=g_col))
  }else{
    if(g_label == "none"){
      gg <- gg + guides(color = guide_legend(title=NULL))
    }else{
      gg <- gg + guides(color = guide_legend(title=g_label))
    }
  }
  
  # remove x axis title
  if(x_label == "none"){
    gg <- gg + theme(axis.title.x = element_blank())
  }else{
    gg <- gg + xlab(x_label)}
  
  if(flip_horizontal == TRUE){
    gg <- gg + coord_flip() +
      theme(legend.position = legend_position)
    gg <- gg + theme(axis.title.x = element_text())
    if(x_label == "none"){
      gg <- gg + theme(axis.title.y = element_blank())
    }else{
      gg <- gg + ylab(x_label)
    }
  }
  
  # gg
  
  return(gg) 
}


## ----plot_pvalues---------------------------------------------------------------
plot_pvalues <- function(
    beta = FALSE,
    gg,
    fit_emm,
    fit_pairs,
    contrast_rows,
    y_pos = NULL
){
  fit_emm_dt <- emm_table(fit_emm)
  pvalue_table_dt <- pvalue_table(fit_pairs)
  
  if(which(names(fit_emm_dt) == "emmean") == 2){
    x_col <- names(fit_emm_dt)[1]
    g_col <- NA
  }else{
    x_col <- names(fit_emm_dt)[1]
    g_col <- names(fit_emm_dt)[2]
  }
  
  
  if(is.na(g_col)){
    #    pd <- position_dodge(width = 0)
    pvalue_table_dt[, x_min_col := group1]
    pvalue_table_dt[, x_max_col := group2]
  }else{
    #    pd <- position_dodge(width = dodge_width)
    pvalue_table_dt[, x_min_col := NA]
    pvalue_table_dt[, x_max_col := NA]
  }
  
  if(is.na(g_col)){
    color_col <- x_col
  }else{
    color_col <- g_col
  }
  
  if(contrast_rows[1] == "all"){
    contrast_rows <- 1:nrow(pvalue_table_dt)
  }
  
  if(sum(contrast_rows) > 0){ # show p-values
    # get x positions for p-values
    # [[3]] is 3rd layer which is means
    if(is.na(pvalue_table_dt[1, x_min_col])){
      gg_data <- cbind(fit_emm_dt,
                       ggplot_build(gg)$data[[3]])
      
      gg_data[, cell := paste(get(x_col), get(g_col), sep=",")]
      match_it_ref <- match(pvalue_table_dt$group1, gg_data$cell)
      
      pvalue_table_dt[, rowa := match(group1, gg_data$cell)]
      pvalue_table_dt[, rowb := match(group2, gg_data$cell)]
      pvalue_table_dt[, x_min_col := gg_data[rowa, x]]
      pvalue_table_dt[, x_max_col := gg_data[rowb, x]]
    }
    
    #    if(is.null(p_pos)){
    p_pos <- 1:length(contrast_rows)
    #   }
    
    # get min/max y
    # conditional: for values that always range between 0 and 1, the max is 1 and min is 0
    # so here, using a dirty workaround for now, we evaluate whether model_contrasts contains 'odds.ratio'
    y_range <- ggplot_build(gg)$layout$panel_params[[1]]$y.range
    max_y <- if(beta == TRUE){1} else y_range[2]
    min_y <- if(beta == TRUE){0} else y_range[1]
    
    increment <- 0.08*(max_y - min_y)
    if(is.null(y_pos)){
      y_position <- max_y + increment*p_pos
    }else{
      y_position <- y_pos
    }
    
    
    gg <- gg + 
      stat_pvalue_manual(pvalue_table_dt[contrast_rows],
                         label = "p",
                         y.position = y_position,
                         xmin = "x_min_col",
                         xmax = "x_max_col",
                         size = 5, 
                         tip.length = 0.01)
    
    # make sure ylim includes p-value
    y_hi <- max(y_position)
    y_lo <- 0 - 0.05*(max_y) 
    # used to be min_y - 0.05*(max_y - min_y), now it's properly scaled, same for all plots
    gg <- gg + coord_cartesian(ylim = c(y_lo, y_hi))
  }
  
  
  
  # gg
  
  return(gg)
}


## ----ggplot_the_response--------------------------------------------------------
ggplot_the_response <- function(
    beta = FALSE, # for layout of p-values for value in (0,1)
    fit, # model fit from lm, lmer, nlme, glmmTMB
    fit_emm,
    fit_pairs,
    wrap_col = NULL,
    x_label = "none",
    y_label = "Response (units)",
    g_label = NULL,
    dots = "sina1", # default
    dodge_width = 0.8, # changing this to 0 didn't make a difference
    adjust = 0.5,
    contrast_rows = "all",
    y_pos = NULL,
    palette = pal_okabe_ito,
    legend_position = "top",
    flip_horizontal = FALSE,
    group_lines = FALSE  
){
  
  # base plot
  gg <- response_plot(
    fit = fit,
    fit_emm = fit_emm,
    x_label = x_label,
    y_label = y_label,
    g_label = g_label,
    dots = dots,
    dodge_width = dodge_width,
    adjust = adjust,
    palette = palette,
    legend_position = legend_position,
    flip_horizontal = flip_horizontal,
    group_lines = group_lines    
  )
  
  # add p-values
  if(contrast_rows[1] != "none" |
     contrast_rows[1] == 0 ){
    gg <- plot_pvalues(
      beta = beta,
      gg = gg,
      fit_emm = fit_emm,
      fit_pairs = fit_pairs,
      contrast_rows = contrast_rows,
      y_pos = y_pos
    )
  }
  
  gg
  return(gg)
}

## ----ggplot_the_effects---------------------------------------------------------
ggplot_the_effects <- function(fit,
                               fit_pairs,
                               contrast_rows = "all",
                               show_p = TRUE,
                               p_position = "top", # "top", "right", "align"
                               nudge_y = 0.33,
                               effect_label = "Effect (units)",
                               effect_x_lim = "auto"){
  
  if(!is.data.frame(fit_pairs)){
    stop("contrast object (fit_pairs) not converted to a data frame. Convert using fit_pairs <- summary(infer = TRUE)")
  }
  
  if(contrast_rows[1] == "all"){
    contrast_rows <- 1:nrow(fit_pairs)
  }
  
  fit_pairs_dt <- effects_table(fit_pairs[contrast_rows,])
  fit_pairs_dt[, contrast := fct_rev(contrast)]
  
  if("lower.CL" %not_in% colnames(fit_pairs_dt)){
    stop("contrast object (fit_pairs) does not include confidence interval columns. Recompute using fit_pairs <- summary(infer = TRUE)")
  }
  
  if("ratio" %in% names(fit_pairs)){
    nill_null <- 1.0
    setnames(fit_pairs_dt, old = "ratio", new = "estimate")
    effect_type <- "ratio"
  }else{
    nill_null <- 0.0
  }
  
  neg_space <- 0.4
  p <- nrow(fit_pairs_dt)
  min_bound <- min(fit_pairs_dt[, lower.CL])
  max_bound <- max(fit_pairs_dt[, upper.CL])
  x_lo <- min_bound - 0.1*abs(min_bound)
  x_hi <- max_bound + 0.1*abs(max_bound)
  if(x_lo > -neg_space*x_hi){
    x_lo <- -neg_space*x_hi
  }
  if(x_hi < -neg_space*x_lo){
    x_hi <- -neg_space*x_lo
  }
  
  if(effect_x_lim[1] == "auto"){
    x_lims <- c(x_lo, x_hi)
  }else
  {
    x_lims <- effect_x_lim
  }
  
  gg <- ggplot(data = fit_pairs_dt,
               aes(y = contrast,
                   x = estimate)) +
    geom_errorbar(aes(xmin=lower.CL, 
                      xmax=upper.CL),
                  width=0, 
                  color="black") +
    geom_point(size = 3) +
    
    geom_vline(xintercept = nill_null, linetype = 2) +
    theme_pubr() +
    xlab(effect_label) +
    theme(axis.title.y = element_blank(),
          axis.text.x = element_blank()) +
    scale_x_continuous(position="top") +
    NULL
  
  gg <- gg + coord_cartesian(xlim = x_lims)
  
  
  nudge_x <- 0
  hjust_i <- 0.5
  if(show_p == TRUE){
    if(p_position == "top"){
      xpos <- fit_pairs_dt[, estimate]
      hjust_i <- 0.5
    }
    if(p_position == "right"){
      xpos <- fit_pairs_dt[, upper.CL] + (x_hi - x_lo)/100
      nudge_y <- 0.0
      hjust_i <- 0
    }
    gg <- gg +
      # p-values
      annotate(geom = "text",
               label = fit_pairs_dt[, p],
               x = xpos + nudge_x,
               y = p:1 + nudge_y,
               size = 2.5,
               hjust = hjust_i)    
  }
  
  # gg
  
  return(gg)
}


## ----ggplot_the_model-----------------------------------------------------------
ggplot_the_model <- function(fit,
                             fit_emm,
                             fit_pairs,
                             x_label = "none",
                             y_label = "Response (units)",
                             g_label = NULL,
                             effect_label = "Effect (units)",
                             dots = "sina",
                             dodge_width = 0.8,
                             adjust = 0.5,
                             contrast_rows = "all",
                             show_p = TRUE,
                             p_position = "top",
                             nudge_y = 0.33,
                             y_pos = NULL,
                             palette = pal_okabe_ito,
                             legend_position = "bottom",
                             flip_horizontal = FALSE,
                             rel_heights = c(1,1),
                             group_lines = FALSE,
                             effect_x_lim = "auto"){
  
  gg1 <- response_plot(
    fit = fit,
    fit_emm = fit_emm,
    x_label = x_label,
    y_label = y_label,
    g_label = g_label,
    dots = dots,
    dodge_width = dodge_width,
    adjust = adjust,
    palette = palette,
    legend_position = legend_position,
    flip_horizontal = flip_horizontal,
    group_lines = group_lines
  )
  
  gg2 <- ggplot_the_effects(fit = fit,
                            fit_pairs = fit_pairs,
                            contrast_rows = contrast_rows,
                            show_p = show_p,
                            p_position = p_position,
                            nudge_y = nudge_y,
                            effect_label = effect_label,
                            effect_x_lim = effect_x_lim)
  plot_grid(gg2,
            gg1,
            nrow=2,
            align = "v",
            axis = "lr",
            rel_heights = rel_heights)
}


## ----plot_treatments------------------------------------------------------------
plot_treatments <- function(gg,
                            x_levels,
                            text_size = 5){
  
  x_pos <- ggplot_build(gg)$data[[2]][, "x"]
  # x_lim <- c(min(ggplot_build(gg)$data[[1]][, "xmin"]),
  #            max(ggplot_build(gg)$data[[1]][, "xmax"]))
  x_range <- ggplot_build(gg)$layout$panel_params[[1]]$x.range
  x_lim <- x_range
  
  y_pos <- 0
  p <- nrow(x_levels)
  y_sep <- 1
  y_breaks <- y_pos + y_sep*(p-1):0
  y_lim <- c(y_pos, y_pos + y_sep*(p-1)) +
    + c(-y_sep/2, y_sep/2)
  
  gg_levels <- ggplot(data = data.table(x_pos, y_pos),
                      aes(x = x_pos,
                          y = y_pos))
  
  for(i in 1:nrow(x_levels)){
    gg_levels <- gg_levels +
      annotate(geom = "text",
               label = x_levels[i,],
               x = x_pos,
               y = y_breaks[i],
               size = text_size)
  }
  gg_levels <- gg_levels +
    coord_cartesian(ylim = y_lim,
                    xlim = x_lim) +
    theme(axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_blank(),
          panel.background = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "pt")) +
    scale_y_continuous(breaks = y_breaks,
                       labels = row.names(x_levels)) +
    scale_x_discrete(expand = expansion(mult = c(0, 0))) +
    NULL 
  
  return(gg_levels)
}




## ----ggplot_the_treatments------------------------------------------------------
ggplot_the_treatments <- function(
    gg,
    x_levels,
    text_size = 5,
    show_groups = FALSE,
    rel_heights = c(1, 0.1)
){
  gg_levels <- plot_treatments(gg,
                               x_levels,
                               text_size)
  
  if(show_groups == FALSE){
    gg <- gg +
      theme(axis.text.x= element_blank())
  }
  
  gg_out <- plot_grid(
    gg,
    gg_levels,
    nrow = 2,
    rel_heights = rel_heights,
    align = "v",
    axis = "lr"
  )
  return(gg_out)
}


## ----vertical-brackets-functions, echo=FALSE------------------------------------
bracketsGrob <- function(...){
  l <- list(...)
  e <- new.env()
  e$l <- l
  grid:::recordGrob(  {
    do.call(grid.brackets, l)
  }, e)
}


geom_bracket <- function(x,
                         y,
                         yend,
                         tip.length = 0.01,
                         line.size = 0.4,
                         label = "",
                         text.size = 4,
                         text.nudge_x = 0.02,
                         text.hjust = 0.5,
                         text.vjust = 0.5,
                         text.angle = 0,
                         parse = FALSE,
                         ...){
  # x_range <- layer_scales(gg)$x$get_limits()
  # if(!is.numeric(x_range)){
  #   x_range <- ggplot_build(gg)$layout$panel_scales_x[[1]]$range_c$range
  # }
  # x_scale <- diff(x_range)
  x_scale <- 1
  
  list(geom_segment(x = x,
                    y = y,
                    xend = x,
                    yend = yend,
                    size = line.size,
                    ...),
       geom_segment(x = x,
                    y = y,
                    xend = x - tip.length*x_scale,
                    yend = y,
                    size = line.size,
                    ...),
       geom_segment(x = x,
                    y = yend,
                    xend = x - tip.length*x_scale,
                    yend = yend,
                    size = line.size,
                    ...),
       annotate(geom = "text",
                x = x + text.nudge_x*x_scale,
                y = y  + text.vjust*(yend - y),
                label = label,
                size = text.size,
                hjust = text.hjust,
                angle = text.angle,
                parse = parse,
                ...)
  )
}
geom_vbracket <- function(gg,
                          x,
                          y,
                          tip.length = 0.02,
                          width = 1,
                          color = "black"){
  x_range <- layer_scales(gg)$x$get_limits()
  if(!is.numeric(x_range)){
    x_range <- ggplot_build(gg)$layout$panel_scales_x[[1]]$range_c$range
  }
  x_scale <- diff(x_range)
  gg <- gg +
    geom_segment(x = x,
                 y = y[1],
                 xend = x,
                 yend = y[2],
                 size = 0.1*width,
                 color = color) +
    geom_segment(x = x,
                 y = y[1],
                 xend = x - tip.length*x_scale,
                 yend = y[1],
                 size = 0.1*width,
                 color = color) +
    geom_segment(x = x,
                 y = y[2],
                 xend = x - tip.length*x_scale,
                 yend = y[2],
                 size = 0.1*width,
                 color = color)
  return(gg)
}

geom_pbracket <- function(gg,
                          x,
                          y,
                          line.tip = 0.02,
                          line.width = 1,
                          label = NULL,
                          text.size = 4,
                          text.hjust = 1,
                          text.vjust = 0.5,
                          parse = FALSE,
                          color = "black"){
  
  x_range <- layer_scales(gg)$x$get_limits()
  if(!is.numeric(x_range)){
    x_range <- ggplot_build(gg)$layout$panel_scales_x[[1]]$range_c$range
  }
  nudge_x <- diff(x_range)*0.01
  
  gg <- geom_vbracket(gg,x,y,line.tip,line.width,color) +
    
    annotate(geom = "text",
             x = x + nudge_x*text.hjust,
             y = y[1]  + text.vjust*diff(y),
             label = label,
             size = text.size,
             color = color,
             hjust = ifelse(text.hjust==1,0,1),
             parse = parse)
  return(gg)
}


geom_vbracket_1 <- function(gg,
                            x,
                            y,
                            xend,
                            yend,
                            type = 4,
                            tip.length = 0.02,
                            width = 1,
                            color = "black"){
  #  y_range <- ggplot_build(gg)$layout$panel_scales_y[[1]]$range$range
  #  y_range <- layer_scales(gg)$y$range$range
  y_range <- layer_scales(gg)$y$get_limits()
  y_diff <- diff(y_range)
  y_diff_expand <- y_diff*1.05
  y_limits <- y_range + c(- 0.5*(y_diff_expand - y_diff),
                          0.5*(y_diff_expand - y_diff))
  x_range <- layer_scales(gg)$x$get_limits()
  if(!is.numeric(x_range)){
    x_range <- ggplot_build(gg)$layout$panel_scales_x[[1]]$range_c$range
  }
  x1 <- (x - x_range[1])/diff(x_range)
  x2 <- (xend - x_range[1])/diff(x_range)
  y1 <- (y - y_limits[1])/y_diff_expand
  y2 <- (yend - y_limits[1])/y_diff_expand
  
  b1 <- bracketsGrob(x1, y1, x2, y2,
                     type = type,
                     h = tip.length,
                     lwd = width,
                     col = color)
  return(b1)
}



## -------------------------------------------------------------------------------
testy <- function(fit_emm){
  exp1c_m1_emm_dt <- fit_emm %>%
    summary()
  b <- exp1c_m1_emm_dt[2, "emmean"]
  gg <- ggplot(data = exp1c_m1_emm_dt,
               aes(x = hac,
                   y = emmean)) +
    geom_point() +
    geom_hline(yintercept = b) +
    geom_hline(yintercept = exp1c_m1_emm_dt[1, "emmean"]) +
    geom_segment(x = 1,
                 y = exp1c_m1_emm_dt[1, "emmean"],
                 xend = 2,
                 yend = exp1c_m1_emm_dt[4, "emmean"]) +
    annotate(geom = "text",
             x = 1.5,
             y = 30,
             label = "test") +
    theme_pubr()
  return(gg)
}




# png("ITS_NMDS_hull.png", width=800, height=600)
# plot_ordination(ps_css, ps_ord_nmds_bray_css, color="landuse") +
#   geom_point(size=1.8) +
#   # stat_ellipse(aes(fill=landuse), geom = "polygon", type = "norm", alpha=0.2, linetype=0) +
#   geom_mark_hull(aes(fill = landuse)) +
#   scale_color_manual(values=vb_colors) + 
#   scale_fill_manual(values=vb_colors) + 
#   theme_minimal(base_size=16) +
#   theme(legend.position = "none")
# dev.off()
pal_okabe_ito_NMDS <- pal_okabe_ito[c(6,2,3,4)]

myNMDS <- function(df, NMDS_res, info = codes_info){
  # works for input (df) of community dataframe and distance/dissimilarity matrix
  
  data_scores <- as.data.frame(NMDS_res[["points"]]) %>%
    rownames_to_column(var = "code") %>%
    left_join(codes_info) %>%
    # to have FF on top for points
    arrange(desc(landuse)) %>%
    # to have order of layering FO < FF < GL < AF for hulls
    mutate(landuse = factor(landuse, levels = c("Forest", "Food_forest",
                                                "Grassland", "Arable_field")))
  
    # species scores only exist if input was community matrix (so not of class "dist")
  if(is.data.frame(df)){
    species_scores <- as.data.frame(scores(NMDS_res, "species"))
    species_scores$species <- rownames(species_scores)
  }
    
  # For the polygons/hulls
  hull_data <- data.frame()
  for(k in c("Food_forest", "Forest", "Grassland", "Arable_field")){
    x <- data_scores[data_scores$landuse == k, ][chull(data_scores[data_scores$landuse ==
                                                                     k, c("MDS1", "MDS2")]), ]
    hull_data <- rbind(hull_data, x)
  }
  
  # For the centroids
  centroids <- data_scores %>% group_by(landuse) %>% summarize(mean1=mean(MDS1), mean2=mean(MDS2))
  
  p <- ggplot(data = data_scores, aes(x = MDS1, y = MDS2,
                                      colour = landuse)) +
    geom_point(size = 1.8) + # before hull to have dots below hulls
    # geom_text_repel(data=data_scores,aes(x=MDS1,y=MDS2,label=code),size=2, hjust=0, vjust=0) +
    geom_mark_hull(aes(fill = landuse), expand = unit(2.5,  "mm"),
                   concavity = 1.5) + # expand: border of hull around dots
    geom_point(data= centroids, aes(x = mean1, y=mean2, bg=landuse), size = 2.2, pch = 21, stroke=0.8, col="black") +
    annotate(geom = "label", 
             x = min(data_scores$MDS1) + 
               ((max(data_scores$MDS1) - min(data_scores$MDS1))/4), 
             y = max(data_scores$MDS2) + 
               ((max(data_scores$MDS2) - min(data_scores$MDS2))/4), 
             size = 5,
             label = paste("Stress: ", 
                           round(NMDS_res$stress, digits = 3))) +
    theme_minimal(base_size = 16) +
    theme(legend.position = "none",
          panel.grid.major = element_line(size = 0.5), 
          panel.grid.minor = element_line(size = 0.5),
          panel.border = element_rect(colour = "black", 
                                      fill = NA, size = 0.5)) +
    # to make the plot a bit lower than the lowest point and higher than the stress label
    scale_y_continuous(limits = c(min(data_scores$MDS2) - 
                                    ((max(data_scores$MDS2) - min(data_scores$MDS2))/12),
                                  max(data_scores$MDS2) + 
                                    ((max(data_scores$MDS2) - min(data_scores$MDS2))/3))) + 
    scale_colour_manual(values = pal_okabe_ito_NMDS) + 
    scale_fill_manual(values = pal_okabe_ito_NMDS)
}

# 
# # Round of numbers in dataframe -------------------------------------------
# 
# round_df <- function(df, digits) {
#   nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
#   
#   df[,nums] <- round(df[,nums], digits = digits)
#   
#   (df)
# }
