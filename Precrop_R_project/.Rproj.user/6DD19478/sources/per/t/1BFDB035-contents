# Script with functions used in the analysis

# plot model assumptions
plot.mod.vg <- function(mod){
  par(mfrow = c(1,2))
  plot1 <- plot(fitted(mod), resid(mod), xlab = "Fitted values", ylab = "Residuals", main = "Tukey-Anscombe plot")
  plot2 <- car::qqPlot(resid(mod), dist = "norm", mean = mean(resid(mod)), sd = sd(resid(mod)),xlab = "Theoretical quantiles", ylab = "Empirical quantiles", main = "Q-Q plot of residuals")
}

# tidy ANOVA table
create_anova_table <- function(mod){
  t.anova <- Anova(mod) %>% 
    rownames_to_column(var = "Variable") %>%
    tibble()
  
  if(colnames(t.anova)[4] == "Pr(>Chisq)") {
    t.anova <- t.anova %>% 
      mutate(Variable = str_replace_all(Variable, "trt:pre_crop_long", "G x P:"), 
             Variable = str_replace_all(Variable, "trt", "Genotype:"), 
             Variable = str_replace_all(Variable, "pre_crop_long", "Precrop:"),
             across(where(is.numeric), ~ as.character(signif(., 2))),
             across(everything(), ~ replace_na(., "")))
    
  } else {
    t.anova <- t.anova %>% 
      mutate(Variable = str_replace_all(Variable, "trt:pre_crop_long", "G x P:"),
             Variable = str_replace_all(Variable, "trt", "Genotype:"),
             Variable = str_replace_all(Variable, "pre_crop_long", "Precrop:"),
             across(where(is.numeric), ~ as.character(signif(., 2))),
             across(everything(), ~ replace_na(., ""))) 
    
  }

  return(t.anova)
}

# tidy ANOVA table complementation
create_anova_table_comp <- function(mod){
  t.anova <- Anova(mod) %>% 
    rownames_to_column(var = "Variable") %>%
    tibble()
  
  if(colnames(t.anova)[4] == "Pr(>Chisq)") {
    t.anova <- t.anova %>% 
      mutate(Variable = str_replace_all(Variable, "trt_comp", "Treatment:"), 
             across(where(is.numeric), ~ as.character(signif(., 2))),
             across(everything(), ~ replace_na(., "")))
    
  } else {
    t.anova <- t.anova %>% 
      mutate(Variable = str_replace_all(Variable, "trt_comp", "Treatment:"),
             across(where(is.numeric), ~ as.character(signif(., 2))),
             across(everything(), ~ replace_na(., ""))) 
    
  }
  
  return(t.anova)
}

# tidy ANOVA table sterilization
create_anova_table_ster <- function(mod){
  t.anova <- Anova(mod) %>% 
    rownames_to_column(var = "Variable") %>%
    tibble()
  
  if(colnames(t.anova)[4] == "Pr(>Chisq)") {
    t.anova <- t.anova %>% 
      mutate(Variable = str_replace_all(Variable, "trt:pre_x_in", "G x S:"), 
             Variable = str_replace_all(Variable, "trt", "Genotype:"), 
             Variable = str_replace_all(Variable, "pre_x_in", "Soil cond.:"),
             across(where(is.numeric), ~ as.character(signif(., 2))),
             across(everything(), ~ replace_na(., "")))
    
  } else {
    t.anova <- t.anova %>% 
      mutate(Variable = str_replace_all(Variable, "trt:pre_x_in", "G x S:"),
             Variable = str_replace_all(Variable, "trt", "Genotype:"),
             Variable = str_replace_all(Variable, "pre_x_in", "Soil cond.:"),
             across(where(is.numeric), ~ as.character(signif(., 2))),
             across(everything(), ~ replace_na(., ""))) 
    
  }
  
  return(t.anova)
}

# tidy ANOVA table precrop
create_anova_table_pre <- function(mod){
  t.anova <- Anova(mod) %>% 
    rownames_to_column(var = "Variable") %>%
    tibble()
  
  if(colnames(t.anova)[4] == "Pr(>Chisq)") {
    t.anova <- t.anova %>% 
      mutate(Variable = str_replace_all(Variable, "trt:precrop", "G x P:"), 
             Variable = str_replace_all(Variable, "trt", "Genotype:"), 
             Variable = str_replace_all(Variable, "precrop", "Precrop.:"),
             across(where(is.numeric), ~ as.character(signif(., 2))),
             across(everything(), ~ replace_na(., "")))
    
  } else {
    t.anova <- t.anova %>% 
      mutate(Variable = str_replace_all(Variable, "trt:precrop", "G x P:"),
             Variable = str_replace_all(Variable, "trt", "Genotype:"),
             Variable = str_replace_all(Variable, "precrop", "Precrop:"),
             across(where(is.numeric), ~ as.character(signif(., 2))),
             across(everything(), ~ replace_na(., ""))) 
    
  }
  
  return(t.anova)
}

# create ANOVA label (to print in ggplot as text for lm and gls models)
create_anova_label <- function(mod){
  t.anova <- Anova(mod) %>% 
    rownames_to_column(var = "Variable") %>%
    tibble()
  
  if(colnames(t.anova)[4] == "Pr(>Chisq)") {
    t.anova <- t.anova %>% 
      mutate(Variable = str_replace_all(Variable, "trt:pre_crop_long", "G x P:"), 
             Variable = str_replace_all(Variable, "trt", "Genotype:"), 
             Variable = str_replace_all(Variable, "pre_crop_long", "Precrop:"),
             p_new = if_else(`Pr(>Chisq)` < 0.001, paste0("p < 0.001"), 
                             paste0("p = ", sprintf("%.3f", `Pr(>Chisq)`))),
             across(where(is.numeric), ~ as.character(signif(., 2))),
             across(everything(), ~ replace_na(., "")))
    
    t.anova_red <- t.anova
    an_title <- "**ANOVA (GLS)**<Br>"
  } else {
    t.anova <- t.anova %>% 
      mutate(Variable = str_replace_all(Variable, "trt:pre_crop_long", "G x P:"),
             Variable = str_replace_all(Variable, "trt", "Genotype:"),
             Variable = str_replace_all(Variable, "pre_crop_long", "Precrop:"),
             p_new = if_else(`Pr(>F)` < 0.001, paste0("p < 0.001 "), 
                             paste0("p = ", sprintf("%.3f", `Pr(>F)`))),
             across(where(is.numeric), ~ as.character(signif(., 2))),
             across(everything(), ~ replace_na(., ""))) 
    
    t.anova_red <- t.anova %>% 
      filter(Variable != "Residuals")
    an_title <- "**ANOVA**<Br>"
  }
  
  lab_anova <- paste(an_title, 
                     t.anova_red$Variable[1], t.anova_red$p_new[1], "<Br>", 
                     t.anova_red$Variable[2], t.anova_red$p_new[2], "<Br>",
                     t.anova_red$Variable[3], t.anova_red$p_new[3])
  
    return(lab_anova)
}

# create ANOVA label for sub experiment "sterilization" (to print in ggplot as text for lm and gls models)
create_anova_label_ster <- function(mod){
  t.anova <- Anova(mod) %>% 
    rownames_to_column(var = "Variable") %>%
    tibble()
  
  if(colnames(t.anova)[4] == "Pr(>Chisq)") {
    t.anova <- t.anova %>% 
      mutate(Variable = str_replace_all(Variable, "trt:pre_x_in", "G x S:"), 
             Variable = str_replace_all(Variable, "trt", "Genotype:"), 
             Variable = str_replace_all(Variable, "pre_x_in", "Soil cond.:"),
             p_new = if_else(`Pr(>Chisq)` < 0.001, paste0("p < 0.001"), 
                             paste0("p = ", sprintf("%.3f", `Pr(>Chisq)`))),
             across(where(is.numeric), ~ as.character(signif(., 2))),
             across(everything(), ~ replace_na(., "")))
    
    t.anova_red <- t.anova
    an_title <- "**ANOVA (GLS)**<Br>"
  } else {
    t.anova <- t.anova %>% 
      mutate(Variable = str_replace_all(Variable, "trt:pre_x_in", "G x S:"),
             Variable = str_replace_all(Variable, "trt", "Genotype:"),
             Variable = str_replace_all(Variable, "pre_x_in", "Soil cond.:"),
             p_new = if_else(`Pr(>F)` < 0.001, paste0("p < 0.001 "), 
                             paste0("p = ", sprintf("%.3f", `Pr(>F)`))),
             across(where(is.numeric), ~ as.character(signif(., 2))),
             across(everything(), ~ replace_na(., ""))) 
    
    t.anova_red <- t.anova %>% 
      filter(Variable != "Residuals")
    an_title <- "**ANOVA**<Br>"
  }
  
  lab_anova <- paste(an_title, 
                     t.anova_red$Variable[1], t.anova_red$p_new[1], "<Br>", 
                     t.anova_red$Variable[2], t.anova_red$p_new[2], "<Br>",
                     t.anova_red$Variable[3], t.anova_red$p_new[3])
  
  return(lab_anova)
}

# create ANOVA label precrop (to print in ggplot as text for lm and gls models)
create_anova_label_pre <- function(mod){
  t.anova <- Anova(mod) %>% 
    rownames_to_column(var = "Variable") %>%
    tibble()
  
  if(colnames(t.anova)[4] == "Pr(>Chisq)") {
    t.anova <- t.anova %>% 
      mutate(Variable = str_replace_all(Variable, "trt:precrop", "G x P:"), 
             Variable = str_replace_all(Variable, "trt", "Genotype:"), 
             Variable = str_replace_all(Variable, "precrop", "Precrop:"),
             p_new = if_else(`Pr(>Chisq)` < 0.001, paste0("p < 0.001"), 
                             paste0("p = ", sprintf("%.3f", `Pr(>Chisq)`))),
             across(where(is.numeric), ~ as.character(signif(., 2))),
             across(everything(), ~ replace_na(., "")))
    
    t.anova_red <- t.anova
    an_title <- "**ANOVA (GLS)**<Br>"
  } else {
    t.anova <- t.anova %>% 
      mutate(Variable = str_replace_all(Variable, "trt:precrop", "G x P:"),
             Variable = str_replace_all(Variable, "trt", "Genotype:"),
             Variable = str_replace_all(Variable, "precrop", "Precrop:"),
             p_new = if_else(`Pr(>F)` < 0.001, paste0("p < 0.001 "), 
                             paste0("p = ", sprintf("%.3f", `Pr(>F)`))),
             across(where(is.numeric), ~ as.character(signif(., 2))),
             across(everything(), ~ replace_na(., ""))) 
    
    t.anova_red <- t.anova %>% 
      filter(Variable != "Residuals")
    an_title <- "**ANOVA**<Br>"
  }
  
  lab_anova <- paste(an_title, 
                     t.anova_red$Variable[1], t.anova_red$p_new[1], "<Br>", 
                     t.anova_red$Variable[2], t.anova_red$p_new[2], "<Br>",
                     t.anova_red$Variable[3], t.anova_red$p_new[3])
  
  return(lab_anova)
}

# Plot effects of BXs on tolerating precrop legacies
ggplot_precrop <- function(data, x, y, fill, tab.emm, lab.anova, ...) {
  ggplot(data, aes(x = {{x}}, y = {{y}}, fill = {{fill}})) +
    geom_boxplot(position = position_dodge(width = 0.85), alpha = 0.4,
                 color = "grey35", outlier.colour = NA) + 
    geom_quasirandom(dodge.width = 0.85, size = 2, shape = 21,
                     width = 0.085, alpha = 0.7, color = "black",
                     show.legend = FALSE) +
    stat_summary(fun = mean, geom = "point",
                 shape = 18, size = 2.75,
                 position = position_dodge(width = 0.85),
                 show.legend = FALSE) +
    stat_summary(fun.data = mean_se, geom = "errorbar",
                 size = 1, width = 0.4,
                 position = position_dodge(width = 0.85)) +
    scale_fill_manual(name = "Genotype", 
                      labels = fill_labels, 
                      values = fill_values) + 
    scale_x_discrete("Precrop", labels = x_labs,
                     guide = guide_axis(n.dodge = 2)) +
    geom_text(data = tab.emm,
              aes(y = max.y, group = 0.5,
                  label = if_else(p.value < 0.001, paste0("p < 0.001 "),
                                  paste0("p = ", sprintf("%.3f", p.value))))) +
    theme(axis.text.x = element_text(face = "italic"),
          plot.margin = margin(t = 3, r = 8, b = 3, l = 3, unit = "pt"))  +
    geom_richtext(aes(x = Inf, y = Inf, label = lab.anova),
                              stat = "unique",
                              fill = "white", label.color = "white", 
                              label.padding = grid::unit(c(0.5, 1.5, 1.5, 0.75), "lines"),
                              hjust = 0, vjust = 1.0, size = 3.5,
                              label.r = unit(0, "lines")) +
    coord_cartesian(clip = "off")
}

# Plot effects of BXs on tolerating precrop legacies (without dodging axis titles)
ggplot_precrop2 <- function(data, x, y, fill, tab.emm, lab.anova, ...) {
  ggplot(data, aes(x = {{x}}, y = {{y}}, fill = {{fill}})) +
    geom_boxplot(position = position_dodge(width = 0.85), alpha = 0.4,
                 color = "grey35", outlier.colour = NA) + 
    geom_quasirandom(dodge.width = 0.85, size = 2, shape = 21,
                     width = 0.085, alpha = 0.7, color = "black",
                     show.legend = FALSE) +
    stat_summary(fun = mean, geom = "point",
                 shape = 18, size = 2.75,
                 position = position_dodge(width = 0.85),
                 show.legend = FALSE) +
    stat_summary(fun.data = mean_se, geom = "errorbar",
                 size = 1, width = 0.4,
                 position = position_dodge(width = 0.85)) +
    scale_fill_manual(name = "Genotype", 
                      labels = fill_labels, 
                      values = fill_values) + 
    scale_x_discrete(labels = x_labs) +
    geom_text(data = tab.emm,
              aes(y = max.y, group = 0.5,
                  label = if_else(p.value < 0.001, paste0("p < 0.001 "),
                                  paste0("p = ", sprintf("%.3f", p.value)))),
                  size = 3.5) +
    theme(axis.text.x = element_text(face = "italic"),
          plot.margin = margin(t = 3, r = 8, b = 3, l = 3, unit = "pt"))  +
    geom_richtext(aes(x = Inf, y = Inf, label = lab.anova),
                  stat = "unique",
                  fill = "white", label.color = "white", 
                  label.padding = grid::unit(c(0.5, 1.5, 1.5, 0.75), "lines"),
                  hjust = 0, vjust = 1.0, size = 3.5,
                  label.r = unit(0, "lines")) +
    coord_cartesian(clip = "off")
}

# Plot effects of BXs on tolerating precrop legacies for maize_line vs. soil 
ggplot_precrop3 <- function(data, x, y, fill, tab.emm = NULL, lab.anova = NULL, ...) {
  ggplot(data, aes(x = {{x}}, y = {{y}}, fill = {{fill}})) +
    geom_boxplot(position = position_dodge(width = 0.85), alpha = 0.4,
                 color = "grey35", outlier.colour = NA) + 
    geom_quasirandom(dodge.width = 0.85, size = 2, shape = 21,
                     width = 0.085, alpha = 0.7, color = "black",
                     show.legend = FALSE) +
    stat_summary(fun = mean, geom = "point",
                 shape = 18, size = 2.75,
                 position = position_dodge(width = 0.85),
                 show.legend = FALSE) +
    stat_summary(fun.data = mean_se, geom = "errorbar",
                 size = 1, width = 0.4,
                 position = position_dodge(width = 0.85)) +
    scale_fill_manual(name = "Genotype", 
                      labels = fill_labels, 
                      values = fill_values) + 
    scale_x_discrete(labels = x_labs) +
    geom_text(data = tab.emm,
              aes(y = max.y, group = 0.5,
                  label = if_else(p.value < 0.001, paste0("p < 0.001 "),
                                  paste0("p = ", sprintf("%.3f", p.value))))) +
    theme(axis.text.x = element_text(face = "italic"),
          plot.margin = margin(t = 3, r = 8, b = 3, l = 3, unit = "pt"))  +
    coord_cartesian(clip = "off") +
    facet_wrap(vars(soil, maize_line), labeller = facet_labels) +
    theme(strip.text = element_textbox(face = "bold"))
}

# add anova to plot for maize_line and soil
add_anova_tab_lin_soil <- function(tab, x, y, soil, maize_line,  ...) {
  geom_richtext(data = data.frame(soil = {{soil}}, maize_line = {{maize_line}}), 
                mapping =  aes(x = {{x}}, y = {{y}}, label = {{tab}}),
                stat = "unique",
                fill = "white", label.color = "white", 
                label.padding = grid::unit(rep(0.1, 3), "lines"),
                hjust = 0,
                vjust = 1, size = 2.75,
                label.r = unit(0, "lines"))
}

