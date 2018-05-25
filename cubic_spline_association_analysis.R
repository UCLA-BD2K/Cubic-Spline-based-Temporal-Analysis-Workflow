# rm(list = ls(all = T));

library(tidyverse); library(stringr); # leading packages
require(splines); library(factoextra); library(NbClust); library(mclust);

setwd(...);

#####################################################################################################
### load data
#####################################################################################################
hwbw = read_csv("hwbw_ratio.csv");
mb_dta = read_csv("2-21-18 Elution MB modification.csv");
so3h_dta = read_csv("2-21-18 SO3H modification.csv");
so2h_dta = read_csv("2-21-18 SO2H modification.csv");

#####################################################################################################
### filter data
#####################################################################################################
dta_ftr =
  mb_dta %>%
  add_column(cys_type = "reversible") %>%
  bind_rows(so3h_dta %>% 
              add_column(cys_type = "so3h"), 
            so2h_dta %>% 
              add_column(cys_type = "so2h")) %>% 
  select(protein_ptm = cys, iso_over_control = Ratio.H.L.normalizedj, sample = Sample, cys_type) %>% 
  mutate(sample = sample %>% str_replace("7_3_2", "7_1")) %>% 
  separate(sample, into = c("timepoint", "replicate"), sep = "_");

ptm_ftr = 
  dta_ftr %>% 
  spread(key = replicate, value = iso_over_control) %>% 
  mutate(valid_rep_num = apply(!is.na(cbind(`1`,`2`,`3`,`4`)), MARGIN = 1, FUN = sum));

ptm_ftr = 
  ptm_ftr %>% 
  select(protein_ptm, timepoint, valid_rep_num, cys_type) %>% 
  spread(key = timepoint, value = valid_rep_num) %>% 
  mutate(num_more_than_two = ((`10`>1) + (`1`>1) + (`3`>1) + (`5`>1) + (`7`>1) + (`14`>1))) %>% 
  mutate(num_more_than_one = ((`10`>0) + (`1`>0) + (`3`>0) + (`5`>0) + (`7`>0) + (`14`>0))) %>% 
  filter(num_more_than_one >= 4)

dta_ftr =   
  ptm_ftr %>% 
  select(protein_ptm, cys_type) %>% 
  left_join(dta_ftr,
            by = c("protein_ptm", "cys_type"));

#####################################################################################################
### Data De-noising via cubic spline
summary_cubic_spline = function(input_timepoint, input_concentration, input_seq_range){ # cubic spline function 
  inner_tidy = tibble(timepoint = input_timepoint, concentration = input_concentration);
  fit = lm(concentration ~ bs(timepoint, degree = 3), data = inner_tidy);
  
  fit_summary = (summary(fit) %>% .$r.squared);
  fit_p_value = anova(fit)$`Pr(>F)`[1];
  output = paste(predict(fit, data.frame(timepoint = input_seq_range)), collapse = "_");
  
  output_list = list("output" = output, "fit_summary" = fit_summary, "fit_p_value" = fit_p_value);
  
  return(output_list);
}
#####################################################################################################
# sep_num = 50;
seq_range = c(1,3,5,7,10,14);

dta_modeling = 
  dta_ftr %>% 
  mutate(timepoint = parse_double(timepoint)) %>% 
  group_by(protein_ptm, cys_type) %>% 
  dplyr::summarize(models = (summary_cubic_spline(timepoint, iso_over_control, seq_range) %>% .$output),
                   r_squared = (summary_cubic_spline(timepoint, iso_over_control, seq_range) %>% .$fit_summary),
                   p_value = (summary_cubic_spline(timepoint, iso_over_control, seq_range) %>% .$fit_p_value));

dta_modeling = 
  dta_modeling %>% 
  separate(models, into = paste("SP", seq_range, sep = "_"), sep = "_", convert = T) %>% 
  tbl_df();

#####################################################################################################
### clustering
wssplot <- function(data, nc=30, seed=1234){
  wss <- (nrow(data)-1)*sum(apply(data,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss[i] <- sum(kmeans(data, centers=i, nstart = 100)$withinss)}
  plot(1:nc, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")}

sort_collapse = function(a_list){
  output =
    a_list %>% 
    sort(decreasing = T) %>% 
    paste(collapse = ":")
  
  return(output);
}

white2red_palette = colorRampPalette(c("white", "red3"));

multiplot = function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
#####################################################################################################
dta_kmeans = # k-means data
  dta_modeling;

# find the best number of cluster
# wssplot(concentration_kmeans %>% select(SP_0:SP_14));
elbow_result = 
  fviz_nbclust(dta_kmeans %>% 
                 select(SP_1:SP_14),
               kmeans,
               method = "wss",
               k.max = 25,
               nstart = 25,
               nboot = 100) +
  geom_vline(xintercept = 4.5, linetype = 2) + 
  geom_vline(xintercept = 7.5, linetype = 2) + 
  labs(subtitle = "Elbow method");

# perform kmeans 
num_cluster = 6;
fit.km <- kmeans(dta_kmeans %>% 
                   select(SP_1:SP_14), 
                 centers = num_cluster, nstart = 100, iter.max = 100);

### dynamic plots
dynamic_plot = list();
for(each_cluster in c(1:num_cluster)){
  # each_cluster = 1
  
  # data per cluster
  concentration_each_cluster = 
    dta_kmeans %>% 
    filter(fit.km$cluster == each_cluster) %>% 
    gather(key = "sps_timepoint", value = "spline_concentration", SP_1:SP_14) %>% 
    separate(sps_timepoint, into = c("sps", "timepoint"), sep = "_") %>%
    mutate(timepoint = parse_double(timepoint));
  
  concentration_each_cluster_bound =  # library(tolerance); nonparametric tolerance intervals provides 1-sided or 2-sided nonparametric 
    concentration_each_cluster %>%    # (i.e., distribution-free) tolerance intervals for any continuous data set.
    group_by(timepoint) %>% 
    dplyr::summarize(upper_bound = tolerance::nptol.int(spline_concentration, P = 0.95) %>% .$`1-sided.upper`,
                     lower_bound = tolerance::nptol.int(spline_concentration, P = 0.95) %>% .$`1-sided.lower`);
  
  concentration_each_cluster_bound = 
    dplyr::tibble(timepoint = loess.smooth(concentration_each_cluster_bound$timepoint, concentration_each_cluster_bound$upper_bound,family = "gaussian") %>% .$x,
                  upper_bound = loess.smooth(concentration_each_cluster_bound$timepoint, concentration_each_cluster_bound$upper_bound, family = "gaussian") %>% .$y,
                  lower_bound = loess.smooth(concentration_each_cluster_bound$timepoint, concentration_each_cluster_bound$lower_bound, family = "gaussian") %>% .$y,
                  centers = loess.smooth(concentration_each_cluster_bound$timepoint, fit.km$centers[each_cluster,], family = "gaussian") %>% .$y);
  
  # ggplot
  each_dynamic_plot = 
    concentration_each_cluster %>% 
    ggplot(aes(x = timepoint, y = spline_concentration, group = protein_ptm))+ # interaction(strain,metabolite))) + 
    geom_hline(yintercept = 0, color = "black", linetype = 2) +
    geom_ribbon(data = concentration_each_cluster_bound,
                aes(x = timepoint, ymin = lower_bound, ymax = upper_bound), inherit.aes = F,
                alpha = 0.5, fill = "lightskyblue1") +
    geom_line(colour = "red2", size = 0.5, alpha = 0.4) +
    geom_line(colour = "blue3", size = 2.5, aes(x = timepoint, y = centers),
              inherit.aes = F, data = concentration_each_cluster_bound) +
    geom_line(colour = "blue", size = 0.5, aes(x = timepoint, y = upper_bound),
              inherit.aes = F, data = concentration_each_cluster_bound) +
    geom_line(colour = "blue", size = 0.5, aes(x = timepoint, y = lower_bound),
              inherit.aes = F, data = concentration_each_cluster_bound) +
    scale_x_continuous(name = "Timepoint", limits = c(0,14),
                       breaks = c(0, 1, 3, 5, 7, 14)) + 
    # scale_y_continuous(name = "Scaled and Modeled Metabolite Concentration", limits = c(-2,2), 
    #                    breaks = c(-2:2)) +
    scale_y_continuous(name = "", limits = c(0.5, 2), 
                       breaks = c(-1.5, -0.5, 0.5, 1.5)) +
    labs(title = paste("Cluster ", each_cluster, " (n = ", 
                       sum(fit.km$cluster == each_cluster), ")")) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(size = 2, colour = "black"),
          axis.text = element_text(size = 15, colour = "black"),
          axis.title = element_text(size = 15),
          plot.title = element_text(size = 15, hjust = 0.5),
          axis.ticks.length = unit(0.25, "cm"));
  
  dynamic_plot[[each_cluster]] = each_dynamic_plot;
}

pdf(paste("5-04-18 Fig3 Dynamics Plot-all-fix.pdf", sep = ""), width = 12, height = 8);  
multiplot(dynamic_plot[[1]], dynamic_plot[[2]],
          dynamic_plot[[3]], dynamic_plot[[4]], 
          dynamic_plot[[5]], dynamic_plot[[6]],
          cols = 2);
dev.off();

#####################################################################################################
### association study
#####################################################################################################
hwbw_model = 
  hwbw %>% 
  mutate(hwbw_1 = ( `HW1 (g)`/ `BW1 (g)`),
         hwbw_2 = ( `BW2 (g)_1`/ `BW2 (g)`),
         hwbw_3 = ( `BW3 (g)_1`/ `BW3 (g)`)) %>% 
  select(Timepoint, Group, hwbw_1:hwbw_3) %>% 
  gather(key = "hwbw_num", value = "hwbw_val", hwbw_1:hwbw_3) %>% 
  spread(key = "Group", value = "hwbw_val") %>% 
  mutate(ratio_iso_over_control = ISO/Sham) %>% 
  select(Timepoint, ratio_iso_over_control) %>% 
  mutate(ratio_iso_over_control = log2(ratio_iso_over_control))
  
hwbw_model_tmp=
  summary_cubic_spline(hwbw_model$Timepoint, hwbw_model$ratio_iso_over_control, seq_range) %>% 
  .$output %>% 
  str_split("_") %>% 
  .[[1]] %>% 
  parse_double();

### phenotype plot
plot(hwbw_model$Timepoint, hwbw_model$ratio_iso_over_control)#, ylim = c(0.5, 2))
lines(seq_range, hwbw_model_tmp)
phe_data = tibble(timepoint = seq_range, hwbw_model_value = hwbw_model_tmp)

pheno_group_plot =
    hwbw_model %>% 
    ggplot(aes(x = Timepoint, y = ratio_iso_over_control)) +
    geom_point(size = 1.5, alpha = 1) +
  geom_line(aes(x = timepoint, y = hwbw_model_value), inherit.aes = F, data = phe_data, size = 1.5, alpha = 1) +
  scale_x_continuous(name = "Timepoint", limits = c(1,14),
                     breaks = c(1, 3, 5, 7, 10, 14)) + 
  scale_y_continuous(name = "log2(Fold Change of HW/BW in ISO Over CTRL)", limits = c(-0.2, 0.8)) + 
  theme(#panel.grid.major = element_blank(), 
    #panel.grid.minor = element_blank(),
    #panel.background = element_blank(), 
    # legend.position = c(.5,.5),
    legend.title = element_text(size = 10),
    axis.line = element_line(size = 2, colour = "black"),
    axis.text = element_text(size = 15, colour = "black"),
    axis.title = element_text(size = 15),
    plot.title = element_text(size = 15, hjust = 0.5),
    axis.ticks.length = unit(0.25, "cm"));

pdf(paste("Fig4 Phenotype.pdf", sep = ""), width = 7, height = 5, useDingbats = F);  
pheno_group_plot
dev.off()

#####################################################################################################
### extract and display molecules positively and negatively correlated with phenotypic characteristics
#####################################################################################################
dta_corr = 
  dta_modeling %>% 
  mutate(SP_1 = log2(SP_1),
         SP_3 = log2(SP_3),
         SP_5 = log2(SP_5),
         SP_7 = log2(SP_7),
         SP_10 = log2(SP_10),
         SP_14 = log2(SP_14)) %>% 
  gather(key = "timepoint", value = "abundance_ratio", SP_1:SP_14) %>% 
  left_join(tibble(hwbw_man =  hwbw_model_tmp,
                   timepoint = paste("SP", seq_range, sep = "_")),
            by = "timepoint");
  
dta_corr_tmp = 
  dta_corr %>% 
  group_by(protein_ptm, cys_type) %>% 
  summarize(cor_value = cor(abundance_ratio, hwbw_man),
            cor_p_value = 
              cor.test(abundance_ratio, hwbw_man) %>% 
              .$p.value) %>% 
  mutate(cor_adjust_p_value = 
           p.adjust(cor_p_value, "fdr", n = 504));

dta_corr_result = 
  dta_corr_tmp %>% 
  filter(cor_p_value < 0.01)

positives = 
dta_corr_result %>% 
  filter(cor_value > 0)

negatives = 
dta_corr_result %>% 
  filter(cor_value < 0)

positives_val = 
  positives %>% 
  select(protein_ptm, cys_type) %>% 
  left_join(dta_corr,
            by = c("protein_ptm", "cys_type"));

negatives_val = 
  negatives %>% 
  select(protein_ptm, cys_type) %>% 
  left_join(dta_corr,
            by = c("protein_ptm", "cys_type"));

negatives_val
num_cluster = 1;
fit.km <- kmeans(negatives_val %>% select(protein_ptm, timepoint, abundance_ratio) %>% 
                   spread(key = "timepoint", value = "abundance_ratio") %>% 
                   tbl_df() %>% 
                   select(SP_1, SP_3, SP_5, SP_7, SP_10, SP_14),
                 centers = num_cluster)

positives %>% write_csv("positives.csv")
negatives %>% write_csv("negatives.csv")
# tmp =
#   dta_corr %>% 
#   filter(protein_ptm == negative$protein_ptm)
# 
# plot(hwbw_model$Timepoint, hwbw_model$ratio_iso_over_control, ylim = c(0.5, 2))
# lines(seq_range, hwbw_model_tmp)
# lines(seq_range, tmp$abundance_ratio)

each_cluster =
  negatives_val %>% 
  #positives_val %>% 
  select(protein_ptm, timepoint, abundance_ratio) %>% 
  separate(timepoint, into = c("sps", "timepoint"), sep = "_") %>%
  mutate(timepoint = parse_double(timepoint));

each_cluster_bound =  # library(tolerance); nonparametric tolerance intervals provides 1-sided or 2-sided nonparametric 
  each_cluster %>%    # (i.e., distribution-free) tolerance intervals for any continuous data set.
  group_by(timepoint) %>% 
  dplyr::summarize(upper_bound = tolerance::nptol.int(abundance_ratio, P = 0.9) %>% .$`1-sided.upper`,
                   lower_bound = tolerance::nptol.int(abundance_ratio, P = 0.9) %>% .$`1-sided.lower`);

each_cluster_bound = 
  dplyr::tibble(timepoint = loess.smooth(each_cluster_bound$timepoint, each_cluster_bound$upper_bound,family = "gaussian") %>% .$x,
                upper_bound = loess.smooth(each_cluster_bound$timepoint, each_cluster_bound$upper_bound, family = "gaussian") %>% .$y,
                lower_bound = loess.smooth(each_cluster_bound$timepoint, each_cluster_bound$lower_bound, family = "gaussian") %>% .$y,
                centers = loess.smooth(each_cluster_bound$timepoint, fit.km$centers, family = "gaussian") %>% .$y);

### display 
each_dynamic_plot = 
  each_cluster %>% 
  ggplot(aes(x = timepoint, y = abundance_ratio, group = protein_ptm))+ # interaction(strain,metabolite))) + 
  geom_hline(yintercept = 0, color = "black", linetype = 2) +
  geom_ribbon(data = each_cluster_bound,
              aes(x = timepoint, ymin = lower_bound, ymax = upper_bound), inherit.aes = F,
              alpha = 0.5, fill = "lightskyblue1") +
  geom_line(colour = "red2", size = 0.5, alpha = 0.4) +
  geom_line(colour = "blue3", size = 2.5, aes(x = timepoint, y = centers),
            inherit.aes = F, data = each_cluster_bound) +
  geom_line(colour = "blue", size = 0.5, aes(x = timepoint, y = upper_bound),
            inherit.aes = F, data = each_cluster_bound) +
  geom_line(colour = "blue", size = 0.5, aes(x = timepoint, y = lower_bound),
            inherit.aes = F, data = each_cluster_bound) +
  scale_x_continuous(name = "Timepoint", limits = c(0,14),
                     breaks = c(0, 1, 3, 5, 7, 14)) + 
  # scale_y_continuous(name = "Scaled and Modeled Metabolite Concentration", limits = c(-2,2), 
  #                    breaks = c(-2:2)) +
  scale_y_continuous(name = "", limits = c(-0.7,1.2), 
                     breaks = c(-1, -0.5,0, 0.5, 1, 1.5)) +
  # labs(title = paste("Cluster ", each_cluster, " (n = ", 
  #                    sum(fit.km$cluster == each_cluster), ")")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(size = 2, colour = "black"),
        axis.text = element_text(size = 15, colour = "black"),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 15, hjust = 0.5),
        axis.ticks.length = unit(0.25, "cm"));

pdf(paste("positive.pdf", sep = ""), width = 7, height = 5, useDingbats = F);  
pdf(paste("negative.pdf", sep = ""), width = 7, height = 5, useDingbats = F);  
each_dynamic_plot
dev.off()
