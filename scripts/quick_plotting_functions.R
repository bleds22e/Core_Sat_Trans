# Plotting functions for quick looks at data
# Feb 21, 2017

#########################
# LIBRARIES

library(dplyr)
library(ggplot2)

#########################
# Plotting functions

rank_abundance <- function(samples, reads, sp = c('PP', 'DM', 'DO'), cut_off = 0.001) {
  # create rank-abundance plot of mean OTUs w/ error bars for a given species
  
  filtered_samples <- samples[samples$species %in% sp,]
  joined_reads <- inner_join(filtered_samples, reads, by = "Sample")
  mean_sd <- joined_reads %>% group_by(OTU.ID) %>%
    summarise_at(vars(Reads), funs(mean, sd)) %>%
    filter(mean >= cut_off) %>%
    arrange(desc(mean))
  
  ggplot(data = mean_sd) +
    geom_col(aes(x = reorder(OTU.ID, desc(mean)), y = mean)) +
    geom_errorbar(aes(
      x = OTU.ID,
      ymin = mean - sd,
      ymax = mean + sd
    )) +
    labs(x = "OTU.ID", y = "Mean Reads") +
    theme_bw() +
    theme(axis.text.x = element_text(
      angle = 45,
      size = 8,
      hjust = 1
    ))
}

fresh_vs_trap <- function(samples, reads, sp = c('PP', 'DM', 'DO'), cut_off = 0.005) {
  # plot fresh scat vs trap scat for a given species
  
  filtered_samples <- samples[samples$species %in% sp,]
  joined_reads <- inner_join(filtered_samples, reads, by = "Sample")
  mean_sd_by_type <-
    joined_reads %>% group_by(OTU.ID, sample_type) %>%
    summarise_at(vars(Reads), funs(mean, sd)) %>%
    filter(mean >= cut_off) %>%
    arrange(desc(mean))
  
  ggplot(data = mean_sd_by_type, aes(
    x = reorder(OTU.ID, desc(mean)),
    y = mean,
    fill = sample_type
  )) +
    geom_bar(position = "dodge", stat = "identity") +
    geom_errorbar(aes(
      x = OTU.ID,
      ymin = mean - sd,
      ymax = mean + sd
    ), position = position_dodge()) +
    labs(
      x = "OTU.ID",
      y = "Mean Reads"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(
      angle = 45,
      size = 8,
      hjust = 1
    ))
}

control_vs_krat <- function(samples, reads, sp = c('PP', 'DM', 'DO'), sample_type, cut_off = 0.005) {
  # plot reads by species, sample type, and plot type
  filtered <- samples[samples$species %in% sp, ]
  filtered <- filtered[filtered$sample_type %in% sample_type, ]
  
  filtered$plot_type <- filtered$plot
  for (i in 1:length(filtered$plot_type)) {
    if (filtered$plot[i] %in% c(4, 11, 14, 17)) {
      filtered$plot_type[i] = 'control'
    } else {
      filtered$plot_type[i] = 'krat_exclosure'
    }
  }
  
  joined_reads <- inner_join(filtered, reads, by = "Sample")
    
  mean_sd_plot_type <-
    joined_reads %>% group_by(OTU.ID, plot_type) %>%
    summarise_at(vars(Reads), funs(mean, sd)) %>%
    filter(mean >= cut_off) %>%
    arrange(desc(mean))
  
  ggplot(data = mean_sd_plot_type, aes(
    x = reorder(OTU.ID, desc(mean)),
    y = mean,
    fill = plot_type
  )) +
    geom_bar(position = "dodge", stat = "identity") +
    geom_errorbar(aes(
      x = OTU.ID,
      ymin = mean - sd,
      ymax = mean + sd
    ), position = position_dodge()) +
    labs(x = "OTU.ID", y = "Mean Reads") +
    theme_bw() +
    theme(axis.text.x = element_text(
      angle = 45,
      size = 8,
      hjust = 1
    ))
}

sum_by_family <- function(taxa, samples, reads, sp = c('PP', 'DM', 'DO'), cut_off = 0.001) {
  # plot of all reads by (plant) family
  family <- select(taxa, OTU.ID, Family) %>% group_by(Family)
  otu_reads <- reads %>% group_by(OTU.ID)
  filtered_samples <- samples[samples$species %in% sp, ]
  filtered <- inner_join(filtered_samples, otu_reads, by = "Sample")
  family_sum <- inner_join(family, filtered, by = "OTU.ID") %>%
    summarise(sum = sum(Reads)) %>%
    arrange(desc(sum)) %>%
    filter(sum > cut_off)
  
  ggplot(data = family_sum, aes(x = reorder(Family, desc(sum)), y = sum)) +
    geom_col() +
    labs(x = "Plant Family", y = "Total Reads", title = "Reads by Family (> 0.001)") +
    geom_text(aes(label = round((sum/sum(sum)), 3)), vjust = -0.2) +
    theme_bw() +
    theme(axis.text.x = element_text(
      angle = 45,
      size = 8,
      hjust = 1
    ))
}

trap_vs_fresh_indv <- function(samples, reads, sp = c('PP', 'DM', 'DO'), cut_off = 0.01) {
  # plot trap vs fresh samples by individual animal
  filtered_samples <- samples[samples$species %in% sp, ]
  pit_tags <- samples %>% group_by(PIT_tag) %>%
    summarise(num_samples = n_distinct(sample_type)) %>%
    filter(num_samples > 1)
  both_types <- semi_join(samples, pit_tags, by = "PIT_tag")
  add_reads <- inner_join(both_types, reads, by = "Sample") %>%
    tidyr::unite(col = individual, species, PIT_tag, sep = "_") %>%
    filter(Reads > cut_off)
  
  ggplot(data = add_reads, aes(
    x = reorder(OTU.ID, desc(Reads)),
    y = Reads,
    fill = sample_type
  )) +
    geom_bar(position = "dodge", stat = "identity") +
    facet_wrap( ~ individual) +
    labs(x = "OTU.ID", y = "Reads") +
    theme_bw() +
    theme(axis.text.x = element_text(
      angle = 45,
      size = 5,
      hjust = 1
    ))
}