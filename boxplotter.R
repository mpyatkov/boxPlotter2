## loading R environment
current_dir=getwd()
setwd("/projectnb/wax-es/routines/RENV_443")
source("renv/activate.R")
setwd(current_dir)

library(Rsamtools)
library(GenomicRanges)
library(dplyr)
library(readxl)
library(purrr)
library(stringr)
library(tidyr)
library(cowplot)
library(grid)
library(gridExtra)
library(readr)
library(writexl)
library(ggplot2)
library(patchwork)
library(ggpubr)## for representing stats
library(valr) ## for bam coverage
library(vroom)
library(future)
library(furrr)

CUSTOM_COLORS <- c(
  "#66C2E0",  # Light Blue  
  "#FFB347",  # Light Orange  
  "#77DD77",  # Light Green  
  "#FF6B6B",  # Light Red  
  "#C4A484",  # Light Brown  
  "#C5A5CF",  # Light Purple  
  "#BEBEBE",  # Light Gray  
  "#F9B5D0",  # Light Pink  
  "#7FDBFF",  # Light Cyan  
  "#D4D962",  # Light Olive  
  "#FFEA70",  # Light Gold  
  "#D4E6FF",   # Very Light Blue  
  "#d62728",
  "#17becf" 
)

# ========================================
# CONFIG RELATED FUNCTIONS
# ========================================
parse_config <- function(config_path) {
  groups <- NULL
  comparisons <- NULL
  regions <- NULL
  
  # Get all sheet names
  sheet_names <- excel_sheets(config_path)
  
  # Process each sheet
  for (s in sheet_names) {
    if (str_starts(s, "GROUP")) {
      groups <- read_excel(config_path, col_names = TRUE, sheet = s) %>% 
        select(group_id = 1, group_desc = 2, sample_id = 3)
    } else if (str_starts(s, "COMPARISON")) {
      comparisons <- read_excel(config_path, col_names = TRUE, sheet = s) %>% 
        select(num = 1, group_treat= 2, group_ctrl = 3)
    } else {
      tmp <- read_excel(config_path, col_names = TRUE, sheet = s) %>% 
        select(region_group_id = 1, region_group_name = 2, region_id = 3, chrom = 4, start = 5, end = 6)
      if (is.null(regions)) {
        regions <- tmp
      } else {
        regions <- rbind(regions, tmp)
      }
    }
  }
  list(groups = groups, comparisons = comparisons, regions = regions)
}

add_gr_order_and_colors <- function(config_groups, colors) {
  
  groupid_order <- config_groups %>% 
    select(group_id) %>% 
    mutate(tmp_order = row_number()) %>%
    group_by(group_id) %>% 
    summarise(tmp_order = min(tmp_order)) %>%
    ungroup() %>% 
    arrange(tmp_order) %>% 
    select(-tmp_order) %>% 
    mutate(color_order = row_number())
  
  tmp <- left_join(groupid_order, 
                   tibble(color = colors) %>% mutate(color_order = row_number()), 
                   join_by(color_order)) 
  
  left_join(config_groups, tmp)
}

# Function to load and process configuration
load_and_process_config <- function(config_path, cache_path, progress_reporter = NULL) {
  #config_path <- "./CONFIG_PAIRWISE_BOXPLOT.xlsx"
  #cache_path <- "/projectnb/wax-es/CHIPSEQ_BAM_CACHE/"
  
  CONFIG <- parse_config(config_path)
  CONFIG$groups <- add_gr_order_and_colors(CONFIG$groups, CUSTOM_COLORS)
  
  if (!is.null(progress_reporter)) {
    progress_reporter(value = 1, 
             message = "Loading config step 1 of 1")
  } 
  
  return(list(
    config = CONFIG,
    cache_path = cache_path
  ))
}


# ========================================
# FRIP RELATED FUNCTIONS
# ========================================

# combine all MACS2 narrow peaks and merge them
get_peaks_union <- function(cache_path, samples) {
  lf <- list.files(path = cache_path, 
                   pattern = "narrow_MACS2_peaks.xls", 
                   recursive = T, full.names = T) %>% 
    keep(~str_detect(.,str_c(samples, collapse = "|")))

  peakcaller_xls <- vroom(lf, id = "path", col_names = T, comment = "#", show_col_types = F) %>%
    select(path, seqnames = chr, start, end) %>% 
    mutate(sample_id = str_extract(basename(path), "G[[:alnum:]]+_?M[[:alnum:]]+")) %>%
    select(-path)
  
  peakcaller_xls_union <- makeGRangesFromDataFrame(peakcaller_xls, keep.extra.columns = F) %>% 
    GenomicRanges::reduce(., min.gapwidth = 1L) %>% 
    as_tibble() %>% select(chrom = 1, start = 2, end = 3) %>% 
    mutate(chrom = as.character(chrom))
  
  rm(peakcaller_xls)
  gc()
  
  peakcaller_xls_union
}


calculate_frip <- function(cache_path, samples, peaks_union, progress_reporter = NULL) {
  # for each fragment file need to find overlap with a MACS2 peaks union
  lfrag <- list.files(path = cache_path, pattern = "fragments\\.bed", recursive = T, full.names = T) %>% 
    keep(~str_detect(.,str_c(samples, collapse = "|")))
  
  sid_fragments_in_peaks_union <- function(frag_path, peaks_union) {
    
    fname_prefix <- str_extract(basename(frag_path), "G[[:alnum:]]+_?M[[:alnum:]]+")
    fragm <- read_bed(frag_path, sort = FALSE)
    bedcov <- bed_coverage(peaks_union, fragm)
    n_frag_in_union <- sum(bedcov$.ints)
    total_frag <- nrow(fragm)
    
    rm(list = c("bedcov", "fragm"))
    gc()

    tibble(sample_id = fname_prefix, 
           bam_total_n_fragments = total_frag, 
           n_fragments_in_peaks_union = n_frag_in_union, 
           nfrag_union_ratio = round(n_frag_in_union/total_frag,3))
  }
  
  # imap(lfrag, \(file_path, ix) {
  #   
  #   if (!is.null(progress_reporter)) {
  #     progress_reporter(value = ix/length(lfrag), 
  #                       message = paste("Calculating FRiP", ix, "of", length(lfrag)))
  #   } 
  #   
  #   sid_fragments_in_peaks_union(file_path, peaks_union)}) %>% 
  
  future_map(lfrag, \(file_path) {
    sid_fragments_in_peaks_union(file_path, peaks_union)
  }, .progress = FALSE, .options = furrr_options(seed = TRUE)) %>% 
    
    list_rbind() %>% 
    mutate(rippm = round(min(n_fragments_in_peaks_union)/n_fragments_in_peaks_union, 3),
           frip = round(n_fragments_in_peaks_union/mean(n_fragments_in_peaks_union),3))
}


# Function to get peak union and calculate FRIP (wrapper for functions above)
get_peak_union_and_frip <- function(config_data, progress_reporter = NULL) {
  samples_description <- config_data$config$groups %>% 
    select(sample_id, status = group_desc)
  
  peakcaller_xls_union <- get_peaks_union(
    config_data$cache_path, 
    config_data$config$groups$sample_id
  )
  
  rippm1 <- calculate_frip(
    config_data$cache_path, 
    config_data$config$groups$sample_id, 
    peakcaller_xls_union,
    progress_reporter
  )
  
  rippm_description <- left_join(samples_description, rippm1, join_by(sample_id)) %>% 
    mutate(bam_path = str_glue("{config_data$cache_path}/{sample_id}/bam/{sample_id}_sorted_filtered.bam"))
  
  return(list(
    peak_union = peakcaller_xls_union,
    rippm_description = rippm_description
  ))
}



# ========================================
# REGIONS RELATED FUNCTIONS
# ========================================

calc_reads_in_region <- function(regions, bam_paths, progress_reporter = NULL) {
  
  region <- makeGRangesFromDataFrame(regions, keep.extra.columns = F,
                                     seqnames.field = "seqnames",
                                     start.field = "start",
                                     end.field = "end")
  param <- ScanBamParam(which = region)
  
  ## bam counts
  # imap(bam_paths, \(bam_file, ix) {
  #   
  #   if (!is.null(progress_reporter)) {
  #     progress_reporter(value = ix/length(bam_paths), 
  #              message = paste("Counting reads in regions of interest", ix, "of", length(bam_paths)))
  #   } 
  #   
  #   countBam(bam_file, param = param)
  # })
  future_map(bam_paths, \(bam_file) {
    countBam(bam_file, param = param)
  }, .progress = FALSE, .options = furrr_options(seed = TRUE)) %>% list_rbind()
}

# Function to process regions and calculate read counts
process_regions_and_counts <- function(config_data, peak_data, progress_reporter = NULL) {
  all_regions <- config_data$config$regions %>% 
    mutate(dataset_id = str_glue("{region_group_id}_{region_group_name}")) %>% 
    select(region_group_id, region_id, seqnames = chrom, start, end, dataset_id)
  
  bam_counts <- calc_reads_in_region(all_regions, peak_data$rippm_description$bam_path, progress_reporter)
  
  counts_norm <- bam_counts %>% 
    mutate(sample_id = str_replace(file, "_sorted_filtered.bam", "")) %>% 
    left_join(., peak_data$rippm_description) %>% 
    left_join(., all_regions %>% select(seqnames, start, end, dataset_id, region_id, region_group_id)) %>% 
    left_join(., config_data$config$groups %>% select(sample_id, group_id, color_order), join_by(sample_id)) %>% 
    mutate(
      records_width_norm = 1000 * records / width,
      counts_width_frip_norm = records_width_norm / frip,
      log2data = ifelse(counts_width_frip_norm == 0, 0, log2(counts_width_frip_norm))
    ) %>% 
    mutate(
      log2data_group_avg = sum(log2data)/n(), 
      .by = c("seqnames", "start", "end", "status")
    )
  
  return(list(
    all_regions = all_regions,
    bam_counts = bam_counts,
    counts_norm = counts_norm
  ))
}

# ========================================
# PLOTS RELATED FUNCTIONS
# ========================================

create_pairwise_plots_for_one_dataset <- function(data, limits, title,
                                                  color_palette = CUSTOM_COLORS[1:2]) {
  # Helper function for ordered factors
  ordered_factor <- function(x,d) reorder(x, as.numeric(factor(d$group_order)))
  
  # Common theme and aesthetics
  common_theme <- theme_pubr() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          plot.title = element_text(hjust = 0.5, size = 10),
          panel.grid.major.y = element_line())
  common_y_scale <- scale_y_continuous(limits = limits, breaks = scales::pretty_breaks(10))
  
  
  individual_plot <- ggplot(data, aes(x = ordered_factor(sample_id, data), y = log2data, 
                                      fill = ordered_factor(status_up, data))) +
    geom_boxplot() + scale_fill_manual(name = "Group", values = color_palette) +
    common_theme + ggtitle("Individual") + xlab("") + ylab("log2 FRiP")+
    common_y_scale
  
  all_points_plot <- ggplot(data, aes(x = ordered_factor(status_up, data), y = log2data, 
                                      fill = ordered_factor(status_up, data))) +
    geom_boxplot() + scale_fill_manual(name = "Group", values = color_palette) +
    common_theme + ggtitle("All points") + xlab("") + ylab("")+
    common_y_scale
  
  avg_data <- data %>% select(status_up, group_id, log2data_group_avg, group_order) %>% distinct()
  avg_plot <- ggplot(avg_data, aes(x = ordered_factor(status_up, avg_data), y = log2data_group_avg, 
                                   fill = ordered_factor(status_up, avg_data))) +
    geom_boxplot() + scale_fill_manual(name = "Group", values = color_palette) +
    common_theme + ggtitle("Average by\nsamples in group") + xlab("") + ylab("")+
    common_y_scale
  
  list(individual = individual_plot,
       avg = avg_plot,
       all_points = all_points_plot)
}

generate_individual_sample_plots <- function(counts_norm, preset_colors, output_dir = "./") {
  all_samples_together <- counts_norm %>% 
    group_by(dataset_id) %>% 
    group_map(\(data, info) {
      ds_id <- info$dataset_id
      ggplot(data, aes(x = reorder(sample_id, as.numeric(factor(status))), y = log2data, fill = status)) +
        geom_boxplot() +
        scale_fill_manual(name = "Group name", values = preset_colors) +
        theme_pubr(legend = "right") +
        scale_y_continuous(
          breaks = scales::pretty_breaks(10), 
          limits = c(min(counts_norm$log2data), max(counts_norm$log2data))) +
        theme(
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          panel.grid.major.y = element_line()
        ) +
        ggtitle(ds_id) +
        ylab("log2 counts normalized by region width\n and FRIP") +
        xlab("")
    })
  
  output_path <- fs::path(output_dir, "Individual_samples.pdf")
  all_samples_together %>% 
    marrangeGrob(ncol = 1, nrow = 2) %>% 
    ggsave(filename = output_path, width = 14, height = 10)
}

# Function to generate group aggregated plots
generate_group_aggregated_plots <- function(counts_norm, preset_colors, output_dir = "./") {
  by_groups <- counts_norm %>% 
    group_by(dataset_id) %>% 
    group_map(\(data, info) {
      ds_id <- info$dataset_id
      pp3 <- ggplot(data, aes(x = reorder(status, as.numeric(factor(color_order))), y = log2data, fill = status)) +
        geom_boxplot() +
        scale_fill_manual(name = "Group name", values = preset_colors) +
        theme_pubr(legend = "right") +
        scale_y_continuous(
          breaks = scales::pretty_breaks(10), 
          limits = c(min(counts_norm$log2data), max(counts_norm$log2data))) +
        theme(
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          panel.grid.major.y = element_line()
        ) +
        ggtitle("Keep all points for all samples in group for each region") +
        ylab("log2 counts normalized by region width\n and FRIP") +
        xlab("")
      
      pp4 <- ggplot(
        data %>% select(status, group_id, log2data_group_avg, color_order) %>% distinct(),
        aes(x = reorder(status, as.numeric(factor(color_order))), y = log2data_group_avg, fill = status)
      ) +
        geom_boxplot() +
        theme_pubr(legend = "right") +
        scale_fill_manual(name = "Group name", values = preset_colors) +
        scale_y_continuous(
          breaks = scales::pretty_breaks(10), 
          limits = c(min(counts_norm$log2data), max(counts_norm$log2data))) +
        theme(
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          panel.grid.major.y = element_line()
        ) +
        ggtitle("Average of all samples in group for each region") +
        ylab("log2 counts normalized by region width\n and FRIP") +
        xlab("")
      
      cowplot::plot_grid(
        pp3 + pp4 + plot_layout(guides = "collect", axis_titles = "collect_y") + 
          plot_annotation(title = ds_id)
      )
    })
  
  output_path <- fs::path(output_dir, "Individual_samples_aggregated.pdf")
  
  by_groups %>% 
    marrangeGrob(ncol = 1, nrow = 2) %>% 
    ggsave(filename = output_path, width = 14, height = 10)
}

# Function to generate pairwise comparison plots
generate_pairwise_comparison_plots <- function(counts_norm, config, output_dir = "./") {
  comparison_groups <- crossing(
    tibble(dataset_id = sort(unique(counts_norm$dataset_id))),
    config$comparisons
  )
  
  comparison_groups %>% 
    group_by(dataset_id) %>% 
    group_walk(\(data, info) {
      ds_id <- info$dataset_id
      output_path <- fs::path(output_dir, str_c(ds_id,".pdf"))
      
      pmap(data, \(num, group_treat, group_ctrl) {
        d <- counts_norm %>% 
          filter(dataset_id == ds_id & (group_id == group_treat | group_id == group_ctrl)) %>% 
          mutate(
            status_up = paste(status, sprintf("(%s)", group_id)),
            group_id = factor(group_id, levels = sort(unique(group_id))))
        
        # Calculate stats
        stats <- stats_for_one_dataset_id(d, num, group_ctrl, group_treat, ds_id) %>%
          pluck("short") 
        
        tab_title <- "*Statistics based on linear scale \n(no log2 transformation as on plots)" #%>% strwrap(30) %>% paste(collapse = "\n")
        gg_stats <- ggtexttable(stats, rows = NULL, theme = ttheme("classic", base_size = 8))  %>%
          tab_add_footnote(text = tab_title, size = 8, face = "italic") %>% 
          tab_add_hline(at.row = 5:6, row.side = "top", linewidth = 4)
        
        
        limits <- range(d$log2data, na.rm = TRUE)
        title <- str_glue("{num}. ({group_treat} vs {group_ctrl}) ({ds_id})")
        
        treat_color <- config$groups %>% filter(group_id == group_treat) %>% pull(color) %>% first()
        ctrl_color <- config$groups %>% filter(group_id == group_ctrl) %>% pull(color) %>% first()
        group_colors <- c(treat_color, ctrl_color)
        
        ## add group_order column to keep order for boxplots as treatment first, control second
        d <- d %>% mutate(group_order = case_when(group_id == group_treat ~ 1,
                                                  .default = 2))
        
        boxplots <- create_pairwise_plots_for_one_dataset(d, limits, title, group_colors)
        
        tmp <- boxplots$individual + boxplots$avg + boxplots$all_points 
        tmp <- tmp + gg_stats +
          plot_layout(axis_titles = 'collect', guides = 'collect', widths = c(2, 1, 1, 1)) &
          theme(legend.position = "bottom") &
          plot_annotation(title = title)
        
        cowplot::plot_grid(tmp)
      }) %>% 
        marrangeGrob(nrow = 3, ncol = 1) %>%
        ggsave(filename = output_path, width = 13, height = 16)
    })
}



# Function to generate all plots
generate_all_plots <- function(region_data, config_data, output_dir = "./", progress_reporter = NULL) {
  
  if (!is.null(progress_reporter)) {
    progress_reporter(value = 1, 
             message = "Generating plots 1 of 1")
  } 
  
  # Set factor levels and colors
  region_data$counts_norm$dataset_id <- factor(
    region_data$counts_norm$dataset_id,
    levels = sort(unique(region_data$counts_norm$dataset_id)))
  
  region_data$counts_norm$status <- factor(
    region_data$counts_norm$status,
    levels = unique(config_data$config$groups$group_desc))
  
  preset_colors <- config_data$config$groups %>%
    select(color, color_order) %>%
    distinct() %>%
    arrange(color_order) %>%
    pull(color)
  
  # Generate individual sample plots
  generate_individual_sample_plots(region_data$counts_norm, preset_colors, output_dir)
  
  # Generate group aggregated plots
  generate_group_aggregated_plots(region_data$counts_norm, preset_colors, output_dir)
  
  # Generate pairwise comparison plots
  generate_pairwise_comparison_plots(region_data$counts_norm, config_data$config, output_dir)
  
  
  
  return("plots_generated")
  
}

# ========================================
# REPORTS RELATED FUNCTIONS
# ========================================

stats_for_one_dataset_id <- function(df, num, ctrl_id, treat_id,  dataset_id){
  
  ### ALL DATA
  group_c <- df %>% filter(group_id == ctrl_id & dataset_id == dataset_id)
  group_t <- df %>% filter(group_id == treat_id & dataset_id == dataset_id)
  
  group_ctrl_flat <- group_c$counts_width_frip_norm
  group_treat_flat <- group_t$counts_width_frip_norm
  
  ### AVERAGE
  ### We want make analysis on not log2 data
  group_ctrl_avg <- group_c %>% 
    summarise(avg_counts_width_frip_norm = mean(counts_width_frip_norm), 
              .by = c("seqnames","start","end","status")) %>% 
    pull(avg_counts_width_frip_norm)
  
  group_treat_avg <- group_t %>% 
    summarise(avg_counts_width_frip_norm = mean(counts_width_frip_norm), 
              .by = c("seqnames","start","end","status")) %>% 
    pull(avg_counts_width_frip_norm)
  
  flat_wilcox_p <- wilcox.test(group_ctrl_flat, group_treat_flat, exact = FALSE)$p.value
  flat_ks_p <- ks.test(group_ctrl_flat, group_treat_flat)$p.value
  flat_ttest_p <- t.test(group_ctrl_flat, group_treat_flat)$p.value
  
  mean_wilcox_p <- wilcox.test(group_ctrl_avg, group_treat_avg, exact = FALSE)$p.value
  mean_ks_p <- ks.test(group_ctrl_avg, group_treat_avg)$p.value
  mean_ttest_p <- t.test(group_ctrl_avg, group_treat_avg)$p.value
  
  linear_avg_median_ratio <- round(median(group_treat_avg)/median(group_ctrl_avg),3)
  linear_flat_median_ratio <- round(median(group_treat_flat)/median(group_ctrl_flat),3)
  
  wide_res <- data.frame(
    dataset_id = dataset_id,
    treat = treat_id,
    ctrl = ctrl_id,
    compar_num = num,
    wilcox_flat = flat_wilcox_p,
    wilcox_mean = mean_wilcox_p,
    ks_flat = flat_ks_p,
    ks_mean = mean_ks_p,
    ttest_flat = flat_ttest_p,
    ttest_mean = mean_ttest_p,
    linear_avg_median_ratio = linear_avg_median_ratio,
    linear_flat_median_ratio = linear_flat_median_ratio
  )
  
  ratio_row <- data.frame(test = str_glue("Median ratio\n(linear scale)\n({treat_id}/{ctrl_id})"), 
                          average = as.character(linear_avg_median_ratio), 
                          all_points = as.character(linear_flat_median_ratio))
  
  short_res <- data.frame(
    test =c("wilcox", "KS", "t.test"),
    average = c(mean_wilcox_p,mean_ks_p,mean_ttest_p),
    all_points = c(flat_wilcox_p,flat_ks_p,flat_ttest_p)
  ) %>% mutate(all_points = case_when(all_points > 0.05 ~ "NS",
                                      .default = sprintf("%.2E", all_points)),
               average = case_when(average > 0.05 ~ "NS",
                                   .default = sprintf("%.2E", average))) %>% 
    bind_rows(ratio_row)
  
  list(short = short_res,
       wide = wide_res)
}

# Helper function to create comparison stats
create_comparison_stats <- function(counts_norm, config) {
  comparison_groups <- crossing(
    tibble(dataset_id = sort(unique(counts_norm$dataset_id))),
    config$comparisons
  )
  
  comparison_groups %>% 
    group_by(dataset_id) %>% 
    group_map(\(data, info) {
      ds_id <- info$dataset_id
      pmap(data, \(num, group_treat, group_ctrl) {
        d <- counts_norm %>% filter(dataset_id == ds_id & (group_id == group_treat | group_id == group_ctrl))
        stats_for_one_dataset_id(d, num, group_ctrl, group_treat, ds_id) %>% pluck("wide")
      }) %>% list_rbind()
    }) %>% 
    list_rbind() %>% 
    arrange(dataset_id, compar_num) %>% 
    select(
      dataset_id,
      treatment = treat,
      control = ctrl,
      comparison_num = compar_num,
      `Wilcox (average)` = wilcox_mean,
      `Wilcox (all points)` = wilcox_flat,
      `KS (average)` = ks_mean,
      `KS (all points)` = ks_flat,
      `t.test (average)` = ttest_mean,
      `t.test (all points)` = ttest_flat,
      `linear median ratio (average)` = linear_avg_median_ratio,
      `linear median ratio (all points)` = linear_flat_median_ratio
    )
}

# Function to generate report outputs
generate_report_outputs <- function(region_data, config_data, rippm_description, output_dir = "./", progress_reporter = NULL) {
  
  if (!is.null(progress_reporter)) {
    progress_reporter(value = 1, 
             message = "Generating reports 1 of 1")
  } 
  
  # Create pivoted data frames
  counts_norm_pivoted_raw <- create_pivoted_data(region_data$counts_norm, "records", config_data$config)
  counts_norm_pivoted_widthnorm <- create_pivoted_data(region_data$counts_norm, "records_width_norm", config_data$config)
  counts_norm_pivoted_frip <- create_pivoted_data(region_data$counts_norm, "counts_width_frip_norm", config_data$config)
  counts_norm_pivoted <- create_pivoted_data(region_data$counts_norm, "log2data", config_data$config)
  
  # Create comparison stats
  comparison_stats <- create_comparison_stats(region_data$counts_norm, config_data$config)
  
  # Create individual FRIPs data
  individual_frips <- rippm_description %>%
    mutate(average = round(mean(n_fragments_in_peaks_union), 1)) %>%
    select(
      sample_id,
      group_id = status,
      fragment_count = bam_total_n_fragments,
      fragment_in_peak_count = n_fragments_in_peaks_union,
      fragment_in_peak_ratio = nfrag_union_ratio,
      fragment_in_peak_average = average,
      normalization_factor = frip
    )
  
  output_path <- fs::path(output_dir, "Counts_before_and_after_normalization_FRIP.xlsx")
  
  # Write to Excel
  write_xlsx(
    list(
      RAW_COUNTS = counts_norm_pivoted_raw,
      RAW_AFTER_WIDTH_NORM = counts_norm_pivoted_widthnorm,
      AFTER_FRIP_NORM = counts_norm_pivoted_frip,
      FINAL_LOG2_AFTER_RIPPM_NORM = counts_norm_pivoted,
      INDIVIDUAL_FRIPS = individual_frips,
      COMPARISON_STATS = comparison_stats
    ),
    path = output_path,
    col_names = TRUE
  )
  
  return("reports_generated")
}

# Helper function to create pivoted data
create_pivoted_data <- function(counts_norm, value_col, config) {
  counts_norm %>% 
    distinct() %>% 
    select(region_id, dataset_id, chr = seqnames, start, end, width, sample_id, status, counts = !!value_col) %>%
    left_join(., config$groups %>% mutate(ord = row_number()) %>% select(sample_id, ord), join_by(sample_id)) %>% 
    mutate(caption = str_glue("{sample_id} ({status})")) %>%
    arrange(ord) %>%
    select(-c(status, sample_id, ord)) %>%
    mutate(counts = round(counts,3)) %>% 
    pivot_wider(names_from = caption, values_from = counts) %>% 
    arrange(dataset_id, chr, start, end)
}


# Main function to run the entire analysis pipeline
run_analysis_pipeline <- function(config_path, cache_path, progress_reporter) {
  # Load and process configuration
  config_data <- load_and_process_config(config_path = config_path, cache_path = cache_path, progress_reporter)
  
  # Get peak union and calculate FRIP
  peak_data <- get_peak_union_and_frip(config_data, progress_reporter)
  
  # Process regions and calculate read counts
  region_data <- process_regions_and_counts(config_data, peak_data, progress_reporter)
  
  # Generate all plots
  generate_all_plots(region_data, config_data, progress_reporter = progress_reporter)
  
  # Generate report outputs
  generate_report_outputs(region_data, config_data, peak_data$rippm_description, progress_reporter = progress_reporter)
}

# ========================================
# MAIN
# ========================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript script.R path_to_config_file.xlsx")
}

CONFIG_FILE_PATH <- args[1]
CACHE_PATH <- "/projectnb/wax-es/CHIPSEQ_BAM_CACHE/"
PROGRESS_FUNCTION <- function(value,message) {print(message)}

NCPU <- ifelse(is.na(args[2]), 1, as.numeric(args[2]))
print(paste0("NCPU: ", NCPU))
plan(multicore, workers = NCPU) 

run_analysis_pipeline(CONFIG_FILE_PATH, CACHE_PATH, PROGRESS_FUNCTION)



