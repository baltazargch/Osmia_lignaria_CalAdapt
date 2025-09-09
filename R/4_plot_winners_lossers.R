library(tidyverse)
library(scales)
library(gt)

res <-  read_csv('outputs/final_maxent/csv/final_models_metrics.csv')

prres <- read_csv('outputs/models/maxent/log_models_tunning.csv')

table(res$fc)

res$npseu %>% table()
nrow(res)

prres$or.10p.avg %>% summary()
res$auc_train %>% summary

res$current_area %>% summary

out <- res %>% 
  select(current_area:last_col()) %>% 
  mutate(
    across(!contains('current|species'), \(x) (x - current_area)/current_area*100)
  )

scenarios <- grep("^\\d{4}-\\d{4}_ssp\\d+$", names(res), value = TRUE)

# summary table: (current_area - scenario) for each scenario
summ_diff <- map_dfr(scenarios, function(scn) {
  diffs <- (res[[scn]] - res$current_area)/res$current_area*100
  s <- summary(diffs)  # Min., 1st Qu., Median, Mean, 3rd Qu., Max (and NA's)
  tibble(
    scenario = scn,
    Min      = unname(s["Min."]),
    Q1       = unname(s["1st Qu."]),
    Median   = unname(s["Median"]),
    Mean     = unname(s["Mean"]),
    Q3       = unname(s["3rd Qu."]),
    Max      = unname(s["Max."]),
    NAs      = if ("NA's" %in% names(s)) unname(s["NA's"]) else 0L
  )
})

summ_diff %>% select(-NAs)

bin_breaks <- c(-Inf,-75,-50,-25,-10,10,25,50,75, Inf)
bin_labels <- c("≤ -75%","-75% to -50%","-50% to -25%","-25% to -10%",
                "-10% to 10%","10% to 25%","25% to 50%","50% to 75%","≥ 75%")

prep_counts <- function(df, col){
  df %>%
    mutate(change_bin = cut(.data[[col]], breaks = bin_breaks,
                            labels = bin_labels, ordered_result = TRUE)) %>%
    mutate(change_bin = forcats::fct_expand(change_bin, bin_labels)) %>%
    count(change_bin, name = "n_species", .drop = FALSE) %>% 
    mutate(Case = ifelse(str_detect(change_bin, '-'), 'Loss', 'Gain')) %>% 
    mutate(Case = ifelse(str_detect(change_bin, '-10% to'), 'No Change', Case)) %>% 
    mutate(Case = factor(Case, levels = c('Loss', 'No Change', 'Gain')))
}

ssp585_2015 <- prep_counts(out, "2075-2100_ssp245") 
ssp585_2075 <- prep_counts(out, "2075-2100_ssp585")

ssp585a <- ssp585_2015 %>% 
  filter(!is.na(Case)) %>% 
  ggplot(aes(x = change_bin, y = n_species, fill = Case)) +
  geom_col(position = "dodge") +
  # geom_vline(xintercept = 5, linetype = "dashed", color = "gray12") +
  envalysis::theme_publish(base_size = 14) +
  ggsci::scale_fill_npg(alpha = 0.75)+
  labs(x = "% Change bin", y = "Number of species",  title = "SSP245 2075-2100") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ssp585b <- ssp585_2075 %>% 
  filter(!is.na(Case)) %>% 
  ggplot(aes(x = change_bin, y = n_species, fill = Case)) +
  geom_col(position = "dodge") +
  # geom_vline(xintercept = 5, linetype = "dashed", color = "gray12") +
  envalysis::theme_publish(base_size = 14) +
  ggsci::scale_fill_npg(alpha = 0.75)+
  labs(x = "% Change bin", y = "Number of species",  title = "SSP585 2075-2100") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

library(patchwork)

ssp585a / ssp585b + 
  plot_layout(guides = 'collect') &
  plot_annotation(tag_levels = 'a') &
  theme(legend.position = 'bottom')

ggsave(filename = 'figures/win_los.png', 
       dpi=600,
       width = 14, height = 16, scale=1.5, units = 'cm')
