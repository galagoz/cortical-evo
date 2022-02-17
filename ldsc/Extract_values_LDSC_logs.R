
library(tidyverse)

files <- list.files("/data/clusterfs/lag/users/gokala/enigma-evol/dti/data/munged/ldsc", full.names = T, include.dirs = F, pattern = ".log")

intercepts <- tibble("Region" = as.character(), 
                     "Total_h2" = as.numeric(),
                     "Total_h2_sterr" = as.numeric(),
                     "Lambda_GC" = as.numeric(),
                     "Mean_chisq" = as.numeric(),
                     "LDSC_intercept" = as.numeric(),
                     "LDSC_int_sterr" = as.numeric())

for (i in 1:length(files)) {
  results <- read_lines(files[i], skip = 25, n_max = 4)
  values <- as.numeric(unlist(str_extract_all(string = results, 
                                              pattern = "(?<=[:punct:]{1}\\s?)[:digit:]{1}\\.[:digit:]{2,4}")))
  info1 <- str_split(files[i], pattern = "/")
  Region = unlist(sapply(strsplit(info1[[1]][12], split = "_allChr_"), `[`, 1))
  to_add <- tibble("Region" = as.character(Region),
         "Total_h2" = as.numeric(values[1]),
         "Total_h2_sterr" = as.numeric(values[2]),
         "Lambda_GC" = as.numeric(values[3]),
         "Mean_chisq" = as.numeric(values[4]),
         "LDSC_intercept" = as.numeric(values[5]),
         "LDSC_int_sterr" = as.numeric(values[6]))
  intercepts <<- rbind(intercepts, to_add)
}

#write_csv(intercepts, "/data/clusterfs/lag/users/gokala/enigma-evol/dti/data/munged/ldsc/DTI_LDSC_intercepts.csv")

# Plot total h2

results_left = intercepts[grep(intercepts$Region, pattern = "left"),]
results_right = intercepts[grep(intercepts$Region, pattern = "right"),]
results_global = intercepts[-c(grep(intercepts$Region, pattern = "left"), grep(intercepts$Region, pattern = "right")),]

# plot global meaures
y_max1 <- max(results_global[, 2])
y_axis_max1 <- y_max1 + results_global[results_global$Total_h2 == y_max1, 3] + 0.02

pWM_global = ggplot(data = results_global, mapping = aes(Region, Total_h2)) +
  geom_bar(stat = "identity", position = "dodge", fill = "#ED6B06") + ##00786A
  geom_linerange(position = position_dodge(width = 0.9), aes(ymin = Total_h2 - Total_h2_sterr, ymax = Total_h2 + Total_h2_sterr)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "bottom") +
    labs(
    x = "Region",
    y = expression(paste("Total SNP-", italic("h"^{
      2
    }))),
    title = "Global FA measures"
    )

# plot left hemispheric meaures
y_max1 <- max(results_left[, 2])
y_axis_max1 <- y_max1 + results_leftl[results_left$Total_h2 == y_max1, 3] + 0.02

pWM_left = ggplot(data = results_left, mapping = aes(Region, Total_h2)) +
  geom_bar(stat = "identity", position = "dodge", fill = "#ED6B06") + ##00786A
  geom_linerange(position = position_dodge(width = 0.9), aes(ymin = Total_h2 - Total_h2_sterr, ymax = Total_h2 + Total_h2_sterr)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "bottom") +
  labs(
    x = "Region",
    y = expression(paste("Total SNP-", italic("h"^{
      2
    }))),
    title = "Left hemispheric FA measures"
    )

# plot right hemispheric meaures
y_max1 <- max(results_right[, 2])
y_axis_max1 <- y_max1 + results_right[results_right$Total_h2 == y_max1, 3] + 0.02

pWM_right = ggplot(data = results_right, mapping = aes(Region, Total_h2)) +
  geom_bar(stat = "identity", position = "dodge", fill = "#ED6B06") + ##00786A
  geom_linerange(position = position_dodge(width = 0.9), aes(ymin = Total_h2 - Total_h2_sterr, ymax = Total_h2 + Total_h2_sterr)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "bottom") +
  labs(
    x = "Region",
    y = expression(paste("Total SNP-", italic("h"^{
      2
    }))),
    title = "Right hemispheric FA measures"
    )

plot_list = list()
plot_list = plot_grid(pWM_global, pWM_left, pWM_right,
                      labels = c("A", "B", "C"),
                      ncol = 1, nrow = 3)

ggsave("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/results/partitioned_heritability/final_results/dti/WM_total_h2_barplots.pdf", plot_list)
