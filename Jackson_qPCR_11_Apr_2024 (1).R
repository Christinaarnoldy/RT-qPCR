### qPCR analysis
### Rob Jackson (Van Doorslaer Lab)
### Created 11 Apr 2024
### Updated 12 Apr 2024

### RStudio 2023.12.1+402 "Ocean Storm"
### R version 4.3.2 "Eye Holes"

### Housekeeping
### Install packages as needed via install.packages(), BiocManager::install(), or remotes::install_github

### Set working directory to the R script location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

### Tutorial: https://liz-is.github.io/qpcr-analysis-with-r/01-import-data/index.html

### Load packages
library(tidyverse)
library(readxl)
library(ggprism)
library(cowplot)
library(ggpubr)
library(ggrepel)
library(dplyr)  

### Import data
qpcr_data <- read_excel("CA_HTK_Raft_1.xlsx", sheet = "Results", skip = 52, n_max = 235)
head(qpcr_data)
tail(qpcr_data)
tidy_data <- separate(qpcr_data, "Well Position", into = c("row", "column"), sep = 1)
head(tidy_data)
tail(tidy_data)

ggplot(tidy_data, aes(x = reorder(column, as.numeric(column)), y = forcats::fct_rev(row), fill = `Target Name`, label = `Sample Name`)) +
  geom_tile(colour = "black") +
  geom_text()

host_data <- filter(tidy_data, `Target Name` != "E6*" & `Target Name` != "E1^E4" & `Sample Name` != "NTC")

summarised_host_data <- host_data %>%
  mutate(CT = as.numeric(CT)) %>%
  group_by(`Sample Name`, `Target Name`) %>%
  summarise(mean_Ct = mean(CT))
summarised_host_data

ggplot(summarised_host_data, aes(x = `Sample Name`, y = mean_Ct, colour = `Target Name`)) +
  geom_point() + theme_prism(axis_text_angle = 45) + 
  scale_y_continuous(guide = "prism_minor", limits = c(0, 40), expand = c(0, 0))

C7_data <- summarised_host_data %>%
  filter(`Target Name` == "CEACAM7") %>%
  rename("C7_Ct" = "mean_Ct")
C7_data

ref_data <- summarised_host_data %>%
  filter(`Target Name` == "TBP") %>%
  rename("TBP_Ct" = "mean_Ct")
ref_data

combined_data <- left_join(C7_data, ref_data, by = "Sample Name")
combined_data
combined_data <- mutate(combined_data, delta_Ct = TBP_Ct - mean_Ct)

ggplot(combined_data, aes(x = `Sample Name`, y = delta_Ct)) +
  geom_point() + theme_prism(axis_text_angle = 45) + 
  scale_y_continuous(guide = "prism_minor", limits = c(0, 6), expand = c(0, 0))

mean_control <- filter(combined_data, `Sample Name` == "2_T2Y-DMSO") %>% pull(delta_Ct)
mean_control

combined_data <- combined_data %>% 
  mutate(delta_delta_Ct = mean_control - delta_Ct)

ggplot(combined_data, aes(x = `Sample Name`, y = delta_delta_Ct)) +
  geom_point() + theme_prism(axis_text_angle = 45) + 
  scale_y_continuous(guide = "prism_minor")

combined_data <- combined_data %>%
  mutate(rel_conc = 2^-delta_delta_Ct)

combined_data <- combined_data %>%
  mutate(sample_num = str_extract(`Sample Name`, "^[^_]+(?=_)")) %>%
  mutate(sample_group = str_extract(`Sample Name`, "(?<=_)[^_]+(?=-)")) %>%
  mutate(treatment = str_extract(`Sample Name`, "(?<=-).*"))
combined_data

ggplot(combined_data, aes(x = reorder(`Sample Name`, as.numeric(sample_num)), y = rel_conc)) +
  geom_point() + theme_prism(axis_text_angle = 45) + 
  scale_y_continuous(guide = "prism_minor", limits = c(0, 50), expand = c(0, 0)) + ylab("CEACAM7 mRNA FC") +
  geom_hline(yintercept = 1, linetype = 2, linewidth = 1) + xlab("Sample")

p1 <- ggplot(combined_data, aes(fill = treatment, x = sample_group, y = rel_conc)) +
  geom_bar(position = "dodge", stat = "identity", color = "black") + theme_prism(axis_text_angle = 0, border= T) + scale_fill_prism(palette = "floral") +
  scale_y_continuous(guide = "prism_minor", limits = c(0, 50), expand = c(0, 0)) + ylab("CEACAM7 mRNA FC") +
  geom_hline(yintercept = 1, linetype = 2) + xlab("Sample") + guides(colour = guide_legend(position = "inside")) +
  theme(legend.position = c(0.9, 0.8)) +
  coord_cartesian(clip = "off")
p1



IVL_data <- summarised_host_data %>%
  filter(`Target Name` == "IVL")%>%
  rename("IVL_Ct" = "mean_Ct")
IVL_data

combined_data <- left_join(IVL_data, ref_data, by = "Sample Name")
combined_data
combined_data <- mutate(combined_data, delta_Ct = TBP_Ct - mean_Ct)
combined_data

ggplot(combined_data, aes(x = `Sample Name`, y = delta_Ct)) +
  geom_point() + theme_prism(axis_text_angle = 45) + 
  scale_y_continuous(guide = "prism_minor", limits = c(0, 10), expand = c(0, 0))

mean_control <- filter(combined_data, `Sample Name` == "2_T2Y-DMSO") %>% pull(delta_Ct)
mean_control

combined_data <- combined_data %>% 
  mutate(delta_delta_Ct = mean_control - delta_Ct)

ggplot(combined_data, aes(x = `Sample Name`, y = delta_delta_Ct)) +
  geom_point() + theme_prism(axis_text_angle = 45) + 
  scale_y_continuous(guide = "prism_minor")

combined_data <- combined_data %>%
  mutate(rel_conc = 2^-delta_delta_Ct)

combined_data <- combined_data %>%
  mutate(sample_num = str_extract(`Sample Name`, "^[^_]+(?=_)")) %>%
  mutate(sample_group = str_extract(`Sample Name`, "(?<=_)[^_]+(?=-)")) %>%
  mutate(treatment = str_extract(`Sample Name`, "(?<=-).*"))
combined_data

ggplot(combined_data, aes(x = reorder(`Sample Name`, as.numeric(sample_num)), y = rel_conc)) +
  geom_point() + theme_prism(axis_text_angle = 45) + 
  scale_y_continuous(guide = "prism_minor", limits = c(0, 25), expand = c(0, 0)) + ylab("IVL mRNA FC") +
  geom_hline(yintercept = 1, linetype = 2, linewidth = 1) + xlab("Sample")

p2 <- ggplot(combined_data, aes(fill = treatment, x = sample_group, y = rel_conc)) +
  geom_bar(position = "dodge", stat = "identity", color = "black") + theme_prism(axis_text_angle = 0, border= T) + scale_fill_prism(palette = "floral") +
  scale_y_continuous(guide = "prism_minor", limits = c(0, 25), expand = c(0, 0)) + ylab("IVL mRNA FC") +
  geom_hline(yintercept = 1, linetype = 2) + xlab("Sample") + guides(colour = guide_legend(position = "inside")) +
  theme(legend.position = c(0.1, 0.8)) +
  coord_cartesian(clip = "off")
p2


PPARA_data <- summarised_host_data %>%
  filter(`Target Name` == "PPARa") %>%
  rename("PPARA_Ct" = "mean_Ct")
PPARA_data

combined_data <- left_join(PPARA_data, ref_data, by = "Sample Name")
combined_data
combined_data <- mutate(combined_data, delta_Ct = TBP_Ct - mean_Ct)
combined_data

ggplot(combined_data, aes(x = `Sample Name`, y = delta_Ct)) +
  geom_point() + theme_prism(axis_text_angle = 45) + 
  scale_y_continuous(guide = "prism_minor", limits = c(-4, 0), expand = c(0, 0))

mean_control <- filter(combined_data, `Sample Name` == "2_T2Y-DMSO") %>% pull(delta_Ct)
mean_control

combined_data <- combined_data %>% 
  mutate(delta_delta_Ct = mean_control - delta_Ct)

ggplot(combined_data, aes(x = `Sample Name`, y = delta_delta_Ct)) +
  geom_point() + theme_prism(axis_text_angle = 45) + 
  scale_y_continuous(guide = "prism_minor")

combined_data <- combined_data %>%
  mutate(rel_conc = 2^-delta_delta_Ct)

combined_data <- combined_data %>%
  mutate(sample_num = str_extract(`Sample Name`, "^[^_]+(?=_)")) %>%
  mutate(sample_group = str_extract(`Sample Name`, "(?<=_)[^_]+(?=-)")) %>%
  mutate(treatment = str_extract(`Sample Name`, "(?<=-).*"))
combined_data

ggplot(combined_data, aes(x = reorder(`Sample Name`, as.numeric(sample_num)), y = rel_conc)) +
  geom_point() + theme_prism(axis_text_angle = 45) + 
  scale_y_continuous(guide = "prism_minor", limits = c(0, 10), expand = c(0, 0)) + ylab("PPARA mRNA FC") +
  geom_hline(yintercept = 1, linetype = 2, linewidth = 1) + xlab("Sample")

p3 <- ggplot(combined_data, aes(fill = treatment, x = sample_group, y = rel_conc)) +
  geom_bar(position = "dodge", stat = "identity", color = "black") + theme_prism(axis_text_angle = 0, border= T) + scale_fill_prism(palette = "floral") +
  scale_y_continuous(guide = "prism_minor", limits = c(0, 10), breaks = seq(0, 10, 2), minor_breaks = seq(0, 10, 1), expand = c(0, 0)) + ylab("PPARA mRNA FC") +
  geom_hline(yintercept = 1, linetype = 2) + xlab("Sample") + guides(colour = guide_legend(position = "inside")) +
  theme(legend.position = c(0.9, 0.85)) +
  coord_cartesian(clip = "off")
p3



viral_data <- filter(tidy_data, `Target Name` != "IVL" & `Target Name` != "CEACAM7" & `Target Name` != "PPARa" & 
                       `Sample Name` != "NTC" &
                       `Sample Name` != "1_T2Y-GW" &
                       `Sample Name` != "2_T2Y-DMSO" &
                       `Sample Name` != "3_T2Y-Feno" &
                       `Sample Name` != "4_T3Y-GW" &
                       `Sample Name` != "5_T3Y-DMSO" &
                       `Sample Name` != "6_T3Y-Feno")
viral_data

summarised_viral_data <- viral_data %>%
  mutate(CT = as.numeric(CT)) %>%
  group_by(`Sample Name`, `Target Name`) %>%
  summarise(mean_Ct = mean(CT))
summarised_viral_data

ggplot(summarised_viral_data, aes(x = `Sample Name`, y = mean_Ct, colour = `Target Name`)) +
  geom_point() + theme_prism(axis_text_angle = 45) + 
  scale_y_continuous(guide = "prism_minor", limits = c(0, 40), expand = c(0, 0))

E6_data <- summarised_viral_data %>%
  filter(`Target Name` == "E6*") %>%
  rename("E6_Ct" = "mean_Ct")
E6_data

ref_data <- summarised_viral_data %>%
  filter(`Target Name` == "TBP") %>%
  rename("TBP_Ct" = "mean_Ct")
ref_data

combined_data <- left_join(E6_data, ref_data, by = "Sample Name")
combined_data
combined_data <- mutate(combined_data, delta_Ct = TBP_Ct - E6_Ct)
combined_data

ggplot(combined_data, aes(x = `Sample Name`, y = delta_Ct)) +
  geom_point() + theme_prism(axis_text_angle = 45) + 
  scale_y_continuous(guide = "prism_minor", limits = c(0, 7), expand = c(0, 0))

mean_control <- filter(combined_data, `Sample Name` == "8_T216-DMSO") %>% pull(delta_Ct)
mean_control

combined_data <- combined_data %>% 
  mutate(delta_delta_Ct = mean_control - delta_Ct)

ggplot(combined_data, aes(x = `Sample Name`, y = delta_delta_Ct)) +
  geom_point() + theme_prism(axis_text_angle = 45) + 
  scale_y_continuous(guide = "prism_minor")

combined_data <- combined_data %>%
  mutate(rel_conc = 2^-delta_delta_Ct)

combined_data <- combined_data %>%
  mutate(sample_num = str_extract(`Sample Name`, "^[^_]+(?=_)")) %>%
  mutate(sample_group = str_extract(`Sample Name`, "(?<=_)[^_]+(?=-)")) %>%
  mutate(treatment = str_extract(`Sample Name`, "(?<=-).*"))
combined_data

ggplot(combined_data, aes(x = reorder(`Sample Name`, as.numeric(sample_num)), y = rel_conc)) +
  geom_point() + theme_prism(axis_text_angle = 45) + 
  scale_y_continuous(guide = "prism_minor", limits = c(0, 2), expand = c(0, 0)) + ylab("E6* mRNA FC") +
  geom_hline(yintercept = 1, linetype = 2, linewidth = 1) + xlab("Sample")

p4 <- ggplot(combined_data, aes(fill = treatment, x = sample_group, y = rel_conc)) +
  geom_bar(position = "dodge", stat = "identity", color = "black") + theme_prism(axis_text_angle = 0, border= T) + scale_fill_prism(palette = "floral") +
  scale_y_continuous(guide = "prism_minor", limits = c(0, 2), expand = c(0, 0)) + ylab("E6* mRNA FC") +
  geom_hline(yintercept = 1, linetype = 2) + xlab("Sample") + guides(colour = guide_legend(position = "inside")) +
  theme(legend.position = c(0.2, 0.8)) +
  coord_cartesian(clip = "off")
p4




E1E4_data <- summarised_viral_data %>%
  filter(`Target Name` == "E1^E4") %>%
rename("E1E4_Ct" = "mean_Ct")
E1E4_data

combined_data <- left_join(E1E4_data, ref_data, by = "Sample Name")
combined_data
combined_data <- mutate(combined_data, delta_Ct = TBP_Ct - mean_Ct)
combined_data

ggplot(combined_data, aes(x = `Sample Name`, y = delta_Ct)) +
  geom_point() + theme_prism(axis_text_angle = 45) + 
  scale_y_continuous(guide = "prism_minor", limits = c(0, 7), expand = c(0, 0))

mean_control <- filter(combined_data, `Sample Name` == "8_T216-DMSO") %>% pull(delta_Ct)
mean_control

combined_data <- combined_data %>% 
  mutate(delta_delta_Ct = mean_control - delta_Ct)

ggplot(combined_data, aes(x = `Sample Name`, y = delta_delta_Ct)) +
  geom_point() + theme_prism(axis_text_angle = 45) + 
  scale_y_continuous(guide = "prism_minor")

combined_data <- combined_data %>%
  mutate(rel_conc = 2^-delta_delta_Ct)

combined_data <- combined_data %>%
  mutate(sample_num = str_extract(`Sample Name`, "^[^_]+(?=_)")) %>%
  mutate(sample_group = str_extract(`Sample Name`, "(?<=_)[^_]+(?=-)")) %>%
  mutate(treatment = str_extract(`Sample Name`, "(?<=-).*"))
combined_data

ggplot(combined_data, aes(x = reorder(`Sample Name`, as.numeric(sample_num)), y = rel_conc)) +
  geom_point() + theme_prism(axis_text_angle = 45) + 
  scale_y_continuous(guide = "prism_minor", limits = c(0, 2), expand = c(0, 0)) + ylab("E1^E4 mRNA FC") +
  geom_hline(yintercept = 1, linetype = 2, linewidth = 1) + xlab("Sample")

p5 <- ggplot(combined_data, aes(fill = treatment, x = sample_group, y = rel_conc)) +
  geom_bar(position = "dodge", stat = "identity", color = "black") + theme_prism(axis_text_angle = 0, border= T) + scale_fill_prism(palette = "floral") +
  scale_y_continuous(guide = "prism_minor", limits = c(0, 2), expand = c(0, 0)) + ylab("E1^E4 mRNA FC") +
  geom_hline(yintercept = 1, linetype = 2) + xlab("Sample") + guides(colour = guide_legend(position = "inside")) +
  theme(legend.position = c(0.2, 0.8)) +
  coord_cartesian(clip = "off")
p5

plot_all <- plot_grid(p1+theme_prism(base_size = 8)+theme(legend.position="none"),
          p2+theme_prism(base_size = 8)+theme(legend.position="none"),
          p3+theme_prism(base_size = 8)+theme(legend.position="none"),
          p4+theme_prism(base_size = 8)+theme(legend.position="none"),
          p5+theme_prism(base_size = 8)+theme(legend.position=c(1.5, 0.5)),
          labels = c("A", "B", "C", "D", "E"),
          label_size = 12)
plot_all
save_plot("CA_HTK_Raft_1.pdf", plot_all)






### 2nd set
qpcr_data_2 <- read_excel("CA_HTK_Raft_2_corrected.xlsx", sheet = "Results", skip = 52, n_max = 276) #n_max = total number after the skip
head(qpcr_data_2)
tail(qpcr_data_2)

#split well position into row and column
tidy_data_2 <- separate(qpcr_data_2, "Well Position", into = c("row", "column"), sep = 1)
head(tidy_data_2)
tail(tidy_data_2)

#plot the plate grid
ggplot(tidy_data_2, aes(x = reorder(column, as.numeric(column)), y = forcats::fct_rev(row), fill = `Target Name`, label = `Sample Name`)) +
  geom_tile(colour = "black") +
  geom_text()

#choose just the host genes
host_data_2 <- filter(tidy_data_2, `Target Name` != "E6*" & `Target Name` != "E1^E4" & `Sample Name` != "NTC")

#make CT a number
#group by sample and target name
#calculate the mean CT
summarised_host_data_2 <- host_data_2 %>%
  mutate(CT = as.numeric(CT)) %>%
  group_by(`Sample Name`, `Target Name`) %>%
  summarise(mean_Ct = mean(CT))
summarised_host_data_2

#plot the mean CTs of host genes
ggplot(summarised_host_data_2, aes(x = reorder(`Sample Name`, as.numeric(`Sample Name`)), y = mean_Ct, colour = `Target Name`)) +
  geom_point() + theme_prism(axis_text_angle = 45) + 
  scale_y_continuous(guide = "prism_minor", limits = c(0, 40), expand = c(0, 0)) + xlab("Sample")

#focus on ADD3, pull all the runs, rename the mean_Ct to ADD3_Ct
ADD3_data_2 <- summarised_host_data_2 %>%
  filter(`Target Name` == "ADD3") %>%
  arrange(as.numeric(`Sample Name`)) %>%
  rename("ADD3_Ct" = "mean_Ct")
ADD3_data_2

#do the same for the reference gene TBP
ref_data_2 <- summarised_host_data_2 %>%
  filter(`Target Name` == "TBP") %>%
  rename("TBP_Ct" = "mean_Ct") %>%
  arrange(as.numeric(`Sample Name`))
ref_data_2

#merge the ADD3 and TBP tables, calc dCT
combined_data_2 <- left_join(ADD3_data_2, ref_data_2, by = "Sample Name")
combined_data_2
combined_data_2 <- mutate(combined_data_2, ADD3_delta_Ct = TBP_Ct - ADD3_Ct)
combined_data_2


ggplot(combined_data_2, aes(x = reorder(`Sample Name`, as.numeric(`Sample Name`)), y = ADD3_delta_Ct)) +
  geom_point() + theme_prism(axis_text_angle = 0) + 
  scale_y_continuous(guide = "prism_minor", limits = c(0, 6), expand = c(0, 0)) + xlab("Sample")

#choose the sample to normalize on, 2.1 in this case
ADD3_mean_control <- filter(combined_data_2, `Sample Name` == "2.1") %>% pull(ADD3_delta_Ct)
ADD3_mean_control

#2^deltadeltaCT to get the fold change
combined_data_2 <- combined_data_2 %>% 
  mutate(ADD3_delta_delta_Ct = ADD3_mean_control - ADD3_delta_Ct)

ggplot(combined_data_2, aes(x = reorder(`Sample Name`, as.numeric(`Sample Name`)), y = ADD3_delta_delta_Ct)) +
  geom_point() + theme_prism(axis_text_angle = 0) + 
  scale_y_continuous(guide = "prism_minor") + xlab("Sample")

combined_data_2 <- combined_data_2 %>%
  mutate(ADD3_rel_conc = 2^-ADD3_delta_delta_Ct)

ggplot(combined_data_2, aes(x = reorder(`Sample Name`, as.numeric(`Sample Name`)), y = ADD3_rel_conc)) +
  geom_point() + theme_prism(axis_text_angle = 0) + 
  scale_y_continuous(guide = "prism_minor", limits = c(0, 3), expand = c(0, 0)) + ylab("ADD3 mRNA FC") +
  geom_hline(yintercept = 1, linetype = 2, linewidth = 1) + xlab("Sample")

combined_data_2$sample_group = c(rep("T2Y", 3), rep("T3Y", 3), "T216_1", "T216_2", "T216_1", "T216_2", "T216_1", "T216_2", "T316_1", "T316_2", "T316_1", "T316_2", "T316_1", "T316_2")
combined_data_2$treatment =c("GW", "DMSO", "Feno", "GW", "DMSO", "Feno", "GW", "GW", "DMSO", "DMSO", "Feno", "Feno", "GW", "GW", "DMSO", "DMSO", "Feno", "Feno")

p1 <- ggplot(combined_data_2, aes(fill = treatment, x = sample_group, y = ADD3_rel_conc)) +
  geom_bar(position = "dodge", stat = "identity", color = "black") + theme_prism(axis_text_angle = 0, border= T) + scale_fill_prism(palette = "floral") +
  scale_y_continuous(guide = "prism_minor", limits = c(0, 3), expand = c(0, 0)) + ylab("ADD3 mRNA FC") +
  geom_hline(yintercept = 1, linetype = 2) + xlab("Sample") + guides(colour = guide_legend(position = "inside")) +
  theme(legend.position = c(0.8, 0.8)) +
  coord_cartesian(clip = "off")
p1 ### doesn't match CA slides?



#do it again for LYPD2
LYPD2_data <- summarised_host_data_2 %>%
  filter(`Target Name` == "LYPD2") %>%
  arrange(as.numeric(`Sample Name`)) %>%
  rename("LYPD2_Ct" = "mean_Ct")
LYPD2_data

combined_data_2 <- left_join(combined_data_2, LYPD2_data, by = "Sample Name")
#View(combined_data_2)

combined_data_2 <- mutate(combined_data_2, LYPD2_delta_Ct = TBP_Ct - LYPD2_Ct)
#View(combined_data_2)

ggplot(combined_data_2, aes(x = reorder(`Sample Name`, as.numeric(`Sample Name`)), y = LYPD2_delta_Ct)) +
  geom_point() + theme_prism(axis_text_angle = 0) + 
  scale_y_continuous(guide = "prism_minor", limits = c(0, 7), expand = c(0, 0)) + xlab("Sample")

LYPD2_mean_control <- filter(combined_data_2, `Sample Name` == "2.1") %>% pull(LYPD2_delta_Ct)
LYPD2_mean_control

combined_data_2 <- combined_data_2 %>% 
  mutate(LYPD2_delta_delta_Ct = LYPD2_mean_control - LYPD2_delta_Ct)

ggplot(combined_data_2, aes(x = reorder(`Sample Name`, as.numeric(`Sample Name`)), y = LYPD2_delta_delta_Ct)) +
  geom_point() + theme_prism(axis_text_angle = 0) + 
  scale_y_continuous(guide = "prism_minor") + xlab("Sample")

combined_data_2 <- combined_data_2 %>%
  mutate(LYPD2_rel_conc = 2^-LYPD2_delta_delta_Ct)

ggplot(combined_data_2, aes(x = reorder(`Sample Name`, as.numeric(`Sample Name`)), y = LYPD2_rel_conc)) +
  geom_point() + theme_prism(axis_text_angle = 0) + 
  scale_y_continuous(guide = "prism_minor", limits = c(0, 8), expand = c(0, 0)) + ylab("LYPD2 mRNA FC") +
  geom_hline(yintercept = 1, linetype = 2, linewidth = 1) + xlab("Sample")

p2 <- ggplot(combined_data_2, aes(fill = treatment, x = sample_group, y = LYPD2_rel_conc)) +
  geom_bar(position = "dodge", stat = "identity", color = "black") + theme_prism(axis_text_angle = 0, border= T) + scale_fill_prism(palette = "floral") +
  scale_y_continuous(guide = "prism_minor", limits = c(0, 8), expand = c(0, 0)) + ylab("LYPD2 mRNA FC") +
  geom_hline(yintercept = 1, linetype = 2) + xlab("Sample") + guides(colour = guide_legend(position = "inside")) +
  theme(legend.position = c(0.8, 0.8)) +
  coord_cartesian(clip = "off")
p2 


plot_all_2 <- plot_grid(p1+theme_prism(base_size = 8)+theme(legend.position=c(0.9, 0.9)),
                      p2+theme_prism(base_size = 8)+theme(legend.position="none"),
                      labels = c("A", "B"),
                      label_size = 12)
save_plot("CA_HTK_Raft_2.pdf", plot_all_2)


#are ADD3 and LYPD2 inversely correlated? (they should be)
ggplot(combined_data_2, aes(x=-ADD3_delta_delta_Ct, y=-LYPD2_delta_delta_Ct)) + geom_point() + geom_smooth(method=lm, fullrange = T) +
  stat_cor(aes(label=..rr.label..)) + geom_text_repel(aes(label=`Sample Name`)) + theme_prism(border = T) +
  scale_y_continuous(limits = c(-4, 4), guide = "prism_minor", expand = c(0, 0)) + scale_x_continuous(limits = c(-1, 1.5), guide = "prism_minor", expand = c(0, 0)) +coord_cartesian(clip = "off") + 
  geom_hline(yintercept=0) + geom_vline(xintercept=0)
#log2FC



PPARA_data_2 <- summarised_host_data_2 %>%
  filter(`Target Name` == "PPARa") %>%
  arrange(as.numeric(`Sample Name`)) %>%
  rename("PPARA_Ct" = "mean_Ct")
PPARA_data_2

combined_data_2 <- left_join(combined_data_2, PPARA_data_2, by = "Sample Name")
#View(combined_data_2)

combined_data_2 <- mutate(combined_data_2, PPARA_delta_Ct = TBP_Ct - PPARA_Ct)

ggplot(combined_data_2, aes(x = reorder(`Sample Name`, as.numeric(`Sample Name`)), y = PPARA_delta_Ct)) +
  geom_point() + theme_prism(axis_text_angle = 0) + 
  scale_y_continuous(guide = "prism_minor", limits = c(-2, 2), expand = c(0, 0)) + xlab("Sample")

PPARA_mean_control <- filter(combined_data_2, `Sample Name` == "8.2") %>% pull(PPARA_delta_Ct)
PPARA_mean_control

combined_data_2 <- combined_data_2 %>% 
  mutate(PPARA_delta_delta_Ct = PPARA_mean_control - PPARA_delta_Ct)

ggplot(combined_data_2, aes(x = reorder(`Sample Name`, as.numeric(`Sample Name`)), y = PPARA_delta_delta_Ct)) +
  geom_point() + theme_prism(axis_text_angle = 0) + 
  scale_y_continuous(guide = "prism_minor") + xlab("Sample")

combined_data_2 <- combined_data_2 %>%
  mutate(PPARA_rel_conc = 2^-PPARA_delta_delta_Ct)

ggplot(combined_data_2, aes(x = reorder(`Sample Name`, as.numeric(`Sample Name`)), y = PPARA_rel_conc)) +
  geom_point() + theme_prism(axis_text_angle = 0) + 
  scale_y_continuous(guide = "prism_minor", limits = c(0, 2), expand = c(0, 0)) + ylab("PPARA mRNA FC") +
  geom_hline(yintercept = 1, linetype = 2, linewidth = 1) + xlab("Sample")

p3 <- ggplot(combined_data_2, aes(fill = treatment, x = sample_group, y = PPARA_rel_conc)) +
  geom_bar(position = "dodge", stat = "identity", color = "black") + theme_prism(axis_text_angle = 0, border= T) + scale_fill_prism(palette = "floral") +
  scale_y_continuous(guide = "prism_minor", limits = c(0, 2), expand = c(0, 0)) + ylab("PPARA mRNA FC") +
  geom_hline(yintercept = 1, linetype = 2) + xlab("Sample") + guides(colour = guide_legend(position = "inside")) +
  theme(legend.position = c(0.8, 0.8)) +
  coord_cartesian(clip = "off")
p3 #looks good



C7_data_2 <- summarised_host_data_2 %>%
  filter(`Target Name` == "CEACAM7") %>%
  arrange(as.numeric(`Sample Name`)) %>%
  rename("C7_Ct" = "mean_Ct")
C7_data_2

combined_data_2 <- left_join(combined_data_2, C7_data_2, by = "Sample Name")
#View(combined_data_2)

combined_data_2 <- mutate(combined_data_2, C7_delta_Ct = TBP_Ct - C7_Ct)

ggplot(combined_data_2, aes(x = reorder(`Sample Name`, as.numeric(`Sample Name`)), y = C7_delta_Ct)) +
  geom_point() + theme_prism(axis_text_angle = 0) + 
  scale_y_continuous(guide = "prism_minor", limits = c(-1, 6), expand = c(0, 0)) + xlab("Sample")

C7_mean_control <- filter(combined_data_2, `Sample Name` == "8.2") %>% pull(C7_delta_Ct)
C7_mean_control

combined_data_2 <- combined_data_2 %>% 
  mutate(C7_delta_delta_Ct = C7_mean_control - C7_delta_Ct)

ggplot(combined_data_2, aes(x = reorder(`Sample Name`, as.numeric(`Sample Name`)), y = C7_delta_delta_Ct)) +
  geom_point() + theme_prism(axis_text_angle = 0) + 
  scale_y_continuous(guide = "prism_minor") + xlab("Sample")

combined_data_2 <- combined_data_2 %>%
  mutate(C7_rel_conc = 2^-C7_delta_delta_Ct)

ggplot(combined_data_2, aes(x = reorder(`Sample Name`, as.numeric(`Sample Name`)), y = C7_rel_conc)) +
  geom_point() + theme_prism(axis_text_angle = 0) + 
  scale_y_continuous(guide = "prism_minor", limits = c(0, 2), expand = c(0, 0)) + ylab("CEACAM7 mRNA FC") +
  geom_hline(yintercept = 1, linetype = 2, linewidth = 1) + xlab("Sample")

p4 <- ggplot(combined_data_2, aes(fill = treatment, x = sample_group, y = C7_rel_conc)) +
  geom_bar(position = "dodge", stat = "identity", color = "black") + theme_prism(axis_text_angle = 0, border= T) + scale_fill_prism(palette = "floral") +
  scale_y_continuous(guide = "prism_minor", limits = c(0, 2), expand = c(0, 0)) + ylab("CEACAM7 mRNA FC") +
  geom_hline(yintercept = 1, linetype = 2) + xlab("Sample") + guides(colour = guide_legend(position = "inside")) +
  theme(legend.position = c(0.8, 0.8)) +
  coord_cartesian(clip = "off")
p4 #looks good             



IVL_data_2 <- summarised_host_data_2 %>%
  filter(`Target Name` == "IVL") %>%
  arrange(as.numeric(`Sample Name`)) %>%
  rename("IVL_Ct" = "mean_Ct")
IVL_data_2

combined_data_2 <- left_join(combined_data_2, IVL_data_2, by = "Sample Name")
#View(combined_data_2)

combined_data_2 <- mutate(combined_data_2, IVL_delta_Ct = TBP_Ct - IVL_Ct)

ggplot(combined_data_2, aes(x = reorder(`Sample Name`, as.numeric(`Sample Name`)), y = IVL_delta_Ct)) +
  geom_point() + theme_prism(axis_text_angle = 0) + 
  scale_y_continuous(guide = "prism_minor", limits = c(-1, 6), expand = c(0, 0)) + xlab("Sample")

IVL_mean_control <- filter(combined_data_2, `Sample Name` == "8.2") %>% pull(IVL_delta_Ct)
IVL_mean_control

combined_data_2 <- combined_data_2 %>% 
  mutate(IVL_delta_delta_Ct = IVL_mean_control - IVL_delta_Ct)

ggplot(combined_data_2, aes(x = reorder(`Sample Name`, as.numeric(`Sample Name`)), y = IVL_delta_delta_Ct)) +
  geom_point() + theme_prism(axis_text_angle = 0) + 
  scale_y_continuous(guide = "prism_minor") + xlab("Sample")

combined_data_2 <- combined_data_2 %>%
  mutate(IVL_rel_conc = 2^-IVL_delta_delta_Ct)

ggplot(combined_data_2, aes(x = reorder(`Sample Name`, as.numeric(`Sample Name`)), y = IVL_rel_conc)) +
  geom_point() + theme_prism(axis_text_angle = 0) + 
  scale_y_continuous(guide = "prism_minor", limits = c(0, 2), expand = c(0, 0)) + ylab("IVL mRNA FC") +
  geom_hline(yintercept = 1, linetype = 2, linewidth = 1) + xlab("Sample")

p5 <- ggplot(combined_data_2, aes(fill = treatment, x = sample_group, y = IVL_rel_conc)) +
  geom_bar(position = "dodge", stat = "identity", color = "black") + theme_prism(axis_text_angle = 0, border= T) + scale_fill_prism(palette = "floral") +
  scale_y_continuous(guide = "prism_minor", limits = c(0, 2), expand = c(0, 0)) + ylab("IVL mRNA FC") +
  geom_hline(yintercept = 1, linetype = 2) + xlab("Sample") + guides(colour = guide_legend(position = "inside")) +
  theme(legend.position = c(0.9, 0.8)) +
  coord_cartesian(clip = "off")
p5 #looks good                       


plot_all <- plot_grid(p1+theme_prism(base_size = 8)+theme(legend.position="none"),
                      p2+theme_prism(base_size = 8)+theme(legend.position="none"),
                      p3+theme_prism(base_size = 8)+theme(legend.position="none"),
                      p4+theme_prism(base_size = 8)+theme(legend.position="none"),
                      p5+theme_prism(base_size = 8)+theme(legend.position=c(1.5, 0.5)),
                      labels = c("A", "B", "C", "D", "E"),
                      label_size = 12)
plot_all
save_plot("CA_HTK_Raft_2_host_only.pdf", plot_all)







#add in the data from the first run, side by side


#focus on PPARA, pull all the runs, rename the mean_Ct to PPARA_Ct
PPARA_data_all <- bind_rows(PPARA_data, PPARA_data_2)
PPARA_data_all


#do the same for the reference gene TBP
ref_data_all <- bind_rows(ref_data, ref_data_2)
ref_data_all

#merge the PPARA and TBP tables, calc dCT
combined_data_all <- left_join(PPARA_data_all, ref_data_all, by = "Sample Name")
combined_data_all
combined_data_all <- mutate(combined_data_all, PPARA_delta_Ct = TBP_Ct - PPARA_Ct)
combined_data_all


ggplot(combined_data_all, aes(x = reorder(`Sample Name`, as.numeric(`Sample Name`)), y = PPARA_delta_Ct)) +
  geom_point() + theme_prism(axis_text_angle = 0) + 
  scale_y_continuous(guide = "prism_minor", limits = c(-4, 1), expand = c(0, 0)) + xlab("Sample")

#choose the sample to normalize on, 2.1 in this case
PPARA_mean_control <- filter(combined_data_all, `Sample Name` == "2_T2Y-DMSO") %>% pull(PPARA_delta_Ct)
PPARA_mean_control

#2^deltadeltaCT to get the fold change
combined_data_all <- combined_data_all %>% 
  mutate(PPARA_delta_delta_Ct = PPARA_mean_control - PPARA_delta_Ct)

ggplot(combined_data_all, aes(x = reorder(`Sample Name`, as.numeric(`Sample Name`)), y = PPARA_delta_delta_Ct)) +
  geom_point() + theme_prism(axis_text_angle = 0) + 
  scale_y_continuous(guide = "prism_minor") + xlab("Sample")

combined_data_all <- combined_data_all %>%
  mutate(PPARA_rel_conc = 2^-PPARA_delta_delta_Ct)

ggplot(combined_data_all, aes(x = reorder(`Sample Name`, as.numeric(`Sample Name`)), y = PPARA_rel_conc)) +
  geom_point() + theme_prism(axis_text_angle = 0) + 
  scale_y_continuous(guide = "prism_minor", limits = c(0, 10), expand = c(0, 0)) + ylab("PPARA mRNA FC") +
  geom_hline(yintercept = 1, linetype = 2, linewidth = 1) + xlab("Sample")


#View(combined_data_all)
combined_data_all$sample_group = c(rep("T316", 3), rep("T2Y", 3), rep("T3Y", 3),  rep("T216", 3), rep("T216_2", 3), rep("T316_2", 3))
combined_data_all$treatment = c("GW", "DMSO", "Feno", "GW", "DMSO", "Feno", "GW", "DMSO", "Feno", "GW", "DMSO", "Feno", "GW", "DMSO", "Feno", "GW", "DMSO", "Feno")


p3 <- ggplot(combined_data_all, aes(fill = treatment, x = sample_group, y = PPARA_rel_conc)) +
  geom_bar(position = "dodge", stat = "identity", color = "black") + theme_prism(axis_text_angle = 0, border= T) + scale_fill_prism(palette = "floral") +
  scale_y_continuous(guide = "prism_minor", limits = c(0, 10), expand = c(0, 0)) + ylab("PPARA mRNA FC") +
  geom_hline(yintercept = 1, linetype = 2) + xlab("Sample") + guides(colour = guide_legend(position = "inside")) +
  theme(legend.position = c(0.8, 0.8)) +
  coord_cartesian(clip = "off")
p3 





#Ceacam7
C7_data_all <- bind_rows(C7_data, C7_data_2)
C7_data_all

#merge the PPARA and TBP tables, calc dCT
combined_data_all <- left_join(combined_data_all, C7_data_all, by = "Sample Name")
#View(combined_data_all)
combined_data_all <- mutate(combined_data_all, C7_delta_Ct = TBP_Ct - C7_Ct)
#View(combined_data_all)


ggplot(combined_data_all, aes(x = reorder(`Sample Name`, as.numeric(`Sample Name`)), y = C7_delta_Ct)) +
  geom_point() + theme_prism(axis_text_angle = 0) + 
  scale_y_continuous(guide = "prism_minor", limits = c(-4, 10), expand = c(0, 0)) + xlab("Sample")

#choose the sample to normalize on, 2.1 in this case
C7_mean_control <- filter(combined_data_all, `Sample Name` == "2_T2Y-DMSO") %>% pull(C7_delta_Ct)
C7_mean_control

#2^deltadeltaCT to get the fold change
combined_data_all <- combined_data_all %>% 
  mutate(C7_delta_delta_Ct = C7_mean_control - C7_delta_Ct)

ggplot(combined_data_all, aes(x = reorder(`Sample Name`, as.numeric(`Sample Name`)), y = C7_delta_delta_Ct)) +
  geom_point() + theme_prism(axis_text_angle = 0) + 
  scale_y_continuous(guide = "prism_minor") + xlab("Sample")

combined_data_all <- combined_data_all %>%
  mutate(C7_rel_conc = 2^-C7_delta_delta_Ct)

ggplot(combined_data_all, aes(x = reorder(`Sample Name`, as.numeric(`Sample Name`)), y = C7_rel_conc)) +
  geom_point() + theme_prism(axis_text_angle = 0) + 
  scale_y_continuous(guide = "prism_minor", limits = c(0, 50), expand = c(0, 0)) + ylab("C7 mRNA FC") +
  geom_hline(yintercept = 1, linetype = 2, linewidth = 1) + xlab("Sample")

p4 <- ggplot(combined_data_all, aes(fill = treatment, x = sample_group, y = C7_rel_conc)) +
  geom_bar(position = "dodge", stat = "identity", color = "black") + theme_prism(axis_text_angle = 0, border= T) + scale_fill_prism(palette = "floral") +
  scale_y_continuous(guide = "prism_minor", limits = c(0, 50), expand = c(0, 0)) + ylab("C7 mRNA FC") +
  geom_hline(yintercept = 1, linetype = 2) + xlab("Sample") + guides(colour = guide_legend(position = "inside")) +
  theme(legend.position = c(0.8, 0.8)) +
  coord_cartesian(clip = "off")
p4 


#IVL
IVL_data_all <- bind_rows(IVL_data, IVL_data_2)
IVL_data_all

#merge the IVL and TBP tables, calc dCT
combined_data_all <- left_join(combined_data_all, IVL_data_all, by = "Sample Name")
#View(combined_data_all)
combined_data_all <- mutate(combined_data_all, IVL_delta_Ct = TBP_Ct - IVL_Ct)
#View(combined_data_all)


ggplot(combined_data_all, aes(x = reorder(`Sample Name`, as.numeric(`Sample Name`)), y = IVL_delta_Ct)) +
  geom_point() + theme_prism(axis_text_angle = 0) + 
  scale_y_continuous(guide = "prism_minor", limits = c(-4, 10), expand = c(0, 0)) + xlab("Sample")

#choose the sample to normalize on, 2.1 in this case
IVL_mean_control <- filter(combined_data_all, `Sample Name` == "2_T2Y-DMSO") %>% pull(IVL_delta_Ct)
IVL_mean_control

#2^deltadeltaCT to get the fold change
combined_data_all <- combined_data_all %>% 
  mutate(IVL_delta_delta_Ct = IVL_mean_control - IVL_delta_Ct)

ggplot(combined_data_all, aes(x = reorder(`Sample Name`, as.numeric(`Sample Name`)), y = IVL_delta_delta_Ct)) +
  geom_point() + theme_prism(axis_text_angle = 0) + 
  scale_y_continuous(guide = "prism_minor") + xlab("Sample")

combined_data_all <- combined_data_all %>%
  mutate(IVL_rel_conc = 2^-IVL_delta_delta_Ct)

ggplot(combined_data_all, aes(x = reorder(`Sample Name`, as.numeric(`Sample Name`)), y = IVL_rel_conc)) +
  geom_point() + theme_prism(axis_text_angle = 0) + 
  scale_y_continuous(guide = "prism_minor", limits = c(0, 25), expand = c(0, 0)) + ylab("IVL mRNA FC") +
  geom_hline(yintercept = 1, linetype = 2, linewidth = 1) + xlab("Sample")

p5 <- ggplot(combined_data_all, aes(fill = treatment, x = sample_group, y = IVL_rel_conc)) +
  geom_bar(position = "dodge", stat = "identity", color = "black") + theme_prism(axis_text_angle = 0, border= T) + scale_fill_prism(palette = "floral") +
  scale_y_continuous(guide = "prism_minor", limits = c(0, 25), expand = c(0, 0)) + ylab("IVL mRNA FC") +
  geom_hline(yintercept = 1, linetype = 2) + xlab("Sample") + guides(colour = guide_legend(position = "inside")) +
  theme(legend.position = c(0.8, 0.8)) +
  coord_cartesian(clip = "off")
p5 


#cowplot
plot_all <- plot_grid(p1+theme_prism(base_size = 8)+theme(legend.position="none"),
                      p2+theme_prism(base_size = 8)+theme(legend.position="none"),
                      p3+theme_prism(base_size = 8)+theme(legend.position="none"),
                      p4+theme_prism(base_size = 8)+theme(legend.position="none"),
                      p5+theme_prism(base_size = 8)+theme(legend.position=c(1.5, 0.5)),
                      labels = c("A", "B", "C", "D", "E"),
                      label_size = 12)
plot_all
save_plot("CA_HTK_Raft_1and2_host_only.pdf", plot_all, base_width = 9, base_height = 6)



#average the values from the first and second run for T216 and T316







#viral, for raft exp 2
viral_data_2 <- filter(tidy_data_2, `Target Name` != "IVL" & `Target Name` != "CEACAM7" & `Target Name` != "PPARa" & `Target Name` != "ADD3" & `Target Name` != "LYPD2" & 
                       `Sample Name` != "NTC" &
                       `Sample Name` != "1.1" &
                       `Sample Name` != "2.1" &
                       `Sample Name` != "3.1" &
                       `Sample Name` != "4.1" &
                       `Sample Name` != "5.1" &
                       `Sample Name` != "6.1" &
                       `Sample Name` != "7.1" &
                       `Sample Name` != "8.1" &
                       `Sample Name` != "9.1" &
                       `Sample Name` != "10.1" &
                       `Sample Name` != "11.1" &
                       `Sample Name` != "12.1")
viral_data_2

summarised_viral_data_2 <- viral_data_2 %>%
  mutate(CT = as.numeric(CT)) %>%
  group_by(`Sample Name`, `Target Name`) %>%
  summarise(mean_Ct = mean(CT))
summarised_viral_data_2

ggplot(summarised_viral_data_2, aes(x = `Sample Name`, y = mean_Ct, colour = `Target Name`)) +
  geom_point() + theme_prism(axis_text_angle = 45) + 
  scale_y_continuous(guide = "prism_minor", limits = c(0, 40), expand = c(0, 0))



## E6

E6_data_2 <- summarised_viral_data_2 %>%
  filter(`Target Name` == "E6*") %>%
  arrange(as.numeric(`Sample Name`)) %>%
  rename("E6_Ct" = "mean_Ct")
E6_data_2

ref_data_2 <- summarised_viral_data_2 %>%
  filter(`Target Name` == "TBP") %>%
  rename("TBP_Ct" = "mean_Ct") %>%
  arrange(as.numeric(`Sample Name`))
ref_data_2

combined_data <- left_join(E6_data_2, ref_data_2, by = "Sample Name")
combined_data
combined_data <- mutate(combined_data, E6_delta_Ct = TBP_Ct - E6_Ct)
combined_data

ggplot(combined_data, aes(x = reorder(`Sample Name`, as.numeric(`Sample Name`)), y = E6_delta_Ct)) +
  geom_point() + theme_prism(axis_text_angle = 0) + 
  scale_y_continuous(guide = "prism_minor", limits = c(0, 7), expand = c(0, 0)) + xlab("Sample")

E6_mean_control <- filter(combined_data, `Sample Name` == "8.2") %>% pull(E6_delta_Ct)
E6_mean_control

combined_data <- combined_data %>% 
  mutate(E6_delta_delta_Ct = E6_mean_control - E6_delta_Ct)

ggplot(combined_data, aes(x = reorder(`Sample Name`, as.numeric(`Sample Name`)), y = E6_delta_delta_Ct)) +
  geom_point() + theme_prism(axis_text_angle = 0) + 
  scale_y_continuous(guide = "prism_minor") + xlab("Sample")

combined_data <- combined_data %>%
  mutate(E6_rel_conc = 2^-E6_delta_delta_Ct)

ggplot(combined_data, aes(x = reorder(`Sample Name`, as.numeric(`Sample Name`)), y = E6_rel_conc)) +
  geom_point() + theme_prism(axis_text_angle = 0) + 
  scale_y_continuous(guide = "prism_minor", limits = c(0, 2), expand = c(0, 0)) + ylab("E6* mRNA FC") +
  geom_hline(yintercept = 1, linetype = 2, linewidth = 1) + xlab("Sample")

combined_data$sample_group = c(rep("T216", 3), rep("T316", 3))
combined_data$treatment =c("GW", "DMSO", "Feno", "GW", "DMSO", "Feno")

p1 <- ggplot(combined_data, aes(fill = treatment, x = sample_group, y = E6_rel_conc)) +
  geom_bar(position = position_dodge(0.6), stat = "identity", color = "black", width = 0.5) + theme_prism(axis_text_angle = 0, border= T) + scale_fill_prism(palette = "floral") +
  scale_y_continuous(guide = "prism_minor", limits = c(0, 2), expand = c(0, 0)) + ylab("E6* mRNA FC") +
  geom_hline(yintercept = 1, linetype = 2) + xlab("Sample") + guides(colour = guide_legend(position = "inside")) +
  theme(legend.position = c(0.8, 0.8)) +
  coord_cartesian(clip = "off")
p1 



## E1^E4

E1E4_data_2 <- summarised_viral_data_2 %>%
  filter(`Target Name` == "E1^E4") %>%
  arrange(as.numeric(`Sample Name`)) %>%
  rename("E1E4_Ct" = "mean_Ct")
E1E4_data_2

combined_data <- left_join(E1E4_data_2, ref_data_2, by = "Sample Name")
combined_data
combined_data <- mutate(combined_data, E1E4_delta_Ct = TBP_Ct - E1E4_Ct)
combined_data

ggplot(combined_data, aes(x = reorder(`Sample Name`, as.numeric(`Sample Name`)), y = E1E4_delta_Ct)) +
  geom_point() + theme_prism(axis_text_angle = 0) + 
  scale_y_continuous(guide = "prism_minor", limits = c(0, 7), expand = c(0, 0)) + xlab("Sample")

E1E4_mean_control <- filter(combined_data, `Sample Name` == "8.2") %>% pull(E1E4_delta_Ct)
E1E4_mean_control

combined_data <- combined_data %>% 
  mutate(E1E4_delta_delta_Ct = E1E4_mean_control - E1E4_delta_Ct)

ggplot(combined_data, aes(x = reorder(`Sample Name`, as.numeric(`Sample Name`)), y = E1E4_delta_delta_Ct)) +
  geom_point() + theme_prism(axis_text_angle = 0) + 
  scale_y_continuous(guide = "prism_minor") + xlab("Sample")

combined_data <- combined_data %>%
  mutate(E1E4_rel_conc = 2^-E1E4_delta_delta_Ct)

ggplot(combined_data, aes(x = reorder(`Sample Name`, as.numeric(`Sample Name`)), y = E1E4_rel_conc)) +
  geom_point() + theme_prism(axis_text_angle = 0) + 
  scale_y_continuous(guide = "prism_minor", limits = c(0, 2), expand = c(0, 0)) + ylab("E1^E4 mRNA FC") +
  geom_hline(yintercept = 1, linetype = 2, linewidth = 1) + xlab("Sample")

combined_data$sample_group = c(rep("T216", 3), rep("T316", 3))
combined_data$treatment =c("GW", "DMSO", "Feno", "GW", "DMSO", "Feno")

p2 <- ggplot(combined_data, aes(fill = treatment, x = sample_group, y = E1E4_rel_conc)) +
  geom_bar(position = position_dodge(0.6), stat = "identity", color = "black", width = 0.5) + theme_prism(axis_text_angle = 0, border= T) + scale_fill_prism(palette = "floral") +
  scale_y_continuous(guide = "prism_minor", limits = c(0, 2), expand = c(0, 0)) + ylab("E1^E4 mRNA FC") +
  geom_hline(yintercept = 1, linetype = 2) + xlab("Sample") + guides(colour = guide_legend(position = "inside")) +
  theme(legend.position = c(0.2, 0.8)) +
  coord_cartesian(clip = "off")
p2 ### looks good



#combined viral reads
#focus on E6, pull all the runs, rename the mean_Ct to E6_Ct
E6_data_all <- bind_rows(E6_data, E6_data_2)
E6_data_all


#do the same for the reference gene TBP
ref_data_all <- bind_rows(ref_data, ref_data_2)
ref_data_all

#merge the E6 and TBP tables, calc dCT
combined_data_all <- left_join(E6_data_all, ref_data_all, by = "Sample Name")
combined_data_all
combined_data_all <- mutate(combined_data_all, E6_delta_Ct = TBP_Ct - E6_Ct)
combined_data_all


ggplot(combined_data_all, aes(x = reorder(`Sample Name`, as.numeric(`Sample Name`)), y = E6_delta_Ct)) +
  geom_point() + theme_prism(axis_text_angle = 0) + 
  scale_y_continuous(guide = "prism_minor", limits = c(-4, 10), expand = c(0, 0)) + xlab("Sample")

#choose the sample to normalize on, 8_T216-DMSO in this case
E6_mean_control <- filter(combined_data_all, `Sample Name` == "8_T216-DMSO") %>% pull(E6_delta_Ct)
E6_mean_control

#delta delta
combined_data_all <- combined_data_all %>% 
  mutate(E6_delta_delta_Ct = E6_mean_control - E6_delta_Ct)

ggplot(combined_data_all, aes(x = reorder(`Sample Name`, as.numeric(`Sample Name`)), y = E6_delta_delta_Ct)) +
  geom_point() + theme_prism(axis_text_angle = 0) + 
  scale_y_continuous(guide = "prism_minor") + xlab("Sample")

#2^deltadeltaCT to get the fold change
combined_data_all <- combined_data_all %>%
  mutate(E6_rel_conc = 2^-E6_delta_delta_Ct)

ggplot(combined_data_all, aes(x = reorder(`Sample Name`, as.numeric(`Sample Name`)), y = E6_rel_conc)) +
  geom_point() + theme_prism(axis_text_angle = 0) + 
  scale_y_continuous(guide = "prism_minor", limits = c(0, 2), expand = c(0, 0)) + ylab("E6 mRNA FC") +
  geom_hline(yintercept = 1, linetype = 2, linewidth = 1) + xlab("Sample")

View(combined_data_all)
combined_data_all$cells = c(rep("T316", 3), rep("T216", 6), rep("T316", 3))
combined_data_all$treatment = c("GW", "DMSO", "Feno", "GW", "DMSO", "Feno", "GW", "DMSO", "Feno", "GW", "DMSO", "Feno")
combined_data_all$cells_treatment= c("T316_GW", "T316_DMSO", "T316_Feno", "T216_GW", "T216_DMSO", "T216_Feno", "T216_GW", "T216_DMSO", "T216_Feno", "T316_GW", "T316_DMSO", "T316_Feno")


#sample_alias
combined_data_all_practice <- combined_data_all
combined_data_all_practice$sample_alias = c("T316_1", "T316_1", "T316_1", "T216_1", "T216_1", "T216_1", "T216_2", "T216_2", "T216_2", "T316_2", "T316_2", "T316_2")
view(combined_data_all_practice)

#plot rep 1 and 2 side by side
p5 <- ggplot(combined_data_all_practice, aes(fill = treatment, x = sample_alias, y = E6_rel_conc)) +
  geom_bar(position = position_dodge(0.6), stat = "identity", color = "black", width = 0.5) + theme_prism(axis_text_angle = 0, border= T) + scale_fill_prism(palette = "floral") +
  scale_y_continuous(guide = "prism_minor", limits = c(0, 2), expand = c(0, 0)) + ylab("E6* mRNA FC") +
  geom_hline(yintercept = 1, linetype = 2) + xlab("Sample") + guides(colour = guide_legend(position = "inside")) +
  theme(legend.position = c(0.2, 0.8)) +
  coord_cartesian(clip = "off")
p5



# now for stats, using dplyr
# mean
avgs <- combined_data_all %>%
  # Specify group indicator, column, function
  group_by(cells_treatment) %>%
  # Calculate the mean of the "Frequency" column for each group
  summarise_at(vars(E6_rel_conc),
               list(groupMean = mean))

# standard deviation
sds <- combined_data_all %>%
  # Specify group indicator, column, function
  group_by(cells_treatment) %>%
  # Calculate the mean of the "Frequency" column for each group
  summarise_at(vars(E6_rel_conc),
               list(groupSd = sd))

#combined stats
stats <- left_join(avgs, sds, by = "cells_treatment")

#combine stats with combined_data_all
combined_data_all <- left_join(combined_data_all, stats, by = "cells_treatment")


p3 <- ggplot(combined_data_all, aes(fill = treatment, x = cells_treatment, y = groupMean)) +
  geom_bar(position = "dodge", stat = "identity", color = "black") + geom_point(aes(y = E6_rel_conc)) + theme_prism(axis_text_angle = 0, border= T) + scale_fill_prism(palette = "floral") +
  scale_y_continuous(guide = "prism_minor", limits = c(0, 2), expand = c(0, 0)) + ylab("E6 mRNA FC") +
  geom_hline(yintercept = 1, linetype = 2) + xlab("Sample") + guides(colour = guide_legend(position = "inside")) +
  geom_errorbar(aes(ymin = groupMean-groupSd, ymax = groupMean+groupSd),
                position = position_dodge(0.9), width = .3) +
  theme(legend.position.inside = c(0.8, 0.8)) +
  coord_cartesian(clip = "off")
p3 







#focus on E1^E4, pull all the runs, rename the mean_Ct to E1E4_Ct
E1E4_data_all <- bind_rows(E1E4_data, E1E4_data_2)
E1E4_data_all


#do the same for the reference gene TBP
ref_data_all <- bind_rows(ref_data, ref_data_2)
ref_data_all

#merge the E1E4 and TBP tables, calc dCT
combined_data_all <- left_join(E1E4_data_all, ref_data_all, by = "Sample Name")
combined_data_all
combined_data_all <- mutate(combined_data_all, E1E4_delta_Ct = TBP_Ct - E1E4_Ct)
combined_data_all


ggplot(combined_data_all, aes(x = reorder(`Sample Name`, as.numeric(`Sample Name`)), y = E1E4_delta_Ct)) +
  geom_point() + theme_prism(axis_text_angle = 0) + 
  scale_y_continuous(guide = "prism_minor", limits = c(-4, 10), expand = c(0, 0)) + xlab("Sample")

#choose the sample to normalize on, 8_T216-DMSO in this case
E1E4_mean_control <- filter(combined_data_all, `Sample Name` == "8_T216-DMSO") %>% pull(E1E4_delta_Ct)
E1E4_mean_control

#delta delta
combined_data_all <- combined_data_all %>% 
  mutate(E1E4_delta_delta_Ct = E1E4_mean_control - E1E4_delta_Ct)

ggplot(combined_data_all, aes(x = reorder(`Sample Name`, as.numeric(`Sample Name`)), y = E1E4_delta_delta_Ct)) +
  geom_point() + theme_prism(axis_text_angle = 0) + 
  scale_y_continuous(guide = "prism_minor") + xlab("Sample")

#2^deltadeltaCT to get the fold change
combined_data_all <- combined_data_all %>%
  mutate(E1E4_rel_conc = 2^-E1E4_delta_delta_Ct)

ggplot(combined_data_all, aes(x = reorder(`Sample Name`, as.numeric(`Sample Name`)), y = E1E4_rel_conc)) +
  geom_point() + theme_prism(axis_text_angle = 0) + 
  scale_y_continuous(guide = "prism_minor", limits = c(0, 2), expand = c(0, 0)) + ylab("E1E4 mRNA FC") +
  geom_hline(yintercept = 1, linetype = 2, linewidth = 1) + xlab("Sample")

View(combined_data_all)
combined_data_all$cells = c(rep("T316", 3), rep("T216", 6), rep("T316", 3))
combined_data_all$treatment = c("GW", "DMSO", "Feno", "GW", "DMSO", "Feno", "GW", "DMSO", "Feno", "GW", "DMSO", "Feno")
combined_data_all$cells_treatment= c("T316_GW", "T316_DMSO", "T316_Feno", "T216_GW", "T216_DMSO", "T216_Feno", "T216_GW", "T216_DMSO", "T216_Feno", "T316_GW", "T316_DMSO", "T316_Feno")

#sample_alias
combined_data_all_practice <- combined_data_all
combined_data_all_practice$sample_alias = c("T316_1", "T316_1", "T316_1", "T216_1", "T216_1", "T216_1", "T216_2", "T216_2", "T216_2", "T316_2", "T316_2", "T316_2")
view(combined_data_all_practice)

#plot rep 1 and 2 side by side
p6 <- ggplot(combined_data_all_practice, aes(fill = treatment, x = sample_alias, y = E1E4_rel_conc)) +
  geom_bar(position = position_dodge(0.6), stat = "identity", color = "black", width = 0.5) + theme_prism(axis_text_angle = 0, border= T) + scale_fill_prism(palette = "floral") +
  scale_y_continuous(guide = "prism_minor", limits = c(0, 2), expand = c(0, 0)) + ylab("E1^E4 mRNA FC") +
  geom_hline(yintercept = 1, linetype = 2) + xlab("Sample") + guides(colour = guide_legend(position = "inside")) +
  theme(legend.position = c(0.2, 0.8)) +
  coord_cartesian(clip = "off")
p6

# now for stats, using dplyr
# mean
avgs <- combined_data_all %>%
  # Specify group indicator, column, function
  group_by(cells_treatment) %>%
  # Calculate the mean of the "Frequency" column for each group
  summarise_at(vars(E1E4_rel_conc),
               list(groupMean = mean))

# standard deviation
sds <- combined_data_all %>%
  # Specify group indicator, column, function
  group_by(cells_treatment) %>%
  # Calculate the mean of the "Frequency" column for each group
  summarise_at(vars(E1E4_rel_conc),
               list(groupSd = sd))

#combined stats
stats <- left_join(avgs, sds, by = "cells_treatment")

#combine stats with combined_data_all
combined_data_all <- left_join(combined_data_all, stats, by = "cells_treatment")


p4 <- ggplot(combined_data_all, aes(fill = treatment, x = cells_treatment, y = groupMean)) +
  geom_bar(position = "dodge", stat = "identity", color = "black") + geom_point(aes(y = E1E4_rel_conc)) + theme_prism(axis_text_angle = 0, border= T) + scale_fill_prism(palette = "floral") +
  scale_y_continuous(guide = "prism_minor", limits = c(0, 2.5), expand = c(0, 0)) + ylab("E1E4 mRNA FC") +
  geom_hline(yintercept = 1, linetype = 2) + xlab("Sample") + guides(colour = guide_legend(position = "inside")) +
  geom_errorbar(aes(ymin = groupMean-groupSd, ymax = groupMean+groupSd),
                position = position_dodge(0.9), width = .3) +
  theme(legend.position.inside = c(0.8, 0.8)) +
  coord_cartesian(clip = "off")
p4


#cowplot
plot_all <- plot_grid(p3+theme_prism(base_size = 8)+theme(legend.position="none"),
                      p4+theme_prism(base_size = 8)+theme(legend.position="none"),
                      p5+theme_prism(base_size = 8)+theme(legend.position="none"),
                      p6+theme_prism(base_size = 8)+theme(legend.position.inside=c(1.5, 0.5)),
                      labels = c("A", "B", "C", "D"),
                      label_size = 12)
plot_all
save_plot("CA_HTK_Raft_1and2_viral.pdf", plot_all, base_width = 9, base_height = 6)
