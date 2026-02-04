library(hrbrthemes)
library(viridis)
library(tidyverse)
library(ggpubr)
library(moments)
library(rstatix)
library(ggalluvial)
library(effsize)
library(patchwork)
library(biomaRt)
library(GOSemSim)
library(org.At.tair.db) 


# load function -----------------------------------------------------------


theme_Publication <- function(base_size=14) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#fdb462","grey","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#fdb462","grey","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}







# load data -----------------------------------------------------------
# Set working directory to script location (run from RStudio or source with chdir=TRUE)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

load('data/plasticity.rdata')
All_genes_corrected_CV2_table <- read_csv("data/All_genes_corrected_CV2_table.csv")

standardized_noise <- All_genes_corrected_CV2_table %>% 
  pivot_longer(-Gene, values_to = 'var', names_to = 'time') %>% 
  group_by(Gene) %>% 
  summarize(BioVar = median(var, na.rm = TRUE))



paralogs <- read_csv("data/paralogs.csv")
duplication_events <- read_csv("data/Young.csv")
duplication_events <- duplication_events %>% 
  mutate(ID = substr(ID,1,9))


gene_age <-  read_table("data/gene_age.csv") %>% filter(species == "Arabidopsis") %>% 
  select(-species)




df <- paralogs %>% select(RCC_scA,RCC_scB,`Duplication type1`,Ka,Ks) %>% 
  #dplyr::left_join(TF_hierarchy_bottom_up, by = c('RCC_scA' = 'name')) %>% ## network properties
  #dplyr::rename(A_hierarchy = level) %>% 
  #dplyr::left_join(TF_hierarchy_bottom_up, by = c('RCC_scB' = 'name')) %>% 
  #dplyr::rename(B_hierarchy = level) %>% 
  #mutate(A_hierarchy = ifelse(is.na(A_hierarchy), 0,A_hierarchy)) %>% 
  #mutate(B_hierarchy = ifelse(is.na(B_hierarchy), 0,B_hierarchy)) %>% 
  dplyr::left_join(standardized_noise, by = c('RCC_scA' = 'Gene')) %>% 
  dplyr::rename( A_BioVar = BioVar) %>% 
  dplyr::left_join(standardized_noise, by = c('RCC_scB' = 'Gene')) %>% 
  dplyr::rename( B_BioVar = BioVar) %>% 
  left_join(duplication_events %>%  select(ID, `Duplication Status`), by = c('RCC_scA' = 'ID')) %>% 
  dplyr::rename(A_status = `Duplication Status`) %>% 
  left_join(duplication_events %>%  select(ID, `Duplication Status`), by = c('RCC_scB' = 'ID')) %>% 
  dplyr::rename(B_status = `Duplication Status`) %>% 
  left_join(gene_age, by = c('RCC_scA' = 'GeneID')) %>% 
  dplyr::rename(A_Clade = Clade, A_exp_mean = ExpressionMean) %>% 
  left_join(gene_age, by = c('RCC_scB' = 'GeneID')) %>% 
  dplyr::rename(B_Clade = Clade, B_exp_mean = ExpressionMean) 


df_with_NEV <- df %>% filter(is.na(A_BioVar) == FALSE & is.na(B_BioVar) == FALSE)






# Fig 1 -------------------------------------------------------------------

p1_1 <- standardized_noise %>% 
  left_join(duplication_events, by = c('Gene' = 'ID')) %>% 
  dplyr::rename(status = `Duplication Status`) %>% 
  mutate(status = case_when(
    status %in% c('Early WGDs', 'Recent WGDs') ~ 'WGD',
    status == 'Non WGDs' ~ 'SSD',
    TRUE ~ status
  )) %>% 
  filter(!is.na(status)) %>% 
  ggplot(aes(x = BioVar, color = status, fill = status)) +
  geom_density(alpha = .1, linewidth = 0.5) +
  theme_Publication()+
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = 'top',
    legend.title = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text = element_text(size = 15),
    axis.title.y = element_text(size = 15, face = "bold", vjust = 1.5), 
    axis.title.x = element_text(size = 15, face = "bold")
  ) +
  labs(
    x = 'Normalized expression variability (NEV)',
    y = 'Density'
  ) +
  scale_fill_manual(values = c('grey', '#C08599', '#406D96')) +
  scale_color_manual(values = c('grey', '#C08599', '#406D96'))

group_A <- standardized_noise %>% 
  left_join(duplication_events, by = c('Gene' = 'ID')) %>% dplyr::rename(status = `Duplication Status`) %>% 
  mutate(status = ifelse(status == 'Early WGDs'| status == 'Recent WGDs', 'WGD',status)) %>% 
  mutate(status = ifelse(status == 'Non WGDs', 'SSD', status)) %>% 
  filter(is.na(status) == FALSE) %>% 
  filter(status == 'SSD')

group_B <- standardized_noise %>% 
  left_join(duplication_events, by = c('Gene' = 'ID')) %>% dplyr::rename(status = `Duplication Status`) %>% 
  mutate(status = ifelse(status == 'Early WGDs'| status == 'Recent WGDs', 'WGD',status)) %>% 
  mutate(status = ifelse(status == 'Non WGDs', 'SSD', status)) %>% 
  filter(is.na(status) == FALSE) %>% 
  filter(status == 'Singletons')

group_C <- standardized_noise %>% 
  left_join(duplication_events, by = c('Gene' = 'ID')) %>% dplyr::rename(status = `Duplication Status`) %>% 
  mutate(status = ifelse(status == 'Early WGDs'| status == 'Recent WGDs', 'WGD',status)) %>% 
  mutate(status = ifelse(status == 'Non WGDs', 'SSD', status)) %>% 
  filter(is.na(status) == FALSE) %>% 
  filter(status == 'WGD')



cliff.delta(group_A$BioVar, group_B$BioVar)
cliff.delta(group_C$BioVar, group_B$BioVar)
cliff.delta(group_A$BioVar, group_C$BioVar)






p1_2 <- df_with_NEV %>% 
  dplyr::mutate(`Duplicate (max)` = ifelse(A_BioVar > B_BioVar, A_BioVar, B_BioVar),
                `Duplicate (min)` = ifelse(A_BioVar > B_BioVar, B_BioVar, A_BioVar)) %>% 
  select(`Duplication type1`,`Duplicate (max)`, `Duplicate (min)`) %>% 
  pivot_longer(-`Duplication type1`, values_to = 'expression_variability', names_to = 'group') %>% 
  bind_rows(standardized_noise %>% 
              left_join(duplication_events, by = c('Gene' = 'ID')) %>% dplyr::rename(`Duplication type1` = `Duplication Status`, expression_variability = BioVar) %>% 
              filter(`Duplication type1` == 'Singletons') %>% 
              mutate(group = `Duplication type1`) %>% 
              select(`Duplication type1`,group, expression_variability)) %>% 
  mutate(duplication = ifelse(`Duplication type1` %in% c('alpha','beta','gamma'), 'WGD','SSD')) %>% 
  mutate(duplication = ifelse(`Duplication type1` == 'Singletons', 'Singletons', duplication)) %>%
  arrange(group) %>% 
  mutate(group=fct_reorder(group,expression_variability)) %>% 
  ggplot(aes(x = expression_variability, color = group ))+
  geom_density(aes(fill = group), alpha = .1, size = 0.5)+
  theme(legend.title = element_blank())+
  theme_Publication()+
  scale_colour_Publication() +
  theme(
    legend.position="top",
    legend.title = element_blank(),
    plot.title = element_text(size=40),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15)) +
  xlab('Normalized expression variability (NEV)')+
  ylab('Density')+
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )








summary <- standardized_noise %>% 
  left_join(duplication_events, by = c('Gene' = 'ID')) %>% dplyr::rename(status = `Duplication Status`) %>% 
  left_join(responsiveness %>% mutate(GeneID = str_sub(GeneID, start = 1, end = 9)), by = c('Gene' = 'GeneID')) #%>% 




p1_3 <- summary %>%  
  ggplot(aes(BioVar, responsiveness))+
  geom_point(fill = 'darkolivegreen', alpha = .5, size = 2, color ='white', shape = 21)+
  theme_Publication()+
  ylab('Environmental responsiveness')+
  xlab('Normalized expression variability (NEV)')+
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15))+
  #geom_hex(bins = 120) +
  geom_smooth(method = 'lm', color = 'red')+
  theme(legend.position = 'top',
        legend.text = element_blank())+
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.001, size = 5, method = "pearson", label.x = -3, label.y = 150000)+
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )




p1_1 + p1_2 + p1_3 + plot_annotation(tag_levels = 'A')

ggsave(filename = "fig_1.pdf")






# fig_2 -------------------------------------------------------------------



count_data <- df_with_NEV %>% 
  dplyr::mutate(high = ifelse(A_BioVar > B_BioVar, A_BioVar, B_BioVar),
                low = ifelse(A_BioVar > B_BioVar, B_BioVar, A_BioVar)) %>%
  mutate(A = ifelse(high == A_BioVar, RCC_scA, RCC_scB),
         B = ifelse(low == A_BioVar, RCC_scA, RCC_scB)) %>% 
  mutate(duplication = ifelse(`Duplication type1` %in% c('alpha','beta','gamma'), 'WGD','SSD')) %>% 
  mutate(group = row_number()) %>% 
  dplyr::select(RCC_scA, RCC_scB, high, low,duplication, group, A, B,`Duplication type1`) %>%
  mutate(diff = high - low) %>% 
  arrange(diff) %>% 
  slice_max(diff, prop = 0.05) %>%
  group_by(duplication) %>%
  summarise(n = n()) %>%
  ungroup()


p2_1 <- df_with_NEV %>% 
  dplyr::mutate(high = ifelse(A_BioVar > B_BioVar, A_BioVar, B_BioVar),
                low = ifelse(A_BioVar > B_BioVar, B_BioVar, A_BioVar)) %>%
  mutate(A = ifelse(high == A_BioVar, RCC_scA, RCC_scB),
         B = ifelse(low == A_BioVar, RCC_scA, RCC_scB)) %>% 
  mutate(duplication = ifelse(`Duplication type1` %in% c('alpha','beta','gamma'), 'WGD','SSD')) %>% 
  mutate(group = row_number()) %>% 
  dplyr::select(RCC_scA, RCC_scB, high, low,duplication, group, A, B,`Duplication type1`) %>%
  mutate(diff = high - low) %>% 
  arrange(diff) %>% 
  slice_max(diff, prop = 0.05) %>%
  pivot_longer(c(high,low),names_to = 'status',values_to = 'Expression variability') %>% 
  ggplot(aes(x = status, y = `Expression variability`, group = group))+geom_line(alpha = 0.9, size = 0.5, color = 'grey')+
  geom_point()+
  #geom_boxplot()+
  facet_wrap(~duplication)+
  geom_text(
    data = count_data, 
    aes(x = "high", y = Inf, label = paste0("n = ", n)), 
    inherit.aes = FALSE,  # Crucial: prevents inheritance of group aesthetic
    hjust = -0.2, vjust = 1.2, size = 5, color = "black"
  ) +
  scale_x_discrete("") +
  theme_Publication() +
  ylab('Normalized expression variability (NEV)')+
  theme(legend.position = "none", 
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 12),
        panel.grid = element_blank(),
        #axis.text.x = element_blank(),
        strip.text = element_text(size=20),
        axis.line.y = element_line(size = .5),
        legend.title = element_blank(),
        strip.background = element_blank())




#GO pairwise semantic comparison 


GO_term <- read_delim("data/Ath_gene2go",
                      delim = "\t", escape_double = FALSE,
                      col_names = FALSE, trim_ws = TRUE)




# BiocManager::install("GOSemSim")
# BiocManager::install("org.At.tair.db")




hsGO <- godata('org.At.tair.db', ont = "BP",  keytype = "TAIR")


data <- df_with_NEV %>% 
  dplyr::mutate(max = ifelse(A_BioVar > B_BioVar, A_BioVar, B_BioVar),
                min = ifelse(A_BioVar > B_BioVar, B_BioVar, A_BioVar)) %>%
  mutate(A = ifelse(max == A_BioVar, RCC_scA, RCC_scB),
         B = ifelse(min == A_BioVar, RCC_scA, RCC_scB)) %>% 
  mutate(duplication = ifelse(`Duplication type1` %in% c('alpha','beta','gamma'), 'WGD','SSD')) %>% 
  mutate(group = row_number()) %>% 
  dplyr::select(RCC_scA, RCC_scB, max, min,duplication, group, A, B,`Duplication type1`) %>%
  #left_join(Ath_gene2go, by = c('A' = 'X1')) %>% 
  #rename(A_GO = X2) %>% 
  #left_join(Ath_gene2go, by = c('B' = 'X1')) %>% 
  #rename(B_GO = X2) %>% 
  mutate(diff = max - min) %>% 
  arrange(diff) %>% 
  slice_max(diff, prop = 0.05) %>% 
  dplyr::select(RCC_scA,RCC_scB,duplication) %>% 
  left_join(GO_term, by = c("RCC_scA" = "X1")) %>% 
  dplyr::rename(A_go = `X2`) %>% 
  # mutate(
  #   A_go = case_when(
  #     is.na(A_go) | A_go == "" ~ NA_character_,
  #     TRUE ~ {
  #       result <- str_replace_all(A_go, c("^" = "c(\"",
  #                                             ", " = "\", \"",
  #                                             "$" = "\")"))
  #       str_replace_all(result, "\\\\", "") # Remove backslashes
  #     }
  #   )
  # ) %>% 
  left_join(GO_term, by = c("RCC_scB" = "X1")) %>% 
  dplyr::rename(B_go = X2) 




parse_go_terms <- function(go_string) {
  if (is.na(go_string) || go_string == "") {
    return(character(0)) # Return empty vector for NA/empty input
  }
  # Split by comma, allowing for spaces around the comma
  terms <- strsplit(go_string, ",\\s*")[[1]]
  # Trim leading/trailing whitespace from each term
  terms <- str_trim(terms)
  # Filter out any empty strings that might result (e.g., trailing comma)
  terms <- terms[terms != ""]
  return(terms)
}



data$pair_id <- rownames(data)
data_clean <- data %>%
  filter(!is.na(A_go) & A_go != "" & !is.na(B_go) & B_go != "")

list_A_go <- lapply(data_clean$A_go, parse_go_terms)
list_B_go <- lapply(data_clean$B_go, parse_go_terms)


calculate_similarity <- function(go_A, go_B, sem_data, measure, combine) {
  # Error handling: Return NA if either GO vector is empty
  if (length(go_A) == 0 || length(go_B) == 0) {
    return(NA)
  }
  
  # Calculate semantic similarity using mgoSim
  similarity <- tryCatch({
    mgoSim(go_A, go_B, semData = sem_data, measure = measure, combine = combine)
  }, error = function(e) {
    # If there's an error, return NA
    message(paste("Error calculating similarity:", e))
    return(NA)
  })
  return(similarity)
}

sim_scores_BP <- mapply(
  calculate_similarity,
  list_A_go,
  list_B_go,
  MoreArgs = list(sem_data = hsGO, measure = "Jiang", combine = "BMA")
)


data_clean$sim <- sim_scores_BP


p2_2 <- data_clean %>% 
  ggplot(aes( sim, color = duplication, fill = duplication))+
  geom_density(alpha = .4, size = 1)+
  theme_Publication()+
  theme(axis.text.y = element_blank(), legend.title = element_blank())+
  theme(legend.position = 'top',
        plot.title = element_text(size=25),
        axis.text = element_text(size = 15),
        #axis.text.y = element_blank(),
        axis.title = element_text(size = 15)
  ) +xlab('Semantic Similarity')+ylab('Density')+
  scale_fill_manual(values = c('#C08599','#406D96'))+
  scale_color_manual(values = c('#C08599','#406D96'))



p2_1 + p2_2 
ggsave(filename = 'fig_2.pdf')




# fig_3 -------------------------------------------------------------------



N = 3
symnum.args <- list(
  cutpoints = c(0, 0.0001/N, 0.001/N, 0.01/N, 0.05/N, 1),
  symbols = c("<0.0001", "<0.001", "<0.01", "<0.05", "ns")
)


library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

pdat <- df_with_NEV %>% 
  mutate(duplication = ifelse(`Duplication type1` %in% c("alpha","beta","gamma"), "WGD", "SSD")) %>% 
  pivot_longer(c(RCC_scA, RCC_scB), names_to = "type", values_to = "gene_name") %>% 
  dplyr::select(gene_name, duplication, type) %>% 
  left_join(standardized_noise, by = c("gene_name" = "Gene")) %>% 
  left_join(gene_age, by = c("gene_name" = "GeneID")) %>% 
  filter(!is.na(Clade)) %>% 
  mutate(Clade = as.numeric(Clade))

# breaks from the actual data
clade_breaks <- sort(unique(pdat$Clade))

p3_2 <- ggplot(pdat, aes(x = Clade, y = BioVar, color = duplication)) +
  geom_point(shape = 16, alpha = 0.2, size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(aes(color = duplication),
           size = 4,
           label.x = min(clade_breaks),
           show.legend = FALSE) +
  theme_Publication() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 18),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  ylab("NEV") +
  xlab("Clade") +
  scale_color_manual(values = c(SSD = "#C08599", WGD = "#406D96")) +
  scale_x_continuous(breaks = clade_breaks)




p3_3 <- df_with_NEV %>%
  mutate(
    duplication = ifelse(`Duplication type1` %in% c("alpha","beta","gamma"), "WGD", "SSD"),
    duplication = factor(duplication, levels = c("SSD", "WGD")),
    diff_biovar = A_BioVar - B_BioVar,
    age_diff    = A_Clade - B_Clade
  ) %>%
  mutate(a = diff_biovar * age_diff) %>%
  filter(a != 0) %>%
  ggplot(aes(age_diff, diff_biovar, color = duplication)) +
  geom_point(alpha = 0.4, size = 2) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1.2, show.legend = FALSE) +
  stat_cor(aes(color = duplication),
           label.x = -10, size = 4,
           show.legend = FALSE) +
  theme_Publication() +
  scale_color_manual(values = c(SSD = "#C08599", WGD = "#406D96")) +
  xlab(expression(Delta*Clade == Clade[A] - Clade[B])) +
  ylab(expression(Delta*NEV == NEV[A] - NEV[B]))+
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 4))) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.key = element_blank(),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 18),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )




(p3_2 | p3_3) + plot_annotation(tag_levels = 'A')

ggsave(filename = 'fig_3.pdf')
