# Load libraries
library(tidyverse)
library(ProteoBayes)
library(DAPAR)
library(cp4p)
library(mi4p)
library(yardstick)

# Set random generator
set.seed(17)

# Load peptide-level data
peptides <- read.delim("DATA/Arabido_UPS/peptides.txt")

# Data preprocessing
db_ARATH <- peptides %>% 
  # Select columns of interest
  select(Sequence, Leading.razor.protein, starts_with("Intensity.")) %>% 
  # Reshape data in long format
  gather(key = "Condition", value = "Intensity", 
         -c(Sequence, Leading.razor.protein)) %>% 
  # Transform the key into Group and Sample
  separate(col = Condition, into = c(NA, "Group","Sample")) %>% 
  # Rename columns to match ProteoBayes requirements
  rename(Peptide = Sequence, 
         Protein = Leading.razor.protein,
         Output = Intensity) %>%
  # Remove reverse and contaminant proteins, remove iRT. 
  filter(!str_detect(Protein, "CON") & !str_detect(Protein, "REV") &
           !str_detect(Protein, "iRT") & Output != 0) %>% 
  # Log2-transform intensity values
  # If Intensity = 0, the intensity value is missing, 
  # so it should be recoded as NA.
  mutate(Output = log2(Output))
  #mutate(Output = if_else(Output !=0, log2(Output), NA))

# Group = recode(Group,
#                "Point1" = "0.05fmol", "Point2" = "0.25fmol", 
#                "Point3" = "0.5fmol", "Point4" = "1.25fmol", 
#                "Point5" = "2.5fmol", "Point6" = "5fmol", 
#                "Point7" = "10fmol")



# ProteoBayes - Univariate setting
Start_time <- Sys.time()
res_uni_ARATH <- posterior_mean(data = db_ARATH,
                                mu_0 = db_ARATH %>% 
                                  group_by(Group) %>% 
                                  mutate(mu_0 = mean(Output)) %>% pull(mu_0))
End_time <- Sys.time()
duration_uni_ARATH <- End_time - Start_time

diff_uni_ARATH <- identify_diff(res_uni_ARATH)

# Moderated t-test - using DAPAR
# Note, limma and DAPAR use data in wide format.
db_ARATH_limma <- db_ARATH %>% 
  # Merge Group and Sample columns to denote biological samples analysed
  unite(col = "Cond_Rep", c(Group,Sample),sep = "_") %>% 
  # Reshape data in wide format
  spread(key = "Cond_Rep", value = "Output")

# Create quantitative data matrix to match DAPAR requirements
qARATH <- db_ARATH_limma %>%
  # Tibble to dataframe to add row names to be propagated in results tables
  as.data.frame(row.names = pull(db_ARATH_limma %>% select(Peptide))) %>% 
  select(-c(Peptide,Protein)) 

# Create design dataframe to match DAPAR requirements
metadata <- data.frame(Sample.name = colnames(qARATH),
                       Condition = rep(c("Point1", "Point2", "Point3",
                                         "Point4", "Point5", "Point6",
                                         "Point7"), rep(3,7)),
#                       Condition = rep(c("0.05fmol", "0.25fmol", "0.5fmol",
#                            "1.25fmol", "2.5fmol", "5fmol","10fmol"), rep(3,7)),
                       Bio.Rep = 1:ncol(qARATH))
                       
# Moderated t-test - OnevsOne setting
Start_time <- Sys.time()
#res_DAPAR <- limmaCompleteTest(qData = as.matrix(qARATH),
#                               sTab = metadata,
#                               comp.type = "OnevsOne")
# Issue with the formatting of the output of DAPAR function, 
# Therefore, using the one in mi4p.
res_limma <- mi4p::limmaCompleteTest.mod(qData = as.matrix(qARATH),
                              sTab = metadata,
                              comp.type = "OnevsOne")$res.l
End_time <- Sys.time()
duration_limma <- End_time - Start_time

res_limma$P_Value %>% rownames_to_column() %>%
  #as_tibble %>% 
  mutate(across(starts_with("Point"), 
         ~ cp4p::adjust.p(p = ., alpha = 0.01)$adjp %>% 
           select(adjusted.p) %>% pull, .names = "{col}_adj")) %>% head

res_limma_adj <- res_limma$P_Value %>% 
  rownames_to_column() %>% 
  as_tibble %>% 
  mutate(across(starts_with("Point"), 
                ~ cp4p::adjust.p(p = ., alpha = 0.01)$adjp %>% 
                  select(adjusted.p) %>% pull))
res_limma_adj %>% rownames_to_column()

res_limma_FC <- res_limma$logFC %>% 
  rownames_to_column(var = "Peptide") %>% 
  gather(key = "Comparison", 
         value = "log2FC", -Peptide) %>% 
  separate(col = Comparison, 
           into = c("Comparison", NA), 
           sep = "_logFC")

# -- Results -- #

db_results <- full_join(y = res_limma_adj %>% 
                        #y = res_limma$P_Value %>% 
                        #  rownames_to_column() %>% 
                          gather(key = "Comparison", 
                                 value = "pval", -rowname) %>% 
                          separate(col = Comparison, 
                                   into = c("Group1vsGroup2", NA), 
                                   sep = "_pval"),
                        x = diff_uni_ARATH %>% 
                          select(Peptide, Group, Group2, Distinct) %>%
                          unite(col = "Group1vsGroup2", c(Group,Group2),
                                sep = "_vs_"),
                        by = join_by(Peptide == rowname, Group1vsGroup2)) %>% 
  rename(Comparison = Group1vsGroup2) %>% 
  # Keep the 21 comparisons of interest
  filter(Comparison %in% c("Point1_vs_Point2", "Point1_vs_Point3", 
                               "Point1_vs_Point4", "Point1_vs_Point5",
                               "Point1_vs_Point6", "Point1_vs_Point7",
                               "Point2_vs_Point3", "Point2_vs_Point4", 
                               "Point2_vs_Point5", "Point2_vs_Point6", 
                               "Point2_vs_Point7",
                               "Point3_vs_Point4", "Point3_vs_Point5",
                               "Point3_vs_Point6", "Point3_vs_Point7",
                               "Point4_vs_Point5", "Point4_vs_Point6",
                               "Point4_vs_Point7",
                               "Point5_vs_Point6", "Point5_vs_Point7",
                               "Point6_vs_Point7")) %>% 
  # Add Protein Column
  left_join(y = .,
            x = db_ARATH %>%
               select(Peptide,Protein) %>%
               distinct,
             by = "Peptide") %>% 
  # Add column to flag differential peptides as in limma
  # and ground truth
  mutate(Limma = if_else(pval < 0.05, true = T, false = F, missing = NA),
         Truth = if_else(str_detect(Protein, "UPS"), true = T, false = F)
         ) %>%
  # Formatting
  rename(ProteoBayes = Distinct) %>% 
  select(-pval) %>% 
  filter(!(is.na(ProteoBayes) & is.na(Limma))) 

# True_diff = case_when(
#   # Non diff
#   !Truth ~ 0,
#   # Point 1
#   Truth & Comparison == "Point1_vs_Point2" ~ log2(0.25/0.05),
#   Truth & Comparison == "Point1_vs_Point3" ~ log2(0.5/0.05),
#   Truth & Comparison == "Point1_vs_Point4" ~ log2(1.25/0.05),
#   Truth & Comparison == "Point1_vs_Point5" ~ log2(2.5/0.05),
#   Truth & Comparison == "Point1_vs_Point6" ~ log2(5/0.05),
#   Truth & Comparison == "Point1_vs_Point7" ~ log2(10/0.05),
#   # Point 2
#   Truth & Comparison == "Point2_vs_Point3" ~ log2(0.5/0.25),
#   Truth & Comparison == "Point2_vs_Point4" ~ log2(1.25/0.25),
#   Truth & Comparison == "Point2_vs_Point5" ~ log2(2.5/0.25),
#   Truth & Comparison == "Point2_vs_Point6" ~ log2(5/0.25),
#   Truth & Comparison == "Point2_vs_Point7" ~ log2(10/0.25),
#   # Point 3
#   Truth & Comparison == "Point3_vs_Point4" ~ log2(1.25/0.5),
#   Truth & Comparison == "Point3_vs_Point5" ~ log2(2.5/0.5),
#   Truth & Comparison == "Point3_vs_Point6" ~ log2(5/0.5),
#   Truth & Comparison == "Point3_vs_Point7" ~ log2(10/0.5),
#   # Point 4
#   Truth & Comparison == "Point4_vs_Point5" ~ log2(2.5/1.25),
#   Truth & Comparison == "Point4_vs_Point6" ~ log2(5/1.25),
#   Truth & Comparison == "Point4_vs_Point7" ~ log2(10/1.25),
#   # Point 5
#   Truth & Comparison == "Point5_vs_Point6" ~ log2(5/2.5),
#   Truth & Comparison == "Point5_vs_Point7" ~ log2(10/2.5),
#   # Point 6
#   Truth & Comparison == "Point5_vs_Point7" ~ log2(10/5)
# )


db_results %>%
  group_by(Comparison, Truth, ProteoBayes, Limma) %>% 
  summarise(Count = n())

db_results %>% 
  group_by(ProteoBayes, Limma) %>% 
  summarise(Count = n())

# Perfomance indicators for ProteoBayes
perf_PB <- db_results %>% 
  as.data.frame() %>% 
  mutate(Truth = as.factor(Truth),
         ProteoBayes = as.factor(ProteoBayes)) %>% 
  group_by(Comparison) %>%
  group_modify(.f = ~ conf_mat(data = .x,
                               truth = Truth,
                               estimate = ProteoBayes) %>%
                 tidy %>%
                 mutate(name = recode(name,
                                      "cell_1_1" = "TN",
                                      "cell_2_1" = "FP",
                                      "cell_1_2" = "FN",
                                      "cell_2_2" = "TP"))) %>%
  spread(key = "name", value = "value") %>%
  mutate(Sensitivity = round(TP/(TP + FN)*100, digits = 1),
         Specificity = round(TN/(TN + FP)*100, digits = 1),
         Accuracy = round((TP + TN)/(TP + TN + FP + FN)*100, digits = 1),
         F1Score = round(2*TP/(2*TP+FP+FN)*100, digits = 1))

# Perfomance indicators for Limma
perf_LM <- db_results %>% 
  as.data.frame() %>% 
  mutate(Truth = as.factor(Truth),
         Limma = as.factor(Limma)) %>% 
  group_by(Comparison) %>%
  group_modify(.f = ~ conf_mat(data = .x,
                               truth = Truth,
                               estimate = Limma) %>%
                 tidy %>%
                 mutate(name = recode(name,
                                      "cell_1_1" = "TN",
                                      "cell_2_1" = "FP",
                                      "cell_1_2" = "FN",
                                      "cell_2_2" = "TP"))) %>%
  spread(key = "name", value = "value") %>%
  mutate(Sensitivity = round(TP/(TP + FN)*100, digits = 1),
         Specificity = round(TN/(TN + FP)*100, digits = 1),
         Accuracy = round((TP + TN)/(TP + TN + FP + FN)*100, digits = 1),
         F1Score = round(2*TP/(2*TP+FP+FN)*100, digits = 1))

db_perf_adj <- bind_rows(ProteoBayes = perf_PB, Limma = perf_LM, .id = "Method") %>% 
  arrange(Comparison) %>% relocate(Comparison)

db_ARATH %>% 
  mutate(Is_diff = if_else(str_detect(Protein, "UPS"), true = T, false = F)) %>% 
  group_by(Group, Is_diff) %>% 
  summarise(Count = n())

# Save results
# save(res_uni_ARATH, 
#      file = "DATA/Arabido_UPS/Arabido_UPS_PB_res")
# save(diff_uni_ARATH,
#      file = "DATA/Arabido_UPS/Arabido_UPS_PB_diff")
# save(res_limma,
#      file = "DATA/Arabido_UPS/Arabido_UPS_LM_res")
# save(db_perf,
#      file = "DATA/Arabido_UPS/Arabido_UPS_performance")

# Load results
load(file = "DATA/Arabido_UPS/Arabido_UPS_PB_res")
load(file = "DATA/Arabido_UPS/Arabido_UPS_PB_diff")
load(file = "DATA/Arabido_UPS/Arabido_UPS_LM_res")
load(file = "DATA/Arabido_UPS/Arabido_UPS_performance")


sample_post <- sample_distrib(res_uni_ARATH)

test_post <- sample_post %>% slice(1:100) %>%  
  group_by(Peptide) %>% 
  reframe(pivot_wider(., names_from = Group, values_from = Sample)) %>% rbind

# ----------------------------------------------------------- #
load(file = "Arabido_UPS_1of3inall_ImpMLE_long")
load(file = "Arabido_UPS_1of3inall_ImpMLE_res_dapar")
head(data.lg)

db_ARATH <- data.lg %>% 
  rename(Draw = Imp.Draw, Peptide = ID, Output = Intensity)
db_ARATH %>% head

# - Univariate approach - #

start_uni <- Sys.time()
res_uni_ARATH <- posterior_mean(db_ARATH)
end_uni <- Sys.time()
duration_uni <- end_uni - start_uni #15 sec env

start_sample_uni <- Sys.time()
sample_uni_ARATH <- sample_distrib(res_uni_ARATH)
end_sample_uni <- Sys.time()
duration_sample_uni <- end_sample_uni - start_sample_uni #38 sec env

# Add Protein to results
res_uni_ARATH <- right_join(x = res_uni_ARATH,
                            y = db_ARATH %>% 
                              select(Peptide,Protein) %>% 
                              distinct,
                            by = "Peptide")

sample_uni_ARATH <- right_join(x = sample_uni_ARATH,
                         y = db_ARATH %>% 
                           select(Peptide,Protein) %>% 
                           distinct,
                         by = "Peptide")

# - Multivariate approach - #
start_multi <- Sys.time()
res_multi_ARATH <- db_ARATH %>% 
  filter(Protein %in% (db_ARATH %>% 
                          select(Protein) %>% 
                          distinct %>% 
                          pull(Protein))) %>% 
  group_by(Protein) %>% 
  group_modify(.f = ~ multi_posterior_mean(.x, nu_0 = 30,
                                           vectorised = T))
end_multi <- Sys.time()
duration_multi <- end_multi - start_multi #56 min -> 9min

start_sample_multi <- Sys.time()
sample_multi_ARATH <- sample_distrib(res_multi_ARATH)
end_sample_multi <- Sys.time()
duration_sample_multi <- end_sample_multi - start_sample_multi

# Error in sample_distrib(res_multi_ARATH) : 
# The 'nu' parameter is too small compared to the number of Peptides 
# (nu < P - 1). Consider changing prior value for 'nu' or decreasing the number 
# of Peptides.

# nu = nu + 50 (mutate)

# Apply multivariate framework on each protein to benefit from intra-protein 
# correlation.
# However, when there is only 1 peptide in the protein group, there is no 
# sense in using a multivariate framework. 
# BUT we have imputed data --> How to manage the univariate framework ? 
# --> Uses more
