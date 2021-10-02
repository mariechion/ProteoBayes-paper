# Load functions and packages
source("bayes_proteo.R")
library(gridExtra)


# Load ARATH+UPS quantification data (MaxQuant export)
peptidesMQ <- read.delim("DATA/Arabido_UPS/peptides.txt")
# Data pre-filtering
peptidesMQ[,grep(pattern = "Intensity.",colnames(peptidesMQ))][peptidesMQ[,grep(pattern = "Intensity.",
                                                                                colnames(peptidesMQ))]==0]=NA
peptidesMQ_clean<-peptidesMQ[which(apply(is.na(peptidesMQ[,grep(pattern = "Intensity.Point1",colnames(peptidesMQ))]),1,sum)<3
                                   & apply(is.na(peptidesMQ[,grep(pattern = "Intensity.Point2",colnames(peptidesMQ))]),1,sum)<3
                                   & apply(is.na(peptidesMQ[,grep(pattern = "Intensity.Point3",colnames(peptidesMQ))]),1,sum)<3
                                   & apply(is.na(peptidesMQ[,grep(pattern = "Intensity.Point4",colnames(peptidesMQ))]),1,sum)<3
                                   & apply(is.na(peptidesMQ[,grep(pattern = "Intensity.Point5",colnames(peptidesMQ))]),1,sum)<3
                                   & apply(is.na(peptidesMQ[,grep(pattern = "Intensity.Point6",colnames(peptidesMQ))]),1,sum)<3
                                   & apply(is.na(peptidesMQ[,grep(pattern = "Intensity.Point7",colnames(peptidesMQ))]),1,sum)<3),]
peptidesMQ_clean<-subset(peptidesMQ_clean, subset = peptidesMQ_clean$Reverse!="+" & 
                           peptidesMQ_clean$Potential.contaminant!="+")

# ----------------------------------------------------------- #
# Load ARATH+UPS data imputed with MLE
load(file="DATA/Arabido_UPS/Arabido_UPS_1of3inall_ImpMLE")
dim(data.pept.imp)

# Array to list
data.imp.list <- lapply(seq(dim(data.pept.imp)[3]), function(aaa){
  data.pept.imp[,,aaa]
})
names(data.imp.list) <- as.character(seq(dim(data.pept.imp)[3]))
# List to long data frame with Draw ID
data.lg <- map(data.imp.list,~data.frame(.)) %>% 
  # Add ID variable using Peptide Sequence and 
  # Protein variable using Leading razor protein from the quantitative data
  map(., ~mutate(.,ID = peptidesMQ_clean$Sequence,
                 Protein = peptidesMQ_clean$Leading.razor.protein)) %>% 
  # Bind rows by Draw ID
  bind_rows(.id = 'Imp.Draw') %>% 
  # Rename sample columns using the quantitative data
  rename_with(.fn= ~colnames(peptidesMQ_clean %>% select(starts_with("Intensity."))),
              .cols = starts_with("X")) %>% 
  # Longer table with sample names
  pivot_longer(cols = starts_with("Intensity."), 
               names_to = "Sample",
               names_prefix = "Intensity.",
               values_to = "Intensity") %>% 
  # Add Group variable
  mutate(., Group = gsub("_.*","",Sample))

#save(data.lg, file = "DATA/Arabido_UPS/Arabido_UPS_1of3inall_ImpMLE_long")
# ----------------------------------------------------------- #
   
# Load post-imputation data in long format
load(file="DATA/Arabido_UPS/Arabido_UPS_1of3inall_ImpMLE_long")

# --- Find peptides of interest for illustration purpose --- #
## List of unique peptides with corresponding protein name
list.ID <- data.lg %>% select(starts_with(c("ID","Protein"))) %>% unique() 
## Print tibble of UPS peptides (differentially expressed)
list.ID %>% slice(grep("ups", list.ID$Protein))
## Print tibble of ARATH peptides (not differentially expressed)
list.ID %>% slice(grep("ARATH", list.ID$Protein))
list.ID %>% slice(intersect(grep("ARATH", list.ID$Protein),which(startsWith(list.ID$ID,"L"))))

## Function to plot the distribution in each condition
## for given peptides and conditions (groups)
plot.qData <- function(qData.lg, pept = NULL, cond = NULL, ...){
  if (is.null(pept)) {pept = levels(as.factor(qData.lg$ID))}
  if (is.null(cond)) {cond = levels(as.factor(qData.lg$Group))}
  df.plot = qData.lg %>% filter(ID %in% pept & Group %in% cond)
  print(df.plot)
  ggplot(df.plot, aes(x = Group, y= Intensity, color = Imp.Draw)) + 
    geom_boxplot() +
    theme_classic() +
    labs(color = "Imputation\nindex")
}

df.plot <- data.lg %>% filter(ID %in% "AALEELVK" & Group %in% c("Point1","Point7"))
plot.qData(data.lg, pept = "AALEELVK", cond = c("Point1","Point7"))
plot.qData(data.lg, pept = "AALEELVK") # UPS
plot.qData(data.lg, pept = "AAAVGANNQAAQSILK") # ARATH
plot.qData(data.lg, pept = "LAAEAYSIFR") # ARATH


load(file="DATA/Arabido_UPS/Arabido_UPS_1of3inall_ImpMLE_res_dapar")
load(file="DATA/Arabido_UPS/Arabido_UPS_1of3inall_ImpMLE_res_mi4limma")
library(cp4p)
df.pval <- data.frame(pval = adjust.p(p = res.mi4limma$P_Value$Point3_vs_Point7_pval)$adjp,
                      Seq = peptidesMQ_clean$Sequence, 
                      Prot = peptidesMQ_clean$Leading.razor.protein)
df.pval %>% arrange(pval.rawp) %>% head
df.pval %>% arrange(pval.rawp) %>% filter(Seq %in% c("AALEELVK","LAAEAYSIFR"))


which(list.ID$ID %in% c("AALEELVK","LAAEAYSIFR"))




# De A-tu

list_ID = data.lg$ID %>% unique()

data = data.lg %>%
  arrange(Imp.Draw,Sample, ID) %>% 
  mutate(Output = Intensity, Draw = Imp.Draw) %>% 
  select(- c(Intensity, Imp.Draw))

dataa = data %>% filter(ID %in% list_ID[c(13731)])
dim = dataa$ID %>% n_distinct()

res = post_mean_diff(
  data = dataa ,
  mu_0 = rep(0, dim), 
  lambda_0 = 1,
  Sigma_0 = diag(1, nrow = dim, ncol = dim),
  nu_0 = 1
)

gg1 = plot_dif(res, c('Point7', 'Point1'), peptide = 1) + xlim(c(-30,30))
gg4 = plot_dif(res, c('Point1', 'Point7'), peptide = 1) + xlim(c(-30,30))
gg5 = plot_dif(res, c('Point7', 'Point3'), peptide = 1) + xlim(c(-30,30))
gg2 = plot_dif(res, c('Point3', 'Point7'), peptide = 1) + xlim(c(-30,30))

gg2 = plot_dif(res, c('Point1'), 42)
gg3 = plot_dif(res, c('Point2'), 42)

grid.arrange(gg1, gg4, gg5, gg2, nrow = 4)

datab = data %>% filter(ID %in% list_ID[1:1000])

dimb = datab$ID %>% n_distinct()

resb = post_mean_diff(
  data = datab,
  mu_0 = rep(0, dimb), 
  lambda_0 = 1,
  Sigma_0 = diag(1, nrow = dimb, ncol = dimb),
  nu_0 = 10
)

gg1b = plot_dif(resb, c('Point7', 'Point1'), peptide = 42) + xlim(c(-30,30))
gg4b = plot_dif(resb, c('Point1', 'Point7'), peptide = 42) + xlim(c(-30,30))
gg5b = plot_dif(resb, c('Point7', 'Point3'), peptide = 42) + xlim(c(-30,30))
gg2b = plot_dif(resb, c('Point3', 'Point7'), peptide = 42) + xlim(c(-30,30))

grid.arrange(gg1b, gg4b, gg5b, gg2b, nrow = 4)



# datac = data %>%
#   arrange(Imp.Draw,Sample, ID) %>% 
#   mutate(Output = Intensity, Draw = Imp.Draw) %>% 
#   select(- c(Intensity, Imp.Draw))

dimc = datac$ID %>% n_distinct()

sub_datac = datac %>%
  filter(ID %in% list_ID[1:100])

resc = post_mean_diff(
  data = datac,
  mu_0 = rep(0, dimc), 
  lambda_0 = 1,
  Sigma_0 = diag(1, nrow = dimc, ncol = dimc),
  nu_0 = 1
)

## Univariate case ##
db = data %>% filter(Draw == 1) %>% select(- Draw)

sub_db = db %>%
  filter(ID %in% list_ID[1])

res_uni = post_mean_diff_uni(
  data = sub_db,
  mu_0 = 0, 
  lambda_0 = 1,
  beta_0 = 1,
  alpha_0 = 50
)

## Rajouter dans la fonction la boucle sur la totalité des peptides
