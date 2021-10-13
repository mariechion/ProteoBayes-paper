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
# Long format for raw (not imputed) data
db.raw <- peptidesMQ_clean %>%
  select(starts_with(c("Sequence","Leading", "Intensity."))) %>% 
  pivot_longer(cols = starts_with("Intensity."), 
               names_to = "Sample",
               names_prefix = "Intensity.",
               values_to = "Intensity") %>% 
  mutate(., Group = gsub("_.*","",Sample)) %>% 
  transmute(ID = Sequence, Protein = Leading.razor.protein,
            Sample, Output = log2(Intensity), Group)


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
  ggplot2::ggplot(df.plot, aes(x = Group, y= Intensity, fill = Imp.Draw)) + 
    geom_boxplot(outlier.shape = NA) +
    theme_classic() +
    labs(fill = "Imputation\nindex")
}

df.plot <- data.lg %>% filter(ID %in% "AALEELVK" & Group %in% c("Point1","Point7"))
plot.qData(data.lg, pept = "AALEELVK", cond = c("Point1","Point7"))
plot.qData(data.lg, pept = "AALEELVK") # UPS
plot.qData(data.lg, pept = "AAAVGANNQAAQSILK") # ARATH
plot.qData(data.lg, pept = "AAAVGANNQAAQSILK", cond = c("Point1","Point7")) # ARATH
plot.qData(data.lg, pept = "LAAEAYSIFR", cond = c("Point1","Point7")) # ARATH


load(file="DATA/Arabido_UPS/Arabido_UPS_1of3inall_ImpMLE_res_dapar")
load(file="DATA/Arabido_UPS/Arabido_UPS_1of3inall_ImpMLE_res_mi4limma")
library(cp4p)
df.pval <- data.frame(pval = adjust.p(p = res.mi4limma$P_Value$Point3_vs_Point7_pval)$adjp,
                      Seq = peptidesMQ_clean$Sequence, 
                      Prot = peptidesMQ_clean$Leading.razor.protein)
df.pval %>% arrange(pval.rawp) %>% head
df.pval %>% arrange(pval.rawp) %>% filter(Seq %in% c("AALEELVK","LAAEAYSIFR"))


which(list.ID$ID %in% c("AALEELVK","LAAEAYSIFR"))

data.plot = data.lg %>% filter(ID %in%  c("AALEELVK","LAAEAYSIFR"))
data.plot

data = data.plot %>%
  arrange(Imp.Draw,Sample, ID) %>% 
  mutate(Output = Intensity, Draw = Imp.Draw) %>% 
  select(- c(Intensity, Imp.Draw))

dim = data$ID %>% n_distinct()

res = post_mean_diff(
  data = data,
  mu_0 = rep(0, 2), 
  lambda_0 = 1,
  Sigma_0 = diag(1, nrow = 2, ncol = 2),
  nu_0 = 1
)

gg1 = plot_dif(res, c('Point7', 'Point1'), peptide = 1) + xlim(c(-30,30))
plot_dif(res, c('Point7', 'Point1'), peptide = 2) + xlim(c(-30,30))


data.lg %>% filter(Protein == "P12081ups|SYHC_HUMAN_UPS") %>% select(ID) %>% unique()
data.lg %>% filter(Protein == "P12081ups|SYHC_HUMAN_UPS") %>% select(ID) %>% n_distinct()
data.lg %>% filter(Protein == "sp|Q94A41|AMY3_ARATH") %>% select(ID) %>% unique()
data.lg %>% filter(Protein == "sp|Q94A41|AMY3_ARATH") %>% select(ID) %>% n_distinct()

data.lg %>% group_by(Protein) %>% summarise(n = n_distinct(ID)) %>% filter(n==9)

data.lg %>% filter(Protein == "P12081ups|SYHC_HUMAN_UPS") %>% select(ID) %>% unique()
data.lg %>% filter(Protein == "sp|F4I893|ILA_ARATH") %>% select(ID) %>% unique()


# P12081ups|SYHC_HUMAN_UPS
plot.qData(data.lg, 
           pept = "AALEELVK",
           cond = c("Point1","Point7"))
plot.qData(data.lg, 
           pept = "ASAELIEEEVAK",
           cond = c("Point1","Point7"))
plot.qData(data.lg, 
           pept = "DQGGELLSLR",
           cond = c("Point1","Point7"))
plot.qData(data.lg, 
           pept = "IFSIVEQR",
           cond = c("Point1","Point7"))
plot.qData(data.lg, 
           pept = "KVPCVGLSIGVER",
           cond = c("Point1","Point7"))
plot.qData(data.lg, 
           pept = "REDLVEEIK",
           cond = c("Point1","Point7"))
plot.qData(data.lg, 
           pept = "RHGAEVIDTPVFELK",
           cond = c("Point1","Point7"))
plot.qData(data.lg, 
           pept = "TICSSVDK",
           cond = c("Point1","Point7"))
plot.qData(data.lg, 
           pept = "VPCVGLSIGVER",
           cond = c("Point1","Point7"))

# sp|F4I893|ILA_ARATH
plot.qData(data.lg, 
           pept = "ALADPNTDVR",
           cond = c("Point1","Point7"))
plot.qData(data.lg, 
           pept = "EVQELAQEAAER",
           cond = c("Point1","Point7"))
plot.qData(data.lg, 
           pept = "LVLPSLLK",
           cond = c("Point1","Point7"))
plot.qData(data.lg, 
           pept = "QSSVELLGDLLFK",
           cond = c("Point1","Point7"))
plot.qData(data.lg, 
           pept = "SPIVSAAAFENLVK",
           cond = c("Point1","Point7"))
plot.qData(data.lg, 
           pept = "TDVSLSVR",
           cond = c("Point1","Point7"))
plot.qData(data.lg, 
           pept = "VLPLIIPILSK",
           cond = c("Point1","Point7"))
plot.qData(data.lg, 
           pept = "VVIDVLSSIVSALHDDSSEVR",
           cond = c("Point1","Point7"))
plot.qData(data.lg, 
           pept = "YALELLPVILPQAR",
           cond = c("Point1","Point7"))

# De A-tu

list_ID = data.lg$ID %>% unique()

data = data.lg %>%
  arrange(Imp.Draw,Sample, ID) %>% 
  mutate(Output = Intensity, Draw = Imp.Draw) %>% 
  select(- c(Intensity, Imp.Draw))

dataa = data %>% filter(ID %in% list_ID[c(13731)])
dim = dataa$ID %>% n_distinct()

nu_0 = rgamma(1,2, 0.1)

res = post_mean_diff(
  data = dataa,
  mu_0 = rep(20, dim), 
  lambda_0 = 1,
  Sigma_0 = diag(1, nrow = dim, ncol = dim),
  nu_0 = 1
)

pept = dataa$ID[1]

plot_dif(res, c('Point1'), peptide = pept)
plot_dif(res, c('Point7'), peptide = pept)

gg1 = plot_dif(res, c('Point7', 'Point1'), peptide = pept) + xlim(c(-10,10))
gg4 = plot_dif(res, c('Point1', 'Point7'), peptide = pept) + xlim(c(-10,10))
gg5 = plot_dif(res, c('Point4', 'Point3'), peptide = pept) + xlim(c(-10,10))
gg2 = plot_dif(res, c('Point3', 'Point4'), peptide = pept) + xlim(c(-10,10))

grid.arrange(gg1, gg4, gg5, gg2, nrow = 4)

## Univariate case ##
db = data %>% filter(Draw == 1) %>% select(- Draw)

sub_db = db %>%
  filter(ID %in% list_ID[1:100])

res_uni = post_mean_diff_uni(
  data = sub_db,
  mu_0 = 20, 
  lambda_0 = 1,
  beta_0 = 0,
  alpha_0 = 10
)

############ Graphsfor the article #################
db_pept_vary = db.raw %>% filter(ID == "WCAVSEHEATK") %>% drop_na()

## Graph illustration of the method 
res_graph1_ups = post_mean_diff_uni(
  data = db_pept_vary,
  mu_0 = 25, 
  lambda_0 = 1,
  beta_0 = 1,
  alpha_0 = 2
)

plot_dif(res_graph1_ups, c('Point1', 'Point7'), peptide = "WCAVSEHEATK") +
  xlim(c(-20,30))

