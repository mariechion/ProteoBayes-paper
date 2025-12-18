# Load functions and packages
library(ProteoBayes)
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
  if (is.null(qData.lg$Intensity)) {qData.lg$Intensity = qData.lg$Output}
  df.plot = qData.lg %>% filter(ID %in% pept & Group %in% cond)
  print(df.plot)
  if (!is.null(qData.lg$Imp.Draw)) {
    ggplot2::ggplot(df.plot, aes(x = Group, y= Intensity, fill = Imp.Draw)) + 
      geom_boxplot(outlier.shape = NA) +
      theme_classic() +
      labs(fill = "Imputation\nindex")
  }
  if (is.null(qData.lg$Imp.Draw)) {
    ggplot2::ggplot(df.plot, aes(x = Group, y= Intensity, fill = Group)) + 
      geom_boxplot(outlier.shape = NA) +
      theme_classic() 
  }
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
