# Packages
library(reshape)
library(tidyverse)

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
  # Add ID variable using Peptide Sequence from the quantitative data
  map(., ~mutate(.,ID = peptidesMQ_clean$Sequence)) %>% 
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

save(data.lg, file = "DATA/Arabido_UPS/Arabido_UPS_1of3inall_ImpMLE_long")
   
  




