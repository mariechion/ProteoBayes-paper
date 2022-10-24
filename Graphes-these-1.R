load("Arabido_raw_data_lg")
load("Arabido_imp_data_lg")

source("bayes_proteo.R")

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
      #ylim(19,29) +
      theme_classic() 
  }
}



############ Graphs for the article #################
db_pept_ups = db.raw %>% filter(ID == "AALEELVK") %>% drop_na()

## Graph illustration of the method 
res_graph1_ups = post_mean_diff_uni(
  data = db_pept_ups,
  mu_0 = 25, 
  lambda_0 = 1,
  beta_0 = 1,
  alpha_0 = 2
)

gg1.1 <- plot_dif(res_graph1_ups, c('Point1', 'Point7'), peptide = "AALEELVK") +
  xlim(c(-8,5)) + ylim(c(0,0.7))
gg2.1 <- plot_dif(res_graph1_ups, c('Point7', 'Point1'), peptide = "AALEELVK") +
  xlim(c(-8,5)) + ylim(c(0,0.7))


db_pept_arath = db.raw %>% filter(ID == "VLPLIIPILSK") %>% drop_na()
## Graph illustration of the method 
res_graph1_arath = post_mean_diff_uni(
  data = db_pept_arath,
  mu_0 = 25, 
  lambda_0 = 1,
  beta_0 = 1,
  alpha_0 = 2
)

gg1.2 <- plot_dif(res_graph1_arath, c('Point1', 'Point7'), peptide = "VLPLIIPILSK") +
  xlim(c(-8,5)) + ylim(c(0,0.7))
gg2.2 <- plot_dif(res_graph1_arath, c('Point7', 'Point1'), peptide = "VLPLIIPILSK") +
  xlim(c(-8,5)) + ylim(c(0,0.7))

cowplot::plot_grid(gg1.1, gg1.2, gg2.1, gg2.2, nrow = 2, ncol =2)
 
#Add boxplots
gg3.1 <- plot.qData(db.raw, 
                    pept = "AALEELVK",
                    cond = c("Point1","Point7")) + ylim(21,29) + theme(legend.position = "none")
gg3.2 <- plot.qData(db.raw, 
                    pept = "VLPLIIPILSK",
                    cond = c("Point1","Point7")) + ylim(21,29) + theme(legend.position = "none")

cowplot::plot_grid(gg1.1, gg1.2, gg2.1, gg2.2, gg3.1, gg3.2, nrow = 3, ncol =2)

png('FIGURES/ch5_graph1.png', res = 600, width = 4800, height = 2400, units = "px")
cowplot::plot_grid(gg1.1, gg1.2, gg3.1, gg3.2, nrow = 2, ncol =2)
dev.off()

png('FIGURES/ch5_graph1-1.png', res = 600, width = 3200, height = 4600, units = "px")
cowplot::plot_grid(gg1.1, gg2.1, gg3.1, nrow = 3, ncol =1)
dev.off()

png('FIGURES/ch5_graph1-2.png', res = 600, width = 3200, height = 4600, units = "px")
cowplot::plot_grid(gg1.2, gg2.2, gg3.2, nrow = 3, ncol =1)
dev.off()

