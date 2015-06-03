setwd("~/nemo_cluster")
Fst_Qst<-read.table("Fst_Qst.txt")
Fst_Qst_conti<-read.table("Fst_Qst_conti.txt")
library("ggplot2")

#On enlève les modèles à l'équilibre, le cas 3 et la fécondité de 10
Fst_Qst <- Fst_Qst[Fst_Qst$Model!="equilibrium",]
Fst_Qst <- Fst_Qst[Fst_Qst$Cas!=3 & Fst_Qst$Fecundity != 10,]
Fst_Qst_conti <- Fst_Qst_conti[Fst_Qst_conti$Model!="equilibrium",]
Fst_Qst_conti <- Fst_Qst_conti[Fst_Qst_conti$Cas!=3 & Fst_Qst_conti$Fecundity != 10,]

#On mets les paramètres en facteurs
Fst_Qst$Cas <- as.factor(Fst_Qst$Cas)
Fst_Qst$Selfing <- as.factor(Fst_Qst$Selfing)
Fst_Qst$QTLs <- as.factor(Fst_Qst$QTLs)
Fst_Qst$Fecundity <- as.factor(Fst_Qst$Fecundity)
Fst_Qst$Selection <- as.factor(Fst_Qst$Selection)
Fst_Qst$Repetition <- as.factor(Fst_Qst$Repetition)
Fst_Qst_conti$Cas <- as.factor(Fst_Qst_conti$Cas)
Fst_Qst_conti$Selfing <- as.factor(Fst_Qst_conti$Selfing)
Fst_Qst_conti$QTLs <- as.factor(Fst_Qst_conti$QTLs)
Fst_Qst_conti$Fecundity <- as.factor(Fst_Qst_conti$Fecundity)
Fst_Qst_conti$Selection <- as.factor(Fst_Qst_conti$Selection)
Fst_Qst_conti$Repetition <- as.factor(Fst_Qst_conti$Repetition)

##### On réordonne les cas 1, 2, 3, 4 en 2, 1, 3, 4 et on les renomme, on réordonne les sélections en 50,20,5,1 et les QTLs en 10,50 #####
Fst_Qst$Cas<-factor(Fst_Qst$Cas,levels=c("2","1","3","4")) ; levels(Fst_Qst$Cas) <- c("A","B","C","D")
Fst_Qst$Selection<-factor(Fst_Qst$Selection,levels=c("50","20","5","1"));levels(Fst_Qst$Selection)<-c("50","20","05","01")
Fst_Qst$QTLs<-factor(Fst_Qst$QTLs,levels=c("10","50"))
Fst_Qst$Selfing<-factor(Fst_Qst$Selfing,levels=c("0","3")) ; levels(Fst_Qst$Selfing) <- c("0","0.3")

Fst_Qst_conti$Cas<-factor(Fst_Qst_conti$Cas,levels=c("2","1","3","4")) ; levels(Fst_Qst_conti$Cas) <- c("A","B","C","D")
Fst_Qst_conti$Selection<-factor(Fst_Qst_conti$Selection,levels=c("50","20","5","1"));levels(Fst_Qst_conti$Selection)<-c("50","20","05","01")
Fst_Qst_conti$QTLs<-factor(Fst_Qst_conti$QTLs,levels=c("10","50"))
Fst_Qst_conti$Selfing<-factor(Fst_Qst_conti$Selfing,levels=c("0","3")) ; levels(Fst_Qst_conti$Selfing) <- c("0","0.3")
Fst_Qst_conti$Model<-factor(Fst_Qst_conti$Model,levels=c("2_conti","2 conti_nul","3_conti","3 conti_nul")) ; levels(Fst_Qst_conti$Model) <- c("2","2 nul","3","3 nul")

Case <- c('A'='Case A','B'='Case B','C'='Case C','D'='Case D'); Fecun <- c("010"="Mean Fecundity = 10","100"="Mean Fecundity = 100","200"="Mean Fecundity = 200")

#On change le système de contraintes par défaut
options(contrasts=c("contr.sum","contr.sum"))

#On fait des plot.design des différentes variables
par(mfrow = c(2,5))
pdf("plot.design de Fst_Qst")
plot.design(Fst_Qst[-c(2,7)])
dev.off()
pdf("plot.design de Fst_Qst_conti")
plot.design(Fst_Qst_conti[-c(2,7)])
dev.off()

anova_glm <- function(i,j) anova(glm(i~(Model+Selfing+QTLs+Fecundity+Selection)^2, data = j,family = quasibinomial(link = "logit")))
coeff_glm <- function(i,j) coefficients(glm(i~(Model+Selfing+QTLs+Fecundity+Selection)^2, data = j,family = quasibinomial(link = "logit")))

anova_df_Fst_Qst <- apply(Fst_Qst[c(8:15,20:22)],MARGIN = 2,FUN = anova_glm,Fst_Qst)
anova_df_Fst_Qst_conti <- apply(Fst_Qst_conti[c(8:15,20:22)],MARGIN = 2,FUN = anova_glm,Fst_Qst_conti)
coeff_df_Fst_Qst <- apply(Fst_Qst[c(8:15,20:22)],MARGIN = 2,FUN = coeff_glm,Fst_Qst)
coeff_df_Fst_Qst_conti <- apply(Fst_Qst_conti[c(8:15,20:22)],MARGIN = 2,FUN = coeff_glm,Fst_Qst_conti)


####Plots des statistiques de différentiation####
#Au niveau Environnement/population
for (f in c("100","200")){
  g <- ggplot(Fst_Qst[which(Fst_Qst$Fecundity==f&Fst_Qst$Model!="2 nul"&Fst_Qst$Model!="3 nul"),]) + facet_grid(Selfing~Model, labeller = label_both) + aes(color = QTLs)+ scale_shape_manual(values = c(1,2,3)) +
    stat_summary(aes(Selection, FstEnv_neutral, shape = "neutral Fst",size = 0.5), fun.y = "mean", geom = "point", position = position_dodge(width = 0.5)) +
    stat_summary(aes(Selection, FstEnv_neutral,group = QTLs, shape = "neutral Fst"), fun.y = "mean", geom = "line", position = position_dodge(width = 0.5),linetype = "dashed") +
    stat_summary(aes(Selection, FstEnv, shape = "Fst",size = 0.5), fun.y = "mean", geom = "point", position = position_dodge(width = 0.5)) +
    stat_summary(aes(Selection, FstEnv,group = QTLs, shape = "Fst"), fun.y = "mean", geom = "line", position = position_dodge(width = 0.5),linetype = "dashed") +
    stat_summary(aes(Selection, QSTenv, shape = "Qst",size = 0.5), fun.y = "mean", geom = "point", position = position_dodge(width = 0.5)) +
    stat_summary(aes(Selection, QSTenv,group = QTLs, shape = "Qst"), fun.y = "mean", geom = "line", position = position_dodge(width = 0.5),linetype = "dashed") +
    guides(shape = guide_legend(override.aes = list(size=2),title = "Statistics"), colour = guide_legend(override.aes = list(size=1),title = "Number of QTLs"), size = FALSE) +
    labs(title = paste("Diallelic mutation model, Environment level, mean fecundity =",f), y = "") +
    theme_bw() + theme(panel.grid=element_blank()) 
  print(g)
  pdf(paste("Diallelic_Differentiation_stat_Env_level_Fec",f,".pdf",sep=""), height = 12, width = 8, paper = "special")
  print(g)
  dev.off()
}

for (f in c("100","200")){
  g <- ggplot(Fst_Qst_conti[which(Fst_Qst_conti$Fecundity==f&Fst_Qst_conti$Model!="2 nul"&Fst_Qst_conti$Model!="3 nul"),]) + facet_grid(Selfing~Model, labeller = label_both) + aes(color = QTLs)+ scale_shape_manual(values = c(1,2,3)) +
    stat_summary(aes(Selection, FstEnv_neutral, shape = "neutral Fst",size = 0.5), fun.y = "mean", geom = "point", position = position_dodge(width = 0.5)) +
    stat_summary(aes(Selection, FstEnv_neutral,group = QTLs, shape = "neutral Fst"), fun.y = "mean", geom = "line", position = position_dodge(width = 0.5),linetype = "dashed") +
    stat_summary(aes(Selection, FstEnv, shape = "Fst",size = 0.5), fun.y = "mean", geom = "point", position = position_dodge(width = 0.5)) +
    stat_summary(aes(Selection, FstEnv,group = QTLs, shape = "Fst"), fun.y = "mean", geom = "line", position = position_dodge(width = 0.5),linetype = "dashed") +
    stat_summary(aes(Selection, QSTenv, shape = "Qst",size = 0.5), fun.y = "mean", geom = "point", position = position_dodge(width = 0.5)) +
    stat_summary(aes(Selection, QSTenv,group = QTLs, shape = "Qst"), fun.y = "mean", geom = "line", position = position_dodge(width = 0.5),linetype = "dashed") +
    guides(shape = guide_legend(override.aes = list(size=2),title = "Statistics"), colour = guide_legend(override.aes = list(size=1),title = "Number of QTLs"), size = FALSE) +
    labs(title = paste("Continuous mutation model, Environment level, mean fecundity =",f), y = "") +
    theme_bw() + theme(panel.grid=element_blank()) 
  pdf(paste("Continuous_Differentiation_stat_Env_level_Fec",f,".pdf",sep=""), height = 12, width = 8, paper = "special")
  print(g)
  dev.off()
}

#Au niveau population
for (f in c("100","200")){
  g <- ggplot(Fst_Qst[which(Fst_Qst$Fecundity==f),]) + facet_grid(Selfing~Model, labeller = label_both) + aes(color = QTLs)+ scale_shape_manual(values = c(1,2,3)) +
    stat_summary(aes(Selection, FstPop_neutral, shape = "neutral Fst",size = 0.5), fun.y = "mean", geom = "point", position = position_dodge(width = 0.5)) +
    stat_summary(aes(Selection, FstPop_neutral,group = QTLs, shape = "neutral Fst"), fun.y = "mean", geom = "line", position = position_dodge(width = 0.5),linetype = "dashed") +
    stat_summary(aes(Selection, FstPop, shape = "Fst",size = 0.5), fun.y = "mean", geom = "point", position = position_dodge(width = 0.5)) +
    stat_summary(aes(Selection, FstPop,group = QTLs, shape = "Fst"), fun.y = "mean", geom = "line", position = position_dodge(width = 0.5),linetype = "dashed") +
    stat_summary(aes(Selection, QSTpop, shape = "Qst",size = 0.5), fun.y = "mean", geom = "point", position = position_dodge(width = 0.5)) +
    stat_summary(aes(Selection, QSTpop,group = QTLs, shape = "Qst"), fun.y = "mean", geom = "line", position = position_dodge(width = 0.5),linetype = "dashed") +
    guides(shape = guide_legend(override.aes = list(size=2),title = "Statistics"), colour = guide_legend(override.aes = list(size=1),title = "Number of QTLs"), size = FALSE) +
    labs(title = paste("Diallelic mutation model, Population level, mean fecundity =",f), y = "") +
    theme_bw() + theme(panel.grid=element_blank()) 
  pdf(paste("Diallelic_Differentiation_stat_Pop_level_Fec",f,".pdf",sep=""), height = 12, width = 16, paper = "special")
  print(g)
  dev.off()
}

for (f in c("100","200")){
  g <- ggplot(Fst_Qst_conti[which(Fst_Qst_conti$Fecundity==f),]) + facet_grid(Selfing~Model, labeller = label_both) + aes(color = QTLs)+ scale_shape_manual(values = c(1,2,3)) +
    stat_summary(aes(Selection, FstPop_neutral, shape = "neutral Fst",size = 0.5), fun.y = "mean", geom = "point", position = position_dodge(width = 0.5)) +
    stat_summary(aes(Selection, FstPop_neutral,group = QTLs, shape = "neutral Fst"), fun.y = "mean", geom = "line", position = position_dodge(width = 0.5),linetype = "dashed") +
    stat_summary(aes(Selection, FstPop, shape = "Fst",size = 0.5), fun.y = "mean", geom = "point", position = position_dodge(width = 0.5)) +
    stat_summary(aes(Selection, FstPop,group = QTLs, shape = "Fst"), fun.y = "mean", geom = "line", position = position_dodge(width = 0.5),linetype = "dashed") +
    stat_summary(aes(Selection, QSTpop, shape = "Qst",size = 0.5), fun.y = "mean", geom = "point", position = position_dodge(width = 0.5)) +
    stat_summary(aes(Selection, QSTpop,group = QTLs, shape = "Qst"), fun.y = "mean", geom = "line", position = position_dodge(width = 0.5),linetype = "dashed") +
    guides(shape = guide_legend(override.aes = list(size=2),title = "Statistics"), colour = guide_legend(override.aes = list(size=1),title = "Number of QTLs"), size = FALSE) +
    labs(title = paste("Continuous mutation model, Population level, mean fecundity =",f), y = "") +
    theme_bw() + theme(panel.grid=element_blank()) 
  pdf(paste("Continuous_Differentiation_stat_Pop_level_Fec",f,".pdf",sep=""), height = 12, width = 16, paper = "special")
  print(g)
  dev.off()
}

#At the patch level
for (f in c("100","200")){
  g <- ggplot(Fst_Qst[which(Fst_Qst$Fecundity==f),]) + facet_grid(Selfing~Model, labeller = label_both) + aes(color = QTLs)+ ylim(0,1)+ scale_shape_manual(values = c(1,2,3)) +
    stat_summary(aes(Selection, FstPatch_neutral, shape = "neutral Fst",size = 0.5), fun.y = "mean", geom = "point", position = position_dodge(width = 0.5)) +
    stat_summary(aes(Selection, FstPatch_neutral,group = QTLs, shape = "neutral Fst"), fun.y = "mean", geom = "line", position = position_dodge(width = 0.5),linetype = "dashed") +
    stat_summary(aes(Selection, FstPatch, shape = "Fst",size = 0.5), fun.y = "mean", geom = "point", position = position_dodge(width = 0.5)) +
    stat_summary(aes(Selection, FstPatch,group = QTLs, shape = "Fst"), fun.y = "mean", geom = "line", position = position_dodge(width = 0.5),linetype = "dashed") +
    stat_summary(aes(Selection, QSTpatch, shape = "Qst",size = 0.5), fun.y = "mean", geom = "point", position = position_dodge(width = 0.5)) +
    stat_summary(aes(Selection, QSTpatch,group = QTLs, shape = "Qst"), fun.y = "mean", geom = "line", position = position_dodge(width = 0.5),linetype = "dashed") +
    guides(shape = guide_legend(override.aes = list(size=2),title = "Statistics"), colour = guide_legend(override.aes = list(size=1),title = "Number of QTLs"), size = FALSE) +
    labs(title = paste("Diallelic mutation model, Patch level, mean fecundity =",f), y = "") +
    theme_bw() + theme(panel.grid=element_blank()) 
  pdf(paste("Diallelic_Differentiation_stat_Patch_level_Fec",f,".pdf",sep=""), height = 12, width = 16, paper = "special")
  print(g)
  dev.off()
}

for (f in c("100","200")){
  g <- ggplot(Fst_Qst_conti[which(Fst_Qst_conti$Fecundity==f),]) + facet_grid(Selfing~Model, labeller = label_both) + aes(color = QTLs)+ylim(0,1)+ scale_shape_manual(values = c(1,2,3)) +
    stat_summary(aes(Selection, FstPatch_neutral, shape = "neutral Fst",size = 0.5), fun.y = "mean", geom = "point", position = position_dodge(width = 0.5)) +
    stat_summary(aes(Selection, FstPatch_neutral,group = QTLs, shape = "neutral Fst"), fun.y = "mean", geom = "line", position = position_dodge(width = 0.5),linetype = "dashed") +
    stat_summary(aes(Selection, FstPatch, shape = "Fst",size = 0.5), fun.y = "mean", geom = "point", position = position_dodge(width = 0.5)) +
    stat_summary(aes(Selection, FstPatch,group = QTLs, shape = "Fst"), fun.y = "mean", geom = "line", position = position_dodge(width = 0.5),linetype = "dashed") +
    stat_summary(aes(Selection, QSTpatch, shape = "Qst",size = 0.5), fun.y = "mean", geom = "point", position = position_dodge(width = 0.5)) +
    stat_summary(aes(Selection, QSTpatch,group = QTLs, shape = "Qst"), fun.y = "mean", geom = "line", position = position_dodge(width = 0.5),linetype = "dashed") +
    guides(shape = guide_legend(override.aes = list(size=2),title = "Statistics"), colour = guide_legend(override.aes = list(size=1),title = "Number of QTLs"), size = FALSE) +
    labs(title = paste("Continuous mutation model, Patch level, mean fecundity =",f), y = "") +
    theme_bw() + theme(panel.grid=element_blank()) 
  pdf(paste("Continuous_Differentiation_stat_Patch_level_Fec",f,".pdf",sep=""), height = 12, width = 16, paper = "special")
  print(g)
  dev.off()
}