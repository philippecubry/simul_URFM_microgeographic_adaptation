#définit le répertoire de travail
setwd("~/nemo_cluster/Sims2")

##### Lecture des fichiers #####
results<- read.table(file= "Number_of_Occupied_Patches.txt")
results$freq_occupied_patches <- results$nb_occupied_patches/24 ; results$nb_unoccupied_patches <- 24 - results$nb_occupied_patches
results$Fecundity <- as.factor(results$Fecundity);results$QTLs <- as.factor(results$QTLs);results$Selfing <- as.factor(results$Selfing);results$Selection <- as.factor(results$Selection)

######Et on fait des GLM######

#On change le système de contraintes par défaut
options(contrasts=c("contr.sum","contr.sum"))

#On fait le GLM sur l'ensemble des données
occupied_glm<-glm(data = results,formula = (nb_occupied_patches/24) ~ (Model+ QTLs +Selfing+Fecundity+Selection), weights = rep(24,times=length(nb_occupied_patches)), family = binomial(link = "logit"),start=NULL)
summary(occupied_glm)
aov_glm_occupied_patches <- anova(occupied_glm,test="Chisq")
drop_glm_occupied_patches <- drop1(occupied_glm, scope = (nb_occupied_patches/24) ~ (Model+ QTLs +Selfing+ Model/Cas+Fecundity+Selection+Model:QTLs+Model:Selfing+Model:Fecundity+Model:Selection+QTLs:Selfing+QTLs:Fecundity+QTLs:Selection+Selfing:Fecundity+Selfing:Selection+Fecundity:Selection+Model:Cas:QTLs+Model:Cas:Selfing), test = "Chisq")
coeff_glm_occupied_patches <- dummy.coef(occupied_glm)

write.table(aov_glm_occupied_patches,file = "aov_glm_occupied_patches.txt")
write.table(drop_glm_occupied_patches,file = "drop_glm_occupied_patches.txt")

for (i in 1:length(coeff_glm_occupied_patches_conti)) {
  write.table(coeff_glm_occupied_patches[i], file = paste("glm_coeff_",names(coeff_glm_occupied_patches[i]),".txt",sep=""))
}

#On fait un GLM sur les modèles "nuls"
occupied_glm_nul<-glm(data = results[which(results$Model=="Overlapping without within population environmental heterogeneity"|results$Model == "Divided without within population environmental heterogeneity"),],
                      formula = (nb_occupied_patches/24) ~ (QTLs+Selfing+Cas+Cas/Model+Fecundity+Selection)^2,
                      weights = rep(24,times=length(nb_occupied_patches)),
                      family = binomial(link = "logit"),start=NULL)
summary(occupied_glm_nul)
aov_glm_occupied_patches_nul <- anova(occupied_glm_nul,test="Chisq")
drop_glm_occupied_patches_nul <- drop1(occupied_glm_nul, scope = ~ QTLs+Selfing+Cas+Cas/Model+Fecundity+Selection, test = "Chisq")
coeff_glm_occupied_patches_nul <- dummy.coef(occupied_glm_nul)

write.table(aov_glm_occupied_patches_nul,file = "aov_glm_occupied_patches_nul.txt")
write.table(drop_glm_occupied_patches_nul,file = "drop_glm_occupied_patches_nul.txt")

for (i in 1:length(coeff_glm_occupied_patches_nul)) {
  write.table(coeff_glm_occupied_patches_nul[i], file = paste("glm_coeff_nul_",names(coeff_glm_occupied_patches_nul[i]),".txt",sep=""))
}

#On fait un GLM sur les modèles avec hétérogénéité
occupied_glm_het<-glm(results[which(results$Model=="Overlapping with within population environmental heterogeneity"|results$Model == "Divided with within population environmental heterogeneity"),],
                      formula = (nb_occupied_patches/24) ~ (Model+ Model/Cas+QTLs+Selfing+Fecundity+Selection)^2,
                      weights = rep(24,times=length(nb_occupied_patches)),
                      family = binomial(link = "logit"),start=NULL)
summary(occupied_glm_het)
aov_glm_occupied_patches_het <- anova(occupied_glm_het,test="Chisq")
drop_glm_occupied_patches_het <- drop1(occupied_glm_het, scope = ~ (Model), test = "Chisq")
drop_glm_occupied_patches_het_more <- drop1(occupied_glm_het, test = "Chisq")
coeff_glm_occupied_patches_het <- dummy.coef(occupied_glm_het)

write.table(aov_glm_occupied_patches_het,file = "aov_glm_occupied_patches_het.txt")
write.table(drop_glm_occupied_patches_het,file = "drop_glm_occupied_patches_het.txt")

for (i in 1:length(coeff_glm_occupied_patches_het)) {
  write.table(coeff_glm_occupied_patches_het[i], file = paste("glm_coeff_het_",names(coeff_glm_occupied_patches_het[i]),".txt",sep=""))
}

#On fait un GLM sur les modèles avec divided
occupied_glm_div<-glm(results[which(results$Model== "Divided without within population environmental heterogeneity"|results$Model == "Divided with within population environmental heterogeneity"),],
                      formula = (nb_occupied_patches/24) ~ (Cas+Cas/Model+QTLs+Selfing+Fecundity+Selection)^2,
                      weights = rep(24,times=length(nb_occupied_patches)),
                      family = binomial(link = "logit"),start=NULL)
summary(occupied_glm_div)
aov_glm_occupied_patches_div <- anova(occupied_glm_div,test="Chisq")
drop_glm_occupied_patches_div <- drop1(occupied_glm_div, scope= ~(Cas + Cas/Model + QTLs + Selfing + Fecundity + Selection), test = "Chisq")
coeff_glm_occupied_patches_div <- dummy.coef(occupied_glm_div)

write.table(aov_glm_occupied_patches_div,file = "aov_glm_occupied_patches_div.txt")
write.table(drop_glm_occupied_patches_div,file = "drop_glm_occupied_patches_div.txt")

for (i in 1:length(coeff_glm_occupied_patches_div)) {
  write.table(coeff_glm_occupied_patches_div[i], file = paste("glm_coeff_div_",names(coeff_glm_occupied_patches_div[i]),".txt",sep=""))
}

#On fait un GLM sur les modèles avec overlapping
occupied_glm_over<-glm(results[which(results$Model== "Overlapping without within population environmental heterogeneity"|results$Model == "Overlapping with within population environmental heterogeneity"),],
                       formula = (nb_occupied_patches/24) ~ (Cas+Cas/Model+QTLs+Selfing+Fecundity+Selection)^2,
                       weights = rep(24,times=length(nb_occupied_patches)),
                       family = binomial(link = "logit"),start=NULL)
summary(occupied_glm_over)
aov_glm_occupied_patches_over <- anova(occupied_glm_over,test="Chisq")
drop_glm_occupied_patches_over <- drop1(occupied_glm_over, ~Model+Model/Cas+QTLs+Selfing+Fecundity+Selection, test = "Chisq")
coeff_glm_occupied_patches_over <- dummy.coef(occupied_glm_over)

write.table(aov_glm_occupied_patches_over,file = "aov_glm_occupied_patches_over.txt")
write.table(drop_glm_occupied_patches_over,file = "drop_glm_occupied_patches_over.txt")

for (i in 1:length(coeff_glm_occupied_patches_over)) {
  write.table(coeff_glm_occupied_patches_over[i], file = paste("glm_coeff_over_",names(coeff_glm_occupied_patches_over[i]),".txt",sep=""))
}