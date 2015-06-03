## On lit les différents fichiers

setwd("~/nemo_cluster")
Fst_neutral <- read.table("FST_neutral_Sims2.txt")
Fst_neutral_conti <- read.table("FST_neutral_Sims2Conti.txt")
Fst_quanti <- read.table("FSTq_Sims2.txt")
Fst_quanti_conti <- read.table("FSTq_Sims2Conti.txt")
Qst_quanti <- read.table("Sims2/QSTs_Sims2.txt")
Qst_quanti_conti <- read.table("Sims2Conti/QSTs_Sims2Conti.txt")
Genic_variances <- read.table('Sims2/Genic_variances_components/main/normal')
Genic_variances_nul <- read.table('Sims2/Genic_variances_components/main/nul')
Genic_variances_conti <- read.table('Sims2Conti/Genic_variances_components/main/normal_conti')
Genic_variances_conti_nul <- read.table('Sims2Conti/Genic_variances_components/main/nul')


## On concatène les tables de variances géniques et on transforme en facteurs certains paramètres
Genic_variances <- rbind(Genic_variances,Genic_variances_nul) ; Genic_variances_conti <- rbind(Genic_variances_conti,Genic_variances_conti_nul)
Genic_variances$Model <- as.factor(Genic_variances$Model) ; Genic_variances_conti$Model <- as.factor(Genic_variances_conti$Model)
Genic_variances$Selection <- as.factor(Genic_variances$Selection) ; Genic_variances_conti$Selection <- as.factor(Genic_variances_conti$Selection)
Genic_variances$Fecundity <- as.factor(Genic_variances$Fecundity) ; Genic_variances_conti$Fecundity <- as.factor(Genic_variances_conti$Fecundity)

## On fait du vide
rm(Genic_variances_conti_nul, Genic_variances_nul)

## On harmonise les noms de variables
names(Fst_neutral)[1] <- "Model" ; names(Fst_neutral_conti)[1] <- "Model"
names(Fst_neutral)[8:11] <- c("FstPop_neutral","FstEnv_neutral","FstPatch_neutral","Fis_neutral")
names(Fst_neutral_conti)[8:11] <- c("FstPop_neutral","FstEnv_neutral","FstPatch_neutral","Fis_neutral")

## On harmonise le nom des paramètres
Fst_neutral$Model <- factor(Fst_neutral$Model,c("2_neutral","2nul_neutral","3_neutral","3nul_neutral","equilibrium_neutral"))
levels(Fst_neutral$Model) <- c("2", "2 nul", "3", "3 nul", "equilibrium")
Fst_neutral_conti$Model <- factor(Fst_neutral_conti$Model,c("2_conti_neutral","2 conti_nul_neutral","3_conti_neutral","3 conti_nul_neutral","conti_equilibrium_neutral"))
levels(Fst_neutral_conti$Model) <- c("2_conti","2 conti_nul","3_conti","3 conti_nul","conti_equilibrium")

Genic_variances$Model <- factor(Genic_variances$Model,c("2","2nul","3","3nul"))
levels(Genic_variances$Model) <- c("2", "2 nul", "3", "3 nul")
Genic_variances$Fecundity <- factor(Genic_variances$Fecundity,c("10","100","200"))
levels(Genic_variances$Fecundity) <- c("010","100","200")
Genic_variances$Selection <- factor(Genic_variances$Selection,c("1","5","20","50"))
levels(Genic_variances$Selection) <- c("01","05","20","50")

Genic_variances_conti$Model <- factor(Genic_variances_conti$Model,c("2_conti","2 conti_nul","3_conti","3 conti_nul"))
levels(Genic_variances_conti$Model) <- c("2_conti","2 conti_nul","3_conti","3 conti_nul")
Genic_variances_conti$Fecundity <- factor(Genic_variances_conti$Fecundity,c("10","100","200"))
levels(Genic_variances_conti$Fecundity) <- c("010","100","200")
Genic_variances_conti$Selection <- factor(Genic_variances_conti$Selection,c("1","5","20","50"))
levels(Genic_variances_conti$Selection) <- c("01","05","20","50")

## On merge tout
Fst_Qst <- merge(Fst_neutral,Fst_quanti) ; Fst_Qst <- merge(Fst_Qst,Qst_quanti) ; Fst_Qst <- merge(Fst_Qst,Genic_variances)
Fst_Qst_conti <- merge(Fst_neutral_conti,Fst_quanti_conti); Fst_Qst_conti <- merge(Fst_Qst_conti,Qst_quanti_conti) ; Fst_Qst_conti <- merge(Fst_Qst_conti,Genic_variances_conti)
summary(Fst_Qst)
summary(Fst_Qst_conti)

write.table(Fst_Qst,'Fst_Qst.txt')
write.table(Fst_Qst_conti,'Fst_Qst_conti.txt')
