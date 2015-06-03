#définit le répertoire de travail
setwd("~/nemo_cluster/Sims2")

#Charge les bibliothèques nécessaires
library("ggplot2", lib.loc="/home/pcubry/R/x86_64-pc-linux-gnu-library/3.1")
library("Hmisc", lib.loc="/home/pcubry/R/x86_64-pc-linux-gnu-library/3.1")

#Définit une fonction basée sur les packages ggplot2 et Hmisc pour le plot de la moyenne et de l'écart-type
stat_sum_df <- function(fun, geom="crossbar", ...) {
  stat_summary(fun.data=fun, geom=geom, width=0.2, ...)
}

#Définit les variables utilisées dans les simulations et le nom des fichiers

model <- c("2","3") #On ne considère plus le modèle 1 qui est un cas intermédiaire entre les autres et présente peu d'intérêt
cas <- c("1","2","3","4")
fecundity <- c("010","100","200")
selfing <- c("0","03")
qtls <- c("10","50")
selection <- c("01","05","20","50")

#Création d'un objet vide qui va accueillir les résultats
results = NULL

#Boucle de lecture des fichiers, pour l'instant on lit le modèle 2
for(m in model) {
  for(s in selfing) {
    for(q in qtls) {
      for(f in fecundity) {
        for(w in selection) {
          for(c in cas) {
            temp <- read.delim(file=(paste("Sim2_Model",m,"_",s,"self_",q,"QTLs_fitperpop","/stats/Sim2_Model",m,"_",s,"self_",q,"QTLs_Cas",c,"_SelInt",w,"_fec",f,".txt",sep = "")),header=T)
            #On ne garde que ce qui nous intéresse
            temp <- temp[,c("replicate","generation","adlt.nbr",
                            #"resid.p1","resid.p2","resid.p3","resid.p4","resid.p5","resid.p6","resid.p7","resid.p8",
                            #"resid.p9","resid.p10","resid.p11","resid.p12","resid.p13","resid.p14","resid.p15","resid.p16",
                            #"resid.p17","resid.p18","resid.p19","resid.p20","resid.p21","resid.p22","resid.p23","resid.p24",
                            "adlt.allnbp","adlt.allnb","adlt.fixlocp","adlt.fixloc","adlt.ho","adlt.hsnei","adlt.htnei","adlt.fis","adlt.fst","adlt.fst.WH",
                            "adlt.q1","adlt.q1.Va","adlt.q1.Vb","adlt.q1.Vp","adlt.q1.Qst","fitness.mean",
                            "adlt.q1.p1","adlt.q1.p2","adlt.q1.p3","adlt.q1.p4","adlt.q1.p5","adlt.q1.p6","adlt.q1.p7","adlt.q1.p8",
                            "adlt.q1.p9","adlt.q1.p10","adlt.q1.p11","adlt.q1.p12","adlt.q1.p13","adlt.q1.p14","adlt.q1.p15","adlt.q1.p16",
                            "adlt.q1.p17","adlt.q1.p18","adlt.q1.p19","adlt.q1.p20","adlt.q1.p21","adlt.q1.p22","adlt.q1.p23","adlt.q1.p24",
                            "adlt.fem.p1","adlt.fem.p2","adlt.fem.p3","adlt.fem.p4","adlt.fem.p5","adlt.fem.p6","adlt.fem.p7","adlt.fem.p8",
                            "adlt.fem.p9","adlt.fem.p10","adlt.fem.p11","adlt.fem.p12","adlt.fem.p13","adlt.fem.p14","adlt.fem.p15","adlt.fem.p16",
                            "adlt.fem.p17","adlt.fem.p18","adlt.fem.p19","adlt.fem.p20","adlt.fem.p21","adlt.fem.p22","adlt.fem.p23","adlt.fem.p24",
                            "adlt.Va.q1.p1","adlt.Va.q1.p2","adlt.Va.q1.p3","adlt.Va.q1.p4","adlt.Va.q1.p5","adlt.Va.q1.p6","adlt.Va.q1.p7","adlt.Va.q1.p8",
                            "adlt.Va.q1.p9","adlt.Va.q1.p10","adlt.Va.q1.p11","adlt.Va.q1.p12","adlt.Va.q1.p13","adlt.Va.q1.p14","adlt.Va.q1.p15","adlt.Va.q1.p16",
                            "adlt.Va.q1.p17","adlt.Va.q1.p18","adlt.Va.q1.p19","adlt.Va.q1.p20","adlt.Va.q1.p21","adlt.Va.q1.p22","adlt.Va.q1.p23","adlt.Va.q1.p24",
                            #"adlt.Vp.q1.p1","adlt.Vp.q1.p2","adlt.Vp.q1.p3","adlt.Vp.q1.p4","adlt.Vp.q1.p5","adlt.Vp.q1.p6","adlt.Vp.q1.p7","adlt.Vp.q1.p8",
                            #"adlt.Vp.q1.p9","adlt.Vp.q1.p10","adlt.Vp.q1.p11","adlt.Vp.q1.p12","adlt.Vp.q1.p13","adlt.Vp.q1.p14","adlt.Vp.q1.p15","adlt.Vp.q1.p16",
                            #"adlt.Vp.q1.p17","adlt.Vp.q1.p18","adlt.Vp.q1.p19","adlt.Vp.q1.p20","adlt.Vp.q1.p21","adlt.Vp.q1.p22","adlt.Vp.q1.p23","adlt.Vp.q1.p24",
                            "adlt.W.avg.p1","adlt.W.avg.p2","adlt.W.avg.p3","adlt.W.avg.p4","adlt.W.avg.p5","adlt.W.avg.p6","adlt.W.avg.p7",
                            "adlt.W.avg.p8","adlt.W.avg.p9","adlt.W.avg.p10","adlt.W.avg.p11","adlt.W.avg.p12","adlt.W.avg.p13","adlt.W.avg.p14","adlt.W.avg.p15",
                            "adlt.W.avg.p16","adlt.W.avg.p17","adlt.W.avg.p18","adlt.W.avg.p19","adlt.W.avg.p20","adlt.W.avg.p21","adlt.W.avg.p22","adlt.W.avg.p23",
                            "adlt.W.avg.p24"
                            )]
            #Définition des noms de paramètres
            Model <- rep(m,length(temp[,1]))
            Selfing <- rep(s,length(temp[,1]))
            QTLs <- rep(q,length(temp[,1]))
            Fecundity <- rep(f,length(temp[,1]))
            Selection <- rep(w,length(temp[,1]))
            Cas <- rep(c,length(temp[,1]))
            #On ajoute les valeurs des paramètres à temp
            temp <- cbind(Model,Cas,Selfing,QTLs,Fecundity,Selection,temp) 
            #et on ajoute temp aux résultats !!!
            results <- rbind(temp,results)
          }  
        }
      }
    }
  }
}

#On transforme certaines variables en facteurs et on réordonne les cas 1, 2, 3, 4 en 2, 1, 3, 4
results$replicate<-as.factor(results$replicate)
results$Cas<-factor(results$Cas,levels=c("2","1","3","4"))

#Définition des Zopts
Z_Model2_Cas1_10QTLs <- c(-13.557282933,-4.519094311,-7.53182385166667,1.5063647703,13.557282933,4.519094311,7.53182385166667,-1.5063647703)
Z_Model2_Cas2_10QTLs <- c(-16.9466036663,-5.6488678888,-9.4147798146,1.8829559629,16.9466036663,5.6488678888,9.4147798146,-1.8829559629)
Z_Model2_Cas3_10QTLs <- c(-6.7786414665,-2.2595471555,-3.76591192583333,0.7531823852,6.7786414665,2.2595471555,3.76591192583333,-0.7531823852)
Z_Model2_Cas4_10QTLs <- c(-3.38932073325,-1.12977357775,-1.88295596291667,0.3765911926,3.38932073325,1.12977357775,1.88295596291667,-0.3765911926)
Z_Model2_Cas1_50QTLs <- c(-80.325348195,-26.775116065,-44.6251934417,8.92503868833333,80.325348195,26.775116065,44.6251934417,-8.92503868833333)
Z_Model2_Cas2_50QTLs <- c(-100.40668524375,-33.46889508125,-55.7814918020833,11.1562983604167,100.40668524375,33.46889508125,55.7814918020833,-11.1562983604167)
Z_Model2_Cas3_50QTLs <- c(-40.1626740975,-13.3875580325,-22.3125967208333,4.46251934416667,40.1626740975,13.3875580325,22.3125967208333,-4.46251934416667)
Z_Model2_Cas4_50QTLs <- c(-20.08133704875,-6.69377901625,-11.1562983604167,2.23125967208333,20.08133704875,6.69377901625,11.1562983604167,-2.23125967208333)

Z_Model3_Cas1_10QTLs <- c(-13.557282933,-9.68377352357143,-5.81026411414286,-1.93675470471429,13.557282933,9.68377352357143,5.81026411414286,1.93675470471429)
Z_Model3_Cas2_10QTLs <- c(-16.94660366625,-12.1047169044643,-7.26283014267857,-2.42094338089286,16.94660366625,12.1047169044643,7.26283014267857,2.42094338089286)
Z_Model3_Cas3_10QTLs <- c(-6.7786414665,-4.84188676178571,-2.90513205707143,-0.968377352357143,6.7786414665,4.84188676178571,2.90513205707143,0.968377352357143)
Z_Model3_Cas4_10QTLs <- c(-3.38932073325,-2.42094338089286,-1.45256602853571,-0.484188676178571,3.38932073325,2.42094338089286,1.45256602853571,0.484188676178571)
Z_Model3_Cas1_50QTLs <- c(-80.325348195,-57.3752487107143,-34.4251492264286,-11.4750497421429,80.325348195,57.3752487107143,34.4251492264286,11.4750497421429)
Z_Model3_Cas2_50QTLs <- c(-100.40668524375,-71.7190608883929,-43.0314365330357,-14.3438121776786,100.40668524375,71.7190608883929,43.0314365330357,14.3438121776786)
Z_Model3_Cas3_50QTLs <- c(-40.1626740975,-28.6876243553571,-17.2125746132143,-5.73752487107143,40.1626740975,28.6876243553571,17.2125746132143,5.73752487107143)
Z_Model3_Cas4_50QTLs <- c(-20.08133704875,-14.3438121776786,-8.60628730660714,-2.86876243553571,20.08133704875,14.3438121776786,8.60628730660714,2.86876243553571)

#On convertit les NaN en 0 dans les valeurs phénotypiques réalisées
results[,c("adlt.q1.p1","adlt.q1.p2","adlt.q1.p3","adlt.q1.p4","adlt.q1.p5","adlt.q1.p6","adlt.q1.p7","adlt.q1.p8",
           "adlt.q1.p9","adlt.q1.p10","adlt.q1.p11","adlt.q1.p12","adlt.q1.p13","adlt.q1.p14","adlt.q1.p15","adlt.q1.p16",
           "adlt.q1.p17","adlt.q1.p18","adlt.q1.p19","adlt.q1.p20","adlt.q1.p21","adlt.q1.p22","adlt.q1.p23","adlt.q1.p24",
           "adlt.Va.q1.p1","adlt.Va.q1.p2","adlt.Va.q1.p3","adlt.Va.q1.p4","adlt.Va.q1.p5","adlt.Va.q1.p6","adlt.Va.q1.p7","adlt.Va.q1.p8",
           "adlt.Va.q1.p9","adlt.Va.q1.p10","adlt.Va.q1.p11","adlt.Va.q1.p12","adlt.Va.q1.p13","adlt.Va.q1.p14","adlt.Va.q1.p15","adlt.Va.q1.p16",
           "adlt.Va.q1.p17","adlt.Va.q1.p18","adlt.Va.q1.p19","adlt.Va.q1.p20","adlt.Va.q1.p21","adlt.Va.q1.p22","adlt.Va.q1.p23","adlt.Va.q1.p24")]<-
  rapply( results[,c("adlt.q1.p1","adlt.q1.p2","adlt.q1.p3","adlt.q1.p4","adlt.q1.p5","adlt.q1.p6","adlt.q1.p7","adlt.q1.p8",
                            "adlt.q1.p9","adlt.q1.p10","adlt.q1.p11","adlt.q1.p12","adlt.q1.p13","adlt.q1.p14","adlt.q1.p15","adlt.q1.p16",
                            "adlt.q1.p17","adlt.q1.p18","adlt.q1.p19","adlt.q1.p20","adlt.q1.p21","adlt.q1.p22","adlt.q1.p23","adlt.q1.p24","adlt.Va.q1.p1","adlt.Va.q1.p2","adlt.Va.q1.p3","adlt.Va.q1.p4","adlt.Va.q1.p5","adlt.Va.q1.p6","adlt.Va.q1.p7","adlt.Va.q1.p8",
                            "adlt.Va.q1.p9","adlt.Va.q1.p10","adlt.Va.q1.p11","adlt.Va.q1.p12","adlt.Va.q1.p13","adlt.Va.q1.p14","adlt.Va.q1.p15","adlt.Va.q1.p16",
                            "adlt.Va.q1.p17","adlt.Va.q1.p18","adlt.Va.q1.p19","adlt.Va.q1.p20","adlt.Va.q1.p21","adlt.Va.q1.p22","adlt.Va.q1.p23","adlt.Va.q1.p24")], f=function(x) ifelse(is.nan(x),0,x), how="replace" )

# #Création des fichiers de données à traiter par ggplot
# for (m in model) {
#   for (f in fecundity) {
#     for (c in cas) {
#       for (q in qtls) {
#         assign(paste("Model",m,"_Fecundity",f,"_Case",c,"_",q,"QTLs",sep =""),ggplot(results[results$Model==m&
#                               results$generation>=1 & results$generation<=200  &
#                               results$Cas==c &
#                               results$Fecundity ==f &
#                               results$QTLs == q
#                               ,]))
#       }
#     }
#   }
# }
# 
# #Création des graphiques moyennes de Q par environnement
# for (m in model) {
#   for (f in fecundity) {
#     for (c in cas) {
#       for (q in qtls) {
#         #Création du graphe Mean Q per environment
#         assign(paste("Model",m,"_Fecundity",f,"_Case",c,"_",q,"QTLs","_MeanQperEnv",sep =""), eval(parse(text=paste("Model",m,"_Fecundity",f,"_Case",c,"_",q,"QTLs",sep=""))) + aes(generation) +
#           labs(y="", title = paste("Model ",m,"/Fecundity ",f,"/Case ",c,"/",q," QTLs - mean Q per environment",sep=""))+
#           facet_grid(Selfing~Selection, labeller = label_both) +
#            ylim(eval(parse(text=paste("Z_Model",m,"_Cas",c,"_",q,"QTLs",sep="")))[1]-1,-eval(parse(text=paste("Z_Model",m,"_Cas",c,"_",q,"QTLs",sep="")))[1]+1) +
#            geom_point(aes(y=((adlt.q1.p1*adlt.fem.p1+ adlt.q1.p2*adlt.fem.p2+adlt.q1.p3*adlt.fem.p3)/(adlt.fem.p1+adlt.fem.p2+adlt.fem.p3)),color="Pop1-Env1"))+
#            geom_point(aes(y=((adlt.q1.p4*adlt.fem.p4+ adlt.q1.p5*adlt.fem.p5+adlt.q1.p6*adlt.fem.p6)/(adlt.fem.p4+adlt.fem.p5+adlt.fem.p6)),color="Pop1-Env2"))+
#            geom_point(aes(y=((adlt.q1.p7*adlt.fem.p7+ adlt.q1.p8*adlt.fem.p8+adlt.q1.p9*adlt.fem.p9)/(adlt.fem.p7+adlt.fem.p8+adlt.fem.p9)),color="Pop2-Env3"))+
#            geom_point(aes(y=((adlt.q1.p10*adlt.fem.p10+ adlt.q1.p11*adlt.fem.p11+adlt.q1.p12*adlt.fem.p12)/(adlt.fem.p10+adlt.fem.p11+adlt.fem.p12)),color="Pop2-Env4"))+
#            geom_point(aes(y=((adlt.q1.p13*adlt.fem.p13+ adlt.q1.p14*adlt.fem.p14+adlt.q1.p15*adlt.fem.p15)/(adlt.fem.p13+adlt.fem.p14+adlt.fem.p15)),color="Pop3-Env5"))+
#            geom_point(aes(y=((adlt.q1.p16*adlt.fem.p16+ adlt.q1.p17*adlt.fem.p17+adlt.q1.p18*adlt.fem.p18)/(adlt.fem.p16+adlt.fem.p17+adlt.fem.p18)),color="Pop3-Env6"))+
#            geom_point(aes(y=((adlt.q1.p19*adlt.fem.p19+ adlt.q1.p20*adlt.fem.p20+adlt.q1.p21*adlt.fem.p21)/(adlt.fem.p19+adlt.fem.p20+adlt.fem.p21)),color="Pop4-Env7"))+
#            geom_point(aes(y=((adlt.q1.p22*adlt.fem.p22+ adlt.q1.p23*adlt.fem.p23+adlt.q1.p24*adlt.fem.p24)/(adlt.fem.p22+adlt.fem.p23+adlt.fem.p24)),color="Pop4-Env8"))+
#            geom_abline(intercept=eval(parse(text=paste("Z_Model",m,"_Cas",c,"_",q,"QTLs",sep="")))[1],slope= 0,size=0.5,aes(colour="Pop1-Env1"))+
#            geom_abline(intercept=eval(parse(text=paste("Z_Model",m,"_Cas",c,"_",q,"QTLs",sep="")))[2],slope= 0,size=0.5,aes(colour="Pop1-Env2"))+
#            geom_abline(intercept=eval(parse(text=paste("Z_Model",m,"_Cas",c,"_",q,"QTLs",sep="")))[3],slope= 0,size=0.5,aes(colour="Pop2-Env3"))+
#            geom_abline(intercept=eval(parse(text=paste("Z_Model",m,"_Cas",c,"_",q,"QTLs",sep="")))[4],slope= 0,size=0.5,aes(colour="Pop2-Env4"))+
#            geom_abline(intercept=eval(parse(text=paste("Z_Model",m,"_Cas",c,"_",q,"QTLs",sep="")))[5],slope= 0,size=0.5,aes(colour="Pop4-Env8"))+
#            geom_abline(intercept=eval(parse(text=paste("Z_Model",m,"_Cas",c,"_",q,"QTLs",sep="")))[6],slope= 0,size=0.5,aes(colour="Pop4-Env7"))+
#            geom_abline(intercept=eval(parse(text=paste("Z_Model",m,"_Cas",c,"_",q,"QTLs",sep="")))[7],slope= 0,size=0.5,aes(colour="Pop3-Env6"))+
#            geom_abline(intercept=eval(parse(text=paste("Z_Model",m,"_Cas",c,"_",q,"QTLs",sep="")))[8],slope= 0,size=0.5,aes(colour="Pop3-Env5"))+
#           theme_classic(base_size = 20, base_family = ""))
#         
#     #Ecriture du pdf
#         pdf(file = paste("Mod",m,"Fec",f,"Cas",c,"QTLs",q,"Qperenv.pdf",sep=""), height = 12, width = 16, paper = "special")
#         print(eval(parse(text=paste("Model",m,"_Fecundity",f,"_Case",c,"_",q,"QTLs","_MeanQperEnv",sep =""))))
#         dev.off()  
#       
#       }
#     }
#   }
# }
# 
# #Création des graphiques moyennes de Q par population
# for (m in model) {
#   for (f in fecundity) {
#     for (c in cas) {
#       for (q in qtls) {
#         #Création du graphe Mean Q per population
#         assign(paste("Model",m,"_Fecundity",f,"_Case",c,"_",q,"QTLs","_MeanQperPop",sep =""), eval(parse(text=paste("Model",m,"_Fecundity",f,"_Case",c,"_",q,"QTLs",sep=""))) + aes(generation) +
#                  labs(y="", title = paste("Model ",m,"/Fecundity ",f,"/Case ",c,"/",q," QTLs - mean Q per population",sep=""))+
#                  facet_grid(Selfing~Selection, labeller = label_both) +
#                  ylim(eval(parse(text=paste("Z_Model",m,"_Cas",c,"_",q,"QTLs",sep="")))[1]-1,-eval(parse(text=paste("Z_Model",m,"_Cas",c,"_",q,"QTLs",sep="")))[1]+1) +
#                  geom_point(aes(y=((adlt.q1.p1*adlt.fem.p1+ adlt.q1.p2*adlt.fem.p2+adlt.q1.p3*adlt.fem.p3+adlt.q1.p4*adlt.fem.p4+ adlt.q1.p5*adlt.fem.p5+adlt.q1.p6*adlt.fem.p6)/
#                                      (adlt.fem.p1+adlt.fem.p2+adlt.fem.p3+adlt.fem.p4+adlt.fem.p5+adlt.fem.p6)),color="Pop1"))+
#                  geom_point(aes(y=((adlt.q1.p7*adlt.fem.p7+ adlt.q1.p8*adlt.fem.p8+adlt.q1.p9*adlt.fem.p9+adlt.q1.p10*adlt.fem.p10+ adlt.q1.p11*adlt.fem.p11+adlt.q1.p12*adlt.fem.p12)/
#                                      (adlt.fem.p7+adlt.fem.p8+adlt.fem.p9+adlt.fem.p10+adlt.fem.p11+adlt.fem.p12)),color="Pop2"))+
#                  geom_point(aes(y=((adlt.q1.p13*adlt.fem.p13+ adlt.q1.p14*adlt.fem.p14+adlt.q1.p15*adlt.fem.p15+adlt.q1.p16*adlt.fem.p16+ adlt.q1.p17*adlt.fem.p17+adlt.q1.p18*adlt.fem.p18)/
#                                      (adlt.fem.p13+adlt.fem.p14+adlt.fem.p15+adlt.fem.p16+adlt.fem.p17+adlt.fem.p18)),color="Pop3"))+
#                  geom_point(aes(y=((adlt.q1.p19*adlt.fem.p19+ adlt.q1.p20*adlt.fem.p20+adlt.q1.p21*adlt.fem.p21+adlt.q1.p22*adlt.fem.p22+ adlt.q1.p23*adlt.fem.p23+adlt.q1.p24*adlt.fem.p24)/
#                                      (adlt.fem.p19+adlt.fem.p20+adlt.fem.p21+adlt.fem.p22+adlt.fem.p23+adlt.fem.p24)),color="Pop4"))+
#                  geom_abline(intercept=eval(parse(text=paste("Z_Model",m,"_Cas",c,"_",q,"QTLs",sep="")))[1],slope= 0,size=0.5,aes(colour="Pop1"))+
#                  geom_abline(intercept=eval(parse(text=paste("Z_Model",m,"_Cas",c,"_",q,"QTLs",sep="")))[2],slope= 0,size=0.5,aes(colour="Pop1"))+
#                  geom_abline(intercept=eval(parse(text=paste("Z_Model",m,"_Cas",c,"_",q,"QTLs",sep="")))[3],slope= 0,size=0.5,aes(colour="Pop2"))+
#                  geom_abline(intercept=eval(parse(text=paste("Z_Model",m,"_Cas",c,"_",q,"QTLs",sep="")))[4],slope= 0,size=0.5,aes(colour="Pop2"))+
#                  geom_abline(intercept=eval(parse(text=paste("Z_Model",m,"_Cas",c,"_",q,"QTLs",sep="")))[5],slope= 0,size=0.5,aes(colour="Pop4"))+
#                  geom_abline(intercept=eval(parse(text=paste("Z_Model",m,"_Cas",c,"_",q,"QTLs",sep="")))[6],slope= 0,size=0.5,aes(colour="Pop4"))+
#                  geom_abline(intercept=eval(parse(text=paste("Z_Model",m,"_Cas",c,"_",q,"QTLs",sep="")))[7],slope= 0,size=0.5,aes(colour="Pop3"))+
#                  geom_abline(intercept=eval(parse(text=paste("Z_Model",m,"_Cas",c,"_",q,"QTLs",sep="")))[8],slope= 0,size=0.5,aes(colour="Pop3"))+
#                  theme_classic(base_size = 20, base_family = ""))
#         
#         #Ecriture du pdf
#         pdf(file = paste("Mod",m,"Fec",f,"Cas",c,"QTLs",q,"Qperpop.pdf",sep=""), height = 12, width = 16, paper = "special")
#         print(eval(parse(text=paste("Model",m,"_Fecundity",f,"_Case",c,"_",q,"QTLs","_MeanQperPop",sep =""))))
#         dev.off()  
#         
#       }
#     }
#   }
# }
# ##

## Représentation des Mean Q par environnement, plot Sélection vs Selfing + Fecundity ##
#Création des fichiers de données à traiter par ggplot
for (m in model) {
  for (c in cas) {
    for (q in qtls) {
      assign(paste("Model",m,"_Case",c,"_",q,"QTLs",sep =""),ggplot(results[results$Model==m&
                                                                              results$generation>=1 & results$generation<=200  &
                                                                              results$Cas==c &
                                                                              results$QTLs == q
                                                                            ,]))
      
    }
  }
}

#Création des graphiques moyennes de Q par environnement
for (m in model) {
  for (c in cas) {
    for (q in qtls) {
      #Création du graphe Mean Q per environment
      assign(paste("Model",m,"_Case",c,"_",q,"QTLs","_MeanQperEnv",sep =""), eval(parse(text=paste("Model",m,"_Case",c,"_",q,"QTLs",sep=""))) + aes(generation) +
               labs(y="", title = paste("Model ",m,"/Case ",c,"/",q," QTLs - mean Q per environment",sep=""))+
               facet_grid(Selfing+Fecundity~Selection, labeller = label_both) +
               ylim(eval(parse(text=paste("Z_Model",m,"_Cas",c,"_",q,"QTLs",sep="")))[1]-1,-eval(parse(text=paste("Z_Model",m,"_Cas",c,"_",q,"QTLs",sep="")))[1]+1) +
               geom_point(aes(y=((adlt.q1.p1*adlt.fem.p1+ adlt.q1.p2*adlt.fem.p2+adlt.q1.p3*adlt.fem.p3)/(adlt.fem.p1+adlt.fem.p2+adlt.fem.p3)),color="Pop1-Env1"))+
               geom_point(aes(y=((adlt.q1.p4*adlt.fem.p4+ adlt.q1.p5*adlt.fem.p5+adlt.q1.p6*adlt.fem.p6)/(adlt.fem.p4+adlt.fem.p5+adlt.fem.p6)),color="Pop1-Env2"))+
               geom_point(aes(y=((adlt.q1.p7*adlt.fem.p7+ adlt.q1.p8*adlt.fem.p8+adlt.q1.p9*adlt.fem.p9)/(adlt.fem.p7+adlt.fem.p8+adlt.fem.p9)),color="Pop2-Env3"))+
               geom_point(aes(y=((adlt.q1.p10*adlt.fem.p10+ adlt.q1.p11*adlt.fem.p11+adlt.q1.p12*adlt.fem.p12)/(adlt.fem.p10+adlt.fem.p11+adlt.fem.p12)),color="Pop2-Env4"))+
               geom_point(aes(y=((adlt.q1.p13*adlt.fem.p13+ adlt.q1.p14*adlt.fem.p14+adlt.q1.p15*adlt.fem.p15)/(adlt.fem.p13+adlt.fem.p14+adlt.fem.p15)),color="Pop3-Env5"))+
               geom_point(aes(y=((adlt.q1.p16*adlt.fem.p16+ adlt.q1.p17*adlt.fem.p17+adlt.q1.p18*adlt.fem.p18)/(adlt.fem.p16+adlt.fem.p17+adlt.fem.p18)),color="Pop3-Env6"))+
               geom_point(aes(y=((adlt.q1.p19*adlt.fem.p19+ adlt.q1.p20*adlt.fem.p20+adlt.q1.p21*adlt.fem.p21)/(adlt.fem.p19+adlt.fem.p20+adlt.fem.p21)),color="Pop4-Env7"))+
               geom_point(aes(y=((adlt.q1.p22*adlt.fem.p22+ adlt.q1.p23*adlt.fem.p23+adlt.q1.p24*adlt.fem.p24)/(adlt.fem.p22+adlt.fem.p23+adlt.fem.p24)),color="Pop4-Env8"))+
               geom_abline(intercept=eval(parse(text=paste("Z_Model",m,"_Cas",c,"_",q,"QTLs",sep="")))[1],slope= 0,size=0.5,aes(colour="Pop1-Env1"))+
               geom_abline(intercept=eval(parse(text=paste("Z_Model",m,"_Cas",c,"_",q,"QTLs",sep="")))[2],slope= 0,size=0.5,aes(colour="Pop1-Env2"))+
               geom_abline(intercept=eval(parse(text=paste("Z_Model",m,"_Cas",c,"_",q,"QTLs",sep="")))[3],slope= 0,size=0.5,aes(colour="Pop2-Env3"))+
               geom_abline(intercept=eval(parse(text=paste("Z_Model",m,"_Cas",c,"_",q,"QTLs",sep="")))[4],slope= 0,size=0.5,aes(colour="Pop2-Env4"))+
               geom_abline(intercept=eval(parse(text=paste("Z_Model",m,"_Cas",c,"_",q,"QTLs",sep="")))[5],slope= 0,size=0.5,aes(colour="Pop4-Env8"))+
               geom_abline(intercept=eval(parse(text=paste("Z_Model",m,"_Cas",c,"_",q,"QTLs",sep="")))[6],slope= 0,size=0.5,aes(colour="Pop4-Env7"))+
               geom_abline(intercept=eval(parse(text=paste("Z_Model",m,"_Cas",c,"_",q,"QTLs",sep="")))[7],slope= 0,size=0.5,aes(colour="Pop3-Env6"))+
               geom_abline(intercept=eval(parse(text=paste("Z_Model",m,"_Cas",c,"_",q,"QTLs",sep="")))[8],slope= 0,size=0.5,aes(colour="Pop3-Env5"))+
               theme_classic(base_size = 20, base_family = ""))
      
      #Ecriture du pdf
      pdf(file = paste("Mod",m,"Cas",c,"QTLs",q,"Qperenv.pdf",sep=""), height = 12, width = 16, paper = "special")
      print(eval(parse(text=paste("Model",m,"_Case",c,"_",q,"QTLs","_MeanQperEnv",sep =""))))
      dev.off()  
      
    }
  }
}

#Création des graphiques moyennes de Q par population
for (m in model) {
    for (c in cas) {
      for (q in qtls) {
        #Création du graphe Mean Q per population
        assign(paste("Model",m,"_Case",c,"_",q,"QTLs","_MeanQperPop",sep =""), eval(parse(text=paste("Model",m,"_Case",c,"_",q,"QTLs",sep=""))) + aes(generation) +
                 labs(y="", title = paste("Model ",m,"/Case ",c,"/",q," QTLs - mean Q per population",sep=""))+
                 facet_grid(Selfing+Fecundity~Selection, labeller = label_both) +
                 ylim(eval(parse(text=paste("Z_Model",m,"_Cas",c,"_",q,"QTLs",sep="")))[1]-1,-eval(parse(text=paste("Z_Model",m,"_Cas",c,"_",q,"QTLs",sep="")))[1]+1) +
                 geom_point(aes(y=((adlt.q1.p1*adlt.fem.p1+ adlt.q1.p2*adlt.fem.p2+adlt.q1.p3*adlt.fem.p3+adlt.q1.p4*adlt.fem.p4+ adlt.q1.p5*adlt.fem.p5+adlt.q1.p6*adlt.fem.p6)/
                                     (adlt.fem.p1+adlt.fem.p2+adlt.fem.p3+adlt.fem.p4+adlt.fem.p5+adlt.fem.p6)),color="Pop1"))+
                 geom_point(aes(y=((adlt.q1.p7*adlt.fem.p7+ adlt.q1.p8*adlt.fem.p8+adlt.q1.p9*adlt.fem.p9+adlt.q1.p10*adlt.fem.p10+ adlt.q1.p11*adlt.fem.p11+adlt.q1.p12*adlt.fem.p12)/
                                     (adlt.fem.p7+adlt.fem.p8+adlt.fem.p9+adlt.fem.p10+adlt.fem.p11+adlt.fem.p12)),color="Pop2"))+
                 geom_point(aes(y=((adlt.q1.p13*adlt.fem.p13+ adlt.q1.p14*adlt.fem.p14+adlt.q1.p15*adlt.fem.p15+adlt.q1.p16*adlt.fem.p16+ adlt.q1.p17*adlt.fem.p17+adlt.q1.p18*adlt.fem.p18)/
                                     (adlt.fem.p13+adlt.fem.p14+adlt.fem.p15+adlt.fem.p16+adlt.fem.p17+adlt.fem.p18)),color="Pop3"))+
                 geom_point(aes(y=((adlt.q1.p19*adlt.fem.p19+ adlt.q1.p20*adlt.fem.p20+adlt.q1.p21*adlt.fem.p21+adlt.q1.p22*adlt.fem.p22+ adlt.q1.p23*adlt.fem.p23+adlt.q1.p24*adlt.fem.p24)/
                                     (adlt.fem.p19+adlt.fem.p20+adlt.fem.p21+adlt.fem.p22+adlt.fem.p23+adlt.fem.p24)),color="Pop4"))+
                 geom_abline(intercept=eval(parse(text=paste("Z_Model",m,"_Cas",c,"_",q,"QTLs",sep="")))[1],slope= 0,size=0.5,aes(colour="Pop1"))+
                 geom_abline(intercept=eval(parse(text=paste("Z_Model",m,"_Cas",c,"_",q,"QTLs",sep="")))[2],slope= 0,size=0.5,aes(colour="Pop1"))+
                 geom_abline(intercept=eval(parse(text=paste("Z_Model",m,"_Cas",c,"_",q,"QTLs",sep="")))[3],slope= 0,size=0.5,aes(colour="Pop2"))+
                 geom_abline(intercept=eval(parse(text=paste("Z_Model",m,"_Cas",c,"_",q,"QTLs",sep="")))[4],slope= 0,size=0.5,aes(colour="Pop2"))+
                 geom_abline(intercept=eval(parse(text=paste("Z_Model",m,"_Cas",c,"_",q,"QTLs",sep="")))[5],slope= 0,size=0.5,aes(colour="Pop4"))+
                 geom_abline(intercept=eval(parse(text=paste("Z_Model",m,"_Cas",c,"_",q,"QTLs",sep="")))[6],slope= 0,size=0.5,aes(colour="Pop4"))+
                 geom_abline(intercept=eval(parse(text=paste("Z_Model",m,"_Cas",c,"_",q,"QTLs",sep="")))[7],slope= 0,size=0.5,aes(colour="Pop3"))+
                 geom_abline(intercept=eval(parse(text=paste("Z_Model",m,"_Cas",c,"_",q,"QTLs",sep="")))[8],slope= 0,size=0.5,aes(colour="Pop3"))+
                 theme_classic(base_size = 20, base_family = ""))
        
        #Ecriture du pdf
        pdf(file = paste("Mod",m,"Cas",c,"QTLs",q,"Qperpop.pdf",sep=""), height = 12, width = 16, paper = "special")
        print(eval(parse(text=paste("Model",m,"_Case",c,"_",q,"QTLs","_MeanQperPop",sep =""))))
        dev.off()  
        
    }
  }
}
##

