#définit le répertoire de travail
setwd("~/nemo_cluster/Sims2Conti")

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
            temp <- read.delim(file=(paste("Sim2_conti_Model",m,"_",s,"self_",q,"QTLs_fitperpop","/stats/Sim2_conti_Model",m,"_",s,"self_",q,"QTLs_Cas",c,"_SelInt",w,"_fec",f,".txt",sep = "")),header=T)
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
Z_Model2_Cas1_10QTLs <- c(-13.49,-4.49666666666667,-7.4944444444,1.49888888888889,13.49,4.49666666666667,7.4944444444,-1.49888888888889)
Z_Model2_Cas2_10QTLs <- c(-16.8625,-5.62083333333333,-9.36805555555556,1.87361111111111,16.8625,5.62083333333333,9.36805555555556,-1.87361111111111)
Z_Model2_Cas3_10QTLs <- c(-6.745,-2.24833333333333,-3.74722222222222,0.749444444444444,6.745,2.24833333333333,3.74722222222222,-0.749444444444444)
Z_Model2_Cas4_10QTLs <- c(-3.3725,-1.1241666667,-1.87361111111111,0.374722222222222,3.3725,1.1241666667,1.87361111111111,-0.374722222222222)
Z_Model2_Cas1_50QTLs <- c(-67.45,-22.4833333333333,-37.4722222222222,7.49444444444444,67.45,22.4833333333333,37.4722222222222,-7.49444444444444)
Z_Model2_Cas2_50QTLs <- c(-84.3125,-28.1041666666667,-46.8402777777778,9.36805555555556,84.3125,28.1041666666667,46.8402777777778,-9.36805555555556)
Z_Model2_Cas3_50QTLs <- c(-33.725,-11.2416666666667,-18.7361111111111,3.74722222222222,33.725,11.2416666666667,18.7361111111111,-3.74722222222222)
Z_Model2_Cas4_50QTLs <- c(-16.8625,-5.62083333333333,-9.36805555555556,1.87361111111111,16.8625,5.62083333333333,9.36805555555556,-1.87361111111111)

Z_Model3_Cas1_10QTLs <- c(-13.49,-9.63571428571429,-5.78142857142857,-1.92714285714286,13.49,9.63571428571429,5.78142857142857,1.92714285714286)
Z_Model3_Cas2_10QTLs <- c(-16.8625,-12.0446428571429,-7.22678571428572,-2.40892857142857,16.8625,12.0446428571429,7.22678571428572,2.40892857142857)
Z_Model3_Cas3_10QTLs <- c(-6.745,-4.81785714285714,-2.89071428571429,-0.963571428571429,6.745,4.81785714285714,2.89071428571429,0.963571428571429)
Z_Model3_Cas4_10QTLs <- c(-3.3725,-2.40892857142857,-1.44535714285714,-0.481785714285714,3.3725,2.40892857142857,1.44535714285714,0.481785714285714)
Z_Model3_Cas1_50QTLs <- c(-67.45,-48.1785714285714,-28.9071428571429,-9.63571428571429,67.45,48.1785714285714,28.9071428571429,9.63571428571429)
Z_Model3_Cas2_50QTLs <- c(-84.3125,-60.2232142857143,-36.1339285714286,-12.0446428571429,84.3125,60.2232142857143,36.1339285714286,12.0446428571429)
Z_Model3_Cas3_50QTLs <- c(-33.725,-24.0892857142857,-14.4535714285714,-4.81785714285714,33.725,24.0892857142857,14.4535714285714,4.81785714285714)
Z_Model3_Cas4_50QTLs <- c(-16.8625,-12.0446428571429,-7.22678571428572,-2.40892857142857,16.8625,12.0446428571429,7.22678571428572,2.40892857142857)

#On convertit les NaN en 0 dans les valeurs phénotypiques réalisées
results[,c("adlt.q1.p1","adlt.q1.p2","adlt.q1.p3","adlt.q1.p4","adlt.q1.p5","adlt.q1.p6","adlt.q1.p7","adlt.q1.p8",
           "adlt.q1.p9","adlt.q1.p10","adlt.q1.p11","adlt.q1.p12","adlt.q1.p13","adlt.q1.p14","adlt.q1.p15","adlt.q1.p16",
           "adlt.q1.p17","adlt.q1.p18","adlt.q1.p19","adlt.q1.p20","adlt.q1.p21","adlt.q1.p22","adlt.q1.p23","adlt.q1.p24",
           "adlt.Va.q1.p1","adlt.Va.q1.p2","adlt.Va.q1.p3","adlt.Va.q1.p4","adlt.Va.q1.p5","adlt.Va.q1.p6","adlt.Va.q1.p7","adlt.Va.q1.p8",
           "adlt.Va.q1.p9","adlt.Va.q1.p10","adlt.Va.q1.p11","adlt.Va.q1.p12","adlt.Va.q1.p13","adlt.Va.q1.p14","adlt.Va.q1.p15","adlt.Va.q1.p16",
           "adlt.Va.q1.p17","adlt.Va.q1.p18","adlt.Va.q1.p19","adlt.Va.q1.p20","adlt.Va.q1.p21","adlt.Va.q1.p22","adlt.Va.q1.p23","adlt.Va.q1.p24")]<-rapply( results[,c("adlt.q1.p1","adlt.q1.p2","adlt.q1.p3","adlt.q1.p4","adlt.q1.p5","adlt.q1.p6","adlt.q1.p7","adlt.q1.p8",
                   "adlt.q1.p9","adlt.q1.p10","adlt.q1.p11","adlt.q1.p12","adlt.q1.p13","adlt.q1.p14","adlt.q1.p15","adlt.q1.p16",
                   "adlt.q1.p17","adlt.q1.p18","adlt.q1.p19","adlt.q1.p20","adlt.q1.p21","adlt.q1.p22","adlt.q1.p23","adlt.q1.p24","adlt.Va.q1.p1","adlt.Va.q1.p2","adlt.Va.q1.p3","adlt.Va.q1.p4","adlt.Va.q1.p5","adlt.Va.q1.p6","adlt.Va.q1.p7","adlt.Va.q1.p8",
                   "adlt.Va.q1.p9","adlt.Va.q1.p10","adlt.Va.q1.p11","adlt.Va.q1.p12","adlt.Va.q1.p13","adlt.Va.q1.p14","adlt.Va.q1.p15","adlt.Va.q1.p16",
                   "adlt.Va.q1.p17","adlt.Va.q1.p18","adlt.Va.q1.p19","adlt.Va.q1.p20","adlt.Va.q1.p21","adlt.Va.q1.p22","adlt.Va.q1.p23","adlt.Va.q1.p24")], f=function(x) ifelse(is.nan(x),0,x), how="replace" )

## Représentation des Mean Q par environnement, plot Sélection vs Selfing + Fecundity ##
#Création des fichiers de données à traiter par ggplot
for (m in model) {
  for (c in cas) {
    for (q in qtls) {
      assign(paste("Model",m,"conti_Case",c,"_",q,"QTLs",sep =""),ggplot(results[results$Model==m&
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
      assign(paste("Model",m,"conti_Case",c,"_",q,"QTLs","_MeanQperEnv",sep =""), eval(parse(text=paste("Model",m,"conti_Case",c,"_",q,"QTLs",sep=""))) + aes(generation) +
               labs(y="", title = paste("Model ",m,"conti /Case ",c,"/",q," QTLs - mean Q per environment",sep=""))+
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
      pdf(file = paste("Mod",m,"contiCas",c,"QTLs",q,"Qperenv.pdf",sep=""), height = 12, width = 16, paper = "special")
      print(eval(parse(text=paste("Model",m,"conti_Case",c,"_",q,"QTLs","_MeanQperEnv",sep =""))))
      dev.off()  
      
    }
  }
}

#Création des graphiques moyennes de Q par population
for (m in model) {
  for (c in cas) {
    for (q in qtls) {
      #Création du graphe Mean Q per population
      assign(paste("Model",m,"conti_Case",c,"_",q,"QTLs","_MeanQperPop",sep =""), eval(parse(text=paste("Model",m,"conti_Case",c,"_",q,"QTLs",sep=""))) + aes(generation) +
               labs(y="", title = paste("Model ",m,"conti /Case ",c,"/",q," QTLs - mean Q per population",sep=""))+
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
      pdf(file = paste("Mod",m,"contiCas",c,"QTLs",q,"Qperpop.pdf",sep=""), height = 12, width = 16, paper = "special")
      print(eval(parse(text=paste("Model",m,"conti_Case",c,"_",q,"QTLs","_MeanQperPop",sep =""))))
      dev.off()  
      
    }
  }
}
##