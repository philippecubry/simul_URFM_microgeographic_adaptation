setwd("~/nemo_cluster/")
library("ggplot2")
QSTsBiallelic<-read.table("Sims2/results_QSTs.txt",dec=",")
FSTqsBiallelic<-read.table("Sims2/FstTableV2.txt")
QSTsConti<-read.table("Sims2Conti/results_QSTs_conti.txt")
FSTqsConti<-read.table("Sims2Conti/FstTableContiV2.txt")
Biallelic <- merge (QSTsBiallelic,FSTqsBiallelic)
Continuous <- merge (QSTsConti,FSTqsConti)


qplot(Selection,Fis,data=Biallelic,colour=factor(Selfing),fill=Model,alpha=0.5,geom="boxplot", outlier.colour="green",main = "Fis vs Selection, Biallelic mutation model")
qplot(Selection,Fis,data=Continuous,colour=factor(Selfing),fill=Model,alpha=0.5,geom="boxplot", outlier.colour="green",main = "Fis vs Selection, Continuous mutation model")

qplot(Selection,FstPop,data=Biallelic,facets = Cas ~ Model + Selfing,colour=factor(QTLs),geom="jitter", alpha=0.7,outlier.colour="green",main = "FSTqPop vs Selection, Biallelic mutation model") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
qplot(Selection,QSTpop,data=Biallelic,facets = Cas ~ Model + Selfing,colour=factor(QTLs),geom="jitter", alpha=0.7,outlier.colour="green",main = "QSTPop vs Selection, Biallelic mutation model") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
qplot(Selection,FstEnv,data=Biallelic,facets = Cas ~ Model + Selfing,colour=factor(QTLs),geom="jitter", alpha=0.7,outlier.colour="green",main = "FSTqEnv vs Selection, Biallelic mutation model") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
qplot(Selection,QSTenv,data=Biallelic,facets = Cas ~ Model + Selfing,colour=factor(QTLs),geom="jitter", alpha=0.7,outlier.colour="green",main = "QSTEnv vs Selection, Biallelic mutation model") + theme(axis.text.x = element_text(angle = 90, hjust = 1))

qplot(Selection,FstPop,data=Continuous,facets = Cas ~ Model + Selfing,colour=factor(QTLs),geom="jitter", alpha=0.7,outlier.colour="green",main = "FSTqPop vs Selection, Continuous mutation model") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
qplot(Selection,QSTpop,data=Continuous,facets = Cas ~ Model + Selfing,colour=factor(QTLs),geom="jitter", alpha=0.7,outlier.colour="green",main = "QSTPop vs Selection, Continuous mutation model") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
qplot(Selection,FstEnv,data=Continuous,facets = Cas ~ Model + Selfing,colour=factor(QTLs),geom="jitter", alpha=0.7,outlier.colour="green",main = "FSTqEnv vs Selection, Continuous mutation model") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
qplot(Selection,QSTenv,data=Continuous,facets = Cas ~ Model + Selfing,colour=factor(QTLs),geom="jitter", alpha=0.7,outlier.colour="green",main = "QSTEnv vs Selection, Continuous mutation model") + theme(axis.text.x = element_text(angle = 90, hjust = 1))


qplot(Selection,FstPop,data=Biallelic,facets = . ~ QTLs,colour=factor(Selfing),geom="boxplot", outlier.colour="green",main = "FSTqPop vs Selection, Biallelic mutation model")
qplot(Selection,QSTpop,data=Biallelic,facets = . ~ QTLs,colour=factor(Selfing),geom="boxplot", outlier.colour="green",main = "QSTpop vs Selection, Biallelic mutation model")
qplot(Selection,FstEnv,data=Biallelic,facets = . ~ QTLs,colour=factor(Selfing),geom="boxplot", outlier.colour="green",main = "FSTqEnv vs Selection, Biallelic mutation model")
qplot(Selection,QSTenv,data=Biallelic,facets = . ~ QTLs,colour=factor(Selfing),geom="boxplot", outlier.colour="green",main = "QSTenv vs Selection, Biallelic mutation model")

qplot(Selection,FstPop,data=Continuous,facets = . ~ QTLs,colour=factor(Selfing),geom="boxplot", outlier.colour="green",main = "FSTqPop vs Selection, Continuous mutation model")
qplot(Selection,QSTpop,data=Continuous,facets = . ~ QTLs,colour=factor(Selfing),geom="boxplot", outlier.colour="green",main = "QSTpop vs Selection, Continuous mutation model")
qplot(Selection,FstEnv,data=Continuous,facets = . ~ QTLs,colour=factor(Selfing),geom="boxplot", outlier.colour="green",main = "FSTqEnv vs Selection, Continuous mutation model")
qplot(Selection,QSTenv,data=Continuous,facets = . ~ QTLs,colour=factor(Selfing),geom="boxplot", outlier.colour="green",main = "QSTenv vs Selection, Continuous mutation model")

qplot(FstEnv,QSTenv,data=Continuous,facets = Model ~ Selfing,colour=factor(QTLs),geom="point", outlier.colour="green",main = "QSTenv vs FstEnv, Continuous mutation model")
qplot(FstPop,QSTpop,data=Continuous,facets = Model ~ Selfing,colour=factor(QTLs),geom="point", outlier.colour="green",main = "QSTpop vs FstPop, Continuous mutation model")

qplot(FstEnv,QSTenv,data=Biallelic,facets = Model ~ Selfing,colour=factor(QTLs),geom="point", outlier.colour="green",main = "QSTenv vs FstEnv, Biallelic mutation model")
qplot(FstPop,QSTpop,data=Biallelic,facets = Model ~ Selfing,colour=factor(QTLs),geom="point", outlier.colour="green",main = "QSTpop vs FstPop, Biallelic mutation model")
