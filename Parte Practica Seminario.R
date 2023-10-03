#### Replicación del código seminario #########
## Algunas librerias necesarias #######
library(zoo)
library(readxl)
library(ggplot2)
library(reshape2)
library(stats)
library(xtable)
library(summarytools)
library(markovchain)
library(sqldf)
#-----------------------------------------------------
########### Carga de las bases de datos ##############
setwd("C:/Users/andre/OneDrive/Escritorio/Andrés/2023-2/Seminario/Bases de datos expo 2")

barrie<-read.table("Barrie.txt", header = TRUE, sep = "\t")
bellevile<-read.table("Bellevile.txt", header = TRUE, sep = "\t")
brantford<-read.table("Brantford.txt", header = TRUE, sep = "\t")
burlington<-read.table("Burlington.txt", header = TRUE, sep = "\t")
chatam<-read.table("Chatham.txt", header = TRUE, sep = "\t")
cornwall<-read.table("Cornwall.txt", header = TRUE, sep = "\t")
dorset<-read.table("Dorset.txt", header = TRUE, sep = "\t")
grandbend<-read.table("GrandBend.txt", header = TRUE, sep = "\t")
guelph<-read.table("Guelph.txt", header = TRUE, sep = "\t")
hamiltonD<-read.table("HamiltonD.txt", header = TRUE, sep = "\t")
hamiltonW<-read.table("HamiltonW.txt", header = TRUE, sep = "\t")
kingston<-read.table("Kingston.txt", header = TRUE, sep = "\t")
kitchner<-read.table("Kitchner.txt", header = TRUE, sep = "\t")
london<-read.table("London.txt", header = TRUE, sep = "\t")
mississauga<-read.table("Mississauga.txt", header = TRUE, sep = "\t")
newmarket<-read.table("Newmarket.txt", header = TRUE, sep = "\t")
northbay<-read.table("Northbay.txt", header = TRUE, sep = "\t")
oakville<-read.table("Oakville.txt", header = TRUE, sep = "\t")
ottawaD<-read.table("OttawaD.txt", header = TRUE, sep = "\t")
parrysound<-read.table("Parrysound.txt", header = TRUE, sep = "\t")
peterbourgh<-read.table("Peterbourgh.txt", header = TRUE, sep = "\t")
portstanley<-read.table("Portstanley.txt", header = TRUE, sep = "\t") 
sarnia<-read.table("Sarnia.txt", header = TRUE, sep = "\t")
stchaterines<-read.table("Stchatherines.txt", header = TRUE, sep = "\t")
sudbury<-read.table("Sudbury.txt", header = TRUE, sep = "\t")
thunderbay<-read.table("Thunderbay.txt", header = TRUE, sep = "\t")
tiverton<-read.table("Tiverton.txt", header = TRUE, sep = "\t")
torontoE<-read.table("TorontoE.txt", header = TRUE, sep = "\t")
torontoW<-read.table("TorontoW.txt", header = TRUE, sep = "\t")
windsorD<-read.table("WindsorD.txt", header = TRUE, sep = "\t")
windsorW<-read.table("WindsorW.txt", header = TRUE, sep = "\t")
#-----------------------------------------------------
########### Barrie ##################
barrie$estado <- cut(barrie$X2, 
                              breaks = c(-Inf, 3.001, 6.0001), 
                              labels = c("1", "2"))
# Crear la cadena de Markov 
cadena_barrie <- markovchainFit(data = barrie$estado)
## Matriz de transición 
cadena_barrie$estimate
## Distribucion estacionaria
Pbarrie <- matrix(c(0.9436759, 0.05632411,
                     0.7037037, 0.29629630),byrow = T, nrow = 2)
xtable(Pbarrie)
markovbarrie <- new("markovchain", states = c("1", "2"), transitionMatrix =
                        Pbarrie,name = "A markovchain Object" )
disesta_barrie <- steadyStates(markovbarrie);disesta_barrie

##Conteo y chi cuadrado.
matriz_nijbarrie <- as.matrix(table(barrie$estado[-length(barrie$estado)], barrie$estado[-1]))
chi2_barrie <- chisq.test(matriz_nijbarrie, p = Pbarrie^n)
##Tiempo medio de persistencia.
timeper_barrie<-round(c(1/(1-cadena_barrie$estimate[1,1]),1/(1-cadena_barrie$estimate[2,2])),2)
#-----------------------------------------------------
########### Bellevile ##################
bellevile$estado <- cut(bellevile$X2, 
                     breaks = c(-Inf, 3.001, 6.0001, Inf), 
                     labels = c("1", "2", "3"))
# Crear la cadena de Markov 
cadena_bellevile <- markovchainFit(data = bellevile$estado)
## Matriz de transición 
cadena_bellevile$estimate
## Distribucion estacionaria
Pbellevile <- matrix(c(cadena_bellevile$estimate[1,1], cadena_bellevile$estimate[1,2],cadena_bellevile$estimate[1,3],
                       cadena_bellevile$estimate[2,1], cadena_bellevile$estimate[2,2],cadena_bellevile$estimate[2,3],
                       cadena_bellevile$estimate[3,1], cadena_bellevile$estimate[3,2],cadena_bellevile$estimate[3,3]),
                     byrow = T, nrow = 3)
xtable(Pbellevile)
markovbellevile <- new("markovchain", states = c("1", "2","3"), transitionMatrix =
                         Pbellevile,name = "A markovchain Object" )
disesta_bellevile <- steadyStates(markovbellevile);disesta_bellevile
##Conteo y chi cuadrado.
matriz_nijbellevile <- as.matrix(table(bellevile$estado[-length(bellevile$estado)], bellevile$estado[-1]))
chi2_bellevile <- chisq.test(matriz_nijbellevile, p = Pbellevile^n)
##Tiempo medio de persistencia.
timeper_bellevile<-round(c(1/(1-cadena_bellevile$estimate[1,1]),1/(1-cadena_bellevile$estimate[2,2]),1/(1-cadena_bellevile$estimate[3,3])),2)
#-----------------------------------------------------
########### Brantford ##################
brantford$estado <- cut(brantford$X3, 
                        breaks = c(-Inf, 3.001, 6.0001, Inf), 
                        labels = c("1", "2", "3"))
# Crear la cadena de Markov 
cadena_brantford <- markovchainFit(data = brantford$estado)
## Matriz de transición 
cadena_brantford$estimate
## Distribucion estacionaria
Pbrantford <- matrix(c(cadena_brantford$estimate[1,1], cadena_brantford$estimate[1,2],cadena_brantford$estimate[1,3],
                       cadena_brantford$estimate[2,1], cadena_brantford$estimate[2,2],cadena_brantford$estimate[2,3],
                       cadena_brantford$estimate[3,1], cadena_brantford$estimate[3,2],cadena_brantford$estimate[3,3]),
                     byrow = T, nrow = 3)
xtable(Pbrantford)
markovbrantford <- new("markovchain", states = c("1", "2","3"), transitionMatrix =
                         Pbrantford,name = "A markovchain Object" )
disesta_brantford <- steadyStates(markovbrantford);disesta_brantford
##Conteo y chi cuadrado.
matriz_nijbrantford <- as.matrix(table(brantford$estado[-length(brantford$estado)], brantford$estado[-1]))
chi2_brantford <- chisq.test(matriz_nijbrantford, p = Pbrantford^n)
##Tiempo medio de persistencia.
timeper_brantford<-round(c(1/(1-cadena_brantford$estimate[1,1]),1/(1-cadena_brantford$estimate[2,2]),1/(1-cadena_brantford$estimate[3,3])),2)
#-----------------------------------------------------
########### Burlington ##################
burlington$estado <- cut(burlington$X2, 
                        breaks = c(-Inf, 3.001, 6.0001, Inf), 
                        labels = c("1", "2", "3"))
# Crear la cadena de Markov 
cadena_burlington <- markovchainFit(data = burlington$estado)
## Matriz de transición 
cadena_burlington$estimate
## Distribucion estacionaria
Pburlington <- matrix(c(cadena_burlington$estimate[1,1], cadena_burlington$estimate[1,2],cadena_burlington$estimate[1,3],
                       cadena_burlington$estimate[2,1], cadena_burlington$estimate[2,2],cadena_burlington$estimate[2,3],
                       cadena_burlington$estimate[3,1], cadena_burlington$estimate[3,2],cadena_burlington$estimate[3,3]),
                     byrow = T, nrow = 3)
xtable(Pburlington)
markovburlington <- new("markovchain", states = c("1", "2","3"), transitionMatrix =
                         Pburlington,name = "A markovchain Object" )
disesta_burlington <- steadyStates(markovburlington);disesta_burlington
##Conteo y chi cuadrado.
matriz_nijburlington <- as.matrix(table(burlington$estado[-length(burlington$estado)], burlington$estado[-1]))
chi2_burlington <- chisq.test(matriz_nijburlington, p = Pburlington^n)
##Tiempo medio de persistencia.
timeper_burlington<-round(c(1/(1-cadena_burlington$estimate[1,1]),1/(1-cadena_burlington$estimate[2,2]),1/(1-cadena_burlington$estimate[3,3])),2)
#-----------------------------------------------------
########### chatham ##################
chatam$estado <- cut(chatam$X2, 
                         breaks = c(-Inf, 3.001, 6.0001, Inf), 
                         labels = c("1", "2", "3"))
# Crear la cadena de Markov 
cadena_chatam <- markovchainFit(data = chatam$estado)
## Matriz de transición 
cadena_chatam$estimate
## Distribucion estacionaria
Pchatam <- matrix(c(cadena_chatam$estimate[1,1], cadena_chatam$estimate[1,2],cadena_chatam$estimate[1,3],
                        cadena_chatam$estimate[2,1], cadena_chatam$estimate[2,2],cadena_chatam$estimate[2,3],
                        cadena_chatam$estimate[3,1], cadena_chatam$estimate[3,2],cadena_chatam$estimate[3,3]),
                      byrow = T, nrow = 3)
xtable(Pchatam)
markovchatam <- new("markovchain", states = c("1", "2","3"), transitionMatrix =
                          Pchatam,name = "A markovchain Object" )
disesta_chatam <- steadyStates(markovchatam);disesta_chatam
##Conteo y chi cuadrado.
matriz_nijchatam <- as.matrix(table(chatam$estado[-length(chatam$estado)], chatam$estado[-1]))
chi2_chatam <- chisq.test(matriz_nijchatam, p = Pchatam^n)
##Tiempo medio de persistencia.
timeper_chatam<-round(c(1/(1-cadena_chatam$estimate[1,1]),1/(1-cadena_chatam$estimate[2,2]),1/(1-cadena_chatam$estimate[3,3])),2)
#-----------------------------------------------------
########### Cornwall ##################
cornwall$estado <- cut(cornwall$X2, 
                     breaks = c(-Inf, 3.001, 6.0001, Inf), 
                     labels = c("1", "2", "3"))
# Crear la cadena de Markov 
cadena_cornwall <- markovchainFit(data = chatam$estado)
## Matriz de transición 
cadena_cornwall$estimate
## Distribucion estacionaria
Pcornwall <- matrix(c(cadena_cornwall$estimate[1,1], cadena_cornwall$estimate[1,2],cadena_cornwall$estimate[1,3],
                    cadena_cornwall$estimate[2,1], cadena_cornwall$estimate[2,2],cadena_cornwall$estimate[2,3],
                    cadena_cornwall$estimate[3,1], cadena_cornwall$estimate[3,2],cadena_cornwall$estimate[3,3]),
                  byrow = T, nrow = 3)
xtable(Pcornwall)
markovcornwall <- new("markovchain", states = c("1", "2","3"), transitionMatrix =
                      Pcornwall,name = "A markovchain Object" )
disesta_cornwall <- steadyStates(markovcornwall);disesta_cornwall
##Conteo y chi cuadrado.
matriz_nijcornwall <- as.matrix(table(cornwall$estado[-length(cornwall$estado)], cornwall$estado[-1]))
chi2_cornwall <- chisq.test(matriz_nijcornwall, p = Pcornwall^n)
##Tiempo medio de persistencia.
timeper_cornwall<-round(c(1/(1-cadena_cornwall$estimate[1,1]),1/(1-cadena_cornwall$estimate[2,2]),1/(1-cadena_cornwall$estimate[3,3])),2)
#-----------------------------------------------------
########### Dorset ##################
dorset$estado <- cut(dorset$X2, 
                       breaks = c(-Inf, 3.001, 6.0001), 
                       labels = c("1", "2"))
# Crear la cadena de Markov 
cadena_dorset <- markovchainFit(data = dorset$estado)
## Matriz de transición 
cadena_dorset$estimate
## Distribucion estacionaria
Pdorset <- matrix(c(cadena_dorset$estimate[1,1], cadena_dorset$estimate[1,2],
                      cadena_dorset$estimate[2,1], cadena_dorset$estimate[2,2]),
                    byrow = T, nrow = 2)
xtable(Pdorset)
markovdorset <- new("markovchain", states = c("1", "2"), transitionMatrix =
                        Pdorset,name = "A markovchain Object" )
disesta_dorset <- steadyStates(markovdorset);disesta_dorset
##Conteo y chi cuadrado.
matriz_nijdorset <- as.matrix(table(dorset$estado[-length(dorset$estado)], dorset$estado[-1]))
chi2_dorset <- chisq.test(matriz_nijdorset, p = Pdorset^n)
##Tiempo medio de persistencia.
timeper_dorset<-round(c(1/(1-cadena_dorset$estimate[1,1]),1/(1-cadena_dorset$estimate[2,2])),2)
#-----------------------------------------------------
########### Grandbend ##################
grandbend$estado <- cut(grandbend$X2, 
                       breaks = c(-Inf, 3.001, 6.0001, Inf), 
                       labels = c("1", "2", "3"))
# Crear la cadena de Markov 
cadena_grandbend <- markovchainFit(data = grandbend$estado)
## Matriz de transición 
cadena_grandbend$estimate
## Distribucion estacionaria
Pgrandbend <- matrix(c(cadena_grandbend$estimate[1,1], cadena_grandbend$estimate[1,2],cadena_grandbend$estimate[1,3],
                      cadena_grandbend$estimate[2,1], cadena_grandbend$estimate[2,2],cadena_grandbend$estimate[2,3],
                      cadena_grandbend$estimate[3,1], cadena_grandbend$estimate[3,2],cadena_grandbend$estimate[3,3]),
                    byrow = T, nrow = 3)
xtable(Pgrandbend)
markovgrandbend <- new("markovchain", states = c("1", "2","3"), transitionMatrix =
                        Pgrandbend,name = "A markovchain Object" )
disesta_grandbend <- steadyStates(markovgrandbend);disesta_grandbend
##Conteo y chi cuadrado.
matriz_nijgrandbend <- as.matrix(table(grandbend$estado[-length(grandbend$estado)], grandbend$estado[-1]))
chi2_grandbend <- chisq.test(matriz_nijgrandbend, p = Pgrandbend^n)
##Tiempo medio de persistencia.
timeper_grandbend<-round(c(1/(1-cadena_grandbend$estimate[1,1]),1/(1-cadena_grandbend$estimate[2,2]),1/(1-cadena_grandbend$estimate[3,3])),2)
#-----------------------------------------------------
########### Guelph ##################
guelph$estado <- cut(guelph$X2, 
                        breaks = c(-Inf, 3.001, 6.0001, Inf), 
                        labels = c("1", "2", "3"))
# Crear la cadena de Markov 
cadena_guelph <- markovchainFit(data = guelph$estado)
## Matriz de transición 
cadena_guelph$estimate
## Distribucion estacionaria
Pguelph <- matrix(c(cadena_guelph$estimate[1,1], cadena_guelph$estimate[1,2],cadena_guelph$estimate[1,3],
                       cadena_guelph$estimate[2,1], cadena_guelph$estimate[2,2],cadena_guelph$estimate[2,3],
                       cadena_guelph$estimate[3,1], cadena_guelph$estimate[3,2],cadena_guelph$estimate[3,3]),
                     byrow = T, nrow = 3)
xtable(Pguelph)
markovguelph <- new("markovchain", states = c("1", "2","3"), transitionMatrix =
                         Pguelph,name = "A markovchain Object" )
disesta_guelph <- steadyStates(markovguelph);disesta_guelph
##Conteo y chi cuadrado.
matriz_nijguelph <- as.matrix(table(guelph$estado[-length(guelph$estado)], guelph$estado[-1]))
chi2_guelph <- chisq.test(matriz_nijguelph, p = Pguelph^n)
##Tiempo medio de persistencia.
timeper_guelph<-round(c(1/(1-cadena_guelph$estimate[1,1]),1/(1-cadena_guelph$estimate[2,2]),1/(1-cadena_guelph$estimate[3,3])),2)
#-----------------------------------------------------
########### HamiltonD ##################
hamiltonD$estado <- cut(hamiltonD$X3, 
                     breaks = c(-Inf, 3.001, 6.0001, Inf), 
                     labels = c("1", "2", "3"))
# Crear la cadena de Markov 
cadena_hamiltonD <- markovchainFit(data = hamiltonD$estado)
## Matriz de transición 
cadena_hamiltonD$estimate
## Distribucion estacionaria
PhamiltonD <- matrix(c(cadena_hamiltonD$estimate[1,1], cadena_hamiltonD$estimate[1,2],cadena_hamiltonD$estimate[1,3],
                    cadena_hamiltonD$estimate[2,1], cadena_hamiltonD$estimate[2,2],cadena_hamiltonD$estimate[2,3],
                    cadena_hamiltonD$estimate[3,1], cadena_hamiltonD$estimate[3,2],cadena_hamiltonD$estimate[3,3]),
                  byrow = T, nrow = 3)
xtable(PhamiltonD)
markovhamiltonD <- new("markovchain", states = c("1", "2","3"), transitionMatrix =
                      PhamiltonD,name = "A markovchain Object" )
disesta_hamiltonD <- steadyStates(markovhamiltonD);disesta_hamiltonD
##Conteo y chi cuadrado.
matriz_nijhamiltonD <- as.matrix(table(hamiltonD$estado[-length(hamiltonD$estado)], hamiltonD$estado[-1]))
chi2_hamiltonD <- chisq.test(matriz_nijhamiltonD, p = PhamiltonD^n)
##Tiempo medio de persistencia.
timeper_hamiltonD<-round(c(1/(1-cadena_hamiltonD$estimate[1,1]),1/(1-cadena_hamiltonD$estimate[2,2]),1/(1-cadena_hamiltonD$estimate[3,3])),2)
#-----------------------------------------------------
########### HamiltonW ##################
hamiltonW$estado <- cut(hamiltonW$X2, 
                        breaks = c(-Inf, 3.001, 6.0001, Inf), 
                        labels = c("1", "2", "3"))
# Crear la cadena de Markov 
cadena_hamiltonW <- markovchainFit(data = hamiltonW$estado)
## Matriz de transición 
cadena_hamiltonW$estimate
## Distribucion estacionaria
PhamiltonW <- matrix(c(cadena_hamiltonW$estimate[1,1], cadena_hamiltonW$estimate[1,2],cadena_hamiltonW$estimate[1,3],
                       cadena_hamiltonW$estimate[2,1], cadena_hamiltonW$estimate[2,2],cadena_hamiltonW$estimate[2,3],
                       cadena_hamiltonW$estimate[3,1], cadena_hamiltonW$estimate[3,2],cadena_hamiltonW$estimate[3,3]),
                     byrow = T, nrow = 3)
xtable(PhamiltonW)
markovhamiltonW <- new("markovchain", states = c("1", "2","3"), transitionMatrix =
                         PhamiltonW,name = "A markovchain Object" )
disesta_hamiltonW <- steadyStates(markovhamiltonW);disesta_hamiltonW
##Conteo y chi cuadrado.
matriz_nijhamiltonW <- as.matrix(table(hamiltonW$estado[-length(hamiltonW$estado)], hamiltonW$estado[-1]))
chi2_hamiltonW <- chisq.test(matriz_nijhamiltonW, p = PhamiltonW^n)
##Tiempo medio de persistencia.
timeper_hamiltonW<-round(c(1/(1-cadena_hamiltonW$estimate[1,1]),1/(1-cadena_hamiltonW$estimate[2,2]),1/(1-cadena_hamiltonW$estimate[3,3])),2)
#-----------------------------------------------------
########### Kingston ##################
kingston$estado <- cut(kingston$X2, 
                        breaks = c(-Inf, 3.001, 6.0001, Inf), 
                        labels = c("1", "2", "3"))
# Crear la cadena de Markov 
cadena_kingston <- markovchainFit(data = kingston$estado)
## Matriz de transición 
cadena_kingston$estimate
## Distribucion estacionaria
Pkingston <- matrix(c(cadena_kingston$estimate[1,1], cadena_kingston$estimate[1,2],cadena_kingston$estimate[1,3],
                       cadena_kingston$estimate[2,1], cadena_kingston$estimate[2,2],cadena_kingston$estimate[2,3],
                       cadena_kingston$estimate[3,1], cadena_kingston$estimate[3,2],cadena_kingston$estimate[3,3]),
                     byrow = T, nrow = 3)
xtable(Pkingston)
markovkingston<- new("markovchain", states = c("1", "2","3"), transitionMatrix =
                         Pkingston,name = "A markovchain Object" )
disesta_kingston <- steadyStates(markovkingston);disesta_kingston
##Conteo y chi cuadrado.
matriz_nijkingston <- as.matrix(table(kingston$estado[-length(kingston$estado)], kingston$estado[-1]))
chi2_kingston <- chisq.test(matriz_nijkingston, p = Pkingston^n)
##Tiempo medio de persistencia.
timeper_kingston<-round(c(1/(1-cadena_kingston$estimate[1,1]),1/(1-cadena_kingston$estimate[2,2]),1/(1-cadena_kingston$estimate[3,3])),2)
#-----------------------------------------------------
########### Kitchner ##################
kitchner$estado <- cut(kitchner$X3, 
                       breaks = c(-Inf, 3.001, 6.0001, Inf), 
                       labels = c("1", "2", "3"))
# Crear la cadena de Markov 
cadena_kitchner <- markovchainFit(data = kitchner$estado)
## Matriz de transición 
cadena_kitchner$estimate
## Distribucion estacionaria
Pkitchner <- matrix(c(cadena_kitchner$estimate[1,1], cadena_kitchner$estimate[1,2],cadena_kitchner$estimate[1,3],
                      cadena_kitchner$estimate[2,1], cadena_kitchner$estimate[2,2],cadena_kitchner$estimate[2,3],
                      cadena_kitchner$estimate[3,1], cadena_kitchner$estimate[3,2],cadena_kitchner$estimate[3,3]),
                    byrow = T, nrow = 3)
xtable(Pkitchner)
markovkitchner<- new("markovchain", states = c("1", "2","3"), transitionMatrix =
                       Pkitchner,name = "A markovchain Object" )
disesta_kitchner <- steadyStates(markovkitchner);disesta_kitchner
##Conteo y chi cuadrado.
matriz_nijkitchner <- as.matrix(table(kitchner$estado[-length(kitchner$estado)], kitchner$estado[-1]))
chi2_kitchner <- chisq.test(matriz_nijkitchner, p = Pkitchner^n)
##Tiempo medio de persistencia.
timeper_kitchner<-round(c(1/(1-cadena_kitchner$estimate[1,1]),1/(1-cadena_kitchner$estimate[2,2]),1/(1-cadena_kitchner$estimate[3,3])),2)
#-----------------------------------------------------
########### London ##################
london$estado <- cut(london$X2, 
                       breaks = c(-Inf, 3.001, 6.0001, Inf), 
                       labels = c("1", "2", "3"))
# Crear la cadena de Markov 
cadena_london <- markovchainFit(data = london$estado)
## Matriz de transición 
cadena_london$estimate
## Distribucion estacionaria
Plondon <- matrix(c(cadena_london$estimate[1,1], cadena_london$estimate[1,2],cadena_london$estimate[1,3],
                      cadena_london$estimate[2,1], cadena_london$estimate[2,2],cadena_london$estimate[2,3],
                      cadena_london$estimate[3,1], cadena_london$estimate[3,2],cadena_london$estimate[3,3]),
                    byrow = T, nrow = 3)
xtable(Plondon)
markovlondon<- new("markovchain", states = c("1", "2","3"), transitionMatrix =
                       Plondon,name = "A markovchain Object" )
disesta_london <- steadyStates(markovlondon);disesta_london
##Conteo y chi cuadrado.
matriz_nijlondon <- as.matrix(table(london$estado[-length(london$estado)], london$estado[-1]))
chi2_london <- chisq.test(matriz_nijlondon, p = Plondon^n)
##Tiempo medio de persistencia.
timeper_london<-round(c(1/(1-cadena_london$estimate[1,1]),1/(1-cadena_london$estimate[2,2]),1/(1-cadena_london$estimate[3,3])),2)
#-----------------------------------------------------
########### Mississauga ##################
mississauga$estado <- cut(mississauga$X2, 
                     breaks = c(-Inf, 3.001, 6.0001, Inf), 
                     labels = c("1", "2", "3"))
# Crear la cadena de Markov 
cadena_mississauga <- markovchainFit(data = mississauga$estado)
## Matriz de transición 
cadena_mississauga$estimate
## Distribucion estacionaria
Pmississauga <- matrix(c(cadena_mississauga$estimate[1,1], cadena_mississauga$estimate[1,2],cadena_mississauga$estimate[1,3],
                    cadena_mississauga$estimate[2,1], cadena_mississauga$estimate[2,2],cadena_mississauga$estimate[2,3],
                    cadena_mississauga$estimate[3,1], cadena_mississauga$estimate[3,2],cadena_mississauga$estimate[3,3]),
                  byrow = T, nrow = 3)
xtable(Pmississauga)
markovmississauga<- new("markovchain", states = c("1", "2","3"), transitionMatrix =
                     Pmississauga,name = "A markovchain Object" )
disesta_mississauga <- steadyStates(markovmississauga);disesta_mississauga
##Conteo y chi cuadrado.
matriz_nijmississauga <- as.matrix(table(mississauga$estado[-length(mississauga$estado)], mississauga$estado[-1]))
chi2_mississauga <- chisq.test(matriz_nijmississauga, p = Pmississauga^n)
##Tiempo medio de persistencia.
timeper_mississauga<-round(c(1/(1-cadena_mississauga$estimate[1,1]),1/(1-cadena_mississauga$estimate[2,2]),1/(1-cadena_mississauga$estimate[3,3])),2)
#-----------------------------------------------------
########### Newmarket ##################
newmarket$estado <- cut(newmarket$X2, 
                          breaks = c(-Inf, 3.001, 6.0001, Inf), 
                          labels = c("1", "2", "3"))
# Crear la cadena de Markov 
cadena_newmarket <- markovchainFit(data = newmarket$estado)
## Matriz de transición 
cadena_newmarket$estimate
## Distribucion estacionaria
Pnewmarket<- matrix(c(cadena_newmarket$estimate[1,1], cadena_newmarket$estimate[1,2],cadena_newmarket$estimate[1,3],
                         cadena_newmarket$estimate[2,1], cadena_newmarket$estimate[2,2],cadena_newmarket$estimate[2,3],
                         cadena_newmarket$estimate[3,1], cadena_newmarket$estimate[3,2],cadena_newmarket$estimate[3,3]),
                       byrow = T, nrow = 3)
xtable(Pnewmarket)
markovnewmarket<- new("markovchain", states = c("1", "2","3"), transitionMatrix =
                          Pnewmarket,name = "A markovchain Object" )
disesta_newmarket <- steadyStates(markovnewmarket);disesta_newmarket
##Conteo y chi cuadrado.
matriz_nijnewmarket <- as.matrix(table(newmarket$estado[-length(newmarket$estado)], newmarket$estado[-1]))
chi2_newmarket <- chisq.test(matriz_nijnewmarket, p = Pnewmarket^n)
##Tiempo medio de persistencia.
timeper_newmarket<-round(c(1/(1-cadena_newmarket$estimate[1,1]),1/(1-cadena_newmarket$estimate[2,2]),1/(1-cadena_newmarket$estimate[3,3])),2)
#-----------------------------------------------------
########### Northbay ##################
northbay$estado <- cut(northbay$X2, 
                        breaks = c(-Inf, 3.001, 6.0001, Inf), 
                        labels = c("1", "2", "3"))
# Crear la cadena de Markov 
cadena_northbay <- markovchainFit(data = northbay$estado)
## Matriz de transición 
cadena_northbay$estimate
## Distribucion estacionaria
Pnorthbay<- matrix(c(cadena_northbay$estimate[1,1], cadena_northbay$estimate[1,2],cadena_northbay$estimate[1,3],
                      cadena_northbay$estimate[2,1], cadena_northbay$estimate[2,2],cadena_northbay$estimate[2,3],
                      cadena_northbay$estimate[3,1], cadena_northbay$estimate[3,2],cadena_northbay$estimate[3,3]),
                    byrow = T, nrow = 3)
xtable(Pnorthbay)
markovnorthbay<- new("markovchain", states = c("1", "2","3"), transitionMatrix =
                        Pnorthbay,name = "A markovchain Object" )
disesta_northbay <- steadyStates(markovnorthbay);disesta_northbay
##Conteo y chi cuadrado.
matriz_nijnorthbay <- as.matrix(table(northbay$estado[-length(northbay$estado)], northbay$estado[-1]))
chi2_northbay <- chisq.test(matriz_nijnorthbay, p = Pnorthbay^n)
##Tiempo medio de persistencia.
timeper_northbay<-round(c(1/(1-cadena_northbay$estimate[1,1]),1/(1-cadena_northbay$estimate[2,2]),1/(1-cadena_northbay$estimate[3,3])),2)
#-----------------------------------------------------
########### Oakville ##################
oakville$estado <- cut(oakville$X2, 
                       breaks = c(-Inf, 3.001, 6.0001, Inf), 
                       labels = c("1", "2", "3"))
# Crear la cadena de Markov 
cadena_oakville <- markovchainFit(data = oakville$estado)
## Matriz de transición 
cadena_oakville$estimate
## Distribucion estacionaria
Poakville<- matrix(c(cadena_oakville$estimate[1,1], cadena_oakville$estimate[1,2],cadena_oakville$estimate[1,3],
                     cadena_oakville$estimate[2,1], cadena_oakville$estimate[2,2],cadena_oakville$estimate[2,3],
                     cadena_oakville$estimate[3,1], cadena_oakville$estimate[3,2],cadena_oakville$estimate[3,3]),
                   byrow = T, nrow = 3)
xtable(Poakville)
markovoakville<- new("markovchain", states = c("1", "2","3"), transitionMatrix =
                       Poakville,name = "A markovchain Object" )
disesta_oakville <- steadyStates(markovoakville);disesta_oakville
##Conteo y chi cuadrado.
matriz_nijoakville <- as.matrix(table(oakville$estado[-length(oakville$estado)], oakville$estado[-1]))
chi2_oakville <- chisq.test(matriz_nijoakville, p = Poakville^n)
##Tiempo medio de persistencia.
timeper_oakville<-round(c(1/(1-cadena_oakville$estimate[1,1]),1/(1-cadena_oakville$estimate[2,2]),1/(1-cadena_oakville$estimate[3,3])),2)
#-----------------------------------------------------
########### OttawaD ##################
ottawaD$estado <- cut(ottawaD$X2, 
                     breaks = c(-Inf, 3.001, 6.0001), 
                     labels = c("1", "2"))
# Crear la cadena de Markov 
cadena_ottawaD <- markovchainFit(data = ottawaD$estado)
## Matriz de transición 
cadena_ottawaD$estimate
## Distribucion estacionaria
PottawaD <- matrix(c(cadena_ottawaD$estimate[1,1], cadena_ottawaD$estimate[1,2],
                    cadena_ottawaD$estimate[2,1], cadena_ottawaD$estimate[2,2]),
                  byrow = T, nrow = 2)
xtable(PottawaD)
markovottawaD <- new("markovchain", states = c("1", "2"), transitionMatrix =
                      PottawaD,name = "A markovchain Object" )
disesta_ottawaD <- steadyStates(markovottawaD);disesta_ottawaD
##Conteo y chi cuadrado.
matriz_nijottawaD<- as.matrix(table(ottawaD$estado[-length(ottawaD$estado)], ottawaD$estado[-1]))
chi2_ottawaD <- chisq.test(matriz_nijottawaD, p = PottawaD^n)
##Tiempo medio de persistencia.
timeper_ottawaD<-round(c(1/(1-cadena_ottawaD$estimate[1,1]),1/(1-cadena_ottawaD$estimate[2,2])),2)
#-----------------------------------------------------
########### Parrysound ##################
parrysound$estado <- cut(parrysound$X2, 
                       breaks = c(-Inf, 3.001, 6.0001, Inf), 
                       labels = c("1", "2", "3"))
# Crear la cadena de Markov 
cadena_parrysound <- markovchainFit(data = parrysound$estado)
## Matriz de transición 
cadena_parrysound$estimate
## Distribucion estacionaria
Pparrysound<- matrix(c(cadena_parrysound$estimate[1,1], cadena_parrysound$estimate[1,2],cadena_parrysound$estimate[1,3],
                     cadena_parrysound$estimate[2,1], cadena_parrysound$estimate[2,2],cadena_parrysound$estimate[2,3],
                     cadena_parrysound$estimate[3,1], cadena_parrysound$estimate[3,2],cadena_parrysound$estimate[3,3]),
                   byrow = T, nrow = 3)
xtable(Pparrysound)
markovparrysound<- new("markovchain", states = c("1", "2","3"), transitionMatrix =
                       Pparrysound,name = "A markovchain Object" )
disesta_parrysound <- steadyStates(markovparrysound);disesta_parrysound
##Conteo y chi cuadrado.
matriz_nijparrysound<- as.matrix(table(parrysound$estado[-length(parrysound$estado)], parrysound$estado[-1]))
chi2_parrysound <- chisq.test(matriz_nijparrysound, p = Pparrysound^n)
##Tiempo medio de persistencia.
timeper_parrysound<-round(c(1/(1-cadena_parrysound$estimate[1,1]),1/(1-cadena_parrysound$estimate[2,2]),1/(1-cadena_parrysound$estimate[3,3])),2)
#-----------------------------------------------------
########### Peterbourgh ##################
peterbourgh$estado <- cut(peterbourgh$X2, 
                         breaks = c(-Inf, 3.001, 6.0001, Inf), 
                         labels = c("1", "2", "3"))
# Crear la cadena de Markov 
cadena_peterbourgh <- markovchainFit(data = peterbourgh$estado)
## Matriz de transición 
cadena_peterbourgh$estimate
## Distribucion estacionaria
Ppeterbourgh<- matrix(c(cadena_peterbourgh$estimate[1,1], cadena_peterbourgh$estimate[1,2],cadena_peterbourgh$estimate[1,3],
                       cadena_peterbourgh$estimate[2,1], cadena_peterbourgh$estimate[2,2],cadena_peterbourgh$estimate[2,3],
                       cadena_peterbourgh$estimate[3,1], cadena_peterbourgh$estimate[3,2],cadena_peterbourgh$estimate[3,3]),
                     byrow = T, nrow = 3)
xtable(Ppeterbourgh)
markovpeterbourgh<- new("markovchain", states = c("1", "2","3"), transitionMatrix =
                         Ppeterbourgh,name = "A markovchain Object" )
disesta_peterbourgh <- steadyStates(markovpeterbourgh);disesta_peterbourgh
##Conteo y chi cuadrado.
matriz_nijpeterbourgh<- as.matrix(table(peterbourgh$estado[-length(peterbourgh$estado)], peterbourgh$estado[-1]))
chi2_peterbourgh <- chisq.test(matriz_nijpeterbourgh, p = Ppeterbourgh^n)
##Tiempo medio de persistencia.
timeper_peterbourgh<-round(c(1/(1-cadena_peterbourgh$estimate[1,1]),1/(1-cadena_peterbourgh$estimate[2,2]),1/(1-cadena_peterbourgh$estimate[3,3])),2)
#-----------------------------------------------------
########### Portstanley ################## 
portstanley$estado <- cut(portstanley$X2, 
                          breaks = c(-Inf, 3.001, 6.0001, Inf), 
                          labels = c("1", "2", "3"))
# Crear la cadena de Markov 
cadena_portstanley <- markovchainFit(data = portstanley$estado)
## Matriz de transición 
cadena_portstanley$estimate
## Distribucion estacionaria
Pportstanley<- matrix(c(cadena_portstanley$estimate[1,1], cadena_portstanley$estimate[1,2],cadena_portstanley$estimate[1,3],
                        cadena_portstanley$estimate[2,1], cadena_portstanley$estimate[2,2],cadena_portstanley$estimate[2,3],
                        cadena_portstanley$estimate[3,1], cadena_portstanley$estimate[3,2],cadena_portstanley$estimate[3,3]),
                      byrow = T, nrow = 3)
xtable(Pportstanley)
markovportstanley<- new("markovchain", states = c("1", "2","3"), transitionMatrix =
                          Pportstanley,name = "A markovchain Object" )
disesta_portstanley <- steadyStates(markovportstanley);disesta_portstanley
##Conteo y chi cuadrado.
matriz_nijportstanley<- as.matrix(table(portstanley$estado[-length(portstanley$estado)], portstanley$estado[-1]))
chi2_portstanley <- chisq.test(matriz_nijportstanley, p = Pportstanley^n)
##Tiempo medio de persistencia.
timeper_portstanley<-round(c(1/(1-cadena_portstanley$estimate[1,1]),1/(1-cadena_portstanley$estimate[2,2]),1/(1-cadena_portstanley$estimate[3,3])),2)
#-----------------------------------------------------
########### Sarnia ##################
sarnia$estado <- cut(sarnia$X3, 
                          breaks = c(-Inf, 3.001, 6.0001, Inf), 
                          labels = c("1", "2", "3"))
# Crear la cadena de Markov 
cadena_sarnia <- markovchainFit(data = sarnia$estado)
## Matriz de transición 
cadena_sarnia$estimate
## Distribucion estacionaria
Psarnia<- matrix(c(cadena_sarnia$estimate[1,1], cadena_sarnia$estimate[1,2],cadena_sarnia$estimate[1,3],
                   cadena_sarnia$estimate[2,1], cadena_sarnia$estimate[2,2],cadena_sarnia$estimate[2,3],
                   cadena_sarnia$estimate[3,1], cadena_sarnia$estimate[3,2],cadena_sarnia$estimate[3,3]),
                      byrow = T, nrow = 3)
xtable(Psarnia)
markovsarnia<- new("markovchain", states = c("1", "2","3"), transitionMatrix =
                          Psarnia,name = "A markovchain Object" )
disesta_sarnia <- steadyStates(markovsarnia);disesta_sarnia
##Conteo y chi cuadrado.
matriz_nijsarnia<- as.matrix(table(sarnia$estado[-length(sarnia$estado)], sarnia$estado[-1]))
chi2_sarnia <- chisq.test(matriz_nijsarnia, p = Psarnia^n)
##Tiempo medio de persistencia.
timeper_sarnia<-round(c(1/(1-cadena_sarnia$estimate[1,1]),1/(1-cadena_sarnia$estimate[2,2]),1/(1-cadena_sarnia$estimate[3,3])),2)
#-----------------------------------------------------
########### Stchaterines ##################
stchaterines$estado <- cut(stchaterines$X2, 
                      breaks = c(-Inf, 3.001, 6.0001), 
                      labels = c("1", "2"))
# Crear la cadena de Markov 
cadena_stchaterines <- markovchainFit(data = stchaterines$estado)
## Matriz de transición 
cadena_stchaterines$estimate
## Distribucion estacionaria
Pstchaterines <- matrix(c(cadena_stchaterines$estimate[1,1], cadena_stchaterines$estimate[1,2],
                          cadena_stchaterines$estimate[2,1], cadena_stchaterines$estimate[2,2]),
                   byrow = T, nrow = 2)
xtable(Pstchaterines)
markovstchaterines <- new("markovchain", states = c("1", "2"), transitionMatrix =
                       Pstchaterines,name = "A markovchain Object" )
disesta_stchaterines<- steadyStates(markovstchaterines);disesta_stchaterines
##Conteo y chi cuadrado.
matriz_nijstchaterines<- as.matrix(table(stchaterines$estado[-length(stchaterines$estado)], stchaterines$estado[-1]))
chi2_stchaterines <- chisq.test(matriz_nijstchaterines, p = Pstchaterines^n)
##Tiempo medio de persistencia.
timeper_stchaterines<-round(c(1/(1-cadena_stchaterines$estimate[1,1]),1/(1-cadena_stchaterines$estimate[2,2])),2)
#-----------------------------------------------------
########### Sudbury ##################
sudbury$estado <- cut(sudbury$X2, 
                      breaks = c(-Inf, 3.001, 6.0001), 
                      labels = c("1", "2"))
# Crear la cadena de Markov 
cadena_sudbury <- markovchainFit(data = sudbury$estado)
## Matriz de transición 
cadena_sudbury$estimate
## Distribucion estacionaria
Psudbury <- matrix(c(cadena_sudbury$estimate[1,1], cadena_sudbury$estimate[1,2],
                     cadena_sudbury$estimate[2,1], cadena_sudbury$estimate[2,2]),
                   byrow = T, nrow = 2)
xtable(Psudbury)
markovsudbury <- new("markovchain", states = c("1", "2"), transitionMatrix =
                       Psudbury,name = "A markovchain Object" )
disesta_sudbury <- steadyStates(markovsudbury);disesta_sudbury
##Conteo y chi cuadrado.
matriz_nijsudbury<- as.matrix(table(sudbury$estado[-length(sudbury$estado)], sudbury$estado[-1]))
chi2_sudbury <- chisq.test(matriz_nijsudbury, p = Psudbury^n)
##Tiempo medio de persistencia.
timeper_sudbury<-round(c(1/(1-cadena_sudbury$estimate[1,1]),1/(1-cadena_sudbury$estimate[2,2])),2)
#-----------------------------------------------------
########### Thunderbay ##################
thunderbay$estado <- cut(thunderbay$X2, 
                          breaks = c(-Inf, 3.001, 6.0001, Inf), 
                          labels = c("1", "2", "3"))
# Crear la cadena de Markov 
cadena_thunderbay <- markovchainFit(data = thunderbay$estado)
## Matriz de transición 
cadena_thunderbay$estimate
## Distribucion estacionaria
Pthunderbay<- matrix(c(cadena_thunderbay$estimate[1,1], cadena_thunderbay$estimate[1,2],cadena_thunderbay$estimate[1,3],
                        cadena_thunderbay$estimate[2,1], cadena_thunderbay$estimate[2,2],cadena_thunderbay$estimate[2,3],
                        cadena_thunderbay$estimate[3,1], cadena_thunderbay$estimate[3,2],cadena_thunderbay$estimate[3,3]),
                      byrow = T, nrow = 3)
xtable(Pthunderbay)
markovthunderbay<- new("markovchain", states = c("1", "2","3"), transitionMatrix =
                          Pthunderbay,name = "A markovchain Object" )
disesta_thunderbay <- steadyStates(markovthunderbay);disesta_thunderbay
##Conteo y chi cuadrado.
matriz_nijthunderbay<- as.matrix(table(thunderbay$estado[-length(thunderbay$estado)], thunderbay$estado[-1]))
chi2_thunderbay <- chisq.test(matriz_nijthunderbay, p = Pthunderbay^n)
##Tiempo medio de persistencia.
timeper_thunderbay<-round(c(1/(1-cadena_thunderbay$estimate[1,1]),1/(1-cadena_thunderbay$estimate[2,2]),1/(1-cadena_thunderbay$estimate[3,3])),2)
#-----------------------------------------------------
########### Tiverton ##################
tiverton$estado <- cut(tiverton$X2, 
                         breaks = c(-Inf, 3.001, 6.0001, Inf), 
                         labels = c("1", "2", "3"))
# Crear la cadena de Markov 
cadena_tiverton <- markovchainFit(data = tiverton$estado)
## Matriz de transición 
cadena_tiverton$estimate
## Distribucion estacionaria
Ptiverton<- matrix(c(cadena_tiverton$estimate[1,1], cadena_tiverton$estimate[1,2],cadena_tiverton$estimate[1,3],
                     cadena_tiverton$estimate[2,1], cadena_tiverton$estimate[2,2],cadena_tiverton$estimate[2,3],
                     cadena_tiverton$estimate[3,1], cadena_tiverton$estimate[3,2],cadena_tiverton$estimate[3,3]),
                     byrow = T, nrow = 3)
xtable(Ptiverton)
markovtiverton<- new("markovchain", states = c("1", "2","3"), transitionMatrix =
                         Ptiverton,name = "A markovchain Object" )
disesta_tiverton <- steadyStates(markovtiverton);disesta_tiverton
##Conteo y chi cuadrado.
matriz_nijtiverton<- as.matrix(table(tiverton$estado[-length(tiverton$estado)], tiverton$estado[-1]))
chi2_tiverton <- chisq.test(matriz_nijtiverton, p = Ptiverton^n)
##Tiempo medio de persistencia.
timeper_tiverton<-round(c(1/(1-cadena_tiverton$estimate[1,1]),1/(1-cadena_tiverton$estimate[2,2]),1/(1-cadena_tiverton$estimate[3,3])),2)
#-----------------------------------------------------
########### TorontoE ##################
torontoE$estado <- cut(torontoE$X2, 
                       breaks = c(-Inf, 3.001, 6.0001, Inf), 
                       labels = c("1", "2", "3"))
# Crear la cadena de Markov 
cadena_torontoE <- markovchainFit(data = torontoE$estado)
## Matriz de transición 
cadena_torontoE$estimate
## Distribucion estacionaria
PtorontoE<- matrix(c(cadena_torontoE$estimate[1,1], cadena_torontoE$estimate[1,2],cadena_torontoE$estimate[1,3],
                     cadena_torontoE$estimate[2,1], cadena_torontoE$estimate[2,2],cadena_torontoE$estimate[2,3],
                     cadena_torontoE$estimate[3,1], cadena_torontoE$estimate[3,2],cadena_torontoE$estimate[3,3]),
                   byrow = T, nrow = 3)
xtable(PtorontoE)
markovtorontoE<- new("markovchain", states = c("1", "2","3"), transitionMatrix =
                       PtorontoE,name = "A markovchain Object" )
disesta_torontoE <- steadyStates(markovtorontoE);disesta_torontoE
##Conteo y chi cuadrado.
matriz_nijtorontoE<- as.matrix(table(torontoE$estado[-length(torontoE$estado)], torontoE$estado[-1]))
chi2_torontoE <- chisq.test(matriz_nijtorontoE, p = PtorontoE^n)
##Tiempo medio de persistencia.
timeper_torontoE<-round(c(1/(1-cadena_torontoE$estimate[1,1]),1/(1-cadena_torontoE$estimate[2,2]),1/(1-cadena_torontoE$estimate[3,3])),2)
#-----------------------------------------------------
########### TorontoW ##################
torontoW$estado <- cut(torontoW$X2, 
                       breaks = c(-Inf, 3.001, 6.0001, Inf), 
                       labels = c("1", "2", "3"))
# Crear la cadena de Markov 
cadena_torontoW <- markovchainFit(data = torontoW$estado)
## Matriz de transición 
cadena_torontoW$estimate
## Distribucion estacionaria
PtorontoW<- matrix(c(cadena_torontoW$estimate[1,1], cadena_torontoW$estimate[1,2],cadena_torontoW$estimate[1,3],
                     cadena_torontoW$estimate[2,1], cadena_torontoW$estimate[2,2],cadena_torontoW$estimate[2,3],
                     cadena_torontoW$estimate[3,1], cadena_torontoW$estimate[3,2],cadena_torontoW$estimate[3,3]),
                   byrow = T, nrow = 3)
xtable(PtorontoW)
markovtorontoW<- new("markovchain", states = c("1", "2","3"), transitionMatrix =
                       PtorontoW,name = "A markovchain Object" )
disesta_torontoW <- steadyStates(markovtorontoW);disesta_torontoW
##Conteo y chi cuadrado.
matriz_nijtorontoW<- as.matrix(table(torontoW$estado[-length(torontoW$estado)], torontoW$estado[-1]))
chi2_torontoW <- chisq.test(matriz_nijtorontoW, p = PtorontoW^n)
##Tiempo medio de persistencia.
timeper_torontoW<-round(c(1/(1-cadena_torontoW$estimate[1,1]),1/(1-cadena_torontoW$estimate[2,2]),1/(1-cadena_torontoW$estimate[3,3])),2)
#-----------------------------------------------------
########### WindsorD##################
windsorD$estado <- cut(windsorD$X3, 
                       breaks = c(-Inf, 3.001, 6.0001, Inf), 
                       labels = c("1", "2", "3"))
# Crear la cadena de Markov 
cadena_windsorD <- markovchainFit(data = windsorD$estado)
## Matriz de transición 
cadena_windsorD$estimate
## Distribucion estacionaria
PwindsorD<- matrix(c(cadena_windsorD$estimate[1,1], cadena_windsorD$estimate[1,2],cadena_windsorD$estimate[1,3],
                     cadena_windsorD$estimate[2,1], cadena_windsorD$estimate[2,2],cadena_windsorD$estimate[2,3],
                     cadena_windsorD$estimate[3,1], cadena_windsorD$estimate[3,2],cadena_windsorD$estimate[3,3]),
                   byrow = T, nrow = 3)
xtable(windsorD)
markovwindsorD<- new("markovchain", states = c("1", "2","3"), transitionMatrix =
                       PwindsorD,name = "A markovchain Object" )
disesta_windsorD <- steadyStates(markovwindsorD);disesta_windsorD
##Conteo y chi cuadrado.
matriz_nijwindsorD<- as.matrix(table(windsorD$estado[-length(windsorD$estado)], windsorD$estado[-1]))
chi2_windsorD <- chisq.test(matriz_nijwindsorD, p = PwindsorD^n)
##Tiempo medio de persistencia.
timeper_windsorD<-round(c(1/(1-cadena_windsorD$estimate[1,1]),1/(1-cadena_windsorD$estimate[2,2]),1/(1-cadena_windsorD$estimate[3,3])),2)
#-----------------------------------------------------
########### WindsorW##################
windsorW$estado <- cut(windsorW$X3, 
                       breaks = c(-Inf, 3.001, 6.0001, Inf), 
                       labels = c("1", "2", "3"))
# Crear la cadena de Markov 
cadena_windsorW <- markovchainFit(data = windsorW$estado)
## Matriz de transición 
cadena_windsorW$estimate
## Distribucion estacionaria
PwindsorW<- matrix(c(cadena_windsorW$estimate[1,1], cadena_windsorW$estimate[1,2],cadena_windsorW$estimate[1,3],
                     cadena_windsorW$estimate[2,1], cadena_windsorW$estimate[2,2],cadena_windsorW$estimate[2,3],
                     cadena_windsorW$estimate[3,1], cadena_windsorW$estimate[3,2],cadena_windsorW$estimate[3,3]),
                   byrow = T, nrow = 3)
xtable(windsorW)
markovwindsorW<- new("markovchain", states = c("1", "2","3"), transitionMatrix =
                       PwindsorW,name = "A markovchain Object" )
disesta_windsorW <- steadyStates(markovwindsorW);disesta_windsorD
##Conteo y chi cuadrado.
matriz_nijwindsorW<- as.matrix(table(windsorW$estado[-length(windsorW$estado)], windsorW$estado[-1]))
chi2_windsorW <- chisq.test(matriz_nijwindsorW, p = PwindsorW^n)
##Tiempo medio de persistencia.
timeper_windsorW<-round(c(1/(1-cadena_windsorW$estimate[1,1]),1/(1-cadena_windsorW$estimate[2,2]),1/(1-cadena_windsorW$estimate[3,3])),2)
#-----------------------------------------------------
########### Creación de tablas########
latlong<-read.table("latlong.txt", header = F, sep = "")[,2:4]
## Probabilidades de transición, distribución estacionaria y tiempo medio de persistencia #####
barrie_NA<-cbind(rbind(Pbarrie, c(0, 0)),c(0,0,0))
dorset_NA<-cbind(rbind(Pdorset, c(0, 0)),c(0,0,0))
ottawaD_NA<-cbind(rbind(PottawaD, c(0, 0)),c(0,0,0))
stcatherines_NA<-cbind(rbind(Pstchaterines, c(0, 0)),c(0,0,0))
sudbury_NA<-cbind(rbind(Psudbury, c(0, 0)),c(0,0,0))
pij<-rbind(c(t(barrie_NA)),c(t(Pbellevile)),c(t(Pbrantford)),c(t(Pburlington)),c(t(Pchatam)),
           c(t(Pcornwall)),c(t(dorset_NA)),c(t(Pgrandbend)),c(t(Pguelph)),c(t(PhamiltonD)),
           c(t(PhamiltonW)),c(t(Pkingston)),c(t(Pkitchner)),c(t(Plondon)),c(t(Pmississauga)),
           c(t(Pnewmarket)),c(t(Pnorthbay)),c(t(Poakville)),c(t(ottawaD_NA)),c(t(Pparrysound)),
           c(t(Ppeterbourgh)),c(t(Pportstanley)),c(t(Psarnia)),c(t(stcatherines_NA)),c(t(sudbury_NA)),#bien
           c(t(Pthunderbay)),c(t(Ptiverton)),c(t(PtorontoE)),c(t(PtorontoW)),c(t(PwindsorD)),
           c(t(PwindsorW)))
pi_ij<-rbind(c(disesta_barrie,0),disesta_bellevile,disesta_brantford,disesta_burlington,disesta_chatam,
             disesta_cornwall,c(disesta_dorset,0),disesta_grandbend,disesta_guelph,disesta_hamiltonD,
             disesta_hamiltonW,disesta_kingston,disesta_kitchner,disesta_london,disesta_mississauga,
             disesta_newmarket,disesta_northbay,disesta_oakville,c(disesta_ottawaD,0),disesta_parrysound,
             disesta_peterbourgh,disesta_portstanley,disesta_sarnia,c(disesta_stchaterines,0),c(disesta_sudbury,0),
             disesta_thunderbay,disesta_tiverton,disesta_torontoE,disesta_torontoW,disesta_windsorD,
             disesta_windsorW)
sij<-rbind(c(timeper_barrie,0),timeper_bellevile,timeper_brantford,timeper_burlington,timeper_chatam,
           timeper_cornwall,c(timeper_dorset,0),timeper_grandbend,timeper_guelph,timeper_hamiltonD,
           timeper_hamiltonW,timeper_kingston,timeper_kitchner,timeper_london,timeper_mississauga,
           timeper_newmarket,timeper_northbay,timeper_oakville,c(timeper_ottawaD,0),timeper_parrysound,
           timeper_peterbourgh,timeper_portstanley,timeper_sarnia,c(timeper_stchaterines,0),c(timeper_sudbury,0),
           timeper_thunderbay,timeper_tiverton,timeper_torontoE,timeper_torontoW,timeper_windsorD,
           timeper_windsorW)
### Tabla de conteos y valor chi cuadrado #####
barrie_NA2<-cbind(rbind(matriz_nijbarrie, c(0, 0)),c(0,0,0))
dorset_NA2<-cbind(rbind(matriz_nijdorset, c(0, 0)),c(0,0,0))
ottawaD_NA2<-cbind(rbind(matriz_nijottawaD, c(0, 0)),c(0,0,0))
stcatherines_NA2<-cbind(rbind(matriz_nijstchaterines, c(0, 0)),c(0,0,0))
sudbury_NA2<-cbind(rbind(matriz_nijsudbury, c(0, 0)),c(0,0,0))
nij<-rbind(c(t(barrie_NA2)),c(t(matriz_nijbellevile)),c(t(matriz_nijbrantford)),c(t(matriz_nijburlington)),c(t(matriz_nijchatam)),
           c(t(matriz_nijcornwall)),c(t(dorset_NA2)),c(t(matriz_nijgrandbend)),c(t(matriz_nijguelph)),c(t(matriz_nijhamiltonD)),
           c(t(matriz_nijhamiltonW)),c(t(matriz_nijkingston)),c(t(matriz_nijkitchner)),c(t(matriz_nijlondon)),c(t(matriz_nijmississauga)),
           c(t(matriz_nijnewmarket)),c(t(matriz_nijnorthbay)),c(t(matriz_nijoakville)),c(t(ottawaD_NA2)),c(t(matriz_nijparrysound)),
           c(t(matriz_nijpeterbourgh)),c(t(matriz_nijportstanley)),c(t(matriz_nijsarnia)),c(t(stcatherines_NA2)),c(t(sudbury_NA2)),#bien
           c(t(matriz_nijthunderbay)),c(t(matriz_nijtiverton)),c(t(matriz_nijtorontoE)),c(t(matriz_nijtorontoW)),c(t(matriz_nijwindsorD)),
           c(t(matriz_nijwindsorW)))
nij2<-cbind(latlong$Estación,nij)
pij2<-cbind(latlong$Estación,round(pij,2))
piijsij<-cbind(latlong$Estación,round(pi_ij,3),round(sij,0))
###---------------------------------------
latlong<-cbind(latlong,pij,pi_ij,sij)
names(latlong)<-c("Estación","latitud","longitud","P11","P12","P13",
                 "P21","P22","P23","P31","P32","P33","Pi11","Pi22","Pi33",
                 "S11","S22","S33")
##############Construcción de los mapas ####################
install.packages("rgeos")
install.packages("maptools")
install.packages("ggmap")
library(rgeos)
library(maptools)
library(ggmap)
library(leaflet)
library(sf)
sp.Canada <- readShapeSpatial('C:/Users/andre/OneDrive/Escritorio/Andrés/2023-2/Seminario/canada geo/gpr_000a06a_e/gpr_000a06a_e.shp',
                              proj4string = CRS("+proj=longlat +datum=WGS84"))
sp.Ontario <- sp.Canada[sp.Canada@data$PRNAME == "Ontario",]
Ontario <- st_as_sf(sp.Ontario)
####################
#install.packages("leaflet")
#devtools::install_github("bhaskarvk/leaflet.extras")
library(leaflet)
library(leaflet.extras)
ontario_map <- leaflet() %>%
  setView(lng = -79.3832, lat = 43.6532, zoom = 6) 

ontario_map <- ontario_map %>%
  addPolygons(data = Ontario, fillColor = "transparent", color = "black",weight = 1,smoothFactor = 0.1)  # Relleno transparente y borde azul

# Agregar los puntos de las estaciones al mapa
ontario_map <- ontario_map %>%
  addCircleMarkers(
    data = latlong,
    lng = ~longitud,
    lat = ~latitud,
    radius = 4,       # Tamaño de los marcadores
    color = "purple",    # Color de los marcadores
    fillOpacity = 4,  # Opacidad de relleno
    popup = ~paste("Estación:", Estación, "<br>Dias:", round(S11, 3))
  ) %>%
  addHeatmap(
    data = latlong,
    lng = ~longitud,
    lat = ~latitud,
    intensity = ~round(S11,3),
    blur = 25,        # Radio del desenfoque del calor
    max = 25,          # Valor máximo de intensidad
    radius = 30       # Radio de influencia del calor
  )#%>%
 # addMarkers(
    #data = latlong,
    #lng = ~longitud,
    #lat = ~latitud,
    #label = ~Estación,
    #labelOptions = labelOptions(direction = "top")
  #)
ontario_map
round(colMeans(latlong[,4:18]),3)













































