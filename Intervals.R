#### Library ####

library(ggplot2)

#### File location ####

setwd("C:/Users/rochi/OneDrive/Ambiente de Trabalho/Tese 2/Resultados/Resultados")

#### Loading data ####

dataV <- read.table("Results-V-rho01-tab2x3-MarginalExemplo.txt")
dataW <- read.table("Results-W-rho01-tab2x3-MarginalExemplo.txt")
dataVtilde <- read.table("Results-Vtilde-rho01-tab2x3-MarginalExemplo.txt")

#### V^2 and Vtilde^2 to V e Vtilde ####

dataV[5:12][dataV[5:12]<0] <- 0
dataV[5:12] <- sqrt(dataV[5:12])

dataVtilde[5:12][dataVtilde[5:12]<0] <- 0
dataVtilde[5:12] <- sqrt(dataVtilde[5:12])

#### probability calculation ####

m=30
mm=30
M=30
ProbC=matrix(0,m,mm)

par(mfrow=c(1,3))

#### V
  
# percentil v
for (i in 1:m) {
  for (j in 1:mm) {
    ProbC[i,j]=sum(dataV$V_Mae[((i-1)*M*mm+(j-1)*M+1):((i-1)*M*mm+(j-1)*M+M)]
                   >=dataV$LimInf_perc[((i-1)*M*mm+(j-1)*M+1):((i-1)*M*mm+(j-1)*M+M)]
                   & dataV$V_Mae[((i-1)*M*mm+(j-1)*M+1):((i-1)*M*mm+(j-1)*M+M)]
                   <=dataV$LimSup_perc[((i-1)*M*mm+(j-1)*M+1):((i-1)*M*mm+(j-1)*M+M)])/30*100
  }
  
}
#apply(ProbC,1,mean)
plot(1:30,apply(ProbC,1,mean), main="Medida V - Método Percentile", xlab="Tabelas Mãe", ylab="Probabilidade de cobertura",
     xlim=c(1,30), ylim=c(0, 100))
  
# bca v
  
for (i in 1:m) {
  for (j in 1:mm) {
    ProbC[i,j]=sum(dataV$V_Mae[((i-1)*M*mm+(j-1)*M+1):((i-1)*M*mm+(j-1)*M+M)]
                   >=dataV$LimInf_bca[((i-1)*M*mm+(j-1)*M+1):((i-1)*M*mm+(j-1)*M+M)]
                   & dataV$V_Mae[((i-1)*M*mm+(j-1)*M+1):((i-1)*M*mm+(j-1)*M+M)]
                   <=dataV$LimSup_bca[((i-1)*M*mm+(j-1)*M+1):((i-1)*M*mm+(j-1)*M+M)])/30*100
  }
  
}
#apply(ProbC,1,mean)
plot(1:30,apply(ProbC,1,mean),  main="Medida V - Método BCa", xlab="Tabelas Mãe", ylab="Probabilidade de cobertura",
     xlim=c(1,30), ylim=c(0, 100)) 
  
# basic v
for (i in 1:m) {
  for (j in 1:mm) {
    ProbC[i,j]=sum(dataV$V_Mae[((i-1)*M*mm+(j-1)*M+1):((i-1)*M*mm+(j-1)*M+M)]
                   >=dataV$LimInf_basic[((i-1)*M*mm+(j-1)*M+1):((i-1)*M*mm+(j-1)*M+M)]
                   & dataV$V_Mae[((i-1)*M*mm+(j-1)*M+1):((i-1)*M*mm+(j-1)*M+M)]
                   <=dataV$LimSup_basic[((i-1)*M*mm+(j-1)*M+1):((i-1)*M*mm+(j-1)*M+M)])/30*100
  }
  
}
#apply(ProbC,1,mean)
plot(1:30,apply(ProbC,1,mean),  main="Medida V - Método Basic", xlab="Tabelas Mãe", ylab="Probabilidade de cobertura",
     xlim=c(1,30), ylim=c(0, 100)) 

#### Vtilde
  
# percentil Vtilde
for (i in 1:m) {
  for (j in 1:mm) {
    ProbC[i,j]=sum(dataVtilde$Vtilde_mae[((i-1)*M*mm+(j-1)*M+1):((i-1)*M*mm+(j-1)*M+M)]
                   >=dataVtilde$LimInf_perc[((i-1)*M*mm+(j-1)*M+1):((i-1)*M*mm+(j-1)*M+M)]
                   & dataVtilde$Vtilde_mae[((i-1)*M*mm+(j-1)*M+1):((i-1)*M*mm+(j-1)*M+M)]
                   <=dataVtilde$LimSup_perc[((i-1)*M*mm+(j-1)*M+1):((i-1)*M*mm+(j-1)*M+M)])/30*100
  }
  
}
#apply(ProbC,1,mean)
plot(1:30,apply(ProbC,1,mean), main="Medida Vtilde - Método Percentile", xlab="Tabelas Mãe", ylab="Probabilidade de cobertura",
     xlim=c(1,30), ylim=c(0, 100)) 
  
# bca Vtilde
for (i in 1:m) {
  for (j in 1:mm) {
    ProbC[i,j]=sum(dataVtilde$Vtilde_mae[((i-1)*M*mm+(j-1)*M+1):((i-1)*M*mm+(j-1)*M+M)]
                   >=dataVtilde$LimInf_bca[((i-1)*M*mm+(j-1)*M+1):((i-1)*M*mm+(j-1)*M+M)]
                   & dataVtilde$Vtilde_mae[((i-1)*M*mm+(j-1)*M+1):((i-1)*M*mm+(j-1)*M+M)]
                   <=dataVtilde$LimSup_bca[((i-1)*M*mm+(j-1)*M+1):((i-1)*M*mm+(j-1)*M+M)])/30*100
  }
  
}
#apply(ProbC,1,mean)
plot(1:30,apply(ProbC,1,mean), main="Medida Vtilde - Método BCa", xlab="Tabelas Mãe", ylab="Probabilidade de cobertura",
     xlim=c(1,30), ylim=c(0, 100))   
  
  
# basic Vtilde
for (i in 1:m) {
  for (j in 1:mm) {
    ProbC[i,j]=sum(dataVtilde$Vtilde_mae[((i-1)*M*mm+(j-1)*M+1):((i-1)*M*mm+(j-1)*M+M)]
                   >=dataVtilde$LimInf_perc[((i-1)*M*mm+(j-1)*M+1):((i-1)*M*mm+(j-1)*M+M)]
                   & dataVtilde$Vtilde_mae[((i-1)*M*mm+(j-1)*M+1):((i-1)*M*mm+(j-1)*M+M)]
                   <=dataVtilde$LimSup_perc[((i-1)*M*mm+(j-1)*M+1):((i-1)*M*mm+(j-1)*M+M)])/30*100
  }
  
}
#apply(ProbC,1,mean)
plot(1:30,apply(ProbC,1,mean), main="Medida Vtilde - Método Basic", xlab="Tabelas Mãe", ylab="Probabilidade de cobertura",
     xlim=c(1,30), ylim=c(0, 100)) 
  

#### W

# percentil w
for (i in 1:m) {
  for (j in 1:mm) {
    ProbC[i,j]=sum(dataW$W_mae[((i-1)*M*mm+(j-1)*M+1):((i-1)*M*mm+(j-1)*M+M)]
                   >=dataW$LimInf_perc[((i-1)*M*mm+(j-1)*M+1):((i-1)*M*mm+(j-1)*M+M)]
                   & dataW$W_mae[((i-1)*M*mm+(j-1)*M+1):((i-1)*M*mm+(j-1)*M+M)]
                   <=dataW$LimSup_perc[((i-1)*M*mm+(j-1)*M+1):((i-1)*M*mm+(j-1)*M+M)])/30*100
  }
  
}
#apply(ProbC,1,mean)
plot(1:30,apply(ProbC,1,mean), main="Medida W - Método Percentile", xlab="Tabelas Mãe", ylab="Probabilidades de cobertura",
     xlim=c(1,30), ylim=c(0, 100)) 

# bca w
for (i in 1:m) {
  for (j in 1:mm) {
    ProbC[i,j]=sum(dataW$W_mae[((i-1)*M*mm+(j-1)*M+1):((i-1)*M*mm+(j-1)*M+M)]
                   >=dataW$LimInf_bca[((i-1)*M*mm+(j-1)*M+1):((i-1)*M*mm+(j-1)*M+M)]
                   & dataW$W_mae[((i-1)*M*mm+(j-1)*M+1):((i-1)*M*mm+(j-1)*M+M)]
                   <=dataW$LimSup_bca[((i-1)*M*mm+(j-1)*M+1):((i-1)*M*mm+(j-1)*M+M)])/30*100
  }
  
}
#apply(ProbC,1,mean)
plot(1:30,apply(ProbC,1,mean), main="Medida W - Método BCa", xlab="Tabelas Mãe", ylab="Probabilidades de cobertura",
     xlim=c(1,30), ylim=c(0, 100))   


# basic w
for (i in 1:m) {
  for (j in 1:mm) {
    ProbC[i,j]=sum(dataW$W_mae[((i-1)*M*mm+(j-1)*M+1):((i-1)*M*mm+(j-1)*M+M)]
                   >=dataW$LimInf_basic[((i-1)*M*mm+(j-1)*M+1):((i-1)*M*mm+(j-1)*M+M)]
                   & dataW$W_mae[((i-1)*M*mm+(j-1)*M+1):((i-1)*M*mm+(j-1)*M+M)]
                   <=dataW$LimSup_basic[((i-1)*M*mm+(j-1)*M+1):((i-1)*M*mm+(j-1)*M+M)])/30*100
  }
  
}
#apply(ProbC,1,mean)
plot(1:30,apply(ProbC,1,mean), main="Medida W - Método Basic", xlab="Tabelas Mãe", ylab="Probabilidades de cobertura",
     xlim=c(1,30), ylim=c(0, 100)) 


#### Amplitude dos intervalos ####
  
V_amp_percent<-c()
V_amp_bca<-c()
V_amp_basic<-c()

Vtilde_amp_percent<-c()
Vtilde_amp_bca<-c()
Vtilde_amp_basic<-c()

W_amp_percent<-c()
W_amp_bca<-c()
W_amp_basic<-c()

for (i in 1:27000) {
  V_amp_percent<-c(V_amp_percent,dataV$LimSup_perc[i]-dataV$LimInf_perc[i])
  V_amp_bca<-c(V_amp_bca,dataV$LimSup_bca[i]-dataV$LimInf_bca[i])
  V_amp_basic<-c(V_amp_basic,dataV$LimSup_basic[i]-dataV$LimInf_basic[i])
  
  Vtilde_amp_percent<-c(Vtilde_amp_percent,dataVtilde$LimSup_perc[i]-dataVtilde$LimInf_perc[i])
  Vtilde_amp_bca<-c(Vtilde_amp_bca,dataVtilde$LimSup_bca[i]-dataVtilde$LimInf_bca[i])
  Vtilde_amp_basic<-c(Vtilde_amp_basic,dataVtilde$LimSup_basic[i]-dataVtilde$LimInf_basic[i])
  
  W_amp_percent<-c(W_amp_percent,dataW$LimSup_perc[i]-dataW$LimInf_perc[i])
  W_amp_bca<-c(W_amp_bca,dataW$LimSup_bca[i]-dataW$LimInf_bca[i])
  W_amp_basic<-c(W_amp_basic,dataW$LimSup_basic[i]-dataW$LimInf_basic[i])
}

Amplitude=matrix(0,nrow=3,ncol=3)
rownames(Amplitude) <- c("V","Vtilde","W")
colnames(Amplitude) <- c("Percentile", "BCa","W")

Amplitude[1,1] <- mean(V_amp_percent)
Amplitude[1,2] <- mean(V_amp_bca)
Amplitude[1,3] <- mean(V_amp_basic)

Amplitude[2,1] <- mean(Vtilde_amp_percent)
Amplitude[2,2] <- mean(Vtilde_amp_bca)
Amplitude[2,3] <- mean(Vtilde_amp_basic)

Amplitude[3,1] <- mean(W_amp_percent)
Amplitude[3,2] <- mean(W_amp_bca)
Amplitude[3,3] <- mean(W_amp_basic)

Amplitude
