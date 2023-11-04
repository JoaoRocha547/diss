rm()

#### Library ####

library(boot)
library(vcd)
library(foreach)
library(doParallel)

#### Functions ####

# Intervalo (min(x), max(x)) e' dividido em nr*gx subintervalos uniformes
Intervalos=function(x,gx,nr)   
{
  Int=c(min(x))
  for (i in 1:(nr*gx)){Int[i+1]=quantile(x,i/(nr*gx))}
  Int[1:(nr*gx+1)]
}

# Intervalo (min(x), max(x)) e' dividido obedecendo a uma "marginal"
Intervalos_marginal=function(x,marginal)   
{
  nm=length(marginal)
  Int=c(min(x))
  for (i in 1:(nm)){Int[i+1]=quantile(x,sum(marginal[1:i]))}
  Int[1:(nm+1)]
}

# Separa o conjunto x de acordo com uma distribuição uniforme
Classes= function(x,gx,nr)
{
  n=length(x)
  xx=NULL
  for (i in 1:n){xx[i]= findInterval(x[i], Intervalos(x,gx,nr))}
  xx
  xx[xx==nr*gx+1]= nr*gx
  aux=sort(unique(xx))
  Categoria=list() ; for (i in 1:nr){Categoria[[i]]=aux[((i-1)*gx+1):(i*gx)]}
  yy=xx
  for (i in 1:nr) {yy[which(xx%in%Categoria[[i]])]=i}
  xx=yy
  xx
}

# Separa o conjunto x de acordo com a "marginal"
Classes_marginal= function(x,marginal)
{
  n=length(x)
  xx=NULL
  for (i in 1:n){xx[i]= findInterval(x[i], Intervalos_marginal(x,marginal))}
  xx[xx==length(marginal)+1]=length(marginal)
  xx
}

# Separa o conjunto x de acordo com uma distribuição normal
# seguindo o paper: Hyeong Chul Jeong, Myoungshic Jhun, Daehak Kim,
# Bootstrap tests for independence in two-way
# contingency tables.Computational Statistics & Data Analysis,
# 48 (2005) 623 – 631. (Section 3.1)
Classes_norm= function(x,nr)
{
  n=length(x)
  pontecorte2=c(min(x),0,max(x))
  pontecorte3=c(min(x),-0.6,0.6,max(x))
  pontecorte4=c(min(x),-0.8,0,0.8,max(x))
  pontecorte5=c(min(x),-0.6,-0.3,0.3,0.6,max(x))
  xx=NULL
  if (nr==2) {
    for (i in 1:n){xx[i]= findInterval(x[i], pontecorte2)}
  }
  if (nr==3) {
    for (i in 1:n){xx[i]= findInterval(x[i], pontecorte3)}
  }
  else if (nr==4) {
    for (i in 1:n){xx[i]= findInterval(x[i], pontecorte4)}
  }
  else if (nr==5) {
    for (i in 1:n){xx[i]= findInterval(x[i], pontecorte5)}
  }
  xx
  xx[xx==nr+1]= nr
  xx
}

### Calculo do coeficiente V de Cramer (cl?ssico)
FuncaoCalcV = function(input, index){
  Input = table(input[index,1],input[index,2])
  InputInicial=Input
  id.row.zeros=which(rowSums(Input) == 0)
  if (length(id.row.zeros)!=0){Input=Input[-c(id.row.zeros),]}
  if (length(Input)==dim(InputInicial)[2])
  {ResultV=0}else{id.col.zeros=which(colSums(Input) == 0)
  if (length(id.col.zeros)!=0){Input=Input[,-id.col.zeros]}
  if(length(Input)==dim(InputInicial)[1]){ResultV=0}else{ResultV=suppressWarnings(assocstats(Input)$cramer)}
  }
  return(ResultV)}

## Calculo medida alternativa W
#
# paper: Tarald O. Kvålseth (2018) An alternative to Cramér's coefficient of
# association, Communications in Statistics - Theory and Methods, 47:23, 5662-5674
# DOI:10.1080/03610926.2017.1400056
FuncaoCalcW = function(input, index){
    Input = table(input[index,1],input[index,2])
    InputInicial=Input
    id.row.zeros=which(rowSums(Input) == 0)
    if (length(id.row.zeros)!=0){Input=Input[-c(id.row.zeros),]}
    if (length(Input)==dim(InputInicial)[2])
    {ResultW=0}else{id.col.zeros=which(colSums(Input) == 0)
    if (length(id.col.zeros)!=0){Input=Input[,-id.col.zeros]}
    if(length(Input)==dim(InputInicial)[1]){ResultW=0}else{
      d=0
      ml=sum(apply(Input/sum(Input),1,sum)^2)
      mc=sum(apply(Input/sum(Input),2,sum)^2)
      w1=sum((Input/sum(Input))^2)
      for (a in 1:nrow(Input)) {
        for (b in 1:ncol(Input)) {
          d = d + (Input[a,b]/sum(Input)-sum(Input[a,]/sum(Input))*sum(Input[,b]/sum(Input)))^2
        }
      }
      ws <- c(ml,mc)
      w2=min(ws)
      w=sqrt(d/(d-w1+w2))
      ResultW=suppressWarnings(w)}
    }
  return(ResultW)}

## Calculo Cramer Vtilde (Cramer corrigido do vi�s)
#
# paper: Bergsma(2013) A bias-correction for Cram�r's V and Tschuprow's T
# Journal of the Korean Statistical Society, Volume 42, Issue 3, 323-328
# DOI: 10.1016/j.jkss.2012.10.002
FuncaoCalcVtilde = function(input, index){
  Input = table(input[index,1],input[index,2])
  InputInicial=Input
  id.row.zeros=which(rowSums(Input) == 0)
  if (length(id.row.zeros)!=0){Input=Input[-c(id.row.zeros),]}
  if (length(Input)==dim(InputInicial)[2])
  {ResultVtilde=0}else{id.col.zeros=which(colSums(Input) == 0)
  if (length(id.col.zeros)!=0){Input=Input[,-id.col.zeros]}
  if(length(Input)==dim(InputInicial)[1]){ResultVtilde=0}else{
    phiTildeMais=max(0,suppressWarnings(assocstats(Input)$chisq_tests[2,1])/sum(Input)
                     -1/(sum(Input)-1)*(dim(Input)[1]-1)*(dim(Input)[2]-1))
    rtilde= dim(Input)[1]-(dim(Input)[1]-1)^2/(sum(Input)-1)
    ctilde= dim(Input)[2]-(dim(Input)[2]-1)^2/(sum(Input)-1)
    ResultVtilde=sqrt(phiTildeMais / min(rtilde-1,ctilde-1))}
  }
  return(ResultVtilde)}

#### Simulation Study part1: Geração de tcs Mães ####

#Número de observações que permitirão gerar cada tc Mãe
NNN=1000  

#tamanho:tabelas 2x3 ou 2x5
nr=2

#TIPOS de marginais:
Marginal="Exemplo"

#Correlação da normal bivariada que origina as tcs
rho=.1 #.5#.9

if (Marginal=="Exemplo")
{
  marginalX=c(0.4,0.6)
  marginalY=c(0.5,0.4,0.1)
  #nr=length(marginalX)  
  nc=length(marginalY)
}

if (Marginal=="Normal")
{
  nc=5  #3  
}  

if (Marginal=="Uniforme")
{
  nc=5  #5  
  gx=5
  gy=5
}

#### Geração de m tcs Mães e correspondentes estatisticas rho, V, W e Vtilde

m=30 #30 # numero de tcs Mães a gerar. Cada tc Mãe contem  NN observações

TC <- list()
TabResultadosTC=matrix(0,nrow=m,ncol=4)
colnames(TabResultadosTC)=c("rho","CramerV","AlternativeW","VTilde")

set.seed(626)
for (i in 1:m) 
{
    NNN=NNN
    x1=rnorm(NNN,0,1)
    x2=rnorm(NNN,0,1)
    a=(sqrt(1+rho)+sqrt(1-rho) )/2
    b=(sqrt(1+rho)-sqrt(1-rho) )/2
    x= a*x1 + b*x2
    y= b*x1 + a*x2
    
    #Construção da tc dependendo da marginal:
    if (Marginal=="Exemplo")
    {
      tc=table(Classes_marginal(x,marginalX),Classes_marginal(y,marginalY))
    }
    if (Marginal=="Normal")
    {
      tc=table(Classes_norm(x,nr),Classes_norm(y,nc))
    }  
    
    if (Marginal=="Uniforme")
    {
      tc=table(Classes(x,gx,nr),Classes(y,gy,nc))
    }
    
      TC[[i]]=tc
      
      #Cálculo da medida V de Cramer:
      V=suppressWarnings(assocstats(tc)$cramer)
      
      #Cálculo da medida alternativa W:
      d=0
      ml=sum(apply(tc/NNN,1,sum)^2)
      mc=sum(apply(tc/NNN,2,sum)^2)
      w1=sum((tc/NNN)^2)
      for (a in 1:nr) {
        for (b in 1:nc) {
          d = d + (tc[a,b]/NNN-sum(tc[a,]/NNN)*sum(tc[,b]/NNN))^2
        }
      }
      ws <- c(ml,mc)
      w2=min(ws)
      W=sqrt(d/(d-w1+w2))
      
      #C�lculo da medida de Cramer corrigida (Vtilde)
      phiTildeMais=max(0,suppressWarnings(assocstats(tc)$chisq_tests[2,1]/sum(tc))
                         -1/(sum(tc)-1)*(dim(tc)[1]-1)*(dim(tc)[2]-1))
      rtilde= dim(tc)[1]-(dim(tc)[1]-1)^2/(sum(tc)-1)
      ctilde= dim(tc)[2]-(dim(tc)[2]-1)^2/(sum(tc)-1)
      Vtilde=sqrt(phiTildeMais/ min(rtilde-1,ctilde-1))
   
      TabResultadosTC[i,1]=cor(x, y, method = c("pearson"))
      TabResultadosTC[i,2]=V
      TabResultadosTC[i,3]=W
      TabResultadosTC[i,4]=Vtilde
  }

par(mfrow=c(1,4))
plot(TabResultadosTC[,1], ylim=c(0,1), ylab="rho", xlab="N�mero de tabelas m�e")
abline(h=0.1, col="black")
plot(TabResultadosTC[,2], ylim=c(0,1), ylab="V de Cram�r", xlab="N�mero de tabelas m�e")
abline(h=0.1, col="black")
plot(TabResultadosTC[,3], ylim=c(0,1), ylab="W de Kv�lseth", xlab="N�mero de tabelas m�e")
abline(h=0.1, col="black")
plot(TabResultadosTC[,4], ylim=c(0,1), ylab="Vtilde", xlab="N�mero de tabelas m�e")
abline(h=0.1, col="black")

#### Settings das simulações bootstrap:

mm=30               #numero de tcs Filhas geradas para cada tc mãe
R=2000 # 1000       #numero de amostras bootstrap geradas para cada tc Filha
                    #com obtenção de 1 ICbootstrap por tipo e tc
M=30  #100          #numero de replicas para calcular probabilidades de cobertura




# Preparação do objeto do R onde irão ser depositados os resultados
# dos IC para, no final da simulação ser transferido para um ficheiro txt
# um ficheiro por cada medida

###V de Cramer:
ResultadosFinaisV=matrix(0,nrow=M*m*mm,ncol=12)
for (i in 1:m){
  ResultadosFinaisV[(M*mm*(i-1)+1):(M*mm*(i-1)+M*mm),1]=rep(i,M*mm)
  ResultadosFinaisV[(M*mm*(i-1)+1):(M*mm*(i-1)+M*mm),3]=rep(1:M,mm)
  ResultadosFinaisV[(M*mm*(i-1)+1):(M*mm*(i-1)+M*mm),4]=rep(TabResultadosTC[i,1],M*mm)
  ResultadosFinaisV[(M*mm*(i-1)+1):(M*mm*(i-1)+M*mm),5]=rep(TabResultadosTC[i,2],M*mm) 
  for (k in 1:mm){ResultadosFinaisV[(M*mm*(i-1)+1+(k-1)*M):(M*mm*(i-1)+M+(k-1)*M),2]=rep(k,M)}
}
colnames(ResultadosFinaisV)=c("tcMae","tcFilha","Replica","rho_Obs", "V_Mae","V0_boot",
                              "LimInf_perc","LimSup_perc",
                              "LimInf_basic","LimSup_basic",
                              "LimInf_bca","LimSup_bca")
#Medida Alternative W
ResultadosFinaisW=ResultadosFinaisV
for (i in 1:m){
  ResultadosFinaisW[(M*mm*(i-1)+1):(M*mm*(i-1)+M*mm),5]=rep(TabResultadosTC[i,3],M*mm) 
}
colnames(ResultadosFinaisW)[c(5,6)]=c("W_mae","W0_boot")

#Medida Alternative Vtilde
ResultadosFinaisVtilde=ResultadosFinaisV
for (i in 1:m){
  ResultadosFinaisVtilde[(M*mm*(i-1)+1):(M*mm*(i-1)+M*mm),5]=rep(TabResultadosTC[i,4],M*mm) 
}
colnames(ResultadosFinaisVtilde)[c(5,6)]=c("Vtilde_mae","Vtilde0_boot")

#### Bootstrap com R amostras bootstrap sobre cada uma das mm tc filhas geradas de cada uma das m tc mae

# Set up parallel processing
cl <- makeCluster(detectCores())
registerDoParallel(cl)
system.time({ 
x<-foreach(i = 1:m, .packages = c("vcd","boot")) %dopar% 
  { # m tcs Mães
    
    ABC=matrix(0,M*mm,21)
    for (k in 1:mm)   # filha k
    {
    zz=rmultinom(n=1, size=50, prob=as.vector(TC[[i]]/NNN))  #filha
    tabDados=as.table(matrix(zz,nrow=dim(TC[[i]])[1]))
    n=sum(tabDados)
    xx=as.data.frame(tabDados)
    NN=dim(tabDados)[1]*dim(tabDados)[2]
    nn=1
    raw.dados=NULL
    for (k1 in 1:NN)
    {
      if (xx$Freq[k1]!=0)
      {
        for (k2 in 1: xx$Freq[k1])
        {
          raw.dados=rbind(raw.dados,xx[k1,1:2])
        }
      }
    }
    
    
    ##
    # Uma boa explicação teórica sobre os IC bootstrap e sua obtenção com o R:
    #https://www.r-bloggers.com/2019/09/understanding-bootstrap-confidence-interval-output-from-the-r-boot-package/
    
    # M réplicas de amostragem bootstrap sobre cada uma das m tc geradas
    # para permitir calcular probabilidades de cobertura:
      
      for (j in 1:M) #para rho=0.1 tirar bca
      {     # M réplicas para cada tc
        
        #Bootstrap para a medida V de Cramer
        set.seed((i-1)*mm*M +(k-1)*M +j+3)
        BootV = boot(raw.dados, FuncaoCalcV, R=R)
        #IntV=boot.ci(BootV, type = c("basic", "perc","bca"))
        #erro BCa
        IntV=boot.ci(BootV, type = c("basic", "perc"))
        
        #Bootstrap para a medida alternativa W
        set.seed((i-1)*mm*M +(k-1)*M +j+3)
        BootW = boot(raw.dados, FuncaoCalcW, R=R)
        #IntW=boot.ci(BootW, type = c("basic", "perc","bca"))
        #erro BCa
        IntW=boot.ci(BootW, type = c("basic", "perc"))
        
        #Bootstrap para a medida alternativa Vtilde
        set.seed((i-1)*mm*M +(k-1)*M +j+3)
        BootVtilde = boot(raw.dados, FuncaoCalcVtilde, R=R)
        #IntVtilde=boot.ci(BootVtilde, type = c("basic", "perc","bca"))
        #erro BCa
        IntVtilde=boot.ci(BootVtilde, type = c("basic", "perc"))
        
        #ABC[(k-1)*M +j,1:21]=c(IntV$t0,IntV$percent[4:5],IntV$basic[4:5],
                               #IntV$bca[4:5],IntW$t0,IntW$percent[4:5],
                               #IntW$basic[4:5],IntW$bca[4:5],IntVtilde$t0,
                               #IntVtilde$percent[4:5],IntVtilde$basic[4:5],
                               #IntVtilde$bca[4:5]) 
        #erro BCa
        ABC[(k-1)*M +j,1:21]=c(IntV$t0,IntV$percent[4:5],IntV$basic[4:5],
                               0,0,IntW$t0,IntW$percent[4:5],
                               IntW$basic[4:5],0,0,IntVtilde$t0,
                               IntVtilde$percent[4:5],IntVtilde$basic[4:5],
                               0,0) 
        
      }  # fim do ciclo j para a replica
    
    }  #fim do ciclo k
    
    return(ABC)    
  }  # fim d2o ciclo i para a obtenção dos resultados de m tcs mãe
})
x
stopCluster(cl)


#Input dos resultados nas tabelas de resultados
for (i in 1:m) {
  ResultadosFinaisV[(M*mm*(i-1)+1):(M*mm*(i-1)+M*mm),6:12] <- x[[i]][,1:7]
  ResultadosFinaisW[(M*mm*(i-1)+1):(M*mm*(i-1)+M*mm),6:12] <- x[[i]][,8:14]
  ResultadosFinaisVtilde[(M*mm*(i-1)+1):(M*mm*(i-1)+M*mm),6:12] <- x[[i]][,15:21]
}    

####Escrita dos resultados finais em ficheiros txt####l,sss'

#Para o V de Cramer
write.table(ResultadosFinaisV,file=paste("Results-V-rho0",rho*10,"-tab",
                                         nr,"x",nc,"-Marginal",Marginal,".txt",sep=""))
#Para a medida alternative W
write.table(ResultadosFinaisW,file=paste("Results-W-rho0",rho*10,"-tab",
                                         nr,"x",nc,"-Marginal",Marginal,".txt",sep=""))

#Para a medida alternative Vtilde
write.table(ResultadosFinaisVtilde,file=paste("Results-Vtilde-rho0",rho*10,"-tab",
                                         nr,"x",nc,"-Marginal",Marginal,".txt",sep=""))





