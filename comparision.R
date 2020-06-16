require(asreml)
require(dae)
require(MASS)

################################### simulation ###############################
column <- factor(rep(1:10,each=6))
row <- factor(rep(1:6,10))
treat <- factor(sample(rep(1:12,5)))
VarianceMatrix <- kronecker(mat.ar1(0.5,10),mat.banded(c(1,-0.3),6,6))
EnvEffect <- t(mvrnorm(1000,matrix(100,nrow=60),VarianceMatrix))
Ug <- rnorm(12,0,0.8)
Genetic <- model.matrix(~treat-1)%*% Ug
Y <- matrix(NA, ncol=1000,nrow=60)
for(i in 1:1000)
Y[,i] <- Genetic + EnvEffect[,i]
SimulEffect <- data.frame(column,row, treat, Y)
summary(SimulEffect)[,1:5]

##############################################################################
####################   Analysis ##############################################
##############################################################################
banded <- function(order, kappa) {
  H <- mat.banded(c(0,1),order,order)
  V <- kappa^H - (mat.J(order)-mat.banded(c(1,1),order,order))
  ## derivative
  dV <- H*(kappa^(H-1))
  return(list(V, dV))
}

CBand2 <- function(order, kappa) {
  H <- mat.banded(c(0,1,rep(1000,order-2)),order,order)
  V <- kappa^H 
  ## derivative
  dV <- H*(kappa^(H-1))
  return(list(V, dV))
}

### for design 

  parametOWN <- matrix(NA,ncol=5,nrow = 1000)
  parametOWN2 <- matrix(NA,ncol=5,nrow = 1000)
  parametBAND <- matrix(NA,ncol=5,nrow = 1000)
  for(i in 989:1000)
  {
    SimulEffect$Y <- SimulEffect[,(i+3)]
    model_op <- asreml(Y ~ 1, random =~ treat,
                       residual = ~ ar1(column):own(row,"banded",0.1,"R"),
                       aom=T,trace=F,data = SimulEffect,maxit=30,gammaPar=TRUE)
    if(model_op$converge == TRUE)
    {
      parametOWN[i,1:4] <- summary(model_op)$varcomp[,1]
      parametOWN[i,5] <- summary(model_op)$loglik
    } 
    model2 <- asreml(Y ~ 1, random =~ treat,
                       residual = ~ ar1(column):own(row,"CBand2",0.1,"R"),
                       aom=T,trace=F,data = SimulEffect,maxit=30,gammaPar=TRUE)
    if(model2$converge == TRUE)
    {
      parametOWN2[i,1:4] <- summary(model2)$varcomp[,1]
      parametOWN2[i,5] <- summary(model2)$loglik
    } 
    model_x <- asreml(Y ~ 1, random =~ treat,
                       residual = ~ ar1(column):corb(row,1),
                       aom=T,trace=F,data = SimulEffect,maxit=30,gammaPar=TRUE)
    if(model_x$converge == TRUE)
    {
      parametBAND[i,1:4] <- summary(model_x)$varcomp[,1]
      parametBAND[i,5] <- summary(model_x)$loglik
    } 
  }
res<-summary(parametBAND)
  summary(parametOWN)
  summary(parametBAND- parametOWN)[1:6,]
  summary(parametBAND- parametOWN2)[1:6,]
 par(mfrow=c(1,2))
   boxplot(parametBAND- parametOWN, xaxt = "n",main="Band-OWn1") 
text((1:5+0.5),-3,pos = 2, cex = 0.9,srt = 30, xpd = TRUE, 
     labels = c("genetic","residual","column cor","row cor","log lik"))
boxplot(parametOWN- parametOWN2, main="OWN1 - OWN2",xaxt = "n") 
text((1:5+0.5),-1.2,pos = 2, cex = 0.9,srt = 30, xpd = TRUE, 
     labels = c("genetic","residual","column cor","row cor","log lik"))

par(mfrow=c(2,2))
boxplot(cbind(parametBAND[,1],parametOWN[,1]), main="treat - 0.8",xaxt="n")
text((1:2+0.5),0.2,pos = 2, cex = 0.9,srt = 30, xpd = TRUE, 
     labels = c("Banded","OWN"))

boxplot(cbind(parametBAND[,2],parametOWN[,2]),main="residuo",xaxt="n")
text((1:2+0.5),0.2,pos = 2, cex = 0.9,srt = 30, xpd = TRUE, 
     labels = c("Banded","OWN"))

boxplot(cbind(parametBAND[,3],parametOWN[,3]),main="column correlation - 0.5",xaxt="n")
text((1:2+0.5),-0.1,pos = 2, cex = 0.9,srt = 30, xpd = TRUE, 
     labels = c("Banded","OWN"))
boxplot(cbind(parametBAND[,4],parametOWN[,4]),main="row correlation - (-0.3)",xaxt="n")
text((1:2+0.5),-0.7,pos = 2, cex = 0.9,srt = 30, xpd = TRUE, 
     labels = c("Banded","OWN"))
