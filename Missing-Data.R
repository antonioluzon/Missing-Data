#Lectura de los datos
library(foreign)
datos<-read.spss("datos.sav",to.data.frame = TRUE)


summary(is.na(datos))


library(miceRanger)
library(mice)
library(VIM)

#AMPUTACIÓN 1
amp1<-amputeData(datos,cols=c("P1","P2","P3","V1","P4"),perc=0.05)
summary(is.na(amp1))

#Proporción de datos faltantes
pMiss <- function(x){sum(is.na(x))/length(x)*100}
apply(amp1,2,pMiss)

#Gráfico para observar el númro de casos faltantes
md.pattern(amp1,rotate.names = TRUE)

#Gráfico para observar el número de casos faltantes
aggr_plot <- aggr(amp1, col=c('blue','red'), numbers=TRUE, 
                  sortVars=TRUE, labels=names(amp1), cex.axis=.7, 
                  gap=3, ylab=c("Histograma de datos faltantes","Pattern"))


imp1<-mice(amp1)

#Comprobación qué variables son utilizadas como predictoras
imp1$predictorMatrix

#Proporción de casos utilizables
p <- md.pairs(amp1)
round(p$mr/(p$mr + p$mm), 3)

#Variables predictoras
quickpred(amp1)

#IMPUTACIÓN 1 
imp1 <- mice(amp1, pred = quickpred(amp1, minpuc = 0.25),
             meth=c("","","","logreg","polyreg",
                    "logreg","polr","logreg"),m=5)

#Convergencia MICE
plot(imp1,c("P1","P2","P3","V1","P4"))

#Puntaje de precisión
V1.na <- is.na(amp1$V1)
fit.V1 <- with(imp1, glm(V1.na ~ P1 + P2 + P3 + Tit))
ps4 <- rep(rowMeans(sapply(fit.V1$analyses, fitted.values)), 6)
xyplot(imp1, V1 ~ ps4 | .imp, pch = c(1, 20), cex = c(0.8, 1.2))

#Completar la imputación y comprobación de casos faltantes en los datos imputados 1
comp1<-complete(imp1)
summary(is.na(comp1))

#AMPUTACIÓN 2
amp2<-amputeData(datos,cols=c("P1","P2","P3","V1","P4"),perc=0.1)
summary(is.na(amp2))

#Porcentaje de datos faltantes
apply(amp2,2,pMiss)

#Gráfico para observar los datos faltantes en la amputación 2
md.pattern(amp2,rotate.names = TRUE)

#Gráfico para observar los datos faltantes en la amputación 2
aggr_plot2 <- aggr(amp2, col=c('blue','red'), numbers=TRUE, 
                   sortVars=TRUE, labels=names(amp2), cex.axis=.7, 
                   gap=3, ylab=c("Histograma de datos faltantes","Pattern"))


#Variables que actuan como predictoras
quickpred(amp2)


#IMPUTACIÓN 2
imp2 <- mice(amp2, pred = quickpred(amp2, minpuc = 0.25),
             meth=c("","","","logreg","polyreg",
                    "logreg","polr","logreg"),m=5)

#Convergencia
plot(imp2,c("P1","P2","P3","V1"))

#Complerar datos imputados 2 y observar casos faltantes
comp2<-complete(imp2)
summary(is.na(comp2))



##COMPARACIÓN DATOS COMPLETOS CON IMPUTADOS MEDIANTE ANALISIS DE DATOS CUALITATIVO


## P1~SEXO

#completos
with(datos,table(Sexo,P1))

round(chisq.test(with(datos,table(Sexo,P1)))$expected,3)

#Test de independencia chi-cuadrado
chisq.test(with(datos,table(Sexo,P1)))

#Diferencia de proporciones
dif_pro<-function(x,alpha=0.05){
  dif_p=x[1,1]/sum(x[1,]-x[2,1]/sum(x[2,]))
  sigma_difp=sqrt(x[1,1]*x[1,2]/sum(x[1,])^3+
                    x[2,1]*x[2,2]/sum(x[2,])^3)
  inf=dif_p-qnorm(1-alpha/2)* sigma_difp
  sup=dif_p+qnorm(1-alpha/2)*sigma_difp
  return(list(dif_p=dif_p,sigma_p=sigma_difp,IC_P=c(inf,sup)))
}

dif_pro(with(datos,table(Sexo,P1)))

#Riesgo Relativo
riesgo<-function(x,alpha=0.05){
  riesgo=(x[1,1]/sum(x[1,]))/(x[2,1]/sum(x[2,]))
  return(list(riesgo=riesgo))
}
riesgo(with(datos,table(Sexo,P1)))


#Cociente de ventajas
odds.ratio<-function(x,alpha=0.05){
  theta<-x[1,1]*x[2,2]/(x[1,2]*x[2,1])
  sigma<-sqrt(sum(1/x))
  Za<-qnorm(1-alpha/2)
  inf<-exp(log(theta)-Za*sigma)
  sup<-exp(log(theta)+Za*sigma)
  IC<-c(inf,sup)
  return(list(theta=theta,sigma=sigma,IC=IC))
  
odds.ratio(with(datos,table(P1,Sexo)))

#Gráfico cociente de ventajas
fourfoldplot(with(datos,table(P1,Sexo)))

## Imputados 1
with(comp1,table(Sexo,P1))

round(chisq.test(with(comp1,table(Sexo,P1)))$expected,3)

chisq.test(with(comp1,table(Sexo,P1)))

dif_pro(with(comp1,table(Sexo,P1)))

riesgo(with(comp1,table(Sexo,P1)))

odds.ratio(with(comp1,table(P1,Sexo)))

fourfoldplot(with(comp1,table(P1,Sexo)))


## Imputados 2
with(comp2,table(Sexo,P1))

round(chisq.test(with(comp2,table(Sexo,P1)))$expected,3)

chisq.test(with(comp2,table(Sexo,P1)))

dif_pro(with(comp2,table(Sexo,P1)))

riesgo(with(comp2,table(Sexo,P1)))

odds.ratio(with(comp2,table(P1,Sexo)))

fourfoldplot(with(comp2,table(P1,Sexo)))



##P2~P3

with(datos,table(P2,P3))

round(chisq.test(with(datos,table(P2,P3)))$expected,2)

#Test de independencia chi-cuadrado
chisq.test(with(datos,table(P2,P3)))

library(vcd)

#Residuos estandarizados
residuos.S<-chisq.test(with(datos,table(P2,P3)))$stdres
mosaic(with(datos,table(P2,P3)),residuals=residuos.S,
       residuals_type='residuos est.',
       gp=shading_hcl,labeling=labeling_residuals, rotate.names=TRUE)


#V de Cramer
library(vcd)
assocstats(with(datos,table(P2,P3)))


## Imputados 1
with(comp1,table(P2,P3))

round(chisq.test(with(comp1,table(P2,P3)))$expected,2)

chisq.test(with(comp1,table(P2,P3)))

residuos.S<-chisq.test(with(comp1,table(P2,P3)))$stdres

mosaic(with(comp1,table(P2,P3)),residuals=residuos.S,
       residuals_type='residuos est.',
       gp=shading_hcl,labeling=labeling_residuals, rotate.names=TRUE)

assocstats(with(comp1,table(P2,P3)))


## Imputados 2
with(comp2,table(P2,P3))

round(chisq.test(with(comp2,table(P2,P3)))$expected,2)

chisq.test(with(comp2,table(P2,P3)))

residuos.S<-chisq.test(with(comp2,table(P2,P3)))$stdres

mosaic(with(comp2,table(P2,P3)),residuals=residuos.S,
       residuals_type='residuos est.',
       gp=shading_hcl,labeling=labeling_residuals, rotate.names=TRUE)

assocstats(with(com2,table(P2,P3)))


## V1~P4

with(datos,table(V1,P4))
round(chisq.test(with(datos,table(V1,P4)))$expected,3)

#Test de independencia chi-cuadrado

chisq.test(with(datos,table(V1,P2)))

#Residuos estandarizados
residuos.S<-chisq.test(with(datos,table(V1,P4)))$stdres
mosaic(with(datos,table(V1,P4)),residuals=residuos.S,
       residuals_type='residuos est.',
       gp=shading_hcl,labeling=labeling_residuals, rotate.names=TRUE)

#Test de tendencia lineal
orig<-with(datos,table(V1,P4))
indorig<-chisq.test(with(datos,table(V1,P4)))
yj<-1:2
xi<-1:11
xy.ij<-expand.grid(xi,yj)
nij<-as.numeric(orig)
orig.ord<-data.frame(x=xy.ij[,1],y=xy.ij[,2],nij=nij)
res.lm<-lm(x~y,data=orig.ord,weights=nij)
res.lm
R2<-summary(res.lm)$r.squared
R2
X2.I<-as.numeric(indorig$statistic)
X2.I
R2<-summary(res.lm)$r.squared
R2
X.I.L<-sum(nij)*R2
X.I.L
X2.L<-X2.I-X.I.L
X2.L


#V de Cramer
assocstats(with(datos,table(V1,P4)))


## Imputados 1
with(comp1,table(V1,P4))

round(chisq.test(with(comp1,table(V1,P4)))$expected,3)

chisq.test(with(comp1,table(V1,P2)))

residuos.S<-chisq.test(with(comp1,table(V1,P4)))$stdres

mosaic(with(comp1,table(V1,P4)),residuals=residuos.S,
       residuals_type='residuos est.',
       gp=shading_hcl,labeling=labeling_residuals, rotate.names=TRUE)

orig<-with(comp1,table(V1,P4))
indorig<-chisq.test(with(comp1,table(V1,P4)))
yj<-1:2
xi<-1:11
xy.ij<-expand.grid(xi,yj)
nij<-as.numeric(orig)
orig.ord<-data.frame(x=xy.ij[,1],y=xy.ij[,2],nij=nij)
res.lm<-lm(x~y,data=orig.ord,weights=nij)
res.lm
R2<-summary(res.lm)$r.squared
R2
X2.I<-as.numeric(indorig$statistic)
X2.I
R2<-summary(res.lm)$r.squared
R2
X.I.L<-sum(nij)*R2
X.I.L
X2.L<-X2.I-X.I.L
X2.L

assocstats(with(comp1,table(V1,P4)))

## Imputados 2
with(comp2,table(V1,P4))

round(chisq.test(with(comp2,table(V1,P4)))$expected,3)

chisq.test(with(comp2,table(V1,P2)))

residuos.S<-chisq.test(with(comp2,table(V1,P4)))$stdres

mosaic(with(comp2,table(V1,P4)),residuals=residuos.S,
       residuals_type='residuos est.',
       gp=shading_hcl,labeling=labeling_residuals, rotate.names=TRUE)

orig<-with(comp2,table(V1,P4))
indorig<-chisq.test(with(comp2,table(V1,P4)))
yj<-1:2
xi<-1:11
xy.ij<-expand.grid(xi,yj)
nij<-as.numeric(orig)
orig.ord<-data.frame(x=xy.ij[,1],y=xy.ij[,2],nij=nij)
res.lm<-lm(x~y,data=orig.ord,weights=nij)
res.lm
R2<-summary(res.lm)$r.squared
R2
X2.I<-as.numeric(indorig$statistic)
X2.I
R2<-summary(res.lm)$r.squared
R2
X.I.L<-sum(nij)*R2
X.I.L
X2.L<-X2.I-X.I.L
X2.L

assocstats(with(comp2,table(V1,P4)))

