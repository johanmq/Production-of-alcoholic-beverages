#rm(list=ls())

#install.packages('tseries')
library(tseries)

#install.packages("astsa")
library(astsa)

#install.packages("stargazer")
library(stargazer)


data = read.csv('valeurs_mensuelles.csv',sep=';')
data = data[349:1,]
plot(x=seq(1990,2019,1/12), y=data$indice, type='l', ylab='Indice de production',xlab='Annee')

#### Question 2

serie <- ts(data$indice, start=0, end=336)
serieD <- (serie - lag(serie,-12)) # On enleve la saisonnalite D=1


a = acf(serie, lag.max = 36)


res_acf <- acf(serieD)

adf.test(serieD) #L'hypothese de non-stat est rejetee a 1%, pas de racine unitaire
pp.test(serieD) #Mieux que DF car par d'hyp sur les residus, H0 rejetee -> stationnaire
kpss.test(serieD, null='Trend') # PAS STATIONNAIRE : H0 = stationnaire
kpss.test(serieD, null='Level') # PAS STATIONNAIRE : H0 = stationnaire

serieD2 = (serieD - lag(serieD,-1)) #on rajoute une difference : d=1, toujours D=1

adf.test(serieD2) #L'hypothese de non-stat est rejetee a 1%, pas de racine unitaire
pp.test(serieD2) #Mieux que DF car par d'hyp sur les residus, H0 rejetee -> stationnaire
kpss.test(serieD2, null='Trend') # STATIONNAIRE : H0 = stationnaire
kpss.test(serieD2, null='Level') # STATIONNAIRE : H0 = stationnaire

###### Question 3 ######

plot(serie, type='l',ylab='Xt')
plot(serieD2, type="l",ylab='Yt') #Semble stationnaire


##### Question 4 ######

res_acf <- acf(serieD2,main='')
# --> Q*=1, q*=1
res_pacf <- pacf(serieD2, lag.max=72, main='')
# --> P*=4, p*=8

# La fonctino Qtest renvoie les p-val des Ljung-Box test (Portmanteau) pour un lag maximal k
Qtest <- function(series, k) {
  t(apply(matrix(1:k), 1, FUN=function(l) {
    pval <- Box.test(series, lag=l, type="Ljung-Box")$p.value
    return(c("lag"=l,"pval"=pval))
  }))
}

signif <- function(estim){ #fonction de test des significations individuelles des coefficients des ordres max
  coef <- estim$fit$coef
  se <- sqrt(diag(estim$fit$var.coef))
  t <- coef/se
  pval <- (1-pnorm(abs(t)))*2
  return(rbind(coef,se,pval))
}

# Les sous-modeles possibles sont tous les SARIMA{12}(P,1,Q)(p,1,q) tels que P<=4 et Q<=1,
# p<=8 et q<=1.
# On choisira un modele siginificatif (les coefficients associes aux ordres les
# plus eleves de retards sont individuellement siginificatives) et les residus
# non autocorreles

Pmax = 4
pmax = 8
Qmax = 1
qmax = 1
lag = 36
m2=c()
i=c() #seulement pour suivre l'avancement

for (P in seq(0,Pmax)){
  for (Q in seq(0,Qmax)){
    for (p in seq(0,pmax)) {
      for (q in seq(0,qmax)) {
        i=rbind(i,c(P,Q,p,q))
        modele=try(sarima(serie,p,1,q,P,1,Q,12,no.constant=FALSE, details = FALSE))
        if(class(modele) == 'try-error'){next}
        if (length(which(Qtest(modele$fit$residuals,lag)[,2]>0.05)) == lag ){ # si le bruit est non correle on garde
          if(length( which( c(signif(modele)[3,p],signif(modele)[3,p+q],signif(modele)[3,p+q+P],signif(modele)[3,p+q+P+Q]) < 0.05) ) == 4){
            m2=rbind(m2,c(P,Q,p,q))
          }
        }
      }
    }
  }
}

stargazer(m2)
#m2
#      [,1] [,2] [,3] [,4]
#[1,]    0    1    8    0
#[2,]    0    1    8    1
#[3,]    2    0    8    0
#[4,]    3    0    8    0
#[5,]    4    0    8    0
#[6,]    4    0    8    1


##### Selectionner le meilleur avec min AIC, BIC

AIC2 = c()
BIC2 = c()

for ( i in seq(1,dim(m2)[1]) ){
  P = m2[i,1]
  Q = m2[i,2]
  p = m2[i,3]
  q = m2[i,4]
  #modele=try(arima(serie, c(p,0,q), seasonal = list(order=c(P,1,Q),period=12)))#, method='CSS')
  modele = try(sarima(serie,p,1,q,P,1,Q,12,no.constant=FALSE, details=FALSE))
  if ( class(modele)=="try-error"){next}#si converge pas
  else{
    #AIC = c(AIC,AIC(modele))
    AIC2 = c(AIC2,modele$AIC)
    BIC2 = c(BIC2,modele$BIC)#c(BIC,BIC(modele))
  }
}

which(AIC2==min(AIC2))
which(BIC2==min(BIC2))
m2[6,] # P=4, Q=0, p=8, q=1 #min AIC
m2[2,] # P=0, Q=1, p=8, q=1 #min BIC

est4 = sarima(serie,8,1,1,4,1,0,12,no.constant=FALSE, details = FALSE) #min AIC
est42 = sarima(serie,8,1,1,0,1,1,12,no.constant=FALSE, details = FALSE) #min BIC


##### Prediction de la derniere annee de production

#pred4
pred4 = sarima.for(serie,n.ahead=10,8,1,1,4,1,0,S=12) #Min AIC
ypred4 = pred4$pred

inf4 = ypred4 -1.96*pred4$se #int conf 95%
sup4 = ypred4 +1.96*pred4$se

plot(x=seq(2017+1/12,2019,1/12),data$indice[326:349],col='blue',type='o',ylim=c(50,150), ylab='Index',xlab='Date (year)') # bleu = vraie serie
lines(x=seq(2017+1/12,2019,1/12), c(rep('',12),ypred4,rep('',2)),type="o",col='red',lwd=2)
lines(x=seq(2017+1/12,2019,1/12),c(rep('',12),as.vector(inf4),rep('',2)),col='green', type='c',lwd=2)
lines(x=seq(2017+1/12,2019,1/12),c(rep('',12),as.vector(sup4),rep('',2)),col='green', type='c',lwd=2)
legend('topleft', legend=c("Index","Confidence interval (95%)","Forecasting"),
       col=c("blue","green","red"),lty=1:2, cex=0.8)
sqrt(mean((data$indice[338:347]-ypred4)**2))


#pred42
pred42 = sarima.for(serie,n.ahead=10,8,1,1,0,1,1,S=12) #Min BIC
ypred42 = pred42$pred

inf42 = ypred42 -1.96*pred42$se #int conf 95%
sup42 = ypred42 +1.96*pred42$se

plot(x=seq(2017+1/12,2019,1/12),data$indice[326:349],col='blue',type='o', ylim=c(50,150), ylab='Index',xlab='Date (year)') # bleu = vraie serie
lines(x=seq(2017+1/12,2019,1/12), c(rep('',12),ypred42,rep('',2)),type="o",col='red',lwd=2)
lines(x=seq(2017+1/12,2019,1/12),c(rep('',12),as.vector(inf42),rep('',2)),col='green', type='c',lwd=2)
lines(x=seq(2017+1/12,2019,1/12),c(rep('',12),as.vector(sup42),rep('',2)),col='green', type='c',lwd=2)
legend('topleft', legend=c("Index","Confidence interval (95%)","Forecasting"),
       col=c("blue","green","red"),lty=1:2, cex=0.8)

sqrt(mean((data$indice[338:347]-ypred42)**2))


## Question 6
# limite : normalite residus

qqnorm(est4$fit$residuals, main='Q-Q Plot pour le modele SARIMA{12}[(8,1,1),(4,1,0)]')
qqline(est4$fit$residuals) #residus normaux

jarque.bera.test(est4$fit$residuals)

## Question 7

serieQ7 <- ts(data$indice, start=0, end=346)
predQ7 = sarima.for(serieQ7,n.ahead=2,8,1,1,4,1,0,S=12)

ypred7 = predQ7$pred

inf7 = ypred7 -1.96*predQ7$se #int conf 95%
sup7 = ypred7 +1.96*predQ7$se

plot(x=seq(2018+6/12,2019+1/12,1/12),data$indice[342:349],col='blue',type='o', ylim=c(50,150), ylab='Index',xlab='Date (year)') # bleu = vraie serie
lines(x=seq(2018+6/12,2019+1/12,1/12), c(rep('',6),ypred7),type="o",col='red',lwd=2)
lines(x=seq(2018+6/12,2019+1/12,1/12),c(rep('',6),as.vector(inf7)),col='green', type='c')
lines(x=seq(2018+6/12,2019+1/12,1/12),c(rep('',6),as.vector(sup7)),col='green', type='c')
legend('topleft', legend=c("Index","Confidence interval (95%)","Forecasting"),
       col=c("blue","green","red"),lty=1:2, cex=0.8)

sqrt(mean((data$indice[348:349]-ypred42)**2))

