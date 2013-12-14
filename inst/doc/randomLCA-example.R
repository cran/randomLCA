### R code from vignette source 'randomLCA-example.Rnw'

###################################################
### code chunk number 1: randomLCA-example.Rnw:13-21
###################################################
library(lattice)
library(xtable)
library(randomLCA)
data(dentistry)
load("dentistry.lca2.outcomes.boot.RData")
load("dentistry.lca2random.outcomes.boot.RData")
set.seed(12345)
options(width=60)


###################################################
### code chunk number 2: randomLCA-example.Rnw:54-60
###################################################
myocardial.lca1 <- randomLCA(myocardial[,1:4],
	freq=myocardial$freq,nclass=1)
myocardial.lca2 <- randomLCA(myocardial[,1:4],
	freq=myocardial$freq,nclass=2,calcSE=TRUE)
myocardial.lca3 <- randomLCA(myocardial[,1:4],
	freq=myocardial$freq,nclass=3)


###################################################
### code chunk number 3: randomLCA-example.Rnw:65-67
###################################################
bic.data <- data.frame(classes=1:3,bic=c(BIC(myocardial.lca1),
	BIC(myocardial.lca2),BIC(myocardial.lca3)))


###################################################
### code chunk number 4: randomLCA-example.Rnw:72-74
###################################################
print(xtable(bic.data, digits = c(0,0,1),caption="BIC by class.",
	label="tab:bic1"),include.rownames=FALSE)


###################################################
### code chunk number 5: randomLCA-example.Rnw:83-84
###################################################
summary(myocardial.lca2)


###################################################
### code chunk number 6: randomLCA-example.Rnw:88-89
###################################################
outcomep.data <- summary(myocardial.lca2)$outcomep


###################################################
### code chunk number 7: randomLCA-example.Rnw:94-96
###################################################
print(xtable(outcomep.data, digits = c(0,3,3,3,3),caption="Outcome Probabilities.",
	label="tab:outcomep1"))


###################################################
### code chunk number 8: randomLCA-example.Rnw:108-110
###################################################
classprobs <- cbind(myocardial.lca2$patterns,myocardial.lca2$classprob)
colnames(classprobs) <- c(names(myocardial)[1:4],"Class 1","Class 2")


###################################################
### code chunk number 9: randomLCA-example.Rnw:117-119
###################################################
print(xtable(classprobs, digits = c(0,0,0,0,0,3,3),caption="Class Probabilities.",
	label="tab:classprob1"),include.rownames=FALSE)


###################################################
### code chunk number 10: randomLCA-example.Rnw:127-130
###################################################
trellis.par.set(col.whitebg())
print(plot(myocardial.lca2,type="l",xlab="Test",ylab="Outcome Probability",
	scales=list(x=list(at=1:4,labels=names(myocardial)[1:4]))))


###################################################
### code chunk number 11: randomLCA-example.Rnw:142-150
###################################################
dentistry.lca1 <- randomLCA(dentistry[,1:5],
	freq=dentistry$freq,nclass=1)
dentistry.lca2 <- randomLCA(dentistry[,1:5],
	freq=dentistry$freq,nclass=2)
dentistry.lca3 <- randomLCA(dentistry[,1:5],
	freq=dentistry$freq,nclass=3,quadpoints=31)
dentistry.lca4 <- randomLCA(dentistry[,1:5],
	freq=dentistry$freq,nclass=4,quadpoints=41)


###################################################
### code chunk number 12: randomLCA-example.Rnw:155-157
###################################################
bic.data <- data.frame(classes=1:4,bic=c(BIC(dentistry.lca1),
	BIC(dentistry.lca2),BIC(dentistry.lca3),BIC(dentistry.lca4)))


###################################################
### code chunk number 13: randomLCA-example.Rnw:163-165
###################################################
print(xtable(bic.data, digits = c(0,0,1),caption="BIC by class.",
	label="tab:bic"),include.rownames=FALSE)


###################################################
### code chunk number 14: randomLCA-example.Rnw:172-174
###################################################
trellis.par.set(col.whitebg())
print(plot(dentistry.lca3,type="l",xlab="Dentist",ylab="Outcome Probability"))


###################################################
### code chunk number 15: randomLCA-example.Rnw:182-184
###################################################
trellis.par.set(col.whitebg())
print(plot(dentistry.lca2,type="l",xlab="Dentist",ylab="Outcome Probability"))


###################################################
### code chunk number 16: randomLCA-example.Rnw:193-194
###################################################
outcome.probs(dentistry.lca2)


###################################################
### code chunk number 17: randomLCA-example.Rnw:199-212
###################################################
probs <- outcome.probs(dentistry.lca2)
# this swaps around the probabilities based on the knowledge that the undiseased class is more common
order <- ifelse(dentistry.lca2$classp[2]>dentistry.lca2$classp[1],1,2)

spec <- NULL
sens <- NULL
for (i in 1:5) {
	sens <- c(sens,sprintf("%3.2f (%3.2f,%3.2f)",probs[[order]]$Outcome[i],probs[[order]]$"2.5 %"[i],probs[[order]]$"97.5 %"[i]))
	spec <- c(spec,sprintf("%3.2f (%3.2f,%3.2f)",1-probs[[3-order]]$Outcome[i],1-probs[[3-order]]$"97.5 %"[i],1-probs[[3-order]]$"2.5 %"[i]))
}
stable <- data.frame(sens,spec)
names(stable) <- c("Sensitivity","Specificity")
row.names(stable) <- paste("V",1:5,sep="")


###################################################
### code chunk number 18: randomLCA-example.Rnw:215-218
###################################################
print(xtable(stable, digits = c(0,2,2),
	caption="Sensitivity and Specificity",
	label="tab:outcomeconfint"),include.rownames=TRUE)


###################################################
### code chunk number 19: randomLCA-example.Rnw:231-244
###################################################
probs <- dentistry.lca2.outcomes.boot
# this swaps around the probabilities based on the knowledge that outcome probabilities are higher in the diseased class
order <- ifelse(probs[[1]]$Outcome[1]>probs[[2]]$Outcome[1],1,2)

spec <- NULL
sens <- NULL
for (i in 1:5) {
	sens <- c(sens,sprintf("%3.2f (%3.2f,%3.2f)",probs[[order]]$Outcome[i],probs[[order]]$"2.5 %"[i],probs[[order]]$"97.5 %"[i]))
	spec <- c(spec,sprintf("%3.2f (%3.2f,%3.2f)",1-probs[[3-order]]$Outcome[i],1-probs[[3-order]]$"97.5 %"[i],1-probs[[3-order]]$"2.5 %"[i]))
}
stable <- data.frame(sens,spec)
names(stable) <- c("Sensitivity","Specificity")
row.names(stable) <- paste("V",1:5,sep="")


###################################################
### code chunk number 20: randomLCA-example.Rnw:247-250
###################################################
print(xtable(stable, digits = c(0,2,2),
	caption="Sensitivity and Specificity",
	label="tab:outcomeconfintboot"),include.rownames=TRUE)


###################################################
### code chunk number 21: randomLCA-example.Rnw:256-260
###################################################
itpr <- ifelse(dentistry.lca2$classp[2]>dentistry.lca2$classp[1],1,2)
ifpr <- 3-itpr
probs <- outcome.probs(dentistry.lca2)
probs <- data.frame(tpr=probs[[itpr]][,1],fpr=probs[[ifpr]][,1])


###################################################
### code chunk number 22: randomLCA-example.Rnw:265-269
###################################################
trellis.par.set(col.whitebg())
print(plot(tpr~fpr,type="p",
	xlab="False Positive Rate\n(1-Specificity)",
	ylab="True Positive Rate (Sensitivity)",data=probs))


###################################################
### code chunk number 23: randomLCA-example.Rnw:300-302
###################################################
dentistry.lca2random <- randomLCA(dentistry[,1:5],freq=dentistry$freq,
	initmodel=dentistry.lca2,nclass=2,random=TRUE,quadpoints=41,probit=TRUE)


###################################################
### code chunk number 24: randomLCA-example.Rnw:307-310
###################################################
dentistry.lca2random1 <- randomLCA(dentistry[,1:5],freq=dentistry$freq,
	initmodel=dentistry.lca2random,nclass=2,random=TRUE,probit=TRUE,
	quadpoints=41,blocksize=5)


###################################################
### code chunk number 25: randomLCA-example.Rnw:317-320
###################################################
dentistry.lca2random2 <- randomLCA(dentistry[,1:5],freq=dentistry$freq,
	initmodel=dentistry.lca2random1,nclass=2,random=TRUE,probit=TRUE,
	blocksize=5,byclass=TRUE,quadpoints=41)


###################################################
### code chunk number 26: randomLCA-example.Rnw:331-333
###################################################
trellis.par.set(col.whitebg())
print(plot(dentistry.lca2random1,graphtype="marginal",type="l",xlab="Dentist",ylab="Marginal Outcome Probability"))


###################################################
### code chunk number 27: randomLCA-example.Rnw:349-362
###################################################
probs <-dentistry.lca2random.outcomes.boot
# this swaps around the probabilities based on the knowledge that outcome probabilities are higher in the diseased class
order <- ifelse(probs[[1]]$Outcome[1]>probs[[2]]$Outcome[1],1,2)

spec <- NULL
sens <- NULL
for (i in 1:5) {
	sens <- c(sens,sprintf("%3.2f (%3.2f,%3.2f)",probs[[order]]$Outcome[i],probs[[order]]$"2.5 %"[i],probs[[order]]$"97.5 %"[i]))
	spec <- c(spec,sprintf("%3.2f (%3.2f,%3.2f)",1-probs[[3-order]]$Outcome[i],1-probs[[3-order]]$"97.5 %"[i],1-probs[[3-order]]$"2.5 %"[i]))
}
stable <- data.frame(sens,spec)
names(stable) <- c("Sensitivity","Specificity")
row.names(stable) <- paste("V",1:5,sep="")


###################################################
### code chunk number 28: randomLCA-example.Rnw:366-369
###################################################
print(xtable(stable, digits = c(0,2,2),
	caption="Sensitivity and Specificity",
	label="tab:outcomeconfintboot2"),include.rownames=TRUE)


###################################################
### code chunk number 29: randomLCA-example.Rnw:374-377
###################################################
obs.data <- data.frame(dentistry.lca2random1$patterns,dentistry.lca2random1$observed,
	dentistry.lca2$fitted,dentistry.lca2random1$fitted)
names(obs.data) <- c("V1","V2","V3","V4","V5","Obs","Exp 2LC","Exp 2LCR")


###################################################
### code chunk number 30: randomLCA-example.Rnw:382-384
###################################################
print(xtable(obs.data, caption="Observed and expected frequencies",
	digits = c(0,rep(0,5),0,1,1),label="tab:obs"),include.rownames=FALSE)


