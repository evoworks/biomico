#setwd("/BioMiCo/")
load("predictionImage.RData")

predictions=round(test.X/rowSums(test.X),3)
write.table(predictions, file="predictions", append=FALSE, quote=TRUE, sep="\t", eol="\n", row.names=TRUE, col.names=TRUE)

est.phi=t(Count.T.G[[20]])  
for(i in 1:25) est.phi[i,]=(est.phi[i,]+0.001)/sum(est.phi[i,]+0.001)
sum(est.phi[1,])  

tr_est.phi=t(est.phi)
write.table(tr_est.phi, file="OTU_pp", append=FALSE, quote=TRUE, sep="\t", eol="\n", row.names=TRUE, col.names=TRUE)

est.rho=t(Count.G.V[[20]])  
for(i in 1:V) est.rho[i,]=(est.rho[i,]+0.001)/sum(est.rho[i,]+0.001)
sum(est.rho[1,]) 

tr_est.rho=t(est.rho)
write.table(tr_est.rho, file="Assemblage_pp", append=FALSE, quote=TRUE, sep="\t", eol="\n", row.names=TRUE, col.names=TRUE)

