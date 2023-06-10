#setwd("/BioMiCo/")

TRAIN.Mat = read.delim(file="train.ix", header=TRUE, sep="\t", row.names=1)
TEST.test = read.delim(file="test.ix", header=TRUE, sep="\t", row.names=1)
ENVS=read.delim(file="env", header=FALSE, sep="\t", row.names=1)

envs=list()
for (i in rownames(ENVS)) envs[[i]]=as.character(ENVS[i,])
# num.communities can be changed  In the original paper we used 100 but we have since determined that 25 is sufficiently big
num.communities = 25

# load package BioMiCo
source('BioMiCoScripts.R')  #this must point to the BioMiCoScripts.R file if not running in this folder


st <- BioMiCo(TRAIN.Mat, envs, unknownEnv=FALSE)

# Estimate source proportions in test data
source.unknownEnv=FALSE
train.results = train.BioMiCo(st, G=num.communities, source.unknownEnv=FALSE, burnin=2500, nrestarts=1, ndraws.per.restart=20, delay=2000, alpha.pi=1, alpha.phi=1, alpha.theta=1, rarefaction_depth=1000, verbosity=1)

V = train.results$V
T = train.results$T
Count.S.V = train.results$Count.S.V
Count.T.G = train.results$Count.T.G
Count.G.V = train.results$Count.G.V
train.draws = train.results$train.draws
alpha.phi.train = train.results$alpha.phi.train
alpha.theta.mat.train = train.results$alpha.theta.mat.train

save(file="trainImage.RData", num.communities, source.unknownEnv, V, T, Count.S.V, Count.T.G, Count.G.V, train.draws, alpha.phi.train, alpha.theta.mat.train, TRAIN.Mat, TEST.test, envs)
#load(file="trainImage.RData")

# Estimate source proportions in test data

test.results = test.BioMiCo(st, train.draws, Count.T.G, Count.G.V, test=TEST.test, test.sample.names=c(test.sample), G=num.communities, source.unknownEnv=FALSE, sink.unknownEnv=FALSE, burnin=100, nrestarts=1, ndraws.per.restart=20, delay=50, alpha.pi=1, alpha.phi=1, alpha.theta=1, rarefaction_depth=1000, verbosity=1)

test.draws = test.results$test.draws
test.X = test.results$test.X


save(file="predictionImage.RData", num.communities, V, T, test.draws, test.X, envs, Count.S.V, Count.T.G, Count.G.V)
