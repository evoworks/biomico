"BioMiCo" <- function(train, envs, unknownEnv=TRUE, rarefaction_depth=1000){
    train <- as.matrix(train)

    # enforce integer data
    if(sum(as.integer(train) != as.numeric(train)) > 0){
        stop('Data must be integral. Consider using "ceiling(datatable)" or ceiling(1000*datatable) to convert floating-point data to integers.')
    }
    
    envs.factors <- factor(unlist(envs))
    envs.list <- sort(unique(levels(unlist(envs.factors))))
    if ( "Unknown" %in% envs.list ){
      stop('One of the environments is labeled as Unknown. Consider renaming it to something else.')
    }
    if (unknownEnv)
      # add "Unknown" environment to every sample's list of environments
      for (sample.id in 1:length(envs) )
        envs[[sample.id]] = c( envs[[sample.id]], "Unknown")
    
    
    ## rarefy samples above maxdepth if requested
    #if(!is.null(rarefaction_depth)) train <- rarefy(train, rarefaction_depth)
    
    ret <- list(train=train, envs=envs)
    class(ret) <- "BioMiCo"
    return(invisible(ret))
}

"predict.BioMiCo" <- function(stobj, test=NULL, G=50, source.unknownEnv=TRUE, sink.unknownEnv=TRUE,
            burnin=100, nrestarts=2, ndraws.per.restart=2, delay=50,
            alpha.pi=1e-1, alpha.phi=1e-2, alpha.theta=1e-2, rarefaction_depth=1000,
            verbosity=1){

  
  library("BioMiCo")
  
  if(!is.null(test)){
    # if test is numeric, cast as a row matrix
    if(class(test) == "numeric" || class(test) == "integer"){
      test <- matrix(test, nrow=1)
    } else {
      test <- as.matrix(test)
    }
    if(sum(as.integer(test) != as.numeric(test)) > 0){
      stop('Data must be integral. Consider using "ceiling(datatable)" or ceiling(1000*datatable) to convert floating-point data to integers.')
    }
  
    train = stobj$train
    envs = stobj$envs
  
    N.train = nrow(train)# number of samples
    cat(paste("number of training samples: ", N.train, "\n", sep=""))
    
    T = ncol(train) # number of taxa
    cat(paste("number of taxa in training data: ", T, "\n", sep=""))  
    Taxa.list = colnames(train)
  
    envs.list <- unique(levels(factor(unlist(envs))))
    V = length(envs.list) # number of envs including the "Unknown" environment
    envs.list.idx = 0:(V-1)
    names(envs.list.idx) = envs.list
  
    cat(paste("number of envs including the Unknown environment: ", V, "\n", sep=""))
    cat( envs.list)
    cat(paste("\n total number of taxa in all samples:", sum(train), "\n", sep=""))

    envs.int = list()
    for (sample.id in 1:length(envs) ){
      envs.int[[sample.id]] = as.vector(envs.list.idx[ envs[[sample.id]] ])
    }

    source.G = G
    sink.G = G
    if(source.unknownEnv)
      source.G = source.G+1
    
    Z = list()
    X = list()
    draws.train= .Call("runGibbs", train, envs.int, source.G, V, T, N.train, source.unknownEnv, Z, X, burnin, nrestarts, ndraws.per.restart, delay, alpha.pi, alpha.phi, alpha.theta, rarefaction_depth, verbosity, 1, N.train, PACKAGE="BioMiCo")

    
    Count.T.G = list()
    Count.G.V = list()
    Count.S.V = list()
    alpha.phi.train = list()
    alpha.theta.mat = matrix(0, nrow=V, ncol=G)
    for (draw in 1:length(draws.train)){
      Count.T.G[[draw]] = matrix(0, nrow=T, ncol=source.G)
      Count.G.V[[draw]] = matrix(0, nrow=source.G, ncol=V)
      colnames(Count.G.V[[draw]]) = envs.list
      Count.S.V[[draw]] = matrix(0, nrow=N.train, ncol=V)
      colnames(Count.S.V[[draw]]) = envs.list
      alpha.theta.mat = alpha.theta.mat + draws.train[[draw]]$ALPHA_THETA_MAT
      alpha.phi.train[[draw]] = draws.train[[draw]]$ALPHA_PHI
      for (n in 1:N.train) # looping over all samples
        for (t in 1:T)# looping over all taxa
          if (length( draws.train[[draw]]$Z[[(n-1)*T + t]]) > 0)  # (train[n,t] > 0)
            for (i in 1:length( draws.train[[draw]]$Z[[(n-1)*T + t]]) ){  #(i in 1:train[n,t]){
              cur.G =  draws.train[[draw]]$Z[[(n-1)*T + t]][i] + 1 # indices in C++ are zero-based
              if (is.character(draws.train[[draw]]$X[[(n-1)*T + t]][i]) )
                cur.envs = draws.train[[draw]]$X[[(n-1)*T + t]][i]
              else
                cur.envs = draws.train[[draw]]$X[[(n-1)*T + t]][i] + 1 # indices in C++ are zero-based
              Count.T.G[[draw]][t, cur.G] = Count.T.G[[draw]][t, cur.G] +1
              Count.G.V[[draw]][cur.G, cur.envs] = Count.G.V[[draw]][cur.G, cur.envs] +1
              Count.S.V[[draw]][n, cur.envs] = Count.S.V[[draw]][n, cur.envs] +1
            }
    }
    alpha.theta.mat = alpha.theta.mat / length(draws.train)


    N.test = nrow(test) # number of test samples
    X.test.samples = matrix(0, nrow=N.test, ncol=V+1)
    rownames(X.test.samples) = rownames(test)
    for (test.sample in rownames(test)){
      cat(paste("predicting for sample ", test.sample, "\n", sep = ""))
      draws.test = .Call("runGibbsOnTest", as.matrix(t(test[test.sample,])), sink.G, V, T, 1, sink.unknownEnv, Count.T.G[[length(draws.train)]], Count.G.V[[length(draws.train)]], burnin, nrestarts, ndraws.per.restart, delay, alpha.pi, alpha.phi, alpha.theta, rarefaction_depth, verbosity, 1, 1, PACKAGE="BioMiCo")
      X.test = matrix(0, nrow=length(draws.test), ncol=V+1) #list()
      for (draw in 1:length(draws.test)){

        for (t in 1:T)# looping over all taxa
          if (length( draws.test[[draw]]$Z[[t]]) > 0)  
            for (i in 1:length( draws.test[[draw]]$X[[t]]) ){  
              if (is.character(draws.test[[draw]]$X[[t]][i]) )
                cur.envs = draws.test[[draw]]$X[[t]][i]
              else
                cur.envs = draws.test[[draw]]$X[[t]][i] + 1 
              X.test[draw, cur.envs] = X.test[draw, cur.envs] + 1
            }
      }
      X.test.samples[test.sample,] = X.test[length(draws.test),]

    
    }
    return(invisible( list(train.draws=draws.train, test.X = X.test.samples, alpha.theta.mat.train=alpha.theta.mat, alpha.phi.train=alpha.phi.train, Count.T.G= Count.T.G, Count.G.V = Count.G.V, Count.S.V = Count.S.V)  ))
    
    
  }
}


"train.BioMiCo" <- function(stobj, G=50, source.unknownEnv=TRUE, sink.unknownEnv=TRUE,
            burnin=100, nrestarts=2, ndraws.per.restart=2, delay=50,
            alpha.pi=1e-1, alpha.phi=1e-2, alpha.theta=1e-2, rarefaction_depth=1000,
            verbosity=1){

  
  library("BioMiCo")
  

  
  train = stobj$train
  envs = stobj$envs
  
  N.train = nrow(train)# number of samples
  cat(paste("number of training samples: ", N.train, "\n", sep=""))
  
  T = ncol(train) # number of taxa
  cat(paste("number of taxa in training data: ", T, "\n", sep=""))  
  Taxa.list = colnames(train)
  
  envs.list <- unique(levels(factor(unlist(envs))))
  V = length(envs.list) # number of envs including the "Unknown" environment
  envs.list.idx = 0:(V-1)
  names(envs.list.idx) = envs.list
  
  cat(paste("number of envs including the Unknown environment: ", V, "\n", sep=""))
  cat( envs.list)
  cat(paste("\n total number of taxa in all samples:", sum(train), "\n", sep=""))

  envs.int = list()
  for (sample.id in 1:length(envs) ){
    envs.int[[sample.id]] = as.vector(envs.list.idx[ envs[[sample.id]] ])
  }

  source.G = G
  sink.G = G
  if(source.unknownEnv)
    source.G = source.G+1
    
  Z = list()
  X = list()
  draws.train= .Call("runGibbs", train, envs.int, source.G, V, T, N.train, source.unknownEnv, Z, X, burnin, nrestarts, ndraws.per.restart, delay, alpha.pi, alpha.phi, alpha.theta, rarefaction_depth, verbosity, 1, N.train, PACKAGE="BioMiCo")

    
  Count.T.G = list()
  Count.G.V = list()
  Count.S.V = list()
  alpha.phi.train = list()
  alpha.theta.mat = matrix(0, nrow=V, ncol=G)
  for (draw in 1:length(draws.train)){
    Count.T.G[[draw]] = matrix(0, nrow=T, ncol=source.G)
    Count.G.V[[draw]] = matrix(0, nrow=source.G, ncol=V)
    colnames(Count.G.V[[draw]]) = envs.list
    Count.S.V[[draw]] = matrix(0, nrow=N.train, ncol=V)
    colnames(Count.S.V[[draw]]) = envs.list
    alpha.theta.mat = alpha.theta.mat + draws.train[[draw]]$ALPHA_THETA_MAT
    alpha.phi.train[[draw]] = draws.train[[draw]]$ALPHA_PHI
    for (n in 1:N.train) # looping over all samples
      for (t in 1:T)# looping over all taxa
        if (length( draws.train[[draw]]$Z[[(n-1)*T + t]]) > 0)  # (train[n,t] > 0)
          for (i in 1:length( draws.train[[draw]]$Z[[(n-1)*T + t]]) ){  #(i in 1:train[n,t]){
            cur.G =  draws.train[[draw]]$Z[[(n-1)*T + t]][i] + 1 # indices in C++ are zero-based
            if (is.character(draws.train[[draw]]$X[[(n-1)*T + t]][i]) )
              cur.envs = draws.train[[draw]]$X[[(n-1)*T + t]][i]
            else
              cur.envs = draws.train[[draw]]$X[[(n-1)*T + t]][i] + 1 # indices in C++ are zero-based
            Count.T.G[[draw]][t, cur.G] = Count.T.G[[draw]][t, cur.G] +1
            Count.G.V[[draw]][cur.G, cur.envs] = Count.G.V[[draw]][cur.G, cur.envs] +1
            Count.S.V[[draw]][n, cur.envs] = Count.S.V[[draw]][n, cur.envs] +1
          }
  }
  alpha.theta.mat = alpha.theta.mat / length(draws.train)

  return(invisible( list(V=V, T=T, train.draws=draws.train, alpha.theta.mat.train=alpha.theta.mat, alpha.phi.train=alpha.phi.train, Count.T.G= Count.T.G, Count.G.V = Count.G.V, Count.S.V = Count.S.V)  ))
  
}

"test.BioMiCo" <- function(stobj, draws.train, Count.T.G, Count.G.V, test=NULL, test.sample.names=NULL, G=50, source.unknownEnv=TRUE, sink.unknownEnv=TRUE,
            burnin=100, nrestarts=2, ndraws.per.restart=2, delay=50,
            alpha.pi=1e-1, alpha.phi=1e-2, alpha.theta=1e-2, rarefaction_depth=1000,
            verbosity=1){

  
  library("BioMiCo")
  
  if(!is.null(test)){
    # if test is numeric, cast as a row matrix
    if(class(test) == "numeric" || class(test) == "integer"){
      test <- matrix(test, nrow=1)
      rownames(test) = test.sample.names
    } else {
      test <- as.matrix(test)
    }
    if(sum(as.integer(test) != as.numeric(test)) > 0){
      stop('Data must be integral. Consider using "ceiling(datatable)" or ceiling(1000*datatable) to convert floating-point data to integers.')
    }
  
    train = stobj$train
    envs = stobj$envs
  
    N.train = nrow(train)# number of samples
    cat(paste("number of training samples: ", N.train, "\n", sep=""))
    
    T = ncol(train) # number of taxa
    cat(paste("number of taxa in training data: ", T, "\n", sep=""))  
    Taxa.list = colnames(train)
  
    envs.list <- unique(levels(factor(unlist(envs))))
    V = length(envs.list) # number of envs including the "Unknown" environment
    envs.list.idx = 0:(V-1)
    names(envs.list.idx) = envs.list
  
    cat(paste("number of envs including the Unknown environment: ", V, "\n", sep=""))
    cat( envs.list)
    cat(paste("\n total number of taxa in all samples:", sum(train), "\n", sep=""))

    envs.int = list()
    for (sample.id in 1:length(envs) ){
      envs.int[[sample.id]] = as.vector(envs.list.idx[ envs[[sample.id]] ])
    }

    source.G = G
    sink.G = G
    if(source.unknownEnv)
      source.G = source.G+1
    

    N.test = nrow(test) # number of test samples
    cat(paste("number of test samples for prediction: ", N.test, "\n", sep=""))
    X.test.samples = matrix(0, nrow=N.test, ncol=V+1) 
    colnames(X.test.samples)= c(envs.list,"unknown")
    rownames(X.test.samples) = rownames(test)
    for (test.sample in rownames(test)){
      cat(paste("predicting for sample ", test.sample, "\n", sep = ""))
      draws.test = .Call("runGibbsOnTest", as.matrix(t(test[test.sample,])), sink.G, V, T, 1, sink.unknownEnv, Count.T.G[[length(draws.train)]], Count.G.V[[length(draws.train)]], burnin, nrestarts, ndraws.per.restart, delay, alpha.pi, alpha.phi, alpha.theta, rarefaction_depth, verbosity, 1, 1, PACKAGE="BioMiCo")
      X.test = matrix(0, nrow=length(draws.test), ncol=V+1) #list()
      for (draw in 1:length(draws.test)){

        for (t in 1:T)# looping over all taxa
          if (length( draws.test[[draw]]$Z[[t]]) > 0)  
            for (i in 1:length( draws.test[[draw]]$X[[t]]) ){  
              if (is.character(draws.test[[draw]]$X[[t]][i]) )
                cur.envs = draws.test[[draw]]$X[[t]][i]
              else
                cur.envs = draws.test[[draw]]$X[[t]][i] + 1 
              X.test[draw, cur.envs] = X.test[draw, cur.envs] + 1
            }
      }
      X.test.samples[test.sample,] = X.test[length(draws.test),]

    
    }
    return(invisible( list(test.draws = draws.test, test.X = X.test.samples)  ))
    
    
  }
}
