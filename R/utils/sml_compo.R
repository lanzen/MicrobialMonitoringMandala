####################################################################################
## Supervised machine learning utility classes
## Tristan Cordier (2018), further developed by Anders Lanzen (2019)
####################################################################################

# we permute the farm to be tested, so that we train a model on all the others farms and test performance. 

sml_compo <- function(otu_table, metadata, index, algo, 
                      optim_overfit = F, classif = F, cross_val = "Locality",
                      numTrees = 300) {

  ## we have to assign everything that ranger will use to the global environment (data.frame, intermediate object like mod...)
  ## even if something is created just before passing it to the function, it won't be found by ranger 
  ## maybe test creating an environment only for "inside" the sml_compo function
  comp <- metadata
  assign("otu_table", otu_table, envir = .GlobalEnv)
  assign("comp", comp, envir = .GlobalEnv)
  assign("index", index, envir = .GlobalEnv)

  ## library
  require(ranger)
  #require(iRF)
  #require(caret)
  #require(kohonen)
  #require(mxnet)
  #require(xgboost)
  require(e1071)
  require(iwrlars)
  #require(neuralnet)
  
  ## random forest and support vector machines works, the others are to be done.. 
  
  OTU <- data.frame(otu_table) 
  COMP <- comp 
  
  # regression by default
  classification <- classif
  
  if (classification) minNode <- 1
  if (!classification) minNode <- 5
  
  ## test of right length between data and compo
  if (dim(OTU)[[1]] != length(COMP[,index])) {
    print("Error: No same size between datasets, check it out...")
    return()
  }
  
  if (numTrees != 300 & optim_overfit)   print("Warning: number of trees cannot (and won't) be specified with optimal overfit")
  

  cross <- unique(COMP[,cross_val]) 
  if (length(cross) == 0) message("the permutation is not set...")
  
  ## ugly vectors for output
  combined1_rf <- c()
  combined2_rf <- c()
  combined1_sm <- c()
  combined2_sm <- c()
  combined1_nn <- c()
  combined2_nn <- c()
  combined1_lm <- c()
  combined2_lm <- c()
  combined1_dn <- c()
  combined2_dn <- c()
  combined1_sv <- c()
  combined2_sv <- c()
  combined1_xg <- c()
  combined2_xg <- c()
  combined1_la <- c()
  combined2_la <- c()
  
  ## gather RMSE for optim and check overfit
  over <- array(NA, c(2,10))
  rownames(over) <- c("training RMSE", "testing RMSE")
  colnames(over) <- c(10,50,100,150,200,250,300,350,400,500)
  
  farm_nam <- c()
  cpt <- 1
  
  for (i in cross)
  {
    # subsetting
    assign("nam_t",paste("test_farm_m", i,sep="_"), envir = .GlobalEnv)
    assign("nam_comp_t", paste("test_farm_c_m", i,sep="_"), envir = .GlobalEnv)
    assign(nam_t, subset(OTU, COMP[,cross_val] ==i), envir = .GlobalEnv)
    assign(nam_comp_t, subset(COMP, COMP[,cross_val] ==i), envir = .GlobalEnv)
    
    assign("nam", paste("train_farm_m", i,sep="_"), envir = .GlobalEnv)
    assign("nam_comp", paste("train_farm_c_m", i,sep="_"), envir = .GlobalEnv)
    assign(nam, subset(OTU, COMP[,cross_val] !=i), envir = .GlobalEnv)
    assign(nam_comp, subset(COMP, COMP[,cross_val] !=i), envir = .GlobalEnv)
    
    
    
    ## if random forest 
    if (algo == "RF")
    {
      ### if optim and check for overfit 
      if (optim_overfit == T) 
      {
        ## gather RMSE for optim and check overfit
        over <- array(NA, c(2,10))
        rownames(over) <- c("training RMSE", "testing RMSE")
        colnames(over) <- c(10,50,100,150,200,250,300,350,400,500)
        for (opt in colnames(over))
        {
          mod <- paste("RF_farm_m", i,sep="_")
          set.seed(1)
          assign(mod,ranger(get(nam_comp)[,index] ~ ., 
                            data=get(nam), mtry=floor(dim(get(nam))[2]/3), 
                            classification = classification, 
                            num.trees = as.numeric(opt), 
                            importance= "impurity", write.forest = T))
          ## prediction for new data
          predict_tr_rf <- predict(get(mod), get(nam_t))
          ## paste in over the RMSE values
          over["training RMSE",opt] <- get(mod)$prediction.error
          over["testing RMSE",opt] <- sqrt(mean((get(nam_comp_t)[,index] - predict_tr_rf$predictions)^2))
        }
        ### then get the minimum testing RMSE to get the optim num.tree param
        best_over <- as.numeric(attributes(which.min(over["training RMSE",]))$names)
        print("num.tree effect on testing RMSE")
        print(over)
        print(paste("Using: ", best_over))
        # maybe export a plot and curves for overfitting check...
        mod <- paste("RF_farm_m", i,sep="_")
        set.seed(1)
        assign(mod,ranger(get(nam_comp)[,index] ~ ., data=get(nam), 
                          mtry=floor(dim(get(nam))[2]/3), classification = classification, 
                          num.trees = best_over, importance= "impurity", write.forest = T))
        ## prediction for new data
        predict_tr_rf <- predict(get(mod), get(nam_t))
        combined1_rf <- c(combined1_rf,predict_tr_rf$predictions)
        combined2_rf <- c(combined2_rf,get(nam_comp_t)[,index])
      } else 
      {
        # ## now the fitting, random forest with Ranger package with default mtry (1/3) for regression according to Breiman
        # similar here the mod objet has to be in the global env. to be fetchable afterwards
        assign("mod", paste("RF_farm_m", i,sep="_"), envir = .GlobalEnv)
        
        ## removing empty OTUs
        OTUtr <- get(nam)[,colSums(get(nam)) != 0]
        # remove those ones for the testing dataset also
        OTUte <- get(nam_t)[,colSums(get(nam)) != 0]
        
        set.seed(1)
        assign(mod, ranger(get(nam_comp)[,index] ~ ., data=OTUtr, mtry=floor(dim(OTUtr)[2]/3), 
                           classification = classification, num.trees = numTrees, 
                           importance= "impurity", write.forest = T, min.node.size = minNode),envir = .GlobalEnv)
        
        #refer <- get(nam_comp)[,index]
        #mod <- ranger(refer ~ ., data=get(nam), mtry=floor(dim(get(nam))[2]/3), classification = classification, num.trees = 300, importance= "impurity", write.forest = T)
        #mod <- ranger(COMP[,"AMBI"] ~ ., data=OTU, mtry=floor(dim(OTU)[2]/3), classification = classification, num.trees = 300, importance= "impurity", write.forest = T)
        #assign(mod, ranger(refer ~ ., data=get(nam), mtry=floor(dim(get(nam))[2]/3), classification = classification, num.trees = 300, importance= "impurity", write.forest = T),envir = .GlobalEnv)
        #print("just after ranger")
        
        ## prediction for new data
        predict_tr_rf <- predict(get(mod), OTUte)
        combined1_rf <- c(combined1_rf,predict_tr_rf$predictions)
        combined2_rf <- c(combined2_rf,get(nam_comp_t)[,index])
      }
    }
    
    ## if random forest with probability of split vector (very experimental)
    if (algo == "RFProb")
    {
      ### if optim and check for overfit 
      if (optim_overfit == T) 
      {
        ## gather RMSE for optim and check overfit
        over <- array(NA, c(2,10))
        rownames(over) <- c("training RMSE", "testing RMSE")
        colnames(over) <- c(10,50,100,150,200,250,300,350,400,500)
        for (opt in colnames(over))
        {
          mod <- paste("RF_farm_m", i,sep="_")
          set.seed(1)
          assign(mod,ranger(get(nam_comp)[,index] ~ ., data=get(nam), mtry=floor(dim(get(nam))[2]/3), 
                            splitrule="extratrees",split.select.weights = splitProb,
                            classification = classification, num.trees = as.numeric(opt), 
                            importance= "impurity", write.forest = T))
          
          ## prediction for new data
          predict_tr_rf <- predict(get(mod), get(nam_t))
          ## paste in over the RMSE values
          over["training RMSE",opt] <- get(mod)$prediction.error
          over["testing RMSE",opt] <- sqrt(mean((get(nam_comp_t)[,index] - predict_tr_rf$predictions)^2))
        }
        ### then get the minimum testing RMSE to get the optim num.tree param
        best_over <- as.numeric(attributes(which.min(over["training RMSE",]))$names)
        print("num.tree effect on testing RMSE")
        print(over)
        print(paste("Using: ", best_over))
        # maybe export a plot and curves for overfitting check...
        mod <- paste("RF_farm_m", i,sep="_")
        set.seed(1)
        assign(mod,ranger(get(nam_comp)[,index] ~ ., data=get(nam), 
                          mtry=floor(dim(get(nam))[2]/3), 
                          splitrule="extratrees",split.select.weights = splitProb,
                          classification = classification, num.trees = best_over, 
                          importance= "impurity", write.forest = T))
        ## prediction for new data
        predict_tr_rf <- predict(get(mod), get(nam_t))
        combined1_rf <- c(combined1_rf,predict_tr_rf$predictions)
        combined2_rf <- c(combined2_rf,get(nam_comp_t)[,index])
      } else 
      {
        # ## now the fitting, random forest with Ranger package with default mtry (1/3) for regression according to Breiman
        # similar here the mod objet has to be in the global env. to be fetchable afterwards
        assign("mod", paste("RF_farm_m", i,sep="_"), envir = .GlobalEnv)
        
        ## removing empty OTUs
        OTUtr <- get(nam)[,colSums(get(nam)) != 0]
        # remove those ones for the testing dataset also
        OTUte <- get(nam_t)[,colSums(get(nam)) != 0]
        
        set.seed(1)
        notus = dim(OTUtr)[2]
        splitProb = c(rep(1/notus,notus-1),1)
        
        assign(mod, ranger(get(nam_comp)[,index] ~ ., data=OTUtr, 
                           mtry=floor(dim(OTUtr)[2]/3), 
                           splitrule="extratrees",split.select.weights = splitProb,
                           classification = classification, num.trees = 300, 
                           importance= "impurity", write.forest = T, 
                           min.node.size = minNode),envir = .GlobalEnv)
        
        #refer <- get(nam_comp)[,index]
        #mod <- ranger(refer ~ ., data=get(nam), mtry=floor(dim(get(nam))[2]/3), classification = classification, num.trees = 300, importance= "impurity", write.forest = T)
        #mod <- ranger(COMP[,"AMBI"] ~ ., data=OTU, mtry=floor(dim(OTU)[2]/3), classification = classification, num.trees = 300, importance= "impurity", write.forest = T)
        #assign(mod, ranger(refer ~ ., data=get(nam), mtry=floor(dim(get(nam))[2]/3), classification = classification, num.trees = 300, importance= "impurity", write.forest = T),envir = .GlobalEnv)
        #print("just after ranger")
        
        ## prediction for new data
        predict_tr_rf <- predict(get(mod), OTUte)
        combined1_rf <- c(combined1_rf,predict_tr_rf$predictions)
        combined2_rf <- c(combined2_rf,get(nam_comp_t)[,index])
      }
    }
    
    ## if random forest with bias correction (experimental..)
    if (algo == "RF_bc")
    {
      ### if optim and check for overfit 
      if (optim_overfit == T) 
      {
        ## gather RMSE for optim and check overfit
        over <- array(NA, c(2,10))
        rownames(over) <- c("training RMSE", "testing RMSE")
        colnames(over) <- c(10,50,100,150,200,250,300,350,400,500)
        for (opt in colnames(over))
        {
          mod <- paste("RF_farm_m", i,sep="_")
          set.seed(1)
          assign(mod,ranger(get(nam_comp)[,index] ~ ., data=get(nam), mtry=floor(dim(get(nam))[2]/3), classification = classification, num.trees = as.numeric(opt), importance= "impurity", write.forest = T))
          ## prediction for new data
          predict_tr_rf <- predict(get(mod), get(nam_t))
          ## paste in over the RMSE values
          over["training RMSE",opt] <- get(mod)$prediction.error
          over["testing RMSE",opt] <- sqrt(mean((get(nam_comp_t)[,index] - predict_tr_rf$predictions)^2))
        }
        ### then get the minimum testing RMSE to get the optim num.tree param
        best_over <- as.numeric(attributes(which.min(over["training RMSE",]))$names)
        print("num.tree effect on testing RMSE")
        print(over)
        print(paste("Using: ", best_over))
        # maybe export a plot and curves for overfitting check...
        mod <- paste("RF_farm_m", i,sep="_")
        set.seed(1)
        assign(mod,ranger(get(nam_comp)[,index] ~ ., data=get(nam), mtry=floor(dim(get(nam))[2]/3), classification = classification, num.trees = best_over, importance= "impurity", write.forest = T))
        ## prediction for new data
        predict_tr_rf <- predict(get(mod), get(nam_t))
        combined1_rf <- c(combined1_rf,predict_tr_rf$predictions)
        combined2_rf <- c(combined2_rf,get(nam_comp_t)[,index])
      } else 
      {
        # ## now the fitting, random forest with Ranger package with default mtry (1/3) for regression according to Breiman
        # similar here the mod objet has to be in the global env. to be fetchable afterwards
        assign("mod", paste("RF_bc_farm_m", i,sep="_"), envir = .GlobalEnv)
        set.seed(1)
        assign(mod, randomForest(get(nam_comp)[,index] ~ ., data=get(nam), mtry=floor(dim(get(nam))[2]/3), ntree = 300, importance= T, corr.bias = T),envir = .GlobalEnv)
        
        #refer <- get(nam_comp)[,index]
        #mod <- ranger(refer ~ ., data=get(nam), mtry=floor(dim(get(nam))[2]/3), classification = classification, num.trees = 300, importance= "impurity", write.forest = T)
        #mod <- ranger(COMP[,"AMBI"] ~ ., data=OTU, mtry=floor(dim(OTU)[2]/3), classification = classification, num.trees = 300, importance= "impurity", write.forest = T)
        #assign(mod, ranger(refer ~ ., data=get(nam), mtry=floor(dim(get(nam))[2]/3), classification = classification, num.trees = 300, importance= "impurity", write.forest = T),envir = .GlobalEnv)
        #print("just after ranger")
        
        ## prediction for new data
        predict_tr_rf <- predict(get(mod), get(nam_t))
        combined1_rf <- c(combined1_rf, as.vector(predict_tr_rf))
        combined2_rf <- c(combined2_rf,get(nam_comp_t)[,index])
      }
    }
    
    ### lasso with lars
    if (algo == "LASSO")
    {
      ### if optim and check for overfit 
      if (optim_overfit == T) 
      {
        ## gather RMSE for optim and check overfit
        over <- array(NA, c(2,10))
        rownames(over) <- c("training RMSE", "testing RMSE")
        colnames(over) <- c(10,50,100,150,200,250,300,350,400,500)
        for (opt in colnames(over))
        {
          mod <- paste("RF_farm_m", i,sep="_")
          set.seed(1)
          assign(mod,ranger(get(nam_comp)[,index] ~ ., data=get(nam), mtry=floor(dim(get(nam))[2]/3), classification = classification, num.trees = as.numeric(opt), importance= "impurity", write.forest = T))
          ## prediction for new data
          predict_tr_rf <- predict(get(mod), get(nam_t))
          ## paste in over the RMSE values
          over["training RMSE",opt] <- get(mod)$prediction.error
          over["testing RMSE",opt] <- sqrt(mean((get(nam_comp_t)[,index] - predict_tr_rf$predictions)^2))
        }
        ### then get the minimum testing RMSE to get the optim num.tree param
        best_over <- as.numeric(attributes(which.min(over["training RMSE",]))$names)
        print("num.tree effect on testing RMSE")
        print(over)
        print(paste("Using: ", best_over))
        # maybe export a plot and curves for overfitting check...
        mod <- paste("RF_farm_m", i,sep="_")
        set.seed(1)
        assign(mod,ranger(get(nam_comp)[,index] ~ ., data=get(nam), mtry=floor(dim(get(nam))[2]/3), classification = classification, num.trees = best_over, importance= "impurity", write.forest = T))
        ## prediction for new data
        predict_tr_rf <- predict(get(mod), get(nam_t))
        combined1_rf <- c(combined1_rf,predict_tr_rf$predictions)
        combined2_rf <- c(combined2_rf,get(nam_comp_t)[,index])
      } else 
      {
        # ## now the fitting, random forest with Ranger package with default mtry (1/3) for regression according to Breiman
        # similar here the mod objet has to be in the global env. to be fetchable afterwards
        assign("mod", paste("LAS_farm_m", i,sep="_"), envir = .GlobalEnv)
        
        ## removing empty OTUs
        OTUtr <- get(nam)[,colSums(get(nam)) != 0]
        # remove those ones for the testing dataset also
        OTUte <- get(nam_t)[,colSums(get(nam)) != 0]
        
        set.seed(1)
        assign(mod, lars(as.matrix(OTUtr), get(nam_comp)[,index], type="lar", use.Gram=FALSE))
        
        #refer <- get(nam_comp)[,index]
        #mod <- ranger(refer ~ ., data=get(nam), mtry=floor(dim(get(nam))[2]/3), classification = classification, num.trees = 300, importance= "impurity", write.forest = T)
        #mod <- ranger(COMP[,"AMBI"] ~ ., data=OTU, mtry=floor(dim(OTU)[2]/3), classification = classification, num.trees = 300, importance= "impurity", write.forest = T)
        #assign(mod, ranger(refer ~ ., data=get(nam), mtry=floor(dim(get(nam))[2]/3), classification = classification, num.trees = 300, importance= "impurity", write.forest = T),envir = .GlobalEnv)
        #print("just after ranger")
        
        ## prediction for new data
        predict_tr_la <- predict(get(mod), as.matrix(OTUte))
        combined1_la <- c(combined1_la, as.vector(predict_tr_la$fit[,max(predict_tr_la$s)]))
        combined2_la <- c(combined2_la,get(nam_comp_t)[,index])
      }
    }
    
    
    
    
    ### xgboost
    if (algo == "XGBOOST")
    {
      ## preparing the data
      
      dtrain <- xgb.DMatrix(get(nam), label = get(nam_comp)[,index])
      dtest <- xgb.DMatrix(get(nam_t), label = get(nam_comp_t)[,index])
      watchlist <- list(eval = dtest, train = dtrain)

      param <- list(max_depth = 6, eta = 0.3, silent = 1, nthread = 4,
                    objective = "reg:linear", eval_metric = "rmse")
      
      
      ### if optim and check for overfit 
      if (optim_overfit == T) 
      {
        ## gather RMSE for optim and check overfit
        over <- array(NA, c(2,10))
        rownames(over) <- c("training RMSE", "testing RMSE")
        colnames(over) <- c(10,50,100,150,200,250,300,350,400,500)
        for (opt in colnames(over))
        {
          mod <- paste("RF_farm_m", i,sep="_")
          set.seed(1)
          assign(mod,ranger(get(nam_comp)[,index] ~ ., data=get(nam), mtry=floor(dim(get(nam))[2]/3), classification = classification, num.trees = as.numeric(opt), importance= "impurity", write.forest = T))
          ## prediction for new data
          predict_tr_rf <- predict(get(mod), get(nam_t))
          ## paste in over the RMSE values
          over["training RMSE",opt] <- get(mod)$prediction.error
          over["testing RMSE",opt] <- sqrt(mean((get(nam_comp_t)[,index] - predict_tr_rf$predictions)^2))
        }
        ### then get the minimum testing RMSE to get the optim num.tree param
        best_over <- as.numeric(attributes(which.min(over["training RMSE",]))$names)
        print("num.tree effect on testing RMSE")
        print(over)
        print(paste("Using: ", best_over))
        # maybe export a plot and curves for overfitting check...
        mod <- paste("RF_farm_m", i,sep="_")
        set.seed(1)
        assign(mod,ranger(get(nam_comp)[,index] ~ ., data=get(nam), mtry=floor(dim(get(nam))[2]/3), classification = classification, num.trees = best_over, importance= "impurity", write.forest = T))
        ## prediction for new data
        predict_tr_rf <- predict(get(mod), get(nam_t))
        combined1_rf <- c(combined1_rf,predict_tr_rf$predictions)
        combined2_rf <- c(combined2_rf,get(nam_comp_t)[,index])
      } else 
      {
        # ## now the fitting, random forest with Ranger package with default mtry (1/3) for regression according to Breiman
        mod <- paste("XGBoost_farm_m", i,sep="_")
        set.seed(1)
        assign(mod,xgb.train(param, dtrain, nrounds = 100, watchlist))
        ## prediction for new data
        predict_tr_xg <- predict(get(mod), get(nam_t))
        combined1_xg <- c(combined1_rf,predict_tr_rf$predictions)
        combined2_xg <- c(combined2_rf,get(nam_comp_t)[,index])
      }
    }
    
    ### if LM with lm function
    if (algo =="LM")
    {
      if (optim_overfit == T) 
      {
        ## gather RMSE for optim and check overfit
        over <- array(NA, c(2,10))
        rownames(over) <- c("training RMSE", "testing RMSE")
        colnames(over) <- c(10,50,100,150,200,250,300,350,400,500)
        for (opt in colnames(over))
        {
          mod <- paste("DeepNet_farm_m", i,sep="_")
          mx.set.seed(1)
          assign(mod,mx.mlp(as.matrix(get(nam)), get(nam_comp)[,index], momentum=0.1, array.layout="rowmajor", learning.rate=0.01,hidden_node=100, out_node=1, dropout=NULL, num.round=opt, activation="tanh", out_activation='rmse', eval.metric=mx.metric.rmse))
          ## prediction on same data for RMSE (training error)
          preds_tr <- predict(get(mod), as.matrix(get(nam)), array.layout="rowmajor") 
          ## prediction for new data
          preds <- predict(get(mod), as.matrix(get(nam_t)), array.layout="rowmajor")
          ## paste in over the RMSE values
          over["training RMSE",opt] <- sqrt(mean((get(nam_comp)[,index] - preds_tr)^2)) 
          over["testing RMSE",opt] <- sqrt(mean((get(nam_comp_t)[,index] - preds)^2))
        }
        ### then get the minimum testing RMSE to get the optim num.tree param
        best_over <- as.numeric(attributes(which.min(over["training RMSE",]))$names)
        print("num.round effect on testing RMSE")
        print(over)
        print(paste("Using: ", best_over))
        # maybe export a plot and curves for overfitting check...
        mod <- paste("DeepNet_farm_m", i,sep="_")
        mx.set.seed(1)
        assign(mod,mx.mlp(as.matrix(get(nam)), get(nam_comp)[,index], momentum=0.1, array.layout="rowmajor", learning.rate=0.01,hidden_node=100, out_node=1, dropout=NULL, num.round=best_over, activation="tanh", out_activation='rmse', eval.metric=mx.metric.rmse))
        ## prediction for new data
        preds <- predict(get(mod), as.matrix(get(nam_t)), array.layout="rowmajor")
        combined1_dn <- c(combined1_dn,as.vector(preds))
        combined2_dn <- c(combined2_dn,get(nam_comp_t)[,index])
      } else 
      {
        mod <- paste("LM_farm_m", i,sep="_")
        set.seed(1)
        
        ## removing empty OTUs
        OTUtr <- get(nam)[,colSums(get(nam)) != 0]
        # remove those ones for the testing dataset also
        OTUte <- get(nam_t)[,colSums(get(nam)) != 0]
        
        
        assign(mod,lm(get(nam_comp)[,index] ~ ., data=OTUtr))
        
        preds <- predict(get(mod), OTUte)
        combined1_lm <- c(combined1_lm,as.vector(preds))
        combined2_lm <- c(combined2_lm,get(nam_comp_t)[,index])
      }
    }
    
    
    ### linear regression 
    # mod <- paste("LM_farm_m", i,sep="_")
    # set.seed(1)
    # assign(mod, lm(get(nam_comp)[,index] ~ ., data=get(nam)))
    # ## predictions for new data
    # predict_tr_lm <- predict(get(mod), get(nam_t))
    # combined1_lm <- c(combined1_lm,predict_tr_lm)
    # combined2_lm <- c(combined2_lm,get(nam_comp_t)[,index])
    
    
    # #### ranger with distance added
    # train <- cbind(get(nam), get(nam_comp)[,"Campaign"])
    # dimnames(train)[[2]][dim(train)[2]] <- "Campaign"
    # mod <- ranger(get(nam_comp)[,index] ~ Campaign/., data=train, mtry=floor(dim(train)[2]/3), importance= "impurity", num.trees = 300, write.forest = T)
    # test <- cbind(get(nam_t),get(nam_comp_t)[,"Campaign"])
    # dimnames(test)[[2]][dim(test)[2]] <- "Campaign"
    # predict_tr_rf <- predict(mod, test)
    # combined1_rf <- c(combined1_rf,predict_tr_rf$predictions)
    # combined2_rf <- c(combined2_rf,get(nam_comp_t)[,index])
    
    ### if SOM 
    if (algo == "SOM")
    {
      ### if optim and check for overfit 
      if (optim_overfit == T) 
      {
        ## gather RMSE for optim and check overfit BUT weird because random search of parameters, so hard to check...
        over <- array(NA, c(2,10))
        rownames(over) <- c("training RMSE", "testing RMSE")
        colnames(over) <- c(10,50,100,150,200,250,300,350,400,500)
        for (opt in colnames(over))
        {
          mod <- paste("SOM_farm_m", i,sep="_")
          set.seed(1)
          assign(mod,train(as.matrix(get(nam)), get(nam_comp)[,index], method = "xyf", trControl = fitControl, tuneLength = 100, rlen = as.numeric(opt)))
          ## prediction for new data
          predict_tr_sm <- predict(get(mod), as.matrix(get(nam_t)))
          ## paste in over the RMSE values
          over["training RMSE",opt] <- min(get(mod)$results[,"RMSE"])
          over["testing RMSE",opt] <- sqrt(mean((get(nam_comp_t)[,index] - predict_tr_sm)^2))
        }
        ### then get the minimum testing RMSE to get the optim num.tree param
        best_over <- as.numeric(attributes(which.min(over["training RMSE",]))$names)
        print("rlen effect on testing RMSE")
        print(over)
        print(paste("Using: ", best_over))
        # maybe export a plot and curves for overfitting check...
        mod <- paste("RF_farm_m", i,sep="_")
        set.seed(1)
        assign(mod,train(as.matrix(get(nam)), get(nam_comp)[,index], method = "xyf", trControl = fitControl, tuneLength = 100, rlen = best_over))
        ## prediction for new data
        predict_tr_sm <- predict(get(mod), get(nam_t))
        combined1_sm <- c(combined1_sm,predict_tr_sm)
        combined2_sm <- c(combined2_sm,get(nam_comp_t)[,index])
      } else 
      {
        # ## now the fitting, SOM
        mod <- paste("SOM_farm_m", i,sep="_")
        
        ## removing empty OTUs
        OTUtr <- get(nam)[,colSums(get(nam)) != 0]
        # remove those ones for the testing dataset also
        OTUte <- get(nam_t)[,colSums(get(nam)) != 0]
        
        
        set.seed(1)
        #assign(mod,train(as.matrix(get(nam)), get(nam_comp)[,index], method = "bdk", trControl = fitControl, tuneLength = 100, rlen = 100))
        #assign(mod,train(as.matrix(get(nam)), get(nam_comp)[,index], method = "bdk", trControl = fitControl, tuneLength = 100, rlen = best_over))
        assign(mod,xyf(as.matrix(OTUtr), get(nam_comp)[,index], grid=somgrid(3, 3, "rectangular"), rlen = 100, alpha = c(0.05, 0.01),
                       radius = c(1, -1), cores = -1, keep.data = TRUE))
                       #xweight = 0.5, contin = T, toroidal = F, n.hood = "circular"))
        
        #assign(mod,train(as.matrix(get(nam)), get(nam_comp)[,index], method = "bdk", trControl = fitControl, tuneLength = 100, rlen = best_over))
        
        ## prediction for new data
        predict_tr_sm <- predict(get(mod), newdata = as.matrix(OTUte), whatmap = 1)
        combined1_sm <- c(combined1_sm,as.vector(predict_tr_sm$predictions[[2]]))
        combined2_sm <- c(combined2_sm,get(nam_comp_t)[,index])
      }
      
    }
    
    
    ## now the fitting, self organizing map
    # mod <- paste("SelfMap_farm_m", i,sep="_")
    # assign(mod,train(as.matrix(get(nam)), get(nam_comp)[,index], method = "bdk", trControl = fitControl, tuneLength = 100, rlen = 100))
    # 
    # ## prediction for new data
    # predict_tr_sm <- predict (get(mod), as.matrix(get(nam_t)))
    # combined1_sm <- c(combined1_sm,predict_tr_sm$prediction)
    # combined2_sm <- c(combined2_sm,get(nam_comp_t)[,index])
    
    ### if DL with mxnet
    if (algo =="DL")
    {
      if (optim_overfit == T) 
      {
        ## gather RMSE for optim and check overfit
        over <- array(NA, c(2,10))
        rownames(over) <- c("training RMSE", "testing RMSE")
        colnames(over) <- c(10,50,100,150,200,250,300,350,400,500)
        for (opt in colnames(over))
        {
          mod <- paste("DeepNet_farm_m", i,sep="_")
          mx.set.seed(1)
          assign(mod,mx.mlp(as.matrix(get(nam)), get(nam_comp)[,index], momentum=0.1, array.layout="rowmajor", learning.rate=0.01,hidden_node=100, out_node=1, dropout=NULL, num.round=opt, activation="tanh", out_activation='rmse', eval.metric=mx.metric.rmse))
          ## prediction on same data for RMSE (training error)
          preds_tr <- predict(get(mod), as.matrix(get(nam)), array.layout="rowmajor") 
          ## prediction for new data
          preds <- predict(get(mod), as.matrix(get(nam_t)), array.layout="rowmajor")
          ## paste in over the RMSE values
          over["training RMSE",opt] <- sqrt(mean((get(nam_comp)[,index] - preds_tr)^2)) 
          over["testing RMSE",opt] <- sqrt(mean((get(nam_comp_t)[,index] - preds)^2))
        }
        ### then get the minimum testing RMSE to get the optim num.tree param
        best_over <- as.numeric(attributes(which.min(over["training RMSE",]))$names)
        print("num.round effect on testing RMSE")
        print(over)
        print(paste("Using: ", best_over))
        # maybe export a plot and curves for overfitting check...
        mod <- paste("DeepNet_farm_m", i,sep="_")
        mx.set.seed(1)
        assign(mod,mx.mlp(as.matrix(get(nam)), get(nam_comp)[,index], momentum=0.1, array.layout="rowmajor", learning.rate=0.01,hidden_node=100, out_node=1, dropout=NULL, num.round=best_over, activation="tanh", out_activation='rmse', eval.metric=mx.metric.rmse))
        ## prediction for new data
        preds <- predict(get(mod), as.matrix(get(nam_t)), array.layout="rowmajor")
        combined1_dn <- c(combined1_dn,as.vector(preds))
        combined2_dn <- c(combined2_dn,get(nam_comp_t)[,index])
      } else 
      {
        mod <- paste("DeepNet_farm_m", i,sep="_")
        mx.set.seed(1)
        assign(mod,mx.mlp(as.matrix(get(nam)), get(nam_comp)[,index], momentum=0.1, array.layout="rowmajor", learning.rate=0.01,hidden_node=100, out_node=1, dropout=NULL, num.round=100, activation="tanh", out_activation='rmse', eval.metric=mx.metric.rmse))
        
        preds <- predict(get(mod), as.matrix(get(nam_t)))
        combined1_dn <- c(combined1_dn,as.vector(preds))
        combined2_dn <- c(combined2_dn,get(nam_comp_t)[,index])
      }
    }
    
    
    ### if SVM with e1071
    if (algo =="SVM")
    {
      if (optim_overfit == T) 
      {
        ## gather RMSE for optim and check overfit
        over <- array(NA, c(2,10))
        rownames(over) <- c("training RMSE", "testing RMSE")
        colnames(over) <- c(10,50,100,150,200,250,300,350,400,500)
        for (opt in colnames(over))
        {
          mod <- paste("DeepNet_farm_m", i,sep="_")
          mx.set.seed(1)
          assign(mod,mx.mlp(as.matrix(get(nam)), get(nam_comp)[,index], momentum=0.1, array.layout="rowmajor", learning.rate=0.01,hidden_node=100, out_node=1, dropout=NULL, num.round=opt, activation="tanh", out_activation='rmse', eval.metric=mx.metric.rmse))
          ## prediction on same data for RMSE (training error)
          preds_tr <- predict(get(mod), as.matrix(get(nam)), array.layout="rowmajor") 
          ## prediction for new data
          preds <- predict(get(mod), as.matrix(get(nam_t)), array.layout="rowmajor")
          ## paste in over the RMSE values
          over["training RMSE",opt] <- sqrt(mean((get(nam_comp)[,index] - preds_tr)^2)) 
          over["testing RMSE",opt] <- sqrt(mean((get(nam_comp_t)[,index] - preds)^2))
        }
        ### then get the minimum testing RMSE to get the optim num.tree param
        best_over <- as.numeric(attributes(which.min(over["training RMSE",]))$names)
        print("num.round effect on testing RMSE")
        print(over)
        print(paste("Using: ", best_over))
        # maybe export a plot and curves for overfitting check...
        mod <- paste("DeepNet_farm_m", i,sep="_")
        mx.set.seed(1)
        assign(mod,mx.mlp(as.matrix(get(nam)), get(nam_comp)[,index], momentum=0.1, array.layout="rowmajor", learning.rate=0.01,hidden_node=100, out_node=1, dropout=NULL, num.round=best_over, activation="tanh", out_activation='rmse', eval.metric=mx.metric.rmse))
        ## prediction for new data
        preds <- predict(get(mod), as.matrix(get(nam_t)), array.layout="rowmajor")
        combined1_dn <- c(combined1_dn,as.vector(preds))
        combined2_dn <- c(combined2_dn,get(nam_comp_t)[,index])
      } else 
      {
        mod <- paste("SVM_farm_m", i,sep="_")
        
        ## removing empty OTUs
        OTUtr <- get(nam)[,colSums(get(nam)) != 0]
        # remove those ones for the testing dataset also
        OTUte <- get(nam_t)[,colSums(get(nam)) != 0]
        
        set.seed(1)
        assign(mod,svm(as.matrix(OTUtr), get(nam_comp)[,index], scale = F, type = "eps-regression", kernel = "sigmoid"))
        
        preds <- predict(get(mod), OTUte)
        combined1_sv <- c(combined1_sv,as.vector(preds))
        combined2_sv <- c(combined2_sv,get(nam_comp_t)[,index])
      }
    }
    
    # gbm 
    #gbm_mod <- paste("GBM_farm_m", i,sep="_")
    #assign(gbm_mod,gbm(get(nam_comp)[,index] ~ ., data=get(nam), n.trees = 10, shrinkage = 0.1, interaction.depth = 1))
    #predict_tr_gbm <- predict(get(gbm_mod), get(nam_t), 10)
    
    ### if SVM with e1071
    if (algo =="NN")
    {
      if (optim_overfit == T) 
      {
        ## gather RMSE for optim and check overfit
        over <- array(NA, c(2,10))
        rownames(over) <- c("training RMSE", "testing RMSE")
        colnames(over) <- c(10,50,100,150,200,250,300,350,400,500)
        for (opt in colnames(over))
        {
          mod <- paste("DeepNet_farm_m", i,sep="_")
          mx.set.seed(1)
          assign(mod,mx.mlp(as.matrix(get(nam)), get(nam_comp)[,index], momentum=0.1, array.layout="rowmajor", learning.rate=0.01,hidden_node=100, out_node=1, dropout=NULL, num.round=opt, activation="tanh", out_activation='rmse', eval.metric=mx.metric.rmse))
          ## prediction on same data for RMSE (training error)
          preds_tr <- predict(get(mod), as.matrix(get(nam)), array.layout="rowmajor") 
          ## prediction for new data
          preds <- predict(get(mod), as.matrix(get(nam_t)), array.layout="rowmajor")
          ## paste in over the RMSE values
          over["training RMSE",opt] <- sqrt(mean((get(nam_comp)[,index] - preds_tr)^2)) 
          over["testing RMSE",opt] <- sqrt(mean((get(nam_comp_t)[,index] - preds)^2))
        }
        ### then get the minimum testing RMSE to get the optim num.tree param
        best_over <- as.numeric(attributes(which.min(over["training RMSE",]))$names)
        print("num.round effect on testing RMSE")
        print(over)
        print(paste("Using: ", best_over))
        # maybe export a plot and curves for overfitting check...
        mod <- paste("DeepNet_farm_m", i,sep="_")
        mx.set.seed(1)
        assign(mod,mx.mlp(as.matrix(get(nam)), get(nam_comp)[,index], momentum=0.1, array.layout="rowmajor", learning.rate=0.01,hidden_node=100, out_node=1, dropout=NULL, num.round=best_over, activation="tanh", out_activation='rmse', eval.metric=mx.metric.rmse))
        ## prediction for new data
        preds <- predict(get(mod), as.matrix(get(nam_t)), array.layout="rowmajor")
        combined1_dn <- c(combined1_dn,as.vector(preds))
        combined2_dn <- c(combined2_dn,get(nam_comp_t)[,index])
      } else 
      {
        mod <- paste("NN_farm_m", i,sep="_")
        
        ## removing empty OTUs
        OTUtr <- get(nam)[,colSums(get(nam)) != 0]
        # remove those ones for the testing dataset also
        OTUte <- get(nam_t)[,colSums(get(nam)) != 0]
        
        formula <- as.formula(paste('get(nam_comp)[,index] ~ ' ,paste(dimnames(OTUtr)[[2]],collapse='+')))

        set.seed(1)
        assign(mod,neuralnet(formula, data = get(nam), hidden=10, err.fct="sse", algorithm = "rprop+"))
        
        predict_tr_nn <- compute(get(mod), OTUte)
        combined1_nn <- c(combined1_nn,predict_tr_nn$net.result[,1])
        combined2_nn <- c(combined2_nn,get(nam_comp_t)[,index])
      }
    }
    
    # # neural net
    # mod <- paste("NN_farm_m", i,sep="_")
    # formula <- as.formula(paste('get(nam_comp)[,index] ~ ' ,paste(dimnames(get(nam))[[2]],collapse='+')))
    # assign(mod,neuralnet(formula, data = get(nam), hidden=10, err.fct="sse"))
    # 
    # predict_tr_nn <- compute(get(mod), get(nam_t))
    # combined1_nn <- c(combined1_nn,predict_tr_nn$net.result[,1])
    # combined2_nn <- c(combined2_nn,get(nam_comp_t)[,index])
    
    # and keep track af farms to samples correspondance
    farm_nam <- c(farm_nam, rep(x = i, dim(get(nam_t))[1]))
    cat(cpt, '/', length(cross), " - ", i, ": done\n")
    cpt <- cpt+1
  }
  
  if (algo =="RF" | algo == "RF_bc" | algo == "RFProb") return(combined1_rf)
  if (algo =="SOM") return(combined1_sm)
  if (algo =="DL") return(combined1_dn)
  if (algo =="SVM") return(combined1_sv)
  if (algo =="LM") return(combined1_lm)
  if (algo =="NN") return(combined1_nn)
  if (algo =="LASSO") return(combined1_la)
}


