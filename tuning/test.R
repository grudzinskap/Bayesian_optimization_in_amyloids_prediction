library(mlrMBO)

opt <- function(seed, lerner_len, measure, n_iter=100, sigma=1, crit='ei', lambda=NULL){
  obj.fun <- makeSingleObjectiveFunction(
    name = "rf",
    fn = function(x) {
      MBOtarget(x[1], x[2], seed, lerner_len, measure)
    },
    par.set = makeParamSet(
      makeIntegerParam("n_", lower = 100, upper = 1000),
      makeIntegerParam("mtry_", lower = 3, upper = 30)
    ),
    minimize = FALSE
  )
  
  des = generateDesign(n = 5, par.set = getParamSet(obj.fun), fun = lhs::randomLHS)
  des$y = apply(des, 1, obj.fun)
  
  surr.km = makeLearner("regr.gausspr", predict.type = "se", kernel = "rbfdot", sigma = sigma)
  #surr.km = makeLearner("regr.km", predict.type = "se", covtype = "matern3_2", 
  #                      control = list(trace = FALSE), nugget.stability = 10e-8)
  control = makeMBOControl()
  control = setMBOControlTermination(control, iters = n_iter)
  if(crit=='ei'){
    control = setMBOControlInfill(control, crit = makeMBOInfillCritEI())
  }
  else{
    control = setMBOControlInfill(control, crit = makeMBOInfillCritCB(cb.lambda = lambda))
    #control = setMBOControlInfill(control, crit = makeMBOInfillCritAdaCB(cb.lambda.start = 10, cb.lambda.end = 1))
    
  }
  
  run <- mbo(obj.fun, design = des, learner = surr.km, control = control, show.info = TRUE)
  #lambda <- 0
return(c(seed, lerner_len, measure, n_iter, sigma, crit, lambda, run$y, run$x$n_, run$x$mtry_))
}

res_ei <- c()
for (seed in c(1)) {
  for (lerner_len in c(10)) {
    for (sigma in c(2)) {
      print(paste('seed', seed))
      print(paste('ll', lerner_len))
      print(paste('sigma', sigma))
      set.seed(seed)
      res <- opt(seed=seed, lerner_len=lerner_len, measure='AUC', 
             n_iter=50, sigma=sigma, crit='ei', lambda=NULL)
      res_ei <- rbind(res_ei, res)
      write.csv(data.frame(res_ei), 'ei_matern.csv')
    }
  }
}


#0.1, 0.5,  1, 1.5, 2, 2.5
res_cb <- c()
for (seed in c(1)) {
  for (lerner_len in c(10)) {
    for (sigma in c(2,5)) {
     for (lambda in c(4)) {
        set.seed(seed)
        print(paste('seed', seed))
        print(paste('ll', lerner_len))
        print(paste('sigma', sigma))
        #print(paste('lambda', lambda))
        res <- opt(seed=seed, lerner_len=lerner_len, measure='AUC', 
                   n_iter=50, sigma=sigma, crit='cb', lambda=lambda)
        res_cb <- rbind(res_cb, res)
        write.csv(data.frame(res_cb), 'cb_missing.csv')
      }
    }
  }
}

res_cb <- c()
for (seed in c(1)) {
  for (lerner_len in c(10)) {
    for (sigma in c(10)) {
      # for (lambda in c(1, 2, 3, 4, 5)) {
      set.seed(seed)
      print(paste('seed', seed))
      print(paste('ll', lerner_len))
      print(paste('sigma', sigma))
      #print(paste('lambda', lambda))
      res <- opt(seed=seed, lerner_len=lerner_len, measure='AUC', 
                 n_iter=50, sigma=sigma, crit='cb', lambda=lambda)
      res_cb <- rbind(res_cb, res)
      write.csv(data.frame(res_cb), 'cb_progressive_m_3_2.csv')
      #}
    }
  }
}
