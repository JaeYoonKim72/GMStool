  DNN_func = function (train_pheno, train_geno, test_geno) {
  set.seed(10)
  tensorflow::tf$random$set_seed(10)
  tf$compat$v1$set_random_seed(10)
  train_geno_mnum = train_geno + 1
  train_geno_mnum = tf$cast(train_geno_mnum, tf$float32) / 3
  test_geno_mnum = test_geno + 1
  test_geno_mnum <- tf$cast(test_geno_mnum, tf$float32) /3
  batchs = round(dim(train_geno_mnum)[1]/20)
  if ((batchs %% 2) == 1){batchs = batchs + 1}
  model <- keras_model_sequential()
  model %>%
    layer_dense(units = 256, activation = "relu", batch_size = batchs, input_shape = dim(train_geno_mnum)[2], kernel_initializer = 'orthogonal', kernel_regularizer=regularizer_l1_l2(l1 = 0.3, l2 = 0.3)) %>%
    layer_dropout(rate = 0.03) %>%
    layer_dense(units = 128, activation = "relu", batch_size = batchs, kernel_initializer = 'orthogonal', kernel_regularizer=regularizer_l1_l2(l1 = 0.2, l2 = 0.2)) %>%
    layer_dropout(rate = 0.02) %>%
    layer_dense(units = 64, activation = "relu", batch_size = batchs, kernel_initializer = 'orthogonal', kernel_regularizer=regularizer_l1_l2(l1 = 0.0, l2 = 0)) %>%
    layer_dropout(rate = 0.01) %>%
    layer_dense(units = 32, activation = "relu", batch_size = batchs, kernel_initializer = 'orthogonal', kernel_regularizer=regularizer_l1_l2(l1 = 0.0, l2 = 0)) %>%
    layer_batch_normalization(batch_size = batchs) %>%
    layer_dense(units = 16, activation = "linear", batch_size = batchs, kernel_initializer = 'orthogonal', kernel_regularizer=regularizer_l1_l2(l1 = 0.0, l2 = 0)) %>%
    layer_dense(units = 1, activation = "linear", batch_size = batchs)
  model %>% compile(loss = "mse", optimizer = optimizer_adamax(lr=0.01, decay = 0.0003), metrics = c(metric_r2_score, metric_cor))
  history <- model %>%
    fit(train_geno_mnum, train_pheno, epochs = 1000, batch_size = batchs, validation_split = 0.2, verbose = 0,
        callbacks = list(callback_early_stopping(patience = 50), callback_reduce_lr_on_plateau(factor = 0.1)))
  #score = model %>% evaluate(test_geno_mnum, test_pheno, batch_size = batchs, verbose=0) %>% unlist
  pred_pheno <- model %>% predict(test_geno_mnum, batch_size = batchs)
  return_value = list("predicted"=pred_pheno, "model"=model, "batch"=batchs)
  return(return_value)
  }

