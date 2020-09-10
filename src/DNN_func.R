DNN_func = function (train_pheno, train_geno, test_geno) {
  train_geno_mnum = train_geno + 1
  train_geno_mnum = tf$cast(train_geno_mnum, tf$float32) / 3
  test_geno_mnum = test_geno + 1
  test_geno_mnum <- tf$cast(test_geno_mnum, tf$float32) /3
  batchs = round(dim(train_geno_mnum)[1]/20)
  if ((batchs %% 2) == 1){batchs = batchs + 1}
  model <- keras_model_sequential()
  model %>%
    layer_dense(units = 256, activation = "relu", input_shape = dim(train_geno_mnum)[2], kernel_initializer = 'orthogonal', kernel_regularizer=regularizer_l1_l2(l1 = 0.3, l2 = 0.3)) %>%
    layer_dropout(rate = 0.03) %>%
    layer_dense(units = 128, activation = "relu", kernel_initializer = 'orthogonal', kernel_regularizer=regularizer_l1_l2(l1 = 0.2, l2 = 0.2)) %>%
    layer_dropout(rate = 0.02) %>%
    layer_dense(units = 64, activation = "relu", kernel_initializer = 'orthogonal', kernel_regularizer=regularizer_l1_l2(l1 = 0.0, l2 = 0)) %>%
    layer_dropout(rate = 0.01) %>%
    layer_dense(units = 32, activation = "relu", kernel_initializer = 'orthogonal', kernel_regularizer=regularizer_l1_l2(l1 = 0.0, l2 = 0)) %>%
    layer_batch_normalization(batch_size = batchs) %>%
    layer_dense(units = 16, activation = "linear", kernel_initializer = 'orthogonal', kernel_regularizer=regularizer_l1_l2(l1 = 0.0, l2 = 0)) %>%
    layer_dense(units = 1, activation = "linear")
  model %>% compile(loss = "mse", optimizer = optimizer_adamax(lr=0.001, decay = 0.0003), metrics = c(metric_r2_score, metric_cor))
  tensorflow::tf$random$set_seed(200)
  tf$compat$v1$set_random_seed(200)
  set.seed(200)
  tensorflow::tf$random$set_seed(200)
  history <- model %>%  
    fit(train_geno_mnum, train_pheno, epochs = 1000, batch_size = batchs, validation_split = 0.2, verbose = 0,
        callbacks = list(callback_early_stopping(patience = 30), callback_reduce_lr_on_plateau(factor = 0.1)))
  train_pred <- model %>% predict(train_geno_mnum, batch_size = batchs)
  val_pred <- model %>% predict(test_geno_mnum, batch_size = batchs)
  return_value = list("val_predicted"=val_pred, "train_predicted"=train_pred, "model"=model, "batch"=batchs)
  return(return_value)
}
