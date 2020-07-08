  CNN_func = function (train_pheno, train_geno, test_geno) {
  set.seed(10)
  tensorflow::tf$random$set_seed(10)
  tf$compat$v1$set_random_seed(10)
  mnum = dim(train_geno)[2]
  train_geno_mnum <- train_geno + 1
  dim(train_geno_mnum) <- c(dim(train_geno_mnum)[1], 1, dim(train_geno_mnum)[2])
  train_geno_mnum <- tf$cast(train_geno_mnum, tf$float32) /3
  train_geno_mnum <- tf$expand_dims(train_geno_mnum, axis=-1L)
  test_geno_mnum = test_geno + 1
  dim(test_geno_mnum) <- c(dim(test_geno_mnum)[1], 1, dim(test_geno_mnum)[2])
  test_geno_mnum <- tf$cast(test_geno_mnum, tf$float32) /3
  test_geno_mnum <- tf$expand_dims(test_geno_mnum, axis=-1L)
  batchs = round(dim(train_geno_mnum)[1]/20)
  if ((batchs %% 2) == 1){batchs = batchs + 1}
  model <- keras_model_sequential() %>%
    layer_conv_2d(filters = 32, kernel_size = c(1,14), strides = c(1,4), batch_size = batchs, padding="same", activation= "relu", kernel_initializer = 'orthogonal', kernel_regularizer=regularizer_l1_l2(l1 = 0.1, l2 = 0.1), input_shape = shape(1, mnum, 1)) %>%
    layer_dropout(rate = 0.2) %>%
    layer_conv_2d(filters = 16, kernel_size = c(1,10), strides = c(1,3), batch_size = batchs, padding="same", activation = "linear", kernel_initializer = 'orthogonal', kernel_regularizer=regularizer_l1_l2(l1 = 0.00, l2 = 0.00)) %>%
    layer_dropout(rate = 0.1) %>%
    layer_conv_2d(filters = 8, kernel_size = c(1,8), strides = c(1,2), batch_size = batchs, padding="same", activation = "linear", kernel_initializer = 'orthogonal', kernel_regularizer=regularizer_l1_l2(l1 = 0.00, l2 = 0.00)) %>%
    layer_max_pooling_2d(pool_size = c(1,2))  %>%
    layer_batch_normalization(batch_size = batchs)
  model %>%
    layer_flatten() %>%
    layer_dense(units = 64, activation = "linear", batch_size = batchs, kernel_initializer = 'orthogonal', kernel_regularizer=regularizer_l1_l2(l1 = 0.00, l2 = 0.00)) %>%
    layer_dense(units = 32, activation = "linear", batch_size = batchs, kernel_initializer = 'orthogonal', kernel_regularizer=regularizer_l1_l2(l1 = 0.00, l2 = 0.00)) %>%
    layer_dense(units = 16, activation = "linear", batch_size = batchs, kernel_initializer = 'orthogonal', kernel_regularizer=regularizer_l1_l2(l1 = 0.00, l2 = 0.00)) %>%
    layer_batch_normalization(batch_size = batchs)%>%
    layer_dense(units = 1, batch_size = batchs, activation = "linear")
   model %>% compile(loss = "mse", optimizer = optimizer_adamax(lr=0.01, decay = 0.0003), metrics = c(metric_r2_score, metric_cor))
  history <- model %>%
    fit(train_geno_mnum, train_pheno, epochs = 1000, batch_size = batchs, validation_split = 0.2, verbose = 0,
        callbacks = list(callback_early_stopping(patience = 50), callback_reduce_lr_on_plateau(factor = 0.1)))
  pred_pheno <- model %>% predict(test_geno_mnum, batch_size = batchs)
  return_value = list("predicted"=pred_pheno, "model"=model, "batch"=batchs)
  return(return_value)
  }

