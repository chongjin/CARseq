context("loglik")
test_that("the analytic gradient is the same as numerical gradient", {
  set.seed(1234)
  H = 4
  n_B = 60
  K = 1
  M = 3
  coefs = c(rep(0, K), rep(1, H*M), overdispersion=10)  # initial value
  cell_type_specific_variables = array(0, dim=c(n_B, H, M))
  cell_type_specific_variables[1:(n_B/3), , 1] = 1
  cell_type_specific_variables[(n_B/3+1):(n_B*2/3), , 2] = 1
  cell_type_specific_variables[(n_B*2/3+1):n_B, , 3] = 1
  other_variables = matrix(runif(n_B * K, min = -1, max = 1), nrow=n_B, ncol=K)
  read_depth = round(runif(n_B, min = 50, max = 100))
  cellular_proportions = matrix(runif(n_B * H), nrow = n_B, ncol = H)
  cellular_proportions = cellular_proportions / rowSums(cellular_proportions)
  counts = round(runif(n_B, min = 100, max = 200))
  grad_negloglik(coefs, cell_type_specific_variables, other_variables, read_depth, cellular_proportions, counts)
  expect_equal(grad_negloglik(coefs, cell_type_specific_variables, other_variables, read_depth, cellular_proportions, counts),
               nloptr::nl.grad(coefs, negloglik, cell_type_specific_variables=cell_type_specific_variables, other_variables=other_variables, read_depth=read_depth, cellular_proportions=cellular_proportions, counts=counts))
})

test_that("the analytic gradient is the same as numerical gradient when there is no cell type-independent variables", {
  set.seed(1234)
  H = 4
  n_B = 60
  K = 0
  M = 3
  coefs = c(rep(0, K), rep(1, H*M), overdispersion=10)  # initial value
  cell_type_specific_variables = array(0, dim=c(n_B, H, M))
  cell_type_specific_variables[1:(n_B/3), , 1] = 1
  cell_type_specific_variables[(n_B/3+1):(n_B*2/3), , 2] = 1
  cell_type_specific_variables[(n_B*2/3+1):n_B, , 3] = 1
  # other_variables = matrix(runif(n_B * K, min = -1, max = 1), nrow=n_B, ncol=K)
  other_variables = NULL
  read_depth = round(runif(n_B, min = 50, max = 100))
  cellular_proportions = matrix(runif(n_B * H), nrow = n_B, ncol = H)
  cellular_proportions = cellular_proportions / rowSums(cellular_proportions)
  counts = round(runif(n_B, min = 100, max = 200))
  grad_negloglik(coefs, cell_type_specific_variables, other_variables, read_depth, cellular_proportions, counts)
  expect_equal(grad_negloglik(coefs, cell_type_specific_variables, other_variables, read_depth, cellular_proportions, counts),
               nloptr::nl.grad(coefs, negloglik, cell_type_specific_variables=cell_type_specific_variables, other_variables=other_variables, read_depth=read_depth, cellular_proportions=cellular_proportions, counts=counts))
})

test_that("the analytic gradient is the same as numerical gradient in a reduced model", {
  set.seed(1234)
  H = 4
  n_B = 60
  K = 1
  M = 3
  coefs = c(rep(0, K), rep(1, H*M), overdispersion=10)  # initial value
  cell_type_specific_variables = array(0, dim=c(n_B, H, M))
  cell_type_specific_variables[1:(n_B/3), , 1] = 1
  cell_type_specific_variables[(n_B/3+1):(n_B*2/3), , 2] = 1
  cell_type_specific_variables[(n_B*2/3+1):n_B, , 3] = 1
  # The reduced model dictates that cell type 1 expresses the same in all 3 groups:
  cell_type_specific_variables[, 1, 1] = 1
  cell_type_specific_variables[, 1, 2:3] = 0
  other_variables = matrix(runif(n_B * K, min = -1, max = 1), nrow=n_B, ncol=K)
  read_depth = round(runif(n_B, min = 50, max = 100))
  cellular_proportions = matrix(runif(n_B * H), nrow = n_B, ncol = H)
  cellular_proportions = cellular_proportions / rowSums(cellular_proportions)
  counts = round(runif(n_B, min = 100, max = 200))
  grad_negloglik(coefs, cell_type_specific_variables, other_variables, read_depth, cellular_proportions, counts)
  expect_equal(grad_negloglik(coefs, cell_type_specific_variables, other_variables, read_depth, cellular_proportions, counts),
               nloptr::nl.grad(coefs, negloglik, cell_type_specific_variables=cell_type_specific_variables, other_variables=other_variables, read_depth=read_depth, cellular_proportions=cellular_proportions, counts=counts))
})

test_that("The results from CARseq and glm.nb+nnls matches when there is no cell type-independent variables", {
  set.seed(1234)
  H = 4
  n_B = 60
  K = 0
  M = 2
  coefs = c(rep(0, K), rep(1, H*M), overdispersion=10)  # initial value
  cell_type_specific_variables = array(0, dim=c(n_B, H, M))
  # case control label
  x = rep(0, n_B)
  x[1:(n_B/2)] = 1
  cell_type_specific_variables[x == 0, , 1] = 1
  cell_type_specific_variables[x == 1, , 2] = 1
  # other_variables = matrix(runif(n_B * K, min = -1, max = 1), nrow=n_B, ncol=K)
  other_variables = NULL
  read_depth = round(runif(n_B, min = 50, max = 100))
  cellular_proportions = matrix(runif(n_B * H), nrow = n_B, ncol = H)
  cellular_proportions = cellular_proportions / rowSums(cellular_proportions)
  counts = round(runif(n_B, min = 100, max = 200))
  
  # CARseq
  overdispersion_theta_init = 50
  theta_max=1e4
  theta_min=1e-3
  x0 = c(rep(1, H*M), overdispersion_theta_init)
  upper = c(rep(Inf, H*M), theta_max)
  lower = c(rep(1e-30, H*M), theta_min)
  
  res_nlminb_full = stats::nlminb(start = x0, 
                                  objective = CARseq::negloglik, gradient = CARseq::grad_negloglik, 
                                  # hessian = hessian_negloglik,
                                  control=list(eval.max=1000),
                                  cell_type_specific_variables = cell_type_specific_variables,
                                  other_variables = NULL,
                                  read_depth = read_depth,
                                  cellular_proportions = cellular_proportions,
                                  counts = counts)
  
  # glm.nb+nnls
  control = list(maxit = 200, trace = FALSE, epsilon = 1e-8)
  
  res_glmnb_full = MASS::glm.nb(counts ~ 0 + cellular_proportions:read_depth:(x==1),
                                method="glm.fit.cons.nonneg", link="identity", control = control)
  
  expect_equal(res_nlminb_full$par[1:(H*M)],
               as.numeric(res_glmnb_full$coefficients),
               tolerance = 1e-4)
})

test_that("The results from CARseq and glm.nb+nnls matches in a reduced model when there is no cell type-independent variables", {
  set.seed(1234)
  H = 4
  n_B = 60
  K = 0
  M = 2
  coefs = c(rep(0, K), rep(1, H*M), overdispersion=10)  # initial value
  cell_type_specific_variables_reduced = array(0, dim=c(n_B, H, M))
  # case control label
  x = rep(0, n_B)
  x[1:(n_B/2)] = 1
  cell_type_specific_variables_reduced[x == 0, , 1] = 1
  cell_type_specific_variables_reduced[x == 1, , 2] = 1
  # For the cell type specified by indices:
  indices = 1  # 1st cell type
  cell_type_specific_variables_reduced[, indices, 1] = 1
  cell_type_specific_variables_reduced[, indices, 2] = 0
  
  # We calculate which (cell type, cell type-specific variable) pairs are not included in the model
  # For future use.
  is_reduced = apply(cell_type_specific_variables_reduced, c(2,3),
                     function(x) identical(x, rep(0, dim(cell_type_specific_variables_reduced)[1])))
  # other_variables = matrix(runif(n_B * K, min = -1, max = 1), nrow=n_B, ncol=K)
  other_variables = NULL
  read_depth = round(runif(n_B, min = 50, max = 100))
  cellular_proportions = matrix(runif(n_B * H), nrow = n_B, ncol = H)
  cellular_proportions = cellular_proportions / rowSums(cellular_proportions)
  counts = round(runif(n_B, min = 100, max = 200))
  
  # CARseq
  overdispersion_theta_init = 50
  theta_max=1e4
  theta_min=1e-3
  x0 = c(rep(1, H*M), overdispersion_theta_init)
  upper = c(rep(Inf, H*M), theta_max)
  lower = c(rep(1e-30, H*M), theta_min)
  
  res_nlminb = stats::nlminb(start = x0, 
                             objective = CARseq::negloglik, gradient = CARseq::grad_negloglik, 
                             # hessian = hessian_negloglik,
                             control=list(eval.max=1000),
                             cell_type_specific_variables = cell_type_specific_variables_reduced,
                             other_variables = NULL,
                             read_depth = read_depth,
                             cellular_proportions = cellular_proportions,
                             counts = counts)
  # transpose and only take the coefficients in reduced model
  res_nlminb_par_reduced = res_nlminb$par
  res_nlminb_par_reduced[which(is_reduced)] = NA
  res_nlminb_par_reduced_mat = 
    matrix(res_nlminb_par_reduced[-length(res_nlminb_par_reduced)], ncol=2)
  res_nlminb_par_reduced = na.omit(as.numeric(t(res_nlminb_par_reduced_mat)))

  # glm.nb+nnls
  control = list(maxit = 200, trace = FALSE, epsilon = 1e-8)
  
  colnames(cellular_proportions) = c("V1", "V2", "V3", "V4")
  res_glmnb = MASS::glm.nb(counts ~ 0 +
                             cellular_proportions[,"V1"]:read_depth +
                             cellular_proportions[,"V2"]:read_depth:(x==1) +
                             cellular_proportions[,"V3"]:read_depth:(x==1) +
                             cellular_proportions[,"V4"]:read_depth:(x==1),
                           method="glm.fit.cons.nonneg", link="identity", control = control)
  
  expect_equal(as.numeric(res_nlminb_par_reduced),
               as.numeric(res_glmnb$coefficients),
               tolerance = 1e-4)
})
