context("loglik")
library(nloptr)
test_that("the analytic gradient is the same as numerical gradient", {
  set.seed(1234)
  H = 4
  n_B = 60
  K = 1
  M = 3
  coefs = c(rep(0, K), rep(1, H*M), overdispersion=10)  # initial value
  cell_type_specific_variables = matrix(0, nrow=n_B, ncol=M)
  cell_type_specific_variables[1:(n_B/3), 1] = 1
  cell_type_specific_variables[(n_B/3+1):(n_B*2/3), 2] = 1
  cell_type_specific_variables[(n_B*2/3+1):n_B, 3] = 1
  other_variables = matrix(runif(n_B * K, min = -1, max = 1), nrow=n_B, ncol=K)
  read_depth = round(runif(n_B, min = 50, max = 100))
  cellular_proportions = matrix(runif(n_B * H), nrow = n_B, ncol = H)
  cellular_proportions = cellular_proportions / rowSums(cellular_proportions)
  counts = round(runif(n_B, min = 100, max = 200))
  grad_loglik(coefs, cell_type_specific_variables, other_variables, read_depth, cellular_proportions, counts)
  expect_equal(grad_loglik(coefs, cell_type_specific_variables, other_variables, read_depth, cellular_proportions, counts),
               nloptr::nl.grad(coefs, loglik, cell_type_specific_variables=cell_type_specific_variables, other_variables=other_variables, read_depth=read_depth, cellular_proportions=cellular_proportions, counts=counts))
})

test_that("the analytic gradient is the same as numerical gradient when there is no cell type-independent variables", {
  set.seed(1234)
  H = 4
  n_B = 60
  K = 0
  M = 3
  coefs = c(rep(0, K), rep(1, H*M), overdispersion=10)  # initial value
  cell_type_specific_variables = matrix(0, nrow=n_B, ncol=M)
  cell_type_specific_variables[1:(n_B/3), 1] = 1
  cell_type_specific_variables[(n_B/3+1):(n_B*2/3), 2] = 1
  cell_type_specific_variables[(n_B*2/3+1):n_B, 3] = 1
  # other_variables = matrix(runif(n_B * K, min = -1, max = 1), nrow=n_B, ncol=K)
  other_variables = NULL
  read_depth = round(runif(n_B, min = 50, max = 100))
  cellular_proportions = matrix(runif(n_B * H), nrow = n_B, ncol = H)
  cellular_proportions = cellular_proportions / rowSums(cellular_proportions)
  counts = round(runif(n_B, min = 100, max = 200))
  grad_loglik(coefs, cell_type_specific_variables, other_variables, read_depth, cellular_proportions, counts)
  expect_equal(grad_loglik(coefs, cell_type_specific_variables, other_variables, read_depth, cellular_proportions, counts),
               nloptr::nl.grad(coefs, loglik, cell_type_specific_variables=cell_type_specific_variables, other_variables=other_variables, read_depth=read_depth, cellular_proportions=cellular_proportions, counts=counts))
})