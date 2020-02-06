
test_that("aof return data.frame",{
  mu1 <- 50
  mu2 <- 50
  rho1 = 0.5
  rho2 = 0.5
  n.obs <- 5
  sigma1 <- 0.1
  sigma2 <- sigma1
  SimTS <- function(n, mu, rho, sigma){
    X.standard <- arima.sim(n, model = list(ar = rho))
    X.standard/sd(X.standard)*sigma + mu
  }
  t.full <- 0:50
  t.break <- 25
  x.full <- c(SimTS(t.break, mu1, rho1, sigma1),
    SimTS(max(t.full)-t.break+1, mu2, rho2, sigma2))
  keep <- sort(sample(1:length(x.full), n.obs))
  TimeBudget <- data.frame(
    name = "A",
    Age = t.full[keep],
    x = x.full[keep])
  AOF <- aof(
    name = TimeBudget$name,
    Age = TimeBudget$Age,
    x = TimeBudget$x)
  expect_equal(
    object = class(AOF),
    expected = "data.frame"
  )
})

test_that("aof return data.frame with 1 obs",{
  mu1 <- 50
  mu2 <- 50
  rho1 = 0.5
  rho2 = 0.5
  n.obs <- 1
  sigma1 <- 0.1
  sigma2 <- sigma1
  SimTS <- function(n, mu, rho, sigma){
    X.standard <- arima.sim(n, model = list(ar = rho))
    X.standard/sd(X.standard)*sigma + mu
  }
  t.full <- 0:50
  t.break <- 25
  x.full <- c(SimTS(t.break, mu1, rho1, sigma1),
              SimTS(max(t.full)-t.break+1, mu2, rho2, sigma2))
  keep <- sort(sample(1:length(x.full), n.obs))
  TimeBudget <- data.frame(
    name = "A",
    Age = t.full[keep],
    x = x.full[keep])
  AOF <- aof(
    name = TimeBudget$name,
    Age = TimeBudget$Age,
    x = TimeBudget$x)
  expect_equal(
    object = class(AOF),
    expected = "data.frame"
  )
})

test_that("aof return data.frame with n.obs < 5 and > 1",{
  mu1 <- 50
  mu2 <- 50
  rho1 = 0.5
  rho2 = 0.5
  n.obs <- 4
  sigma1 <- 0.1
  sigma2 <- sigma1
  SimTS <- function(n, mu, rho, sigma){
    X.standard <- arima.sim(n, model = list(ar = rho))
    X.standard/sd(X.standard)*sigma + mu
  }
  t.full <- 0:50
  t.break <- 25
  x.full <- c(SimTS(t.break, mu1, rho1, sigma1),
              SimTS(max(t.full)-t.break+1, mu2, rho2, sigma2))
  keep <- sort(sample(1:length(x.full), n.obs))
  TimeBudget <- data.frame(
    name = "A",
    Age = t.full[keep],
    x = x.full[keep])
  AOF <- aof(
    name = TimeBudget$name,
    Age = TimeBudget$Age,
    x = TimeBudget$x)
  expect_equal(
    object = class(AOF),
    expected = "data.frame"
  )
})
