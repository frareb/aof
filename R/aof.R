#' aof
#'
#' Ontogenetic shifts in central-place foraging insects
#'
#' @details A breakpoint-based method to detect ontogenetic shifts in
#' univariate time-activity budget series of central-place foraging insects.
#' The method finds a single breakpoint according to the likelihood function.
#' The method was developped with honey bees in order to detect the
#' Age at Onset of Foraging (AOF), but can be used for the detection of other
#' ontongenetic shifts in other central-place foraging insects. For more
#' details, see Requier et al. (2020) Measuring ontogenetic shifts in
#' central-place foraging insects. Journal of Animal Ecology.
#'
#' @param bee The identity of the bee as factor (e.g. "A00103C00020C301", "bee1").
#' @param Age The age of the bee in day as numeric (e.g. 1, 4, 32).
#' @param x The daily activity of the bee at a given age as a numeric value,
#'   for instance (i) the number of the trips per day, (ii) the duration of
#'   the trips per day, or (iii) the time of the trips per day.
#' @return A data.frame with one row containing the aof results.
#' @examples
#' require("bcpa")
#' # Exemple with simulated data:
#' # No change simulated
#' mu1 <- 50  # behavioral change: 25 or 50
#' mu2 <- 50
#' rho1 = 0.5 # intrerval frequency
#' rho2 = 0.5
#' # N low and V low
#' # A single changepoint in a time-series where the parameters change at
#' # some unknown timepoints t* is done by simply sweeping all possible
#' # breaks, and finding the most-likely changepoint according to the
#' # likelihood. This is illustrated below:
#' n.obs <- 5 # no. observations: 5 to 45
#' sigma1 <- 0.1 # variance: 0.1 to 3
#' sigma2 <- sigma1
#' SimTS <- function(n, mu, rho, sigma){
#'   X.standard <- arima.sim(n, model = list(ar = rho))
#'   X.standard/sd(X.standard)*sigma + mu
#' }
#' # create time series with break at 25
#' t.full <- 0:50
#' t.break <- 25
#' x.full <- c(SimTS(t.break, mu1, rho1, sigma1),
#'   SimTS(max(t.full)-t.break+1, mu2, rho2, sigma2))
#' # subsample of observations (n defined above) and estimate
#' keep <- sort(sample(1:length(x.full), n.obs))
#' TimeBudget <- data.frame(
#'   bee = "A",
#'   Age = t.full[keep],
#'   x = x.full[keep])
#' # Running the algorithm
#' AOF <- aof(
#'   bee = TimeBudget$bee,
#'   Age = TimeBudget$Age,
#'   x = TimeBudget$x)
#' print(AOF)
#' # see vignette for more examples
#' @export
aof <- function(bee, Age, x){
  TimeBudget <- data.frame(
    bee = bee,
    Age = Age,
    x = x)

  # Preparation of the output table
  # --------------------------------------------
  AOF <- data.frame(
    bee = rep(
      x = unique(TimeBudget$bee),
      times = length(colnames(TimeBudget)[3:ncol(TimeBudget)])),
    parameter = sort(
      rep(
        x = colnames(TimeBudget)[3:ncol(TimeBudget)],
        times = length(unique(TimeBudget$bee)))),
    AOF = NA, # will be computed later
    behav.flightspan = NA, # will be computed later
    behav.learning = NA, # will be computed later
    behav.foraging = NA) # will be computed later
  # --------------------------------------------

  # start of the algorithm

  # Initial work time of the function
  # --------------------------------------------
  iter <- length(unique(TimeBudget$bee))
  progress.time <- utils::txtProgressBar(min = 0, max = iter, style = 3)
  prog.i <- 1
  # --------------------------------------------

  # For each bee
  for(bee in unique(TimeBudget$bee)){
    select <- (TimeBudget$bee == bee)
    dataselect <- TimeBudget[select, ]
    dataselect <- dataselect[order(dataselect$Age), ]

    # Work time progress of the function, updated after each bee proceed
    # --------------------------------------------
    utils::setTxtProgressBar(progress.time, prog.i)
    prog.i <- prog.i + 1
    # --------------------------------------------

    # For each parameter (i.e. number, duration and time of the trips per day)
    for(parameter in colnames(TimeBudget)[3:ncol(TimeBudget)]){
      x <- dataselect[, 2] #Age
      y <- dataselect[, parameter]
      data.source <- data.frame(x, y)

      # If the number of trip day is egal to one,
      # informing that "no data" is available for computing the funciton of
      # behavioral change
      if (length(x) <= 1){
        AOF$AOF[AOF$bee == bee & AOF$parameter == parameter] <- "no.data"
      }else{

        # If not (>1), and if the number of trip day is less than five,
        # informing that "not enough data" are available for computing the
        # funciton of behavioral change
        if (length(x) < 5){
          AOF$AOF[AOF$bee == bee &
            AOF$parameter == parameter] <- "not.enough.data"
          AOF$behav.flightspan[AOF$bee == bee &
            AOF$parameter==parameter] <- mean(data.source$y)
          # Here is informed the average value of the flight activity
        }else{

          # for all the cases with more than five days of trips, we can
          # computed the BCPA function (Gurarie, E. (2013) bcpa: Behavioral
          # Change Point Analysis of Animal Movement; see also Gurarie et al.
          # 2009, Ecology Letters)
          # We use the GetBestBreak function of the bcpa R-package to detect
          # the change point in the time series of (i) number, (ii) duration
          # and (iii) time of daily trips

          # If error in the model, we doesn't save the change point and inform
          # it as "undetected"
          if(class(tryCatch(
            bcpa::GetBestBreak(y, x,tau=FALSE),
            error = function(e) e))[1] == "simpleError"){
            AOF$AOF[AOF$bee == bee &
              AOF$parameter == parameter] <- "undetected"
            AOF$behav.flightspan[AOF$bee == bee &
              AOF$parameter == parameter] <- mean(data.source$y)
          }else{

            # We then tested the parcimony of the change point prediction
            # using model selection with delta BIC > 2 between the simple
            # model (without change point) and all the six other models
            # predicting scenarii changes in the 3 differents parameters of
            # the time serie. See the GetModels function for methematical
            # details of each model this can be easily obtained using the
            # code '?GetModels'
            BB <- bcpa::GetBestBreak(y, x, tau = FALSE)
            BICsimple.model <- bcpa::GetModels(y, x, BB[1], tau = FALSE)[1, 3]
            BICmin.change.point.models <- min(
              bcpa::GetModels(y, x, BB[1], tau = FALSE)[2:8, 3])

            # If the delta BIC > 2, we save the change point (called AOF)
            # and compute the behaviour before and after it
            mod <- bcpa::GetModels(y, x, BB[1], tau = FALSE)
            mod[, 3] <- ifelse(mod[,3] %in% c("Inf", "-Inf"), NA, mod[, 3])
            if(diff(c(min(mod[2:8, 3], na.rm = T),mod[1,3])) >= 2){
              AOF$AOF[AOF$bee == bee &
                AOF$parameter == parameter] <- BB[2]
              AOF$behav.flightspan[AOF$bee == bee &
                AOF$parameter == parameter] <- mean(data.source$y)
              AOF$behav.learning[AOF$bee == bee &
                AOF$parameter == parameter] <- mean(
                  subset(data.source, data.source$x <= BB[2])$y)
              AOF$behav.foraging[AOF$bee == bee &
                AOF$parameter == parameter] <- mean(
                  subset(data.source, data.source$x > BB[2])$y)
            }else{
              # If the delta BIC < 2, we doesn't save the change point and
              # inform it as "undetected"
              AOF$AOF[AOF$bee == bee &
                AOF$parameter == parameter] <- "undetected"
              AOF$behav.flightspan[AOF$bee == bee &
                AOF$parameter == parameter] <- mean(data.source$y)
            }
          }
        }
      }
    }
  }

  # Work time end of the function
  #-------------------------------------
  close(progress.time)
  #-------------------------------------
  return(AOF)
}



# plot.aof <- function(AOF)
# graphics::plot(
#   x = TimeBudget$Age,
#   y = TimeBudget$x,
#   las = 1,
#   xlab = "Age (day)", ylab = "Time-activity budget",
#   col = NA)
# # graphics::abline(v = t.break, col = 1, lwd = 1, lty = 3)
# goodbreak1 = max(TimeBudget$Age[TimeBudget$Age < 25])
# goodbreak2 = min(TimeBudget$Age[TimeBudget$Age >= 25])
# # graphics::polygon(c(goodbreak1 - 0.5, goodbreak2 + 0.5,
# #           goodbreak2 + 0.5, goodbreak1 - 0.5),
# #         c(0, 0, 100, 100), col = "gray50", border = NA)
# graphics::points(
#   x = TimeBudget$Age,
#   y = TimeBudget$x,
#   pch = 21,
#   col = "gray20",
#   bg = "gray80")
# graphics::abline(v = AOF$AOF, col = "black", lwd = 2, lty = 3)
# graphics::lines(
#   c(min(TimeBudget$Age), AOF$AOF),
#   c(AOF$behav.learning, AOF$behav.learning),
#   col = "black", lwd = 2)
# graphics::lines(
#   c(AOF$AOF, max(TimeBudget$Age)),
#   c(AOF$behav.foraging,AOF$behav.foraging),
#   col = "black", lwd = 2)








