#' aof
#'
#' Ontogenetic shifts in central-place foraging insects
#'
#' @details A breakpoint-based method to detect ontogenetic shifts in
#' univariate time-activity budget series of central-place foraging insects.
#' The method finds a single change point according to the likelihood function.
#' The method was developed with honey bees in order to detect the
#' Age at Onset of Foraging (AOF), but can be used for the detection of other
#' ontogenetic shifts in other central-place foraging insects. For more
#' details, see Requier et al. (2020) Measuring ontogenetic shifts in
#' central-place foraging insects: a case study with honey bees. Journal of
#' Animal Ecology.
#'
#' @param bee The identity of the insect (e.g. a bee) as factor
#'   (e.g. "A00103C00020C301", "bee1").
#' @param Age The age of the insect in day as numeric (e.g. 1, 4, 32).
#' @param x The daily activity of the insect at a given age as a numeric value,
#'   for instance (i) the number of the trips per day, (ii) the duration of
#'   the trips per day, or (iii) the time of the trips per day.
#' @return A data.frame containing the aof results (one row per insect).
#' @examples
#' require("bcpa")
#' # Exemple with simulated data:
#' # A study case with no change simulated (mu1=mu2)
#' mu1 <- 50  # behavioural value at stage 1
#' mu2 <- 50
#' rho1 <- 0.5 # interval frequency at stage 1
#' rho2 <- rho1
#' # Low number of individuals (N, n.obs) and low variance (V, sigma)
#' # create time series from 0 to 50 with a behavioural change at 25
#' t.full <- 0:50
#' t.break <- 25
#' n.obs <- 5 # no. observations randomly selected in the time series: 5 to 45
#' sigma1 <- 0.1 # variance: 0.1 to 3
#' sigma2 <- sigma1
#' SimTS <- function(n, mu, rho, sigma){
#'   X.standard <- arima.sim(n, model = list(ar = rho))
#'   X.standard/sd(X.standard)*sigma + mu
#' }
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

  # For each insect
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

      # If the number of trip day is equal to one,
      # informing that "no data" is available for computing the function of
      # behavioural change
      if (length(x) <= 1){
        AOF$AOF[AOF$bee == bee & AOF$parameter == parameter] <- "no.data"
      }else{

        # If not (>1), and if the number of trip day is less than five,
        # informing that "not enough data" are available for computing the
        # function of behavioural change
        if (length(x) < 5){
          AOF$AOF[AOF$bee == bee &
            AOF$parameter == parameter] <- "not.enough.data"
          AOF$behav.flightspan[AOF$bee == bee &
            AOF$parameter==parameter] <- mean(data.source$y)
          # Here is informed the average value of the flight activity
        }else{

          # for all cases with more than five days of trips, we can
          # computed the GetBestBreak function of the bcpa R-package
          # (Gurarie, E. (2013) bcpa: Behavioral Change Point Analysis of
          # Animal Movement; see also Gurarie et al. 2009, Ecology Letters)
          # We use GetBestBreak to detect a single behavioural change in time
          # series of  of (i) number, (ii) duration, and (iii) time of daily
          # trips.
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
            # the time series. See the GetModels function for methematical
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
