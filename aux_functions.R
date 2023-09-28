# jonashaslbeck@protonmail.com; March 2nd, 2023

# -------------------------------------------------
# -------- What is happening here? ----------------
# -------------------------------------------------

# Some auxiliary functions

# -------------------------------------------------
# -------- Estimate VAR ---------------------------
# -------------------------------------------------

EstimateVAR <- function(data, roundpoints = 4,
                        summaries = FALSE,
                        changescore = FALSE) {

  p <- ncol(data)
  Phi <- matrix(NA, p, p)
  alpha <- rep(NA, p)
  residuals <- matrix(NA, nrow(data)-1, p)
  summary_list<- list()

  for(i in 1:p) {

    # Estimates the change-score version of the var model instead
    if(changescore == TRUE){
      y = data[-1,i] - data[-nrow(data), i]
    } else {
      y <- data[-1, i]
    }


    # predicted cases
    X <- data[-nrow(data), ]

    coefs <- lm(y ~ X)
    Phi[i, ] <- coefs$coefficients[-1]
    alpha[i] <- coefs$coefficients[1]
    residuals[,i] <- coefs$residuals
    if(summaries) summary_list[[i]] <- summary(coefs)
  }

  Psi <- cov(residuals)
  mu <- as.vector(solve(diag(p)-Phi)%*%alpha)

  coefs <- list(round(alpha, roundpoints),
                round(mu,roundpoints),
                round(Phi, roundpoints),
                round(Psi,roundpoints))

  names(coefs) <- c("intercepts", "means","phi","psi")

  if(summaries){
    return(summary_list)
  } else {
    return(coefs)
  }

} # eoF



# -------------------------------------------------
# -------- Function to extract 5 Symptoms ---------
# -------------------------------------------------

# This has been coded up together with Oisin Ryan for this project:
# # https://github.com/ryanoisin/ComputationalTreatment/


getSymptoms <- function(dframe,
                        av_thresh = .5,
                        as_thresh = .4,
                        oldsim = FALSE){
  # in new simulation, avoidance is "V" rather than "AV" as in the original simulation
  if(isTRUE(oldsim)) AV <- "AV" else AV <- "V"

  plist <- PanicModel::detectPanic(dframe$AF)

  n_panic <- plist$n_panic

  # q1 ; skip "limited symptom" condition, only recode n_panic to none, 1-2, more than 2 but not once a day, >1 a day
  q1_panic <- 0
  if(0 < n_panic && n_panic <= 2 ) q1_panic <- 2
  if(n_panic > 2 && n_panic <= 7) q1_panic <- 3
  if(n_panic > 7) q1_panic <- 4
  # NOTE: no concept of limited symptom panic attack, so no "1" criteria
  # possibility: shoehorn in some "fear" threshold to make it 1 (fear above .25?)

  # q2: how distressing were the panic attacks?
  # assuming severity between 0 and 1
  q2_padistress <- NA
  q2_distress_cont <- 0
  if(n_panic == 0) q2_padistress <- 0
  if(n_panic >= 1){
    msev <- mean(plist$panic_stats[,"severity"])
    if(msev <= .25) q2_padistress <- 1
    if(msev > .25 && msev <= .5 ) q2_padistress <- 2
    if(msev > .5 && msev <= .75) q2_padistress <- 3
    if(.75 < msev) q2_padistress <- 4
    q2_distress_cont <- msev
  }


  # q3: worried or felt anxious about / fear about next panic attack?
  # OR: to me, this is related to either AV or just fear, so i will split these into two questions
  # to account for the fact that you have a recovery period, i need to omit panic attacks + recovery period
  pind <- plist$ind_label
  # find out when the indicator turns on or off
  tmp <- which(pind[-1] != pind[-length(pind)]) + 1

  # find end points of the panic attacks
  ends <- tmp[seq(2,length(tmp))]


  # count two hours from each end point if there was a panic attack
  if(n_panic > 0){
    # edge case; if there is a single panic attack that starts at the end of the window and doesnt end
    if(length(tmp) == 1){
      outside_pa <- seq(1:nrow(dframe))[-tmp]
    }else{
      ends + 60*2
      tmp2 <- as.numeric(sapply(ends, function(s){
        seq(s, s + 60*2,1)
      }))

      recovery <- rep(0,length(pind))
      recovery[tmp2] <- 1

      # here are the time poitns which are outside panic attacks & recovery periods
      # note: previous operation can add extra elements to recovery (if PA happens < 2 hours before end)
      outside_pa <- (pind ==0 & recovery[1:length(pind)] == 0)
    }
  } else{
    outside_pa <- seq(1:nrow(dframe))
  }
  fear_mean <- mean(dframe[outside_pa,"AF"])
  avoid_mean <- mean(dframe[outside_pa,AV]) # in new code,

  #
  q3_fear_cont <- fear_mean
  q3_fear <- NA
  # if(fear_mean <= .01) q3_fear <- 0
  #   if(fear_mean <= .25) q3_fear <- 1
  #   if(fear_mean > .25 && fear_mean <= .5 ) q3_fear <- 2
  #   if(fear_mean > .5 && fear_mean <= .75) q3_fear <- 3
  #   if(.75 < fear_mean) q3_fear <- 4
  if(fear_mean <= .001) q3_fear <- 0
  if(fear_mean > .001 && fear_mean <= .005) q3_fear <- 1
  if(fear_mean > .005 && fear_mean <= .01 ) q3_fear <- 2
  if(fear_mean > .01 && fear_mean <= .02) q3_fear <- 3
  if(.02 < fear_mean) q3_fear <- 4

  # q4: any places or situations you avoided or felt afraid of because of fear of panic?
  q4_avoid_cont <- avoid_mean
  q4_avoid <- NA
  if(avoid_mean <= .01) q4_avoid <- 0
  if(avoid_mean > .01 && avoid_mean <= .25) q4_avoid <- 1
  if(avoid_mean > .25 && avoid_mean <= .5 ) q4_avoid <- 2
  if(avoid_mean > .5 && avoid_mean <= .75) q4_avoid <- 3
  if(.75 < avoid_mean) q4_avoid <- 4

  q5_context_cont <- 1- mean(dframe$p_C)
  q5_context <- NA
  if( q5_context_cont <= .908)  q5_context <- 0
  if(.908 < q5_context_cont && q5_context_cont <= .910)  q5_context <- 1
  if(.910 < q5_context_cont && q5_context_cont <= .991)  q5_context <- 2
  if(.991 < q5_context_cont && q5_context_cont <= .992)  q5_context <- 3
  if(.992 <  q5_context_cont)  q5_context <- 4


  # q5 - avoid situations. p_C; range approx .01 to ...?
  # (range unknown - maybe relative to some baseline/ start value?)
  # same problem for fear and avoid; check min values ofhealthy (no AS)
  # vs max values of unhealthy (AS maxed out)

  sumscore <- q1_panic +  q2_padistress +  q3_fear +  q4_avoid + q5_context
  out = matrix(c(q1_panic, q2_padistress, q3_fear, q4_avoid,q5_context, sumscore,
                 n_panic,q2_distress_cont, q3_fear_cont, q4_avoid_cont, q5_context_cont),1, 11)
  colnames(out) <- c("q1_panic", "q2_padistress", "q3_fear", "q4_avoid","q5_context", "sumscore",
                     "n_panic","q2_distress_cont", "q3_fear_cont", "q4_avoid_cont", "q5_context_cont")

  return(out)

} # eoF


