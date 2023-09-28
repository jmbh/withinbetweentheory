# jonashaslbeck@protonmail.com; September 28th, 2023

# -------------------------------------------------
# -------- What is happening here? ----------------
# -------------------------------------------------

# Loading simulation results and plotting Figure 3

# -------------------------------------------------
# -------- Loading Packages -----------------------
# -------------------------------------------------

# library(devtools)
# install_github("jmbh/Panicmodel")
library(PanicModel)

source("aux_functions.R")

library(graphicalVAR)

library(RColorBrewer)
library(scales)
library(qgraph)

# -------------------------------------------------
# -------- Loading Simulated Data -----------------
# -------------------------------------------------

# The script for the simulation can be found here:
# https://github.com/ryanoisin/ComputationalTreatment/tree/main/Simulation

# The simulation output itself is available upon request. The reason
# is that we were not able to upload them due to their size (1.9GB)
# However, we provide the relevant part of the output as an RDS file
# which is the output of the below subsetting of the simulation results

# simDir <- ""

# Collect only baseline data from all 5xx people
l_base_pers <- list()
counter <- 1
for(i in 1:32) {

  batch_i <- readRDS(paste0(simDir, "/IntSim6_Iter", i, ".RDS"))

  for(j in 1:16) {

    l_base_pers[[counter]] <- batch_i[[j]]$out_baseline
    counter <- counter + 1
  }
  print(i)
}

length(l_base_pers)

# saveRDS(l_base_pers, file="Files/BaselineData.RDS")

# Read from processed file
l_base_pers <- readRDS(file="Files/BaselineData.RDS")



# -------------------------------------------------
# -------- Panel 1: Data of one person ------------
# -------------------------------------------------

# ------ Pick one person -----

l_person_i <- l_base_pers[[1]] # we take Person 1
DataP1 <- l_person_i$outmat

head(DataP1)


# -------------------------------------------------
# -------- Panel 2: Within VAR: A+T ---------------
# -------------------------------------------------

# ----- On original time scale ------
data_VAR <- cbind(DataP1$A, DataP1$PT)
N <- nrow(data_VAR)
colnames(data_VAR) <- c("A", "PT")
out1 <- EstimateVAR(data_VAR)
qgraph(out1$phi, edge.labels=TRUE)


# -------------------------------------------------
# -------- Panel 3: Symptom Networks --------------
# -------------------------------------------------
# Get Symptoms from each person in the last week

n_weeks <- 4
week_length <- (nrow(DataP1)-1)/n_weeks

m_symptoms <- matrix(NA, 500, 5)

for(i in 1:500) {

  l_person_i <- l_base_pers[[i]]
  symptoms_lw_i <-  getSymptoms(l_person_i$outmat[(1+(4-1)*week_length):(week_length*4),])
  m_symptoms[i, ] <- symptoms_lw_i[1, 1:5]

  print(i)
} # end for: subj


# Estimate GGM
library(corpcor)
library(qgraph)
ggm_out <- EBICglasso(cor(m_symptoms), n=500)


# -------------------------------------------------
# -------- Panel 4: Scatter: Sumscore vs. init S --
# -------------------------------------------------

# Candidate: Linear effect of A->PT
m_mix <- as.data.frame(matrix(NA, 500, 9))
colnames(m_mix) <- c("sumscore", "VAR11", "VAR22", "VAR12", "VAR21",
                     "Mar11", "Mar22", "Mar12", "Mar21")
l_person_i$outmat$A

for(i in 1:500) {

  l_person_i <- l_base_pers[[i]]
  N <- nrow(l_person_i$outmat)

  # Symptom sum score last week
  symptoms_lw_i <- getSymptoms(l_person_i$outmat[(1+(4-1)*week_length):(week_length*4),])
  m_mix[i, 1] <- symptoms_lw_i[1, 6] #sumscore

  # VAR model
  data_VAR <- cbind(l_person_i$outmat$A, l_person_i$outmat$PT)
  out1 <- EstimateVAR(data_VAR)
  m_mix[i, 2] <- out1$phi[1,1]
  m_mix[i, 3] <- out1$phi[2,2]
  m_mix[i, 4] <- out1$phi[1,2]
  m_mix[i, 5] <- out1$phi[2,1]

  # Marginal auto/cross correlations
  m_mix[i, 6] <- cor(l_person_i$outmat$A[-N], l_person_i$outmat$A[-1])
  m_mix[i, 7] <- cor(l_person_i$outmat$PT[-N], l_person_i$outmat$PT[-1])
  m_mix[i, 8] <- cor(l_person_i$outmat$A[-N], l_person_i$outmat$PT[-1])
  m_mix[i, 9] <- cor(l_person_i$outmat$PT[-N], l_person_i$outmat$A[-1])

  print(i)
} # end for: subj



# -------------------------------------------------
# -------- Plotting it all For 4-panel Figure -----
# -------------------------------------------------


sc <- 0.7
pdf("Figures/Fig_Data_and_Phenomena_Sept28_23.pdf", width = 15*sc, height = 13*sc)


par(mfrow=c(2,2))

cols <- c("grey", "black", "lightblue", "tomato", "orange")
cex_mtext <- 1.2

# ----- Panel A: Data -----

# Some settings
cex_axis <- 0.8

day2inmin <- (60*24*1)
day <- 14 # days, 19+20
SimData_2days <- DataP1[(day2inmin*day):(day2inmin*(day+2)), ]

N_min <- length(SimData_2days$A)

par(mar=c(4.5,4,3,3))
plot.new()
plot.window(xlim = c(0, N_min),  ylim=c(-0.55, 1))
# box()
axis(2, las=2, seq(-0.5, 1, length=7), cex.axis=cex_axis)
axis(1, c(0, 12, 24, 36, 48), at=seq(0, day2inmin*2, length=5), cex.axis=cex_axis)
title(xlab="Hours", line=2.5)
mtext("A. Simulated Data of Person 1 (Days 13-14)", side=3, adj=-.1, line=.8, cex=cex_mtext)

lwd <- 1.75
lines(SimData_2days$S, col = cols[4],  lwd = lwd)
lines(SimData_2days$A, col = cols[1], lwd = lwd)
lines(SimData_2days$X, col = cols[5], lwd = lwd)
lines(SimData_2days$PT, col = cols[2],lwd = lwd)
lines(SimData_2days$E, col = cols[3], lwd = lwd, lty=2)

legend(N_min*.4, 0.8, legend = c("Arousal","Perceived Threat", "Escape Behavior", "Arousal Schema", "Escape Schema"),
       col = cols, bty = "n", cex = 0.95,
       text.col = cols)


# ----- Panel B: Within-person thing -----
lo <- rbind(c(.2, .8),
            c(.8, .2))
qgraph(out1$phi,
       layout = lo,
       edge.labels = TRUE,
       labels = c("A", "PT"),
       fade = F,
       # color = cols[1:2],
       vsize = 13,
       esize = 18,
       asize = 10,
       edge.label.cex = 1.5,
       theme = "colorblind",
       mar = rep(15,4))
mtext("            B. Within-person: VAR model of Person 1", side=3, adj=-.1, line=.8, cex=cex_mtext)

text(-.2, 1.2, "Arousal", cex=1.4, col="black")
text(-.2, -1.2, "Perceived Threat", cex=1.4, col="black")


# ----- Panel C: Between person symptom network -----
symptom_labels <- c("Panic", "Distress", "Fear", "Avoidance", "AvoidanceX")
symptom_labels_short <- c("P", "D", "F","AvB","AvC")
par(mar=c(3.5,4,2,3))
qgraph(ggm_out,
       labels = symptom_labels_short,
       # nodeNames = symptom_labels,
       # legend=TRUE,
       theme = "colorblind",
       vsize = 13,
       esize = 18,
       edge.label.cex = 1.5,
       edge.labels = TRUE,
       mar = rep(5,4))
mtext(".           C. Between-person: Symptom network", side=3, adj=-.1, line=-.2, cex=cex_mtext)

legend(-1.3,0, legend=c("P = Panic",
                           "D = Distress",
                           "F = Fear",
                           "AvB = Avoid Behavior",
                           "AvC = Avoid Contexts"), bty="n")

# ----- Panel D: Scatter mixed -----
par(mar=c(4.5,4,3,3))
plot.new()
plot.window(xlim=c(0.2, 0.8), ylim=c(0, 20))
axis(1)
axis(2, las=2)
points(m_mix$Mar12,  m_mix$sumscore, col=alpha("black", alpha=0.2), pch=20, cex=2)
title(xlab = expression(atop(phantom(0), "Perceived Threat"[t] %->% Arousal[t])))
title(ylab = "Symptom Sum Score")
lm_obj <- lm(m_mix$sumscore~m_mix$Mar12)
abline(lm_obj, lwd=2, lty=2)
mtext("D. Mixed: VAR parameter & Symptom Sum Score", side=3, adj=-.1, line=.8, cex=cex_mtext)


dev.off()








