## -----------------------------------------------------------------------------
library(tdsa)

## -----------------------------------------------------------------------------
# Parameter values for the dynamic equations.
parms = list(
  b = 1,                                          # Per-capita birth rate.
  a = 0.1,                                        # Competition coefficient.
  mu = function(t){0.5 + 1/(1 + exp((t-10)/2))},  # Per-capita loss rate.
  sigma = 0.2                                     # Immigration rate.
)

# Function that returns the dynamic equations.
dynamic_fn = function(t, y, parms){
  b = parms[["b"]]
  a = parms[["a"]]
  sigma = parms[["sigma"]]
  mu = parms[["mu"]](t)
  
  dy = b*y*(1- a*y) - mu*y + sigma
  return( list(dy) )
}

# Initial conditions.
y_0 = 0.37  # Approximate steady-state population before restoration efforts.

# Function that returns the reward integrand.
reward_fn = function(t, y){
  w = 1  # Per-capita rate at which the ecosystem service is provided.
  return( w * y )
}

# Function that returns the terminal payoff.
terminal_fn = function(y){
  v = 1.74  # Ascribed value per individual at the end of the period.
  return( v * y )
}

# Time steps over management period. We discretise it into 1001 time steps
# (so the step size is 0.02).
times = seq(0, 30, length.out=1001)

## ---- results="hide"----------------------------------------------------------
state_sens_out = state_sens(
  model_type = "continuous",
  dynamic_fn = dynamic_fn,
  parms = parms,
  reward_fn = reward_fn,
  terminal_fn = terminal_fn,
  y_0 = y_0,
  times = times,
  verbose = FALSE
)

## -----------------------------------------------------------------------------
str(state_sens_out)

## ---- fig.dim=c(12,4)---------------------------------------------------------
# Set graphical parameters.
par(mfrow=c(1,3), cex=1)
par(mar=c(3.2,3.2,2,2), mgp=c(2,0.7,0), cex.lab=1.2)

# Plot the per-capita unregulated birth and loss rates.
plot(times, parms[["mu"]](times), type="l", lwd=2,
     xlab="Time (year)", ylab="Demographic rate (/year)")
abline(h=parms[["b"]], col="red", lwd=2)
legend("topright", col=c("red", "black"), lwd=2, bty="n",
       legend=c("Birth rate", "Loss rate"))

# Plot the population size.
plot(times, state_sens_out[["state"]][,1], type="l", lwd=2,
     xlab="Time (year)", ylab="Population size y")

# Plot the time-dependent state sensitivity. Peaks at around t=10, which is
# roughly when mu and b intersects, so the population has just become
# self-sustaining.
plot(times, state_sens_out[["tdss"]][,1], type="l", lwd=2,
     xlab="Time (year)", ylab="State sensitivity of y")

## ---- results="hide"----------------------------------------------------------
parm_sens_out = parm_sens(state_sens_out = state_sens_out,
                          verbose = FALSE)

## -----------------------------------------------------------------------------
str(parm_sens_out)

## ---- fig.dim=c(4,4), out.width="33%"-----------------------------------------
# Set graphical parameters.
par(mar=c(3.2,3.2,2,2), mgp=c(2,0.7,0), cex.lab=1.2)

# Plot the parameter sensitivity of sigma.
plot(times, parm_sens_out[["tdps"]][["sigma"]][,1], type="l", lwd=2,
     xlab="Time (year)", ylab="Param. sensitivity of sigma")

## -----------------------------------------------------------------------------
parms = list(
  mu = c(1, 1, 1, 1, 1),  # Mortality rate of a susceptible individual.
  b = matrix(c(1.976537, 1.976537, 0.000000, 0.000000, 0.000000,
               1.976537, 1.976537, 1.976537, 0.000000, 0.000000,
               0.000000, 1.976537, 1.976537, 1.976537, 0.000000,
               0.000000, 0.000000, 1.976537, 1.976537, 1.976537,
               0.000000, 0.000000, 0.000000, 1.976537, 1.976537),
             nrow=5, byrow=T)  # Matrix of transmission coefficients; values chosen so that disease R_0 = 0.9. 
)

parms2 = list(
  B = c(5, 5, 5, 5, 1.02),  # Unregulated per-capita birth rate; low for Species 5.
  a = c(0.8, 0.8, 0.8, 0.8, 0.02),  # Intraspecific competition; chosen so all species have disease-free carrying capacity = 1.
  nu = c(0, 0, 0, 0, 5),  # Disease-induced mortality of an infected individual; only nonzero for Species 5.
  gamma = c(5, 5, 5, 5, 0),  # Recovery rate of an infected individual; zero for species 5.
  sigma = c(0.2, 0, 0, 0, 0)  # Per-capita exogenous spillover rate; only nonzero for Species 1.
)

## -----------------------------------------------------------------------------
dynamic_fn = function(t, y, parms, parms2){
  # To make the lines below easier to read, we "extract" each coefficient from the parameter lists.
  mu = parms[["mu"]]
  b = parms[["b"]]
  
  B = parms2[["B"]]
  a = parms2[["a"]]
  nu = parms2[["nu"]]
  gamma = parms2[["gamma"]]
  sigma = parms2[["sigma"]]

 # To make the lines below easier to read, we "extract" the susceptible and infected parts from the state vector.
  SS = y[1:5 ]
  II = y[6:10]
  
  # Calculate the species population size.
  NN = SS + II
  
  # RHS of the dynamic equations.
  dSS = B * NN * (1 - a*NN) - SS * (sigma + b%*%II) - mu * SS + gamma * II
  dII = SS * (sigma + b%*%II) - (mu + nu + gamma) * II
  
  return( list( c(dSS, dII) ) )
}

## -----------------------------------------------------------------------------
reward_fn = function(t, y){
  # Parameters.
  W = c(0, 0, 0, 0, 1)  # Per-capita rate of contribution to ecosystem service; only nonzero for Species 5.

  # Split the state vector.
  SS = y[1:5]
  II = y[6:10]
  
  # Return the reward integrand.
  NN = SS + II
  return( sum(W * NN) )
}


terminal_fn = function(y){
  # Parameters.
  V = c(0, 0, 0, 0, 1)  # Per-capita terminal payoff; only nonzero for Species 5.

  # Split the state vector.
  SS = y[1:5]
  II = y[6:10]
  
  # Return the terminal payoff.
  NN = SS + II
  return( sum(V*NN) )
}

## -----------------------------------------------------------------------------
SS_0 = (1-parms[["mu"]]/parms2[["B"]])/parms2[["a"]]     # At carrying capacity.
II_0 = c(0, 0, 0, 0, 0)  # Disease-free.
y_0  = c(SS_0, II_0)

## -----------------------------------------------------------------------------
t_0 = 0
t_1 = 5
times = seq(from=t_0, to=t_1, length.out=1001)

## ---- results="hide"----------------------------------------------------------
state_sens_out = state_sens(
  model_type = "continuous",
  dynamic_fn = dynamic_fn,
  parms = parms,
  reward_fn = reward_fn,
  terminal_fn = terminal_fn,
  y_0 = y_0,
  times = times,
  dynamic_fn_arglist = list(parms2 = parms2),
  verbose = FALSE
)

## -----------------------------------------------------------------------------
str(state_sens_out)

## ---- fig.dim=c(12,6)---------------------------------------------------------
# Calculate the derived quantities to be plotted.

# Disease prevalence:
SS = state_sens_out[["state"]][, 1:5]  # Number of susceptible individuals.
II = state_sens_out[["state"]][, 6:10]  # Number of infected individuals.
NN = SS + II  # Species population size.
PP = II / NN  # Disease prevalence.

# State sensitivities:
lambda_SS = state_sens_out[["tdss"]][, 1:5]  # State sensitivities for S_j.
lambda_II = state_sens_out[["tdss"]][, 6:10]  # State sensitivities for I_j.
lambda_NN = (1 - PP) * lambda_SS + PP * lambda_II

# Generate the plots.
palette = c("#42049EFF", "#900DA4FF", "#CC4678FF", "#F1844BFF", "#FCCE25FF")  # Colour palette.

par(mar=c(5,6.7,4,1.7)+0.1, mfrow=c(1,2))  # Set graphical parameters.

# Plot PP.
plot(NA, xlim=c(0, 5), ylim=c(0.00067, 1), log="y", xaxs="i", yaxs="i", yaxt="n",
     main="Infection dynamics",
     xlab="Time (in units of lifespan)",
     ylab="Fraction of population infected (log scale)",
     cex.main=2, cex.lab=1.5, cex.axis=1.5)
axis(side=2, cex.axis=1.5, at=10^((-4):(0)), labels=c(0.0001, 0.001, 0.01, 0.1, 1))
for(i in 1:5){
  lines(times, PP[,i], lwd=3, col=palette[i])
}
legend("topright", legend=paste("Species",1:5), lwd=3, col=palette[1:5], bty="n")

# Plot lambda_NN.
plot(NA, xlim=c(0, 5), ylim=c(0, 0.11), xaxs="i", yaxs="i",
     main="Time-dep. state sensitivities",
     xlab="Time (in units of lifespan)",
     ylab=bquote(atop("Sensitivity of J to sudden decrease", "in species population "~(-lambda[N["j"]]))),
     cex.main=2, cex.lab=1.5, cex.axis=1.5)
for(i in (1:4)){
  lines(times, -lambda_NN[,i], lwd=3, col=palette[i])
}
legend("topright", legend=paste("Species",1:4), lwd=3, col=palette[1:4], bty="n")

## ---- results="hide"----------------------------------------------------------
parm_sens_out = parm_sens(state_sens_out = state_sens_out,
                          numDeriv_arglist = list(method.args=list(r=2)),
                          verbose = FALSE)

## -----------------------------------------------------------------------------
str(parm_sens_out)

## ---- fig.dim=c(12,6)---------------------------------------------------------
# Extract the list containing the parameter sensitivities.
tdps = parm_sens_out[["tdps"]]

# Parameter sensitivities for the mortality rate of a susceptible individual.
kappa_mu = tdps[["mu"]]

# Parameter sensitivities for the forward transmission rate b_{2,1}, b_{3,2}, etc.
# These are given by tdps[["b"]][,2,1], tdps[["b"]][,3,2], etc.
# A more systematic way to extract these elements is to use mapply.
kappa_b = mapply(function(i,j){tdps[["b"]][,i,j]}, 2:5, 1:4)



# Generate the plots.
palette = c("#42049EFF", "#900DA4FF", "#CC4678FF", "#F1844BFF", "#FCCE25FF")  # Colour palette.

par(mar=c(5,6.7,4,1.7)+0.1, mfrow=c(1,2))  # Set graphical parameters.

# Plot kappa_mu.
plot(NA, xlim=c(0, 5), ylim=c(0, 0.11), xaxs="i", yaxs="i",
     main="Time-dep. parm. sensitivities",
     xlab="Time (in units of lifespan)",
     ylab=bquote(atop("Sensitivity of J to brief increase", "in mortality of susceptibles "~(kappa[mu["j"]]))),
     cex.main=2, cex.lab=1.5, cex.axis=1.5)
for(i in 1:4){
  lines(times, kappa_mu[,i], lwd=3, col=palette[i])
}
legend("topright", legend=paste("Species",1:4), lwd=3, col=palette[1:4], bty="n")

# Plot kappa_b
plot(NA, xlim=c(0, 5), ylim=c(0, 0.077), xaxs="i", yaxs="i",
     main="Time-dep. parm. sensitivities",
     xlab="Time (in units of lifespan)",
     ylab=bquote(atop("Sensitivity of J to brief decrease", "in forward transmission "~(-kappa[b["j+1, j"]]))),
     cex.main=2, cex.lab=1.5, cex.axis=1.5)
for(i in 1:4){
  lines(times, -kappa_b[,i], lwd=3, col=palette[i])
}
legend("topright", legend=paste("Species", 1:4, "to", 2:5), lwd=3, col=palette[1:4], bty="n")

