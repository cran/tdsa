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
# (so the step size is 0.03).
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

