################################################################################
## Purpose: CUse g-formula to investigate the impact of education on depression
##          and mortality in the presence of time-varying confounders
## Author: Austin Carter, aucarter@uw.edu
## To Do: 
## - Sample from full uncertainty
## - Add lost-to-follow-up
## - Model non-binary in log space
## - Add nicely formatted tables
## - Look into clustering (multiple observations per individual)
################################################################################

### Setup
library(data.table); library(ggplot2); library(lme4); library(BMA)
rm(list = ls())

### Paths
data.path <- "data/educ_dep_data.csv"

### Code
## Step 1: Prep data
dt <- fread(data.path)

# Turn everything into a number
for(var in names(dt)) {
    if(typeof(dt[[var]]) == "character") {
        dt[[var]] <- as.numeric(dt[[var]])
    }
}

# Remove NAs for non-death rows
dt <- rbind(na.omit(dt[dead == 0]), dt[dead == 1])

# Remove people without an initial observation
dt <- dt[!(person %in% setdiff(unique(dt$person), unique(dt[time == 1]$person)))]

# Reorder for viewing
dt <- dt[order(person, time)]

time.varying.vars <- c("income", "socialSupport", "familyNear", "disabled", 
                       "pain", "depression")
# remove person as it is just an identifier, but we may want to investigate 
# clustering in the future
time.fixed.vars <- setdiff(setdiff(names(dt), time.varying.vars), "person") 

# Add lags for time.varying.vars
for(var in time.varying.vars) {
    dt[, paste0(var, "_lag") := shift(get(var)), by = person]
}

## Step 2: Estimate conditional probability models
# Functions
vectorize.fit <- function(model.fit, dt){
    # Prep a vector of the coefficients from the model fits that can be used
    # as part of a matrix that includes all variables in the input data
    coefs <-  coef(summary(model.fit))[,1]
    vec <- unlist(lapply(names(dt), function(name) {
        if(name %in% names(coefs)) {
            coefs[name]
        } else {
            0
        }
    }))
    vec <- c(coefs[1], vec)
    return(vec)
}

# Fit time-varying models with all possible variables (except death!)
binary.vars <- c("familyNear", "disabled", "pain", "depression", "dead")
dt[disabled > 1, disabled := 1] # Find out what the other disabled values mean!!
dt[pain > 1, pain := 1] # Find out what the other pain values mean!!
alive.dt <- dt[dead == 0]

model.fits <- lapply(time.varying.vars, function(var) {
    print(var)
    formula <- paste(var, "~", 
                     paste(c(time.fixed.vars, 
                             paste0(
                                 setdiff(time.varying.vars, var), 
                                 "_lag")
                             ), 
                           collapse=" + ")
                     )
    if(var %in% binary.vars) {
        fit <- do.call("glm", list(as.formula(formula), 
                                   data = as.name("alive.dt"), 
                                   family = as.name("binomial")))
    } else {
        fit <- do.call("lm", list(as.formula(formula), 
                                  data = as.name("alive.dt")))
    }
    coef.vec <- vectorize.fit(fit, dt)
    return(coef.vec)
})

coef.matrix <- do.call(rbind, model.fits)

# Fit mortality
dead.fit <- glm(dead ~ time + age + sex + childSES1 + childSES2 + 
                    educationAttainment + income_lag + socialSupport_lag + 
                    familyNear_lag + disabled_lag + pain_lag + depression_lag, 
               data = dt, family = binomial(link = "logit"))
dead.vec <- vectorize.fit(dead.fit, dt)

# Combine into coefficient matrix
coef.matrix <- rbind(coef.matrix, dead.vec)
colnames(coef.matrix) <- c("Intercept", names(dt))
rownames(coef.matrix) <- c(time.varying.vars, "dead")

## Step 3: Sample with replacement from data and simulate
# Functions
gen.baseline <- function(n, baseline.dt, intervene) {
    # Generate a baseline matrix
    M.n <- n * nrow(baseline.dt)
    M.ids <- sample(1:nrow(baseline.dt), M.n, replace = T)
    if(intervene == 1) {
        baseline.dt[, educationAttainment := educationAttainment + 1]
    }
    dt.matrix <- as.matrix(baseline.dt)
    M <- cbind(dt.matrix[M.ids,])
    return(M)
}

update.time <- function(M) {
    # Update time varying predictors: time and age (assume 2-year intervals)
    M[,"time"] <- M[,"time"] + 1
    M[, "age"] <- M[,"age"] + 2
    return(M)
}

update.lags <- function( M) {
    # Update lags on time varying predictors
    for(var in time.varying.vars) {
        M[, paste0(var, "_lag")] <- M[, var]
    }
    return(M)
}

gen.draws <- function(var, M, coef.matrix) {
    # Predict 
    # Case 1: log odds, convert to probability, and sample from Bernoulli
    # Case 2: log, sample from a normal distribution, exponetiate
    var.idx <- which(colnames(M) == var)
    # zero.idx <- which(M[, var.idx] == 0)
    # if(length(zero.idx) > 0) {
        coef.vec <- as.vector(coef.matrix[var,])
        temp.M <- cbind(rep(1, nrow(M)), M)
        colnames(temp.M)[1] <- "(Intercept)"
        # Order variables
        # coef.vec <- coef.vec[which(colnames(temp.M) == names(coef.vec))]
        pred <- temp.M %*% coef.vec
        if(var %in% binary.vars) {
            odds <- exp(pred)
            prob <- odds / (1 + odds)
            val <- unlist(lapply(prob, rbinom, n = 1, size = 1))
        } else {
            log.val <- rnorm(1, pred)
            val <- exp(log.val)
        }
        M[, var.idx] <- val
    # }
    return(M)
}

simulate <- function(M, intervene, coef.matrix) {
    # Simulate 20 time periods (40 years)
    total <- nrow(M)
    out.dt <- data.table()
    var.list <- rownames(coef.matrix)
    # Loop through each time period
    for(i in 1:20) {
        # Update lags and time variables
        M <- update.lags(M)
        M <- update.time(M)
        # Generate a predicted value for each variable in order
        for(var in var.list) {
            var.idx <- which(colnames(M) == var)
            M <- gen.draws(var, M, coef.matrix)
        }
        # Pull out dead and add to out table
        dead.dt <- as.data.table(M[M[, "dead"] == 1, ,drop = F])
        if(nrow(dead.dt) > 0) {
            out.dt <- rbind(out.dt, dead.dt)
            M <- M[!(M[, "dead"] == 1), ,drop = F]
        }
        if(nrow(M) == 0) {
            break
        }
        print(paste0("Year ", 2*i, " of ", max(dt$time), "; P(alive) = ", 
                     round(nrow(M) / total, 2)))
    }
    out.dt <- rbind(out.dt, as.data.table(M))
    return(out.dt)
}

# Natural course
n.draws <- 1
M.nat <- gen.baseline(n.draws, dt[time == 1], intervene = 0)
sim.nat <- simulate(M.nat, intervene = 0, coef.matrix)

# Intervention 
M.educ<- gen.baseline(n.draws, dt[time == 1], intervene = 1)
sim.educ <- simulate(M.educ, intervene = 1, coef.matrix)

## Step 4: Concatentate simulation data sets, write results, and
##         plot kaplan meier survival plots
sim.nat[, scenario := "Natural Course"]
sim.educ[, scenario := "Education increase"]
dt[, scenario := "Data"]
sim.dt <- rbindlist(list(sim.nat, sim.educ, dt), use.names = T, fill = T)
out.dt <- sim.dt[, .(time, dead, scenario)]
write.csv(out.dt, "results/sim_results.csv", row.names = F)

# Make kaplan meier survival plots
surv.fit <- survfit(Surv(time, dead) ~ scenario, data=sim.dt, 
                    conf.type = "log-log")
plot(surv.fit,  col=c("blue", "red", "green"),lty=c("dashed", "solid", "solid"), 
     xlab="Years of Observation", ylab="Survival Probability", 
     main="Kaplan Meier survial curve estimates")

legend("topright",
       c("Observed", "Education Intervention", "Natural course"),
       col=c("blue", "red", "green"),
       lty=c("dashed", "solid", "solid"), lwd=rep(3, 2), cex=0.9)



