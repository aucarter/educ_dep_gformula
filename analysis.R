### Setup
library(data.table); library(ggplot2); library(lme4); library(BMA)
rm(list = ls())

### Paths
data.path <- "data/educ_dep_data.csv"

### Code
## Prep data
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

time.varying.vars <- c("income", "socialSupport", "familyNear", "disabled", "pain", "depression")
time.fixed.vars <- setdiff(setdiff(names(dt), time.varying.vars), "person") # remove person as it is just an identifier, but we may want to investigate clustering in the future

# Add lags for time.varying.vars
for(var in time.varying.vars) {
    dt[, paste0(var, "_lag") := shift(get(var)), by = person]
}


## Estimate conditional probability models
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
binary.vars <- c("familyNear", "disabled", "pain", "depression")
alive.dt <- dt[dead == 0]
dt[disabled > 1, disabled := 1] # Find out what the other disabled values mean!!
dt[pain > 1, pain := 1] # Find out what the other pain values mean!!
model.fits <- lapply(time.varying.vars, function(var) {
    print(var)
    formula <- paste(var, "~", paste(c(time.fixed.vars, paste0(setdiff(time.varying.vars, var), "_lag")), collapse=" + "))
    if(var %in% binary.vars) {
        fit <- do.call("glm", list(as.formula(formula), data = as.name("alive.dt"), family = as.name("binomial")))
    } else {
        fit <- do.call("lm", list(as.formula(formula), data = as.name("alive.dt")))
    }
    coef.vec <- vectorize.fit(fit, dt)
    return(coef.vec)
})

coef.matrix <- do.call(rbind, model.fits)

# Fit mortality!
dead.fit <- glm(dead ~ time + age + sex + childSES1 + childSES2 + educationAttainment + income_lag + socialSupport_lag + familyNear_lag + disabled_lag + pain_lag + depression_lag, 
               data = dt, family = binomial(link = "logit"))
dead.vec <- vectorize.fit(dead.fit, dt)

# Combine into coefficient matrix
coef.matrix <- rbind(coef.matrix, dead.vec)
colnames(coef.matrix) <- c("Intercept", names(dt))
rownames(coef.matrix) <- c(time.varying.vars, "dead")

## Step 3 - sample with replacement from data
# Generate a baseline matrix
gen.baseline <- function(n, baseline.dt, intervene) {
    M.n <- n * nrow(baseline.dt)
    M.ids <- sample(baseline.dt$id, M.n, replace = T)
    if(intervene == 2) {
        baseline.dt[, educationAttainment := educationAttainment + 1]
    }
    dt.matrix <- as.matrix(baseline.dt)
    M <- cbind(dt.matrix[M.ids,])
    return(M)
}

# Update time varying predictors: time and age (assume 2-year intervals)
update.time <- function(M) {
    M[,"time"] <- M[,"time"] + 1
    M[, "age"] <- M[,"age"] + 2
    return(M)
}

# Update lags on time varying predictors
update.lags <- function( M) {
    for(var in time.varying.vars) {
        M[, paste0(var, "_lag")] <- M[, var]
    }
    return(M)
}


# Predict log.odds, convert to probability, and sample from Bernoulli
gen.draws <- function(var, M, coef.matrix) {
    var.idx <- which(colnames(M) == var)
    zero.idx <- which(M[, var.idx] == 0)
    if(length(zero.idx) > 0) {
        coef.vec <- as.matrix(coef.matrix[variable == var, 1:(ncol(coef.matrix) - 1)])[1,]
        temp.M <- cbind(rep(1, length(zero.idx)), M[zero.idx,,drop = F])
        colnames(temp.M)[1] <- "(Intercept)"
        # Order variables
        coef.vec <- coef.vec[which(colnames(temp.M) == names(coef.vec))]
        log.odds <- temp.M %*% coef.vec
        odds <- exp(log.odds)
        prob <- odds / (1 + odds)
        draws <- unlist(lapply(prob, rbinom, n = 1, size = 1))
        M[zero.idx, var.idx] <- draws
    }
    return(M)
}



simulate <- function(M, intervene, coef.matrix) {
    total <- nrow(M)
    var.list <- rownames(coef.matrix)
    # Loop through each day
    for(i in 1:max(dt$day)) {
        # Update lags and time variables
        M <- update.lags(M)
        M <- update.time(M)
        # Generate a predicted value for each variable in order
        for(var in var.list) {
            var.idx <- which(colnames(M) == var)
            M <- gen.draws(var, M, coef.matrix)
        }
        # Pull out dead add to out table
        dead.dt <- as.data.table(M[M[, "dead"] == 1, ,drop = F])
        if(nrow(dead.dt) > 0) {
            out.dt <- rbind(out.dt, dead.dt)
            M <- M[!(M[, "dead"] == 1), ,drop = F]
        }
        if(nrow(M) == 0) {
            break
        }
        print(paste0("Day ", i, " of ", max(dt$time), "; P(alive) = ", round(nrow(M) / total, 2)))
    }
    out.dt <- rbind(out.dt, as.data.table(M))
    return(out.dt)
}

baseline <- dt[,min(time), by = person]

# Natural course
n.draws <- 100
M.nat <- gen.baseline(n.draws, dt[time == 1], intervene = 0)
sim.nat <- simulate(M.nat, intervene = 0, coef.matrix)

# Intervention
M.educ<- gen.baseline(n.draws, dt[time == 1], intervene = 1)
sim.educ <- simulate(M.educ, intervene = 1, coef.matrix)


## Step 6 - concatentate intervetion data sets and run Cox model
sim.nat[, scenario := "Natural Course"]
sim.educ[, scenario := "Education increase"]
data.dt[, scenario := "Data"]
sim.dt <- rbindlist(list(sim.nat, sim.educ), use.names = T, fill = T)
out.dt <- sim.dt[, .(time, dead, scenario)]
write.csv(out.dt, "results/sim_results.csv", row.names = F)

# Make kaplan meier survival plots
surv.fit <- survfit(Surv(time, d) ~ scenario, data=sim.dt, conf.type = "log-log")
plot(surv.fit,  col=c("blue", "red", "green", "purple", "orange"),lty=c("dashed", "solid", "solid", "solid", "solid"), xlab="Days",
     ylab="Survival Probability", main="Kaplan Meier survial curve estimates")

legend("topright",
       c("Observed", "Natural course", "No Acute GvHD", "No Chronic GvHD", "No GvHD"),
       col=c("blue", "red", "green", "purple", "orange"),
       lty=c("dashed", "solid", "solid", "solid", "solid"), lwd=rep(3, 2), cex=0.9)



## Step 7: Cox Proportional Hazard Model and Uncertainty intervals
print.hazard <- function(nat.dt, int.dt) {
    coxph(Surv(day, d) ~ gvhd, data = all.p.day.dt)
    rand <- sample(1:137000, 137000, replace = T)
    boot.a <- unlist(lapply(1:137, function(i) {
        idx <- rand[(i*137+1):((i+1)*137)]
        cox.fit <- coxph(Surv(day, d) ~ scenario, data = rbind(nat.dt[idx], int.dt[idx]))
        return(-1*summary(cox.fit)$coefficients[,"coef"])
    }))
    mean <- round(exp(mean(boot.a)), 2)
    se <- 1.96 * sqrt(sd(exp(boot.a)) / 137)
    lower <- round(mean - se, 2)
    upper <- round(mean + se, 2)
    
    # lower <- round(exp(quantile(boot.a, 0.05)), 2)
    # upper <- round(exp(quantile(boot.a, 0.95)), 2)
    print(paste0(mean, " (", lower, ", ", upper, ")"))
}

# Education intervention
print.hazard(sim.nat, sim.educ)


