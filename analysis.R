################################################################################
## Purpose: Use g-formula to investigate the impact of education on depression
##          and mortality in the presence of time-varying confounders
## Author: Austin Carter, aucarter@uw.edu
## To Do: 
## - Direct vs indirect [ ]
## - Add lost-to-follow-up [X]
## - Change time to be age difference [X]
## - Investigate age-specific mortality rate (prevalence of depression) [ ]
## - Incidence of depression (is the onset delayed) median? [ ]
## - Add a tracker for depression onset in simulation [ ]
## - Sample from full uncertainty [ ]
## - Model non-binary in log space [ ]
## - Add nicely formatted tables [ ]
## - Look into clustering (multiple observations per individual) [ ]
##
## Notes:
## - No one dies before 10 years of follow-up
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

# Add a loss-to-follow-up indicator
dt[, ltfu := ifelse(dead == 0 & time == max(time), 1, 0)]

# Add years of follow-up and remove times
dt[, years := age - min(age), by = .(person)]
dt[, time := NULL]

# Reorder for viewing
dt <- dt[order(person, years)]

# Try removing income because of correlation with education
dt[, income := NULL]

time.varying.vars <- intersect(names(dt), 
    c("income", "socialSupport", "familyNear", "disabled", "pain", "depression"))

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
binary.vars <- c("familyNear", "disabled", "pain", "depression", "ltfu", "dead")
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

# Fit loss-to-follow-up
ltfu.fit <- glm(ltfu ~ years + age + sex + childSES1 + childSES2 + 
                    educationAttainment + socialSupport_lag + 
                    familyNear_lag + disabled_lag + pain_lag + depression_lag, 
                data = dt, family = binomial(link = "logit"))
ltfu.vec <- vectorize.fit(ltfu.fit, dt)

# Fit mortality
dead.fit <- glm(dead ~ years + age + sex + childSES1 + childSES2 + 
                    educationAttainment + socialSupport_lag + 
                    familyNear_lag + disabled_lag + pain_lag + depression_lag, 
               data = dt, family = binomial(link = "logit"))
dead.vec <- vectorize.fit(dead.fit, dt)

# Combine into coefficient matrix
coef.matrix <- rbind(coef.matrix, ltfu.vec)
coef.matrix <- rbind(coef.matrix, dead.vec)
colnames(coef.matrix) <- c("Intercept", names(dt))
rownames(coef.matrix) <- c(time.varying.vars, "ltfu", "dead")

## Step 3: Sample with replacement from data and simulate
# Functions
gen.baseline <- function(n, baseline.dt, intervene) {
    # Generate a baseline matrix
    M.n <- n * nrow(baseline.dt)
    M.ids <- sample(1:nrow(baseline.dt), M.n, replace = T)
    if(intervene == 1) {
        # baseline.dt[, educationAttainment := educationAttainment + 1]
        # baseline.dt[educationAttainment > 5, educationAttainment := 5]
        baseline.dt[, educationAttainment := 5]
    }
    dt.matrix <- as.matrix(baseline.dt)
    M <- cbind(dt.matrix[M.ids,])
    return(M)
}

update.time <- function(M) {
    # Update time varying predictors: years and age (assume 2-year intervals)
    M[,"years"] <- M[,"years"] + 1
    M[, "age"] <- M[,"age"] + 1
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
    zero.idx <- which(M[, var.idx] == 0)
    if(length(zero.idx) > 0 | !(var %in% binary.vars)) {
        coef.vec <- as.vector(coef.matrix[var,])
        if(var %in% binary.vars) {
            temp.M <- cbind(rep(1, length(zero.idx)), M[zero.idx,,drop = F])
        } else {
            temp.M <- cbind(rep(1, nrow(M)), M)
        }
        colnames(temp.M)[1] <- "(Intercept)"
        # Order variables
        # coef.vec <- coef.vec[which(colnames(temp.M) == names(coef.vec))]
        pred <- temp.M %*% coef.vec
        if(var %in% binary.vars) {
            odds <- exp(pred)
            prob <- odds / (1 + odds)
            val <- unlist(lapply(prob, rbinom, n = 1, size = 1))
        } else {
            # log.val <- rnorm(1, pred)
            # val <- exp(log.val)
            val <- pred
        }
        if(var %in% binary.vars) {
            M[zero.idx, var.idx] <- val
        } else {
            M[, var.idx] <- val
        }
    }
    return(M)
}

simulate <- function(M, intervene, coef.matrix, n = 40) {
    # Simulate n time periods (40 years)
    total <- nrow(M)
    out.dt <- data.table()
    var.list <- rownames(coef.matrix)
    # Loop through each time period
    for(i in 1:n) {
        # Update lags and time variables
        M <- update.lags(M)
        M <- update.time(M)
        # Generate a predicted value for each variable in order
        for(var in var.list) {
            M <- gen.draws(var, M, coef.matrix)
        }
        # Pull out dead and add to out table
        dead.dt <- as.data.table(M[M[, "dead"] == 1 | M[, "ltfu"] == 1, ,drop = F])
        if(nrow(dead.dt) > 0) {
            out.dt <- rbind(out.dt, dead.dt)
            M <- M[!(M[, "dead"] == 1), ,drop = F]
        }
        if(nrow(M) == 0) {
            break
        }
        print(paste0("Year ", i, " of ", n, "; P(alive) = ", 
                     round(nrow(M) / total, 2)))
    }
    out.dt <- rbind(out.dt, as.data.table(M))
    return(out.dt)
}

# Natural course
n.draws <- 10
M.nat <- gen.baseline(n.draws, dt[years == 0], intervene = 0)
sim.nat <- simulate(M.nat, intervene = 0, coef.matrix)

# Intervention 
M.educ<- gen.baseline(n.draws, dt[years == 0], intervene = 1)
sim.educ <- simulate(M.educ, intervene = 1, coef.matrix)

## Step 4: Concatentate simulation data sets, write results, and
##         plot kaplan meier survival plots
sim.nat[, scenario := "Natural Course"]
sim.educ[, scenario := "Education increase"]
dt[, scenario := "Data"]
sim.dt <- rbindlist(list(sim.nat, sim.educ, dt), use.names = T, fill = T)
out.dt <- sim.dt[, .(years, dead, depression, scenario)]
# write.csv(out.dt, "results/sim_results.csv", row.names = F)

# Make kaplan meier survival plots
surv.fit <- survfit(Surv(age, depression) ~ scenario , data=sim.dt, 
                    conf.type = "log-log")
plot(surv.fit,  col=c("blue", "red", "green"),lty=c("dashed", "solid", "solid"), 
     xlab="Age", ylab="Probability Not Depressed", 
     main="Kaplan Meier depression curve estimates")

legend("bottomleft",
       c("Observed", "Education Intervention", "Natural course"),
       col=c("blue", "red", "green"),
       lty=c("dashed", "solid", "solid"), lwd=rep(3, 2), cex=0.9)


## Calculate summary statistics
sim.dt[, age_group := age - age%%10]
sim.dt[, total := .N, by = .(scenario, age_group)]
summ.dt <- sim.dt[scenario != "Data", lapply(.SD, sum), by = .(age_group, total, scenario), .SDcols = c("depression", "dead")]
summ.dt[, c("dep_rate", "dead_rate") := .(depression / total, dead / total)]

gg <- ggplot(summ.dt) + geom_line(aes(x = age_group, y = dep_rate, color = scenario)) +
    ylim(0, max(summ.dt$dep_rate)) + ggtitle("Depression rate by age for each scenario") +
    theme_classic()
gg

gg <- ggplot(summ.dt) + geom_line(aes(x = age_group, y = dead_rate, color = scenario)) +
    ylim(0, max(summ.dt$dead_rate)) + ggtitle("Mortality rate by age for each scenario") +
    theme_classic()
gg
## End
