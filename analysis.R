### Setup
library(data.table); library(ggplot2);
rm(list = ls())

### Paths
data.path <- "data/educ_dep_data.csv"

### Code

## Read input data and prep person-day data

## Estimate conditional probability models
vectorize.fit <- function(model.fit, dt){
    coefs <-  coef(summary(model.fit))[,1]
    vec <- unlist(lapply(names(dt), function(name) {
        if(name %in% names(coefs)) {
            coefs[name]
        } else {
            0
        }
    }))
    vec <- c(coefs[1], vec)
}

get.pred.vars <- function(model.fit) {
    rownames(attr(summary(model.fit)$terms, "factors"))[1]
}

# Social support

# Family nearby

# Other comorbidities

# Depression


## Step 3 - sample with replacement from data

# Generate a baseline matrix
gen.baseline <- function(n, baseline.dt, intervene) {
    M.n <- n * nrow(baseline.dt)
    M.ids <- sample(baseline.dt$id, M.n, replace = T)
    days.no <- grep("daysno", names(baseline.dt), value = T)
    baseline.dt[, c("day", "d", "gvhd", "platnorm", "relapse", days.no) := 0]
    if(intervene == 2) {
        baseline.dt[, gvhd := 1]
    }
    dt.matrix <- as.matrix(baseline.dt)
    M <- cbind(dt.matrix[M.ids,])
    return(M)
}

# Update time varying predictors: day, daysq, daycu, daycurs1, dacurs2
update.time <- function(M) {
    day <- M[,"day"][1] + 1
    M[,"day"] <- M[,"day"] + 1
    M[, "daysq"] <- M[,"day"]**2
    M[, "daycu"] <- M[,"day"]**3
    M[, "daycurs1"] <- ((day>63)*((M[,"day"]-63)/63)**3)+((day>716)*((M[,"day"]-716)/63)**3)*(350.0-63) -((day>350)*((M[,"day"]-350)/63)**3)*(716-63)/(716-350)
    M[, "daycurs2"] <- ((day>168)*((M[,"day"]-168)/63)**3)+((day>716)*((M[,"day"]-716)/63)**3)*(350-168) -((day>350)*((M[,"day"]-350)/63)**3)*(716-168)/(716-350)
    return(M)
}

# Update lags on each iteration: platnormm1, gvhdm1, relapsem1
lag.list <- c("gvhd", "relapse", "platnorm")
update.lags <- function(var.list, M) {
    for(var in var.list) {
        M[, paste0(var, "m1")] <- M[, var]
    }
    return(M)
}


# Predict log.odds, convert to probability, and sample from Bernoulli
gen.draws <- function(var, M, coef.dt) {
    var.idx <- which(colnames(M) == var)
    zero.idx <- which(M[, var.idx] == 0)
    if(length(zero.idx) > 0) {
        coef.vec <- as.matrix(coef.dt[variable == var, 1:(ncol(coef.dt) - 1)])[1,]
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

# Update cumulative values on each iteration
update.cum <- function(var, M) {
    # Update days no platnorm, relapse, or gvhd
    M[, paste0("daysno", var)] <- M[, paste0("daysno", var)] + (1 - M[, var])
    # Update days platnorm, relapse, or gvhd
    M[, paste0("days", var)] <- M[, paste0("days", var)] + M[, var]
    if(var == "gvhd") {
        # Interaction terms
        M[, "day_gvhd"] <- M[, "gvhd"] * M[, var]
        M[, "daysq_gvhd"] <- M[, "gvhd"] * M[, var]
        M[, "daycu_gvhd"] <- M[, "gvhd"] * M[, var]
    }
    return(M)
}




simulate <- function(M, intervene, coef.dt) {
    total <- nrow(M)
    out.dt <- data.table()
    # Remove censoring for all interventions
    if(intervene > 0) {
        var.list <- setdiff(coef.dt$variable, c("censlost"))
    } else {
        var.list <- coef.dt$variable 
    }
    # Remove gvhd for always, never and acute
    if(intervene %in% 1:3) {
        var.list <- setdiff(var.list, "gvhd")
    }
    # Loop through each day
    for(i in 1:max(dt$day)) {
        # Add back in gvhd after 100 days for acute
        if(intervene == 3 & i > 100) {
            var.list <- c(var.list, "gvhd")
        }
        # remove gvhd after 100 days for chronic
        if(intervene == 4 & i > 100) {
            var.list <- setdiff(var.list, "gvhd")
        }
        # Update lags and time variables
        M <- update.lags(lag.list, M)
        M <- update.time(M)
        # Generate a predicted value for each variable in order
        for(var in var.list) {
            var.idx <- which(colnames(M) == var)
            M <- gen.draws(var, M, coef.dt)
        }
        # Pull out dead or censored and add to out table
        dead.dt <- as.data.table(M[M[, "d"] == 1 | M[, "censlost"] == 1, ,drop = F])
        if(nrow(dead.dt) > 0) {
            out.dt <- rbind(out.dt, dead.dt)
            M <- M[!(M[, "d"] == 1 | M[, "censlost"] == 1), ,drop = F]
        }
        if(nrow(M) == 0) {
            break
        }
        print(paste0("Day ", i, " of ", max(dt$day), "; P(alive) = ", round(nrow(M) / total, 2)))
        # Update cumulative variables
        for (var in lag.list) {
            M <- update.cum(var, M)
        }
    }
    out.dt <- rbind(out.dt, as.data.table(M))
    return(out.dt)
}

# Natural course
n.draws <- 1000
M.nat <- gen.baseline(n.draws, dt[day == 1], intervene = 0)
sim.nat <- simulate(M.nat, intervene = 0, coef.dt)

# Intervention (No GvHD)
M.ngvhd <- gen.baseline(n.draws, dt[day == 1], intervene = 1)
sim.ngvhd <- simulate(M.ngvhd, intervene = 1, coef.dt)

# Intervention (Always GvHD)
M.agvhd<- gen.baseline(n.draws, dt[day == 1], intervene = 2)
sim.agvhd <- simulate(M.agvhd, intervene = 2, coef.dt)

# Prevent Acute
M.acute<- gen.baseline(n.draws, dt[day == 1], intervene = 3)
sim.acute <- simulate(M.acute, intervene = 3, coef.dt)

# Prevent Chronic
M.chronic <- gen.baseline(n.draws, dt[day == 1], intervene = 4)
sim.chronic <- simulate(M.chronic, intervene = 4, coef.dt)


## Step 6 - concatentate intervetion data sets and run Cox model
sim.nat[, scenario := "Natural Course"]
sim.agvhd[, scenario := "Always GvHD"]
sim.ngvhd[, scenario := "No GvHD"]
sim.acute[, scenario := "No Acute GvHD"]
sim.chronic[, scenario := "No Chronic GvHD"]
setnames(data.dt, c("t", "d_dea"), c("day", "d"))
data.dt[, scenario := "Data"]
sim.dt <- rbindlist(list(sim.nat, sim.agvhd, sim.ngvhd, sim.acute, sim.chronic), use.names = T, fill = T)
out.dt <- sim.dt[, .(day, d, scenario)]
write.csv(out.dt, "results/sim_results.csv", row.names = F)

# Make kaplan meier survival plots
surv.fit <- survfit(Surv(day, d) ~ scenario, data=sim.dt, conf.type = "log-log")
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

# Never GvHD
print.hazard(sim.nat, sim.ngvhd)

# No Acute
print.hazard(sim.nat, sim.acute)

# No Chronic
print.hazard(sim.nat, sim.chronic)
© 2019 GitHub, Inc.
Terms
Privacy
Security
Status
Help
Contact GitHub
Pricing
API
Training
Blog
About
