# -----------------------
# LOAD NECESSARY LIBRARIES
# -----------------------
library(tidyverse)
library(dplyr)
library(parallel)

# -----------------------
# GLOBAL PARAMETER INITIALIZATION
# -----------------------

# Define simulation parameters globally
investigate_date <- c(seq(50, 150, by = 10), seq(400, 600, by = 10))
sim_days <- 601
n <- 50000

inc.pd <- 5
rec.pd <- 7
num_int <- 25
prop_rec <- 0.02  # Proportion of recovered individuals
prop_symp <- 1/3
prop_high_risk <- 0
high_risk_multiplier <- 2 #high risk individuals twice likely to get infected
R <- 2
days_reinf <- 30  # Number of days before chance of reinfection
reinf <- 0.5  # Relative probability of reinfection
prop_vax <- 0.2  # Proportion of individuals to vaccinate
days_wane <- 0  # Number of days before immunity starts to wane
wane_rate <- 0.01
beta <- R / rec.pd
vax_days <- seq(1, 600, 15)

infectious_num <- list()
cumulative_case <- list()
new_case_daily <- list()
cumu_case <- 0

# -----------------------
# SOURCE HELPER FUNCTIONS
# -----------------------
source("helper_functions.R")

# -----------------------
# MANUAL CONFIGURATION
# -----------------------
subtask_num <- 1

simulations_per_subtask <- 25

start_sim <- (subtask_num - 1) * simulations_per_subtask + 1
end_sim <- subtask_num * simulations_per_subtask

output_dir <- paste0("squared_subtask_", subtask_num)
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# -----------------------
# DEFINE SIMULATION FUNCTION
# -----------------------

run_simulation <- function(sim_num, subtask_num) {
  
  
  # Susceptible Nodes
  S_nodes <- tibble(
    ID = 1:(n - num_int * 2 - prop_rec * n),
    vax_day = NA_integer_,
    VEs = 0,
    VEp = 0,
    prior_inf = 0,
    day_prior_inf = NA_integer_,
    day_exposed = NA_integer_,
    symp = NA_integer_,
    day_inf = NA_integer_,
    day_rec = NA_integer_,
    history = "1",
    immunity = 0,
    immunity_start = NA_integer_,
    immunity_model = 0,
    dpv = NA_integer_,
    status = "susceptible"
  ) %>%
    mutate(high_risk = rbinom(n(), 1, prop_high_risk))
  
  
  # Exposed Nodes
  E_nodes <- tibble(
    ID = ((n - num_int * 2 - prop_rec * n) + 1):((n - num_int * 2 - prop_rec * n) + num_int),
    vax_day = NA_integer_,
    VEs = 0,
    VEp = 0,
    prior_inf = 0,
    day_prior_inf = NA_integer_,
    day_exposed = sample(-1:1, num_int, replace = TRUE),
    history = "2",
    immunity_model = 0,
    dpv = NA_integer_,
    status = "exposed"
  ) %>%
    mutate(
      day_inf = day_exposed + inc.pd,
      day_rec = day_inf + rec.pd,
      immunity_start = day_exposed,
      immunity = 0.75,
      symp = NA_integer_,
      high_risk = rbinom(n(), 1, prop_high_risk)
    )
  
  # Infectious Nodes
  I_nodes <- tibble(
    ID = ((n - num_int * 2 - prop_rec * n) + num_int + 1):((n - num_int * 2 - prop_rec * n) + num_int * 2),
    vax_day = NA_integer_,
    VEs = 0,
    VEp = 0,
    prior_inf = 0,
    day_prior_inf = NA_integer_,
    day_exposed = 1 - inc.pd,
    history = "2",
    immunity_start = 1 - inc.pd,
    immunity = 0.75,
    immunity_model = 0,
    dpv = NA_integer_,
    status = "infectious"
  ) %>%
    mutate(
      day_inf = 1,
      day_rec = 1 + rec.pd,
      symp = rbinom(n(), 1, prop_symp),
      high_risk = rbinom(n(), 1, prop_high_risk)
    )
  
  # Recovered Nodes
  R_nodes <- tibble(
    ID = ((n - num_int * 2 - prop_rec * n) + num_int * 2 + 1):n,
    vax_day = NA_integer_,
    VEs = 0,
    VEp = 0,
    prior_inf = 1,
    day_prior_inf = sample(seq(-100, -(inc.pd + rec.pd), 1), prop_rec * n, replace = TRUE),
    history = "23",
    immunity_model = 0,
    dpv = NA_integer_,
    status = "recovered"
  ) %>%
    mutate(
      day_exposed = day_prior_inf,
      day_inf = day_prior_inf + inc.pd,
      day_rec = day_inf + rec.pd,
      immunity_start = day_prior_inf,
      immunity = pmax(0.75 - (1 - immunity_start) * wane_rate, 0),
      symp = NA_integer_,
      high_risk = rbinom(n(), 1, prop_high_risk)
    )
  
  # Initialize cumulative cases
  cumu_case <- 0
  cumulative_case <- list()
  new_case_daily <- list()
  infectious_num <- list()
  
  # Simulation loop for each day
  for (t in 1:sim_days) {
    cat(sprintf("Subtask %d, Simulation %d, Day %d\n", subtask_num, sim_num, t))
    
    epidemic_curve <- NULL
    case_daily <- 0
    
    # Vaccinate if current day is one of the "vax days"
    if (t %in% vax_days) {
      S_nodes <- vaccinate(S_nodes, t)
      R_nodes <- vaccinate(R_nodes, t)
      E_nodes <- vaccinate(E_nodes, t)
      
      # Vaccinate infectious nodes with symp == 0
      I_nodes_symp_0 <- I_nodes[I_nodes$symp == 0, ]
      I_nodes_symp_1 <- I_nodes[I_nodes$symp == 1, ]
      I_nodes_symp_0 <- vaccinate(I_nodes_symp_0, t)
      I_nodes <- bind_rows(I_nodes_symp_0, I_nodes_symp_1)
    }
    
    # Update immunity level with waning
    # Only the immunity of susceptible and recovered agents can wane
    S_nodes <- wane(S_nodes, t)
    R_nodes <- wane(R_nodes, t)
    
    # I -> R 
    recover_I <- recover(I_nodes, t)
    I_nodes <- recover_I[[1]]
    R_nodes <- bind_rows(R_nodes, recover_I[[2]])
    
    # E -> I
    new_infectious <- E_to_I(E_nodes, t)
    E_nodes <- new_infectious[[1]]
    I_nodes <- bind_rows(I_nodes, new_infectious[[2]])
    
    # Infect S
    if (nrow(S_nodes) > 0) {
      S_nodes2 <- S_nodes
      for (node in 1:nrow(S_nodes)) {
        inf <- rbinom(1, 1, 
                      min(1, max(0, if (S_nodes$high_risk[node] == 1) {
                        high_risk_multiplier * beta * (1 - (S_nodes$immunity_model[node])^2) * (nrow(I_nodes) / n)
                      } else {
                        beta * (1 - (S_nodes$immunity_model[node])^2) * (nrow(I_nodes) / n)
                      }))
        )
        if (inf == 1) {
          S_nodes2 <- subset(S_nodes2, ID != S_nodes$ID[node])
          new_exposed_susceptible <- S_nodes %>%
            filter(ID == S_nodes$ID[node]) %>%
            mutate(
              day_exposed = t,
              immunity_start = t,
              day_inf = t + inc.pd,
              day_rec = t + inc.pd + rec.pd,
              history = paste0(history, "2"),
              status = "exposed",
              immunity = case_when(
                grepl("2", history) & grepl("4", history) ~ 1,  # Infected given vaccinated
                grepl("2", history) & !grepl("4", history) ~ 0.75,
                TRUE ~ immunity
              )
            )
          case_daily <- case_daily + nrow(new_exposed_susceptible)
          E_nodes <- bind_rows(new_exposed_susceptible, E_nodes)
        }
      }
      S_nodes <- S_nodes2
    }
    
    # Infect R
    if (nrow(R_nodes) > 0) {
      R_nodes2 <- R_nodes
      for (node in 1:nrow(R_nodes)) {
        inf <- rbinom(1, 1, 
                      min(1, max(0, if (R_nodes$high_risk[node] == 1) {
                        high_risk_multiplier * beta * (1 - (R_nodes$immunity_model[node])^2) * (nrow(I_nodes) / n)
                      } else {
                        beta * (1 - (R_nodes$immunity_model[node])^2) * (nrow(I_nodes) / n)
                      }))
        )
        if (inf == 1) {
          R_nodes2 <- subset(R_nodes2, ID != R_nodes$ID[node])
          new_exposed_recover <- R_nodes %>%
            filter(ID == R_nodes$ID[node]) %>%
            mutate(
              day_exposed = t,
              immunity_start = t,
              day_inf = t + inc.pd,
              day_rec = t + inc.pd + rec.pd,
              history = paste0(history, "2"),
              status = "exposed",
              immunity = 1
            )
          case_daily <- case_daily + nrow(new_exposed_recover)
          E_nodes <- bind_rows(new_exposed_recover, E_nodes)
        }
      }
      R_nodes <- R_nodes2
    }
    

    cumu_case <- cumu_case + case_daily
    cumulative_case <- append(cumulative_case, list(cumu_case))
    new_case_daily <- append(new_case_daily, list(case_daily))
    

    infected_num <- nrow(E_nodes) + nrow(I_nodes)
    cat(sprintf("Number of infectious individuals: %d\n", nrow(I_nodes)))
    

    if (t %in% investigate_date) {
      epidemic_curve <- bind_rows(epidemic_curve, E_nodes, I_nodes)
      specific_day_results <- epidemic_curve %>%
        mutate(eventstatus = 1) %>%
        bind_rows(S_nodes %>% mutate(eventstatus = 0)) %>%
        bind_rows(R_nodes %>% mutate(eventstatus = 0)) %>%
        mutate(
          TrialStatus = ifelse(is.na(vax_day), 0, 1),
          RI = ifelse(prior_inf == 1 & eventstatus == 1, 1, 0)
        ) %>%
        rename(PI = prior_inf, Symptomatic = symp) %>%
        dpv(t) %>%
        mutate(
          prop_rec = prop_rec, 
          prop_symp = prop_symp,
          R = R,
          days_reinf = days_reinf,
          reinf = reinf,
          sim_days = sim_days,
          VEs = VEs,
          VEp = VEp,
          prop_vax = prop_vax,
          days_wane = days_wane,
          wane_rate = wane_rate,
          day = t,
          subtask = subtask_num,
          simulation = sim_num,
          high_risk = high_risk
        )
      infectious_num <- append(infectious_num, list(nrow(I_nodes) + nrow(E_nodes)))
      

      filename <- sprintf("squared_simulation_subtask_%d_sim_%d_day_%d.csv", subtask_num, sim_num, t)
      write.csv(specific_day_results, file.path(output_dir, filename), row.names = FALSE)
    }
  }
}

# -----------------------
# MAIN SIMULATION LOOP WITH PARALLELISM
# -----------------------

num_cores <- 25


cl <- makeCluster(num_cores)

clusterExport(cl, varlist = c("run_simulation", "subtask_num", "sim_days", "n",
                              "inc.pd", "rec.pd", "num_int", "prop_rec",
                              "prop_symp", "R", "days_reinf", "reinf",
                              "prop_vax", "days_wane", "wane_rate", "beta",
                              "vax_days", "investigate_date", "output_dir", "prop_high_risk", "high_risk_multiplier"))



clusterEvalQ(cl, {
  library(tidyverse)
  library(dplyr)
  source("helper_functions.R")
})


simulations <- start_sim:end_sim


results <- parLapply(cl, simulations, function(sim_num) {
  run_simulation(sim_num, subtask_num)
})


stopCluster(cl)

cat(sprintf("Subtask %d: All %d simulations completed successfully.\n", subtask_num, simulations_per_subtask))
