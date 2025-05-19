recover <- function(df, t) {
  df_recovered <- df %>%
    filter(t == day_rec) %>%
    mutate(prior_inf = 1,
           day_prior_inf = day_inf,
           history = paste0(history, "3"),
           status = "recovered",
           immunity_model = immunity)
  
  df <- df %>%
    filter(t < day_rec)
  
  list(df, df_recovered)
}

E_to_I <- function(df, t) {
  df_inf <- df %>%
    filter(t == day_inf) %>%
    rowwise() %>%
    mutate(symp = rbinom(1, 1, prop_symp * (1 - VEp)),
           status = "infectious") %>%
    ungroup()
  
  df <- df %>%
    filter(t < day_inf)
  
  list(df, df_inf)
}

vaccinate <- function(df, t) {
  num_vax <- round(prop_vax * nrow(df %>% filter(is.na(vax_day))) / length(vax_days))
  
  new_vax <- df %>%
    filter((is.na(vax_day))|(t-vax_day >= 56)) %>% #those have not received vaccines or 8-weeks away from last vaccination date
    pull(ID) %>%
    sample(size = num_vax)
  
  df <- df %>%
    mutate(vax_day = case_when(ID %in% new_vax ~ as.integer(t),
                               TRUE ~ as.integer(vax_day)),
           VEs = case_when(vax_day == t ~ VEs,
                           TRUE ~ VEs),
           VEp = case_when(vax_day == t ~ VEp,
                           TRUE ~ VEp),
           history = case_when(ID %in% new_vax ~ paste0(history, "4"),
                               TRUE ~ history),
           immunity = case_when(grepl("4", history) & grepl("2", history) ~ 1,
                                grepl("4", history) & !grepl("2", history)  ~ 0.75,
                                TRUE ~ immunity),
           immunity_start = case_when(ID %in% new_vax ~ t,
                                      TRUE ~ immunity_start))
  
  return(df)
}

wane <- function(df, t) {
  df <- df %>%
    rowwise() %>%
    mutate(immunity = case_when(!is.na(immunity_start) & ((t - immunity_start) > days_wane) ~ pmax(0, immunity - wane_rate),
                                TRUE ~ immunity)) %>% 
    mutate(immunity_model = immunity)
  
  return(df)
}

dpv <- function(df, day) {
  df %>%
    mutate(dpv = ifelse(is.na(vax_day), 0, day - vax_day))
}
