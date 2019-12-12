#' Functions created to run Bayesian time series model usinh INLA


#load packages used during this analysis

packages_list <- c("tidyverse", "INLA", "SpatialEpi", "lubridate", "ggpubr")
lapply(packages_list, require, character.only = TRUE)



make_df <- function(path){
   df <- read_csv(path)
   df <- df %>% 
    mutate(time = lubridate::mdy(time),
          ym = format(time, "%Y-%m")) %>% 
    gather(key = "city", value = "value", -time, -ym) %>% 
    group_by(ym, city) %>% 
    summarise(var_sum = mean(value, na.rm = TRUE))
}


#Split data according to city
split_city <- function(data, city){
  ts <- data_frame(ym = seq(from = ymd("2007/01/01"), to = ymd("2018/12/31"), by = "month"))
  ts <- ts %>% mutate(ym = format(ym, "%Y-%m"))
  data <- data %>% 
    filter(CITY == !! city) %>% 
    select(DT_NOTIFIC, CITY) %>% 
    mutate(ym = format(DT_NOTIFIC, "%Y-%m")) %>% 
    group_by(ym, CITY) %>% 
    summarise(total = n())
  ts <-  left_join(ts, data)
  ts <- ts %>% 
    replace_na(list(total = 0)) %>% 
    mutate(time_str = 1:144,
	   time_non_str = 1:144,
           data = ymd(paste0(ym,"-","01"))) %>%
    rename(city = CITY) %>%
    mutate(city = na.omit(unique(city))[1])
}


#create a function to run a inla model
run_inla <- function(data, clima, plot_ci = FALSE) {

  data <- data %>%
    mutate(E = expected(pop,total, 1))
    model <- total ~ max_temp + min_temp + rain + hum +
    
    #define hyperpriors for the structuted random effect (rw1)
    f(time_str, model = "rw1", hyper = list("prec" = list(prior = "loggamma", param = c(0.5, 0.0005)))) +

    #define hyperpriors for the unstructuted random effect (iid)
    f(time_non_str, model = "iid", hyper = list("prec" = list(prior = "loggamma", param = c(0.5, 0.0005))))


  #run model
  model <- inla(model, family = "poisson", E = E, data = data,
                           
                           # Fixed effect priors: beta0, beta1,...,beta3
                           control.fixed = control.fixed(mean.intercept = 0, prec.intercept = 0.0001, 
                                                         mean = 0, prec = 0.0001),
                           control.predictor = list(compute = TRUE, link = 1),
                           control.compute = list(dic = TRUE, cpo = TRUE))
  
  #extract summary statistics
  table <- (exp(model$summary.fixed))
  
  point_estimate <- lapply(model$marginals.random$time_str,function(X){
    marg <- inla.tmarginal(function(x) exp(x), X)
    inla.emarginal(mean, marg)
  })
  
  lower_bond <- lapply(model$marginals.random$time_str,function(X){
    marg <- inla.tmarginal(function(x) exp(x), X)
    inla.hpdmarginal(0.95, marg)[1]
  })
  
   upper_bond <- lapply(model$marginals.random$time_str,function(X){
    marg <- inla.tmarginal(function(x) exp(x), X)
    inla.hpdmarginal(0.95, marg)[2]
  })
  
#Plot climate data

#Create a ggplot2 object to plot
  gg_df <- data.frame(Time = lubridate::ymd(data$data),
                    Estimated = unlist(point_estimate),
                    LB = unlist(lower_bond),
                    UB = unlist(upper_bond))
  
  graph <- ggplot(gg_df, aes(x = Time, y = Estimated)) + geom_line() + geom_hline(yintercept = 1) +
    theme_bw(base_size = 28) + ylab(expression(exp(omega[t])))
  city <- data$city[1]
  graph + ggsave(paste0("results/RR_",city,".tiff"), device = "tiff", width = 18, height = 12)
  if(isTRUE(plot_ci)){
    print(graph + geom_ribbon(aes(ymin = LB, ymax = UB), alpha = 0.5) + ylim(0, 3.5))
    graph <- graph + geom_ribbon(aes(ymin = LB, ymax = UB), alpha = 0.5) +
      ylim(0, 3.5)
    graph + ggsave(paste0("results/RR_",city,".tiff"), device = "tiff", width = 18, height = 12)
    return(list(graph, table, model))
  }
  else{
    print(graph)
    return(list(graph, table, model))
  }
  
  
  return(list(graph, table, model))
}

#compute prop
compute_prop <- function(data, variable){
    variable <- enquo(variable)
    data %>%
        group_by(!! variable) %>%
        summarise(numb = n()) %>%
        mutate(prop = numb/sum(numb))
}

#Compute desc
compute_desc <- function(data, variable, by, stats = "proportion"){
  variable <- enquo(variable)
  group_var <- enquo(by)

  if(stats == "median"){
    result <- data %>%
      group_by(!! group_var) %>%
        summarise(Median = median(!! variable, na.rm = TRUE))
  }
  if(stats == "proportion"){
    result <- data %>%
        group_by(!! group_var, !! variable) %>%
        summarise(numb = n()) %>%
        group_by(!! group_var) %>%
        mutate(Proportion = numb/sum(numb))
  }
  return(result)
}

#plot time series
plot_seas_ts <- function(data, variable, ylab){
    variable <- enquo(variable)
    data <- data %>%
      mutate(month = month(time),
            city = factor(city, levels = c("Porto Velho", "Ariquemes", "Cacoal", "Vilhena"))) 

    ggplot(data, aes(x = as.factor(month), y = !! variable)) +
      geom_boxplot() +
      theme_bw(base_size = 28) +
      facet_wrap(~ city) +
      labs(x = NULL, y = ylab)
}



predict_inla <- function(model, data){
  estimated <- list()
  
  #Combine values of the random effects
  time_component <- model$summary.random$time_str$mean + model$summary.random$time_non_str$mean
  
  #Retrive fixed effects
  intercept <- model$summary.fixed$mean[1]
  beta_max_temp <- model$summary.fixed$mean[2]
  beta_min_temp <- model$summary.fixed$mean[3]
  beta_rain <- model$summary.fixed$mean[4]
  beta_hum <- model$summary.fixed$mean[5]
  
  #Retrive data information
  max_temp <- data %>% pull("max_temp")
  min_temp <- data %>% pull("min_temp")
  rain <- data %>% pull("rain")
  hum <- data %>% pull("hum")
  pop <- data %>% pull("pop")
  total <- data %>% pull("total")
  
  #use a for loop to compute expected case in each month
  #this loop could be removed futher for optimization
  for(i in 1:nrow(data)){
    estimated[i] <- exp(intercept + beta_max_temp * max_temp[i] + beta_min_temp * min_temp[i] +
                          beta_rain * rain[i] + beta_hum * hum[i] + time_component[i]
                          ) * expected(pop[i],total[i],1)
  }

  
  return(unlist(estimated))
}






