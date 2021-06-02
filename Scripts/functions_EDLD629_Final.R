
## n = level 1; 
## j = level 2

simulate_two_level_interaction <-
  function(
    n, ## number level 1 units
    j, ## number level 2 units
    intercept_lv1, ## intercept at level 1
    main_x, ## main effect of x
    main_z, ## main effect of z
    interaction, ## interaction of x and z
    residual_var_sd_lv1, ## standard deviation of residuals at level 1
    random_int_mean_lv2, ## mean of random intercept at level 2
    random_int_sd_lv2, ## standard deviation of random intercept at level 2,
    start_seed = 123
  ){
    
    ## placeholder for school scores
    scores_j <- vector("list", j)
    
    ## for each variable, make list 
    ## as long as school
    ## fill each element with a list of student scores
    
    
    x_lv1 <- vector("list", j) 
    z_lv1 <- vector("list", j) 
    
    ## School distribution (intercept at school level)
    ## Standard deviation of random intercept at school level (?)
    set.seed(start_seed)
    
    a_j <- rnorm(j, random_int_mean_lv2, random_int_sd_lv2)
    
    for(i in 1:j) {
      
      ## make a level1 predictor variable:
      ## set multiple seeds that change each iteration
      ## prevents identical people but makes replicable
      
      set.seed(start_seed+i)
      x_lv1[[i]] <- rnorm(n)
      set.seed(start_seed-i)
      z_lv1[[i]] <- rnorm(n)
      set.seed(-start_seed + i)
      scores_j[[i]] <- 
        rnorm(
          n, 
          intercept_lv1 + 
            interaction*z_lv1[[i]]*x_lv1[[i]] + 
            main_x*x_lv1[[i]] + 
            main_z*z_lv1[[i]] + 
            a_j[i], 
          ## standard deviation of residual variance
          residual_var_sd_lv1
        )
    }
    
    scores_df <- 
      data.frame(
        scid = rep(1:j, each = n),
        x_lv1 = unlist(x_lv1),
        z_lv1 = unlist(z_lv1),
        score = unlist(scores_j)
      )
    
    return(scores_df)
    
  }


# pathpred <- 
#   function(object, ...){
#   ## coerce to "party" object if necessary
#   if(!inherits(object, "party")) object <- as.party(object)
# 
#   ## get standard predictions (response/prob) and collect in data frame
#   rval <- data.frame(response = predict(object, type = "response"))
#   rval$prob <- predict(object, type = "prob")
# 
#   ## get rules for each node
#   rls <- partykit:::.list.rules.party(object)
# 
#   ## get predicted node and select corresponding rule
#   rval$rule <- rls[as.character(predict(object, type = "node"))]
# 
#   return(rval)
# }
# 
# pathpred(tree1)

##### function to extract split points #####

extract_node_splits <- function(glmertree_mod, skip_test = FALSE) {
    
    if (skip_test == FALSE & length(glmertree_mod[["formula"]][[3]][[3]]) > 1)
      {
        message('More than 1 splitting variable... Check names of split variables')
      plot(tree1$tree, ask = F)
      }
    
    else {
      
    ni <- nodeids(glmertree_mod$tree)
    ni_terminal <- nodeids(glmertree_mod$tree, terminal = TRUE)
    ni_inner <- ni[!ni %in% ni_terminal]

    temp_vector <- 
      vapply(
        ni_inner, 
        function(i){
          split_node(node_party(glmertree_mod$tree[[i]]))$breaks
        },
        FUN.VALUE = numeric(1)
      )
    
    return(temp_vector)
    
    # tibble(split_points = temp_vector) %>%
    #   arrange(split_points) %>%
    #   mutate(number_of_splits = 1:length(temp_vector)) %>%
    #   select(number_of_splits, split_points)
    }
}

##### function to extract number of observations per node ###

get_node_counts <- function(glmertree_mod){
  glmertree_mod$data %>%
    transmute(node = .tree) %>% 
    count(node)
}

####### function to get number of nodes ######

get_number_of_nodes <- function(glmertree_mod){
  length(unique(glmertree_mod$data$.tree))
}

###### function for RMSE ######

rmse <- function(model){
  sqrt(sum(residuals(model)^2)/nobs(model))
}


#### function to extract variable ranges ####


##### function for standard error #####

se <- function(vector){
  if (sum(is.na(vector)) > 0) {
    vector <- na.omit(vector)
    message('Missing')
    sd(vector, na.rm = T)/sqrt(length(vector))
  }
  else{
  sd(vector, na.rm = T)/sqrt(length(vector))
  }
}




##### pool variance #####

# pool_variance <- function(df1, df2, sd1, sd2){
#   (df1*sd1^2 + df2*sd2^2)/(df1+df2)
# }

##### get standard errors from each regression coefficient #####
get_standard_errors <- 
  function(glmertree_mod){
    map_dbl(summary(glmertree_mod$tree),~coef(.x)[, "Std. Error"][2])
  }
##### get number of fixed effects from each regression coefficient #####

get_num_fixed_effects <- 
  function(glmertree_mod){
    map_dbl(summary(glmertree_mod$tree), ~nrow(coef(.x)))
  }


### function to get variances ####
get_variances <- 
  function(std_error, n){
    (std_error*sqrt(n))^2
  }

#### helper function for pooling se #####

# make_lag <- function(df){
#   df %>% 
#     mutate(
#       se_lag = lag(se), 
#       n_lag = lag(n)
#     ) %>% 
#     slice(-1)
# }

##### satterwaithe approximation of pooled SE for moderation effects ####

se_of_moderation <- function(glmertree_mod){
  
  temp_list  <- rep(NA, get_number_of_nodes(glmertree_mod)-1)
  
  for (i in 1:(get_number_of_nodes(glmertree_mod)-1)){
    temp_list [[i]] <- 
      satterwaithe_pooled_se(
        variance1 = get_variances(get_standard_errors(glmertree_mod)[[i]], get_node_counts(glmertree_mod)[,2][i]), 
        n1 = get_node_counts(glmertree_mod)[,2][i], 
        variance2 = get_variances(get_standard_errors(glmertree_mod)[[i+1]], get_node_counts(glmertree_mod)[,2][i+1]), 
        n2 = get_node_counts(glmertree_mod)[,2][i+1]
      )
  }
  temp_list 
}

##### function to make many data sets as one list #####



make_many_df <- function(
  intx_values, 
  num_df, 
  n_range, 
  j_range, 
  intercept_lv1, 
  main_x, 
  main_z,
  residual_var_sd_lv1,
  random_int_mean_lv2, 
  random_int_sd_lv2
  ){
  
  
  all_sims <- list(rep(NA, length(intx_values)))
  
  ns_param <- 
    new_quant_param(
      type = "integer",
      range = n_range,
      inclusive = c(TRUE, TRUE),
      trans = NULL,
      label = c(ns_param = "number of level 1 units"),
      finalize = NULL
    )
  
  
  js_param <- 
    new_quant_param(
      type = "integer",
      range = j_range,
      inclusive = c(TRUE, TRUE),
      trans = NULL,
      label = c(js_param = "number of level 2 units"),
      finalize = NULL
    )
  
  parameter_grid <- 
    grid_max_entropy(
      ns_param, 
      js_param, 
      size = num_df,
      original = FALSE)
  
  for (i in seq_along(intx_values)){
    
    intx_value <- intx_values[i]
    
    sim_dat <-
      pmap(
        list(
          parameter_grid$ns_param, 
          parameter_grid$js_param
          #parameter_grid$intx_param
        ), 
        ~simulate_two_level_interaction(
          n = ..1, 
          j = ..2, 
          intercept_lv1 = 10.00, 
          interaction = intx_value,
          main_x = 1.25, 
          main_z = 1.00,
          residual_var_sd_lv1 = 6.00,
          random_int_mean_lv2 = 5, 
          random_int_sd_lv2 = 2.00, 
          start_seed = 123
        )
      )
    
    all_sims[[i]] <- sim_dat
    
  }
  
  names(all_sims) <- intx_values
  
  return(all_sims)
}


##### function to get moderation estimate ####

### the se estimate here is not right. the project uses the standard error
### of the estimated average moderation effect. This would be useful
### to give the option if I publish this.

get_glmertree_moderation <- 
  function(
    glmertree_mod, 
    which_variable, 
    mean_only = T, 
    return_all_differences = F
    ){
    
    fixed_effects <- coef(glmertree_mod)[,which_variable]
    
    names(fixed_effects) <- NULL
    # 
    # test_tibble <- 
    #   data.frame(
    #     fixed_effect = fixed_effects, 
    #     lagged = lag(fixed_effects)
    #   )
    ### could use diff() instead of this method, too ###
    ## diff() gives differences between this and lagged automatically
    differences_in_slopes <- diff(fixed_effects)
    
    if (mean_only == T){return(mean(differences_in_slopes, na.rm = T))}
    
    output1 <- 
      data.frame(
        average_moderation = mean(differences_in_slopes, na.rm = T), 
        average_se_mod = 
          mean(
            se_of_moderation(glmertree_mod = glmertree_mod), na.rm = T
            )
      )
    
    output2<- 
      list(
        average_moderation = mean(differences_in_slopes, na.rm = T), 
        diff_in_slopes = differences_in_slopes, 
        se_moderation = se_of_moderation(glmertree_mod = glmertree_mod)
      )
    
    if (return_all_differences == F){return(output1)}
    
    else {return(output2)}
  }

##### function to plot the fit of models #####

plot_fits <- function(df){
  
  df %>% 
    group_by(model, metric) %>% 
    mutate(avg_perf = mean(performance, na.rm = T)) %>% 
    ungroup() %>% 
    ggplot(
      aes(x = performance)
    ) + 
    geom_density(
      alpha = 0.2, 
      aes(fill = model)
    ) + 
    geom_vline(aes(xintercept = avg_perf, color = model)) +
    facet_wrap(~metric, scales = 'free') +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 30), 
      legend.position = 'bottom')
  
}


##### function to plot the interaction estimate densities #####
## i.e., dfs for interaction estimates and fits & plots for both

plot_intx_ests <- function(df, intx_value){
  
  df %>% 
    group_by(type) %>% 
    mutate(average_effect = mean(estimate, na.rm = T)) %>% 
    ungroup() %>% 
    ggplot(aes(estimate, fill = type)) + 
    geom_density(alpha = 0.4) +
    geom_vline(
      aes(xintercept = average_effect, color = type),
      size = 2
    ) + 
    geom_vline(
      aes(xintercept = intx_value), 
      lty = 2, 
      size = 2
    ) + 
    theme_minimal() + 
    theme(legend.position = 'bottom')
}


##### large function to make all needed output #####

compare_mods <- function(sim_dat, intx_value){
  
  average_glmertree_moderation <- rep(NA, length(sim_dat))
  average_lmer_moderation <- rep(NA, length(sim_dat))
  lmer_t_value <- rep(NA, length(sim_dat))
  glmertree_aic <- rep(NA, length(sim_dat))
  lmer_aic <- rep(NA, length(sim_dat))
  glmertree_bic <- rep(NA, length(sim_dat))
  lmer_bic <- rep(NA, length(sim_dat))
  glmertree_rmse <- rep(NA, length(sim_dat))
  lmer_rmse <- rep(NA, length(sim_dat))
  
  for (i in 1:length(sim_dat)){
    
    
    temp_tree <- 
      lmertree(
        data = sim_dat[[i]], 
        formula = 
          score ~ 
          x_lv1 + z_lv1 |
          (1 | scid) | 
          z_lv1, 
        cluster = scid, 
        bonferroni = T
      )
    
    temp_lmer <- 
      lmer(
        data = sim_dat[[i]], 
        formula = 
          score ~ 
          x_lv1*z_lv1 + (1 | scid)
      )
    
    average_glmertree_moderation[i] <- 
      get_glmertree_moderation(
        temp_tree, 
        which_variable = 'x_lv1', 
        mean_only = T)
    
    average_lmer_moderation[i] <- 
      fixef(temp_lmer)[4]
    
    lmer_t_value[[i]] <- 
      summary(temp_lmer)$coefficients[[12]]
    
    glmertree_aic[i] <- AIC(temp_tree)
    lmer_aic[i] <- AIC(temp_lmer)
    glmertree_bic[i] <- BIC(temp_tree)
    lmer_bic[i] <- BIC(temp_lmer)
    glmertree_rmse[i] <- rmse(temp_tree)
    lmer_rmse[i] <- rmse(temp_lmer)
    
    print(
      paste0('iteration #', 
             i, 
             ' of ', 
             length(sim_dat), 
             ' complete for interaction = ', 
             intx_value,
             sep = '')
      )
  }
  
  fits <- 
    tibble(
      model = 
        c(
          rep('tree', length(sim_dat)), 
          rep('lmer', length(sim_dat)), 
          rep('tree', length(sim_dat)), 
          rep('lmer', length(sim_dat)), 
          rep('tree', length(sim_dat)), 
          rep('lmer', length(sim_dat))
        ), 
      metric = c(
        rep('aic', 2*length(sim_dat)), 
        rep('bic', 2*length(sim_dat)), 
        rep('rmse', 2*length(sim_dat))
      ),
      performance = 
        c(glmertree_aic, 
          lmer_aic, 
          glmertree_bic, 
          lmer_bic, 
          glmertree_rmse, 
          lmer_rmse)
    )
  
  fits_summary <- 
    fits %>% 
    group_by(model, metric) %>% 
    summarize(
      average_performance = mean(performance, na.rm = T), 
      se_performance = sd(performance, na.rm = T)/sqrt(n())
    ) %>% 
    ungroup() %>% 
    suppressMessages()
  
  intx_ests <- 
    tibble(
      index = rep(1:length(sim_dat), 2), 
      type = c(
        rep('lmer', length(sim_dat)),
        rep('tree', length(sim_dat))
      ), 
      estimate = 
        c(average_lmer_moderation, average_glmertree_moderation)
    )
  
  intx_summary <-
    intx_ests %>% 
    group_by(type) %>% 
    summarize(
      average_estimate = mean(estimate, na.rm = T), 
      se_estimate = sd(estimate, na.rm = T)/sqrt(n())
    ) %>% 
    ungroup() %>% 
    suppressMessages()
  
  output <- list(
    interaction_summary = intx_summary, 
    fits_summary = fits_summary,
    interactions = intx_ests,
    fits = fits,
    plot_fits = plot_fits(fits) + 
      labs(caption = paste0('interaction value = ', intx_value)), 
    plot_intx_ests = 
      plot_intx_ests(intx_ests, intx_value) + 
      labs(caption = paste0('interaction value = ', intx_value)))
  
  return(output)
}

##### make a list with all results you want #####

make_final_results <- 
  function(
    intx_values, 
    num_df, 
    n_range, 
    j_range, 
    intercept_lv1, 
    main_x, 
    main_z,
    residual_var_sd_lv1,
    random_int_mean_lv2, 
    random_int_sd_lv2
  ){
    map2(.x = 
           make_many_df(
             intx_values = intx_values, 
             num_df = num_df, 
             n_range = n_range, 
             j_range = j_range, 
             intercept_lv1 = intercept_lv1, 
             main_x = main_x, 
             main_z = main_z,
             residual_var_sd_lv1 = residual_var_sd_lv1,
             random_int_mean_lv2 = random_int_mean_lv2, 
             random_int_sd_lv2 = random_int_sd_lv2), 
         .y = intx_values, 
         ~compare_mods(.x, .y)
    )
  }

##### helper functions to extract pieces of the list made make_final_results():

extract_sim_ests <- 
function(full_list){map(full_list, ~.x[[1]])}

extract_sim_perfomance <- 
  function(full_list){map(full_list, ~.x[[2]])}

extract_sim_fit_plots <- 
  function(full_list){map(full_list, ~.x[[5]])}

extract_sim_ests_density <- 
  function(full_list){map(full_list, ~.x[[6]])}

###### function to create 95% credible (?) intervals #####
## this especially could be improved
c_intervals <- 
  function(
    full_list, 
    which_model, 
    intx_level,
    lower_level = 0.025, 
    upper_level = 0.975
  ){
    
    if (is.character(which_model) != T | is.character(intx_level) != T){
      warning('`which_model` & `intx_level` must be character')
    }
    vec <- 
      full_list[[intx_level]][['interactions']] %>% 
      filter(type == which_model) %>% 
      pull(estimate)
    
    quantile(vec, probs = c(lower_level, upper_level), na.rm = T)
    
  }
