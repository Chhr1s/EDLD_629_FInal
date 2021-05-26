
## n = level 1; 
## j = level 2

simulate_two_level <-
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