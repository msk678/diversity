## function to calculate multidimensional Shannon diversity as the MLE
mdiv_i <- function(df){
  
  # Shannon entropy based on observed joint proportions
  obs_prop <- c(prop.table(table(df))) # find proportions of all variable-level combinations
  obs_prop <- obs_prop[obs_prop != 0] # remove any zeros
  entropy_obs <- sum(-obs_prop*log(obs_prop), na.rm=T) # calculate the entropy
  diversity_obs <- exp(entropy_obs) # transform entropy to diversity
  
  # Shannon entropy based on joint proportion as products of 
  # predicted proportions of individual variables (i.e. assuming independence)
  nvar <- dim(df)[2] # number of variables in the data set
  
  
  # write observed variable combination as a new column into the data frame
  colNames = colnames(df) # variable names in the data set 
  df$comb = apply(df[, colNames, drop = F], MARGIN = 1, FUN = function(i) paste(i, collapse = ""))
  # write corresponding observed proportion into the data frame
  obs_prop2 <- prop.table(table(df$comb))
  L<-dim(df)[1] # number of observations in the data set
  for (i in 1:L){
    df$obs_prop[i] <- obs_prop2[names(obs_prop2)==df$comb[i]]
  }
  ## generate data table with
  # a row for each unique combination of variables that is observed in the data
  # a column for each variable containing the observed proportion of the level  
  # in the specified combination
  # a column specifying the observed proportion of the level combination
  # a column specifying the predicted proportion of the level combination
  # assuming independence
  
  # create a data table with variable levels replaced by variable level proportions
  df_prop <- data.frame(row.names = (1:L))
  for (i in (1:nvar)){
    p_ind <- prop.table(table(df[i]))
    for (j in 1:L){
      df_prop[j,i] <- p_ind[names(p_ind) == as.character(df[j,i])]
    }
  }
  # copy over combination identifier and observed combination proportion from df
  df_prop$comb <- df$comb
  df_prop$obs_prop <- df$obs_prop
  # keep only unique combinations
  df_prop_unique <- df_prop[which(duplicated(df_prop$comb)==F),]
  # calculate predicted joint probability assuming independence
  df_prop_unique$joint <- apply(df_prop_unique[1:nvar],1,prod)
  
  # calculate entropies for each variable level combinations assuming independence
  df_prop_unique$entropy_joint <- -df_prop_unique$joint*log(df_prop_unique$joint)
  
  # sum into overall predicted entropy of the data set
  entropy_pred <- sum(df_prop_unique$entropy_joint)
  diversity_pred <- exp(entropy_pred)
  
  # output the values ofboth predicted and observed entropy and diversity  
  output <- c(entropy_obs,diversity_obs,entropy_pred,diversity_pred)
  return(output)
}

################################################################################

## function to calculate multidimensional Simpson diversity as the MLE
mdiv_Simpson <- function(df){
  
  # Simpson entropy based on observed joint proportions
  obs_prop <- c(prop.table(table(df))) # find proportions of all variable-level combinations
  obs_prop <- obs_prop[obs_prop != 0] # remove any zeros
  simp_obs <- 1-sum(obs_prop**2, na.rm=T) # calculate the simpson index
  diversity_obs <- 1/(1-simp_obs) # transform entropy to diversity
  
  # Shannon entropy based on joint proportion as products of 
  # predicted proportions of individual variables (i.e. assuming independence)
  nvar <- dim(df)[2] # number of variables in the data set
  
  
  # write observed variable combination as a new column into the data frame
  colNames = colnames(df) # variable names in the data set 
  df$comb = apply(df[, colNames, drop = F], MARGIN = 1, FUN = function(i) paste(i, collapse = ""))
  # write corresponding observed proportion into the data frame
  obs_prop2 <- prop.table(table(df$comb))
  L<-dim(df)[1] # number of observations in the data set
  for (i in 1:L){
    df$obs_prop[i] <- obs_prop2[names(obs_prop2)==df$comb[i]]
  }
  ## generate data table with
  # a row for each unique combination of variables that is observed in the data
  # a column for each variable containing the observed proportion of the level  
  # in the specified combination
  # a column specifying the observed proportion of the level combination
  # a column specifying the predicted proportion of the level combination
  # assuming independence
  
  # create a data table with variable levels replaced by variable level proportions
  df_prop <- data.frame(row.names = (1:L))
  for (i in (1:nvar)){
    p_ind <- prop.table(table(df[i]))
    for (j in 1:L){
      df_prop[j,i] <- p_ind[names(p_ind) == as.character(df[j,i])]
    }
  }
  # copy over combination identifier and observed combination proportion from df
  df_prop$comb <- df$comb
  df_prop$obs_prop <- df$obs_prop
  # keep only unique combinations
  df_prop_unique <- df_prop[which(duplicated(df_prop$comb)==F),]
  # calculate predicted joint probability assuming independence
  df_prop_unique$joint <- apply(df_prop_unique[1:nvar],1,prod)
  
  # calculate entropies for each variable level combinations assuming independence
  df_prop_unique$simp_joint <- -df_prop_unique$joint**2
  
  # sum into overall predicted entropy of the data set
  simp_pred <- 1-sum(df_prop_unique$simp_joint)
  diversity_pred <- 1/(1-simp_pred)
  
  # output the values ofboth predicted and observed entropy and diversity  
  output <- c(simp_obs,diversity_obs,simp_pred,diversity_pred)
  return(output)
  
}

################################################################################


Chao_Shen_Entropy <- function(df){
  # require(tidyverse)
  # joint proportions observed:
  obs_prop <- c(prop.table(table(df)))
  obs_prop <- obs_prop[obs_prop != 0]
  
  # number of singleton observations
  obs_singletons <- c(table(df))
  obs_singletons <-obs_singletons[obs_singletons==1]
  f1 <- sum(obs_singletons)
  
  # number of total observations
  n <- dim(df)[1]
  
  # adjusted coverage
  C_hat <- 1-f1/n
  
  # adjusted estimated probabilities
  p_hat <- obs_prop*C_hat
  
  # Chao_Shen estimator for entropy
  CS_H <- -sum(p_hat*log(p_hat)/(1-(1-p_hat)^n), na.rm=T)
  CS_D <- exp(CS_H)
  
  # output results
  output <- c(CS_H, CS_D)
  return(output)
  
}



