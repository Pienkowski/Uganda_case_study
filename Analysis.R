############################
######### Analysis ######### 
############################

######### Steps ######### 
# 1) Set up environment 
# 2) Descriptive statistics
# 3) Simple model without spatial component 
# 4) Model with spatial component 
# 5) Model with imputed datasets 
# 6) Examine the results 
# 7) Look at predicted prevalence 
# 8) Extrapolate across Uganda 


######### 1) Set up environment ######### 
### Libraries ###
library(readxl)
library(plyr)
library(tidyverse)
library(ggplot2)
library(psych)
library(HH)
library(naniar)
library(mice)
library(logisticPCA)
library(lavaan)
library(semTools)
library(blavaan)
library(mirt)
library(semPlot)
library(DataExplorer)
library(sf)
library(corrplot)
library(gtsummary)
library(INLA)
library(gstat)
library(sjmisc)

### Load data ###
DF_stat.list <- readRDS("HH_locations/DF_stat.list3.rds")

### Project as equal earth ### 
DF_stat.list <- lapply(DF_stat.list, function(x) st_transform(x,crs = "+proj=eqearth"))

### Load functions ###
source("functions.R")

### Seed ###
set.seed(100)

######### 2) Descriptive statistics #########
### Number of respondents   
nrow(DF_stat.list[[1]])

### Proportion who pass the diagnostic threshold 
round(prop.table(table(DF_stat.list[[1]]$phq_01))*100, 2)

### Patterns of missing data ###
vis_miss(DF_stat.list[[1]]) 

### Spearman's rank correlation between the explanatory variables. Correlations greater than 0.7 are quite strong. ### 
# Convert ordinal to numeric 
DF_stat_temp <- st_drop_geometry(DF_stat.list[[1]])

# ord_to_num function 
ord_to_num <- function(data, ordered_columns) {
  for (col in ordered_columns) {
    if (is.factor(data[[col]])) {
      data[[col]] <- as.numeric(data[[col]])
    } else {
      warning(paste("Column '", col, "' is not an ordered factor. Skipping.", sep = ""))
    }
  }
  return(data)
}

# Run the function 
ord_f <- c("Health", "EDUCAT", "age_cat")
DF_stat_temp <- ord_to_num(data = DF_stat_temp,  ordered_columns = ord_f)

# Correlation plot (looking for correlations over .7)
cont_vars <- c("FIES_est", "Economic_poverty","Health","age_cat", "EDUCAT","Soci_est.1",  "Alcohol") 
str(DF_stat_temp[cont_vars])
M <- cor(DF_stat_temp[cont_vars], method = "spearman")
corrplot(M, method="number")

### Associations with the response variable ###
DF_stat_temp2 <- DF_stat.list[[1]]

# Setting JAMA theme for gtsummary
theme_gtsummary_journal("jama")
theme_gtsummary_compact(set_theme = TRUE, font_size = 10)

# Create a table
table_DF <- DF_stat_temp2 %>% select(phq_01, FIES_est,Economic_poverty, Health, age_cat,SEX,EDUCAT, Soci_est.1, MSTATUS, Alcohol,SMOKING) %>% as.data.frame()
table_DF  <- subset(table_DF , select = -geometry)

# Create table 
tab2SDME <-  tbl_summary(table_DF, by = phq_01, missing_text='Missing',
                         label = list(FIES_est  ~ "Food insecurity (latent)", 
                                      Economic_poverty ~ "Economic poverty (latent)",
                                      age_cat ~ "Age", 
                                      SEX ~ "Gender",
                                      EDUCAT ~ "Education",
                                      Soci_est.1 ~ "Social support (latent)",
                                      MSTATUS ~ "Marital status",
                                      Alcohol ~ "Alchohol consumption",
                                      SMOKING ~ "Smoking"), 
                         
                         statistic = list(all_continuous() ~ "{mean} ({sd})")) %>%
  italicize_levels() %>% 
  add_p(test.args = all_tests("fisher.test") ~ list(simulate.p.value=TRUE)) 

tab2SDME

######### 3) Simple model without spatial component ######### 
# Model test dataset
model_test <- DF_stat.list[[1]]

### The model without random effect ### 
start_time <- Sys.time()
Global_model_1 <- inla(phq_01 ~ FIES_est + Economic_poverty+ Health+ age_cat+ SEX+ EDUCAT+ Soci_est.1+ MSTATUS+ Alcohol+ SMOKING, 
                       control.compute = list(dic = TRUE, waic = TRUE),          
                       family = "binomial", 
                       data = model_test)

end_time <- Sys.time()
end_time - start_time

# Fixed effects 
model_1_res <-
  data.frame(
    Variable = row.names(Global_model_1$summary.fixed),
    Odds = exp(Global_model_1$summary.fixed$mean),
    Upper = exp(Global_model_1$summary.fixed$`0.025quant`),
    Lower = exp(Global_model_1$summary.fixed$`0.975quant`))
model_1_res


### The model with iid random effect ### 
start_time <- Sys.time()
Global_model_2 <- inla(phq_01 ~ FIES_est + Economic_poverty+ Health+ age_cat+ SEX+ EDUCAT+ Soci_est.1+ MSTATUS+ Alcohol+ SMOKING +
                         f(LOCB, model = "iid"), 
                       control.compute = list(dic = TRUE, waic = TRUE),          
                       family = "binomial", 
                       data = model_test)

end_time <- Sys.time()
end_time - start_time

# Fixed effects 
model_2_res <-
  data.frame(
    Variable = row.names(Global_model_2$summary.fixed),
    Odds = exp(Global_model_2$summary.fixed$mean),
    Upper = exp(Global_model_2$summary.fixed$`0.025quant`),
    Lower = exp(Global_model_2$summary.fixed$`0.975quant`))
model_2_res

# And compare the two models with DICs and WAICs
dic  <- c(Global_model_1$dic$dic, Global_model_2$dic$dic)
waic <- c(Global_model_1$waic$waic, Global_model_2$waic$wai)
Z    <- cbind(dic, waic)
rownames(Z) <- c("Bernoulli GLM",  
                 "Bernoulli GLM + iid")
Z


######### 4) Model with spatial component #########
### Identify coordinates ###
coords <- as.matrix(st_coordinates(model_test))

# max.edge guess (based on X meters)
max_edge_g <- 3000 / 5

# Create test meshes 
MeshA <- inla.mesh.2d(coords, max.edge = c(0.2, 0.5)*max_edge_g); plot(MeshA)
MeshB <- inla.mesh.2d(coords, max.edge = c(0.2, 0.5)*max_edge_g, cutoff = max_edge_g/6); plot(MeshB)
MeshC <- inla.mesh.2d(coords, max.edge = c(0.5, 1)*max_edge_g); plot(MeshC)
MeshD <- inla.mesh.2d(coords, max.edge = c(0.5, 1)*max_edge_g, cutoff = max_edge_g/4); plot(MeshD)
MeshE <- inla.mesh.2d(coords, max.edge = c(1, 1)*max_edge_g); plot(MeshE)
MeshF <- inla.mesh.2d(coords, max.edge = c(1, 1)*max_edge_g, cutoff = max_edge_g/2); plot(MeshF)

# A matrix maps the Gaussian Markov Random Field (GMRF) 
A.est <- inla.spde.make.A(mesh=MeshB,loc=as.matrix(coords)); dim(A.est)

# Create the spatial structure (SPDE object)
spde <- inla.spde2.matern(MeshB, alpha=2)

# Required indexes for the SPDE
iset <- inla.spde.make.index(name = "spatial.field", spde$n.spde)

# Create a stack 
Uganda.stack.est <- inla.stack(data=list(phq_01=model_test$phq_01),
                               A=list(A.est, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                               effects=list(c(iset, list(Intercept=1)),
                                            list(FIES_est=model_test$FIES_est), 
                                            list(Economic_poverty=model_test$Economic_poverty), 
                                            list(Health_Bad=model_test$Health_Bad), # RL is very bad
                                            list(Health_Fair=model_test$Health_Fair),
                                            list(Health_Good=model_test$Health_Good),
                                            list(Health_V.good=model_test$Health_V.good),
                                            list(age_cat_Under.25=model_test$age_cat_Under.25),
                                            list(age_cat_25.to.39=model_test$age_cat_25.to.39),
                                            list(age_cat_Over.55=model_test$age_cat_Over.55),
                                            list(SEX_Female=model_test$SEX_Female), # RL is female
                                            list(Soci_est.1=model_test$Soci_est.1),
                                            list(MSTATUS_Div.wid=model_test$MSTATUS_Div.wid),
                                            list(MSTATUS_Single=model_test$MSTATUS_Single),
                                            list(Alcohol=model_test$Alcohol),
                                            list(SMOKING=model_test$SMOKING)),
                               tag="est")

# Define the formula 
formula <- phq_01 ~ -1 + Intercept + 
  FIES_est +
  Economic_poverty +
  Health_Bad + Health_Fair + Health_Good + Health_V.good + 
  age_cat_Under.25 + age_cat_25.to.39 + age_cat_Over.55 +
  SEX_Female + 
  Soci_est.1 +
  MSTATUS_Div.wid +
  MSTATUS_Single +
  Alcohol +
  SMOKING +
  f(spatial.field, model=spde) 

# Run the model 
Global_model_3 <- inla(formula,
                       data=inla.stack.data(Uganda.stack.est, spde=spde),
                       family="binomial", Ntrials=1,
                       control.predictor=list(A=inla.stack.A(Uganda.stack.est),
                                              compute=TRUE),
                       control.compute=list(dic=TRUE))

# Convert from log-odds to odds 
data.frame(
  Variable = row.names(Global_model_3$summary.fixed),
  Odds = exp(Global_model_3$summary.fixed$mean),
  Upper = exp(Global_model_3$summary.fixed$`0.025quant`),
  Lower = exp(Global_model_3$summary.fixed$`0.975quant`))


# And compare the two models with DICs and WAICs
dic  <- c(Global_model_1$dic$dic, Global_model_3$dic$dic)
waic <- c(Global_model_1$waic$waic, Global_model_3$waic$wai)
Z    <- cbind(dic, waic)
rownames(Z) <- c("Bernoulli GLM",  
                 "Bernoulli GLM + SRF")
Z

# Precision
prec_2 <- Global_model_3$summary.hyperpar
round(prec_2,3)

# Variance
var_2 <- 1/prec_2
round(var_2,3)

######### 5) Model with imputed datasets ######### 
### Function to run the model ###
model_4_func <- function(model_test) {

  # Create a stack 
  Uganda.stack.est <- inla.stack(data=list(phq_01=model_test$phq_01),
                                 A=list(A.est, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                                 effects=list(c(iset, list(Intercept=1)),
                                              list(FIES_est=model_test$FIES_est), 
                                              list(Economic_poverty=model_test$Economic_poverty), 
                                              list(Health_Bad=model_test$Health_Bad), # RL is very bad
                                              list(Health_Fair=model_test$Health_Fair),
                                              list(Health_Good=model_test$Health_Good),
                                              list(Health_V.good=model_test$Health_V.good),
                                              list(age_cat_Under.25=model_test$age_cat_Under.25),
                                              list(age_cat_25.to.39=model_test$age_cat_25.to.39),
                                              list(age_cat_Over.55=model_test$age_cat_Over.55),
                                              list(SEX_Female=model_test$SEX_Female), # RL is female
                                              list(Soci_est.1=model_test$Soci_est.1),
                                              list(MSTATUS_Div.wid=model_test$MSTATUS_Div.wid),
                                              list(MSTATUS_Single=model_test$MSTATUS_Single),
                                              list(Alcohol=model_test$Alcohol),
                                              list(SMOKING=model_test$SMOKING)),
                                 tag="est")

  
  # Run the model 
  Global_model_4 <- inla(formula,
                         data=inla.stack.data(Uganda.stack.est, spde=spde),
                         family="binomial", Ntrials=1,
                         control.predictor=list(A=inla.stack.A(Uganda.stack.est),
                                                compute=TRUE, link = 1),
                         control.compute=list(dic=TRUE))
  return(Global_model_4)
}

### Run the function on the imputed datasets ###
Global_model_4 <- lapply(DF_stat.list, model_4_func)

### Examine the variation in the results between each model ###
# Implement Bayesian model averaging following Gómez-Rubio (2021) in section "12.4.1 Sampling from the imputation model"
test <- Global_model_4
test[[1]] <- NULL # For some reason, 
Global_model_4_mer <- inla.merge(test, rep(1, length(test)))



######### 6) Examine the results  ######### 
### Function to calculate fixed effect results and 95% CI in odds ###
fixed_exp_fun <- function(model, variables) {
  
  fixed_mean <- list()
  fixed_CI <- list()
  
  for (i in seq_along(1:length(variables))){
    fixed_mean[[i]] <- inla.emarginal(exp, model$marginals.fixed[[i]])
    fixed_CI[[i]] <-  inla.qmarginal(c(0.025, 0.975), inla.tmarginal(exp, model$marginals.fixed[[i]]))
  }
  
  fixed_mean_DF <- data.frame(Variable = variables,
                              Odds = do.call("rbind", fixed_mean),
                              Lower = (do.call("rbind", fixed_CI))[,1],
                              Upper = (do.call("rbind", fixed_CI))[,2],
                              Var_t = "Result")
  return(fixed_mean_DF)
}

# Implement function to extract mean and 95% CI in odds
Global_model_1_av <- fixed_exp_fun(model = Global_model_4_mer, variables = c(Global_model_4_mer$names.fixed))
Global_model_1_av

### Create a small coefficient plot ###
# Create function
plot_prep_f <- function(model_2_res) {
  
  # Remove intercept
  model_2_res <- model_2_res %>%  filter(!Variable=='Intercept')

  # Create reference level DF
  Ref_level <- data.frame(
    Variable = c(
      "Health: Very bad",
      "Age: 40-55",
      "Male",
      "Married/polygamous",
      "Non-smoker"
    ),
    Odds = c(rep(1, 5)),
    Upper = c(rep(1, 5)),
    Lower = c(rep(1, 5)),
    Var_t = c(rep("Reference", 5))
  )
  
  # Combine
  model_2_res <- rbind(model_2_res, Ref_level)
  
  # Rename variables
  model_2_res$Variable <- factor(
    model_2_res$Variable,
    levels = c(
      "FIES_est",
      
      "Economic_poverty",
      
      "Health: Very bad",
      "Health_Bad",
      "Health_Fair",
      "Health_Good",
      "Health_V.good",
      
      "age_cat_Under.25",
      "age_cat_25.to.39",
      "Age: 40-55",
      "age_cat_Over.55",
      
      "Male",
      "SEX_Female",
      
      "Soci_est.1",
      
      "Married/polygamous",
      "MSTATUS_Div.wid",
      "MSTATUS_Single",
      
      "Alcohol",
      
      "Non-smoker", 
      "SMOKING"    
    )
  )
  
  model_2_res$Variable <- model_2_res$Variable %>%
    dplyr::recode(
      "FIES_est" = "Food insecurity",
      
      "Economic_poverty" = "Economic poverty",
      
      "Health_Bad" = "Bad",
      "Health_Fair" = "Fair",
      "Health_Good" = "Good",
      "Health_V.good" = "Very good",
      
      "age_cat_Under.25" = "Under 25",
      "age_cat_25.to.39" = "25-39",
      "age_cat_Over.55" = "Over 55",
      
      "SEX_Female" = "Female",
      
      "Soci_est.1" = "Social support (family)",
      
      "MSTATUS_Div.wid" =  "Divorced or widow/er",
      "MSTATUS_Single" = "Never married",
      
      "Alcohol" = "Alchohold consumption frequency",
      
      "SMOKING" = "Smoker"
    
    )
  
  ### Split labels over multiple lines
  model_2_res$Variable_split <-
    as.factor(sapply(
      strwrap(model_2_res$Variable, 30, simplify = FALSE),
      paste,
      collapse = "\n"
    ))
  model_2_res <- with(model_2_res, model_2_res[order(Variable), ])
  model_2_res$Variable_split <-
    factor(model_2_res$Variable_split,
           levels = unique(model_2_res$Variable_split[order(model_2_res$Variable)]))
  return(model_2_res)
}

Global_model_1_av <- plot_prep_f(Global_model_1_av)


### Main plotting function 
plot_do_f <- function(model_res) {

  # Plot the graph
  p.1 <-
    ggplot(model_res,
           aes(
             x = Variable_split,
             y = Odds,
             ymin = Lower,
             ymax = Upper
           )) +
    geom_linerange(fatten = 0.5, size = 0.5) +
    geom_point(aes(shape = Var_t), size = 1) +
    scale_shape_manual(values = c(1, 16)) +
    scale_colour_manual(values = c("white", "black")) +
    ylab("Odds")
  
  p.1 <-
    p.1 + geom_hline(yintercept = 1,
                     color = "#8b0000",
                     size = 0.2)
  p.1 <- p.1 + coord_flip()
  p.1 <-
    p.1 + scale_x_discrete(breaks = c(levels(model_res$Variable_split)),
                           limits = rev(
                             c(
                               levels(model_res$Variable_split)[1],
                               "skip",
                               levels(model_res$Variable_split)[2],
                               "skip",
                               levels(model_res$Variable_split)[3:7],
                               "skip",
                               levels(model_res$Variable_split)[8:11],
                               "skip",
                               levels(model_res$Variable_split)[12:13],
                               "skip",
                               levels(model_res$Variable_split)[14],
                               "skip",
                               levels(model_res$Variable_split)[15:17],
                               "skip",
                               levels(model_res$Variable_split)[18],
                               "skip",
                               levels(model_res$Variable_split)[19:20]
                             )
                           ))  +
    theme_minimal() +
    theme(
      axis.text.y.right = element_blank(),
      axis.ticks.y.right = element_blank(),
      axis.line.y.right = element_blank(),
      axis.title.y = element_blank()
    ) 
  p.1 <- p.1 + theme(text = element_text(size = 8))
  p.1 <- p.1 + theme(legend.position = "none")  +
    theme(plot.margin = unit(c(0, 0, 0, 0), 'lines'))

  # Return the plot
  return(p.1)
}

p_main <- plot_do_f(Global_model_1_av)

# Save
ggsave("Fig_CS.png", plot = p_main,
       width = 3, height = 3.5, dpi = 500,
       bg = "white")

ggsave("Fig_CS.eps", plot = p_main,
       width = 3, height = 3.5, dpi = 500,
       bg = "white")

### Results ###
# Food insecurity 
FI_O <- format(round(Global_model_1_av[1,2],2), nsmall=2); FI_O
FI_L <- format(round(Global_model_1_av[1,3],2), nsmall=2); FI_L
FI_U <- format(round(Global_model_1_av[1,4],2), nsmall=2); FI_U

# Economic poverty 
EP_O <- format(round(Global_model_1_av[2,2],2), nsmall=2); EP_O
EP_L <- format(round(Global_model_1_av[2,3],2), nsmall=2); EP_L
EP_U <- format(round(Global_model_1_av[2,4],2), nsmall=2); EP_U
 
######### 7) Look at predicted prevalence ######### 
###### Function to run the model with the predicted prevalence ###### 
model_4_SB <- function(model_test) {
  
  # Step 1. Make a DF for the data predictor data
  model_test_train <- model_test
  model_test_train$Fit <- "Fit"
  model_test_pred <- model_test
  model_test_pred$phq_01  <- NA
  model_test_pred$Fit <- "Predict"
  
  # Step 2: Combine the original and predictor DF
  model_test_comb <- rbind(model_test_train, model_test_pred)
  
  ### Identify coordinates ###
  coords <- as.matrix(st_coordinates(model_test_comb))
  
  # max.edge guess (based on X meters)
  max_edge_g <- 3000 / 5
  
  # Create test meshes 
  MeshB <- inla.mesh.2d(coords, max.edge = c(0.2, 0.5)*max_edge_g, cutoff = max_edge_g/6); plot(MeshB)

  # A matrix maps the Gaussian Markov Random Field (GMRF) 
  A.est <- inla.spde.make.A(mesh=MeshB,loc=as.matrix(coords)); dim(A.est)
  
  # Create the spatial structure (SPDE object)
  spde <- inla.spde2.matern(MeshB, alpha=2)
  
  # Required indexes for the SPDE
  iset <- inla.spde.make.index(name = "spatial.field", spde$n.spde)
  
  # Create a stack 
  Uganda.stack.est <- inla.stack(data=list(phq_01=model_test_comb$phq_01),
                                 A=list(A.est, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                                 effects=list(c(iset, list(Intercept=1)),
                                              list(FIES_est=model_test_comb$FIES_est), 
                                              list(Economic_poverty=model_test_comb$Economic_poverty), 
                                              list(Health_Bad=model_test_comb$Health_Bad), # RL is very bad
                                              list(Health_Fair=model_test_comb$Health_Fair),
                                              list(Health_Good=model_test_comb$Health_Good),
                                              list(Health_V.good=model_test_comb$Health_V.good),
                                              list(age_cat_Under.25=model_test_comb$age_cat_Under.25),
                                              list(age_cat_25.to.39=model_test_comb$age_cat_25.to.39),
                                              list(age_cat_Over.55=model_test_comb$age_cat_Over.55),
                                              list(SEX_Female=model_test_comb$SEX_Female), # RL is female
                                              list(Soci_est.1=model_test_comb$Soci_est.1),
                                              list(MSTATUS_Div.wid=model_test_comb$MSTATUS_Div.wid),
                                              list(MSTATUS_Single=model_test_comb$MSTATUS_Single),
                                              list(Alcohol=model_test_comb$Alcohol),
                                              list(SMOKING=model_test_comb$SMOKING)),
                                 tag="est")
  
  # Run the model 
  Global_model_4_pred <- inla(formula,
                         data=inla.stack.data(Uganda.stack.est, spde=spde),
                         family="binomial", Ntrials=1,
                         control.predictor=list(A=inla.stack.A(Uganda.stack.est),
                                                compute=TRUE, link = 1),
                         control.compute=list(dic=TRUE))
  
  # Extract the fitted values
  fit_val <-
    cbind(data.frame(Rows = rownames(Global_model_4_pred$summary.fitted.values)),
          data.frame(Global_model_4_pred$summary.fitted.values))
  
  # Subset to the fitted values for just the observations
  fit_val_sub <- fit_val %>%
    filter(str_detect(Rows, "fitted.APredictor"))
  
  # Extract the mean
  model_test_comb$Fit_val <- fit_val_sub[, "mean"]
  
  # The fitted values for the predicted data
  pred <- model_test_comb[which(model_test_comb$Fit ==  "Predict"), ]$Fit_val
  
  # Original data from the fitting data
  orig <- model_test_comb[which(model_test_comb$Fit ==  "Fit"), ]$phq_01
  
  # Calculate the accuracy - following https://remiller1450.github.io/s230f19/caret2.html
  pred_cat <- ifelse(pred > .5, 1, 0)
  overall_accuracy <- sum(orig == pred_cat) / length(orig)
  overall_accuracy * 100
  
  # Return the results
  results <- list(pred_cat, overall_accuracy)
  
  return(results)
}

### Run the function on the imputed datasets ###
Global_model_4_pred <- lapply(DF_stat.list, model_4_SB)


###### Function to run the model with the predicted prevalence under the access restriction scenario ###### 
model_test <- DF_stat.list[[1]]

# Function to append a factor based on quntile groups of a numeric variable  
qunt_fact <- function(DF, var) {
  # Then let's create five levels
  var_level <-  cut(st_drop_geometry(DF[,var][[1]]),
                    breaks = c(-Inf, 
                               quantile(st_drop_geometry(DF[,var][[1]]), probs = seq(1/3, 2/3, 1/3),  na.rm = T ),
                               Inf), 
                    include.lowest = F, right=FALSE)
  
  table(var_level)
  
  # Recode 
  levels(var_level) <- c("Low", "Middle", "High")
  
  # Re-level so "Middle" is the reference level
  relevel(var_level, ref = "Middle")
  
  # As DF
  var_level_DF <- data.frame(var_level)
  
  # Rename column
  names(var_level_DF) <- paste0(var, "_level") 
  
  # Combine 
  DF <- cbind(DF, var_level_DF)
  return(DF)
}


###### Function to run the model with the predicted prevalence ###### 
model_4_scen <- function(model_test) {
  
  # Normalize the FIES_est and Economic_poverty instruments 
  model_test$FIES_est <- (model_test$FIES_est-min(model_test$FIES_est))/(max(model_test$FIES_est)-min(model_test$FIES_est))
  model_test$Economic_poverty <- (model_test$Economic_poverty-min(model_test$Economic_poverty))/(max(model_test$Economic_poverty)-min(model_test$Economic_poverty))
  
  # Step 1. Make a DF for the data predictor data
  model_test_train <- model_test
  model_test_train$Fit <- "Fit"
  model_test_pred <- model_test
  model_test_pred$phq_01  <- NA
  
  # Create the three scenario datasets 
  model_test_BUS <- model_test_pred
  model_test_FAS <- model_test_pred
  model_test_AFS <- model_test_pred
  
  # Tag the three datasets 
  model_test_BUS$Fit <- "BAU"
  model_test_FAS$Fit <- "FAS"
  model_test_AFS$Fit <- "AFS"
  
  # Respondents that are highly forest dependent experience a 25% increase in food insecurity and economic poverty
  model_test_FAS$FIES_est <- ifelse(qunt_fact(model_test_FAS, "Forest_est.1")$Forest_est.1_level == "High", 
                                    (model_test_FAS$FIES_est + 0.25), 
                                    model_test_FAS$FIES_est)
  
  model_test_FAS$Economic_poverty <- ifelse(qunt_fact(model_test_FAS, "Forest_est.1")$Forest_est.1_level == "High", 
                                    (model_test_FAS$Economic_poverty + 0.25), 
                                    model_test_FAS$Economic_poverty)
  
  # Respondents that have the lowest level of land cover have a 40% reduction in poverty 
  model_test_AFS$Economic_poverty <- ifelse(qunt_fact(model_test_AFS, "Land_est")$Land_est_level == "Low", 
                                            (model_test_AFS$Economic_poverty - 0.30), 
                                            model_test_AFS$Economic_poverty)
  
  # Step 2: Combine the original and predictor DF
  model_test_comb <- rbind(model_test_train, model_test_BUS, model_test_FAS, model_test_AFS)
  
  ### Identify coordinates ###
  coords <- as.matrix(st_coordinates(model_test_comb))
  
  # max.edge guess (based on X meters)
  max_edge_g <- 3000 / 5
  
  # Create test meshes 
  MeshB <- inla.mesh.2d(coords, max.edge = c(0.2, 0.5)*max_edge_g, cutoff = max_edge_g/6); plot(MeshB)
  
  # A matrix maps the Gaussian Markov Random Field (GMRF) 
  A.est <- inla.spde.make.A(mesh=MeshB,loc=as.matrix(coords)); dim(A.est)
  
  # Create the spatial structure (SPDE object)
  spde <- inla.spde2.matern(MeshB, alpha=2)
  
  # Required indexes for the SPDE
  iset <- inla.spde.make.index(name = "spatial.field", spde$n.spde)
  
  # Create a stack 
  Uganda.stack.est <- inla.stack(data=list(phq_01=model_test_comb$phq_01),
                                 A=list(A.est, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                                 effects=list(c(iset, list(Intercept=1)),
                                              list(FIES_est=model_test_comb$FIES_est), 
                                              list(Economic_poverty=model_test_comb$Economic_poverty), 
                                              list(Health_Bad=model_test_comb$Health_Bad), # RL is very bad
                                              list(Health_Fair=model_test_comb$Health_Fair),
                                              list(Health_Good=model_test_comb$Health_Good),
                                              list(Health_V.good=model_test_comb$Health_V.good),
                                              list(age_cat_Under.25=model_test_comb$age_cat_Under.25), # RL is 40 to 55
                                              list(age_cat_25.to.39=model_test_comb$age_cat_25.to.39),
                                              list(age_cat_Over.55=model_test_comb$age_cat_Over.55),
                                              list(SEX_Female=model_test_comb$SEX_Female), # RL is male
                                              list(Soci_est.1=model_test_comb$Soci_est.1),
                                              list(MSTATUS_Div.wid=model_test_comb$MSTATUS_Div.wid),
                                              list(MSTATUS_Single=model_test_comb$MSTATUS_Single),
                                              list(Alcohol=model_test_comb$Alcohol),
                                              list(SMOKING=model_test_comb$SMOKING)),
                                 tag="est")
  
  # Run the model 
  Global_model_4_pred <- inla(formula,
                              data=inla.stack.data(Uganda.stack.est, spde=spde),
                              family="binomial", Ntrials=1,
                              control.predictor=list(A=inla.stack.A(Uganda.stack.est),
                                                     compute=TRUE, link = 1),
                              control.compute=list(dic=TRUE))
  
  # Extract the fitted values
  fit_val <-
    cbind(data.frame(Rows = rownames(Global_model_4_pred$summary.fitted.values)),
          data.frame(Global_model_4_pred$summary.fitted.values))
  
  # Subset to the fitted values for just the observations
  fit_val_sub <- fit_val %>%
    filter(str_detect(Rows, "fitted.APredictor"))
  
  # Extract the mean
  model_test_comb$Fit_val <- fit_val_sub[, "mean"]
  
  # Original data from the fitting data
  orig <- model_test_comb[which(model_test_comb$Fit ==  "Fit"), ]$phq_01
  
  # The fitted values for the predicted data 
  # BAU
  pred_BAU <- model_test_comb[which(model_test_comb$Fit ==  "BAU"), ]$Fit_val
  pred_FAS <- model_test_comb[which(model_test_comb$Fit ==  "FAS"), ]$Fit_val
  pred_AFS <- model_test_comb[which(model_test_comb$Fit ==  "AFS"), ]$Fit_val
  
  # Calculate the accuracy - following https://remiller1450.github.io/s230f19/caret2.html
  pred_BAU_cat <- ifelse(pred_BAU > .5, 1, 0)
  pred_FAS_cat <- ifelse(pred_FAS > .5, 1, 0)
  pred_AFS_cat <- ifelse(pred_AFS > .5, 1, 0)
  
  # Overall accuracy 
  overall_accuracy <- (sum(orig == pred_BAU_cat) / length(orig))* 100

  # Return the results
  results <- list(overall_accuracy, pred_BAU_cat, pred_FAS_cat, pred_AFS_cat)
  
  return(results)
}

### Run the function on the imputed datasets ###
Global_model_4_pred <- lapply(DF_stat.list, model_4_scen)

### Extract the change in predicted depression ###
com_list <- list()
for (i in seq_along(1:length(Global_model_4_pred))) {
  com_list[[i]] <- data.frame(Iteration = i, 
                              BAU_sum = sum(Global_model_4_pred[[i]][[2]]),
                              FAS_sum = sum(Global_model_4_pred[[i]][[3]]),
                              AFS_sum = sum(Global_model_4_pred[[i]][[4]]))
  # Change in instances
  com_list[[i]]$BAU_FAS <- com_list[[i]]$FAS_sum - com_list[[i]]$BAU_sum
  com_list[[i]]$BAU_AFS <- com_list[[i]]$AFS_sum - com_list[[i]]$BAU_sum 
  
  # Percentage change 
  com_list[[i]]$BAU_FAS_p <- (com_list[[i]]$BAU_FAS / com_list[[i]]$BAU_sum)*100
  com_list[[i]]$BAU_AFS_p <- (com_list[[i]]$BAU_AFS / com_list[[i]]$BAU_sum)*100
}

######### 7) Extrapolate around Budongo and across Uganda #########
### Download and re-name WDPA data from: https://www.protectedplanet.net/en

### Read in WDPA data ###
WDPA_Uganda_0 <- read_sf('Spatial_data/WDPA_Uganda/WDPA_Uganda_shp_0/WDPA_Uganda_0.shp')
WDPA_Uganda_1 <- read_sf('Spatial_data/WDPA_Uganda/WDPA_Uganda_shp_1/WDPA_Uganda_1.shp')
WDPA_Uganda_2 <- read_sf('Spatial_data/WDPA_Uganda/WDPA_Uganda_shp_2/WDPA_Uganda_2.shp')

# Equal earth #
WDPA_Uganda_0 <- st_transform(WDPA_Uganda_0, crs = "+proj=eqearth")
WDPA_Uganda_1 <- st_transform(WDPA_Uganda_1, crs = "+proj=eqearth")
WDPA_Uganda_2 <- st_transform(WDPA_Uganda_2, crs = "+proj=eqearth")

# Combine #
WDPA_Uganda <- rbind(WDPA_Uganda_0, WDPA_Uganda_1, WDPA_Uganda_2)

# Filter by CFR #
WDPA_Uganda_BuRW <- WDPA_Uganda[WDPA_Uganda[[c("NAME")]] %in% c("Budongo", "Rwensama"), ]
WDPA_Uganda_CFR <- WDPA_Uganda[WDPA_Uganda[[c("DESIG")]] == "Forest Reserve", ]

### Read in settlement data ###
### Download GRID 3 data from: https://data.grid3.org/datasets/GRID3::grid3-uga-settlement-extents-v1-1/explore
Settle_Uganda <- read_sf('Spatial_data/GRID3_Uganda/GRID3_Uganda_Settlement_Extents%2C_Version_01.01..shp')

# Equal earth #
Settle_Uganda <- st_transform(Settle_Uganda, crs = "+proj=eqearth")

# Filter by rural # 
Settle_Uganda_rur <- Settle_Uganda[Settle_Uganda[[c("type")]] %in% c("Hamlet", "Small Settlement Area"), ]

### Intersect within buffer distance ### 
# Create buffer encompassing the maximum distance from CFR
max_dis <- 3407.636 # (m) Taken from previous script 
WDPA_Uganda_BuRW_buf <- st_buffer(WDPA_Uganda_BuRW, max_dis)
WDPA_Uganda_buf <- st_buffer(WDPA_Uganda, max_dis)

# Find those intersecting with the buffer 
Settle_Budongo_sel <- st_intersection(Settle_Uganda_rur, WDPA_Uganda_BuRW_buf)
Settle_Uganda_sel <- st_intersection(Settle_Uganda_rur, WDPA_Uganda_buf)

### Make the calculation ### 
# Based on UN population estimates 
Total_pop <- 48582334
Child_pop <- 25104152
Adult_pop <- 48582334 - 25104152
Adult_per <- Adult_pop/Total_pop

### Number adult pop around Budongo and CFM 
Adul_Bud <- sum(Settle_Budongo_sel$population, na.rm=T) *Adult_per
Adul_CFR <- sum(Settle_Uganda_sel$population, na.rm=T) *Adult_per
round(Adul_Bud)

### Number ###
calc_list <- list()
for (i in seq_along(1:length(Global_model_4_pred))) {
  # Results
  calc_list[[i]] <- data.frame(Iteration = i, 
                               BAU_Bud = Adul_Bud*(com_list[[i]]$BAU_sum / nrow(DF_stat.list[[1]])),
                               FAS_Bud = Adul_Bud*(com_list[[i]]$FAS_sum / nrow(DF_stat.list[[1]])), 
                               AFS_Bud = Adul_Bud*(com_list[[i]]$AFS_sum / nrow(DF_stat.list[[1]])),
                               BAU_CFR = Adul_CFR*(com_list[[i]]$BAU_sum / nrow(DF_stat.list[[1]])),
                               FAS_CFR = Adul_CFR*(com_list[[i]]$FAS_sum / nrow(DF_stat.list[[1]])),
                               AFS_CFR = Adul_CFR*(com_list[[i]]$AFS_sum / nrow(DF_stat.list[[1]])))
}

# As a DF 
calc_df <- do.call(rbind, calc_list)

### Max difference between scenarios - Budongo 
# FAS scenario - Budongo 
min_FAS <- round(min(calc_df$FAS_Bud) - max(calc_df$BAU_Bud)); min_FAS # Minimum plausible distance 
max_FAS <- round(max(calc_df$FAS_Bud) - min(calc_df$BAU_Bud)); max_FAS # Maximum plausible distance 
mean_FAS <- round((mean(calc_df$FAS_Bud) - mean(calc_df$BAU_Bud)) / mean(calc_df$BAU_Bud)*100,1)

# FAS scenario - Budongo 
min_AFS <- round(min(calc_df$AFS_Bud) - max(calc_df$BAU_Bud)); min_AFS # Minimum plausible distance 
max_AFS <- round(max(calc_df$AFS_Bud) - min(calc_df$BAU_Bud)); max_AFS # Maximum plausible distance 
mean_AFS <- round((mean(calc_df$AFS_Bud) - mean(calc_df$BAU_Bud)) / mean(calc_df$BAU_Bud)*100,1)



