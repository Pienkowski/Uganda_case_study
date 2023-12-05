##############################
######### Processing ######### 
##############################

######### Steps ######### 
# 1) Set up environment 
# 2) Explore and process the raw data
# 3) Examine patterns of missing data and perform the multiple imputation
# 4) Post imputation processing 

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
library(fastDummies)
library(DataExplorer)
library(sf)
library(sjmisc)

### Load data ###
DF_anly_1 <- st_as_sf(readRDS("HH_locations/DF_anly.rds"))

### Load functions ###
source("functions.R")

### Seed ###
set.seed(55)

### Split into spatial and non-spatial datasets ###
DF_anly <- st_drop_geometry(DF_anly_1)
DF_anly_Sp <- DF_anly_1[c("REF-PATID", "geometry", "LOCB")]


######### 2) Explore and process the raw data #########
### Identify training and test observations (for later) ### 
train = sort(sample(nrow(DF_anly), nrow(DF_anly)*2/3))

###### Survey data exploration ###### 
### Location ### 
table(DF_anly$LOCB)

### Gender ### 
table(DF_anly$SEX)

### Age ### 
hist(DF_anly$DOBB)

### Education ### 
table(DF_anly$EDUCAT)

### Marital status ### 
table(DF_anly$MSTATUS)

### Health ### 
table(DF_anly$Health)

### Number of children ### 
hist(DF_anly$HOUCHI)

### Number of adults ### 
hist(DF_anly$HOUADU)

### Main source of income ### 
table(DF_anly$ECCOCC)

### Asset ownership ### 
# Subset DF to assets 
assets <- c("BNSGAS","BNSCLOT","BNSGOAT","BNSMEAL","BNSELEC","BNSBED"  ,
            "BNSSOLA","BNSJERR","BNSMATT","BNSTELV","BNSBRICK","BNSMOTO",
            "BNSHOES","BNSBLAK","BNSBACC","BNSBIKE","BNSTANK","BNSSHOE",
            "BNSROOF","BNSRADI" ,"BNSSAUC","BNSBATT","BNSTELL","BNSCHAI",
            "BNSCAR","BNSFRID","BNSCONC","BNSSAV", "BNSSOF","BNSWOD",
            "BNSSOP")
assets_DF <- DF_anly[assets]
# Some of these are removed later 

# Columns to rows 
assets_DF_l<- gather(assets_DF)

# Yes = 1
assets_DF_l$value <- ifelse(assets_DF_l$value == "Yes", 1, 
                            ifelse(assets_DF_l$value == "No", 0,
                                   ifelse(is.na(assets_DF_l$value)== T, NA, 99)))
table(assets_DF_l$value)

# Ownership of each item 
ggplot(data=assets_DF_l, aes(x=key, y=value)) +
  geom_bar(stat="identity") + coord_flip()

### Subjective financial strain ### 
table(DF_anly$Strain)

### Subjective land size ### 
table(DF_anly$Land_sub)

### Estimated land size ###
hist(DF_anly$Land_size) # All data 
hist(DF_anly[which(DF_anly$Land_size <1),]$Land_size) # Under 1 hectare 

### Growing sugar cane? ### 
table(DF_anly$GROWSUG)

### Land size instrument ###
names_land <- c("LANDBIG", "LANDBAS", "LANDSTR", "LANDWELL", "LANDGOOD", "LANDSMALL")

# Subset to land items 
land_DF <- DF_anly[names_land]

# create a DF with the proportions for the first land variable 
land_prop <- data.frame(Variable = names_land[1], cbind(t(prop.table(table(land_DF[1])))))

# For loop appending the proportions for the remaining variables 
for (i in seq_along(2:(length(names_land)))){
  land_prop <- rbind(land_prop, data.frame(Variable = names_land[i+1], cbind(t(prop.table(table(land_DF[i+1]))))))
}

# Rename variable
land_prop$Variable <- fct_recode(land_prop$Variable, "LANDSTR R" = "LANDSTR", "LANDSMALL R" = "LANDSMALL")

# Likert scaled plot 
HH::likert(Variable~.,land_prop, positive.order=TRUE,as.percent = TRUE,
           main="Land instrument",
           xlab="Percentage", ylab="Variable")

# Convert from character to numeric - positively coded items 
land_pos <- c("LANDBIG", "LANDBAS", "LANDWELL", "LANDGOOD")
table(land_DF[land_pos][1])
land_DF[land_pos] <- apply(land_DF[land_pos], 2, agree.disagree.pos)
table(land_DF[land_pos][1])

# Convert from character to numeric - negatively coded items (i.e. reverse code)
land_neg <- c("LANDSTR",  "LANDSMALL")
table(land_DF[land_neg][1])
land_DF[land_neg] <- apply(land_DF[land_neg], 2, agree.disagree.neg)
table(land_DF[land_neg][1])

# Parallel analysis (with training data)
fa.parallel(land_DF[complete.cases(land_DF),][train,] , cor = "poly", fm="wls", fa="fa",   main = "Parallel analysis")

# Factor analysis, with polychoric correlation, WLS estimator, and oblimin rotation 
fa_mod1 <- efaUnrotate(land_DF[train,], estimator = "WLS", nf = 1)

# Factor loadings 
summary(fa_mod1)

# RMSEA
fitmeasures(fa_mod1)["rmsea"]

# Confirmatory analysis with test data 
CFA_1 <- 'factor.1 =~ LANDBIG  + LANDBAS + LANDSTR + LANDWELL + LANDGOOD + LANDSMALL'
CFA_mod1 <- cfa(CFA_1, data = land_DF[-train,])

# Factor loadings 
summary(CFA_mod1)

# RMSEA - this suggests poor model fit
fitmeasures(CFA_mod1)["rmsea"]

# Factor analysis, with polychoric correlation, WLS estimator, and oblimin rotation using all data 
fa_mod1b <- efaUnrotate(land_DF, estimator = "WLS", nf = 1)

# Factor loadings 
summary(fa_mod1b)

# RMSEA - much better than the above
fitmeasures(fa_mod1b)["rmsea"]

### Forest dependency instrument ###
names_fore <- c("FORELOT", "FOREMON", "FOREBUY", "FOREBAD", "FORFOOD", "FORESUR", "FORENO")

# Subset to forest items 
fore_DF <- DF_anly[names_fore]

# create a DF with the proportions for the first forest variable 
fore_prop <- data.frame(Variable = names_fore[1], cbind(t(prop.table(table(fore_DF[1])))))

# For loop appending the proportions for the remaining variables 
for (i in seq_along(2:(length(names_fore)))){
  fore_prop <- rbind(fore_prop, data.frame(Variable = names_fore[i+1], cbind(t(prop.table(table(fore_DF[i+1]))))))
}

# Likert scaled plot 
HH::likert(Variable~.,fore_prop, positive.order=TRUE,as.percent = TRUE,
           main="Forest instrument",
           xlab="Percentage", ylab="Variable")

# Convert from character to numeric - positively coded items 
table(fore_DF[names_fore][1])
fore_DF[names_fore] <- apply(fore_DF[names_fore], 2, agree.disagree.pos)
table(fore_DF[names_fore][1])

# Parallel analysis 
fa.parallel(fore_DF[complete.cases(fore_DF),][train,], cor = "poly", fm="wls", fa="fa",   main = "Parallel analysis")

# Factor analysis, with polychoric correlation, WLS estimator, and oblimin rotation 
fa_mod2 <- efaUnrotate(fore_DF[train,], estimator = "WLS", nf = 2)

# Factor loadings 
summary(fa_mod2)

# RMSEA
fitmeasures(fa_mod2)["rmsea"]

# Confirmatory analysis with test data 
CFA_2 <- '
factor.1 =~ FOREMON + FOREBUY + FOREBAD + FORFOOD + FORESUR + FORENO
factor.2 =~ FORELOT + FOREBAD + FORESUR 
'
CFA_mod2 <- cfa(CFA_2, data = fore_DF[-train,])

# Factor loadings 
summary(CFA_mod2)

# RMSEA
fitmeasures(CFA_mod2)["rmsea"]

### Food Insecurity Experience Scale (FIES) ### 
FIES_vars <- c("FIESONE","FIESTWO","FIESTHR","FIESFOU","FIESFIV","FIESSIX" , "FIESSEV", "FIESEIG")

# Subset to FIES items 
FIES_DF <- DF_anly[FIES_vars]

# Columns to rows 
FIES_DF_l<- gather(FIES_DF)

# Yes = 1
FIES_DF_l$value <- ifelse(FIES_DF_l$value == "Yes", 1, 
                          ifelse(FIES_DF_l$value == "No", 0,
                                 ifelse(is.na(FIES_DF_l$value)== T, NA, 99)))
table(FIES_DF_l$value)

# Response to each item 
ggplot(data=FIES_DF_l, aes(x=factor(key, levels = rev(c(FIES_vars))), y=value)) +
  geom_bar(stat="identity") + coord_flip()

### Social support instrument ###
names_soc <- c("SOCNEED","SOCREAL","SOCEMOT","SOCTHING","SOCTALK","SOCJOY")

# Subset to social items 
soci_DF <- DF_anly[names_soc]

# create a DF with the proportions for the first social variable 
soci_prop <- data.frame(Variable = names_soc[1], cbind(t(prop.table(table(soci_DF[1])))))

# For loop appending the proportions for the remaining variables 
for (i in seq_along(2:(length(names_soc)))){
  soci_prop <- rbind(soci_prop, data.frame(Variable = names_soc[i+1], cbind(t(prop.table(table(soci_DF[i+1]))))))
}

# Likert scaled plot 
HH::likert(Variable~.,soci_prop, positive.order=TRUE,as.percent = TRUE,
           main="Social support instrument",
           xlab="Percentage", ylab="Variable")

# Convert from character to numeric - positively coded items 
table(soci_DF[names_soc][1])
soci_DF[names_soc] <- apply(soci_DF[names_soc], 2, agree.disagree.pos)
table(soci_DF[names_soc][1])

# Parallel analysis 
fa.parallel(soci_DF[complete.cases(soci_DF),][train,], cor = "poly", fm="wls", fa="fa",   main = "Parallel analysis")

# Factor analysis, with polychoric correlation, WLS estimator, and oblimin rotation 
fa_mod3 <- efaUnrotate(soci_DF[train,], estimator = "WLS", nf = 2)

# Factor loadings 
summary(fa_mod3)

# RMSEA
fitmeasures(fa_mod3)["rmsea"]

# The two factor confirmatory model 
CFA_3 <- '
factor.1 =~ SOCNEED + SOCREAL + SOCEMOT  + SOCJOY
factor.2 =~ SOCTALK + SOCTHING + SOCJOY
'
CFA_mod3 <- cfa(CFA_3, data = soci_DF[-train,])

# Factor loadings 
summary(CFA_mod3)

# RMSEA 
fitmeasures(CFA_mod3)["rmsea"]


### Smoking ### 
table(DF_anly$SMOKING)

### Drinking days per week ### 
table(DF_anly$Alcohol)

### PHQ-8 plus 'strong thoughts' ###
names_PHQ8.S <- c("PH9INTERST", "PH9FEEL","PH9TROUBL" , "PH9TIRED","PH9APPETIT","PH9BADABT","PH9CONCEN","PH9MOVING", "STRONG" )
names_PHQ8 <- c("PH9INTERST", "PH9FEEL","PH9TROUBL" , "PH9TIRED","PH9APPETIT","PH9BADABT","PH9CONCEN","PH9MOVING")

# Subset to PHQ8 items & strong thoughts
PHQ8_DF <- DF_anly[names_PHQ8.S]

# How many crossed the diagnostic for referral? (In total, 20 agreed to be referred, and 6 actually went to the hospital.)

# create a DF with the proportions for the first PHQ8 variable 
PHQ8_prop <- data.frame(Variable = names_PHQ8.S[1], cbind(t(prop.table(table(PHQ8_DF[1])))))

# For loop appending the proportions for the remaining variables 
for (i in seq_along(2:(length(names_PHQ8.S)))){
  PHQ8_prop <- rbind(PHQ8_prop, data.frame(Variable = names_PHQ8.S[i+1], cbind(t(prop.table(table(PHQ8_DF[i+1]))))))
}

# Likert scaled plot 
HH::likert(Variable~.,PHQ8_prop, positive.order=TRUE,as.percent = TRUE,
           main="PHQ-8 & strong thoughts instrument",
           xlab="Percentage", ylab="Variable")

# Convert from character to numeric 
PHQ8_DF <- PHQ8_DF[names_PHQ8] # (exclude 'strong thoughts')
table(PHQ8_DF[names_PHQ8][1])
PHQ8_DF[names_PHQ8] <- apply(PHQ8_DF[names_PHQ8], 2, PHQ8.rec.num)
table(PHQ8_DF[names_PHQ8][1])

# Parallel analysis 
fa.parallel(PHQ8_DF[complete.cases(PHQ8_DF),][train,], cor = "poly", fm="wls", fa="fa",   main = "Parallel analysis")

# Factor analysis, with polychoric correlation, WLS estimator, and oblimin rotation 
fa_mod4 <- efaUnrotate(PHQ8_DF[train,], estimator = "WLS", nf = 1)

# Factor loadings 
summary(fa_mod4)

# RMSEA
fitmeasures(fa_mod4)["rmsea"]

# Confirmatory analysis with test data 
CFA_4 <- '
factor.1 =~ PH9INTERST + PH9FEEL + PH9TROUBL + PH9TIRED + PH9APPETIT + PH9BADABT + PH9CONCEN + PH9MOVING
'
CFA_mod4 <- cfa(CFA_4, data = PHQ8_DF[-train,])

# Factor loadings 
summary(CFA_mod4)

# RMSEA - this suggests poor model fit 
fitmeasures(CFA_mod4)["rmsea"]

# Repeat the parallel analysis using all data 
fa.parallel(PHQ8_DF[complete.cases(PHQ8_DF),], cor = "poly", fm="wls", fa="fa",   main = "Parallel analysis")

# Factor analysis, with polychoric correlation, WLS estimator, and oblimin rotation 
fa_mod4b <- efaUnrotate(PHQ8_DF, estimator = "WLS", nf = 1)

# Factor loadings 
summary(fa_mod4b)

# Fit measures - better than the above
fitmeasures(fa_mod4b)["rmsea"]


### Distance to the nearest forest reserve ### 
hist(DF_anly$FR.dist)

### Distance to Budongo ### 
hist(DF_anly$Budongo.dist)


######### 3) Examine patterns of missing data and perform the multiple imputation #########
###### Patterns of missing data ######
# The observations and variables missing data 
vis_miss(DF_anly) + 
  theme(axis.text.x = element_text(angle = 90)) + theme(text = element_text(size=8)) + coord_flip()

# Patterns of missingness between variables
gg_miss_upset(DF_anly, nsets = n_var_miss(DF_anly))

# High-level pattern of missingness 
gg_miss_var(DF_anly, show_pct = TRUE) + theme(text = element_text(size=8))

# Specific percentages of missing data 
round(table(is.na(DF_anly$SOCJOY))[2]/ nrow(DF_anly)*100,1)

### Observations that did not have any land or asset values
DF_anly[names_land][rowSums(is.na(DF_anly[c(names_land, "Land_size", assets)])) == ncol(DF_anly[c(names_land, "Land_size", assets)]), ]

###### Multiple imputation by chain equation ######
### Extract values for multiple imputation ###
# Variables used in the analysis #
stat.vars2 <- c("LOCB","SEX","DOBB", "EDUCAT","MSTATUS","Health", "ECCOCC","BNSGAS",
                "BNSCLOT","BNSGOAT", "BNSMEAL","BNSELEC","BNSBED",
                "BNSSOLA","BNSJERR" , "BNSMATT","BNSTELV","BNSBRICK",
                "BNSMOTO","BNSHOES", "BNSBLAK","BNSBACC","BNSBIKE",
                "BNSTANK","BNSSHOE","BNSROOF","BNSRADI","BNSSAUC",
                "BNSBATT", "BNSTELL", "BNSCHAI","BNSCAR","BNSFRID",
                "BNSCONC","BNSSAV", "BNSSOF","BNSWOD","BNSSOP","Land_sub", "Land_size",
                "GROWSUG","LANDBIG","LANDBAS","LANDSTR","LANDWELL",
                "LANDGOOD","LANDSMALL","FORELOT","FOREMON","FOREBUY","FOREBAD", "FORFOOD","FORESUR","FORENO",
                "LANDFOR","FIESONE","FIESTWO","FIESTHR","FIESFOU","FIESFIV", "FIESSIX","FIESSEV","FIESEIG",
                "SOCNEED", "SOCREAL","SOCEMOT","SOCTHING","SOCTALK","SOCJOY","SMOKING","PH9INTERST","PH9FEEL",
                "PH9TROUBL" , "PH9TIRED","PH9APPETIT","PH9BADABT","PH9CONCEN","PH9MOVING","STRONG",
                "Strain" , "Alcohol","FR.dist",  "C.A.ratio") 
# Drop "Budongo.dist" - highly correlated with DF.dist

# Subset statistical variables 
DF_stat <- DF_anly[stat.vars2]

# Percentage of non-complete cases
round((1-nrow(DF_stat[complete.cases(DF_stat),])/nrow(DF_stat))*100,1)

# Check the variable format
str(DF_stat)

### Assess the method depending on the data type ###
# Factor = Polytomous logistic regression = polyreg
# Binary = Logistic regression = logreg (or pmm)
# Numeric = Predictive mean matching = pmm 
# Ordered factors = Proportional odds model = polr

# All incomplete PHQ-8 responses were removed so the response variable is not imputed

# Binary variables 
binary_vars <- c("SEX", assets, "GROWSUG", "LANDFOR", FIES_vars, "SMOKING")

# Numeric 
numeric_vars <- c("DOBB", "Land_size", "Alcohol", "FR.dist", "C.A.ratio")

# Factors 
factor_vars <- c("LOCB", "MSTATUS", "ECCOCC")

# Ordered factor 
ordered_vars <- c("EDUCAT", "Health", "Land_sub", "Strain", names_land, names_fore, names_soc, names_PHQ8, "STRONG")

# Create dataset containing the variable names and data types 
DF_stat_type <- data.frame(Name = colnames(DF_stat))
DF_stat_type$Method <- ifelse(DF_stat_type$Name %in% binary_vars, "pmm" , # "logreg", 
                              ifelse(DF_stat_type$Name %in% numeric_vars, "pmm", 
                                     ifelse(DF_stat_type$Name %in% factor_vars, "polyreg", 
                                            ifelse(DF_stat_type$Name %in% ordered_vars, "polr", "ERROR"))))
table(DF_stat_type$Method=="ERROR")

### Run the imputation ### 
# 11 imputed datasets 
n.imp = 11

# Perform the imputation 
DF_stat.imp <- mice(data = DF_stat, m = n.imp, meth = DF_stat_type$Method)

# Check for logged events 
head(DF_stat.imp$loggedEvents, 80)

# Summary 
summary(DF_stat.imp)
saveRDS(DF_stat.imp, "HH_locations/DF_stat.imp.rds")

# Incomplete cases
incomplete <- complete.cases(DF_stat)

### Examining imputed data ### 
# Extract the 10 complete imputed datasets 
DF_stat.list <- list()
for (i in seq_along(1:n.imp)){
  DF_stat.list[[i]] <- complete(DF_stat.imp, i)
}


# Examine numeric variables #
stripplot(DF_stat.imp, pch = 20, cex = 1.2)

# Imputed non-continuous variables 
Imp.nc <- c("BNSBACC", "FIESSEV", "ECCOCC", "Land_sub", "LANDBIG", "LANDBAS","LANDSTR","LANDWELL","LANDGOOD","LANDSMALL","FORELOT","FOREMON" ,"FOREBUY","FOREBAD","FORFOOD","FORESUR","FORENO","SOCJOY","Strain" )

# For loop applying the impt_compare function 
for (i in seq_along(1:length(Imp.nc))){
  print(Imp.nc[[i]])
  print(impt_compare(Orignal.DF = DF_stat, Imp.list = DF_stat.list, variable = Imp.nc[i], n.imp = n.imp))
}

table(is.na(DF_stat.list[[1]][assets]))

######### 4) Post imputation processing #########
### Recode the variables ### 
# Binary variables 
binary_vars <- c(assets,  "GROWSUG" , "LANDFOR" , FIES_vars, "SMOKING" )

# Recode the variables
for (i in seq_along(1:length(DF_stat.list))){
  
  # Land instrument - positively coded items 
  print(apply(DF_stat.list[[i]][land_pos], 2, table))
  DF_stat.list[[i]][land_pos] <- apply(DF_stat.list[[i]][land_pos], 2, agree.disagree.pos)
  print(apply(DF_stat.list[[i]][land_pos], 2, table))
  
  # Land instrument - negatively coded items (i.e. reverse code)
  print(apply(DF_stat.list[[i]][land_neg], 2, table))
  DF_stat.list[[i]][land_neg] <- apply(DF_stat.list[[i]][land_neg], 2, agree.disagree.neg)
  print(apply(DF_stat.list[[i]][land_neg], 2, table))
  
  # Forest dependency instrument - positively coded items 
  print(apply(DF_stat.list[[i]][names_fore], 2, table))
  DF_stat.list[[i]][names_fore] <- apply(DF_stat.list[[i]][names_fore], 2, agree.disagree.pos)
  print(apply(DF_stat.list[[i]][names_fore], 2, table))
  
  # Social support - positively coded items 
  print(apply(DF_stat.list[[i]][names_soc], 2, table))
  DF_stat.list[[i]][names_soc] <- apply(DF_stat.list[[i]][names_soc], 2, agree.disagree.pos)
  print(apply(DF_stat.list[[i]][names_soc], 2, table))
  
  # PHQ-8 
  print(apply(DF_stat.list[[i]][names_PHQ8], 2, table))
  DF_stat.list[[i]][names_PHQ8] <- apply(DF_stat.list[[i]][names_PHQ8], 2, PHQ8.rec.num)
  print(apply(DF_stat.list[[i]][names_PHQ8], 2, table))
  
  # Add back in "REF-PATID"
  DF_stat.list[[i]]$`REF-PATID` <- DF_anly$`REF-PATID`
  
  
  # Binary variables to numeric
  print(apply(DF_stat.list[[i]][binary_vars], 2, table))
  DF_stat.list[[i]][binary_vars] <- apply(DF_stat.list[[i]][binary_vars], 2, function(x)
    ifelse(x %in% c("Yes"), 1,
           ifelse(x %in% c("No"), 0,
                  ifelse(is.na(x)== T, NA, 99)))
  )
  print(apply(DF_stat.list[[i]][binary_vars], 2, table))
  
}

# Save the imputed dataset (run from here if the imputed dataset is already created)
saveRDS(DF_stat.list, "HH_locations/DF_stat.list.rds")

######### 5) Construct asset index and extract plausible values #########
# (NB All newly created variables are scaled and centred.)

# Read in the data 
DF_stat.list <- readRDS("HH_locations/DF_stat.list.rds")

### Asset indexes constructed from the asset ownership variables ###
# Excluding some assets 
exclude <- c("BNSGAS" , "BNSFRID",  "BNSTANK", "BNSELEC" , "BNSCAR" , "BNSTELV" ,   "BNSMOTO" , "BNSBATT", "BNSGOAT",  "BNSBACC")
assets_sub <- setdiff(assets, exclude)

# # Conduct the logistic principal component analysis (following https://cran.r-project.org/web/packages/logisticPCA/vignettes/logisticPCA.html)
# # k = number of principal components to return
# # m = value to approximate the saturated model
# # ms = the different approximations to the saturated model m to try
# 
# # Optimise the number of components to extract based on how strongly it is associated with subjective financial strain (1-10 components)
# association_strengh <- data.frame(K = NA, X.Intercept. = NA, Strain.L = NA, Strain.Q = NA, Strain.C = NA)
# for (i in seq_along(1:10)){
#   
#   # Decide which m to use with cross validation 
#   logpca_cv = cv.lpca(DF_stat.list[[1]][assets_sub], ks = i, ms = 1:10)
#   
#   # Fit the logistic PCA using the minimum m 
#   logpca_asset = logisticPCA(DF_stat.list[[1]][assets_sub], k = i, m = which.min(logpca_cv))
#   
#   # Association between the first component score and subjective financial strain
#   strain_asset <- data.frame(Asset = logpca_asset$PCs[,1], Strain = DF_stat.list[[1]]$Strain)
#   association_strengh  <- rbind(association_strengh, data.frame(K = i, t(coef(lm(Asset ~ Strain, strain_asset)))))
# }
# 
# # Select the k with that yielded the highest association with subjective financial strain 
# association_strengh <- association_strengh[c(2:10),] 
# K_op <- association_strengh[association_strengh$Strain.L == max(association_strengh$Strain.L, na.rm = T),]
# K_op$K
# 
# # Decide which m to use with cross validation using k that yielded the strongest estimate 
# logpca_cv = cv.lpca(DF_stat.list[[1]][assets_sub], ks = K_op$K, ms = 1:10)
# plot(logpca_cv)
# 
# # Fit the logistic PCA using the minimum m
# logpca_asset <- logisticPCA(DF_stat.list[[1]][assets_sub], k = K_op$K, m = which.min(logpca_cv))
# 
# # Summary of loadings 
# round(logpca_asset$U, 2)
# cbind(assets_sub, logpca_asset$U < -.1)
# cbind(assets_sub, logpca_asset$U > -.1)
# 
# # Plot against the sum of asset ownership, split into quantiles 
# survey_response_group <- cut(
#   rowSums(DF_stat.list[[1]][assets_sub]),
#   breaks = quantile(rowSums(DF_stat.list[[1]][assets_sub]), c(0, 0.25, 0.5, 0.75, 1)),
#   labels = c("Poorest", "Poor", "Wealthy", "Wealthiest"),
#   right  = FALSE,
#   include.lowest = TRUE)
# plot(logpca_asset, type = "scores") + geom_point(aes(colour = survey_response_group)) 
# 
# # Plot PCA against subjective financial strain (remember the sign of PCA scores can be flipped)
# plot(logpca_asset, type = "scores") + geom_point(aes(colour = DF_stat.list[[1]]$Strain)) 
# 
# # Association between the first component score and subjective financial strain
# strain_asset <- data.frame(Asset = logpca_asset$PCs[,1], Strain = DF_stat.list[[1]]$Strain)
# summary(lm(Asset ~ Strain, strain_asset))

### Specify k and m
# k_opt <- K_op$K
# min_log <- logpca_cv
k_opt <- 9
min_log <- 10

# Repeat this model for each of the ten imputed datasets - extracting the 1st component in each
for (i in seq_along(1:length(DF_stat.list))) {
  # Fit the logistic PCA using the minimum m and optimal k
  logpca_asset <- logisticPCA(DF_stat.list[[i]][assets_sub], k = k_opt, m = min_log)
  DF_stat.list[[i]]$Assets_index <- logpca_asset$PCs[,1]
  
  # Switching the sign of the Asset index since it currently has an unintuitive direction 
  DF_stat.list[[i]]$Assets_index <- DF_stat.list[[i]]$Assets_index*-1
  
  # Scale and centrer 
  DF_stat.list[[i]]$Assets_index <- scale(DF_stat.list[[i]]$Assets_index, center = T, scale = T)
  
  # Create the economic poverty variable - simply the inverse of the asset index 
  DF_stat.list[[i]]$Economic_poverty <- DF_stat.list[[i]]$Assets_index*-1
}

# # Check that switching the sign worked 
# ggplot() + geom_point(aes(x = survey_response_group, y = DF_stat.list[[1]]$Assets_index)) 
# ggplot() + geom_point(aes(x = survey_response_group, y = DF_stat.list[[1]]$Economic_poverty)) 

### FIES scores ###
# Rasch (where item location varies)
FIES_RM_1 <- mirt(DF_stat.list[[1]][FIES_vars], model = 1, "Rasch")

# Two-parameter model (where item location and discrimination vary)
FIES_RM_2 <- mirt(DF_stat.list[[1]][FIES_vars], model = 1, "2PL")

# Two-parameter model (where item location, discrimination, and guessing vary)
FIES_RM_3 <- mirt(DF_stat.list[[1]][FIES_vars], model = 1, "3PL")

# Compare the AIC of the three models 
extract.mirt(FIES_RM_1, "AIC")
extract.mirt(FIES_RM_2, "AIC") # The two-parameter model performs the best 
extract.mirt(FIES_RM_3, "AIC")

# Inspect fit statistics 
round(M2(FIES_RM_2, type = "C2", calcNULL = FALSE),3)

# Item information curve (or 'information and trace lines') and information and SE
plot(FIES_RM_2, type = 'infotrace', facet_items= F)
plot(FIES_RM_2, type = 'infoSE', facet_items= F) 

# Item characteristic curve (or 'item scoring traceline plots')
plot(FIES_RM_2, type = 'trace', facet_items = F) 

# Extract ten sets of plausible values and inspect them 
plau.draws <- length(DF_stat.list)
FIES_PV <- fscores(FIES_RM_2, plausible.draws = plau.draws)
FIES_PV_combined <- data.frame(Iteration = "1", FIES = FIES_PV[[1]] )
for (i in seq_along(2:length(FIES_PV))){
  FIES_PV_combined <- rbind(FIES_PV_combined, data.frame(Iteration = as.factor(i+1), FIES = FIES_PV[[i+1]] ))
}
# Plot
ggplot(FIES_PV_combined, aes(x=FIES, color=Iteration)) +
  geom_density()

# Extract 1 set of plausible values for each of the ten imputed datasets and scale and center
for (i in seq_along(1:length(DF_stat.list))){
  
  # Two-parameter model 
  FIES_RM_2 <- mirt(DF_stat.list[[i]][FIES_vars], model = 1, "2PL")
  
  # Extract 1 draw of plausible values
  FIES_PV <- fscores(FIES_RM_2, plausible.draws = 1)
  
  # Append it to the imputed dataset 
  DF_stat.list[[i]]$FIES_est <- FIES_PV
  
  # Scale and centrer 
  DF_stat.list[[i]]$FIES_est <- scale(DF_stat.list[[i]]$FIES_est, center = T, scale = T)
}

# Check the sign of the latent variable - should be positive
summary(lm(DF_stat.list[[1]]$FIES_est ~ rowSums(DF_stat.list[[1]][FIES_vars])))

# Explore the association between FIES and subjective land size, controlling for child:adult ration
summary(lm(FIES_est ~ as.ordered(DF_stat.list[[1]]$Land_sub) + C.A.ratio, DF_stat.list[[1]]))

### Land size ###
# Graded response model
Land_GRM <- mirt(DF_stat.list[[1]][names_land], model = 1, "graded")

# Inspect fit statistics 
round(M2(Land_GRM, type = "C2", calcNULL = FALSE),3)

# Item information curve (or 'information and trace lines') and information and SE
plot(Land_GRM, type = 'infotrace', facet_items= F)
plot(Land_GRM, type = 'infoSE', facet_items= F) 

# Item characteristic curve (or 'item scoring traceline plots')
plot(Land_GRM, type = 'trace') 

# Extract ten sets of plausible values and inspect them 
Land_PV <- fscores(Land_GRM, plausible.draws = plau.draws)
Land_PV_combined <- data.frame(Iteration = "1", Land = Land_PV[[1]] )
for (i in seq_along(2:length(Land_PV))){
  Land_PV_combined <- rbind(Land_PV_combined, data.frame(Iteration = as.factor(i+1), Land = Land_PV[[i+1]] ))
}
# Plot
ggplot(Land_PV_combined, aes(x=Land, color=Iteration)) +
  geom_density()

# Extract 1 set of plausible values for each of the ten imputed datasets and scale and center
for (i in seq_along(1:length(DF_stat.list))){
  
  # Graded response model 
  Land_GRM <- mirt(DF_stat.list[[i]][names_land], 1, "graded")
  
  # Extract 1 draw of plausible values
  Land_PV <- fscores(Land_GRM, plausible.draws = 1)
  
  # Append it to the imputed dataset 
  DF_stat.list[[i]]$Land_est <- Land_PV
  
  # Scale and centrer 
  DF_stat.list[[i]]$Land_est <- scale(DF_stat.list[[i]]$Land_est, center = T, scale = T)
}

# Check the sign of the latent variable - should be positive
summary(lm(DF_stat.list[[1]]$Land_est ~ rowSums(DF_stat.list[[1]][names_land])))

### Forest dependence ### 
# The factor structure described by model 'CFA_2' above, following the 'mirt' package syntax.
CFA_2_mirt <- '
F1 = 2-7
F2 = 1,4,6
COV = F1*F2'

# Graded response model 
Forest_GRM <- mirt(DF_stat.list[[1]][names_fore], model = CFA_2_mirt, "graded") 

# Inspect fit statistics 
round(M2(Forest_GRM, type = "C2", calcNULL = FALSE),2)

# Item information curve (or 'information and trace lines') and information and SE
plot(Forest_GRM, type = 'infotrace', facet_items= F)
# plot(Forest_GRM, type = 'infoSE', facet_items= F) 

# Item characteristic curve (or 'item scoring traceline plots')
plot(Forest_GRM, type = 'trace', facet_items= T) 

# Extract ten sets of plausible values and inspect them 
Forest_PV <- fscores(Forest_GRM, plausible.draws = plau.draws)
Forest_PV_combined <- data.frame(Iteration = "1", Forest = Forest_PV[[1]] )
for (i in seq_along(2:length(Forest_PV))){
  Forest_PV_combined <- rbind(Forest_PV_combined, data.frame(Iteration = as.factor(i+1), Forest = Forest_PV[[i+1]] ))
}
# Plot
ggplot(Forest_PV_combined, aes(x=Forest.1, color=Iteration)) +
  geom_density()
ggplot(Forest_PV_combined, aes(x=Forest.2, color=Iteration)) +
  geom_density()

# Extract 1 set of plausible values for each of the ten imputed datasets and scale and centrer
for (i in seq_along(1:length(DF_stat.list))){
  
  # Graded response model 
  Forest_GRM <- mirt(DF_stat.list[[i]][names_fore], model = CFA_2_mirt, "graded") 
  
  # Extract 1 draw of plausible values
  Forest_PV <- fscores(Forest_GRM, plausible.draws = 1)
  
  # Append it to the imputed dataset 
  DF_stat.list[[i]]$Forest_est.1 <- Forest_PV[,1] # Extract the first factor 
  DF_stat.list[[i]]$Forest_est.2 <- Forest_PV[,2] # Extract the second factor 
  
  # Scale and centrer 
  DF_stat.list[[i]]$Forest_est.1 <- scale(DF_stat.list[[i]]$Forest_est.1, center = T, scale = T)
  DF_stat.list[[i]]$Forest_est.2 <- scale(DF_stat.list[[i]]$Forest_est.2, center = T, scale = T)
}


# Items 2-7 should be more positively associated with Forest_est.1 and less positively associated with Forest_est.1
summary(lm(DF_stat.list[[1]]$Forest_est.1 ~ rowSums(DF_stat.list[[1]][names_fore[c(2,3,4,5,6,7)]]))) # Specific food and income forest dependence
summary(lm(DF_stat.list[[1]]$Forest_est.2 ~ rowSums(DF_stat.list[[1]][names_fore[c(2,3,4,5,6,7)]])))

# Items 1, 4, & 6 should be more positively associated with Forest_est.2 and less positively associated with Forest_est.1
summary(lm(DF_stat.list[[1]]$Forest_est.1 ~ rowSums(DF_stat.list[[1]][names_fore[c(1,4,6)]])))
summary(lm(DF_stat.list[[1]]$Forest_est.2 ~ rowSums(DF_stat.list[[1]][names_fore[c(1,4,6)]]))) # General forest dependence


### Social support  ### 
# The factor structure described by model 'CFA_3' above, following the 'mirt' package syntax.
CFA_3_mirt <- '
F1 = 1,2,3,6
F2 = 4,5,6
COV = F1*F2'

# Graded response model 
Soci_GRM <- mirt(data = DF_stat.list[[1]][names_soc], model = CFA_3_mirt, "graded")

# Inspect fit statistics 
round(M2(Soci_GRM, type = "C2", calcNULL = FALSE),3)

# Item information curve (or 'information and trace lines') and information and SE
plot(Soci_GRM, type = 'infotrace', facet_items= F)
# plot(Soci_GRM, type = 'infoSE', facet_items= F) 

# Item characteristic curve (or 'item scoring traceline plots')
plot(Soci_GRM, type = 'trace') 

# Extract ten sets of plausible values and inspect them 
Soci_PV <- fscores(Soci_GRM, plausible.draws = plau.draws)
Soci_PV_combined <- data.frame(Iteration = "1", Soci.1 = Soci_PV[[1]][,1], Soci.2 = Soci_PV[[1]][,2])
for (i in seq_along(2:length(Soci_PV))){
  Soci_PV_combined <- rbind(Soci_PV_combined, data.frame(Iteration = as.factor(i+1), Soci.1 = Soci_PV[[i+1]][,1], Soci.2 = Soci_PV[[i+1]][,2] ))
}

# Plot
ggplot(Soci_PV_combined, aes(x=Soci.1, color=Iteration)) +
  geom_density()
ggplot(Soci_PV_combined, aes(x=Soci.2, color=Iteration)) +
  geom_density()

# Extract 1 set of plausible values for each of the ten imputed datasets and scale and center
for (i in seq_along(1:length(DF_stat.list))){
  
  # Graded response model 
  Soci_GRM <- mirt(DF_stat.list[[i]][names_soc], model = CFA_3_mirt, "graded")
  
  # Extract 1 draw of plausible values
  Soci_PV <- fscores(Soci_GRM, plausible.draws = 1)
  
  # Append it to the imputed dataset 
  DF_stat.list[[i]]$Soci_est.1 <- Soci_PV[,1]
  DF_stat.list[[i]]$Soci_est.2 <- Soci_PV[,2]
  
  # Scale and centrer 
  DF_stat.list[[i]]$Soci_est.1 <- scale(DF_stat.list[[i]]$Soci_est.1, center = T, scale = T)
  DF_stat.list[[i]]$Soci_est.2 <- scale(DF_stat.list[[i]]$Soci_est.2, center = T, scale = T)
}

# Check the sign for social support is as expected - should be positive with raw scores 
# Items 1,2,3, & 6 should be more positively associated with Soci_est.1 and less positively associated with Soci_est.2
summary(lm(DF_stat.list[[1]]$Soci_est.1 ~ rowSums(DF_stat.list[[1]][names_soc[c(1:3,6)]]))) # Family and significant other 
summary(lm(DF_stat.list[[1]]$Soci_est.2 ~ rowSums(DF_stat.list[[1]][names_soc[c(1:3,6)]])))

# Items 4,5 and 6 should be more positively associated with Soci_est.2 and less positively associated with Soci_est.1
summary(lm(DF_stat.list[[1]]$Soci_est.1 ~ rowSums(DF_stat.list[[1]][names_soc[c(4,5,6)]]))) 
summary(lm(DF_stat.list[[1]]$Soci_est.2 ~ rowSums(DF_stat.list[[1]][names_soc[c(4,5,6)]]))) # Friends


### PHQ-8 ### 
# PHQ-8 scores of 10 or above: https://www.thelancet.com/journals/lanpub/article/PIIS2468-2667(21)00047-5/fulltext
# response scale ranged from 0 (not at all) to 3 (nearly every day). 

### Function to sum per individual ###
PHQ_func <- function(x) {
  x$sum_phq <- rowSums(x[names_PHQ8]) 
  x$phq_01 <- ifelse(x$sum_phq > 9, 1, 0)  
  return(x)
}
DF_stat.list <- lapply(DF_stat.list, PHQ_func)

# https://pubmed.ncbi.nlm.nih.gov/36264962/#:~:text=The%20pooled%20prevalence%20of%20depression%20was%2030.2%25%20(95%25%20confidence,pandemic%20period%20(48.1%25%20vs.
### Remove scaling ###
as.num <- c("Assets_index","Economic_poverty","FIES_est","Land_est","Forest_est.1","Forest_est.2" ,"Soci_est.1","Soci_est.2")

# Function to convert specified columns to numeric 
colm_to_nume <- function(data, columns) {
  data[, columns] <- lapply(data[, columns], as.numeric)
  return(data)
}

# Apply function 
DF_stat.list <- lapply(DF_stat.list, colm_to_nume, as.num)

### Add geometry back to dataset ###
geom_fnc <- function(indf){
  indf$geometry <- DF_anly_Sp$geometry
  out_df <- st_as_sf(indf)
  return(out_df)
}
DF_stat.list2 <- lapply(DF_stat.list, geom_fnc)

### Subset to variables of interest ###
DF_stat.list3 <- lapply(DF_stat.list2, function(x) return(x[c("phq_01", "FIES_est", "Economic_poverty", 
                                                        "Health","DOBB","SEX", "EDUCAT","Soci_est.1",
                                                        "MSTATUS",  "Alcohol", "SMOKING", "LOCB", "Land_est", "Forest_est.1")]))

# Convert ordinal to factors 
ord_fact_func <- function(input, cont_vars_ord) {
  output <- st_drop_geometry(input)
  output[cont_vars_ord] <- lapply(output[cont_vars_ord],  function(x) factor(x , ordered = FALSE ))
  st_geometry(output) <- input$geom
  return(output)
}

# Convert ordinal to factors 
fact <- c("Health", "EDUCAT")
DF_stat.list3 <- lapply(DF_stat.list3, ord_fact_func, fact)

# Convert ordinal to age category 
cat_age <- function(input) {
  input$age_cat <- as.factor(ifelse(input$DOBB  < 25, "Under 25",
                                ifelse(input$DOBB >= 25 & input$DOBB <= 39, "25 to 39",
                                       ifelse(input$DOBB >= 40 & input$DOBB <= 55, "40 to 55", "Over 55"))))
  return(input)
}
DF_stat.list3  <- lapply(DF_stat.list3, cat_age)

# Create dummy variables
dummy_function <- function(model_test) {
  model_test <- cbind(model_test, to_dummy(model_test$Health, suffix = "label"))
  model_test <- cbind(model_test, to_dummy(model_test$SEX, suffix = "label"))
  model_test <- cbind(model_test, to_dummy(model_test$EDUCAT, suffix = "label"))
  model_test <- cbind(model_test, to_dummy(model_test$MSTATUS, suffix = "label"))
  model_test <- cbind(model_test, to_dummy(model_test$age_cat, suffix = "label"))
  
  return(model_test)
  
}

# Create dummy variables 
DF_stat.list3 <- lapply(DF_stat.list3, dummy_function)

# Save the final analysis dataset 
saveRDS(DF_stat.list3, "HH_locations/DF_stat.list3.rds")

