#### Functions script ####
# Packages 
library(tidyr)


### A function that collapses the paired question into a single variable ###
# DF = The data frame 
# var1 = The name of the variable corresponding to the first part of the two-part question 
# var2 = The name of the variable corresponding to the second part of the question, if the first was answered positively  
# var3 = The name of the variable corresponding to the second part of the question, if the first was answered negatively 
# name = The name of the new variable 
collapse_func <- function(DF, var1, var2, var3, name) {
  
  # Merge the 2nd and 3rd responses 
  DF <- DF %>% unite(new.variable, var2:var3, remove = T, na.rm = T)
  
  # Where there are gaps, use the first response part 
  DF$new.variable <- ifelse(DF$new.variable == "", DF[,var1], DF$new.variable)
  
  # Rename the variable
  names(DF)[names(DF) == "new.variable"] <- name
  
  # Drop the first part of health 
  DF <- DF[ , !(names(DF) %in% c(var1))]
  return(DF)
}

### Recode positively worded agree or disagree items ###
agree.disagree.pos <- function(DF, variable) {
  out <- ifelse(DF[variable] == "Agree.lot", 5, 
                ifelse(DF[variable] == "Agree.lit", 4,
                       ifelse(DF[variable] == "Middle", 3, 
                              ifelse(DF[variable] == "Disagree.lit", 2, 
                                     ifelse(DF[variable] == "Disagree.lot", 1, 
                                            ifelse(DF[variable] %in% c("Refused", "DN", NA), NA,  99))))))
  return(out)
}


### Recode negatively worded agree or disagree items ###
agree.disagree.neg <- function(DF, variable) {
  out <- ifelse(DF[variable] == "Agree.lot", 1, 
                ifelse(DF[variable] == "Agree.lit", 2,
                       ifelse(DF[variable] == "Middle", 3, 
                              ifelse(DF[variable] == "Disagree.lit", 4, 
                                     ifelse(DF[variable] == "Disagree.lot", 5, 
                                            ifelse(DF[variable] %in% c("Refused", "DN", NA), NA,  99))))))
  return(out)
}

### PHQ-8 recode ###
PHQ8.rec <- function(X) {
  out <- ifelse(X == "0", "Not at all", 
                ifelse(X == "1", "Few days",
                       ifelse(X == "2", "More than half the days", 
                              ifelse(X == "3", "Nearly every day", 
                                            ifelse(X %in% c("Refused", "DN", NA), NA,  "99")))))
  return(out)
}

### Education recode function ###
education.rec.num <- function(X) {
  out <- ifelse(X == "No.educ", 1, 
                ifelse(X == "Start.pri", 2,
                       ifelse(X == "Fin.prim", 3, 
                              ifelse(X == "Start.sec", 4, 
                                     ifelse(X == "Fin.sec", 5, 
                                            ifelse(X == "Bey.sec", 6, 
                                            ifelse(X %in% c("Refused", "DN", NA), NA,  99)))))))
  return(out)
}

### Health recode function ###
health.rec.num <- function(X) {
  out <- ifelse(X == "V.bad", 1, 
                ifelse(X == "Bad", 2,
                       ifelse(X == "Fair", 3, 
                              ifelse(X == "Good", 4, 
                                     ifelse(X == "V.good", 5, 
                                            ifelse(X %in% c("Refused", "DN", NA), NA,  99))))))
  return(out)
}

### Subjective land size recode function ###
subj.land.rec.num <- function(X) {
  out <- ifelse(X == "Very.small", 1, 
                ifelse(X == "Small", 2,
                       ifelse(X == "Middle", 3, 
                              ifelse(X == "Large", 4, 
                                     ifelse(X == "Very.large", 5, 
                                            ifelse(X %in% c("Refused", "DN", NA), NA,  99))))))
  return(out)
}

### PHQ-8 recode  - numeric ###
PHQ8.rec.num <- function(X) {
  out <- ifelse(X == "Not at all", 0, 
                ifelse(X ==  "Few days", 1,
                       ifelse(X == "More than half the days", 2, 
                              ifelse(X == "Nearly every day", 3, 
                                     ifelse(X %in% c("Refused", "DN", NA), NA,  99)))))
  return(out)
}


### FIES recode - numeric ###
FIES.rec.num <- function(X) {
  out <- ifelse(X == "Yes", 1, 
                ifelse(X ==  "No", 0,
                       ifelse(X %in% c("Refused", "DN", NA), NA,  99)))
  return(out)
}


### Comparing original vs imputed valued 
impt_compare <- function(Orignal.DF, Imp.list, variable, n.imp){
  
  # Original DF variable
  table_comp <- data.frame(Name = "Original",  data.frame(rbind(table(Orignal.DF[variable]))))
  
  # Loop appending the imputed DF variable 
  for (i in seq_along(1:n.imp)){
    table_comp <- rbind(table_comp,  data.frame(Name = paste0("Iteration ", i),  data.frame(rbind(table(Imp.list[[i]][, variable])))))
  }
  return(table_comp)
}
