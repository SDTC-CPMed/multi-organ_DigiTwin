library(missForest)
library(glmmTMB)
library(MASS)
library(DHARMa)
library(ggplot2)
#library(mixlm)
library(readxl)

#Patients = read.csv('Lupus_data_supplement.csv', row.names = 1)
#Patients = read.csv('Data S16.csv', row.names = 1)

# Read data
Patients_clinical = read_excel('Lupus_data_supplement.xlsx', sheet = 'clinical')
Patients_proteins = read_excel('Lupus_data_supplement.xlsx', sheet = 'proteins')
Patients_info = read_excel('Lupus_data_supplement.xlsx', sheet = 'patients')

# Remove IL17 and CSF as mentioned in manuscript
Patients_proteins = subset(Patients_proteins, select = -c(IL17,CSF))

# Merge all the files
Patients <- merge(Patients_clinical,Patients_proteins,by=c('Inclusion_number', 'Visit_number'))
Patients <- merge(Patients,Patients_info,by=c('Inclusion_number', 'Visit_number'))

# Restructure the file
Patients$Cortisone_dosage <- as.numeric(Patients$Cortisone_dosage)


summary_ACR = list()

for(i in c(1,6,7,9,10)){
  print(i)
  set.seed(123456789)
  
  Patients_Features = Patients[c(4:25,33:49,51:65)]
  Patients_Features_logged = Patients_Features
  Patients_Features_logged[,24:39] = log(Patients_Features[,24:39]+1)
  Patients_Features_imputed = missForest(Patients_Features_logged[,c(3:41,53:54)])$ximp
  
  
  Patients_proteins = Patients_Features[i+41]
  Patients_proteins = cbind(Patients_proteins, Patients_Features_imputed)
  colnames(Patients_proteins)[1] = 'y'
  Patients_proteins = Patients_proteins[!is.na(Patients_proteins[,1]),]
  Patients_proteins[2:42] = scale(Patients_proteins[2:42])
  
  
  Patients_proteins.bin <- glm(y ~ ., family = binomial(link='logit'), data = Patients_proteins)
  Patients_proteins.AIC = stepAIC(Patients_proteins.bin, direction = 'backward', k =10, steps = 1, trace = 0)
  pvals = summary(Patients_proteins.AIC)$coefficients[,4]
  while(max(pvals[2:length(pvals)]) > 0.05){
    Patients_proteins.AIC = stepAIC(Patients_proteins.AIC, direction = 'backward', k =10, steps = 1, trace = 0)
    pvals = summary(Patients_proteins.AIC)$coefficients[,4]
  }

  
  simulationOutput <- simulateResiduals(fittedModel = Patients_proteins.AIC, plot = T)
  testZeroInflation(simulationOutput)
  plot(residuals(simulationOutput)) #should be Uniform
  
  
  print(mean(Patients_proteins$y)/var(Patients_proteins$y)) #NB model is good choice if variance is greater than mean
  print(summary(Patients_proteins.AIC))
  
  plot(Patients_proteins.AIC$fitted.values, Patients_proteins$y)
  
  summary_ACR = append(summary_ACR, list(summary(Patients_proteins.AIC)))

}


#SLICC.ACR.DI.skadeindex.
ACR = 7
summary_SLICC = list()
simulationOutput_SLICC = list()
plot_SLICC = list()
for(status in c(0,1)){
  set.seed(1)
  Patients_Subset = Patients[which(Patients[59] == status),]
  Patients_Features = Patients_Subset[c(4:49, 51,52, 64:65)]
  
  Patients_Features_logged = Patients_Features
  Patients_Features_logged[,31:46] = log(Patients_Features[,31:46]+1)
  Patients_Features_imputed = missForest(Patients_Features_logged[,c(3:50)])$ximp
  
  Patients_proteins = Patients_Features$SDI
  Patients_proteins = cbind('SDI' = Patients_proteins, Patients_Features_imputed)
  Patients_proteins = Patients_proteins[!is.na(Patients_proteins[,1]),]
  Patients_proteins[2:49] = scale(Patients_proteins[2:49])
  
  #Patients_proteins.nb <- glm.nb(SDI ~ ., data = Patients_proteins, control = glm.control(maxit = 100))
  Patients_proteins.nb <- glm.nb(SDI ~ ., data = Patients_proteins)
  #Patients_proteins.AIC = stepAIC(Patients_proteins.nb, direction = 'both', k =3.8415, trace = 0)
  Patients_proteins.AIC = stepAIC(Patients_proteins.nb, direction = 'backward', k =10, steps = 1, trace = 0)
  pvals = summary(Patients_proteins.AIC)$coefficients[,4]
  while(max(pvals[2:length(pvals)]) > 0.05){
    Patients_proteins.AIC = stepAIC(Patients_proteins.AIC, direction = 'backward', k =10, steps = 1, trace = 0)
    pvals = summary(Patients_proteins.AIC)$coefficients[,4]
  }
  
  
  simulationOutput <- simulateResiduals(fittedModel = Patients_proteins.AIC, plot = T)
  testZeroInflation(simulationOutput)
  plot(residuals(simulationOutput)) #should be Uniform
  
  
  print(mean(Patients_proteins$SDI)/var(Patients_proteins$SDI)) #NB model is good choice if variance is greater than mean
  print(summary(Patients_proteins.AIC))
  
  if (status == 0){
    plot_data = data.frame(fitted = Patients_proteins.AIC$fitted.values, 'SDI' = Patients_proteins$SDI)
    col = rep(0, length(Patients_proteins$SDI))
  }
  else{
    plot_data1 = data.frame(fitted = Patients_proteins.AIC$fitted.values, 'SDI' = Patients_proteins$SDI)
    plot_data1 = rbind(plot_data, plot_data1)
    plot_data1$col = c(col, rep(1, length(Patients_proteins$SDI)))
    p = ggplot(plot_data1, aes(fitted,  SDI)) +
      geom_point(aes(colour = factor(col))) + 
      ggtitle(paste('ACR',as.character(ACR)))
    plot_SLICC = append(plot_SLICC, list(p))
  }
  
  summary_SLICC = append(summary_SLICC, list(summary(Patients_proteins.AIC)))
  simulationOutput_SLICC = append(simulationOutput_SLICC, list(simulationOutput))
  
}


#SLEDAI
summary_SLEDAI = list()
simulationOutput_SLEDAI = list()
plot_SLEDAI = list()
for(status in c(0,1)){
  set.seed(1)
  Patients_Subset = Patients[which(Patients[59] == status),]
  Patients_Features = Patients_Subset[c(4:49, 51,52, 64:65)]
  
  Patients_Features_logged = Patients_Features
  Patients_Features_logged[,31:46] = log(Patients_Features[,31:46]+1)
  Patients_Features_imputed = missForest(Patients_Features_logged[,c(3:50)])$ximp
  
  Patients_proteins = Patients_Features$SLEDAI
  Patients_proteins = cbind('SLEDAI' = Patients_proteins, Patients_Features_imputed)
  Patients_proteins = Patients_proteins[!is.na(Patients_proteins[,1]),]
  Patients_proteins[2:49] = scale(Patients_proteins[2:49])
  
  Patients_proteins.zi <- glmmTMB(SLEDAI ~ Weight + BP_systolic + BP_diastolic + Hb +
                                    LCP + TCP + Neutrophils + Eosinophils + Basophils + 
                                    Lymphocytes + Monocytes + Plasma_creatinine + SR +
                                    CRP + C3 + C4 + Urine_leukocytes + Azatioprin +
                                    HQ + MMF + MTX + Rituximab + Sirolimus +
                                    Cortisone_dosage + Smoker + `TNF-alpha` + `IL-6` +
                                    `IL-10` + `IL-27` + `IL1-beta` + `IFN-gamma` + OSM +
                                    `IL1-alpha` + `IL-4` + `IL-2` + `IL-15` + HGF +
                                    FAS + `CD40-L` + CD40 + TGFB1 + Caucasian + Woman +
                                    age + diagnosis_age,
                                  ziformula = ~ 0, family = "nbinom1", data = Patients_proteins)

  Patients_proteins.AIC = stepAIC(Patients_proteins.zi, direction = 'backward', k =10, steps = 1, trace = 0)
  pvals = summary(Patients_proteins.AIC)$coefficients$cond[,4]
  while(max(pvals[2:length(pvals)]) > 0.05){
    Patients_proteins.AIC = stepAIC(Patients_proteins.AIC, direction = 'backward', k =10, steps = 1, trace = 0)
    pvals = summary(Patients_proteins.AIC)$coefficients$cond[,4]
  }
  
  simulationOutput <- simulateResiduals(fittedModel = Patients_proteins.AIC, plot = T)
  testZeroInflation(simulationOutput)
  plot(residuals(simulationOutput)) #should be Uniform
  
  print(mean(Patients_proteins$SLEDAI)/var(Patients_proteins$SLEDAI)) #NB model is good choice if variance is greater than mean
  print(summary(Patients_proteins.AIC))
  
  if (status == 0){
    plot_data = data.frame(fitted = fitted(Patients_proteins.AIC), 'KlinisktSLEDAI' = Patients_proteins$SLEDAI)
    col = rep(0, length(Patients_proteins$SLEDAI))
  }
  else{
    plot_data1 = data.frame(fitted = fitted(Patients_proteins.AIC), 'KlinisktSLEDAI' = Patients_proteins$SLEDAI)
    plot_data1 = rbind(plot_data, plot_data1)
    plot_data1$col = c(col, rep(1, length(Patients_proteins$SLEDAI)))
    p = ggplot(plot_data1, aes(fitted,  SLEDAI)) +
      geom_point(aes(colour = factor(col))) + 
      ggtitle(paste('ACR',as.character(ACR)))
    plot_SLEDAI = append(plot_SLEDAI, list(p))
  }
  summary_SLEDAI = append(summary_SLEDAI, list(summary(Patients_proteins.AIC)))
  simulationOutput_SLEDAI = append(simulationOutput_SLEDAI, list(simulationOutput))
  
}


library("xlsx")

write.xlsx(summary_SLEDAI[[1]]$coefficients$cond, 'SLEDAI_results.xlsx', sheetName = "ACR-7 = 0", 
           col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx(summary_SLEDAI[[2]]$coefficients$cond, 'SLEDAI_results.xlsx', sheetName = "ACR-7 = 1", 
           col.names = TRUE, row.names = TRUE, append = TRUE)


write.xlsx(summary_SLICC[[1]]$coefficients, 'SLICC_results.xlsx', sheetName = "ACR-7 = 0", 
           col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx(summary_SLICC[[2]]$coefficients, 'SLICC_results.xlsx', sheetName = "ACR-7 = 1", 
           col.names = TRUE, row.names = TRUE, append = TRUE)


write.xlsx(summary_ACR[[1]]$coefficients, 'ACR_results.xlsx', sheetName = "ACR1", 
           col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx(summary_ACR[[2]]$coefficients, 'ACR_results.xlsx', sheetName = "ACR6", 
           col.names = TRUE, row.names = TRUE, append = TRUE)
write.xlsx(summary_ACR[[3]]$coefficients, 'ACR_results.xlsx', sheetName = "ACR7", 
           col.names = TRUE, row.names = TRUE, append = TRUE)
write.xlsx(summary_ACR[[4]]$coefficients, 'ACR_results.xlsx', sheetName = "ACR9", 
           col.names = TRUE, row.names = TRUE, append = TRUE)
write.xlsx(summary_ACR[[5]]$coefficients, 'ACR_results.xlsx', sheetName = "ACR10", 
           col.names = TRUE, row.names = TRUE, append = TRUE)

