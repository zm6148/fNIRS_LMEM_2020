library(ggplot2)
library(reshape2)
library(dplyr)
library(lme4)
library(scales)
library(zoo)

# First remove all objects
rm(list=ls())

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Load the files
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source(paste(dirname(rstudioapi::getSourceEditorContext()$path), '/functions/get_pvals.R', sep = ""))

full_data <- read.csv(file=sprintf(paste(dirname(rstudioapi::getSourceEditorContext()$path), '/LMEM R data exp2/R_data_speech_speechoppo_HbO_HbR.csv', sep = ""), header=FALSE, sep=","))

# build data frame
names(full_data) <- c("SID", "hemisphere", "cortical_structure", "roi_code", "masker_configuration", 
                      "L_R_hand", "R_ear_PTA", "L_ear_PTA",
                      "fnirs_data", "fnirs_data_RC_HbO",
                      "HRF_HbO_1", "HRF_HbO_2", 
                      "HRF_HbO_first_d_1", "HRF_HbO_first_d_2", 
                      "HRF_HbO_second_d_1", "HRF_HbO_second_d_2", 
                      "blk_num", "d_prime", "time",
                      "fnirs_data_RC_HbR", "HRF_HbR_1", "HRF_HbR_2", 
                      "HRF_HbR_first_d_1", "HRF_HbR_first_d_2", 
                      "HRF_HbR_second_d_1", "HRF_HbR_second_d_2",
                      "HbO_HbR")


# convert to factors
full_data$SID <- as.factor(full_data$SID)
full_data$hemisphere <- as.factor(full_data$hemisphere)
full_data$cortical_structure <- as.factor(full_data$cortical_structure)
full_data$roi_code <- as.factor(full_data$roi_code)
full_data$masker_configuration <- as.factor(full_data$masker_configuration)
full_data$L_R_hand <- as.factor(full_data$L_R_hand)
full_data$HbO_HbR <- as.factor(full_data$HbO_HbR)

# rescale to between -1 and 1 for model fit
full_data$HRF_HbO_first_d_1 <- rescale(full_data$HRF_HbO_first_d_1, to = c(-1, 1))
full_data$HRF_HbO_first_d_2 <- rescale(full_data$HRF_HbO_first_d_2, to = c(-1, 1))
full_data$HRF_HbO_second_d_1 <- rescale(full_data$HRF_HbO_second_d_1, to = c(-1, 1))
full_data$HRF_HbO_second_d_2 <- rescale(full_data$HRF_HbO_second_d_2, to = c(-1, 1))
full_data$HRF_HbR_first_d_1 <- rescale(full_data$HRF_HbR_first_d_1, to = c(-1, 1))
full_data$HRF_HbR_first_d_2 <- rescale(full_data$HRF_HbR_first_d_2, to = c(-1, 1))
full_data$HRF_HbR_second_d_1 <- rescale(full_data$HRF_HbR_second_d_1, to = c(-1, 1))
full_data$HRF_HbR_second_d_2 <- rescale(full_data$HRF_HbR_second_d_2, to = c(-1, 1))


# combine 2 masker_configuration and build model
# HbO
full_data$new_HRF_HbO = full_data$HRF_HbO_1 +full_data$HRF_HbO_2
full_data$new_HRF_HbO_first_d = full_data$HRF_HbO_first_d_1 + full_data$HRF_HbO_first_d_2
full_data$new_HRF_HbO_second_d = full_data$HRF_HbO_second_d_1 +  full_data$HRF_HbO_second_d_2
# HbR
full_data$new_HRF_HbR = full_data$HRF_HbR_1 +full_data$HRF_HbR_2
full_data$new_HRF_HbR_first_d = full_data$HRF_HbR_first_d_1 + full_data$HRF_HbR_first_d_2
full_data$new_HRF_HbR_second_d = full_data$HRF_HbR_second_d_1 +  full_data$HRF_HbR_second_d_2


LMEM_model <- lmer(fnirs_data ~
                         # Default effects
                         new_HRF_HbO +
                         new_HRF_HbO_first_d + 
                         new_HRF_HbR +
                         new_HRF_HbR_first_d + 
                         blk_num +
                         fnirs_data_RC_HbO +
                         fnirs_data_RC_HbR +
                         hemisphere +
                         cortical_structure +
                         masker_configuration +
                         #L_R_hand +
                         R_ear_PTA +
                         L_ear_PTA +
                         # two-way interactions
                         cortical_structure:masker_configuration +
                         hemisphere:masker_configuration +
                         hemisphere:cortical_structure +
                         new_HRF_HbO:masker_configuration +
                         new_HRF_HbO:cortical_structure +
                         new_HRF_HbO:hemisphere +
                         new_HRF_HbO_first_d:masker_configuration + 
                         new_HRF_HbO_first_d:cortical_structure + 
                         new_HRF_HbO_first_d:hemisphere + 
                         new_HRF_HbR:masker_configuration +
                         new_HRF_HbR:cortical_structure +
                         new_HRF_HbR:hemisphere +
                         new_HRF_HbR_first_d:masker_configuration + 
                         new_HRF_HbR_first_d:cortical_structure +
                         new_HRF_HbR_first_d:hemisphere +
                         # Random effects
                         #-----------------------------------#
                         (1 + masker_configuration + cortical_structure + hemisphere|SID),
                         #-----------------------------------#
                       data=full_data)

#########################################################################################
## Table 2
#########################################################################################
coefs_1 <- get_pvals(LMEM_model)
coefs_1

#########################################################################################
## figure 2
## plot model predicted trace over real recording
## roi_code 1: right STG; 2: left STG; 3: right cIFS; 4: left cIFS
## masker_configuration 1: speech ; 2: speechoppo
#########################################################################################
# add fitted valued to data frame
full_data$predicted <- fitted(LMEM_model)

# summarize first
df_LMEM_summary <- full_data %>%
  ungroup() %>%
  group_by(hemisphere, cortical_structure, roi_code, masker_configuration, HbO_HbR, time) %>%
  summarise(num_listeners = n(),
            fnirs_data.se = sd(fnirs_data, na.rm = TRUE)/sqrt(n()),
            fnirs_data = round(mean(fnirs_data, na.rm = TRUE),4),
            predicted = mean(predicted))

df_LMEM_summary_HbO = df_LMEM_summary[df_LMEM_summary$HbO_HbR == 1, ] 
df_LMEM_summary_HbR = df_LMEM_summary[df_LMEM_summary$HbO_HbR == 2, ] 

# Plot model prediction on top of the real data
# HbO
ggplot(df_LMEM_summary_HbO)+
  aes(x = time, y = fnirs_data, color = masker_configuration)+
  geom_ribbon(aes(ymin = fnirs_data - fnirs_data.se,
                  ymax = fnirs_data + fnirs_data.se,
                  fill = masker_configuration),
              alpha = 0.5, color = NA)+
  # actual data
  geom_line(aes(y = fnirs_data, group = masker_configuration), size = 1.5)+
  # model prediction
  geom_line(aes(y = predicted, group = masker_configuration), color = "white", size = 1.5)+
  geom_line(aes(y = predicted, group = masker_configuration), size = 0.75, linetype = "dashed")+
  ylim(-1.6, 1.6) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.margin = margin(2,.5,2,.5, "cm")) +
  ggtitle("LMEM fits and raw recordings (HbO)") + 
  facet_grid(. ~ roi_code)

# HbR
ggplot(df_LMEM_summary_HbR)+
  aes(x = time, y = fnirs_data, color = masker_configuration)+
  geom_ribbon(aes(ymin = fnirs_data - fnirs_data.se,
                  ymax = fnirs_data + fnirs_data.se,
                  fill = masker_configuration),
              alpha = 0.5, color = NA)+
  # actual data
  geom_line(aes(y = fnirs_data, group = masker_configuration), size = 1.5)+
  # model prediction
  geom_line(aes(y = predicted, group = masker_configuration), color = "white", size = 1.5)+
  geom_line(aes(y = predicted, group = masker_configuration), size = 0.75, linetype = "dashed")+
  ylim(-1.6, 1.6) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.margin = margin(2,.5,2,.5, "cm")) +
  ggtitle("LMEM fits and raw recordings (HbR)") + 
  facet_grid(. ~ roi_code)

# model goodness of fit measured in R2 
model_fit_r2 = cor(full_data$fnirs_data, full_data$predicted)^2
print(model_fit_r2)

#########################################################################################
## calculate the percetage of reference channel attributed to the
## total activation levels in the LMEM as measured using area under the curve (auc)
#########################################################################################
df_auc <- full_data %>%
  ungroup() %>%
  group_by(time) %>%
  summarise(num_listeners = n(),
            fnirs_data.se = sd(fnirs_data, na.rm = TRUE)/sqrt(n()),
            fnirs_data = round(mean(fnirs_data, na.rm = TRUE),4),
            predicted = mean(predicted),
            fnirs_data_RC_HbO = mean(fnirs_data_RC_HbO),
            fnirs_data_RC_HbR = mean(fnirs_data_RC_HbR),
            fnirs_data_RC = (fnirs_data_RC_HbO +fnirs_data_RC_HbR)/2)

time_axis = abs(df_auc$time)
recording = abs(df_auc$fnirs_data)
LMEM = abs(df_auc$predicted)
RC = abs(df_auc$fnirs_data_RC)
id = order(time_axis)

AUC_LMEM= sum(diff(time_axis[id])*rollmean(LMEM[id],2))
AUC_RC = sum(diff(time_axis[id])*rollmean(RC[id],2))
AUC_percentage = AUC_RC/AUC_LMEM
print(AUC_percentage)

#########################################################################################
## build data for figure 3
## remove reference channel data and find HbO peak during task block
## saved in csv for plot using matlab later
#########################################################################################
full_data_S_RC <- full_data
# remove reference from prediction
full_data_S_RC$predicted <- full_data_S_RC$predicted - full_data$fnirs_data_RC_HbO *  coefs_1$Estimate[7]
full_data_S_RC$predicted <- full_data_S_RC$predicted - full_data$fnirs_data_RC_HbR * coefs_1$Estimate[8]

# calcualte correlation between Speech detection sensitivity (d') and Peak of HbO response
# summarize first
d_prime_vs_HbO_peak <- full_data_S_RC%>%
  ungroup() %>%
  group_by(SID, hemisphere, cortical_structure, roi_code, masker_configuration, d_prime, HbO_HbR, time) %>%
  summarise(num_listeners = n(),
            fnirs_data.se = sd(fnirs_data, na.rm = TRUE)/sqrt(n()),
            fnirs_data = round(mean(fnirs_data, na.rm = TRUE),4),
            predicted = mean(predicted))


## loop through all 14 subjets and record max
##roi_code 1: right STG; 2: left STG; 3: right cIFS; 4: left cIFS
## masker_configuration 1: speech ; 2: speechoppo
## masker_configuration 1 speech, HbO, 4 region
# right STG
peak_R_STG_speech_HbO <- vector()
d_prime_speech_HbO <- vector()
R_STG_speech_HbO <- vector()
for (subject_ID in c(1:14)){
  
  d_prime_vs_HbO_peak_HbO <- d_prime_vs_HbO_peak[d_prime_vs_HbO_peak$HbO_HbR == 1 &
                                                         d_prime_vs_HbO_peak$SID == subject_ID &
                                                         d_prime_vs_HbO_peak$masker_configuration == 1 &
                                                         d_prime_vs_HbO_peak$roi_code == 3, ] 
  
  peak_R_STG_speech_HbO <- c(peak_R_STG_speech_HbO, abs(max(d_prime_vs_HbO_peak_HbO$predicted)))
  d_prime_speech_HbO <- c(d_prime_speech_HbO , d_prime_vs_HbO_peak_HbO$d_prime[1])
  R_STG_speech_HbO <- rbind(R_STG_speech_HbO, d_prime_vs_HbO_peak_HbO$predicted)
}



# left STG
peak_L_STG_speech_HbO <- vector()
d_prime_speech_HbO <- vector()
L_STG_speech_HbO <- vector()
for (subject_ID in c(1:14)){
  
  d_prime_vs_HbO_peak_HbO <- d_prime_vs_HbO_peak[d_prime_vs_HbO_peak$HbO_HbR == 1 &
                                                         d_prime_vs_HbO_peak$SID == subject_ID &
                                                         d_prime_vs_HbO_peak$masker_configuration == 1 &
                                                         d_prime_vs_HbO_peak$roi_code == 4, ] 
  
  peak_L_STG_speech_HbO <- c(peak_L_STG_speech_HbO, abs(max(d_prime_vs_HbO_peak_HbO$predicted)))
  d_prime_speech_HbO <- c(d_prime_speech_HbO , d_prime_vs_HbO_peak_HbO$d_prime[1])
  L_STG_speech_HbO <- rbind(L_STG_speech_HbO, d_prime_vs_HbO_peak_HbO$predicted)
  
}


# right cIFS
peak_R_cIFS_speech_HbO <- vector()
d_prime_speech_HbO <- vector()
R_cIFS_speech_HbO <- vector()
for (subject_ID in c(1:14)){
  
  d_prime_vs_HbO_peak_HbO <- d_prime_vs_HbO_peak[d_prime_vs_HbO_peak$HbO_HbR == 1 &
                                                         d_prime_vs_HbO_peak$SID == subject_ID &
                                                         d_prime_vs_HbO_peak$masker_configuration == 1 &
                                                         d_prime_vs_HbO_peak$roi_code == 1, ] 
  
  peak_R_cIFS_speech_HbO <- c(peak_R_cIFS_speech_HbO, abs(max(d_prime_vs_HbO_peak_HbO$predicted)))
  d_prime_speech_HbO <- c(d_prime_speech_HbO , d_prime_vs_HbO_peak_HbO$d_prime[1])
  R_cIFS_speech_HbO <- rbind(R_cIFS_speech_HbO, d_prime_vs_HbO_peak_HbO$predicted)
  
}

# left cIFS
peak_L_cIFS_speech_HbO <- vector()
d_prime_speech_HbO <- vector()
L_cIFS_speech_HbO <- vector()
for (subject_ID in c(1:14)){
  
  d_prime_vs_HbO_peak_HbO <- d_prime_vs_HbO_peak[d_prime_vs_HbO_peak$HbO_HbR == 1 &
                                                         d_prime_vs_HbO_peak$SID == subject_ID &
                                                         d_prime_vs_HbO_peak$masker_configuration == 1 &
                                                         d_prime_vs_HbO_peak$roi_code == 2, ] 
  
  peak_L_cIFS_speech_HbO <- c(peak_L_cIFS_speech_HbO, abs(max(d_prime_vs_HbO_peak_HbO$predicted)))
  d_prime_speech_HbO<- c(d_prime_speech_HbO , d_prime_vs_HbO_peak_HbO$d_prime[1])
  L_cIFS_speech_HbO <- rbind(L_cIFS_speech_HbO, d_prime_vs_HbO_peak_HbO$predicted)
  
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# masker_configuration 2 speech, HbO, 4 region
# right STG
peak_R_STG_speechoppo_HbO <- vector()
d_prime_speechoppo_HbO <- vector()
R_STG_speechoppo_HbO <- vector()
for (subject_ID in c(1:14)){
  
  d_prime_vs_HbO_peak_HbO <- d_prime_vs_HbO_peak[d_prime_vs_HbO_peak$HbO_HbR == 1 &
                                                         d_prime_vs_HbO_peak$SID == subject_ID &
                                                         d_prime_vs_HbO_peak$masker_configuration == 2 &
                                                         d_prime_vs_HbO_peak$roi_code == 3, ] 
  
  peak_R_STG_speechoppo_HbO <- c(peak_R_STG_speechoppo_HbO, abs(max(d_prime_vs_HbO_peak_HbO$predicted)))
  d_prime_speechoppo_HbO <- c(d_prime_speechoppo_HbO , d_prime_vs_HbO_peak_HbO$d_prime[1])
  R_STG_speechoppo_HbO <- rbind(R_STG_speechoppo_HbO, d_prime_vs_HbO_peak_HbO$predicted)
  
}


# left STG
peak_L_STG_speechoppo_HbO <- vector()
d_prime_speechoppo_HbO <- vector()
L_STG_speechoppo_HbO <- vector()
for (subject_ID in c(1:14)){
  
  d_prime_vs_HbO_peak_HbO <- d_prime_vs_HbO_peak[d_prime_vs_HbO_peak$HbO_HbR == 1 &
                                                         d_prime_vs_HbO_peak$SID == subject_ID &
                                                         d_prime_vs_HbO_peak$masker_configuration == 2 &
                                                         d_prime_vs_HbO_peak$roi_code == 4, ] 
  
  peak_L_STG_speechoppo_HbO <- c(peak_L_STG_speechoppo_HbO, abs(max(d_prime_vs_HbO_peak_HbO$predicted)))
  d_prime_speechoppo_HbO <- c(d_prime_speechoppo_HbO , d_prime_vs_HbO_peak_HbO$d_prime[1])
  L_STG_speechoppo_HbO <- rbind(L_STG_speechoppo_HbO, d_prime_vs_HbO_peak_HbO$predicted)
  
}



# right cIFS
peak_R_cIFS_speechoppo_HbO <- vector()
d_prime_speechoppo_HbO <- vector()
R_cIFS_speechoppo_HbO <-vector()
for (subject_ID in c(1:14)){
  
  d_prime_vs_HbO_peak_HbO <- d_prime_vs_HbO_peak[d_prime_vs_HbO_peak$HbO_HbR == 1 &
                                                         d_prime_vs_HbO_peak$SID == subject_ID &
                                                         d_prime_vs_HbO_peak$masker_configuration == 2 &
                                                         d_prime_vs_HbO_peak$roi_code == 1, ] 
  
  peak_R_cIFS_speechoppo_HbO <- c(peak_R_cIFS_speechoppo_HbO, abs(max(d_prime_vs_HbO_peak_HbO$predicted)))
  d_prime_speechoppo_HbO <- c(d_prime_speechoppo_HbO , d_prime_vs_HbO_peak_HbO$d_prime[1])
  R_cIFS_speechoppo_HbO <- rbind(R_cIFS_speechoppo_HbO, d_prime_vs_HbO_peak_HbO$predicted)
  
}



# left cIFS
peak_L_cIFS_speechoppo_HbO <- vector()
d_prime_speechoppo_HbO <- vector()
L_cIFS_speechoppo_HbO <- vector()
for (subject_ID in c(1:14)){
  
  d_prime_vs_HbO_peak_HbO <- d_prime_vs_HbO_peak[d_prime_vs_HbO_peak$HbO_HbR == 1 &
                                                         d_prime_vs_HbO_peak$SID == subject_ID &
                                                         d_prime_vs_HbO_peak$masker_configuration == 2 &
                                                         d_prime_vs_HbO_peak$roi_code == 2, ] 
  
  peak_L_cIFS_speechoppo_HbO <- c(peak_L_cIFS_speechoppo_HbO, abs(max(d_prime_vs_HbO_peak_HbO$predicted)))
  d_prime_speechoppo_HbO <- c(d_prime_speechoppo_HbO , d_prime_vs_HbO_peak_HbO$d_prime[1])
  L_cIFS_speechoppo_HbO <- rbind(L_cIFS_speechoppo_HbO, d_prime_vs_HbO_peak_HbO$predicted)
  
}



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# masker_configuration 1 speech, HbR, 4 region
# right STG
peak_R_STG_speech_HbR <- vector()
d_prime_speech_HbR<- vector()
R_STG_speech_HbR <- vector()
for (subject_ID in c(1:14)){
  
  d_prime_vs_HbO_peak_HbR <- d_prime_vs_HbO_peak[d_prime_vs_HbO_peak$HbO_HbR == 2 &
                                                         d_prime_vs_HbO_peak$SID == subject_ID &
                                                         d_prime_vs_HbO_peak$masker_configuration == 1 &
                                                         d_prime_vs_HbO_peak$roi_code == 3, ] 
  
  peak_R_STG_speech_HbR <- c(peak_R_STG_speech_HbR, abs(max(d_prime_vs_HbO_peak_HbR$predicted)))
  d_prime_speech_HbR <- c(d_prime_speech_HbR, d_prime_vs_HbO_peak_HbR$d_prime[1])
  R_STG_speech_HbR <- rbind(R_STG_speech_HbR, d_prime_vs_HbO_peak_HbR$predicted)
  
}


# left STG
peak_L_STG_speech_HbR <- vector()
d_prime_speech_HbR <- vector
L_STG_speech_HbR <- vector()
for (subject_ID in c(1:14)){
  
  d_prime_vs_HbO_peak_HbR <- d_prime_vs_HbO_peak[d_prime_vs_HbO_peak$HbO_HbR == 2 &
                                                         d_prime_vs_HbO_peak$SID == subject_ID &
                                                         d_prime_vs_HbO_peak$masker_configuration == 1 &
                                                         d_prime_vs_HbO_peak$roi_code == 4, ] 
  
  peak_L_STG_speech_HbR <- c(peak_L_STG_speech_HbR, abs(max(d_prime_vs_HbO_peak_HbR$predicted)))
  d_prime_speech_HbR <- c(d_prime_speech_HbR, d_prime_vs_HbO_peak_HbR$d_prime[1])
  L_STG_speech_HbR <- rbind(L_STG_speech_HbR, d_prime_vs_HbO_peak_HbR$predicted)
  
}



# right cIFS
peak_R_cIFS_speech_HbR <- vector()
d_prime_speech_HbR <- vector()
R_cIFS_speech_HbR <- vector()
for (subject_ID in c(1:14)){
  
  d_prime_vs_HbO_peak_HbR <- d_prime_vs_HbO_peak[d_prime_vs_HbO_peak$HbO_HbR == 2 &
                                                         d_prime_vs_HbO_peak$SID == subject_ID &
                                                         d_prime_vs_HbO_peak$masker_configuration == 1 &
                                                         d_prime_vs_HbO_peak$roi_code == 1, ] 
  
  peak_R_cIFS_speech_HbR <- c(peak_R_cIFS_speech_HbR, abs(max(d_prime_vs_HbO_peak_HbR$predicted)))
  d_prime_speech_HbR <- c(d_prime_speech_HbR , d_prime_vs_HbO_peak_HbR$d_prime[1])
  R_cIFS_speech_HbR <- rbind(R_cIFS_speech_HbR, d_prime_vs_HbO_peak_HbR$predicted)
  
}


# left cIFS
peak_L_cIFS_speech_HbR <- vector()
d_prime_speech_HbR <- vector()
L_cIFS_speech_HbR <- vector()
for (subject_ID in c(1:14)){
  
  d_prime_vs_HbO_peak_HbR <- d_prime_vs_HbO_peak[d_prime_vs_HbO_peak$HbO_HbR == 2 &
                                                         d_prime_vs_HbO_peak$SID == subject_ID &
                                                         d_prime_vs_HbO_peak$masker_configuration == 1 &
                                                         d_prime_vs_HbO_peak$roi_code == 2, ] 
  
  peak_L_cIFS_speech_HbR <- c(peak_L_cIFS_speech_HbR, abs(max(d_prime_vs_HbO_peak_HbR$predicted)))
  d_prime_speech_HbR <- c(d_prime_speech_HbR , d_prime_vs_HbO_peak_HbR$d_prime[1])
  L_cIFS_speech_HbR <- rbind(L_cIFS_speech_HbR, d_prime_vs_HbO_peak_HbR$predicted)
  
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# masker_configuration 2 speechoppo, HbO, 4 region
# right STG
peak_R_STG_speechoppo_HbR <- vector()
d_prime_speechoppo_HbR <- vector()
R_STG_speechoppo_HbR <- vector()
for (subject_ID in c(1:14)){
  
  d_prime_vs_HbO_peak_HbR <- d_prime_vs_HbO_peak[d_prime_vs_HbO_peak$HbO_HbR == 2 &
                                                         d_prime_vs_HbO_peak$SID == subject_ID &
                                                         d_prime_vs_HbO_peak$masker_configuration == 2 &
                                                         d_prime_vs_HbO_peak$roi_code == 3, ] 
  
  peak_R_STG_speechoppo_HbR <- c(peak_R_STG_speechoppo_HbR, abs(max(d_prime_vs_HbO_peak_HbR$predicted)))
  d_prime_speechoppo_HbR <- c(d_prime_speechoppo_HbR, d_prime_vs_HbO_peak_HbR$d_prime[1])
  R_STG_speechoppo_HbR <- rbind(R_STG_speechoppo_HbR, d_prime_vs_HbO_peak_HbR$predicted)
  
}


# left STG
peak_L_STG_speechoppo_HbR <- vector()
d_prime_speechoppo_HbR <- vector()
L_STG_speechoppo_HbR <- vector()
for (subject_ID in c(1:14)){
  
  d_prime_vs_HbO_peak_HbR <- d_prime_vs_HbO_peak[d_prime_vs_HbO_peak$HbO_HbR == 2 &
                                                         d_prime_vs_HbO_peak$SID == subject_ID &
                                                         d_prime_vs_HbO_peak$masker_configuration == 2 &
                                                         d_prime_vs_HbO_peak$roi_code == 4, ] 
  
  peak_L_STG_speechoppo_HbR <- c(peak_L_STG_speechoppo_HbR, abs(max(d_prime_vs_HbO_peak_HbR$predicted)))
  d_prime_speechoppo_HbR <- c(d_prime_speechoppo_HbR , d_prime_vs_HbO_peak_HbR$d_prime[1])
  L_STG_speechoppo_HbR <- rbind(L_STG_speechoppo_HbR, d_prime_vs_HbO_peak_HbR$predicted)
  
}

# right cIFS
peak_R_cIFS_speechoppo_HbR <- vector()
d_prime_speechoppo_HbR <- vector()
R_cIFS_speechoppo_HbR <- vector()
for (subject_ID in c(1:14)){
  
  d_prime_vs_HbO_peak_HbR <- d_prime_vs_HbO_peak[d_prime_vs_HbO_peak$HbO_HbR == 2 &
                                                         d_prime_vs_HbO_peak$SID == subject_ID &
                                                         d_prime_vs_HbO_peak$masker_configuration == 2 &
                                                         d_prime_vs_HbO_peak$roi_code == 1, ] 
  
  peak_R_cIFS_speechoppo_HbR <- c(peak_R_cIFS_speechoppo_HbR, abs(max(d_prime_vs_HbO_peak_HbR$predicted)))
  d_prime_speechoppo_HbR <- c(d_prime_speechoppo_HbR , d_prime_vs_HbO_peak_HbR$d_prime[1])
  R_cIFS_speechoppo_HbR <- rbind(R_cIFS_speechoppo_HbR, d_prime_vs_HbO_peak_HbR$predicted)
  
}



# left cIFS
peak_L_cIFS_speechoppo_HbR <- vector()
d_prime_speechoppo_HbR <- vector()
L_cIFS_speechoppo_HbR <- vector()
for (subject_ID in c(1:14)){
  
  d_prime_vs_HbO_peak_HbR <- d_prime_vs_HbO_peak[d_prime_vs_HbO_peak$HbO_HbR == 2 &
                                                         d_prime_vs_HbO_peak$SID == subject_ID &
                                                         d_prime_vs_HbO_peak$masker_configuration == 2 &
                                                         d_prime_vs_HbO_peak$roi_code == 2, ] 
  
  peak_L_cIFS_speechoppo_HbR <- c(peak_L_cIFS_speechoppo_HbR, abs(max(d_prime_vs_HbO_peak_HbR$predicted)))
  d_prime_speechoppo_HbR <- c(d_prime_speechoppo_HbR , d_prime_vs_HbO_peak_HbR$d_prime[1])
  L_cIFS_speechoppo_HbR <- rbind(L_cIFS_speechoppo_HbR, d_prime_vs_HbO_peak_HbR$predicted)
  
}


## save data for analysis
#trace data
setwd(paste(dirname(rstudioapi::getSourceEditorContext()$path), '/model only trace exp2', sep = ""))
write.csv(L_cIFS_speech_HbO, file = "L_cIFS_speech_HbO.csv")
write.csv(L_cIFS_speech_HbR, file = "L_cIFS_speech_HbR.csv")
write.csv(L_STG_speech_HbO, file = "L_STG_speech_HbO.csv")
write.csv(L_STG_speech_HbR, file = "L_STG_speech_HbR.csv")
write.csv(R_cIFS_speech_HbO, file = "R_cIFS_speech_HbO.csv")
write.csv(R_cIFS_speech_HbR, file = "R_cIFS_speech_HbR.csv")
write.csv(R_STG_speech_HbO, file = "R_STG_speech_HbO.csv")
write.csv(R_STG_speech_HbR, file = "R_STG_speech_HbR.csv")
write.csv(L_cIFS_speechoppo_HbO, file = "L_cIFS_speechoppo_HbO.csv")
write.csv(L_cIFS_speechoppo_HbR, file = "L_cIFS_speechoppo_HbR.csv")
write.csv(L_STG_speechoppo_HbO, file = "L_STG_speechoppo_HbO.csv")
write.csv(L_STG_speechoppo_HbR, file = "L_STG_speechoppo_HbR.csv")
write.csv(R_cIFS_speechoppo_HbO, file = "R_cIFS_speechoppo_HbO.csv")
write.csv(R_cIFS_speechoppo_HbR, file = "R_cIFS_speechoppo_HbR.csv")
write.csv(R_STG_speechoppo_HbO, file = "R_STG_speechoppo_HbO.csv")
write.csv(R_STG_speechoppo_HbR, file = "R_STG_speechoppo_HbR.csv")
# correlation data
setwd(paste(dirname(rstudioapi::getSourceEditorContext()$path), '/behavioral and hemodynamic responses data exp2', sep = ""))
write.csv(rbind(peak_R_STG_speech_HbO, d_prime_speech_HbO), file = "R_STG_speech_HbO_peak_d.csv")
write.csv(rbind(peak_L_STG_speech_HbO, d_prime_speech_HbO), file = "L_STG_speech_HbO_peak_d.csv")
write.csv(rbind(peak_R_cIFS_speech_HbO, d_prime_speech_HbO), file = "R_cIFS_speech_HbO_peak_d.csv")
write.csv(rbind(peak_L_cIFS_speech_HbO, d_prime_speech_HbO), file = "L_cIFS_speech_HbO_peak_d.csv")
write.csv(rbind(peak_R_STG_speechoppo_HbO, d_prime_speechoppo_HbO), file = "R_STG_speechoppo_HbO_peak_d.csv")
write.csv(rbind(peak_L_STG_speechoppo_HbO, d_prime_speechoppo_HbO), file = "L_STG_speechoppo_HbO_peak_d.csv")
write.csv(rbind(peak_R_cIFS_speechoppo_HbO, d_prime_speechoppo_HbO), file = "R_cIFS_speechoppo_HbO_peak_d.csv")
write.csv(rbind(peak_L_cIFS_speechoppo_HbO, d_prime_speechoppo_HbO), file = "L_cIFS_speechoppo_HbO_peak_d.csv")
write.csv(rbind(peak_R_STG_speech_HbR, d_prime_speech_HbR), file = "R_STG_speech_HbR_peak_d.csv")
write.csv(rbind(peak_L_STG_speech_HbR, d_prime_speech_HbR), file = "L_STG_speech_HbR_peak_d.csv")
write.csv(rbind(peak_R_cIFS_speech_HbR, d_prime_speech_HbR), file = "R_cIFS_speech_HbR_peak_d.csv")
write.csv(rbind(peak_L_cIFS_speech_HbR, d_prime_speech_HbR), file = "L_cIFS_speech_HbR_peak_d.csv")
write.csv(rbind(peak_R_STG_speechoppo_HbR, d_prime_speechoppo_HbR), file = "R_STG_speechoppo_HbR_peak_d.csv")
write.csv(rbind(peak_L_STG_speechoppo_HbR, d_prime_speechoppo_HbR), file = "L_STG_speechoppo_HbR_peak_d.csv")
write.csv(rbind(peak_R_cIFS_speechoppo_HbR, d_prime_speechoppo_HbR), file = "R_cIFS_speechoppo_HbR_peak_d.csv")
write.csv(rbind(peak_L_cIFS_speechoppo_HbR, d_prime_speechoppo_HbR), file = "L_cIFS_speechoppo_HbR_peak_d.csv")