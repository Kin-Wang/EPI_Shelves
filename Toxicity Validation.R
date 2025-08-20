# Toxicity Validation 

# Last modified by: PW
# Date modified: 6/14/2024


# Set WD
# rm(list=ls()[ls() %in% c("a", "b", "check")])

# PC
if (print(Sys.info()["login"])=="wangp3") {
  datapath = "G:/OBCD/Toxicity_Validation/Data/"
  outpath = "G:/OBCD/Toxicity_Validation/Output/"
}

# Mac
if (print(Sys.info()["user"])=="pwang") {
  datapath = "/Users/pwang/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/OBCD/Toxicity_Validation/Data/"
  outpath = "/Users/pwang/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/OBCD/Toxicity_Validation/Output/"
}


# Check WD is correct
datapath
outpath

# Library
library(haven)
library(dplyr)
library(naniar)
library(lubridate)
library(gtsummary)
library(officer)
library(flextable)
library(irr)
library(irrCAC)
library(zoo)
library(tidyr)
library(ggplot2)
library(epiR)


#### Prepare Data for Renal, Hepatic, Anemia, Neutro, Thrombo, and Neuropathy ####
## Renal ##
## Load Data
kpwa_icd <- haven::read_sas(paste0(datapath, "kpwa_comorb_toxic_icd_2024_04_19.sas7bdat"))
kpwa_lab <- haven::read_sas(paste0(datapath, "kpwa_comorbid_lab_2023_02_22.sas7bdat"))
kpnc_icd <- haven::read_sas(paste0(datapath, "kpnc_comorb_tox_icd_2024_04_10.sas7bdat"))
kpnc_lab <- haven::read_sas(paste0(datapath, "kpnc_comorb_tox_lab_2023_04_21.sas7bdat"))
kpwa_combo <- haven::read_sas(paste0(datapath, "kpwa_toxicity_08282023.sas7bdat"))
obcd_cohort <- haven::read_sas(paste0(datapath, "obcd_cohort_15dec23.sas7bdat")) # OBCD Cohort

## Check
# lower case all variable names
names(obcd_cohort) <- tolower(names(obcd_cohort))
names(kpwa_icd) <- tolower(names(kpwa_icd))
names(kpwa_lab) <- tolower(names(kpwa_lab))
names(kpnc_icd) <- tolower(names(kpnc_icd))
names(kpnc_lab) <- tolower(names(kpnc_lab))
names(kpwa_combo) <- tolower(names(kpwa_combo))

# Check dates
anyNA(obcd_cohort$dx_dt) # F
anyNA(kpwa_icd$co_tox_dx_code_date) # T
sum(is.na(kpwa_icd$co_tox_dx_code_date)) # 728
anyNA(kpwa_lab$co_tox_lab_date) # F
anyNA(kpnc_icd$co_tox_dx_code_date) # F
anyNA(kpnc_lab$co_tox_lab_date) # F
sum(is.na(obcd_cohort$first_chemo_dt[obcd_cohort$any_chemo ==1])) 

obcd_cohort <- obcd_cohort[!(obcd_cohort$any_chemo == 1 & is.na(obcd_cohort$first_chemo_dt)),] # delete patients who had chemo but no date
table(obcd_cohort$surgery_d, useNA = "ifany")
table(obcd_cohort$surgery_dt[obcd_cohort$surgery_d == 0], useNA = "ifany")

# Check
glimpse(obcd_cohort) 
glimpse(kpwa_icd) 
glimpse(kpwa_lab) 
glimpse(kpnc_icd) 
glimpse(kpnc_lab) 
glimpse(kpwa_combo) 
table(kpwa_combo$chart_reviewed, useNA = "ifany") 
sum(kpwa_combo$studyid %in% obcd_cohort$studyid) 
sum(kpwa_combo$studyid %in% kpwa_icd$studyid) 
length(unique(kpwa_icd$studyid[!kpwa_icd$studyid %in% obcd_cohort$studyid])) 
length(unique(kpnc_icd$studyid[!kpnc_icd$studyid %in% obcd_cohort$studyid])) 
length(unique(kpwa_lab$studyid[!kpwa_lab$studyid %in% obcd_cohort$studyid])) 
length(unique(kpnc_lab$studyid[!kpnc_lab$studyid %in% obcd_cohort$studyid])) 
sum(is.na(obcd_cohort$final_dx_wt)) 
sum(is.na(obcd_cohort$final_2nd_wt)) 

obcd_cohort[is.na(obcd_cohort$surgery_d), c("studyid", "dx_dt", "surgery_d", "surgery_dt", "any_chemo", "first_chemo_dt")] %>% print(n=Inf)

# Check missing
vis_miss(obcd_cohort) 
vis_miss(kpwa_icd, warn_large_data = FALSE) 
vis_miss(kpwa_lab, warn_large_data = FALSE) 
vis_miss(kpnc_icd, warn_large_data = FALSE) 
vis_miss(kpnc_lab, warn_large_data = FALSE) 
vis_miss(kpwa_combo) 
vis_miss(kpwa_combo[kpwa_combo$chart_reviewed == "Y",]) 
vis_miss(kpwa_combo[kpwa_combo$chart_reviewed == "N",]) 
vis_miss(kpwa_combo[!is.na(kpwa_combo$elig_final) & kpwa_combo$elig_final == 1,]) 

## Merge Data 
icd_data <-  kpwa_icd %>% bind_rows(kpnc_icd) %>%
  left_join(obcd_cohort, by = join_by(studyid, siteid)) %>%
  filter(!is.na(dx_dt)) %>% # patients in icd but not in obcd cohort were excluded
  filter(!is.na(co_tox_dx_code_date)) # patients with missing icd code date
glimpse(icd_data)

lab_data <- kpwa_lab %>% bind_rows(kpnc_lab) %>%
  left_join(obcd_cohort, by = join_by(studyid, siteid)) %>%
  filter(!is.na(dx_dt)) # patients in lab but not in obcd cohort were excluded
glimpse(lab_data)

chart_data <- kpwa_combo %>%
  inner_join(obcd_cohort, by = join_by(studyid)) %>%
  filter(any_chemo == 1)  # for chart review, only patients with chemo were included
glimpse(chart_data) # only chemo side effects; no need to create time frames
table(chart_data$toxicityname, useNA = "ifany")

## Clean Data
length(unique(icd_data$studyid)) 
length(unique(icd_data$studyid[icd_data$any_chemo==1])) 
anyNA(icd_data$co_tox_dx_code_date) # F
table(icd_data$kidney_flag_icd, useNA = "ifany") 
length(unique(lab_data$studyid)) 
length(unique(lab_data$studyid[lab_data$any_chemo==1])) 
anyNA(lab_data$co_tox_lab_date) #F
sum(is.na(lab_data$final_dx_wt)) 
sum(is.na(lab_data$final_2nd_wt)) 
sum(is.na(lab_data$co_tox_lab_result[lab_data$co_tox_lab_name == 11 ])) 
table(lab_data$co_tox_lab_unit[lab_data$co_tox_lab_name==11], useNA = "ifany") 
as.data.frame(lab_data[lab_data$co_tox_lab_name==11 & lab_data$co_tox_lab_unit %in% c(NA, 11), c("studyid", "co_tox_lab_date", "co_tox_lab_result", "co_tox_lab_unit")]) 
as.data.frame(lab_data[lab_data$co_tox_lab_name==11 & lab_data$co_tox_lab_result>20, c("studyid", "co_tox_lab_date", "co_tox_lab_result", "co_tox_lab_unit")]) 
lab_data$co_tox_lab_result[lab_data$studyid == "marked-out id" & lab_data$co_tox_lab_date == ymd("marked-out date")] <- round(121/88.402, digits = 2) # transform umol/L to mg/dL
## there is only mg/dL for serum creatinine now


## Adding Group
# group
# Neither surgery nor chemo
# Only chemo
# Only surgery
# Surgery after chemo
# Surgery before chemo
icd_data <- icd_data %>%
  mutate(group = ifelse(is.na(surgery_d), NA,
                        ifelse(surgery_d == 0 & any_chemo == 0, "Neither",
                               ifelse(surgery_d == 0 & any_chemo == 1, "Only chemo",
                                      ifelse(is.na(surgery_dt), NA,
                                             ifelse(any_chemo == 0, "Only surgery",
                                                    ifelse(surgery_dt < first_chemo_dt, "Surgery before chemo", "Surgery after chemo")))))))
table(icd_data$group, icd_data$any_chemo, useNA = "ifany")
table(icd_data$group, icd_data$surgery_d, useNA = "ifany")


lab_data <- lab_data %>%
  mutate(group = ifelse(is.na(surgery_d), NA,
                        ifelse(surgery_d == 0 & any_chemo == 0, "Neither",
                               ifelse(surgery_d == 0 & any_chemo == 1, "Only chemo",
                                      ifelse(is.na(surgery_dt), NA,
                                             ifelse(any_chemo == 0, "Only surgery",
                                                    ifelse(surgery_dt < first_chemo_dt, "Surgery before chemo", "Surgery after chemo")))))))
table(lab_data$group, lab_data$any_chemo, useNA = "ifany")
table(lab_data$group, lab_data$surgery_d, useNA = "ifany")


obcd_group <- obcd_cohort %>%
  mutate(group = ifelse(is.na(surgery_d), NA,
                        ifelse(surgery_d == 0 & any_chemo == 0, "Neither",
                               ifelse(surgery_d == 0 & any_chemo == 1, "Only chemo",
                                      ifelse(is.na(surgery_dt), NA,
                                             ifelse(any_chemo == 0, "Only surgery",
                                                    ifelse(surgery_dt < first_chemo_dt, "Surgery before chemo", "Surgery after chemo"))))))
         
  )


a <- chart_data[is.na(chart_data$toxicitydate), c("studyid", "toxicityname", "toxicitydate", "dx_dt", "surgery_d", "surgery_dt", "first_chemo_dt")] %>% print(n=Inf, width = Inf)
b <- chart_data[!is.na(chart_data$toxicitydate) & (chart_data$toxicitydate < chart_data$first_chemo_dt), c("studyid", "toxicityid", "toxicityname", "toxicitydate", "dx_dt", "surgery_d", "surgery_dt", "first_chemo_dt", "last_chemo_dt")] %>% print(n=Inf, width = Inf)
c <- chart_data[!is.na(chart_data$toxicitydate) & (chart_data$toxicitydate == chart_data$first_chemo_dt), c("studyid", "toxicityname", "toxicitydate", "dx_dt", "surgery_d", "surgery_dt", "first_chemo_dt")] %>% print(n=Inf, width = Inf)
# exclude 2 patients without toxicity name
chart_data <- chart_data[chart_data$toxicityname != "",]
# exclude 1 patient with date < first chemo and renal toxicity
chart_data <- chart_data[chart_data$toxicityid != "marked-out id",]
# if chemo start at July 2008
table(chart_data$surgery_d, useNA = "ifany")

# export data
# saveRDS(icd_data, file = paste0(datapath, "icd.rds"))
# saveRDS(lab_data, file = paste0(datapath, "lab.rds"))
# saveRDS(chart_data, file = paste0(datapath, "chart.rds"))
# saveRDS(obcd_group, file = paste0(datapath, "obcd.rds"))
# import data
# icd_data <- readRDS(file = paste0(datapath, "icd.rds"))
# lab_data <- readRDS(file = paste0(datapath, "lab.rds"))
# chart_data <- readRDS(file = paste0(datapath, "chart.rds"))
# obcd_group <- readRDS(file = paste0(datapath, "obcd.rds"))


## Adding Time Frame
# time frame for surgery before chemo
# 1: <= date of cancer diagnosis --> pre-existing conditions (1 year)
# 2: <= date of surgery
# 3: <= date of first chemo
# 4: <= end of chemo
icd_data_before <- icd_data %>%
  filter(group == "Surgery before chemo") %>%
  mutate(time_frame = ifelse(co_tox_dx_code_date <= dx_dt, 1,
                             ifelse(co_tox_dx_code_date <= surgery_dt, 2,
                                    ifelse(co_tox_dx_code_date <= first_chemo_dt, 3,
                                           ifelse(co_tox_dx_code_date <= last_chemo_dt, 4, NA))))) %>%
  filter(!is.na(time_frame)) # remove NAs in time frame

table(icd_data_before$time_frame, useNA = "ifany")
table(icd_data_before$kidney_flag_icd, useNA = "ifany") 

lab_data_before <- lab_data %>%
  filter(group == "Surgery before chemo") %>%
  mutate(time_frame = ifelse(co_tox_lab_date <= dx_dt, 1,
                             ifelse(co_tox_lab_date <= surgery_dt, 2,
                                    ifelse(co_tox_lab_date <= first_chemo_dt, 3,
                                           ifelse(co_tox_lab_date <= last_chemo_dt, 4, NA))))) %>%
  filter(!is.na(time_frame)) # remove NAs in time frame

table(lab_data_before$time_frame, useNA = "ifany")
table(lab_data_before$co_tox_lab_name, useNA = "ifany") # 11 = Creatinine_s (serum), 12 = Creatinine_u (urinary), 13 = Creatinine_u24 (urinary, 24 hour), 14 = GFR_NAA_MD, 15 = GFR_AA_MD
sum(is.na(lab_data_before$final_dx_wt)) 
sum(is.na(lab_data_before$final_2nd_wt)) 

# time frame = 1 or 2, use final weight; time frame = 3 or 4, use second weight
lab_data_before$diffyears <- as.numeric(trunc((lab_data_before$co_tox_lab_date - lab_data_before$dx_dt)/365))
lab_data_before$age <- lab_data_before$dx_age + lab_data_before$diffyears

lab_data_before <- lab_data_before %>%
  mutate(CrCl = ifelse(co_tox_lab_name != 11, NA,
                       ifelse(is.na(co_tox_lab_result), NA,
                              ifelse(time_frame < 3 & is.na(final_dx_wt), NA,
                                     ifelse(time_frame > 2 & is.na(final_2nd_wt), NA,
                                            ifelse(time_frame < 3, (140 - age)*final_dx_wt*0.85/(72*co_tox_lab_result),
                                                   ifelse(time_frame > 2, (140 - age)*final_2nd_wt*0.85/(72*co_tox_lab_result), NA))))))) %>%
  mutate(CrCl_15 = ifelse(is.na(CrCl), NA, ifelse(CrCl < 15, 1, 0)),
         CrCl_30 = ifelse(is.na(CrCl), NA, ifelse(CrCl < 30, 1, 0)),
         CrCl_60 = ifelse(is.na(CrCl), NA, ifelse(CrCl < 60, 1, 0)))
summary(lab_data_before$CrCl)
table(lab_data_before$CrCl_15[lab_data_before$co_tox_lab_name == 11], lab_data_before$time_frame[lab_data_before$co_tox_lab_name == 11], useNA = "ifany")

chart_before <- chart_data %>% filter(surgery_d == 1) %>%
  filter(!is.na(surgery_dt)) %>%
  filter(surgery_dt < first_chemo_dt) %>%
  filter(toxicitydate > first_chemo_dt & toxicitydate <= last_chemo_dt) # limit to last timeframe


## Summarize Kidney Toxicity
icd_data_before_sum <- icd_data_before %>% 
  select(studyid, time_frame, kidney_flag_icd) %>%
  group_by(studyid, time_frame) %>%
  summarise(kidney_n = sum(kidney_flag_icd))
icd_data_before_sum <- icd_data_before_sum %>% mutate(kidney_flag = ifelse(kidney_n == 0, 0, 1))
# kidney_flag is individual level indicator by time frame; if patient had at least one kidney icd by time frame, that patient would have 1 in this variable.
icd_data_before_sum$siteid <- ifelse(icd_data_before_sum$studyid>"marked-out id", 2, 1)
table(icd_data_before_sum$kidney_n) # check
icd_data_before_sum[icd_data_before_sum$kidney_n >= 1,] # check

lab_data_before_sum <- lab_data_before %>% 
  filter(co_tox_lab_name == 11) %>%
  select(studyid, time_frame, CrCl, CrCl_15, CrCl_30, CrCl_60) %>%
  group_by(studyid, time_frame) %>%
  summarise(CrCl15 = sum(CrCl_15, na.rm = T),
            CrCl30 = sum(CrCl_30, na.rm = T),
            CrCl60 = sum(CrCl_60, na.rm = T),
            NAs = sum(is.na(CrCl)), N= n())
sum((lab_data_before_sum$NAs == lab_data_before_sum$N))

lab_data_before_sum <- lab_data_before_sum %>% 
  mutate(kidney_flag_15 = ifelse(NAs == N, NA, ifelse(CrCl15 == 0, 0, 1)),
         kidney_flag_30 = ifelse(NAs == N, NA, ifelse(CrCl30 == 0, 0, 1)),
         kidney_flag_60 = ifelse(NAs == N, NA, ifelse(CrCl60 == 0, 0, 1))) %>%
  filter(!is.na(kidney_flag_15)) # remove NAs
# kidney_flag is individual level indicator by time frame; if patient had at least one CrCl below range by time frame, that patient would have 1 in this variable.
lab_data_before_sum$siteid <- ifelse(lab_data_before_sum$studyid>"marked-out id", 2, 1)
table(lab_data_before_sum$CrCl15, useNA = "ifany") # check
table(lab_data_before_sum$CrCl30, useNA = "ifany") # check
table(lab_data_before_sum$CrCl60, useNA = "ifany") # check
table(lab_data_before_sum$kidney_flag_15, useNA = "ifany")
table(lab_data_before_sum$kidney_flag_30, useNA = "ifany")
table(lab_data_before_sum$kidney_flag_60, useNA = "ifany")

chart_before_sum <- chart_before %>%
  mutate(kidney_d = ifelse(toxicityname == "Renal toxicity", 1, 0)) %>%
  group_by(studyid) %>%
  summarise(kidney_n = sum(kidney_d))
chart_before_sum <- chart_before_sum %>% mutate(kidney_tox = ifelse(kidney_n == 0, 0, 1))
table(chart_before_sum$kidney_n, useNA = "ifany")  

## Merge ICD and Lab
lab_data_before_sum$missing_lab <- rep(0)
icd_data_before_sum$missing_icd <- rep(0)
lab_id_1 <- lab_data_before_sum$studyid[lab_data_before_sum$time_frame==1]
lab_id_2 <- lab_data_before_sum$studyid[lab_data_before_sum$time_frame==2]
lab_id_3 <- lab_data_before_sum$studyid[lab_data_before_sum$time_frame==3]
lab_id_4 <- lab_data_before_sum$studyid[lab_data_before_sum$time_frame==4]
icd_id_1 <- icd_data_before_sum$studyid[icd_data_before_sum$time_frame==1]
icd_id_2 <- icd_data_before_sum$studyid[icd_data_before_sum$time_frame==2]
icd_id_3 <- icd_data_before_sum$studyid[icd_data_before_sum$time_frame==3]
icd_id_4 <- icd_data_before_sum$studyid[icd_data_before_sum$time_frame==4]
sum(!(lab_id_1 %in% icd_id_1)) 
sum(!(lab_id_2 %in% icd_id_2)) 
sum(!(lab_id_3 %in% icd_id_3)) 
sum(!(lab_id_4 %in% icd_id_4)) 
sum(!(icd_id_1 %in% lab_id_1)) 
sum(!(icd_id_2 %in% lab_id_2)) 
sum(!(icd_id_3 %in% lab_id_3)) 
sum(!(icd_id_4 %in% lab_id_4)) 

a <- data.frame(
  studyid = c(icd_id_1[!(icd_id_1 %in% lab_id_1)],
              icd_id_2[!(icd_id_2 %in% lab_id_2)],
              icd_id_3[!(icd_id_3 %in% lab_id_3)],
              icd_id_4[!(icd_id_4 %in% lab_id_4)]),
  time_frame = c(rep(1, times = 709),
                 rep(2, times = 1587),
                 rep(3, times = 168),
                 rep(4, times = 152)),
  missing_lab = rep(1), 
  kidney_flag_15 = rep(0),
  kidney_flag_30 = rep(0),
  kidney_flag_60 = rep(0)
)

lab_data_before_sum <- rbind(lab_data_before_sum, a)
lab_data_before_sum$siteid <- ifelse(lab_data_before_sum$studyid>"marked-out id", 2, 1)
anyNA(lab_data_before_sum$studyid)
anyNA(lab_data_before_sum$time_frame)
anyNA(lab_data_before_sum$kidney_flag_15)
anyNA(lab_data_before_sum$CrCl15)

b <- data.frame(
  studyid = c(lab_id_1[!(lab_id_1 %in% icd_id_1)],
              lab_id_2[!(lab_id_2 %in% icd_id_2)],
              lab_id_3[!(lab_id_3 %in% icd_id_3)],
              lab_id_4[!(lab_id_4 %in% icd_id_4)]),
  time_frame = c(rep(1, times = 3045),
                 rep(2, times = 3485),
                 rep(3, times = 7507),
                 rep(4, times = 5631)),
  missing_icd = rep(1),
  kidney_flag = rep(0)
)
icd_data_before_sum <- rbind(icd_data_before_sum, b)
icd_data_before_sum$siteid <- ifelse(icd_data_before_sum$studyid>"marked-out id", 2, 1)
anyNA(icd_data_before_sum$studyid)
anyNA(icd_data_before_sum$kidney_flag)
anyNA(icd_data_before_sum$kidney_n)

icd_data_before_sum <- icd_data_before_sum %>% arrange(studyid, time_frame)
lab_data_before_sum <- lab_data_before_sum %>% arrange(studyid, time_frame)
length(unique(lab_data_before_sum$studyid)) 
length(unique(icd_data_before_sum$studyid)) 
identical(icd_data_before_sum$studyid, lab_data_before_sum$studyid) # T
identical(icd_data_before_sum$time_frame, lab_data_before_sum$time_frame) # T

# merge icd and lab
before_sum <- icd_data_before_sum %>%
  inner_join(lab_data_before_sum[, names(lab_data_before_sum)[c(1:10, 12)]], by = join_by(studyid, time_frame))
adju_cohort <- obcd_group %>% filter(group == "Surgery before chemo") 
adju_cohort <- coredata(adju_cohort)[rep(seq(nrow(adju_cohort)), each = 4), ]
adju_cohort$time_frame <- rep(c(1:4), 11859)
before_sum <- before_sum[, names(before_sum)[c(1:4, 6:15)]] %>%
  right_join(adju_cohort, by = join_by(studyid, time_frame))
glimpse(before_sum)
table(before_sum$kidney_flag, useNA = "ifany")
before_sum$kidney_flag[is.na(before_sum$kidney_flag)] <- 0
table(before_sum$kidney_flag_15, useNA = "ifany")
before_sum$kidney_flag_15[is.na(before_sum$kidney_flag_15)] <- 0
table(before_sum$kidney_flag_30, useNA = "ifany")
before_sum$kidney_flag_30[is.na(before_sum$kidney_flag_30)] <- 0
table(before_sum$kidney_flag_60, useNA = "ifany")
before_sum$kidney_flag_60[is.na(before_sum$kidney_flag_60)] <- 0
table(before_sum$missing_icd, useNA = "ifany")
before_sum$missing_icd[is.na(before_sum$missing_icd)] <- 1
table(before_sum$missing_lab, useNA = "ifany")
before_sum$missing_lab[is.na(before_sum$missing_lab)] <- 1

# create renal_tox_indicator
before_sum <- before_sum %>%
  mutate(renal_tox_indicator = ifelse(kidney_flag == 1 & kidney_flag_60 == 1, "Both",
                                      ifelse(kidney_flag == 1 & kidney_flag_60 == 0, "Only ICD",
                                             ifelse(kidney_flag == 0 & kidney_flag_60 == 1, "Only CrCl60", "Neither"))))

# subset before_sum and export
adjuvant_renal <- before_sum[, c("studyid", "time_frame", "kidney_flag", "kidney_flag_15", "kidney_flag_30", "kidney_flag_60",
                                 "missing_icd", "missing_lab", "dx_dt", "siteid", "first_chemo_dt", "last_chemo_dt",
                                 "surgery_d", "surgery_dt", "group", "renal_tox_indicator")]
write.csv(adjuvant_renal, file = paste0(datapath, "adjuvant_renal.csv"),
          row.names = F)

## Merge IOCD, Lab, and Chart
length(unique(chart_before_sum$studyid)) 
length(unique(before_sum$studyid[before_sum$time_frame==4 & before_sum$siteid==2])) 
tox_sum <- chart_before_sum[, c("studyid", "kidney_tox")] %>%
  right_join(before_sum[before_sum$time_frame == 4 & before_sum$siteid==2,], by = join_by(studyid))
glimpse(tox_sum) 
anyNA(tox_sum$kidney_tox)
tox_sum$kidney_tox[is.na(tox_sum$kidney_tox)] <- 0


# export data
# saveRDS(icd_data_before_sum, file = paste0(datapath, "icd_before_sum.rds"))
# saveRDS(lab_data_before_sum, file = paste0(datapath, "lab_before_sum.rds"))
# saveRDS(before_sum, file = paste0(datapath, "before_sum.rds"))
# saveRDS(tox_sum, file = paste0(datapath, "tox_sum.rds"))
# saveRDS(adju_cohort, file = paste0(datapath, "adju_cohort.rds"))
# import data
# icd_data_before_sum <- readRDS(file = paste0(datapath, "icd_before_sum.rds"))
# lab_data_before_sum <- readRDS(file = paste0(datapath, "lab_before_sum.rds"))
# before_sum <- readRDS(file = paste0(datapath, "before_sum.rds"))
# tox_sum <- readRDS(file = paste0(datapath, "tox_sum.rds"))


## Hepatic ##
## Load Dataset
icd_data <- readRDS(file = paste0(datapath, "icd.rds"))
lab_data <- readRDS(file = paste0(datapath, "lab.rds"))
chart_data <- readRDS(file = paste0(datapath, "chart.rds"))
obcd_group <- readRDS(file = paste0(datapath, "obcd.rds"))

## Check and Clean
length(unique(icd_data$studyid)) 
length(unique(icd_data$studyid[icd_data$any_chemo==1])) 
anyNA(icd_data$co_tox_dx_code_date) # F
table(icd_data$liver_flag_icd, useNA = "ifany") 
length(unique(lab_data$studyid)) 
length(unique(lab_data$studyid[lab_data$any_chemo==1])) 
anyNA(lab_data$co_tox_lab_date) #F
anyNA(lab_data$co_tox_lab_name) #F
sum(is.na(lab_data$co_tox_lab_result[lab_data$co_tox_lab_name == 3 ])) # 1 NAs for bili dir
sum(is.na(lab_data$co_tox_lab_result[lab_data$co_tox_lab_name == 5 ])) # 1 NAs for bili tot
sum(is.na(lab_data$co_tox_lab_result[lab_data$co_tox_lab_name == 6 ])) # 5 NAs for AST
sum(is.na(lab_data$co_tox_lab_result[lab_data$co_tox_lab_name == 7 ])) # 0 NAs for ALT
sum(is.na(lab_data$co_tox_lab_result[lab_data$co_tox_lab_name == 8 ])) # 1 NAs for ALP

table(lab_data$co_tox_lab_unit[lab_data$co_tox_lab_name==3], useNA = "ifany") # all mg/dL
table(lab_data$co_tox_lab_unit[lab_data$co_tox_lab_name==4], useNA = "ifany") # all mg/dL
table(lab_data$co_tox_lab_unit[lab_data$co_tox_lab_name==5], useNA = "ifany") # 6 NA, 1 unknown, other mg/dL
table(lab_data$co_tox_lab_unit[lab_data$co_tox_lab_name==6], useNA = "ifany") # 21 NA, 1 unknown, 103 IU/L, other U/L
table(lab_data$co_tox_lab_unit[lab_data$co_tox_lab_name==7], useNA = "ifany") # 23 NA, 1 unknown, 104 IU/L, other U/L
table(lab_data$co_tox_lab_unit[lab_data$co_tox_lab_name==8], useNA = "ifany") # 3 NA, 31 IU/L, other U/L
# for IU and U: they are assumed to be same

# bili
summary(lab_data$co_tox_lab_result[lab_data$co_tox_lab_name==5])
as.data.frame(lab_data[lab_data$co_tox_lab_name==5 & lab_data$co_tox_lab_unit %in% c(NA, 11), c("studyid", "co_tox_lab_date", "co_tox_lab_result", "co_tox_lab_unit")]) # all should be mg/dL
as.data.frame(lab_data[lab_data$co_tox_lab_name==5 & lab_data$co_tox_lab_result<0, c("studyid", "co_tox_lab_date", "co_tox_lab_result", "co_tox_lab_unit")]) # result should be 0.1
lab_data$co_tox_lab_result[lab_data$co_tox_lab_result == -0.1] <- 0.1
# AST
summary(lab_data$co_tox_lab_result[lab_data$co_tox_lab_name==6])
as.data.frame(lab_data[lab_data$co_tox_lab_name==6 & lab_data$co_tox_lab_unit %in% c(NA, 11), c("studyid", "co_tox_lab_date", "co_tox_lab_result", "co_tox_lab_unit")]) # all should be U/L
# ALT
summary(lab_data$co_tox_lab_result[lab_data$co_tox_lab_name==7])
as.data.frame(lab_data[lab_data$co_tox_lab_name==7 & lab_data$co_tox_lab_unit %in% c(NA, 11), c("studyid", "co_tox_lab_date", "co_tox_lab_result", "co_tox_lab_unit")]) # all should be U/L
# ALP
summary(lab_data$co_tox_lab_result[lab_data$co_tox_lab_name==8])
as.data.frame(lab_data[lab_data$co_tox_lab_name==8 & lab_data$co_tox_lab_unit %in% c(NA, 11), c("studyid", "co_tox_lab_date", "co_tox_lab_result", "co_tox_lab_unit")]) # all should be U/L
# for bili: all are assumed to be mg/dL
# for AST, ALT, ALP: all are assumed to be U/L

# Normal High
sum(is.na(lab_data$co_tox_lab_normal_high[lab_data$co_tox_lab_name == 6 ])) # 2568 NAs for AST
sum(is.na(lab_data$co_tox_lab_normal_high[lab_data$co_tox_lab_name == 7 ])) # 2547 NAs for ALT
sum(is.na(lab_data$co_tox_lab_normal_high[lab_data$co_tox_lab_name == 8 ])) # 2347 NAs for ALP
sum(is.na(lab_data$co_tox_lab_normal_high[lab_data$co_tox_lab_name == 6 & lab_data$siteid == 2])) # 2567 NAs for AST
sum(is.na(lab_data$co_tox_lab_normal_high[lab_data$co_tox_lab_name == 7 & lab_data$siteid == 2])) # 2545 NAs for ALT
sum(is.na(lab_data$co_tox_lab_normal_high[lab_data$co_tox_lab_name == 8 & lab_data$siteid == 2])) # 2347 NAs for ALP
nrow(lab_data[lab_data$co_tox_lab_name == 6,]) 
nrow(lab_data[lab_data$co_tox_lab_name == 7,]) 
nrow(lab_data[lab_data$co_tox_lab_name == 8,]) 
length(unique(lab_data$studyid[is.na(lab_data$co_tox_lab_normal_high) & lab_data$co_tox_lab_name == 6]))
# 1122 unique patients had NA for AST
length(unique(lab_data$studyid[is.na(lab_data$co_tox_lab_normal_high) & lab_data$co_tox_lab_name == 7]))
# 1104 unique patients had NA for ALT
length(unique(lab_data$studyid[is.na(lab_data$co_tox_lab_normal_high) & lab_data$co_tox_lab_name == 8]))
# 1024 unique patients had NA for ALP

# AST: Normal High
summary(lab_data$co_tox_lab_normal_high[lab_data$co_tox_lab_name == 6 & lab_data$co_tox_lab_unit == 6 & !is.na(lab_data$co_tox_lab_unit)])
# Min: 35, 1st: 36, 3rd: 40, Max: 59, Median: 40, 2546 NAs
summary(lab_data$co_tox_lab_normal_high[lab_data$co_tox_lab_name == 6 & lab_data$co_tox_lab_unit %in% c(NA, 11)])
# all 37, 21 NAs
summary(lab_data$co_tox_lab_normal_high[lab_data$co_tox_lab_name == 6 & lab_data$co_tox_lab_unit == 9 & !is.na(lab_data$co_tox_lab_unit)])
# all 40, 1 NAs


# ALT: Normal High
summary(lab_data$co_tox_lab_normal_high[lab_data$co_tox_lab_name == 7 & lab_data$co_tox_lab_unit == 6 & !is.na(lab_data$co_tox_lab_unit)])
# Min: 29, 1st: 36, 3rd: 41, Max: 66, Median: 41, 2523 NAs
summary(lab_data$co_tox_lab_normal_high[lab_data$co_tox_lab_name == 7 & lab_data$co_tox_lab_unit %in% c(NA, 11)])
# all 78, 23 NAs
summary(lab_data$co_tox_lab_normal_high[lab_data$co_tox_lab_name == 7 & lab_data$co_tox_lab_unit == 9 & !is.na(lab_data$co_tox_lab_unit)])
# all 50, 1 NAs

# ALP: Normal High
summary(lab_data$co_tox_lab_normal_high[lab_data$co_tox_lab_name == 8 & lab_data$co_tox_lab_unit == 6 & !is.na(lab_data$co_tox_lab_unit)])
# Min: 110, 1st: 117, 3rd: 125, Max: 153, Median: 117, 2344 NAs
summary(lab_data$co_tox_lab_normal_high[lab_data$co_tox_lab_name == 8 & lab_data$co_tox_lab_unit %in% c(NA, 11)])
# 3 NAs
summary(lab_data$co_tox_lab_normal_high[lab_data$co_tox_lab_name == 8 & lab_data$co_tox_lab_unit == 9 & !is.na(lab_data$co_tox_lab_unit)])
# Min: 110, 1st: 110, Median: 125, 3rd: 130, Max: 140


# check chart data
a <- chart_data[is.na(chart_data$toxicitydate), c("studyid", "toxicityname", "toxicitydate", "dx_dt", "surgery_d", "surgery_dt", "first_chemo_dt")] %>% print(n=Inf, width = Inf)
b <- chart_data[!is.na(chart_data$toxicitydate) & (chart_data$toxicitydate < chart_data$first_chemo_dt), c("studyid", "toxicityid", "toxicityname", "toxicitydate", "dx_dt", "surgery_d", "surgery_dt", "first_chemo_dt", "last_chemo_dt")] %>% print(n=Inf, width = Inf)
c <- chart_data[!is.na(chart_data$toxicitydate) & (chart_data$toxicitydate == chart_data$first_chemo_dt), c("studyid", "toxicityname", "toxicitydate", "dx_dt", "surgery_d", "surgery_dt", "first_chemo_dt")] %>% print(n=Inf, width = Inf)


## Adding Time Frame
# time frame for surgery before chemo
# 1: <= date of cancer diagnosis --> pre-existing conditions (1 year)
# 2: <= date of surgery
# 3: <= date of first chemo
# 4: <= end of chemo
icd_data_before <- icd_data %>%
  filter(group == "Surgery before chemo") %>%
  mutate(time_frame = ifelse(co_tox_dx_code_date <= dx_dt, 1,
                             ifelse(co_tox_dx_code_date <= surgery_dt, 2,
                                    ifelse(co_tox_dx_code_date <= first_chemo_dt, 3,
                                           ifelse(co_tox_dx_code_date <= last_chemo_dt, 4, NA))))) %>%
  filter(!is.na(time_frame)) # remove NAs in time frame

table(icd_data_before$time_frame, useNA = "ifany")
table(icd_data_before$liver_flag_icd, useNA = "ifany") 

lab_data_before <- lab_data %>%
  filter(group == "Surgery before chemo") %>%
  mutate(time_frame = ifelse(co_tox_lab_date <= dx_dt, 1,
                             ifelse(co_tox_lab_date <= surgery_dt, 2,
                                    ifelse(co_tox_lab_date <= first_chemo_dt, 3,
                                           ifelse(co_tox_lab_date <= last_chemo_dt, 4, NA))))) %>%
  filter(!is.na(time_frame)) # remove NAs in time frame

table(lab_data_before$time_frame, useNA = "ifany")
table(lab_data_before$co_tox_lab_name, useNA = "ifany") 
# 3=bili dir; 4=bili indir; 5=bili total; 6=AST; 7=ALT; 8=ALP

chart_before <- chart_data %>% filter(surgery_d == 1) %>%
  filter(!is.na(surgery_dt)) %>%
  filter(surgery_dt < first_chemo_dt) %>%
  filter(toxicitydate > first_chemo_dt & toxicitydate <= last_chemo_dt) # limit to last timeframe


## Check 
# check normal high in lab_data_before
sum(is.na(lab_data_before$co_tox_lab_normal_high[lab_data_before$co_tox_lab_name == 6 ])) # 428 NAs for AST
sum(is.na(lab_data_before$co_tox_lab_normal_high[lab_data_before$co_tox_lab_name == 7 ])) # 428 NAs for ALT
sum(is.na(lab_data_before$co_tox_lab_normal_high[lab_data_before$co_tox_lab_name == 8 ])) # 423 NAs for ALP
sum(is.na(lab_data_before$co_tox_lab_normal_high[lab_data_before$co_tox_lab_name == 6 & lab_data_before$siteid == 2])) # 428 NAs for AST
sum(is.na(lab_data_before$co_tox_lab_normal_high[lab_data_before$co_tox_lab_name == 7 & lab_data_before$siteid == 2])) # 428 NAs for AST
sum(is.na(lab_data_before$co_tox_lab_normal_high[lab_data_before$co_tox_lab_name == 8 & lab_data_before$siteid == 2])) # 423 NAs for AST
nrow(lab_data_before[lab_data_before$co_tox_lab_name == 6,]) 
nrow(lab_data_before[lab_data_before$co_tox_lab_name == 7,]) 
nrow(lab_data_before[lab_data_before$co_tox_lab_name == 8,]) 
length(unique(lab_data_before$studyid[is.na(lab_data_before$co_tox_lab_normal_high) & lab_data_before$co_tox_lab_name == 6]))
# 95 unique patients had NA for AST
length(unique(lab_data_before$studyid[is.na(lab_data_before$co_tox_lab_normal_high) & lab_data_before$co_tox_lab_name == 7]))
# 95 unique patients had NA for ALT
length(unique(lab_data_before$studyid[is.na(lab_data_before$co_tox_lab_normal_high) & lab_data_before$co_tox_lab_name == 8]))
# 95 unique patients had NA for ALP
# all NAs in KPWA
lab_data_before %>% filter(co_tox_lab_name %in% 6:8) %>%
  group_by(time_frame, co_tox_lab_name) %>%
  summarise(NAs = sum(is.na(co_tox_lab_normal_high)), 
            NAs_patient = length(unique(studyid[is.na(co_tox_lab_normal_high)])),
            N = n(),
            N_patient = length(unique(studyid)))
# most of NAs in timeframe=4

lab_kpwa_sum <- lab_data_before %>% filter(co_tox_lab_name %in% 6:8 & siteid == 2) %>%
  group_by(time_frame, co_tox_lab_name) %>%
  summarise(NAs = sum(is.na(co_tox_lab_normal_high)), 
            N = n(),
            NAs_patient = length(unique(studyid[is.na(co_tox_lab_normal_high)])),
            N_patient = length(unique(studyid))) %>%
  mutate(Percent_N = NAs/N*100, Percent_patient = NAs_patient/N_patient*100)

# summary of normal high in KPWA
summary(lab_data_before$co_tox_lab_normal_high[lab_data_before$co_tox_lab_name == 6 & lab_data_before$siteid==2])
summary(lab_data_before$co_tox_lab_normal_high[lab_data_before$co_tox_lab_name == 7 & lab_data_before$siteid==2])
summary(lab_data_before$co_tox_lab_normal_high[lab_data_before$co_tox_lab_name == 8 & lab_data_before$siteid==2])

# check missing in AST, ALT, ALP for KPWA
sum(is.na(lab_data_before$co_tox_lab_result[lab_data_before$co_tox_lab_name == 6 & lab_data_before$siteid==2])) # 0 NAs for AST
sum(is.na(lab_data_before$co_tox_lab_result[lab_data_before$co_tox_lab_name == 7 & lab_data_before$siteid==2])) # 0 NAs for ALT
sum(is.na(lab_data_before$co_tox_lab_result[lab_data_before$co_tox_lab_name == 8 & lab_data_before$siteid==2])) # 0 NAs for ALP

# patients with any missingness in normal high
a <- unique(lab_data_before$studyid[is.na(lab_data_before$co_tox_lab_normal_high) & lab_data_before$co_tox_lab_name %in% 6:8])
lab_missing <- lab_data_before %>% filter(studyid %in% a & co_tox_lab_name %in% 6:8) %>%
  group_by(studyid, time_frame, co_tox_lab_name) %>%
  mutate(count = row_number())
length(unique(lab_missing$studyid))
lab_missing_sum <- lab_missing %>% group_by(studyid, time_frame, co_tox_lab_name) %>%
  summarise(N= n(), NAs = sum(is.na(co_tox_lab_normal_high))) %>%
  mutate(missing_ind = ifelse(N == NAs, "full", ifelse(NAs == 0, "none", "partial")))
b <- lab_missing_sum %>% group_by(time_frame, co_tox_lab_name, missing_ind) %>%
  summarize(N_missing = n()) %>% 
  filter(missing_ind != "none")
b$N_patient <- c(258, 258, 196, 196, 167, 167, 309, 291, 298,
                 666, 666, 583, 583, 657, 657, 715, 715,
                 628, 628, 700, 700)
b$percent_missing <- b$N_missing/b$N_patient*100
b %>% print(n = Inf)

summary(lab_missing$co_tox_lab_result[lab_missing$co_tox_lab_name == 6 & is.na(lab_missing$co_tox_lab_normal_high)])
summary(lab_missing$co_tox_lab_result[lab_missing$co_tox_lab_name == 7 & is.na(lab_missing$co_tox_lab_normal_high)])
summary(lab_missing$co_tox_lab_result[lab_missing$co_tox_lab_name == 8 & is.na(lab_missing$co_tox_lab_normal_high)])
lab_missing[lab_missing$co_tox_lab_name == 6 & lab_missing$co_tox_lab_result >=70 & is.na(lab_missing$co_tox_lab_normal_high), c("studyid", "time_frame", "co_tox_lab_result", "co_tox_lab_normal_high")]
lab_missing[lab_missing$co_tox_lab_name == 7 & lab_missing$co_tox_lab_result >=80 & is.na(lab_missing$co_tox_lab_normal_high), c("studyid", "time_frame", "co_tox_lab_result", "co_tox_lab_normal_high")] %>% print(n = Inf)
lab_missing[lab_missing$co_tox_lab_name == 8 & lab_missing$co_tox_lab_result >=275 & is.na(lab_missing$co_tox_lab_normal_high), c("studyid", "time_frame", "co_tox_lab_result", "co_tox_lab_normal_high")]

c <- lab_missing_sum[lab_missing_sum$studyid %in% c(201073, 205197, 206272, 200561, 202505, 203743, 204389, 204668, 206174), ] %>% print(n = Inf)
write.csv(c, file = paste0(outpath, "lab_missing_ind.csv"))
lab_missing[lab_missing$studyid == 201073, c("time_frame", "co_tox_lab_name", "co_tox_lab_result", "co_tox_lab_normal_high", "co_tox_lab_date")] %>% print(n = Inf, width = Inf)
lab_missing[lab_missing$studyid %in% c(201073, 205197, 206272, 200561, 202505, 203743, 204389, 204668, 206174) & !is.na(lab_missing$co_tox_lab_normal_high), c("studyid", "time_frame", "co_tox_lab_name", "co_tox_lab_result", "co_tox_lab_normal_high", "co_tox_lab_date")] %>% print(n = Inf, width = Inf)
lab_missing[lab_missing$studyid == 200561 & !is.na(lab_missing$co_tox_lab_normal_high) & lab_missing$co_tox_lab_name %in% 6:8 & lab_missing$time_frame == 4, c("studyid", "time_frame", "co_tox_lab_name", "co_tox_lab_result", "co_tox_lab_normal_high", "co_tox_lab_date")] %>% print(n = Inf, width = Inf)
lab_missing[lab_missing$studyid == 201073 & !is.na(lab_missing$co_tox_lab_normal_high) & lab_missing$co_tox_lab_name %in% 6:8 & lab_missing$time_frame == 4, c("studyid", "time_frame", "co_tox_lab_name", "co_tox_lab_result", "co_tox_lab_normal_high", "co_tox_lab_date")] %>% print(n = Inf, width = Inf)
lab_missing[lab_missing$studyid == 202505 & !is.na(lab_missing$co_tox_lab_normal_high) & lab_missing$co_tox_lab_name %in% 6:8 & lab_missing$time_frame == 4, c("studyid", "time_frame", "co_tox_lab_name", "co_tox_lab_result", "co_tox_lab_normal_high", "co_tox_lab_date")] %>% print(n = Inf, width = Inf)
lab_missing[lab_missing$studyid == 203743 & !is.na(lab_missing$co_tox_lab_normal_high) & lab_missing$co_tox_lab_name %in% 6:8 & lab_missing$time_frame == 4, c("studyid", "time_frame", "co_tox_lab_name", "co_tox_lab_result", "co_tox_lab_normal_high", "co_tox_lab_date")] %>% print(n = Inf, width = Inf)
lab_missing[lab_missing$studyid == 204668 & !is.na(lab_missing$co_tox_lab_normal_high) & lab_missing$co_tox_lab_name %in% 6:8 & lab_missing$time_frame == 4, c("studyid", "time_frame", "co_tox_lab_name", "co_tox_lab_result", "co_tox_lab_normal_high", "co_tox_lab_date")] %>% print(n = Inf, width = Inf)
lab_missing[lab_missing$studyid == 206174 & !is.na(lab_missing$co_tox_lab_normal_high) & lab_missing$co_tox_lab_name %in% 6:8 & lab_missing$time_frame == 4, c("studyid", "time_frame", "co_tox_lab_name", "co_tox_lab_result", "co_tox_lab_normal_high", "co_tox_lab_date")] %>% print(n = Inf, width = Inf)


## Imputate Normal High 
lab_data_before <- lab_data_before %>%
  mutate(normal_high_missing = ifelse(is.na(co_tox_lab_normal_high), 1, 0)) # normal high missingness indicator

lab_data_before_liver <- lab_data_before %>%
  filter(co_tox_lab_name %in% 3:8) %>% # only hepatic indicators
  filter(!is.na(co_tox_lab_result)) # remove NAs in lab result

# patients with partial missingness of normal high in specific time frame
lab_missing_sum %>% filter(missing_ind == "partial") %>% print(n = Inf)
d <- unique(lab_missing_sum$studyid[lab_missing_sum$missing_ind == "partial"])
e <- lab_missing[lab_missing$studyid %in% d & !is.na(lab_missing$co_tox_lab_normal_high) & lab_missing$co_tox_lab_name %in% 6:8, c("studyid", "time_frame", "co_tox_lab_name", "co_tox_lab_result", "co_tox_lab_normal_high", "co_tox_lab_date")] %>% print(n = Inf, width = Inf)
f <- e %>% group_by(studyid, time_frame, co_tox_lab_name) %>%
  summarise(min = min(co_tox_lab_normal_high), 
            median = median(co_tox_lab_normal_high), 
            max = max(co_tox_lab_normal_high), n = n()) %>%
  mutate(Ind = ifelse(min == max, 1, 0))
lab_missing[lab_missing$studyid == "marked-out id" & lab_missing$co_tox_lab_name %in% 6:8 & lab_missing$time_frame == 4, c("studyid", "time_frame", "co_tox_lab_name", "co_tox_lab_result", "co_tox_lab_normal_high", "co_tox_lab_date")] %>% print(n = Inf, width = Inf)
## timeframe=4, lastest AST=35, latest ALT=40, latest ALP=125
lab_data_before_liver$co_tox_lab_normal_high[lab_data_before_liver$studyid == "marked-out id" & lab_data_before_liver$co_tox_lab_name == 6 & lab_data_before_liver$time_frame == 4 & is.na(lab_data_before_liver$co_tox_lab_normal_high)] <- 35
lab_data_before_liver$co_tox_lab_normal_high[lab_data_before_liver$studyid == "marked-out id" & lab_data_before_liver$co_tox_lab_name == 7 & lab_data_before_liver$time_frame == 4 & is.na(lab_data_before_liver$co_tox_lab_normal_high)] <- 40
lab_data_before_liver$co_tox_lab_normal_high[lab_data_before_liver$studyid == "marked-out id" & lab_data_before_liver$co_tox_lab_name == 8 & lab_data_before_liver$time_frame == 4 & is.na(lab_data_before_liver$co_tox_lab_normal_high)] <- 125
lab_data_before_liver[lab_data_before_liver$studyid=="marked-out id" & lab_data_before_liver$time_frame == 4 & lab_data_before_liver$co_tox_lab_name %in% 6:8, "co_tox_lab_normal_high"] %>% print(n = Inf)
## other patients had same normal high values
g <- d[d!="marked-out id"]
for (i in 1:length(g)){
  for (j in 6:8) {
    for (k in 1:4){
      t <- length(lab_data_before_liver$co_tox_lab_normal_high[lab_data_before_liver$studyid == g[i] & lab_data_before_liver$co_tox_lab_name == j & lab_data_before_liver$time_frame == k & is.na(lab_data_before_liver$co_tox_lab_normal_high)])
      s <- length(f$min[f$studyid == g[i] & f$co_tox_lab_name == j & f$time_frame == k])
      if (t != 0 & s != 0) {
        lab_data_before_liver$co_tox_lab_normal_high[lab_data_before_liver$studyid == g[i] & lab_data_before_liver$co_tox_lab_name == j & lab_data_before_liver$time_frame == k & is.na(lab_data_before_liver$co_tox_lab_normal_high)] <- f$min[f$studyid == g[i] & f$co_tox_lab_name == j & f$time_frame == k]
      }
    }
  }
}
lab_data_before_liver[lab_data_before_liver$studyid=="marked-out id" & lab_data_before_liver$time_frame == 1 & lab_data_before_liver$co_tox_lab_name %in% 6:8, "co_tox_lab_normal_high"] %>% print(n = Inf)
lab_data_before_liver[lab_data_before_liver$studyid=="marked-out id" & lab_data_before_liver$time_frame == 3 & lab_data_before_liver$co_tox_lab_name %in% 6:8, "co_tox_lab_normal_high"] %>% print(n = Inf)
# min normal high for AST in KPWA: 35
# min normal high for ALT in KPWA: 40
# min normal high for ALP in KPWA: 110
lab_data_before_liver$co_tox_lab_normal_high[lab_data_before_liver$co_tox_lab_name == 6 & is.na(lab_data_before_liver$co_tox_lab_normal_high)] <- 35
lab_data_before_liver$co_tox_lab_normal_high[lab_data_before_liver$co_tox_lab_name == 7 & is.na(lab_data_before_liver$co_tox_lab_normal_high)] <- 40
lab_data_before_liver$co_tox_lab_normal_high[lab_data_before_liver$co_tox_lab_name == 8 & is.na(lab_data_before_liver$co_tox_lab_normal_high)] <- 110

# export data
# saveRDS(lab_data_before_liver, file = paste0(datapath, "lab_before_liver.rds"))


## Flag Liver
anyNA(lab_data_before_liver$co_tox_lab_result) # F
anyNA(lab_data_before_liver$co_tox_lab_normal_high[lab_data_before_liver$co_tox_lab_name %in% 6:8]) # F

lab_data_before_liver <- lab_data_before_liver %>%
  mutate(
    liver_flag_lab = ifelse(co_tox_lab_name == 5 & co_tox_lab_result > 1.5, 1,
                            ifelse(co_tox_lab_name == 6 & co_tox_lab_result > co_tox_lab_normal_high*3, 1,
                                   ifelse(co_tox_lab_name == 7 & co_tox_lab_result > co_tox_lab_normal_high*3, 1,
                                          ifelse(co_tox_lab_name == 8 & co_tox_lab_result > co_tox_lab_normal_high*2.5, 1, 0)))),
    liver_flag_lab1 = ifelse(co_tox_lab_name == 5 & co_tox_lab_result > 3, 1,
                             ifelse(co_tox_lab_name == 6 & co_tox_lab_result > co_tox_lab_normal_high*5, 1,
                                    ifelse(co_tox_lab_name == 7 & co_tox_lab_result > co_tox_lab_normal_high*5, 1,
                                           ifelse(co_tox_lab_name == 8 & co_tox_lab_result > co_tox_lab_normal_high*5, 1, 0)))),
    liver_flag_lab2 = ifelse(co_tox_lab_name == 5 & co_tox_lab_result >= 1.2, 1,
                             ifelse(co_tox_lab_name == 6 & co_tox_lab_result >= co_tox_lab_normal_high*2, 1,
                                    ifelse(co_tox_lab_name == 7 & co_tox_lab_result >= co_tox_lab_normal_high*2, 1,
                                           ifelse(co_tox_lab_name == 8 & co_tox_lab_result > co_tox_lab_normal_high*2.5, 1, 0))))
    
  )

## Summarize Hepatic Toxicity
icd_data_before_sum_liver <- icd_data_before %>% 
  select(studyid, time_frame, liver_flag_icd) %>%
  group_by(studyid, time_frame) %>%
  summarise(liver_icd_n = sum(liver_flag_icd))
icd_data_before_sum_liver <- icd_data_before_sum_liver %>% mutate(liver_flag_icd = ifelse(liver_icd_n == 0, 0, 1))
icd_data_before_sum_liver$siteid <- ifelse(icd_data_before_sum_liver$studyid>"marked-out id", 2, 1)
table(icd_data_before_sum_liver$liver_icd_n, useNA = "ifany") # check
table(icd_data_before_sum_liver$liver_flag_icd, useNA = "ifany") # check
icd_data_before_sum_liver[icd_data_before_sum_liver$liver_icd_n >= 1,] # check

lab_data_before_sum_liver <- lab_data_before_liver %>%
  select(studyid, time_frame, liver_flag_lab, liver_flag_lab1, liver_flag_lab2) %>%
  group_by(studyid, time_frame) %>%
  summarise(
    liver_lab_n = sum(liver_flag_lab),
    liver_lab1_n = sum(liver_flag_lab1),
    liver_lab2_n = sum(liver_flag_lab2)
  ) %>%
  mutate(
    liver_flag_lab = ifelse(liver_lab_n == 0, 0, 1),
    liver_flag_lab1 = ifelse(liver_lab1_n == 0, 0, 1),
    liver_flag_lab2 = ifelse(liver_lab2_n == 0, 0, 1)
  )
lab_data_before_sum_liver$siteid <- ifelse(lab_data_before_sum_liver$studyid>"marked-out id", 2, 1)
table(lab_data_before_sum_liver$liver_lab_n, useNA = "ifany") # check
table(lab_data_before_sum_liver$liver_lab1_n, useNA = "ifany") # check
table(lab_data_before_sum_liver$liver_lab2_n, useNA = "ifany") # check
table(lab_data_before_sum_liver$liver_flag_lab, useNA = "ifany")
table(lab_data_before_sum_liver$liver_flag_lab1, useNA = "ifany")
table(lab_data_before_sum_liver$liver_flag_lab2, useNA = "ifany")


chart_before_sum_liver <- chart_before %>%
  mutate(liver_d = ifelse(toxicityname == "Hepatic toxicity", 1, 0)) %>%
  group_by(studyid) %>%
  summarise(liver_n = sum(liver_d))
chart_before_sum_liver <- chart_before_sum_liver %>% mutate(liver_tox = ifelse(liver_n == 0, 0, 1))
table(chart_before_sum_liver$liver_tox, useNA = "ifany")
table(chart_before_sum_liver$liver_n, useNA = "ifany")


## Merge ICD and Lab
lab_data_before_sum_liver$missing_lab <- rep(0)
icd_data_before_sum_liver$missing_icd <- rep(0)
lab_id_1 <- lab_data_before_sum_liver$studyid[lab_data_before_sum_liver$time_frame==1]
lab_id_2 <- lab_data_before_sum_liver$studyid[lab_data_before_sum_liver$time_frame==2]
lab_id_3 <- lab_data_before_sum_liver$studyid[lab_data_before_sum_liver$time_frame==3]
lab_id_4 <- lab_data_before_sum_liver$studyid[lab_data_before_sum_liver$time_frame==4]
icd_id_1 <- icd_data_before_sum_liver$studyid[icd_data_before_sum_liver$time_frame==1]
icd_id_2 <- icd_data_before_sum_liver$studyid[icd_data_before_sum_liver$time_frame==2]
icd_id_3 <- icd_data_before_sum_liver$studyid[icd_data_before_sum_liver$time_frame==3]
icd_id_4 <- icd_data_before_sum_liver$studyid[icd_data_before_sum_liver$time_frame==4]
sum(!(lab_id_1 %in% icd_id_1)) 
sum(!(lab_id_2 %in% icd_id_2)) 
sum(!(lab_id_3 %in% icd_id_3)) 
sum(!(lab_id_4 %in% icd_id_4)) 
sum(!(icd_id_1 %in% lab_id_1)) 
sum(!(icd_id_2 %in% lab_id_2)) 
sum(!(icd_id_3 %in% lab_id_3)) 
sum(!(icd_id_4 %in% lab_id_4)) 

a <- data.frame(
  studyid = c(icd_id_1[!(icd_id_1 %in% lab_id_1)],
              icd_id_2[!(icd_id_2 %in% lab_id_2)],
              icd_id_3[!(icd_id_3 %in% lab_id_3)],
              icd_id_4[!(icd_id_4 %in% lab_id_4)]),
  time_frame = c(rep(1, times = 1468),
                 rep(2, times = 2389),
                 rep(3, times = 162),
                 rep(4, times = 179)),
  missing_lab = rep(1), 
  liver_flag_lab = rep(0),
  liver_flag_lab1 = rep(0),
  liver_flag_lab2 = rep(0)
)

lab_data_before_sum_liver <- rbind(lab_data_before_sum_liver, a)
lab_data_before_sum_liver$siteid <- ifelse(lab_data_before_sum_liver$studyid>"marked-out id", 2, 1)
anyNA(lab_data_before_sum_liver$studyid)
anyNA(lab_data_before_sum_liver$time_frame)
anyNA(lab_data_before_sum_liver$liver_flag_lab)

b <- data.frame(
  studyid = c(lab_id_1[!(lab_id_1 %in% icd_id_1)],
              lab_id_2[!(lab_id_2 %in% icd_id_2)],
              lab_id_3[!(lab_id_3 %in% icd_id_3)],
              lab_id_4[!(lab_id_4 %in% icd_id_4)]),
  time_frame = c(rep(1, times = 2614),
                 rep(2, times = 2839),
                 rep(3, times = 7623),
                 rep(4, times = 5641)),
  missing_icd = rep(1),
  liver_flag_icd = rep(0)
)
icd_data_before_sum_liver <- rbind(icd_data_before_sum_liver, b)
icd_data_before_sum_liver$siteid <- ifelse(icd_data_before_sum_liver$studyid>"marked-out id", 2, 1)
anyNA(icd_data_before_sum_liver$studyid)
anyNA(icd_data_before_sum_liver$liver_flag_icd)

icd_data_before_sum_liver <- icd_data_before_sum_liver %>% arrange(studyid, time_frame)
lab_data_before_sum_liver <- lab_data_before_sum_liver %>% arrange(studyid, time_frame)
length(unique(lab_data_before_sum_liver$studyid))
length(unique(icd_data_before_sum_liver$studyid))
identical(icd_data_before_sum_liver$studyid, lab_data_before_sum_liver$studyid) # T
identical(icd_data_before_sum_liver$time_frame, lab_data_before_sum_liver$time_frame) # T

# merge icd and lab
before_sum_liver <- icd_data_before_sum_liver %>%
  inner_join(lab_data_before_sum_liver[, names(lab_data_before_sum_liver)[c(1:8, 10)]], by = join_by(studyid, time_frame))
adju_cohort <- readRDS(file = paste0(datapath, "adju_cohort.rds"))
before_sum_liver <- before_sum_liver[, names(before_sum_liver)[c(1:4, 6:13)]] %>%
  right_join(adju_cohort, by = join_by(studyid, time_frame))
glimpse(before_sum_liver)
table(before_sum_liver$liver_flag_icd, useNA = "ifany")
before_sum_liver$liver_flag_icd[is.na(before_sum_liver$liver_flag_icd)] <- 0
table(before_sum_liver$liver_flag_lab, useNA = "ifany")
before_sum_liver$liver_flag_lab[is.na(before_sum_liver$liver_flag_lab)] <- 0
table(before_sum_liver$liver_flag_lab1, useNA = "ifany")
before_sum_liver$liver_flag_lab1[is.na(before_sum_liver$liver_flag_lab1)] <- 0
table(before_sum_liver$liver_flag_lab2, useNA = "ifany")
before_sum_liver$liver_flag_lab2[is.na(before_sum_liver$liver_flag_lab2)] <- 0
table(before_sum_liver$missing_icd, useNA = "ifany")
before_sum_liver$missing_icd[is.na(before_sum_liver$missing_icd)] <- 1
table(before_sum_liver$missing_lab, useNA = "ifany")
before_sum_liver$missing_lab[is.na(before_sum_liver$missing_lab)] <- 1

# create hep_tox_indicator
before_sum_liver <- before_sum_liver %>%
  mutate(hep_tox_indicator = ifelse(liver_flag_icd == 1 & liver_flag_lab == 1, "Both",
                                    ifelse(liver_flag_icd == 1 & liver_flag_lab == 0, "Only ICD",
                                           ifelse(liver_flag_icd == 0 & liver_flag_lab == 1, "Only Grade 2", "Neither"))))

# subset before_sum_liver and export
adjuvant_hepatic <- before_sum_liver[, c("studyid", "time_frame", "liver_flag_icd", "liver_flag_lab", "liver_flag_lab1", "liver_flag_lab2",
                                         "missing_icd", "missing_lab", "dx_dt", "siteid", "first_chemo_dt", "last_chemo_dt",
                                         "surgery_d", "surgery_dt", "group", "hep_tox_indicator")]
write.csv(adjuvant_hepatic, file = paste0(datapath, "adjuvant_hepatic.csv"),
          row.names = F)


## Merge IOCD, Lab, and Chart
length(unique(chart_before_sum_liver$studyid)) 
length(unique(before_sum_liver$studyid[before_sum_liver$time_frame==4 & before_sum_liver$siteid==2])) 
tox_sum_liver <- chart_before_sum_liver[, c("studyid", "liver_tox")] %>%
  right_join(before_sum_liver[before_sum_liver$time_frame == 4 & before_sum_liver$siteid==2,], by = join_by(studyid))
glimpse(tox_sum_liver) 
anyNA(tox_sum_liver$liver_tox)
tox_sum_liver$liver_tox[is.na(tox_sum_liver$liver_tox)] <- 0


# export data
# saveRDS(icd_data_before_sum_liver, file = paste0(datapath, "icd_before_sum_liver.rds"))
# saveRDS(lab_data_before_sum_liver, file = paste0(datapath, "lab_before_sum_liver.rds"))
# saveRDS(before_sum_liver, file = paste0(datapath, "before_sum_liver.rds"))
# saveRDS(tox_sum_liver, file = paste0(datapath, "tox_sum_liver.rds"))
# saveRDS(chart_before_sum_liver, file = paste0(datapath, "chart_before_sum_liver.rds"))
# icd_data_before_sum_liver <- readRDS(file = paste0(datapath, "icd_before_sum_liver.rds"))
# lab_data_before_sum_liver <- readRDS(file = paste0(datapath, "lab_before_sum_liver.rds"))
# before_sum_liver <- readRDS(file = paste0(datapath, "before_sum_liver.rds"))


## Anemia ##
## Load Dataset
icd_data <- readRDS(file = paste0(datapath, "icd.rds"))
lab_data <- readRDS(file = paste0(datapath, "lab.rds"))
chart_data <- readRDS(file = paste0(datapath, "chart.rds"))
obcd_group <- readRDS(file = paste0(datapath, "obcd.rds"))

## Check and Clean 
length(unique(icd_data$studyid)) 
length(unique(icd_data$studyid[icd_data$any_chemo==1])) 
anyNA(icd_data$co_tox_dx_code_date) # F
table(icd_data$anemia_flag_icd, useNA = "ifany") 
length(unique(lab_data$studyid)) 
length(unique(lab_data$studyid[lab_data$any_chemo==1])) 
anyNA(lab_data$co_tox_lab_date) #F
anyNA(lab_data$co_tox_lab_name) #F
sum(is.na(lab_data$co_tox_lab_result[lab_data$co_tox_lab_name == 17 ])) # no NA

table(lab_data$co_tox_lab_unit[lab_data$co_tox_lab_name==17], useNA = "ifany") # mg/dL, unknown, g/24
summary(lab_data$co_tox_lab_result[lab_data$co_tox_lab_name == 17 ]) 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.30   10.80   12.20   11.97   13.30   24.10
summary(lab_data$co_tox_lab_result[lab_data$co_tox_lab_name == 17 & lab_data$co_tox_lab_unit == 11]) 
summary(lab_data$co_tox_lab_result[lab_data$co_tox_lab_name == 17 & lab_data$co_tox_lab_unit == 12]) # similar to 11
lab_data$co_tox_lab_result[lab_data$co_tox_lab_name == 17 & lab_data$co_tox_lab_unit == 5] 
# if transform to g/dL, then all values need to be divided by 1000. does not make sense. 
lab_data$co_tox_lab_result[lab_data$co_tox_lab_name == 17 & is.na(lab_data$co_tox_lab_unit)]
# all units should be g/dL

# Normal Low
sum(is.na(lab_data$co_tox_lab_normal_low[lab_data$co_tox_lab_name == 17 ])) # 3917 NAs for anemia
sum(is.na(lab_data$co_tox_lab_normal_low[lab_data$co_tox_lab_name == 17 & lab_data$siteid == 2])) # 3882 NAs for anemia
nrow(lab_data[lab_data$co_tox_lab_name == 17,]) # 750356 rows
length(unique(lab_data$studyid[is.na(lab_data$co_tox_lab_normal_low) & lab_data$co_tox_lab_name == 17]))
# 1507 unique patients had NA for anemia


# check chart data
a <- chart_data[is.na(chart_data$toxicitydate), c("studyid", "toxicityname", "toxicitydate", "dx_dt", "surgery_d", "surgery_dt", "first_chemo_dt")] %>% print(n=Inf, width = Inf)
# no anemia
b <- chart_data[!is.na(chart_data$toxicitydate) & (chart_data$toxicitydate < chart_data$first_chemo_dt), c("studyid", "toxicityid", "toxicityname", "toxicitydate", "dx_dt", "surgery_d", "surgery_dt", "first_chemo_dt", "last_chemo_dt")] %>% print(n=Inf, width = Inf)
# 3 anemia before first chemo
# exclude those records
c <- chart_data[!is.na(chart_data$toxicitydate) & (chart_data$toxicitydate == chart_data$first_chemo_dt), c("studyid", "toxicityid", "toxicityname", "toxicitydate", "dx_dt", "surgery_d", "surgery_dt", "first_chemo_dt")] %>% print(n=Inf, width = Inf)
# 1 anemia on the first chemo date
# exclude this record

chart_data <- chart_data[!(chart_data$toxicityid %in% c("marked-out id", "marked-out id", "marked-out id", "marked-out id")),]

## Adding Time Frame
# time frame for surgery before chemo
# 1: <= date of cancer diagnosis --> pre-existing conditions (1 year)
# 2: <= date of surgery
# 3: <= date of first chemo
# 4: <= end of chemo
icd_data_before <- icd_data %>%
  filter(group == "Surgery before chemo") %>%
  mutate(time_frame = ifelse(co_tox_dx_code_date <= dx_dt, 1,
                             ifelse(co_tox_dx_code_date <= surgery_dt, 2,
                                    ifelse(co_tox_dx_code_date <= first_chemo_dt, 3,
                                           ifelse(co_tox_dx_code_date <= last_chemo_dt, 4, NA))))) %>%
  filter(!is.na(time_frame)) # remove NAs in time frame

table(icd_data_before$time_frame, useNA = "ifany")
table(icd_data_before$anemia_flag_icd, useNA = "ifany") 

lab_data_before <- lab_data %>%
  filter(group == "Surgery before chemo") %>%
  mutate(time_frame = ifelse(co_tox_lab_date <= dx_dt, 1,
                             ifelse(co_tox_lab_date <= surgery_dt, 2,
                                    ifelse(co_tox_lab_date <= first_chemo_dt, 3,
                                           ifelse(co_tox_lab_date <= last_chemo_dt, 4, NA))))) %>%
  filter(!is.na(time_frame)) # remove NAs in time frame

table(lab_data_before$time_frame, useNA = "ifany")
table(lab_data_before$co_tox_lab_name, useNA = "ifany") 


chart_before <- chart_data %>% filter(surgery_d == 1) %>%
  filter(!is.na(surgery_dt)) %>%
  filter(surgery_dt < first_chemo_dt) %>%
  filter(toxicitydate > first_chemo_dt & toxicitydate <= last_chemo_dt) # limit to last timeframe

## check
# check normal low in lab_data_before
sum(is.na(lab_data_before$co_tox_lab_normal_low[lab_data_before$co_tox_lab_name == 17 ])) # 520 NAs for anemia
sum(is.na(lab_data_before$co_tox_lab_normal_low[lab_data_before$co_tox_lab_name == 17 & lab_data_before$siteid == 2])) # 516 NAs for anemia
nrow(lab_data_before[lab_data_before$co_tox_lab_name == 17,]) # 149796
length(unique(lab_data_before$studyid[is.na(lab_data_before$co_tox_lab_normal_low) & lab_data_before$co_tox_lab_name == 17]))
# 111 unique patients had NA for anemia
lab_data_before %>% filter(co_tox_lab_name == 17) %>%
  group_by(time_frame, co_tox_lab_name) %>%
  summarise(NAs = sum(is.na(co_tox_lab_normal_low)), 
            NAs_patient = length(unique(studyid[is.na(co_tox_lab_normal_low)])),
            N = n(),
            N_patient = length(unique(studyid)))
# most of NAs in timeframe=4

lab_kpwa_sum <- lab_data_before %>% filter(co_tox_lab_name ==17 & siteid == 2) %>%
  group_by(time_frame, co_tox_lab_name) %>%
  summarise(NAs = sum(is.na(co_tox_lab_normal_low)), 
            N = n(),
            NAs_patient = length(unique(studyid[is.na(co_tox_lab_normal_low)])),
            N_patient = length(unique(studyid))) %>%
  mutate(Percent_N = NAs/N*100, Percent_patient = NAs_patient/N_patient*100)

# summary of normal low in KPNC
summary(lab_data_before$co_tox_lab_normal_low[lab_data_before$co_tox_lab_name == 17 & lab_data_before$siteid==1 & !is.na(lab_data_before$co_tox_lab_result)])
#     Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   11.00   11.00   11.50   11.32   11.50   13.00       4
# summary of normal low in KPWA
summary(lab_data_before$co_tox_lab_normal_low[lab_data_before$co_tox_lab_name == 17 & lab_data_before$siteid==2])
#     Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   10.00   11.40   11.40   11.15   11.40   11.60     516 

# change normal low <=10 to NA
lab_data_before$co_tox_lab_normal_low[lab_data_before$co_tox_lab_name==17 & lab_data_before$co_tox_lab_normal_low <=10] <- NA

# check missing in anemia for KPWA
sum(is.na(lab_data_before$co_tox_lab_result[lab_data_before$co_tox_lab_name == 17 & lab_data_before$siteid==2])) # 0 NAs

# patients with any missingness in normal low
a <- unique(lab_data_before$studyid[is.na(lab_data_before$co_tox_lab_normal_low) & lab_data_before$co_tox_lab_name ==17])
lab_missing <- lab_data_before %>% filter(studyid %in% a & co_tox_lab_name ==17) %>%
  group_by(studyid, time_frame, co_tox_lab_name) %>%
  mutate(count = row_number())
length(unique(lab_missing$studyid))
lab_missing_sum <- lab_missing %>% group_by(studyid, time_frame, co_tox_lab_name) %>%
  summarise(N= n(), NAs = sum(is.na(co_tox_lab_normal_low))) %>%
  mutate(missing_ind = ifelse(N == NAs, "full", ifelse(NAs == 0, "none", "partial")))
b <- lab_missing_sum %>% group_by(time_frame, co_tox_lab_name, missing_ind) %>%
  summarize(N_missing = n()) %>% 
  filter(missing_ind != "none")

summary(lab_missing$co_tox_lab_result[lab_missing$co_tox_lab_name == 17 & is.na(lab_missing$co_tox_lab_normal_low)])


## Imputate Normal Low
lab_data_before <- lab_data_before %>%
  mutate(normal_low_missing = ifelse(is.na(co_tox_lab_normal_low), 1, 0)) # normal low missingness indicator

lab_data_before_anemia <- lab_data_before %>%
  filter(co_tox_lab_name == 17) %>% # only anemia
  filter(!is.na(co_tox_lab_result)) # remove NAs in lab result

# patients with partial missingness of normal low in specific time frame
lab_missing_sum %>% filter(missing_ind == "partial") %>% print(n = Inf)
d <- unique(lab_missing_sum$studyid[lab_missing_sum$missing_ind == "partial"])
e <- lab_missing[lab_missing$studyid %in% d & !is.na(lab_missing$co_tox_lab_normal_low) & lab_missing$co_tox_lab_name ==17, c("studyid", "time_frame", "co_tox_lab_name", "co_tox_lab_result", "co_tox_lab_normal_low", "co_tox_lab_date")] %>% print(n = Inf, width = Inf)
f <- e %>% group_by(studyid, time_frame, co_tox_lab_name) %>%
  summarise(min = min(co_tox_lab_normal_low), 
            median = median(co_tox_lab_normal_low), 
            max = max(co_tox_lab_normal_low), n = n()) %>%
  mutate(Ind = ifelse(min == max, 1, 0))
f %>% print(n=Inf)
f[f$Ind==0,] %>% print(n=Inf)
## they had different normal low and partial missing;
lab_missing[lab_missing$studyid %in% c("marked-out id", "marked-out id", "marked-out id", "marked-out id", "marked-out id", "marked-out id", "marked-out id") & lab_missing$co_tox_lab_name ==17, c("studyid", "time_frame", "co_tox_lab_name", "co_tox_lab_result", "co_tox_lab_normal_low", "co_tox_lab_date")] %>% print(n = Inf, width = Inf)
lab_data_before_anemia$co_tox_lab_normal_low[lab_data_before_anemia$studyid == "marked-out id" & lab_data_before_anemia$co_tox_lab_name == 17 & lab_data_before_anemia$time_frame == 4 & is.na(lab_data_before_anemia$co_tox_lab_normal_low)] <- 11.4
lab_data_before_anemia$co_tox_lab_normal_low[lab_data_before_anemia$studyid == "marked-out id" & lab_data_before_anemia$co_tox_lab_name == 17 & lab_data_before_anemia$time_frame == 1 & is.na(lab_data_before_anemia$co_tox_lab_normal_low)] <- 11.3
lab_data_before_anemia$co_tox_lab_normal_low[lab_data_before_anemia$studyid == "marked-out id" & lab_data_before_anemia$co_tox_lab_name == 17 & lab_data_before_anemia$time_frame == 3 & is.na(lab_data_before_anemia$co_tox_lab_normal_low)] <- 11
## other patients had same normal low values
g <- d
for (i in 1:length(g)){
  for (j in 17) {
    for (k in 1:4){
      t <- length(lab_data_before_anemia$co_tox_lab_normal_low[lab_data_before_anemia$studyid == g[i] & lab_data_before_anemia$co_tox_lab_name == j & lab_data_before_anemia$time_frame == k & is.na(lab_data_before_anemia$co_tox_lab_normal_low)])
      s <- length(f$max[f$studyid == g[i] & f$co_tox_lab_name == j & f$time_frame == k])
      if (t != 0 & s != 0) {
        lab_data_before_anemia$co_tox_lab_normal_low[lab_data_before_anemia$studyid == g[i] & lab_data_before_anemia$co_tox_lab_name == j & lab_data_before_anemia$time_frame == k & is.na(lab_data_before_anemia$co_tox_lab_normal_low)] <- f$max[f$studyid == g[i] & f$co_tox_lab_name == j & f$time_frame == k]
      }
    }
  }
}
lab_data_before_anemia[lab_data_before_anemia$studyid=="marked-out id" & lab_data_before_anemia$time_frame == 3 & lab_data_before_anemia$co_tox_lab_name == 17, "co_tox_lab_normal_low"] %>% print(n = Inf)
# summary of normal low in KPNC
summary(lab_data_before_anemia$co_tox_lab_normal_low[lab_data_before_anemia$co_tox_lab_name == 17 & lab_data_before_anemia$siteid==1])
#     Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   11.00   11.00   11.50   11.32   11.50   13.00 
# summary of normal low in KPWA
summary(lab_data_before_anemia$co_tox_lab_normal_low[lab_data_before_anemia$co_tox_lab_name == 17 & lab_data_before_anemia$siteid==2])
#     Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    11.3    11.4    11.4    11.4    11.4    11.6    1911 
# median normal low in KPWA: 11.4
lab_data_before_anemia$co_tox_lab_normal_low[lab_data_before_anemia$co_tox_lab_name == 17 & is.na(lab_data_before_anemia$co_tox_lab_normal_low)] <- 11.4
sum(is.na(lab_data_before_anemia$co_tox_lab_normal_low)) # 0


## Flag Anemia 
lab_data_before_anemia <- lab_data_before_anemia %>%
  mutate(
    anemia_G3_lab = ifelse(co_tox_lab_name == 17 & co_tox_lab_result < 8, 1, 0),
    anemia_G2_lab = ifelse(co_tox_lab_name == 17 & co_tox_lab_result < 10, 1, 0),
    anemia_G1_lab = ifelse(co_tox_lab_name == 17 & co_tox_lab_result < co_tox_lab_normal_low, 1, 0)
  )

## Summarize Anemia Toxicity 
icd_data_before_sum_anemia <- icd_data_before %>% 
  select(studyid, time_frame, anemia_flag_icd) %>%
  group_by(studyid, time_frame) %>%
  summarise(anemia_icd_n = sum(anemia_flag_icd))
icd_data_before_sum_anemia <- icd_data_before_sum_anemia %>% mutate(anemia_flag_icd = ifelse(anemia_icd_n == 0, 0, 1))
icd_data_before_sum_anemia$siteid <- ifelse(icd_data_before_sum_anemia$studyid>"marked-out id", 2, 1)
table(icd_data_before_sum_anemia$anemia_icd_n, useNA = "ifany") # check
table(icd_data_before_sum_anemia$anemia_flag_icd, useNA = "ifany") # check
icd_data_before_sum_anemia[icd_data_before_sum_anemia$anemia_icd_n >= 1,] # check

lab_data_before_sum_anemia <- lab_data_before_anemia %>%
  select(studyid, time_frame, anemia_G3_lab, anemia_G2_lab, anemia_G1_lab) %>%
  group_by(studyid, time_frame) %>%
  summarise(
    anemia_G3_n = sum(anemia_G3_lab),
    anemia_G2_n = sum(anemia_G2_lab),
    anemia_G1_n = sum(anemia_G1_lab)
  ) %>%
  mutate(
    anemia_G3_lab = ifelse(anemia_G3_n == 0, 0, 1),
    anemia_G2_lab = ifelse(anemia_G2_n == 0, 0, 1),
    anemia_G1_lab = ifelse(anemia_G1_n == 0, 0, 1)
  )
lab_data_before_sum_anemia$siteid <- ifelse(lab_data_before_sum_anemia$studyid>"marked-out id", 2, 1)
table(lab_data_before_sum_anemia$anemia_G3_n, useNA = "ifany") # check
table(lab_data_before_sum_anemia$anemia_G2_n, useNA = "ifany") # check
table(lab_data_before_sum_anemia$anemia_G1_n, useNA = "ifany") # check
table(lab_data_before_sum_anemia$anemia_G3_lab, useNA = "ifany")
table(lab_data_before_sum_anemia$anemia_G2_lab, useNA = "ifany")
table(lab_data_before_sum_anemia$anemia_G1_lab, useNA = "ifany")


chart_before_sum_anemia <- chart_before %>%
  mutate(anemia_d = ifelse(toxicityname == "Anemia", 1, 0)) %>%
  group_by(studyid) %>%
  summarise(anemia_n = sum(anemia_d))
chart_before_sum_anemia <- chart_before_sum_anemia %>% mutate(anemia_tox = ifelse(anemia_n == 0, 0, 1))
table(chart_before_sum_anemia$anemia_tox, useNA = "ifany") # 1=185, 0=684
table(chart_before_sum_anemia$anemia_n, useNA = "ifany")


## Merge ICD and Lab
lab_data_before_sum_anemia$missing_lab <- rep(0)
icd_data_before_sum_anemia$missing_icd <- rep(0)
lab_id_1 <- lab_data_before_sum_anemia$studyid[lab_data_before_sum_anemia$time_frame==1]
lab_id_2 <- lab_data_before_sum_anemia$studyid[lab_data_before_sum_anemia$time_frame==2]
lab_id_3 <- lab_data_before_sum_anemia$studyid[lab_data_before_sum_anemia$time_frame==3]
lab_id_4 <- lab_data_before_sum_anemia$studyid[lab_data_before_sum_anemia$time_frame==4]
icd_id_1 <- icd_data_before_sum_anemia$studyid[icd_data_before_sum_anemia$time_frame==1]
icd_id_2 <- icd_data_before_sum_anemia$studyid[icd_data_before_sum_anemia$time_frame==2]
icd_id_3 <- icd_data_before_sum_anemia$studyid[icd_data_before_sum_anemia$time_frame==3]
icd_id_4 <- icd_data_before_sum_anemia$studyid[icd_data_before_sum_anemia$time_frame==4]
sum(!(lab_id_1 %in% icd_id_1)) 
sum(!(lab_id_2 %in% icd_id_2)) 
sum(!(lab_id_3 %in% icd_id_3)) 
sum(!(lab_id_4 %in% icd_id_4)) 
sum(!(icd_id_1 %in% lab_id_1)) 
sum(!(icd_id_2 %in% lab_id_2)) 
sum(!(icd_id_3 %in% lab_id_3)) 
sum(!(icd_id_4 %in% lab_id_4)) 

a <- data.frame(
  studyid = c(icd_id_1[!(icd_id_1 %in% lab_id_1)],
              icd_id_2[!(icd_id_2 %in% lab_id_2)],
              icd_id_3[!(icd_id_3 %in% lab_id_3)],
              icd_id_4[!(icd_id_4 %in% lab_id_4)]),
  time_frame = c(rep(1, times = 1416),
                 rep(2, times = 1262),
                 rep(3, times = 99),
                 rep(4, times = 32)),
  missing_lab = rep(1), 
  anemia_G3_lab = rep(0),
  anemia_G2_lab = rep(0),
  anemia_G1_lab = rep(0)
)

lab_data_before_sum_anemia <- rbind(lab_data_before_sum_anemia, a)
lab_data_before_sum_anemia$siteid <- ifelse(lab_data_before_sum_anemia$studyid>"marked-out id", 2, 1)
anyNA(lab_data_before_sum_anemia$studyid)
anyNA(lab_data_before_sum_anemia$time_frame)
anyNA(lab_data_before_sum_anemia$anemia_G3_lab)

b <- data.frame(
  studyid = c(lab_id_1[!(lab_id_1 %in% icd_id_1)],
              lab_id_2[!(lab_id_2 %in% icd_id_2)],
              lab_id_3[!(lab_id_3 %in% icd_id_3)],
              lab_id_4[!(lab_id_4 %in% icd_id_4)]),
  time_frame = c(rep(1, times = 3050),
                 rep(2, times = 4347),
                 rep(3, times = 7828),
                 rep(4, times = 5885)),
  missing_icd = rep(1),
  anemia_flag_icd = rep(0)
)
icd_data_before_sum_anemia <- rbind(icd_data_before_sum_anemia, b)
icd_data_before_sum_anemia$siteid <- ifelse(icd_data_before_sum_anemia$studyid>"marked-out id", 2, 1)
anyNA(icd_data_before_sum_anemia$studyid)
anyNA(icd_data_before_sum_anemia$anemia_flag_icd)

icd_data_before_sum_anemia <- icd_data_before_sum_anemia %>% arrange(studyid, time_frame)
lab_data_before_sum_anemia <- lab_data_before_sum_anemia %>% arrange(studyid, time_frame)
length(unique(lab_data_before_sum_anemia$studyid))
length(unique(icd_data_before_sum_anemia$studyid))
identical(icd_data_before_sum_anemia$studyid, lab_data_before_sum_anemia$studyid) # T
identical(icd_data_before_sum_anemia$time_frame, lab_data_before_sum_anemia$time_frame) # T

# merge icd and lab
before_sum_anemia <- icd_data_before_sum_anemia %>%
  inner_join(lab_data_before_sum_anemia[, names(lab_data_before_sum_anemia)[c(1:8, 10)]], by = join_by(studyid, time_frame))
adju_cohort <- readRDS(file = paste0(datapath, "adju_cohort.rds"))
before_sum_anemia <- before_sum_anemia[, names(before_sum_anemia)[c(1:4, 6:13)]] %>%
  right_join(adju_cohort, by = join_by(studyid, time_frame))
glimpse(before_sum_anemia)
table(before_sum_anemia$anemia_flag_icd, useNA = "ifany")
before_sum_anemia$anemia_flag_icd[is.na(before_sum_anemia$anemia_flag_icd)] <- 0
table(before_sum_anemia$anemia_G3_lab, useNA = "ifany")
before_sum_anemia$anemia_G3_lab[is.na(before_sum_anemia$anemia_G3_lab)] <- 0
table(before_sum_anemia$anemia_G2_lab, useNA = "ifany")
before_sum_anemia$anemia_G2_lab[is.na(before_sum_anemia$anemia_G2_lab)] <- 0
table(before_sum_anemia$anemia_G1_lab, useNA = "ifany")
before_sum_anemia$anemia_G1_lab[is.na(before_sum_anemia$anemia_G1_lab)] <- 0
table(before_sum_anemia$missing_icd, useNA = "ifany")
before_sum_anemia$missing_icd[is.na(before_sum_anemia$missing_icd)] <- 1
table(before_sum_anemia$missing_lab, useNA = "ifany")
before_sum_anemia$missing_lab[is.na(before_sum_anemia$missing_lab)] <- 1

# create anemia_tox_indicator
before_sum_anemia <- before_sum_anemia %>%
  mutate(anemia_tox_indicator = ifelse(anemia_flag_icd == 1 & anemia_G2_lab == 1, "Both",
                                       ifelse(anemia_flag_icd == 1 & anemia_G2_lab == 0, "Only ICD",
                                              ifelse(anemia_flag_icd == 0 & anemia_G2_lab == 1, "Only Grade 2", "Neither"))))

# subset before_sum_anemia and export
adjuvant_anemia <- before_sum_anemia[, c("studyid", "time_frame", "anemia_flag_icd", "anemia_G3_lab", "anemia_G2_lab", "anemia_G1_lab",
                                         "missing_icd", "missing_lab", "dx_dt", "siteid", "first_chemo_dt", "last_chemo_dt",
                                         "surgery_d", "surgery_dt", "group", "anemia_tox_indicator")]
write.csv(adjuvant_anemia, file = paste0(datapath, "adjuvant_anemia.csv"),
          row.names = F)


## Merge IOCD, Lab, and Chart
tox_sum_anemia <- chart_before_sum_anemia[, c("studyid", "anemia_tox")] %>%
  right_join(before_sum_anemia[before_sum_anemia$time_frame == 4 & before_sum_anemia$siteid==2,], by = join_by(studyid))
glimpse(tox_sum_anemia) 
tox_sum_anemia$anemia_tox[is.na(tox_sum_anemia$anemia_tox)] <- 0


# export data
# saveRDS(icd_data_before_sum_anemia, file = paste0(datapath, "icd_before_sum_anemia.rds"))
# saveRDS(lab_data_before_anemia, file = paste0(datapath, "lab_before_anemia.rds"))
# saveRDS(lab_data_before_sum_anemia, file = paste0(datapath, "lab_before_sum_anemia.rds"))
# saveRDS(before_sum_anemia, file = paste0(datapath, "before_sum_anemia.rds"))
# saveRDS(tox_sum_anemia, file = paste0(datapath, "tox_sum_anemia.rds"))
# saveRDS(chart_before_sum_anemia, file = paste0(datapath, "chart_before_sum_anemia.rds"))
# before_sum_anemia <- readRDS(file = paste0(datapath, "before_sum_anemia.rds"))


## Neuropathy ##
## Load Dataset
icd_data <- readRDS(file = paste0(datapath, "icd.rds"))
chart_data <- readRDS(file = paste0(datapath, "chart.rds"))
obcd_group <- readRDS(file = paste0(datapath, "obcd.rds"))

## Check and Clean
length(unique(icd_data$studyid)) 
length(unique(icd_data$studyid[icd_data$any_chemo==1])) 
anyNA(icd_data$co_tox_dx_code_date) # F
table(icd_data$neuropathy_flag_icd, useNA = "ifany") 
table(icd_data$neuropathy_flag_icd_reg, useNA = "ifany") 
table(icd_data$neuropathy_flag_icd_se, useNA = "ifany") 


# check chart data
a <- chart_data[is.na(chart_data$toxicitydate), c("studyid", "toxicityname", "toxicitydate", "dx_dt", "surgery_d", "surgery_dt", "first_chemo_dt")] %>% print(n=Inf, width = Inf)
## we have missing toxicity date for neuropathy
b <- chart_data[!is.na(chart_data$toxicitydate) & (chart_data$toxicitydate < chart_data$first_chemo_dt), c("studyid", "toxicityid", "toxicityname", "toxicitydate", "dx_dt", "surgery_d", "surgery_dt", "first_chemo_dt", "last_chemo_dt")] %>% print(n=Inf, width = Inf)
## ## we have neuropathy date before first chemo date
c <- chart_data[!is.na(chart_data$toxicitydate) & (chart_data$toxicitydate == chart_data$first_chemo_dt), c("studyid", "toxicityname", "toxicitydate", "dx_dt", "surgery_d", "surgery_dt", "first_chemo_dt")] %>% print(n=Inf, width = Inf)
## we have neuropathy date same as first chemo date
## exclude those records

a <- chart_data[!(is.na(chart_data$toxicitydate) & chart_data$toxicityname == "Neuropathy"),]
chart_data <- a
a <- chart_data[!((chart_data$toxicitydate <= chart_data$first_chemo_dt) & chart_data$toxicityname == "Neuropathy"),]
chart_data <- a

## Adding Time Frame 
# time frame for surgery before chemo
# 1: <= date of cancer diagnosis --> pre-existing conditions (1 year)
# 2: <= date of surgery
# 3: <= date of first chemo
# 4: <= end of chemo
icd_data_before <- icd_data %>%
  filter(group == "Surgery before chemo") %>%
  mutate(time_frame = ifelse(co_tox_dx_code_date <= dx_dt, 1,
                             ifelse(co_tox_dx_code_date <= surgery_dt, 2,
                                    ifelse(co_tox_dx_code_date <= first_chemo_dt, 3,
                                           ifelse(co_tox_dx_code_date <= last_chemo_dt, 4, NA))))) %>%
  filter(!is.na(time_frame)) # remove NAs in time frame

table(icd_data_before$time_frame, useNA = "ifany")
table(icd_data_before$neuropathy_flag_icd, useNA = "ifany") 
table(icd_data_before$neuropathy_flag_icd_reg, useNA = "ifany") 
table(icd_data_before$neuropathy_flag_icd_se, useNA = "ifany")


chart_before <- chart_data %>% filter(surgery_d == 1) %>%
  filter(!is.na(surgery_dt)) %>%
  filter(surgery_dt < first_chemo_dt) %>%
  filter(toxicitydate > first_chemo_dt & toxicitydate <= last_chemo_dt) # limit to last timeframe



## Summarize Neuropathy Toxicity
icd_data_before_sum_neuro <- icd_data_before %>% 
  select(studyid, time_frame, neuropathy_flag_icd, neuropathy_flag_icd_reg, neuropathy_flag_icd_se) %>%
  group_by(studyid, time_frame) %>%
  summarise(neuropathy_icd_n = sum(neuropathy_flag_icd),
            neuropathy_reg_icd_n = sum(neuropathy_flag_icd_reg),
            neuropathy_se_icd_n = sum(neuropathy_flag_icd_se))
icd_data_before_sum_neuro <- icd_data_before_sum_neuro %>% 
  mutate(neuropathy_flag_icd = ifelse(neuropathy_icd_n == 0, 0, 1),
         neuropathy_flag_icd_reg = ifelse(neuropathy_reg_icd_n == 0, 0, 1),
         neuropathy_flag_icd_se = ifelse(neuropathy_se_icd_n == 0, 0, 1))
icd_data_before_sum_neuro$siteid <- ifelse(icd_data_before_sum_neuro$studyid>"marked-out id", 2, 1)
table(icd_data_before_sum_neuro$neuropathy_icd_n, useNA = "ifany") # check
table(icd_data_before_sum_neuro$neuropathy_flag_icd, useNA = "ifany") # check
icd_data_before_sum_neuro[icd_data_before_sum_neuro$neuropathy_icd_n >= 1,] # check

chart_before_sum_neuro <- chart_before %>%
  mutate(neuro_d = ifelse(toxicityname == "Neuropathy", 1, 0)) %>%
  group_by(studyid) %>%
  summarise(neuro_n = sum(neuro_d))
chart_before_sum_neuro <- chart_before_sum_neuro %>% 
  mutate(neuropathy_tox = ifelse(neuro_n == 0, 0, 1))
table(chart_before_sum_neuro$neuropathy_tox, useNA = "ifany") # 1: 389; 0: 480
table(chart_before_sum_neuro$neuro_n, useNA = "ifany")

# find the date of first reported neuropathy
icd_neuro_date <- icd_data_before %>%
  filter(time_frame == 4 & neuropathy_flag_icd == 1) %>%
  select(studyid, co_tox_dx_code_date, dx_dt, first_chemo_dt, last_chemo_dt) %>%
  arrange(studyid, co_tox_dx_code_date) %>%
  group_by(studyid) %>%
  mutate(Count = row_number()) %>%
  filter(Count == 1) %>%
  select(studyid, co_tox_dx_code_date, dx_dt, first_chemo_dt, last_chemo_dt)
names(icd_neuro_date) <- c("studyid", "first_neuropathy_icd_date", "dx_icd", "first_chemo_icd", "last_chemo_icd")

chart_neuro_date <- chart_before %>%
  filter(toxicityname == "Neuropathy") %>%
  select(studyid, toxicitydate, dx_dt, first_chemo_dt, last_chemo_dt) %>%
  arrange(studyid, toxicitydate) %>%
  group_by(studyid) %>%
  mutate(Count = row_number()) %>%
  filter(Count == 1) %>%
  select(studyid, toxicitydate, dx_dt, first_chemo_dt, last_chemo_dt)
names(chart_neuro_date) <- c("studyid", "first_neuropathy_tox_date", "dx_chart", "first_chemo_chart", "last_chemo_chart")


## Merge ICD and Chart
length(unique(chart_before_sum_neuro$studyid)) 
length(unique(icd_data_before_sum_neuro$studyid[icd_data_before_sum_neuro$time_frame==4 & icd_data_before_sum_neuro$siteid==2])) 
sum(icd_data_before_sum_neuro$studyid[icd_data_before_sum_neuro$time_frame==4 & icd_data_before_sum_neuro$siteid==2] %in% chart_before_sum_neuro$studyid) 
tox_sum_neuro <- chart_before_sum_neuro[, c("studyid", "neuropathy_tox")] %>%
  full_join(icd_data_before_sum_neuro[icd_data_before_sum_neuro$time_frame == 4 & icd_data_before_sum_neuro$siteid==2,], by = join_by(studyid)) 
glimpse(tox_sum_neuro) 

tox_sum_neuro <- tox_sum_neuro %>%
  left_join(icd_neuro_date, by = join_by(studyid)) %>%
  left_join(chart_neuro_date, by = join_by(studyid)) 

anyNA(tox_sum_neuro$neuropathy_tox)
tox_sum_neuro$neuropathy_tox[is.na(tox_sum_neuro$neuropathy_tox)] <- 0

tox_sum_neuro$neuropathy_flag_icd[is.na(tox_sum_neuro$neuropathy_flag_icd)] <- 0
tox_sum_neuro$neuropathy_flag_icd_reg[is.na(tox_sum_neuro$neuropathy_flag_icd_reg)] <- 0
tox_sum_neuro$neuropathy_flag_icd_se[is.na(tox_sum_neuro$neuropathy_flag_icd_se)] <- 0
tox_sum_neuro$gap_between_icd_and_tox <- as.numeric(tox_sum_neuro$first_neuropathy_icd_date - tox_sum_neuro$first_neuropathy_tox_date)
tox_sum_neuro$gap_between_icd_and_first_chemo <- as.numeric(tox_sum_neuro$first_neuropathy_icd_date - tox_sum_neuro$first_chemo_icd)
tox_sum_neuro$gap_between_chart_and_first_chemo <- as.numeric(tox_sum_neuro$first_neuropathy_tox_date - tox_sum_neuro$first_chemo_chart)

# export data
# saveRDS(icd_data_before_sum_neuro, file = paste0(datapath, "icd_before_sum_neuro.rds"))
# saveRDS(tox_sum_neuro, file = paste0(datapath, "tox_sum_neuro.rds"))
# saveRDS(chart_before_sum_neuro, file = paste0(datapath, "chart_before_sum_neuro.rds"))



## Neutropenia ##
## Load Dataset
icd_data <- readRDS(file = paste0(datapath, "icd.rds"))
lab_data <- readRDS(file = paste0(datapath, "lab.rds"))
chart_data <- readRDS(file = paste0(datapath, "chart.rds"))
obcd_group <- readRDS(file = paste0(datapath, "obcd.rds"))

## Check and Clean 
length(unique(icd_data$studyid)) 
length(unique(icd_data$studyid[icd_data$any_chemo==1])) 
anyNA(icd_data$co_tox_dx_code_date) # F
table(icd_data$neutropenia_flag_icd, useNA = "ifany") 
length(unique(lab_data$studyid)) 
length(unique(lab_data$studyid[lab_data$any_chemo==1])) 
anyNA(lab_data$co_tox_lab_date) #F
anyNA(lab_data$co_tox_lab_name) #F
sum(is.na(lab_data$co_tox_lab_result[lab_data$co_tox_lab_name == 1 ])) # no NA
sum(is.na(lab_data$co_tox_lab_result[lab_data$co_tox_lab_name == 2 ])) # 1366 NA

table(lab_data$co_tox_lab_unit[lab_data$co_tox_lab_name==1], useNA = "ifany") # /ul, ul, k/ul, percent, unknown
# /ul     ul   k/ul   percent   unknown 
# 108    190 448655   1011         1
table(lab_data$co_tox_lab_unit[lab_data$co_tox_lab_name==2], useNA = "ifany") # /ul, ul, k/ul, NA
#     /ul     ul   k/ul   <NA> 
#   36033  10837 682465     26 
table(lab_data$co_tox_lab_unit[lab_data$co_tox_lab_name==21], useNA = "ifany") # ?
summary(lab_data$co_tox_lab_result[lab_data$co_tox_lab_name == 1 ]) 
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.00     2.70     3.90     6.95     5.60    32634.00 
summary(lab_data$co_tox_lab_result[lab_data$co_tox_lab_name == 21 ]) # >100 as NAs
# Min.  1st Qu.  Median    Mean  3rd Qu.    Max.    NA's 
# 0.00   54.00   63.00    62.58   72.00    187.00    1362 
summary(lab_data$co_tox_lab_result[lab_data$co_tox_lab_name == 1 & lab_data$co_tox_lab_unit == 1]) 
summary(lab_data$co_tox_lab_result[lab_data$co_tox_lab_name == 1 & lab_data$co_tox_lab_unit == 2]) # similar to 1
summary(lab_data$co_tox_lab_result[lab_data$co_tox_lab_name == 1 & lab_data$co_tox_lab_unit == 3]) # different
summary(lab_data$co_tox_lab_result[lab_data$co_tox_lab_name == 1 & lab_data$co_tox_lab_unit == 4]) # different
summary(lab_data$co_tox_lab_result[lab_data$co_tox_lab_name == 1 & lab_data$co_tox_lab_unit == 11]) # different
# 1 K/uL = 1000 uL = 1000 mm^3
sum(lab_data$co_tox_lab_result[lab_data$co_tox_lab_name == 1] == 0 & !is.na(lab_data$co_tox_lab_result[lab_data$co_tox_lab_name == 1])) # 562; 564

# Normal Low
sum(is.na(lab_data$co_tox_lab_normal_low[lab_data$co_tox_lab_name == 1 ])) # 3708 NAs
sum(is.na(lab_data$co_tox_lab_normal_low[lab_data$co_tox_lab_name == 1 & lab_data$siteid == 2])) # 3707 NAs for anemia
summary(lab_data$co_tox_lab_normal_low[lab_data$co_tox_lab_name == 1 & lab_data$co_tox_lab_unit == 1]) # 1000 /ul
summary(lab_data$co_tox_lab_normal_low[lab_data$co_tox_lab_name == 1 & lab_data$co_tox_lab_unit == 2]) # 1000 ul
summary(lab_data$co_tox_lab_normal_low[lab_data$co_tox_lab_name == 1 & lab_data$co_tox_lab_unit == 3]) # k/ul
#     Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   1.800   1.800   2.100   1.991   2.100   2.100    3138
summary(lab_data$co_tox_lab_normal_low[lab_data$co_tox_lab_name == 1 & lab_data$co_tox_lab_unit == 4]) # percent
#     Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   38.00   40.00   40.00   39.87   40.00   40.00     569 
summary(lab_data$co_tox_lab_normal_low[lab_data$co_tox_lab_name == 1 & lab_data$co_tox_lab_unit == 11]) # NA

sum(is.na(lab_data$co_tox_lab_normal_low[lab_data$co_tox_lab_name == 2 ])) # 5261 NAs
sum(is.na(lab_data$co_tox_lab_normal_low[lab_data$co_tox_lab_name == 2 & lab_data$siteid == 2])) # 5260 NAs for anemia
summary(lab_data$co_tox_lab_normal_low[lab_data$co_tox_lab_name == 2 & lab_data$co_tox_lab_unit == 1]) # /ul
#     Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#       0    4000    4000    3997    4000    4000      26
lab_data[lab_data$co_tox_lab_name == 2 & lab_data$co_tox_lab_unit == 1 & lab_data$co_tox_lab_normal_low<4000 & lab_data$group == "Surgery before chemo",
         c("studyid", "co_tox_lab_date", "co_tox_lab_result", "co_tox_lab_normal_low", "dx_dt", "group")] %>% print(n = Inf)
# normal low < 4000: only 2 in adjuvant and their dx date <2004; so no need to worry
summary(lab_data$co_tox_lab_normal_low[lab_data$co_tox_lab_name == 2 & lab_data$co_tox_lab_unit == 2]) # ul
#     Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    3400    3400    3700    3657    3700    6000    1377
summary(lab_data$co_tox_lab_normal_low[lab_data$co_tox_lab_name == 2 & lab_data$co_tox_lab_unit == 3]) # k/ul
#     Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#     3.5     3.5     3.5     3.6     3.7     4.5    3911
summary(lab_data$co_tox_lab_normal_low[lab_data$co_tox_lab_name == 2 & is.na(lab_data$co_tox_lab_unit)]) # 4; should be k/ul

sum(is.na(lab_data$co_tox_lab_normal_low[lab_data$co_tox_lab_name == 21 ])) # 6 NAs
sum(is.na(lab_data$co_tox_lab_normal_low[lab_data$co_tox_lab_name == 21 & lab_data$siteid == 2])) # 0 NAs
summary(lab_data$co_tox_lab_normal_low[lab_data$co_tox_lab_name == 21]) # percent?
#     Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   25.00   41.00   41.00   42.46   42.00   50.00       6 

length(unique(lab_data$studyid[is.na(lab_data$co_tox_lab_normal_low) & lab_data$co_tox_lab_name %in% c(1, 2, 21)]))
# 1995 unique patients had NA


## deal with percentage: result
a <- lab_data[(lab_data$co_tox_lab_name == 1 & lab_data$co_tox_lab_unit == 4) | (lab_data$co_tox_lab_name == 21 & lab_data$co_tox_lab_result <= 100), 
              c("studyid", "co_tox_lab_date", "co_tox_lab_result", "co_tox_lab_name", "co_tox_lab_unit")] %>% print(width = Inf)
a <- a %>% group_by(studyid, co_tox_lab_date) %>% mutate(Count = row_number()) %>% filter(Count == 1) %>% select(studyid, co_tox_lab_date, co_tox_lab_result, co_tox_lab_name, co_tox_lab_unit) %>% print(width = Inf)
summary(a$co_tox_lab_result)
b <- lab_data[lab_data$studyid %in% a$studyid & lab_data$co_tox_lab_name == 2 & lab_data$co_tox_lab_date %in% a$co_tox_lab_date, 
              c("studyid", "co_tox_lab_date", "co_tox_lab_result", "co_tox_lab_name", "co_tox_lab_unit")] %>% print(width = Inf)
table(b$co_tox_lab_unit, useNA = "ifany")
summary(b$co_tox_lab_result)
summary(b$co_tox_lab_result[b$co_tox_lab_unit==1])
#     Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#     200    5900    7500    7914    9400   31800       1 
summary(b$co_tox_lab_result[b$co_tox_lab_unit==3])
#     Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   0.000   4.800   6.400   7.164   8.400 359.500       1 
b[is.na(b$co_tox_lab_unit),] # K/uL
b$co_tox_lab_result[is.na(b$co_tox_lab_unit) | b$co_tox_lab_unit == 3] <- b$co_tox_lab_result[is.na(b$co_tox_lab_unit) | b$co_tox_lab_unit == 3]*1000
b <- b[, 1:3]
names(b) <- c("studyid", "co_tox_lab_date", "WBC")
b <- b %>% group_by(studyid, co_tox_lab_date) %>% mutate(Count = row_number()) %>% filter(Count == 1) %>% select(studyid, co_tox_lab_date, WBC)
a <- a %>% left_join(b, by = join_by(studyid, co_tox_lab_date))
summary(a$co_tox_lab_result)
summary(a$WBC)
a$co_tox_lab_result <- a$co_tox_lab_result*a$WBC/100
a <- a[, c("studyid", "co_tox_lab_date", "co_tox_lab_result", "co_tox_lab_name", "co_tox_lab_unit")]
a$co_tox_lab_result <- round(a$co_tox_lab_result)
names(a) <- c("studyid", "co_tox_lab_date", "lab_by_percent", "co_tox_lab_name", "co_tox_lab_unit")
lab_data <- lab_data %>% left_join(a, by = join_by(studyid, co_tox_lab_name, co_tox_lab_unit, co_tox_lab_date))
glimpse(lab_data)
sum(!is.na(lab_data$lab_by_percent))
lab_data <- lab_data %>% mutate(co_tox_lab_result = ifelse(is.na(lab_by_percent), co_tox_lab_result, lab_by_percent))
lab_data[!is.na(lab_data$lab_by_percent), c("studyid", "co_tox_lab_date", "co_tox_lab_result", "lab_by_percent")]

## deal with percentage: normal low
c <- lab_data[(lab_data$co_tox_lab_name == 1 & lab_data$co_tox_lab_unit == 4) | (lab_data$co_tox_lab_name == 21 & lab_data$co_tox_lab_normal_low <= 100), 
              c("studyid", "co_tox_lab_date", "co_tox_lab_name", "co_tox_lab_unit", "co_tox_lab_normal_low")] %>% print(width = Inf)
c <- c %>% group_by(studyid, co_tox_lab_date) %>% mutate(Count = row_number()) %>% filter(Count == 1) %>% select(studyid, co_tox_lab_date, co_tox_lab_name, co_tox_lab_unit, co_tox_lab_normal_low) %>% print(width = Inf)
summary(c$co_tox_lab_normal_low)
b <- lab_data[lab_data$studyid %in% a$studyid & lab_data$co_tox_lab_name == 2 & lab_data$co_tox_lab_date %in% a$co_tox_lab_date, 
              c("studyid", "co_tox_lab_date", "co_tox_lab_result", "co_tox_lab_name", "co_tox_lab_unit")] %>% print(width = Inf)
table(b$co_tox_lab_unit, useNA = "ifany")
summary(b$co_tox_lab_result)
summary(b$co_tox_lab_result[b$co_tox_lab_unit==1])
summary(b$co_tox_lab_result[b$co_tox_lab_unit==3])
b[is.na(b$co_tox_lab_unit),] # K/uL
b$co_tox_lab_result[is.na(b$co_tox_lab_unit) | b$co_tox_lab_unit == 3] <- b$co_tox_lab_result[is.na(b$co_tox_lab_unit) | b$co_tox_lab_unit == 3]*1000
b <- b[, 1:3]
names(b) <- c("studyid", "co_tox_lab_date", "WBC")
b <- b %>% group_by(studyid, co_tox_lab_date) %>% mutate(Count = row_number()) %>% filter(Count == 1) %>% select(studyid, co_tox_lab_date, WBC)
c <- c %>% left_join(b, by = join_by(studyid, co_tox_lab_date))
c$co_tox_lab_normal_low <- c$co_tox_lab_normal_low*c$WBC/100
c <- c[, c("studyid", "co_tox_lab_date", "co_tox_lab_name", "co_tox_lab_unit", "co_tox_lab_normal_low")]
c$co_tox_lab_normal_low <- round(c$co_tox_lab_normal_low)
names(c) <- c("studyid", "co_tox_lab_date", "co_tox_lab_name", "co_tox_lab_unit", "normal_low_by_percent")
lab_data <- lab_data %>% left_join(c, by = join_by(studyid, co_tox_lab_name, co_tox_lab_unit, co_tox_lab_date))
glimpse(lab_data)
sum(!is.na(lab_data$normal_low_by_percent))
lab_data <- lab_data %>% mutate(co_tox_lab_normal_low = ifelse(is.na(normal_low_by_percent), co_tox_lab_normal_low, normal_low_by_percent))
lab_data[!is.na(lab_data$lab_by_percent), c("studyid", "co_tox_lab_date", "co_tox_lab_result", "lab_by_percent", "normal_low_by_percent")]
lab_data[!is.na(lab_data$normal_low_by_percent), c("studyid", "co_tox_lab_date", "co_tox_lab_result", "lab_by_percent", "normal_low_by_percent")]


## converse K/uL and unknown to uL
lab_data$co_tox_lab_result[lab_data$co_tox_lab_name == 1 & lab_data$co_tox_lab_unit %in% c(3, 11)]
length(lab_data$co_tox_lab_result[lab_data$co_tox_lab_name == 1 & lab_data$co_tox_lab_unit %in% c(3, 11)])
summary(lab_data$co_tox_lab_result[lab_data$co_tox_lab_name == 1 & lab_data$co_tox_lab_unit %in% c(3)])
summary(lab_data$co_tox_lab_result[lab_data$co_tox_lab_name == 1 & lab_data$co_tox_lab_unit %in% c(11)])
lab_data$co_tox_lab_result[lab_data$co_tox_lab_name == 1 & lab_data$co_tox_lab_unit %in% c(3, 11)] <- lab_data$co_tox_lab_result[lab_data$co_tox_lab_name == 1 & lab_data$co_tox_lab_unit %in% c(3, 11)]*1000

lab_data$co_tox_lab_normal_low[lab_data$co_tox_lab_name == 1 & lab_data$co_tox_lab_unit %in% c(3, 11)]
length(lab_data$co_tox_lab_normal_low[lab_data$co_tox_lab_name == 1 & lab_data$co_tox_lab_unit %in% c(3, 11)])
summary(lab_data$co_tox_lab_normal_low[lab_data$co_tox_lab_name == 1 & lab_data$co_tox_lab_unit %in% c(3, 11)])
lab_data$co_tox_lab_normal_low[lab_data$co_tox_lab_name == 1 & lab_data$co_tox_lab_unit %in% c(3, 11)] <- lab_data$co_tox_lab_normal_low[lab_data$co_tox_lab_name == 1 & lab_data$co_tox_lab_unit %in% c(3, 11)]*1000


## deal with 0 as lab results
a <- lab_data[lab_data$co_tox_lab_name %in% c(1, 21) & lab_data$co_tox_lab_result == 0 & !is.na(lab_data$co_tox_lab_result),]
length(unique(a$studyid)) # 348 unique patients; 349
table(a$co_tox_lab_unit, useNA = "ifany") # 1 /uL, 5 uL, 556->558 K/uL, 315 (percent)
# for K/uL, 0 may be caused by rounding
a[, c("studyid", "co_tox_lab_date", "co_tox_lab_unit")] %>% print(n = Inf)
a[a$co_tox_lab_unit!=3 & a$co_tox_lab_unit!=16, ] %>% print(width = Inf) # no chemo date
a[a$co_tox_lab_unit==16, ] %>% print(width = Inf) # some have no chemo date

b <- lab_data[lab_data$studyid %in% a$studyid & lab_data$co_tox_lab_name %in% c(1, 21) & !is.na(lab_data$co_tox_lab_result),]
anyNA(b$co_tox_lab_result)
table(b$group, useNA = "ifany")
c <- b %>% group_by(studyid) %>%
  summarise(Min = min(co_tox_lab_result),
            Q1 = quantile(co_tox_lab_result, 0.25),
            Median = median(co_tox_lab_result),
            Mean = mean(co_tox_lab_result),
            Q3 = quantile(co_tox_lab_result, 0.75),
            Max = max(co_tox_lab_result),
            N = n()) %>%
  arrange(Q1, Median, Q3, N) %>%
  print(n = Inf)
d <- b %>%
  filter(group == "Surgery before chemo") %>%
  mutate(time_frame = ifelse(co_tox_lab_date <= dx_dt, 1,
                             ifelse(co_tox_lab_date <= surgery_dt, 2,
                                    ifelse(co_tox_lab_date <= first_chemo_dt, 3,
                                           ifelse(co_tox_lab_date <= last_chemo_dt, 4, NA))))) %>%
  filter(!is.na(time_frame)) # remove NAs in time frame
e <- d %>% group_by(studyid, time_frame) %>%
  summarise(Min = min(co_tox_lab_result),
            Q1 = quantile(co_tox_lab_result, 0.25),
            Median = median(co_tox_lab_result),
            Mean = mean(co_tox_lab_result),
            Q3 = quantile(co_tox_lab_result, 0.75),
            Max = max(co_tox_lab_result),
            N = n()) %>%
  arrange(Min, Q1, Median, Q3, N) %>%
  print(n = Inf)

b[b$studyid== 206523, c("studyid", "co_tox_lab_date", "co_tox_lab_result", "dx_dt", "first_chemo_dt")] %>% print(width = Inf) %>%   print(n = Inf)

## check ANC <100 for unit = 4 or lab_name = 21
lab_data[lab_data$co_tox_lab_name == 1 & lab_data$co_tox_lab_result < 100 & lab_data$co_tox_lab_unit == 4, 
         c("studyid", "co_tox_lab_date", "co_tox_lab_result", "co_tox_lab_normal_low", "dx_dt", "first_chemo_dt", "group")] %>% print(width = Inf)
# only 4 in adjuvant group
lab_data[lab_data$co_tox_lab_name == 1 & lab_data$co_tox_lab_result < 100 & lab_data$co_tox_lab_unit == 4 & lab_data$group %in% c("Surgery before chemo") & !is.na(lab_data$group), 
         c("studyid", "co_tox_lab_date", "co_tox_lab_result", "co_tox_lab_normal_low", "dx_dt", "first_chemo_dt", "group")] %>% print(width = Inf)
# why: low values or errors in WBC?

t <- lab_data %>% filter(group == "Surgery before chemo" & co_tox_lab_name %in% 21 & co_tox_lab_result < 100) %>%
  select(studyid, co_tox_lab_date, co_tox_lab_result, co_tox_lab_normal_low, dx_dt, first_chemo_dt, group)
length(unique(t$studyid))
# 912 unique patients in adjuvant group; nrow=1720
summary(t$co_tox_lab_result)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00   12.00   34.00   38.65   63.00   99.00
summary(t$co_tox_lab_normal_low)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0   210.0   405.0   485.6   650.0  4400.0 

# check normal low < 1500 (G2)
k <- lab_data %>% filter(group == "Surgery before chemo" & co_tox_lab_name %in% c(1,21) & co_tox_lab_normal_low < 1500) %>%
  select(studyid, co_tox_lab_name, co_tox_lab_date, co_tox_lab_result, co_tox_lab_normal_low, dx_dt, first_chemo_dt, group)
k <- unique(k)
length(unique(k$studyid))
nrow(k)
summary(k$co_tox_lab_normal_low)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   0     984    1189    1125    1353    1476 
# if a patient had low WBC, then using WBC to calculate normal low is incorrect; we need to revise normal low
summary(k$co_tox_lab_normal_low[k$co_tox_lab_name==1])
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 240.0  1000.0  1000.0   987.4  1080.0  1440.0 
summary(k$co_tox_lab_normal_low[k$co_tox_lab_name==21])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   0     984    1189    1126    1353    1476 

# check NA
anyNA(lab_data$co_tox_lab_result[lab_data$co_tox_lab_name %in% c(1, 21)])
lab_data$co_tox_lab_name[lab_data$co_tox_lab_name==21] <- 1


# check chart data
a <- chart_data[is.na(chart_data$toxicitydate), c("studyid", "toxicityname", "toxicitydate", "dx_dt", "surgery_d", "surgery_dt", "first_chemo_dt")] %>% print(n=Inf, width = Inf)
# no neutropenia
b <- chart_data[!is.na(chart_data$toxicitydate) & (chart_data$toxicitydate < chart_data$first_chemo_dt), c("studyid", "toxicityid", "toxicityname", "toxicitydate", "dx_dt", "surgery_d", "surgery_dt", "first_chemo_dt", "last_chemo_dt")] %>% print(n=Inf, width = Inf)
# 1 neutropenia before first chemo
# exclude those records
c <- chart_data[!is.na(chart_data$toxicitydate) & (chart_data$toxicitydate == chart_data$first_chemo_dt), c("studyid", "toxicityid", "toxicityname", "toxicitydate", "dx_dt", "surgery_d", "surgery_dt", "first_chemo_dt")] %>% print(n=Inf, width = Inf)
# 1 neutropenia on the first chemo date
# exclude this record

chart_data <- chart_data[!(chart_data$toxicityid %in% c("marked-out id", "marked-out id")),]

## Adding Time Frame
# time frame for surgery before chemo
# 1: <= date of cancer diagnosis --> pre-existing conditions (1 year)
# 2: <= date of surgery
# 3: <= date of first chemo
# 4: <= end of chemo
icd_data_before <- icd_data %>%
  filter(group == "Surgery before chemo") %>%
  mutate(time_frame = ifelse(co_tox_dx_code_date <= dx_dt, 1,
                             ifelse(co_tox_dx_code_date <= surgery_dt, 2,
                                    ifelse(co_tox_dx_code_date <= first_chemo_dt, 3,
                                           ifelse(co_tox_dx_code_date <= last_chemo_dt, 4, NA))))) %>%
  filter(!is.na(time_frame)) # remove NAs in time frame

table(icd_data_before$time_frame, useNA = "ifany")
table(icd_data_before$neutropenia_flag_icd, useNA = "ifany") 

lab_data_before <- lab_data %>%
  filter(group == "Surgery before chemo") %>%
  mutate(time_frame = ifelse(co_tox_lab_date <= dx_dt, 1,
                             ifelse(co_tox_lab_date <= surgery_dt, 2,
                                    ifelse(co_tox_lab_date <= first_chemo_dt, 3,
                                           ifelse(co_tox_lab_date <= last_chemo_dt, 4, NA))))) %>%
  filter(!is.na(time_frame)) # remove NAs in time frame

table(lab_data_before$time_frame, useNA = "ifany")
table(lab_data_before$co_tox_lab_name, useNA = "ifany") 


chart_before <- chart_data %>% filter(surgery_d == 1) %>%
  filter(!is.na(surgery_dt)) %>%
  filter(surgery_dt < first_chemo_dt) %>%
  filter(toxicitydate > first_chemo_dt & toxicitydate <= last_chemo_dt) # limit to last timeframe

## check 
# check
summary(lab_data_before$co_tox_lab_normal_low[lab_data_before$co_tox_lab_name == 1 & lab_data_before$siteid==1 & lab_data_before$time_frame==4])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 41    1800    2100    2397    2352  101514       3 
table(year(lab_data_before$co_tox_lab_date[lab_data_before$co_tox_lab_name == 1 & lab_data_before$siteid==1 & lab_data_before$co_tox_lab_normal_low <=1500 & lab_data_before$time_frame==4]))
# 2006 to 2021
summary(lab_data_before$co_tox_lab_normal_low[lab_data_before$co_tox_lab_name == 1 & lab_data_before$siteid==2 & lab_data_before$time_frame==4])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 240    2000    2000    2003    2000    8360     380 
table(year(lab_data_before$co_tox_lab_date[lab_data_before$co_tox_lab_name == 1 & lab_data_before$siteid==2 & lab_data_before$co_tox_lab_normal_low <=1500 & lab_data_before$time_frame==4]))
# 2005 to 2016

# check normal low in lab_data_before
sum(is.na(lab_data_before$co_tox_lab_normal_low[lab_data_before$co_tox_lab_name == 1 ])) # 492 NAs for neutropenia
sum(is.na(lab_data_before$co_tox_lab_normal_low[lab_data_before$co_tox_lab_name == 1 & lab_data_before$siteid == 2])) # 489 NAs for neutropenia
nrow(lab_data_before[lab_data_before$co_tox_lab_name == 1,]) 
length(unique(lab_data_before$studyid[is.na(lab_data_before$co_tox_lab_normal_low) & lab_data_before$co_tox_lab_name == 1]))
# 105 unique patients had NA for neutropenia
lab_data_before %>% filter(co_tox_lab_name == 1) %>%
  group_by(time_frame, co_tox_lab_name) %>%
  summarise(NAs = sum(is.na(co_tox_lab_normal_low)), 
            NAs_patient = length(unique(studyid[is.na(co_tox_lab_normal_low)])),
            N = n(),
            N_patient = length(unique(studyid)))
# most of NAs in timeframe=4

lab_kpwa_sum <- lab_data_before %>% filter(co_tox_lab_name ==17 & siteid == 2) %>%
  group_by(time_frame, co_tox_lab_name) %>%
  summarise(NAs = sum(is.na(co_tox_lab_normal_low)), 
            N = n(),
            NAs_patient = length(unique(studyid[is.na(co_tox_lab_normal_low)])),
            N_patient = length(unique(studyid))) %>%
  mutate(Percent_N = NAs/N*100, Percent_patient = NAs_patient/N_patient*100)

# summary of normal low in KPNC
summary(lab_data_before$co_tox_lab_normal_low[lab_data_before$co_tox_lab_name == 1 & lab_data_before$siteid==1 & !is.na(lab_data_before$co_tox_lab_result)])
#     Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#      41    1800    2100    2454    2542  101514       2 
lab_data_before[lab_data_before$co_tox_lab_normal_low<100 & lab_data_before$co_tox_lab_name == 1 & lab_data_before$siteid==1 & !is.na(lab_data_before$co_tox_lab_result), c("studyid", "co_tox_lab_date", "co_tox_lab_result", "co_tox_lab_name", "co_tox_lab_unit", "co_tox_lab_normal_low")]
# summary of normal low in KPWA
summary(lab_data_before$co_tox_lab_normal_low[lab_data_before$co_tox_lab_name == 1 & lab_data_before$siteid==2])
#     Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    240    2000    2000    2004    2000    8360     489 
summary(lab_data_before$co_tox_lab_normal_low[lab_data_before$co_tox_lab_name == 1 & lab_data_before$siteid==2 & lab_data_before$time_frame==4])
length(unique(lab_data_before$studyid[lab_data_before$co_tox_lab_name == 1 & lab_data_before$time_frame==4 & lab_data_before$co_tox_lab_normal_low<=1500 & lab_data_before$siteid==2]))
# 13
length(unique(lab_data_before$studyid[lab_data_before$co_tox_lab_name == 1 & lab_data_before$time_frame==4 & lab_data_before$co_tox_lab_normal_low<2000 & lab_data_before$siteid==2]))
# 15
length(unique(lab_data_before$studyid[lab_data_before$co_tox_lab_name == 1 & lab_data_before$time_frame==4 & lab_data_before$co_tox_lab_normal_low>2000 & lab_data_before$siteid==2]))
# 5 

# change normal low <=1500 to NA
lab_data_before$co_tox_lab_normal_low[lab_data_before$co_tox_lab_name==1 & lab_data_before$co_tox_lab_normal_low <=1500] <- NA

# change normal low >3000 to NA
summary(lab_data_before$co_tox_lab_normal_low[lab_data_before$co_tox_lab_name==1 & lab_data_before$co_tox_lab_normal_low >3000])
lab_data_before[!is.na(lab_data_before$co_tox_lab_normal_low) & lab_data_before$co_tox_lab_normal_low >3000 & lab_data_before$co_tox_lab_name ==1 & lab_data_before$time_frame==4, c("studyid", "co_tox_lab_name", "co_tox_lab_result", "co_tox_lab_normal_low", "co_tox_lab_date")] %>% print(n = Inf, width = Inf)
lab_data_before$co_tox_lab_normal_low[lab_data_before$co_tox_lab_name==1 & lab_data_before$co_tox_lab_normal_low >3000] <- NA

# check missing in neutropenia for KPWA
sum(is.na(lab_data_before$co_tox_lab_result[lab_data_before$co_tox_lab_name == 1 & lab_data_before$siteid==2])) # 0 NAs

# patients with any missingness in normal low
a <- unique(lab_data_before$studyid[is.na(lab_data_before$co_tox_lab_normal_low) & lab_data_before$co_tox_lab_name ==17])
lab_missing <- lab_data_before %>% filter(studyid %in% a & co_tox_lab_name ==1) %>%
  group_by(studyid, time_frame, co_tox_lab_name) %>%
  mutate(count = row_number())
length(unique(lab_missing$studyid)) # 111
lab_missing_sum <- lab_missing %>% group_by(studyid, time_frame, co_tox_lab_name) %>%
  summarise(N= n(), NAs = sum(is.na(co_tox_lab_normal_low))) %>%
  mutate(missing_ind = ifelse(N == NAs, "full", ifelse(NAs == 0, "none", "partial")))
b <- lab_missing_sum %>% group_by(time_frame, co_tox_lab_name, missing_ind) %>%
  summarize(N_missing = n()) %>% 
  filter(missing_ind != "none")

summary(lab_missing$co_tox_lab_result[lab_missing$co_tox_lab_name == 1 & is.na(lab_missing$co_tox_lab_normal_low)])

## Imputate Normal Low
lab_data_before <- lab_data_before %>%
  mutate(normal_low_missing = ifelse(is.na(co_tox_lab_normal_low), 1, 0)) # normal low missingness indicator

lab_data_before_neutropenia <- lab_data_before %>%
  filter(co_tox_lab_name == 1) %>% # only neutropenia
  filter(!is.na(co_tox_lab_result)) # remove NAs in lab result

# patients with partial missingness of normal low in specific time frame
lab_missing_sum %>% filter(missing_ind == "partial") %>% print(n = Inf)
d <- unique(lab_missing_sum$studyid[lab_missing_sum$missing_ind == "partial"])
e <- lab_missing[lab_missing$studyid %in% d & !is.na(lab_missing$co_tox_lab_normal_low) & lab_missing$co_tox_lab_name ==1, c("studyid", "time_frame", "co_tox_lab_name", "co_tox_lab_result", "co_tox_lab_normal_low", "co_tox_lab_date")] %>% print(n = Inf, width = Inf)
f <- e %>% group_by(studyid, time_frame, co_tox_lab_name) %>%
  summarise(min = min(co_tox_lab_normal_low), 
            median = median(co_tox_lab_normal_low), 
            max = max(co_tox_lab_normal_low), n = n()) %>%
  mutate(Ind = ifelse(min == max, 1, 0))
f %>% print(n=Inf)
f[f$Ind==0,] %>% print(n=Inf)
## they had different normal low and partial missing;
lab_missing[lab_missing$studyid %in% c("marked-out id", "marked-out id", "marked-out id", "marked-out id", "marked-out id", "marked-out id", "marked-out id", "marked-out id") & lab_missing$co_tox_lab_name ==1, c("studyid", "time_frame", "co_tox_lab_name", "co_tox_lab_result", "co_tox_lab_normal_low", "co_tox_lab_date")] %>% print(n = Inf, width = Inf)
lab_data_before_neutropenia$co_tox_lab_normal_low[lab_data_before_neutropenia$studyid == "marked-out id" & lab_data_before_neutropenia$co_tox_lab_name == 1 & lab_data_before_neutropenia$time_frame == 3 & is.na(lab_data_before_neutropenia$co_tox_lab_normal_low)] <- 1900
lab_data_before_neutropenia$co_tox_lab_normal_low[lab_data_before_neutropenia$studyid == "marked-out id" & lab_data_before_neutropenia$co_tox_lab_name == 1 & lab_data_before_neutropenia$time_frame == 4 & is.na(lab_data_before_neutropenia$co_tox_lab_normal_low)] <- 2000
lab_data_before_neutropenia$co_tox_lab_normal_low[lab_data_before_neutropenia$studyid == "marked-out id" & lab_data_before_neutropenia$co_tox_lab_name == 1 & lab_data_before_neutropenia$time_frame == 4 & is.na(lab_data_before_neutropenia$co_tox_lab_normal_low)] <- 1900
lab_data_before_neutropenia$co_tox_lab_normal_low[lab_data_before_neutropenia$studyid == "marked-out id" & lab_data_before_neutropenia$co_tox_lab_name == 1 & lab_data_before_neutropenia$time_frame == 2 & is.na(lab_data_before_neutropenia$co_tox_lab_normal_low)] <- 1900
lab_data_before_neutropenia$co_tox_lab_normal_low[lab_data_before_neutropenia$studyid == "marked-out id" & lab_data_before_neutropenia$co_tox_lab_name == 1 & lab_data_before_neutropenia$time_frame == 4 & is.na(lab_data_before_neutropenia$co_tox_lab_normal_low)] <- 2000
lab_data_before_neutropenia$co_tox_lab_normal_low[lab_data_before_neutropenia$studyid == "marked-out id" & lab_data_before_neutropenia$co_tox_lab_name == 1 & lab_data_before_neutropenia$time_frame == 4 & is.na(lab_data_before_neutropenia$co_tox_lab_normal_low)] <- 2000
lab_data_before_neutropenia$co_tox_lab_normal_low[lab_data_before_neutropenia$studyid == "marked-out id" & lab_data_before_neutropenia$co_tox_lab_name == 1 & lab_data_before_neutropenia$time_frame == 4 & is.na(lab_data_before_neutropenia$co_tox_lab_normal_low)] <- 1900
## other patients had same normal low values
g <- d
for (i in 1:length(g)){
  for (j in 1) {
    for (k in 1:4){
      t <- length(lab_data_before_neutropenia$co_tox_lab_normal_low[lab_data_before_neutropenia$studyid == g[i] & lab_data_before_neutropenia$co_tox_lab_name == j & lab_data_before_neutropenia$time_frame == k & is.na(lab_data_before_neutropenia$co_tox_lab_normal_low)])
      s <- length(f$max[f$studyid == g[i] & f$co_tox_lab_name == j & f$time_frame == k])
      if (t != 0 & s != 0) {
        lab_data_before_neutropenia$co_tox_lab_normal_low[lab_data_before_neutropenia$studyid == g[i] & lab_data_before_neutropenia$co_tox_lab_name == j & lab_data_before_neutropenia$time_frame == k & is.na(lab_data_before_neutropenia$co_tox_lab_normal_low)] <- f$max[f$studyid == g[i] & f$co_tox_lab_name == j & f$time_frame == k]
      }
    }
  }
}
# summary of normal low in KPNC
summary(lab_data_before_neutropenia$co_tox_lab_normal_low[lab_data_before_neutropenia$co_tox_lab_name == 1 & lab_data_before_neutropenia$siteid==1])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    1512    1845    2100    2101    2100    3000   59072
# summary of normal low in KPWA
summary(lab_data_before_neutropenia$co_tox_lab_normal_low[lab_data_before_neutropenia$co_tox_lab_name == 1 & lab_data_before_neutropenia$siteid==2])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    1520    2000    2000    1999    2000    2520     283 
# median normal low in KPNC: 2100
lab_data_before_neutropenia$co_tox_lab_normal_low[lab_data_before_neutropenia$siteid==1 & lab_data_before_neutropenia$co_tox_lab_name == 1 & is.na(lab_data_before_neutropenia$co_tox_lab_normal_low)] <- 2100
# median normal low in KPWA: 2000
lab_data_before_neutropenia$co_tox_lab_normal_low[lab_data_before_neutropenia$siteid==2 & lab_data_before_neutropenia$co_tox_lab_name == 1 & is.na(lab_data_before_neutropenia$co_tox_lab_normal_low)] <- 2000
sum(is.na(lab_data_before_neutropenia$co_tox_lab_normal_low)) # 0


## Flag Neutropenia 
lab_data_before_neutropenia <- lab_data_before_neutropenia %>%
  mutate(
    neutropenia_G3_lab = ifelse(co_tox_lab_name == 1 & co_tox_lab_result < 1000, 1, 0),
    neutropenia_G2_lab = ifelse(co_tox_lab_name == 1 & co_tox_lab_result < 1500, 1, 0),
    neutropenia_G1_lab = ifelse(co_tox_lab_name == 1 & co_tox_lab_result < co_tox_lab_normal_low, 1, 0)
  )

## Summarize Neutropenia Toxicity
icd_data_before_sum_neutropenia <- icd_data_before %>% 
  select(studyid, time_frame, neutropenia_flag_icd) %>%
  group_by(studyid, time_frame) %>%
  summarise(neutropenia_icd_n = sum(neutropenia_flag_icd))
icd_data_before_sum_neutropenia <- icd_data_before_sum_neutropenia %>% mutate(neutropenia_flag_icd = ifelse(neutropenia_icd_n == 0, 0, 1))
icd_data_before_sum_neutropenia$siteid <- ifelse(icd_data_before_sum_neutropenia$studyid>"marked-out id", 2, 1)
table(icd_data_before_sum_neutropenia$neutropenia_icd_n, useNA = "ifany") # check
table(icd_data_before_sum_neutropenia$neutropenia_flag_icd, useNA = "ifany") # check
icd_data_before_sum_neutropenia[icd_data_before_sum_neutropenia$neutropenia_icd_n >= 1,] # check

lab_data_before_sum_neutropenia <- lab_data_before_neutropenia %>%
  select(studyid, time_frame, neutropenia_G3_lab, neutropenia_G2_lab, neutropenia_G1_lab) %>%
  group_by(studyid, time_frame) %>%
  summarise(
    neutropenia_G3_n = sum(neutropenia_G3_lab),
    neutropenia_G2_n = sum(neutropenia_G2_lab),
    neutropenia_G1_n = sum(neutropenia_G1_lab)
  ) %>%
  mutate(
    neutropenia_G3_lab = ifelse(neutropenia_G3_n == 0, 0, 1),
    neutropenia_G2_lab = ifelse(neutropenia_G2_n == 0, 0, 1),
    neutropenia_G1_lab = ifelse(neutropenia_G1_n == 0, 0, 1)
  )
lab_data_before_sum_neutropenia$siteid <- ifelse(lab_data_before_sum_neutropenia$studyid>"marked-out id", 2, 1)
table(lab_data_before_sum_neutropenia$neutropenia_G3_n, useNA = "ifany") # check
table(lab_data_before_sum_neutropenia$neutropenia_G2_n, useNA = "ifany") # check
table(lab_data_before_sum_neutropenia$neutropenia_G1_n, useNA = "ifany") # check
table(lab_data_before_sum_neutropenia$neutropenia_G3_lab, useNA = "ifany")
table(lab_data_before_sum_neutropenia$neutropenia_G2_lab, useNA = "ifany")
table(lab_data_before_sum_neutropenia$neutropenia_G1_lab, useNA = "ifany")


chart_before_sum_neutropenia <- chart_before %>%
  mutate(neutropenia_d = ifelse(toxicityname == "Neutropenia/Leukopenia", 1, 0)) %>%
  group_by(studyid) %>%
  summarise(neutropenia_n = sum(neutropenia_d))
chart_before_sum_neutropenia <- chart_before_sum_neutropenia %>% mutate(neutropenia_tox = ifelse(neutropenia_n == 0, 0, 1))
table(chart_before_sum_neutropenia$neutropenia_tox, useNA = "ifany") # 1=248, 0=621
table(chart_before_sum_neutropenia$neutropenia_n, useNA = "ifany")


## Merge ICD and Lab 
lab_data_before_sum_neutropenia$missing_lab <- rep(0)
icd_data_before_sum_neutropenia$missing_icd <- rep(0)
lab_id_1 <- lab_data_before_sum_neutropenia$studyid[lab_data_before_sum_neutropenia$time_frame==1]
lab_id_2 <- lab_data_before_sum_neutropenia$studyid[lab_data_before_sum_neutropenia$time_frame==2]
lab_id_3 <- lab_data_before_sum_neutropenia$studyid[lab_data_before_sum_neutropenia$time_frame==3]
lab_id_4 <- lab_data_before_sum_neutropenia$studyid[lab_data_before_sum_neutropenia$time_frame==4]
icd_id_1 <- icd_data_before_sum_neutropenia$studyid[icd_data_before_sum_neutropenia$time_frame==1]
icd_id_2 <- icd_data_before_sum_neutropenia$studyid[icd_data_before_sum_neutropenia$time_frame==2]
icd_id_3 <- icd_data_before_sum_neutropenia$studyid[icd_data_before_sum_neutropenia$time_frame==3]
icd_id_4 <- icd_data_before_sum_neutropenia$studyid[icd_data_before_sum_neutropenia$time_frame==4]
sum(!(lab_id_1 %in% icd_id_1)) 
sum(!(lab_id_2 %in% icd_id_2))
sum(!(lab_id_3 %in% icd_id_3)) 
sum(!(lab_id_4 %in% icd_id_4)) 
sum(!(icd_id_1 %in% lab_id_1)) 
sum(!(icd_id_2 %in% lab_id_2)) 
sum(!(icd_id_3 %in% lab_id_3)) 
sum(!(icd_id_4 %in% lab_id_4)) 

a <- data.frame(
  studyid = c(icd_id_1[!(icd_id_1 %in% lab_id_1)],
              icd_id_2[!(icd_id_2 %in% lab_id_2)],
              icd_id_3[!(icd_id_3 %in% lab_id_3)],
              icd_id_4[!(icd_id_4 %in% lab_id_4)]),
  time_frame = c(rep(1, times = 2562),
                 rep(2, times = 3076),
                 rep(3, times = 188),
                 rep(4, times = 140)),
  missing_lab = rep(1), 
  neutropenia_G3_lab = rep(0),
  neutropenia_G2_lab = rep(0),
  neutropenia_G1_lab = rep(0)
)

lab_data_before_sum_neutropenia <- rbind(lab_data_before_sum_neutropenia, a)
lab_data_before_sum_neutropenia$siteid <- ifelse(lab_data_before_sum_neutropenia$studyid>"marked-out id", 2, 1)
anyNA(lab_data_before_sum_neutropenia$studyid)
anyNA(lab_data_before_sum_neutropenia$time_frame)
anyNA(lab_data_before_sum_neutropenia$neutropenia_G3_lab)

b <- data.frame(
  studyid = c(lab_id_1[!(lab_id_1 %in% icd_id_1)],
              lab_id_2[!(lab_id_2 %in% icd_id_2)],
              lab_id_3[!(lab_id_3 %in% icd_id_3)],
              lab_id_4[!(lab_id_4 %in% icd_id_4)]),
  time_frame = c(rep(1, times = 1628),
                 rep(2, times = 1701),
                 rep(3, times = 7605),
                 rep(4, times = 5765)),
  missing_icd = rep(1),
  neutropenia_flag_icd = rep(0)
)
icd_data_before_sum_neutropenia <- rbind(icd_data_before_sum_neutropenia, b)
icd_data_before_sum_neutropenia$siteid <- ifelse(icd_data_before_sum_neutropenia$studyid>"marked-out id", 2, 1)
anyNA(icd_data_before_sum_neutropenia$studyid)
anyNA(icd_data_before_sum_neutropenia$neutropenia_flag_icd)

icd_data_before_sum_neutropenia <- icd_data_before_sum_neutropenia %>% arrange(studyid, time_frame)
lab_data_before_sum_neutropenia <- lab_data_before_sum_neutropenia %>% arrange(studyid, time_frame)
length(unique(lab_data_before_sum_neutropenia$studyid))
length(unique(icd_data_before_sum_neutropenia$studyid))
identical(icd_data_before_sum_neutropenia$studyid, lab_data_before_sum_neutropenia$studyid) # T
identical(icd_data_before_sum_neutropenia$time_frame, lab_data_before_sum_neutropenia$time_frame) # T

# merge icd and lab
before_sum_neutropenia <- icd_data_before_sum_neutropenia %>%
  inner_join(lab_data_before_sum_neutropenia[, names(lab_data_before_sum_neutropenia)[c(1:8, 10)]], by = join_by(studyid, time_frame))
adju_cohort <- readRDS(file = paste0(datapath, "adju_cohort.rds"))
before_sum_neutropenia <- before_sum_neutropenia[, names(before_sum_neutropenia)[c(1:4, 6:13)]] %>%
  right_join(adju_cohort, by = join_by(studyid, time_frame))
glimpse(before_sum_neutropenia)
table(before_sum_neutropenia$neutropenia_flag_icd, useNA = "ifany")
before_sum_neutropenia$neutropenia_flag_icd[is.na(before_sum_neutropenia$neutropenia_flag_icd)] <- 0
table(before_sum_neutropenia$neutropenia_G3_lab, useNA = "ifany")
before_sum_neutropenia$neutropenia_G3_lab[is.na(before_sum_neutropenia$neutropenia_G3_lab)] <- 0
table(before_sum_neutropenia$neutropenia_G2_lab, useNA = "ifany")
before_sum_neutropenia$neutropenia_G2_lab[is.na(before_sum_neutropenia$neutropenia_G2_lab)] <- 0
table(before_sum_neutropenia$neutropenia_G1_lab, useNA = "ifany")
before_sum_neutropenia$neutropenia_G1_lab[is.na(before_sum_neutropenia$neutropenia_G1_lab)] <- 0
table(before_sum_neutropenia$missing_icd, useNA = "ifany")
before_sum_neutropenia$missing_icd[is.na(before_sum_neutropenia$missing_icd)] <- 1
table(before_sum_neutropenia$missing_lab, useNA = "ifany")
before_sum_neutropenia$missing_lab[is.na(before_sum_neutropenia$missing_lab)] <- 1

# create hep_tox_indicator
before_sum_neutropenia <- before_sum_neutropenia %>%
  mutate(neutropenia_tox_indicator_G3 = ifelse(neutropenia_flag_icd == 1 & neutropenia_G3_lab == 1, "Both",
                                               ifelse(neutropenia_flag_icd == 1 & neutropenia_G3_lab == 0, "Only ICD",
                                                      ifelse(neutropenia_flag_icd == 0 & neutropenia_G3_lab == 1, "Only Grade 3", "Neither"))),
         neutropenia_tox_indicator_G2 = ifelse(neutropenia_flag_icd == 1 & neutropenia_G2_lab == 1, "Both",
                                               ifelse(neutropenia_flag_icd == 1 & neutropenia_G2_lab == 0, "Only ICD",
                                                      ifelse(neutropenia_flag_icd == 0 & neutropenia_G2_lab == 1, "Only Grade 2", "Neither"))))

# subset before_sum_neutropenia and export
adjuvant_neutropenia <- before_sum_neutropenia[, c("studyid", "time_frame", "neutropenia_flag_icd", "neutropenia_G3_lab", "neutropenia_G2_lab", "neutropenia_G1_lab",
                                                   "missing_icd", "missing_lab", "dx_dt", "siteid", "first_chemo_dt", "last_chemo_dt",
                                                   "surgery_d", "surgery_dt", "group", "neutropenia_tox_indicator_G3", "neutropenia_tox_indicator_G2")]
write.csv(adjuvant_neutropenia, file = paste0(datapath, "adjuvant_neutropenia.csv"),
          row.names = F)


## Merge IOCD, Lab, and Chart 
tox_sum_neutropenia <- chart_before_sum_neutropenia[, c("studyid", "neutropenia_tox")] %>%
  right_join(before_sum_neutropenia[before_sum_neutropenia$time_frame == 4 & before_sum_neutropenia$siteid==2,], by = join_by(studyid))
glimpse(tox_sum_neutropenia) 
tox_sum_neutropenia$neutropenia_tox[is.na(tox_sum_neutropenia$neutropenia_tox)] <- 0


# export data
# saveRDS(icd_data_before_sum_neutropenia, file = paste0(datapath, "icd_before_sum_neutropenia.rds"))
# saveRDS(lab_data_before_neutropenia, file = paste0(datapath, "lab_before_neutropenia.rds"))
# saveRDS(lab_data_before_sum_neutropenia, file = paste0(datapath, "lab_before_sum_neutropenia.rds"))
# saveRDS(before_sum_neutropenia, file = paste0(datapath, "before_sum_neutropenia.rds"))
# saveRDS(tox_sum_neutropenia, file = paste0(datapath, "tox_sum_neutropenia.rds"))
# saveRDS(chart_before_sum_neutropenia, file = paste0(datapath, "chart_before_sum_neutropenia.rds"))
# before_sum_neutropenia <- readRDS(file = paste0(datapath, "before_sum_neutropenia.rds"))




## Thrombocytopenia ##
## Load Dataset 
icd_data <- readRDS(file = paste0(datapath, "icd.rds"))
lab_data <- readRDS(file = paste0(datapath, "lab.rds"))
chart_data <- readRDS(file = paste0(datapath, "chart.rds"))
obcd_group <- readRDS(file = paste0(datapath, "obcd.rds"))

## Check and Clean
length(unique(icd_data$studyid)) 
length(unique(icd_data$studyid[icd_data$any_chemo==1])) 
anyNA(icd_data$co_tox_dx_code_date) # F
table(icd_data$thrombocytopenia_flag_icd, useNA = "ifany") 
length(unique(lab_data$studyid)) 
length(unique(lab_data$studyid[lab_data$any_chemo==1])) 
anyNA(lab_data$co_tox_lab_date) #F
anyNA(lab_data$co_tox_lab_name) #F
sum(is.na(lab_data$co_tox_lab_result[lab_data$co_tox_lab_name == 16 ])) # 104 NA

table(lab_data$co_tox_lab_unit[lab_data$co_tox_lab_name==16], useNA = "ifany") # k/ul, unknown, NA
summary(lab_data$co_tox_lab_result[lab_data$co_tox_lab_name == 16 ]) 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.0   194.0   243.0   249.8   296.0  1823.0     104
summary(lab_data$co_tox_lab_result[lab_data$co_tox_lab_name == 16 & lab_data$co_tox_lab_unit == 3]) 
summary(lab_data$co_tox_lab_result[lab_data$co_tox_lab_name == 16 & lab_data$co_tox_lab_unit == 11]) # similar to 11
# 1 K/uL = 1000 uL = 1000 mm^3
lab_data$co_tox_lab_result[lab_data$co_tox_lab_name == 16 & lab_data$co_tox_lab_unit == 11] 
# should be k/uL
lab_data$co_tox_lab_result[lab_data$co_tox_lab_name == 16 & is.na(lab_data$co_tox_lab_unit)]
# should be k/uL

sum(lab_data$co_tox_lab_result[lab_data$co_tox_lab_name == 16] == 0 & !is.na(lab_data$co_tox_lab_result[lab_data$co_tox_lab_name == 16])) # 7

# Normal Low
sum(is.na(lab_data$co_tox_lab_normal_low[lab_data$co_tox_lab_name == 16 ])) # 3900 NAs
sum(is.na(lab_data$co_tox_lab_normal_low[lab_data$co_tox_lab_name == 16 & lab_data$siteid == 2])) # 3899 NAs 
summary(lab_data$co_tox_lab_normal_low[lab_data$co_tox_lab_name == 16 & lab_data$co_tox_lab_unit == 3]) 
#     Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#     140     140     140     140     140     260    3901 
lab_data[!is.na(lab_data$co_tox_lab_normal_low) & lab_data$co_tox_lab_normal_low >200 & lab_data$co_tox_lab_name ==16, c("studyid", "co_tox_lab_name", "co_tox_lab_result", "co_tox_lab_normal_low", "co_tox_lab_date")] %>% print(n = Inf, width = Inf)
# studyid co_tox_lab_name co_tox_lab_result co_tox_lab_normal_low co_tox_lab_date
# 206058              16               914                   260 1997-01-16 
summary(lab_data$co_tox_lab_normal_low[lab_data$co_tox_lab_name == 16 & lab_data$co_tox_lab_unit == 11]) 
#     Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#     154     154     154     154     154     154      23
summary(lab_data$co_tox_lab_normal_low[lab_data$co_tox_lab_name == 16 & is.na(lab_data$co_tox_lab_unit)]) # 140

length(unique(lab_data$studyid[is.na(lab_data$co_tox_lab_normal_low) & lab_data$co_tox_lab_name == 16]))
# 1487 unique patients had NA


## deal with 0 as lab results
a <- lab_data[lab_data$co_tox_lab_name == 16 & lab_data$co_tox_lab_result == 0,]
length(unique(a$studyid)) # 6 unique patients
table(a$co_tox_lab_unit) # 7 K/uL
# for K/uL, 0 may be caused by rounding


## converse K/uL and unknown to uL
lab_data$co_tox_lab_result[lab_data$co_tox_lab_name == 16 & !is.na(lab_data$co_tox_lab_result)] <- lab_data$co_tox_lab_result[lab_data$co_tox_lab_name == 16 & !is.na(lab_data$co_tox_lab_result)]*1000
lab_data$co_tox_lab_normal_low[lab_data$co_tox_lab_name == 16 & !is.na(lab_data$co_tox_lab_normal_low)] <- lab_data$co_tox_lab_normal_low[lab_data$co_tox_lab_name == 16 & !is.na(lab_data$co_tox_lab_normal_low)]*1000


# check chart data
a <- chart_data[is.na(chart_data$toxicitydate), c("studyid", "toxicityid","toxicityname", "toxicitydate", "dx_dt", "surgery_d", "surgery_dt", "first_chemo_dt")] %>% print(n=Inf, width = Inf)
# 1 Thrombocytopenia
# exclude those records
b <- chart_data[!is.na(chart_data$toxicitydate) & (chart_data$toxicitydate < chart_data$first_chemo_dt), c("studyid", "toxicityid", "toxicityname", "toxicitydate", "dx_dt", "surgery_d", "surgery_dt", "first_chemo_dt", "last_chemo_dt")] %>% print(n=Inf, width = Inf)
# no
c <- chart_data[!is.na(chart_data$toxicitydate) & (chart_data$toxicitydate == chart_data$first_chemo_dt), c("studyid", "toxicityid", "toxicityname", "toxicitydate", "dx_dt", "surgery_d", "surgery_dt", "first_chemo_dt")] %>% print(n=Inf, width = Inf)
# no

chart_data <- chart_data[!(chart_data$toxicityid %in% c("marked-out id")),]

## Adding Time Frame 
# time frame for surgery before chemo
# 1: <= date of cancer diagnosis --> pre-existing conditions (1 year)
# 2: <= date of surgery
# 3: <= date of first chemo
# 4: <= end of chemo
icd_data_before <- icd_data %>%
  filter(group == "Surgery before chemo") %>%
  mutate(time_frame = ifelse(co_tox_dx_code_date <= dx_dt, 1,
                             ifelse(co_tox_dx_code_date <= surgery_dt, 2,
                                    ifelse(co_tox_dx_code_date <= first_chemo_dt, 3,
                                           ifelse(co_tox_dx_code_date <= last_chemo_dt, 4, NA))))) %>%
  filter(!is.na(time_frame)) # remove NAs in time frame

table(icd_data_before$time_frame, useNA = "ifany")
table(icd_data_before$thrombocytopenia_flag_icd, useNA = "ifany") 

lab_data_before <- lab_data %>%
  filter(group == "Surgery before chemo") %>%
  mutate(time_frame = ifelse(co_tox_lab_date <= dx_dt, 1,
                             ifelse(co_tox_lab_date <= surgery_dt, 2,
                                    ifelse(co_tox_lab_date <= first_chemo_dt, 3,
                                           ifelse(co_tox_lab_date <= last_chemo_dt, 4, NA))))) %>%
  filter(!is.na(time_frame)) # remove NAs in time frame

table(lab_data_before$time_frame, useNA = "ifany")
table(lab_data_before$co_tox_lab_name, useNA = "ifany") 


chart_before <- chart_data %>% filter(surgery_d == 1) %>%
  filter(!is.na(surgery_dt)) %>%
  filter(surgery_dt < first_chemo_dt) %>%
  filter(toxicitydate > first_chemo_dt & toxicitydate <= last_chemo_dt) # limit to last timeframe

## che
# check normal low in lab_data_before
sum(is.na(lab_data_before$co_tox_lab_normal_low[lab_data_before$co_tox_lab_name == 16 ])) 
sum(is.na(lab_data_before$co_tox_lab_normal_low[lab_data_before$co_tox_lab_name == 16 & lab_data_before$siteid == 2])) 
nrow(lab_data_before[lab_data_before$co_tox_lab_name == 16,]) 
length(unique(lab_data_before$studyid[is.na(lab_data_before$co_tox_lab_normal_low) & lab_data_before$co_tox_lab_name == 16]))
# 107 unique patients had NA
# all NAs in KPWA
lab_data_before %>% filter(co_tox_lab_name == 16) %>%
  group_by(time_frame, co_tox_lab_name) %>%
  summarise(NAs = sum(is.na(co_tox_lab_normal_low)), 
            NAs_patient = length(unique(studyid[is.na(co_tox_lab_normal_low)])),
            N = n(),
            N_patient = length(unique(studyid)))
# most of NAs in timeframe=4

lab_kpwa_sum <- lab_data_before %>% filter(co_tox_lab_name ==16 & siteid == 2) %>%
  group_by(time_frame, co_tox_lab_name) %>%
  summarise(NAs = sum(is.na(co_tox_lab_normal_low)), 
            N = n(),
            NAs_patient = length(unique(studyid[is.na(co_tox_lab_normal_low)])),
            N_patient = length(unique(studyid))) %>%
  mutate(Percent_N = NAs/N*100, Percent_patient = NAs_patient/N_patient*100)

# summary of normal low in KPWA
summary(lab_data_before$co_tox_lab_normal_low[lab_data_before$co_tox_lab_name == 16 & lab_data_before$siteid==2])
#     Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  140000  140000  140000  140097  140000  150000     516

# check missing in thrombo for KPWA
sum(is.na(lab_data_before$co_tox_lab_result[lab_data_before$co_tox_lab_name == 16 & lab_data_before$siteid==2])) # 3 NAs

# patients with any missingness in normal low
a <- unique(lab_data_before$studyid[is.na(lab_data_before$co_tox_lab_normal_low) & lab_data_before$co_tox_lab_name ==16])
lab_missing <- lab_data_before %>% filter(studyid %in% a & co_tox_lab_name ==16) %>%
  group_by(studyid, time_frame, co_tox_lab_name) %>%
  mutate(count = row_number())
length(unique(lab_missing$studyid))
lab_missing_sum <- lab_missing %>% group_by(studyid, time_frame, co_tox_lab_name) %>%
  summarise(N= n(), NAs = sum(is.na(co_tox_lab_normal_low))) %>%
  mutate(missing_ind = ifelse(N == NAs, "full", ifelse(NAs == 0, "none", "partial")))
b <- lab_missing_sum %>% group_by(time_frame, co_tox_lab_name, missing_ind) %>%
  summarize(N_missing = n()) %>% 
  filter(missing_ind != "none")

summary(lab_missing$co_tox_lab_result[lab_missing$co_tox_lab_name == 16 & is.na(lab_missing$co_tox_lab_normal_low)])

## Imputate Normal Low
lab_data_before <- lab_data_before %>%
  mutate(normal_low_missing = ifelse(is.na(co_tox_lab_normal_low), 1, 0)) # normal low missingness indicator

lab_data_before_thrombocytopenia <- lab_data_before %>%
  filter(co_tox_lab_name == 16) %>% # only thrombocytopenia
  filter(!is.na(co_tox_lab_result)) # remove NAs in lab result

# patients with partial missingness of normal low in specific time frame
lab_missing_sum %>% filter(missing_ind == "partial") %>% print(n = Inf)
d <- unique(lab_missing_sum$studyid[lab_missing_sum$missing_ind == "partial"])
e <- lab_missing[lab_missing$studyid %in% d & !is.na(lab_missing$co_tox_lab_normal_low) & lab_missing$co_tox_lab_name ==16, c("studyid", "time_frame", "co_tox_lab_name", "co_tox_lab_result", "co_tox_lab_normal_low", "co_tox_lab_date")] %>% print(n = Inf, width = Inf)
f <- e %>% group_by(studyid, time_frame, co_tox_lab_name) %>%
  summarise(min = min(co_tox_lab_normal_low), 
            median = median(co_tox_lab_normal_low), 
            max = max(co_tox_lab_normal_low), n = n()) %>%
  mutate(Ind = ifelse(min == max, 1, 0))
f %>% print(n=Inf)
f[f$Ind==0,] %>% print(n=Inf)
## they had different normal low and partial missing;
lab_missing[lab_missing$studyid %in% c("marked-out id", "marked-out id", "marked-out id", "marked-out id", "marked-out id") & lab_missing$co_tox_lab_name ==16, c("studyid", "time_frame", "co_tox_lab_name", "co_tox_lab_result", "co_tox_lab_normal_low", "co_tox_lab_date")] %>% print(n = Inf, width = Inf)
lab_data_before_thrombocytopenia$co_tox_lab_normal_low[lab_data_before_thrombocytopenia$studyid == "marked-out id" & lab_data_before_thrombocytopenia$co_tox_lab_name == 16 & lab_data_before_thrombocytopenia$time_frame == 4 & is.na(lab_data_before_thrombocytopenia$co_tox_lab_normal_low)] <- 140000
lab_data_before_thrombocytopenia$co_tox_lab_normal_low[lab_data_before_thrombocytopenia$studyid == "marked-out id" & lab_data_before_thrombocytopenia$co_tox_lab_name == 16 & lab_data_before_thrombocytopenia$time_frame == 1 & is.na(lab_data_before_thrombocytopenia$co_tox_lab_normal_low)] <- 150000
## other patients had same normal low values
g <- d
for (i in 1:length(g)){
  for (j in 16) {
    for (k in 1:4){
      t <- length(lab_data_before_thrombocytopenia$co_tox_lab_normal_low[lab_data_before_thrombocytopenia$studyid == g[i] & lab_data_before_thrombocytopenia$co_tox_lab_name == j & lab_data_before_thrombocytopenia$time_frame == k & is.na(lab_data_before_thrombocytopenia$co_tox_lab_normal_low)])
      s <- length(f$max[f$studyid == g[i] & f$co_tox_lab_name == j & f$time_frame == k])
      if (t != 0 & s != 0) {
        lab_data_before_thrombocytopenia$co_tox_lab_normal_low[lab_data_before_thrombocytopenia$studyid == g[i] & lab_data_before_thrombocytopenia$co_tox_lab_name == j & lab_data_before_thrombocytopenia$time_frame == k & is.na(lab_data_before_thrombocytopenia$co_tox_lab_normal_low)] <- f$max[f$studyid == g[i] & f$co_tox_lab_name == j & f$time_frame == k]
      }
    }
  }
}
lab_data_before_thrombocytopenia[lab_data_before_thrombocytopenia$studyid=="marked-out id" & lab_data_before_thrombocytopenia$time_frame == 4 & lab_data_before_thrombocytopenia$co_tox_lab_name == 16, "co_tox_lab_normal_low"] %>% print(n = Inf)
# summary of normal low in KPNC
summary(lab_data_before_thrombocytopenia$co_tox_lab_normal_low[lab_data_before_thrombocytopenia$co_tox_lab_name == 16 & lab_data_before_thrombocytopenia$siteid==1]) # 140000
# summary of normal low in KPWA
summary(lab_data_before_thrombocytopenia$co_tox_lab_normal_low[lab_data_before_thrombocytopenia$co_tox_lab_name == 16 & lab_data_before_thrombocytopenia$siteid==2])
#     Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  140000  140000  140000  140098  140000  150000     274 
# median normal low in KPWA: 140000
lab_data_before_thrombocytopenia$co_tox_lab_normal_low[lab_data_before_thrombocytopenia$co_tox_lab_name == 16 & is.na(lab_data_before_thrombocytopenia$co_tox_lab_normal_low)] <- 140000
sum(is.na(lab_data_before_thrombocytopenia$co_tox_lab_normal_low)) # 0
summary(lab_data_before_thrombocytopenia$co_tox_lab_normal_low[lab_data_before_thrombocytopenia$co_tox_lab_name == 16 & lab_data_before_thrombocytopenia$siteid==2 & lab_data_before_thrombocytopenia$time_frame==4])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 140000  140000  140000  140071  140000  150000

## Flag thrombocytopenia
lab_data_before_thrombocytopenia <- lab_data_before_thrombocytopenia %>%
  mutate(
    thrombocytopenia_G3_lab = ifelse(co_tox_lab_name == 16 & co_tox_lab_result < 50000, 1, 0),
    thrombocytopenia_G2_lab = ifelse(co_tox_lab_name == 16 & co_tox_lab_result < 75000, 1, 0),
    thrombocytopenia_G1_lab = ifelse(co_tox_lab_name == 16 & co_tox_lab_result < co_tox_lab_normal_low, 1, 0)
  )

## Summarize thrombocytopenia 
icd_data_before_sum_thrombocytopenia <- icd_data_before %>% 
  select(studyid, time_frame, thrombocytopenia_flag_icd) %>%
  group_by(studyid, time_frame) %>%
  summarise(thrombocytopenia_icd_n = sum(thrombocytopenia_flag_icd))
icd_data_before_sum_thrombocytopenia <- icd_data_before_sum_thrombocytopenia %>% mutate(thrombocytopenia_flag_icd = ifelse(thrombocytopenia_icd_n == 0, 0, 1))
icd_data_before_sum_thrombocytopenia$siteid <- ifelse(icd_data_before_sum_thrombocytopenia$studyid>"marked-out id", 2, 1)
table(icd_data_before_sum_thrombocytopenia$thrombocytopenia_icd_n, useNA = "ifany") # check
table(icd_data_before_sum_thrombocytopenia$thrombocytopenia_flag_icd, useNA = "ifany") # check
icd_data_before_sum_thrombocytopenia[icd_data_before_sum_thrombocytopenia$thrombocytopenia_icd_n >= 1,] # check

lab_data_before_sum_thrombocytopenia <- lab_data_before_thrombocytopenia %>%
  select(studyid, time_frame, thrombocytopenia_G3_lab, thrombocytopenia_G2_lab, thrombocytopenia_G1_lab) %>%
  group_by(studyid, time_frame) %>%
  summarise(
    thrombocytopenia_G3_n = sum(thrombocytopenia_G3_lab),
    thrombocytopenia_G2_n = sum(thrombocytopenia_G2_lab),
    thrombocytopenia_G1_n = sum(thrombocytopenia_G1_lab)
  ) %>%
  mutate(
    thrombocytopenia_G3_lab = ifelse(thrombocytopenia_G3_n == 0, 0, 1),
    thrombocytopenia_G2_lab = ifelse(thrombocytopenia_G2_n == 0, 0, 1),
    thrombocytopenia_G1_lab = ifelse(thrombocytopenia_G1_n == 0, 0, 1)
  )
lab_data_before_sum_thrombocytopenia$siteid <- ifelse(lab_data_before_sum_thrombocytopenia$studyid>"marked-out id", 2, 1)
table(lab_data_before_sum_thrombocytopenia$thrombocytopenia_G3_n, useNA = "ifany") # check
table(lab_data_before_sum_thrombocytopenia$thrombocytopenia_G2_n, useNA = "ifany") # check
table(lab_data_before_sum_thrombocytopenia$thrombocytopenia_G1_n, useNA = "ifany") # check
table(lab_data_before_sum_thrombocytopenia$thrombocytopenia_G3_lab, useNA = "ifany")
table(lab_data_before_sum_thrombocytopenia$thrombocytopenia_G2_lab, useNA = "ifany")
table(lab_data_before_sum_thrombocytopenia$thrombocytopenia_G1_lab, useNA = "ifany")


chart_before_sum_thrombocytopenia <- chart_before %>%
  mutate(thrombocytopenia_d = ifelse(toxicityname == "Thrombocytopenia", 1, 0)) %>%
  group_by(studyid) %>%
  summarise(thrombocytopenia_n = sum(thrombocytopenia_d))
chart_before_sum_thrombocytopenia <- chart_before_sum_thrombocytopenia %>% mutate(thrombocytopenia_tox = ifelse(thrombocytopenia_n == 0, 0, 1))
table(chart_before_sum_thrombocytopenia$thrombocytopenia_tox, useNA = "ifany") # 1=21, 0=848
table(chart_before_sum_thrombocytopenia$thrombocytopenia_n, useNA = "ifany")


## Merge ICD and Lab 
lab_data_before_sum_thrombocytopenia$missing_lab <- rep(0)
icd_data_before_sum_thrombocytopenia$missing_icd <- rep(0)
lab_id_1 <- lab_data_before_sum_thrombocytopenia$studyid[lab_data_before_sum_thrombocytopenia$time_frame==1]
lab_id_2 <- lab_data_before_sum_thrombocytopenia$studyid[lab_data_before_sum_thrombocytopenia$time_frame==2]
lab_id_3 <- lab_data_before_sum_thrombocytopenia$studyid[lab_data_before_sum_thrombocytopenia$time_frame==3]
lab_id_4 <- lab_data_before_sum_thrombocytopenia$studyid[lab_data_before_sum_thrombocytopenia$time_frame==4]
icd_id_1 <- icd_data_before_sum_thrombocytopenia$studyid[icd_data_before_sum_thrombocytopenia$time_frame==1]
icd_id_2 <- icd_data_before_sum_thrombocytopenia$studyid[icd_data_before_sum_thrombocytopenia$time_frame==2]
icd_id_3 <- icd_data_before_sum_thrombocytopenia$studyid[icd_data_before_sum_thrombocytopenia$time_frame==3]
icd_id_4 <- icd_data_before_sum_thrombocytopenia$studyid[icd_data_before_sum_thrombocytopenia$time_frame==4]
sum(!(lab_id_1 %in% icd_id_1)) 
sum(!(lab_id_2 %in% icd_id_2)) 
sum(!(lab_id_3 %in% icd_id_3)) 
sum(!(lab_id_4 %in% icd_id_4)) 
sum(!(icd_id_1 %in% lab_id_1)) 
sum(!(icd_id_2 %in% lab_id_2)) 
sum(!(icd_id_3 %in% lab_id_3)) 
sum(!(icd_id_4 %in% lab_id_4))

a <- data.frame(
  studyid = c(icd_id_1[!(icd_id_1 %in% lab_id_1)],
              icd_id_2[!(icd_id_2 %in% lab_id_2)],
              icd_id_3[!(icd_id_3 %in% lab_id_3)],
              icd_id_4[!(icd_id_4 %in% lab_id_4)]),
  time_frame = c(rep(1, times = 1439),
                 rep(2, times = 1282),
                 rep(3, times = 99),
                 rep(4, times = 31)),
  missing_lab = rep(1), 
  thrombocytopenia_G3_lab = rep(0),
  thrombocytopenia_G2_lab = rep(0),
  thrombocytopenia_G1_lab = rep(0)
)

lab_data_before_sum_thrombocytopenia <- rbind(lab_data_before_sum_thrombocytopenia, a)
lab_data_before_sum_thrombocytopenia$siteid <- ifelse(lab_data_before_sum_thrombocytopenia$studyid>"marked-out id", 2, 1)
anyNA(lab_data_before_sum_thrombocytopenia$studyid)
anyNA(lab_data_before_sum_thrombocytopenia$time_frame)
anyNA(lab_data_before_sum_thrombocytopenia$thrombocytopenia_G3_lab)

b <- data.frame(
  studyid = c(lab_id_1[!(lab_id_1 %in% icd_id_1)],
              lab_id_2[!(lab_id_2 %in% icd_id_2)],
              lab_id_3[!(lab_id_3 %in% icd_id_3)],
              lab_id_4[!(lab_id_4 %in% icd_id_4)]),
  time_frame = c(rep(1, times = 2996),
                 rep(2, times = 4317),
                 rep(3, times = 7824),
                 rep(4, times = 5885)),
  missing_icd = rep(1),
  thrombocytopenia_flag_icd = rep(0)
)
icd_data_before_sum_thrombocytopenia <- rbind(icd_data_before_sum_thrombocytopenia, b)
icd_data_before_sum_thrombocytopenia$siteid <- ifelse(icd_data_before_sum_thrombocytopenia$studyid>"marked-out id", 2, 1)
anyNA(icd_data_before_sum_thrombocytopenia$studyid)
anyNA(icd_data_before_sum_thrombocytopenia$thrombocytopenia_flag_icd)

icd_data_before_sum_thrombocytopenia <- icd_data_before_sum_thrombocytopenia %>% arrange(studyid, time_frame)
lab_data_before_sum_thrombocytopenia <- lab_data_before_sum_thrombocytopenia %>% arrange(studyid, time_frame)
length(unique(lab_data_before_sum_thrombocytopenia$studyid))
length(unique(icd_data_before_sum_thrombocytopenia$studyid))
identical(icd_data_before_sum_thrombocytopenia$studyid, lab_data_before_sum_thrombocytopenia$studyid) # T
identical(icd_data_before_sum_thrombocytopenia$time_frame, lab_data_before_sum_thrombocytopenia$time_frame) # T

# merge icd and lab
before_sum_thrombocytopenia <- icd_data_before_sum_thrombocytopenia %>%
  inner_join(lab_data_before_sum_thrombocytopenia[, names(lab_data_before_sum_thrombocytopenia)[c(1:8, 10)]], by = join_by(studyid, time_frame))
adju_cohort <- readRDS(file = paste0(datapath, "adju_cohort.rds"))
before_sum_thrombocytopenia <- before_sum_thrombocytopenia[, names(before_sum_thrombocytopenia)[c(1:4, 6:13)]] %>%
  right_join(adju_cohort, by = join_by(studyid, time_frame))
glimpse(before_sum_thrombocytopenia)
table(before_sum_thrombocytopenia$thrombocytopenia_flag_icd, useNA = "ifany")
before_sum_thrombocytopenia$thrombocytopenia_flag_icd[is.na(before_sum_thrombocytopenia$thrombocytopenia_flag_icd)] <- 0
table(before_sum_thrombocytopenia$thrombocytopenia_G3_lab, useNA = "ifany")
before_sum_thrombocytopenia$thrombocytopenia_G3_lab[is.na(before_sum_thrombocytopenia$thrombocytopenia_G3_lab)] <- 0
table(before_sum_thrombocytopenia$thrombocytopenia_G2_lab, useNA = "ifany")
before_sum_thrombocytopenia$thrombocytopenia_G2_lab[is.na(before_sum_thrombocytopenia$thrombocytopenia_G2_lab)] <- 0
table(before_sum_thrombocytopenia$thrombocytopenia_G1_lab, useNA = "ifany")
before_sum_thrombocytopenia$thrombocytopenia_G1_lab[is.na(before_sum_thrombocytopenia$thrombocytopenia_G1_lab)] <- 0
table(before_sum_thrombocytopenia$missing_icd, useNA = "ifany")
before_sum_thrombocytopenia$missing_icd[is.na(before_sum_thrombocytopenia$missing_icd)] <- 1
table(before_sum_thrombocytopenia$missing_lab, useNA = "ifany")
before_sum_thrombocytopenia$missing_lab[is.na(before_sum_thrombocytopenia$missing_lab)] <- 1

# create thrombocytopenia_tox_indicator
before_sum_thrombocytopenia <- before_sum_thrombocytopenia %>%
  mutate(thrombocytopenia_tox_indicator = ifelse(thrombocytopenia_flag_icd == 1 & thrombocytopenia_G2_lab == 1, "Both",
                                                 ifelse(thrombocytopenia_flag_icd == 1 & thrombocytopenia_G2_lab == 0, "Only ICD",
                                                        ifelse(thrombocytopenia_flag_icd == 0 & thrombocytopenia_G2_lab == 1, "Only Grade 2", "Neither"))))

# subset before_sum_thrombocytopenia and export
adjuvant_thrombocytopenia <- before_sum_thrombocytopenia[, c("studyid", "time_frame", "thrombocytopenia_flag_icd", "thrombocytopenia_G3_lab", "thrombocytopenia_G2_lab", "thrombocytopenia_G1_lab",
                                                             "missing_icd", "missing_lab", "dx_dt", "siteid", "first_chemo_dt", "last_chemo_dt",
                                                             "surgery_d", "surgery_dt", "group", "thrombocytopenia_tox_indicator")]
write.csv(adjuvant_thrombocytopenia, file = paste0(datapath, "adjuvant_thrombocytopenia.csv"),
          row.names = F)


## Merge IOCD, Lab, and Chart
tox_sum_thrombocytopenia <- chart_before_sum_thrombocytopenia[, c("studyid", "thrombocytopenia_tox")] %>%
  right_join(before_sum_thrombocytopenia[before_sum_thrombocytopenia$time_frame == 4 & before_sum_thrombocytopenia$siteid==2,], by = join_by(studyid))
glimpse(tox_sum_thrombocytopenia) 
tox_sum_thrombocytopenia$thrombocytopenia_tox[is.na(tox_sum_thrombocytopenia$thrombocytopenia_tox)] <- 0


# export data
# saveRDS(icd_data_before_sum_thrombocytopenia, file = paste0(datapath, "icd_before_sum_thrombocytopenia.rds"))
# saveRDS(lab_data_before_thrombocytopenia, file = paste0(datapath, "lab_before_thrombocytopenia.rds"))
# saveRDS(lab_data_before_sum_thrombocytopenia, file = paste0(datapath, "lab_before_sum_thrombocytopenia.rds"))
# saveRDS(before_sum_thrombocytopenia, file = paste0(datapath, "before_sum_thrombocytopenia.rds"))
# saveRDS(tox_sum_thrombocytopenia, file = paste0(datapath, "tox_sum_thrombocytopenia.rds"))
# saveRDS(chart_before_sum_thrombocytopenia, file = paste0(datapath, "chart_before_sum_thrombocytopenia.rds"))
# before_sum_thrombocytopenia <- readRDS(file = paste0(datapath, "before_sum_thrombocytopenia.rds"))



#### Load Dataset ####
tox_sum_renal <- readRDS(file = paste0(datapath, "tox_sum.rds"))
tox_sum_renal <- tox_sum_renal[, c(1:9, 12:40)]
names(tox_sum_renal) <- c(names(tox_sum_renal)[1:3], "kidney_icd_n", "kidney_flag_icd", "missing_icd_kidney", 
                             paste0("kidney_", rep(c("G3", "G2", "G1"), 2), rep(c("_n", "_lab"), each = 3)), 
                             "missing_lab_kidney", names(tox_sum_renal)[14:38])
tox_sum_liver <- readRDS(file = paste0(datapath, "tox_sum_liver.rds"))
names(tox_sum_liver) <- c(names(tox_sum_liver)[1:5], "missing_icd_liver", 
                             paste0("liver_", rep(c("G2", "G3", "G1"), 2), rep(c("_n", "_lab"), each = 3)), 
                             "missing_lab_liver", names(tox_sum_liver)[14:38])
tox_sum_anemia <- readRDS(file = paste0(datapath, "tox_sum_anemia.rds"))
names(tox_sum_anemia)[6] <- "missing_icd_anemia"
names(tox_sum_anemia)[13] <- "missing_lab_anemia"
tox_sum_neutropenia <- readRDS(file = paste0(datapath, "tox_sum_neutropenia.rds"))
names(tox_sum_neutropenia)[6] <- "missing_icd_neutropenia"
names(tox_sum_neutropenia)[13] <- "missing_lab_neutropenia"
tox_sum_thrombocytopenia <- readRDS(file = paste0(datapath, "tox_sum_thrombocytopenia.rds"))
names(tox_sum_thrombocytopenia)[6] <- "missing_icd_thrombocytopenia"
names(tox_sum_thrombocytopenia)[13] <- "missing_lab_thrombocytopenia"
icd_data_tox_sum_neuro <- readRDS(file = paste0(datapath, "tox_sum_neuro.rds"))
names(icd_data_tox_sum_neuro)

#### merge ####
adjuvant_tox <- tox_sum_renal %>%
  left_join(tox_sum_liver[, names(tox_sum_liver)[c(1:13, 38)]], by = join_by(studyid, time_frame)) %>%
  left_join(tox_sum_anemia[, names(tox_sum_anemia)[c(1:13, 38)]], by = join_by(studyid, time_frame)) %>%
  left_join(tox_sum_neutropenia[, names(tox_sum_neutropenia)[c(1:13, 38, 39)]], by = join_by(studyid, time_frame)) %>%
  left_join(tox_sum_thrombocytopenia[, names(tox_sum_thrombocytopenia)[c(1:13, 38)]], by = join_by(studyid, time_frame)) %>%
  left_join(icd_data_tox_sum_neuro[, names(icd_data_tox_sum_neuro)[c(1:4, 7)]], by = join_by(studyid, time_frame))

#### check ####
vis_miss(adjuvant_tox, warn_large_data = FALSE)
adjuvant_tox$kidney_tox[is.na(adjuvant_tox$kidney_tox)] <- 0
# neuropathy
adjuvant_tox$missing_icd_neuropathy <- ifelse(is.na(adjuvant_tox$neuropathy_flag_icd), 1, 0)
adjuvant_tox$neuropathy_flag_icd[is.na(adjuvant_tox$neuropathy_flag_icd)] <- 0
unique(adjuvant_tox$time_frame)

#### merging BMI and RCT indicator into the data ####
new_data <- haven::read_sas(paste0(datapath, "obcd_cohort_03jan24.sas7bdat")) # OBCD Cohort
glimpse(new_data)
new_data <- new_data[, c("studyid", "dx_bmi", "chemo_bmi", "final_rct")]
adjuvant_tox <- adjuvant_tox %>%
  left_join(new_data, by = join_by(studyid))
glimpse(adjuvant_tox)

# export data
# saveRDS(adjuvant_tox, file = paste0(datapath, "adjuvant_tox.rds"))
# adjuvant_tox <- readRDS(file = paste0(datapath, "adjuvant_tox.rds"))

#### further exclusion ####
# dx_dt < 2004
adjuvant_tox[year(adjuvant_tox$dx_dt) < 2004, c("studyid", "time_frame", "dx_dt", "kidney_G3_lab", "liver_G3_lab", "anemia_G3_lab", "neutropenia_G3_lab", "thrombocytopenia_G3_lab")] %>% print(n = Inf, width = Inf) 
## 20 unique patients: 1 anemia, 1 thrombo
adjuvant_tox <- adjuvant_tox[year(adjuvant_tox$dx_dt) >= 2004,]
## 870 unique patients
adjuvant_tox[year(adjuvant_tox$dx_dt) < 2005, c("studyid", "time_frame", "dx_dt", "kidney_G3_lab", "liver_G3_lab", "anemia_G3_lab", "neutropenia_G3_lab", "thrombocytopenia_G3_lab")] %>% print(n = Inf, width = Inf) 
adjuvant_tox[year(adjuvant_tox$dx_dt) > 2015, c("studyid", "time_frame", "dx_dt", "kidney_G3_lab", "liver_G3_lab", "anemia_G3_lab", "neutropenia_G3_lab", "thrombocytopenia_G3_lab")] %>% print(n = Inf, width = Inf) 


# year(first_chemo_dt) > 2019
# none

# no RCT
adjuvant_tox[adjuvant_tox$final_rct == 1, c("studyid", "time_frame", "dx_dt", "kidney_G3_lab", "liver_G3_lab", "anemia_G3_lab", "neutropenia_G3_lab", "thrombocytopenia_G3_lab")] %>% print(n = Inf, width = Inf) 
## 16 patients: 1 anemia, 2 neutropenia
adjuvant_tox <- adjuvant_tox[adjuvant_tox$final_rct == 0,]
## 854 unique patients

#### New Variables ####
# check categories 
unique(adjuvant_tox$cycle1_dose_reduction90) # 0, 1, NA
unique(adjuvant_tox$stage_d) # 1, 2, 3
unique(adjuvant_tox$erpr_model) # 1, 2, NA
unique(adjuvant_tox$her2_model2) # 1, 2, NA
unique(adjuvant_tox$surgery_type_model) # 1, 2, NA
unique(adjuvant_tox$income_quart) # 1, 2, 3, 4, NA
unique(adjuvant_tox$race_eth_model) # 0, 1, 2, 3, 4, 5, NA
unique(adjuvant_tox$dosing_exclude) # 0, 1

# class of variables
adjuvant_tox$cycle1_dose_reduction90 <- factor(adjuvant_tox$cycle1_dose_reduction90, levels = 0:1, labels = c("N", "Y"))
adjuvant_tox$stage_d <- as.factor(adjuvant_tox$stage_d)
adjuvant_tox$erpr_model <- as.factor(adjuvant_tox$erpr_model)
adjuvant_tox$her2_model2 <- as.factor(adjuvant_tox$her2_model2)
adjuvant_tox$surgery_type_model <- as.factor(adjuvant_tox$surgery_type_model)
adjuvant_tox$income_quart <- as.factor(adjuvant_tox$income_quart)

# new variables
adjuvant_tox$race_new <- ifelse(!is.na(adjuvant_tox$race_eth_model) & adjuvant_tox$race_eth_model > 3, 4, adjuvant_tox$race_eth_model)
adjuvant_tox$race_new <- factor(adjuvant_tox$race_new, levels = 0:4, labels = c("White", "Black", "Asian", "Hispanic", "Other"))
adjuvant_tox$age_g <- ifelse(adjuvant_tox$dx_age < 50, "<50",
                                     ifelse(adjuvant_tox$dx_age < 65, "50-64",
                                            ifelse(adjuvant_tox$dx_age < 75, "65-74", "75+")))
adjuvant_tox$dx_year_grp <- ifelse(adjuvant_tox$dx_dt < as.Date("2012-05-01"), "Jan 2004 - April 2012", "May 2012 - Dec 2019")
adjuvant_tox$dx_bmi_grp <- ifelse(is.na(adjuvant_tox$dx_bmi), NA,
                                          ifelse(adjuvant_tox$dx_bmi < 18.5, 1,
                                                 ifelse(adjuvant_tox$dx_bmi < 25, 2,
                                                        ifelse(adjuvant_tox$dx_bmi < 30, 3,
                                                               ifelse(adjuvant_tox$dx_bmi < 35, 4,
                                                                      ifelse(adjuvant_tox$dx_bmi < 40, 5, 6)))))) 
adjuvant_tox$dx_bmi_grp <- factor(adjuvant_tox$dx_bmi_grp, levels = 1:6, labels = c("<18.5", "18.5 - <25", "25 - <30", "30 - <35", "35 - <40", "40+"))
adjuvant_tox$chemo_bmi_grp <- ifelse(is.na(adjuvant_tox$chemo_bmi), NA,
                                             ifelse(adjuvant_tox$chemo_bmi < 18.5, 1,
                                                    ifelse(adjuvant_tox$chemo_bmi < 25, 2,
                                                           ifelse(adjuvant_tox$chemo_bmi < 30, 3,
                                                                  ifelse(adjuvant_tox$chemo_bmi < 35, 4,
                                                                         ifelse(adjuvant_tox$chemo_bmi < 40, 5, 6)))))) 
adjuvant_tox$chemo_bmi_grp <- factor(adjuvant_tox$chemo_bmi_grp, levels = 1:6, labels = c("<18.5", "18.5 - <25", "25 - <30", "30 - <35", "35 - <40", "40+"))
adjuvant_tox$dx_year_grp_more <- ifelse(adjuvant_tox$dx_dt < as.Date("2008-01-01"), "Jan 2004 - Dec 2007",
                                                ifelse(adjuvant_tox$dx_dt < as.Date("2012-05-01"), "Jan 2008 - April 2012", "May 2012 - Dec 2015"))

anyNA(adjuvant_tox$neuropathy_flag_icd)
anyNA(adjuvant_tox$neuropathy_tox)
adjuvant_tox$neuropathy_tox[is.na(adjuvant_tox$neuropathy_tox)] <- 0

#### either ICD or lab ####
adjuvant_tox <- adjuvant_tox %>% 
  mutate(anemia_either_G1 = ifelse(anemia_flag_icd==1, 1,
                                   ifelse(anemia_G1_lab==1, 1, 0)),
         anemia_either_G2 = ifelse(anemia_flag_icd==1, 1,
                                   ifelse(anemia_G2_lab==1, 1, 0)),
         anemia_either_G3 = ifelse(anemia_flag_icd==1, 1,
                                   ifelse(anemia_G3_lab==1, 1, 0)),
         neutropenia_either_G1 = ifelse(neutropenia_flag_icd==1, 1,
                                        ifelse(neutropenia_G1_lab==1, 1, 0)),
         neutropenia_either_G2 = ifelse(neutropenia_flag_icd==1, 1,
                                        ifelse(neutropenia_G2_lab==1, 1, 0)),
         neutropenia_either_G3 = ifelse(neutropenia_flag_icd==1, 1,
                                        ifelse(neutropenia_G3_lab==1, 1, 0)),
         thrombocytopenia_either_G1 = ifelse(thrombocytopenia_flag_icd==1, 1,
                                             ifelse(thrombocytopenia_G1_lab==1, 1, 0)),
         thrombocytopenia_either_G2 = ifelse(thrombocytopenia_flag_icd==1, 1,
                                             ifelse(thrombocytopenia_G2_lab==1, 1, 0)),
         thrombocytopenia_either_G3 = ifelse(thrombocytopenia_flag_icd==1, 1,
                                             ifelse(thrombocytopenia_G3_lab==1, 1, 0))
         
  )
# export data
# saveRDS(adjuvant_tox, file = paste0(datapath, "adjuvant_tox_final.rds"))
# adjuvant_tox <- readRDS(file = paste0(datapath, "adjuvant_tox_final.rds"))

#### descriptive ####
a <- adjuvant_tox %>% ungroup() %>%  
  select(age_g, race_new, dx_bmi_grp, chemo_bmi_grp, income_quart, dx_year_grp, dx_year_grp_more, stage_d, 
         erpr_model, her2_model2, surgery_type_model) %>%
  tbl_summary(
    statistic = list(all_categorical() ~ "{n} ({p}%)"),
    missing_text = "(Missing)",
    digits = list(all_categorical() ~ c(0, 1))) 

# save as word doc
filename_1 = "table_1_tox.docx"

save_as_docx(`descriptive table` = as_flex_table(a), 
             path = paste0(outpath, filename_1))

a <- adjuvant_tox %>% ungroup() %>% 
  filter(missing_lab_anemia == 0 | missing_lab_neutropenia == 0 |
           missing_lab_thrombocytopenia == 0) %>%
  select(age_g, race_new, dx_bmi_grp, chemo_bmi_grp, income_quart, dx_year_grp, dx_year_grp_more, stage_d, 
         erpr_model, her2_model2, surgery_type_model) %>%
  tbl_summary(
    statistic = list(all_categorical() ~ "{n} ({p}%)"),
    missing_text = "(Missing)",
    digits = list(all_categorical() ~ c(0, 1))) 

# save as word doc
filename_1 = "table_1_1_tox.docx"

save_as_docx(`descriptive table` = as_flex_table(a), 
             path = paste0(outpath, filename_1))

a <- adjuvant_tox %>% ungroup() %>% 
  filter(missing_lab_anemia == 0 & missing_lab_neutropenia == 0 &
           missing_lab_thrombocytopenia == 0) %>%
  select(age_g, race_new, dx_bmi_grp, chemo_bmi_grp, income_quart, dx_year_grp, dx_year_grp_more, stage_d, 
         erpr_model, her2_model2, surgery_type_model) %>%
  tbl_summary(
    statistic = list(all_categorical() ~ "{n} ({p}%)"),
    missing_text = "(Missing)",
    digits = list(all_categorical() ~ c(0, 1))) 

# save as word doc
filename_1 = "table_1_2_tox.docx"

save_as_docx(`descriptive table` = as_flex_table(a), 
             path = paste0(outpath, filename_1))

a <- adjuvant_tox %>% ungroup() %>% 
  filter(missing_lab_anemia == 1 | missing_lab_kidney == 1 |
           missing_lab_liver == 1 | missing_lab_neutropenia == 1 |
           missing_lab_thrombocytopenia == 1) %>%
  select(age_g, race_new, dx_bmi_grp, chemo_bmi_grp, income_quart, dx_year_grp, dx_year_grp_more, stage_d, 
         erpr_model, her2_model2, surgery_type_model) %>%
  tbl_summary(
    statistic = list(all_categorical() ~ "{n} ({p}%)"),
    missing_text = "(Missing)",
    digits = list(all_categorical() ~ c(0, 1))) 

# save as word doc
filename_1 = "table_1_3_tox.docx"

save_as_docx(`descriptive table` = as_flex_table(a), 
             path = paste0(outpath, filename_1))


# who was missing neutrophils 
a <- adjuvant_tox %>% ungroup() %>% 
  filter(missing_lab_neutropenia == 1) %>%
  select(age_g, race_new, dx_bmi_grp, chemo_bmi_grp, income_quart, dx_year_grp, dx_year_grp_more, stage_d, 
         erpr_model, her2_model2, surgery_type_model) %>%
  tbl_summary(
    statistic = list(all_categorical() ~ "{n} ({p}%)"),
    missing_text = "(Missing)",
    digits = list(all_categorical() ~ c(0, 1))) 

# save as word doc
filename_1 = "table_1_neutrophil_tox.docx"

save_as_docx(`descriptive table` = as_flex_table(a), 
             path = paste0(outpath, filename_1))

a <- adjuvant_tox %>% ungroup() %>% 
  filter(missing_lab_neutropenia == 0) %>%
  select(age_g, race_new, dx_bmi_grp, chemo_bmi_grp, income_quart, dx_year_grp, dx_year_grp_more, stage_d, 
         erpr_model, her2_model2, surgery_type_model) %>%
  tbl_summary(
    statistic = list(all_categorical() ~ "{n} ({p}%)"),
    missing_text = "(Missing)",
    digits = list(all_categorical() ~ c(0, 1))) 

# save as word doc
filename_1 = "table_1_neutrophil_tox2.docx"

save_as_docx(`descriptive table` = as_flex_table(a), 
             path = paste0(outpath, filename_1))

#### Missing Lab ####
adjuvant_tox %>% 
  summarise(renal = 854 - sum(missing_lab_kidney),
            liver = 854 - sum(missing_lab_liver),
            anemia = 854 - sum(missing_lab_anemia),
            neutropenia = 854 - sum(missing_lab_neutropenia),
            thrombo = 854 - sum(missing_lab_thrombocytopenia))
adjuvant_tox %>% 
  summarise(renal = (854 - sum(missing_lab_kidney))/854*100,
            liver = (854 - sum(missing_lab_liver))/854*100,
            anemia = (854 - sum(missing_lab_anemia))/854*100,
            neutropenia = (854 - sum(missing_lab_neutropenia))/854*100,
            thrombo = (854 - sum(missing_lab_thrombocytopenia))/854*100)

#### validation ####
## anemia
a <- adjuvant_tox %>% # Ref = chart
  mutate(Vali_G3 = ifelse(anemia_tox == 0 & anemia_G3_lab == 0, "TN",
                         ifelse(anemia_tox == 1 & anemia_G3_lab == 1, "TP",
                                ifelse(anemia_tox == 0 & anemia_G3_lab == 1, "FP", "FN"))),
         Vali_G2 = ifelse(anemia_tox == 0 & anemia_G2_lab == 0, "TN",
                         ifelse(anemia_tox == 1 & anemia_G2_lab == 1, "TP",
                                ifelse(anemia_tox == 0 & anemia_G2_lab == 1, "FP", "FN"))),
         Vali_G1 = ifelse(anemia_tox == 0 & anemia_G1_lab == 0, "TN",
                         ifelse(anemia_tox == 1 & anemia_G1_lab == 1, "TP",
                                ifelse(anemia_tox == 0 & anemia_G1_lab == 1, "FP", "FN"))),
         Vali_icd = ifelse(anemia_flag_icd == 0 & anemia_tox == 0, "TN",
                           ifelse(anemia_flag_icd == 1 & anemia_tox == 1, "TP",
                                  ifelse(anemia_flag_icd == 0 & anemia_tox == 1, "FN", "FP"))))
b <- a %>% group_by(Vali_G3) %>% summarise(n = n())
c <- a %>% group_by(Vali_G2) %>% summarise(n = n())
d <- a %>% group_by(Vali_G1) %>% summarise(n = n())
e <- a %>% group_by(Vali_icd) %>% summarise(n = n())
## ICD
f <- as.table(matrix(c(136,80,44,594), nrow = 2, byrow = TRUE))
epi.tests(f, conf.level = 0.95)
## G1
f <- as.table(matrix(c(163,374,17,300), nrow = 2, byrow = TRUE))
epi.tests(f, conf.level = 0.95)
## G2
f <- as.table(matrix(c(112,104,68,570), nrow = 2, byrow = TRUE))
epi.tests(f, conf.level = 0.95)
## G3
f <- as.table(matrix(c(23,5,157,669), nrow = 2, byrow = TRUE))
epi.tests(f, conf.level = 0.95)
### agreement
a <- agree(adjuvant_tox[, c("anemia_flag_icd", "anemia_G1_lab")])
a
b <- 1.96*sqrt((a$value/100)*(1 - a$value/100)/854)
c((a$value/100) - b, (a$value/100) + b)*100
a <- agree(adjuvant_tox[, c("anemia_flag_icd", "anemia_G2_lab")])
a
b <- 1.96*sqrt((a$value/100)*(1 - a$value/100)/854)
c((a$value/100) - b, (a$value/100) + b)*100
a <- agree(adjuvant_tox[, c("anemia_flag_icd", "anemia_G3_lab")])
a
b <- 1.96*sqrt((a$value/100)*(1 - a$value/100)/854)
c((a$value/100) - b, (a$value/100) + b)*100
### kappa
a <- kappa2(adjuvant_tox[, c("anemia_flag_icd", "anemia_G1_lab")])
a
a$value + qnorm(c(0.025, 0.975))*(a$value/a$statistic)
a <- kappa2(adjuvant_tox[, c("anemia_flag_icd", "anemia_G2_lab")])
a
a$value + qnorm(c(0.025, 0.975))*(a$value/a$statistic)
a <- kappa2(adjuvant_tox[, c("anemia_flag_icd", "anemia_G3_lab")])
a
a$value + qnorm(c(0.025, 0.975))*(a$value/a$statistic)

## neutropenia
a <- adjuvant_tox %>% # Ref = chart
  mutate(Vali_G3 = ifelse(neutropenia_tox == 0 & neutropenia_G3_lab == 0, "TN",
                         ifelse(neutropenia_tox == 1 & neutropenia_G3_lab == 1, "TP",
                                ifelse(neutropenia_tox == 0 & neutropenia_G3_lab == 1, "FP", "FN"))),
         Vali_G2 = ifelse(neutropenia_tox == 0 & neutropenia_G2_lab == 0, "TN",
                         ifelse(neutropenia_tox == 1 & neutropenia_G2_lab == 1, "TP",
                                ifelse(neutropenia_tox == 0 & neutropenia_G2_lab == 1, "FP", "FN"))),
         Vali_G1 = ifelse(neutropenia_tox == 0 & neutropenia_G1_lab == 0, "TN",
                          ifelse(neutropenia_tox == 1 & neutropenia_G1_lab == 1, "TP",
                                 ifelse(neutropenia_tox == 0 & neutropenia_G1_lab == 1, "FP", "FN"))),
         Vali_icd = ifelse(neutropenia_flag_icd == 0 & neutropenia_tox == 0, "TN",
                           ifelse(neutropenia_flag_icd == 1 & neutropenia_tox == 1, "TP",
                                  ifelse(neutropenia_flag_icd == 0 & neutropenia_tox == 1, "FN", "FP"))))
b <- a %>% group_by(Vali_G3) %>% summarise(n = n())
c <- a %>% group_by(Vali_G2) %>% summarise(n = n())
d <- a %>% group_by(Vali_G1) %>% summarise(n = n())
e <- a %>% group_by(Vali_icd) %>% summarise(n = n())
## ICD
f <- as.table(matrix(c(160,70,77,547), nrow = 2, byrow = TRUE))
epi.tests(f, conf.level = 0.95)
## G1
f <- as.table(matrix(c(163,193,74,424), nrow = 2, byrow = TRUE))
epi.tests(f, conf.level = 0.95)
## G2
f <- as.table(matrix(c(150,129,87,488), nrow = 2, byrow = TRUE))
epi.tests(f, conf.level = 0.95)
## G3
f <- as.table(matrix(c(117,69,120,548), nrow = 2, byrow = TRUE))
epi.tests(f, conf.level = 0.95)
### agreement
a <- agree(adjuvant_tox[, c("neutropenia_flag_icd", "neutropenia_G1_lab")])
a
b <- 1.96*sqrt((a$value/100)*(1 - a$value/100)/854)
c((a$value/100) - b, (a$value/100) + b)*100
a <- agree(adjuvant_tox[, c("neutropenia_flag_icd", "neutropenia_G2_lab")])
a
b <- 1.96*sqrt((a$value/100)*(1 - a$value/100)/854)
c((a$value/100) - b, (a$value/100) + b)*100
a <- agree(adjuvant_tox[, c("neutropenia_flag_icd", "neutropenia_G3_lab")])
a
b <- 1.96*sqrt((a$value/100)*(1 - a$value/100)/854)
c((a$value/100) - b, (a$value/100) + b)*100
### kappa
a <- kappa2(adjuvant_tox[, c("neutropenia_flag_icd", "neutropenia_G1_lab")])
a
a$value + qnorm(c(0.025, 0.975))*(a$value/a$statistic)
a <- kappa2(adjuvant_tox[, c("neutropenia_flag_icd", "neutropenia_G2_lab")])
a
a$value + qnorm(c(0.025, 0.975))*(a$value/a$statistic)
a <- kappa2(adjuvant_tox[, c("neutropenia_flag_icd", "neutropenia_G3_lab")])
a
a$value + qnorm(c(0.025, 0.975))*(a$value/a$statistic)

## thrombocytopenia
a <- adjuvant_tox %>% # Ref = chart
  mutate(Vali_G3 = ifelse(thrombocytopenia_tox == 0 & thrombocytopenia_G3_lab == 0, "TN",
                         ifelse(thrombocytopenia_tox == 1 & thrombocytopenia_G3_lab == 1, "TP",
                                ifelse(thrombocytopenia_tox == 0 & thrombocytopenia_G3_lab == 1, "FP", "FN"))),
         Vali_G2 = ifelse(thrombocytopenia_tox == 0 & thrombocytopenia_G2_lab == 0, "TN",
                         ifelse(thrombocytopenia_tox == 1 & thrombocytopenia_G2_lab == 1, "TP",
                                ifelse(thrombocytopenia_tox == 0 & thrombocytopenia_G2_lab == 1, "FP", "FN"))),
         Vali_G1 = ifelse(thrombocytopenia_tox == 0 & thrombocytopenia_G1_lab == 0, "TN",
                          ifelse(thrombocytopenia_tox == 1 & thrombocytopenia_G1_lab == 1, "TP",
                                 ifelse(thrombocytopenia_tox == 0 & thrombocytopenia_G1_lab == 1, "FP", "FN"))),
         Vali_icd = ifelse(thrombocytopenia_flag_icd == 0 & thrombocytopenia_tox == 0, "TN",
                           ifelse(thrombocytopenia_flag_icd == 1 & thrombocytopenia_tox == 1, "TP",
                                  ifelse(thrombocytopenia_flag_icd == 0 & thrombocytopenia_tox == 1, "FN", "FP"))))
b <- a %>% group_by(Vali_G3) %>% summarise(n = n())
c <- a %>% group_by(Vali_G2) %>% summarise(n = n())
d <- a %>% group_by(Vali_G1) %>% summarise(n = n())
e <- a %>% group_by(Vali_icd) %>% summarise(n = n())
## ICD
f <- as.table(matrix(c(13,9,8,824), nrow = 2, byrow = TRUE))
epi.tests(f, conf.level = 0.95)
## G1
f <- as.table(matrix(c(20,119,1,714), nrow = 2, byrow = TRUE))
epi.tests(f, conf.level = 0.95)
## G2
f <- as.table(matrix(c(10,9,11,824), nrow = 2, byrow = TRUE))
epi.tests(f, conf.level = 0.95)
## G3
f <- as.table(matrix(c(4,2,17,831), nrow = 2, byrow = TRUE))
epi.tests(f, conf.level = 0.95)
### agreement
a <- agree(adjuvant_tox[, c("thrombocytopenia_flag_icd", "thrombocytopenia_G1_lab")])
a
b <- 1.96*sqrt((a$value/100)*(1 - a$value/100)/854)
c((a$value/100) - b, (a$value/100) + b)*100
a <- agree(adjuvant_tox[, c("thrombocytopenia_flag_icd", "thrombocytopenia_G2_lab")])
a
b <- 1.96*sqrt((a$value/100)*(1 - a$value/100)/854)
c((a$value/100) - b, (a$value/100) + b)*100
a <- agree(adjuvant_tox[, c("thrombocytopenia_flag_icd", "thrombocytopenia_G3_lab")])
a
b <- 1.96*sqrt((a$value/100)*(1 - a$value/100)/854)
c((a$value/100) - b, (a$value/100) + b)*100
### kappa
a <- kappa2(adjuvant_tox[, c("thrombocytopenia_flag_icd", "thrombocytopenia_G1_lab")])
a
a$value + qnorm(c(0.025, 0.975))*(a$value/a$statistic)
a <- kappa2(adjuvant_tox[, c("thrombocytopenia_flag_icd", "thrombocytopenia_G2_lab")])
a
a$value + qnorm(c(0.025, 0.975))*(a$value/a$statistic)
a <- kappa2(adjuvant_tox[, c("thrombocytopenia_flag_icd", "thrombocytopenia_G3_lab")])
a
a$value + qnorm(c(0.025, 0.975))*(a$value/a$statistic)


## anemia: either or
a <- adjuvant_tox %>% # Ref = chart
  mutate(Vali_G3 = ifelse(anemia_tox == 0 & anemia_either_G3  == 0, "TN",
                          ifelse(anemia_tox == 1 & anemia_either_G3  == 1, "TP",
                                 ifelse(anemia_tox == 0 & anemia_either_G3  == 1, "FP", "FN"))),
         Vali_G2 = ifelse(anemia_tox == 0 & anemia_either_G2  == 0, "TN",
                          ifelse(anemia_tox == 1 & anemia_either_G2 == 1, "TP",
                                 ifelse(anemia_tox == 0 & anemia_either_G2 == 1, "FP", "FN"))),
         Vali_G1 = ifelse(anemia_tox == 0 & anemia_either_G1 == 0, "TN",
                          ifelse(anemia_tox == 1 & anemia_either_G1 == 1, "TP",
                                 ifelse(anemia_tox == 0 & anemia_either_G1 == 1, "FP", "FN"))))
b <- a %>% group_by(Vali_G3) %>% summarise(n = n())
c <- a %>% group_by(Vali_G2) %>% summarise(n = n())
d <- a %>% group_by(Vali_G1) %>% summarise(n = n())
## G1
f <- as.table(matrix(c(173,383,7,291), nrow = 2, byrow = TRUE))
epi.tests(f, conf.level = 0.95)
## G2
f <- as.table(matrix(c(151,149,29,525), nrow = 2, byrow = TRUE))
epi.tests(f, conf.level = 0.95)
## G3
f <- as.table(matrix(c(138,82,42,592), nrow = 2, byrow = TRUE))
epi.tests(f, conf.level = 0.95)

## neutropenia: either or
a <- adjuvant_tox %>% # Ref = chart
  mutate(Vali_G3 = ifelse(neutropenia_tox == 0 & neutropenia_either_G3  == 0, "TN",
                          ifelse(neutropenia_tox == 1 & neutropenia_either_G3  == 1, "TP",
                                 ifelse(neutropenia_tox == 0 & neutropenia_either_G3  == 1, "FP", "FN"))),
         Vali_G2 = ifelse(neutropenia_tox == 0 & neutropenia_either_G2  == 0, "TN",
                          ifelse(neutropenia_tox == 1 & neutropenia_either_G2 == 1, "TP",
                                 ifelse(neutropenia_tox == 0 & neutropenia_either_G2 == 1, "FP", "FN"))),
         Vali_G1 = ifelse(neutropenia_tox == 0 & neutropenia_either_G1  == 0, "TN",
                          ifelse(neutropenia_tox == 1 & neutropenia_either_G1 == 1, "TP",
                                 ifelse(neutropenia_tox == 0 & neutropenia_either_G1 == 1, "FP", "FN"))))
b <- a %>% group_by(Vali_G3) %>% summarise(n = n())
c <- a %>% group_by(Vali_G2) %>% summarise(n = n())
d <- a %>% group_by(Vali_G1) %>% summarise(n = n())
## G1
f <- as.table(matrix(c(203,228,34,389), nrow = 2, byrow = TRUE))
epi.tests(f, conf.level = 0.95)
## G2
f <- as.table(matrix(c(199,171,38,446), nrow = 2, byrow = TRUE))
epi.tests(f, conf.level = 0.95)
## G3
f <- as.table(matrix(c(188,118,49,499), nrow = 2, byrow = TRUE))
epi.tests(f, conf.level = 0.95)

## thrombocytopenia: either or
a <- adjuvant_tox %>% # Ref = chart
  mutate(Vali_G3 = ifelse(thrombocytopenia_tox == 0 & thrombocytopenia_either_G3  == 0, "TN",
                          ifelse(thrombocytopenia_tox == 1 & thrombocytopenia_either_G3  == 1, "TP",
                                 ifelse(thrombocytopenia_tox == 0 & thrombocytopenia_either_G3  == 1, "FP", "FN"))),
         Vali_G2 = ifelse(thrombocytopenia_tox == 0 & thrombocytopenia_either_G2  == 0, "TN",
                          ifelse(thrombocytopenia_tox == 1 & thrombocytopenia_either_G2 == 1, "TP",
                                 ifelse(thrombocytopenia_tox == 0 & thrombocytopenia_either_G2 == 1, "FP", "FN"))),
         Vali_G1 = ifelse(thrombocytopenia_tox == 0 & thrombocytopenia_either_G1  == 0, "TN",
                          ifelse(thrombocytopenia_tox == 1 & thrombocytopenia_either_G1 == 1, "TP",
                                 ifelse(thrombocytopenia_tox == 0 & thrombocytopenia_either_G1 == 1, "FP", "FN"))))
b <- a %>% group_by(Vali_G3) %>% summarise(n = n())
c <- a %>% group_by(Vali_G2) %>% summarise(n = n())
d <- a %>% group_by(Vali_G1) %>% summarise(n = n())
## G1
f <- as.table(matrix(c(20,121,1,712), nrow = 2, byrow = TRUE))
epi.tests(f, conf.level = 0.95)
## G2
f <- as.table(matrix(c(16,16,5,817), nrow = 2, byrow = TRUE))
epi.tests(f, conf.level = 0.95)
## G3
f <- as.table(matrix(c(13,10,8,823), nrow = 2, byrow = TRUE))
epi.tests(f, conf.level = 0.95)

#### Sensitivity ####
## anemia
a <- adjuvant_tox[adjuvant_tox$missing_lab_anemia==0,] %>% # Ref = chart
  mutate(Vali_G3 = ifelse(anemia_tox == 0 & anemia_G3_lab == 0, "TN",
                          ifelse(anemia_tox == 1 & anemia_G3_lab == 1, "TP",
                                 ifelse(anemia_tox == 0 & anemia_G3_lab == 1, "FP", "FN"))),
         Vali_G2 = ifelse(anemia_tox == 0 & anemia_G2_lab == 0, "TN",
                          ifelse(anemia_tox == 1 & anemia_G2_lab == 1, "TP",
                                 ifelse(anemia_tox == 0 & anemia_G2_lab == 1, "FP", "FN"))),
         Vali_G1 = ifelse(anemia_tox == 0 & anemia_G1_lab == 0, "TN",
                          ifelse(anemia_tox == 1 & anemia_G1_lab == 1, "TP",
                                 ifelse(anemia_tox == 0 & anemia_G1_lab == 1, "FP", "FN"))),
         Vali_icd = ifelse(anemia_flag_icd == 0 & anemia_tox == 0, "TN",
                           ifelse(anemia_flag_icd == 1 & anemia_tox == 1, "TP",
                                  ifelse(anemia_flag_icd == 0 & anemia_tox == 1, "FN", "FP"))))
b <- a %>% group_by(Vali_G3) %>% summarise(n = n())
c <- a %>% group_by(Vali_G2) %>% summarise(n = n())
d <- a %>% group_by(Vali_G1) %>% summarise(n = n())
e <- a %>% group_by(Vali_icd) %>% summarise(n = n())
## ICD
f <- as.table(matrix(c(132,78,43,560), nrow = 2, byrow = TRUE))
epi.tests(f, conf.level = 0.95)
## G1
f <- as.table(matrix(c(163,374,12,264), nrow = 2, byrow = TRUE))
epi.tests(f, conf.level = 0.95)
## G2
f <- as.table(matrix(c(112,104,63,534), nrow = 2, byrow = TRUE))
epi.tests(f, conf.level = 0.95)
## G3
f <- as.table(matrix(c(23,5,152,633), nrow = 2, byrow = TRUE))
epi.tests(f, conf.level = 0.95)
### agreement
a <- agree(adjuvant_tox[adjuvant_tox$missing_lab_anemia==0, c("anemia_flag_icd", "anemia_G1_lab")])
a
b <- 1.96*sqrt((a$value/100)*(1 - a$value/100)/813)
c((a$value/100) - b, (a$value/100) + b)*100
a <- agree(adjuvant_tox[adjuvant_tox$missing_lab_anemia==0, c("anemia_flag_icd", "anemia_G2_lab")])
a
b <- 1.96*sqrt((a$value/100)*(1 - a$value/100)/813)
c((a$value/100) - b, (a$value/100) + b)*100
a <- agree(adjuvant_tox[adjuvant_tox$missing_lab_anemia==0, c("anemia_flag_icd", "anemia_G3_lab")])
a
b <- 1.96*sqrt((a$value/100)*(1 - a$value/100)/813)
c((a$value/100) - b, (a$value/100) + b)*100
### kappa
a <- kappa2(adjuvant_tox[adjuvant_tox$missing_lab_anemia==0, c("anemia_flag_icd", "anemia_G1_lab")])
a
a$value + qnorm(c(0.025, 0.975))*(a$value/a$statistic)
a <- kappa2(adjuvant_tox[adjuvant_tox$missing_lab_anemia==0, c("anemia_flag_icd", "anemia_G2_lab")])
a
a$value + qnorm(c(0.025, 0.975))*(a$value/a$statistic)
a <- kappa2(adjuvant_tox[adjuvant_tox$missing_lab_anemia==0, c("anemia_flag_icd", "anemia_G3_lab")])
a
a$value + qnorm(c(0.025, 0.975))*(a$value/a$statistic)

## neutropenia
a <- adjuvant_tox[adjuvant_tox$missing_lab_neutropenia==0,] %>% # Ref = chart
  mutate(Vali_G3 = ifelse(neutropenia_tox == 0 & neutropenia_G3_lab == 0, "TN",
                          ifelse(neutropenia_tox == 1 & neutropenia_G3_lab == 1, "TP",
                                 ifelse(neutropenia_tox == 0 & neutropenia_G3_lab == 1, "FP", "FN"))),
         Vali_G2 = ifelse(neutropenia_tox == 0 & neutropenia_G2_lab == 0, "TN",
                          ifelse(neutropenia_tox == 1 & neutropenia_G2_lab == 1, "TP",
                                 ifelse(neutropenia_tox == 0 & neutropenia_G2_lab == 1, "FP", "FN"))),
         Vali_G1 = ifelse(neutropenia_tox == 0 & neutropenia_G1_lab == 0, "TN",
                          ifelse(neutropenia_tox == 1 & neutropenia_G1_lab == 1, "TP",
                                 ifelse(neutropenia_tox == 0 & neutropenia_G1_lab == 1, "FP", "FN"))),
         Vali_icd = ifelse(neutropenia_flag_icd == 0 & neutropenia_tox == 0, "TN",
                           ifelse(neutropenia_flag_icd == 1 & neutropenia_tox == 1, "TP",
                                  ifelse(neutropenia_flag_icd == 0 & neutropenia_tox == 1, "FN", "FP"))))
b <- a %>% group_by(Vali_G3) %>% summarise(n = n())
c <- a %>% group_by(Vali_G2) %>% summarise(n = n())
d <- a %>% group_by(Vali_G1) %>% summarise(n = n())
e <- a %>% group_by(Vali_icd) %>% summarise(n = n())
## ICD
f <- as.table(matrix(c(134,57,53,369), nrow = 2, byrow = TRUE))
epi.tests(f, conf.level = 0.95)
## G1
f <- as.table(matrix(c(163,193,24,233), nrow = 2, byrow = TRUE))
epi.tests(f, conf.level = 0.95)
## G2
f <- as.table(matrix(c(150,129,37,297), nrow = 2, byrow = TRUE))
epi.tests(f, conf.level = 0.95)
## G3
f <- as.table(matrix(c(117,69,70,357), nrow = 2, byrow = TRUE))
epi.tests(f, conf.level = 0.95)
### agreement
a <- agree(adjuvant_tox[adjuvant_tox$missing_lab_neutropenia==0, c("neutropenia_flag_icd", "neutropenia_G1_lab")])
a
b <- 1.96*sqrt((a$value/100)*(1 - a$value/100)/613)
c((a$value/100) - b, (a$value/100) + b)*100
a <- agree(adjuvant_tox[adjuvant_tox$missing_lab_neutropenia==0, c("neutropenia_flag_icd", "neutropenia_G2_lab")])
a
b <- 1.96*sqrt((a$value/100)*(1 - a$value/100)/613)
c((a$value/100) - b, (a$value/100) + b)*100
a <- agree(adjuvant_tox[adjuvant_tox$missing_lab_neutropenia==0, c("neutropenia_flag_icd", "neutropenia_G3_lab")])
a
b <- 1.96*sqrt((a$value/100)*(1 - a$value/100)/613)
c((a$value/100) - b, (a$value/100) + b)*100
### kappa
a <- kappa2(adjuvant_tox[adjuvant_tox$missing_lab_neutropenia==0, c("neutropenia_flag_icd", "neutropenia_G1_lab")])
a
a$value + qnorm(c(0.025, 0.975))*(a$value/a$statistic)
a <- kappa2(adjuvant_tox[adjuvant_tox$missing_lab_neutropenia==0, c("neutropenia_flag_icd", "neutropenia_G2_lab")])
a
a$value + qnorm(c(0.025, 0.975))*(a$value/a$statistic)
a <- kappa2(adjuvant_tox[adjuvant_tox$missing_lab_neutropenia==0, c("neutropenia_flag_icd", "neutropenia_G3_lab")])
a
a$value + qnorm(c(0.025, 0.975))*(a$value/a$statistic)


## neutropenia: after 2007
a <- adjuvant_tox[adjuvant_tox$dx_year_grp_more!="Jan 2004 - Dec 2007",] %>% # Ref = chart
  mutate(Vali_G3 = ifelse(neutropenia_tox == 0 & neutropenia_G3_lab == 0, "TN",
                          ifelse(neutropenia_tox == 1 & neutropenia_G3_lab == 1, "TP",
                                 ifelse(neutropenia_tox == 0 & neutropenia_G3_lab == 1, "FP", "FN"))),
         Vali_G2 = ifelse(neutropenia_tox == 0 & neutropenia_G2_lab == 0, "TN",
                          ifelse(neutropenia_tox == 1 & neutropenia_G2_lab == 1, "TP",
                                 ifelse(neutropenia_tox == 0 & neutropenia_G2_lab == 1, "FP", "FN"))),
         Vali_G1 = ifelse(neutropenia_tox == 0 & neutropenia_G1_lab == 0, "TN",
                          ifelse(neutropenia_tox == 1 & neutropenia_G1_lab == 1, "TP",
                                 ifelse(neutropenia_tox == 0 & neutropenia_G1_lab == 1, "FP", "FN"))),
         Vali_icd = ifelse(neutropenia_flag_icd == 0 & neutropenia_tox == 0, "TN",
                           ifelse(neutropenia_flag_icd == 1 & neutropenia_tox == 1, "TP",
                                  ifelse(neutropenia_flag_icd == 0 & neutropenia_tox == 1, "FN", "FP"))))
b <- a %>% group_by(Vali_G3) %>% summarise(n = n())
c <- a %>% group_by(Vali_G2) %>% summarise(n = n())
d <- a %>% group_by(Vali_G1) %>% summarise(n = n())
e <- a %>% group_by(Vali_icd) %>% summarise(n = n())
## ICD
f <- as.table(matrix(c(124,62,44,335), nrow = 2, byrow = TRUE))
epi.tests(f, conf.level = 0.95)
## G1
f <- as.table(matrix(c(141,160,27,237), nrow = 2, byrow = TRUE))
epi.tests(f, conf.level = 0.95)
## G2
f <- as.table(matrix(c(129,108,39,289), nrow = 2, byrow = TRUE))
epi.tests(f, conf.level = 0.95)
## G3
f <- as.table(matrix(c(99,58,69,339), nrow = 2, byrow = TRUE))
epi.tests(f, conf.level = 0.95)
### agreement
a <- agree(adjuvant_tox[adjuvant_tox$dx_year_grp_more!="Jan 2004 - Dec 2007", c("neutropenia_flag_icd", "neutropenia_G1_lab")])
a
b <- 1.96*sqrt((a$value/100)*(1 - a$value/100)/565)
c((a$value/100) - b, (a$value/100) + b)*100
a <- agree(adjuvant_tox[adjuvant_tox$dx_year_grp_more!="Jan 2004 - Dec 2007", c("neutropenia_flag_icd", "neutropenia_G2_lab")])
a
b <- 1.96*sqrt((a$value/100)*(1 - a$value/100)/565)
c((a$value/100) - b, (a$value/100) + b)*100
a <- agree(adjuvant_tox[adjuvant_tox$dx_year_grp_more!="Jan 2004 - Dec 2007", c("neutropenia_flag_icd", "neutropenia_G3_lab")])
a
b <- 1.96*sqrt((a$value/100)*(1 - a$value/100)/565)
c((a$value/100) - b, (a$value/100) + b)*100
### kappa
a <- kappa2(adjuvant_tox[adjuvant_tox$dx_year_grp_more!="Jan 2004 - Dec 2007", c("neutropenia_flag_icd", "neutropenia_G1_lab")])
a
a$value + qnorm(c(0.025, 0.975))*(a$value/a$statistic)
a <- kappa2(adjuvant_tox[adjuvant_tox$dx_year_grp_more!="Jan 2004 - Dec 2007", c("neutropenia_flag_icd", "neutropenia_G2_lab")])
a
a$value + qnorm(c(0.025, 0.975))*(a$value/a$statistic)
a <- kappa2(adjuvant_tox[adjuvant_tox$dx_year_grp_more!="Jan 2004 - Dec 2007", c("neutropenia_flag_icd", "neutropenia_G3_lab")])
a
a$value + qnorm(c(0.025, 0.975))*(a$value/a$statistic)

## thrombocytopenia
a <- adjuvant_tox[adjuvant_tox$missing_lab_thrombocytopenia==0,] %>% # Ref = chart
  mutate(Vali_G3 = ifelse(thrombocytopenia_tox == 0 & thrombocytopenia_G3_lab == 0, "TN",
                          ifelse(thrombocytopenia_tox == 1 & thrombocytopenia_G3_lab == 1, "TP",
                                 ifelse(thrombocytopenia_tox == 0 & thrombocytopenia_G3_lab == 1, "FP", "FN"))),
         Vali_G2 = ifelse(thrombocytopenia_tox == 0 & thrombocytopenia_G2_lab == 0, "TN",
                          ifelse(thrombocytopenia_tox == 1 & thrombocytopenia_G2_lab == 1, "TP",
                                 ifelse(thrombocytopenia_tox == 0 & thrombocytopenia_G2_lab == 1, "FP", "FN"))),
         Vali_G1 = ifelse(thrombocytopenia_tox == 0 & thrombocytopenia_G1_lab == 0, "TN",
                          ifelse(thrombocytopenia_tox == 1 & thrombocytopenia_G1_lab == 1, "TP",
                                 ifelse(thrombocytopenia_tox == 0 & thrombocytopenia_G1_lab == 1, "FP", "FN"))),
         Vali_icd = ifelse(thrombocytopenia_flag_icd == 0 & thrombocytopenia_tox == 0, "TN",
                           ifelse(thrombocytopenia_flag_icd == 1 & thrombocytopenia_tox == 1, "TP",
                                  ifelse(thrombocytopenia_flag_icd == 0 & thrombocytopenia_tox == 1, "FN", "FP"))))
b <- a %>% group_by(Vali_G3) %>% summarise(n = n())
c <- a %>% group_by(Vali_G2) %>% summarise(n = n())
d <- a %>% group_by(Vali_G1) %>% summarise(n = n())
e <- a %>% group_by(Vali_icd) %>% summarise(n = n())
## ICD
f <- as.table(matrix(c(13,9,8,784), nrow = 2, byrow = TRUE))
epi.tests(f, conf.level = 0.95)
## G1
f <- as.table(matrix(c(20,119,1,674), nrow = 2, byrow = TRUE))
epi.tests(f, conf.level = 0.95)
## G2
f <- as.table(matrix(c(10,9,11,784), nrow = 2, byrow = TRUE))
epi.tests(f, conf.level = 0.95)
## G3
f <- as.table(matrix(c(4,2,17,791), nrow = 2, byrow = TRUE))
epi.tests(f, conf.level = 0.95)
### agreement
a <- agree(adjuvant_tox[adjuvant_tox$missing_lab_thrombocytopenia==0, c("thrombocytopenia_flag_icd", "thrombocytopenia_G1_lab")])
a
b <- 1.96*sqrt((a$value/100)*(1 - a$value/100)/814)
c((a$value/100) - b, (a$value/100) + b)*100
a <- agree(adjuvant_tox[adjuvant_tox$missing_lab_thrombocytopenia==0, c("thrombocytopenia_flag_icd", "thrombocytopenia_G2_lab")])
a
b <- 1.96*sqrt((a$value/100)*(1 - a$value/100)/814)
c((a$value/100) - b, (a$value/100) + b)*100
a <- agree(adjuvant_tox[adjuvant_tox$missing_lab_thrombocytopenia==0, c("thrombocytopenia_flag_icd", "thrombocytopenia_G3_lab")])
a
b <- 1.96*sqrt((a$value/100)*(1 - a$value/100)/814)
c((a$value/100) - b, (a$value/100) + b)*100
### kappa
a <- kappa2(adjuvant_tox[adjuvant_tox$missing_lab_thrombocytopenia==0, c("thrombocytopenia_flag_icd", "thrombocytopenia_G1_lab")])
a
a$value + qnorm(c(0.025, 0.975))*(a$value/a$statistic)
a <- kappa2(adjuvant_tox[adjuvant_tox$missing_lab_thrombocytopenia==0, c("thrombocytopenia_flag_icd", "thrombocytopenia_G2_lab")])
a
a$value + qnorm(c(0.025, 0.975))*(a$value/a$statistic)
a <- kappa2(adjuvant_tox[adjuvant_tox$missing_lab_thrombocytopenia==0, c("thrombocytopenia_flag_icd", "thrombocytopenia_G3_lab")])
a
a$value + qnorm(c(0.025, 0.975))*(a$value/a$statistic)

## anemia: either or
a <- adjuvant_tox[adjuvant_tox$missing_lab_anemia==0,] %>% # Ref = chart
  mutate(Vali_G3 = ifelse(anemia_tox == 0 & anemia_either_G3  == 0, "TN",
                          ifelse(anemia_tox == 1 & anemia_either_G3  == 1, "TP",
                                 ifelse(anemia_tox == 0 & anemia_either_G3  == 1, "FP", "FN"))),
         Vali_G2 = ifelse(anemia_tox == 0 & anemia_either_G2  == 0, "TN",
                          ifelse(anemia_tox == 1 & anemia_either_G2 == 1, "TP",
                                 ifelse(anemia_tox == 0 & anemia_either_G2 == 1, "FP", "FN"))),
         Vali_G1 = ifelse(anemia_tox == 0 & anemia_either_G1 == 0, "TN",
                          ifelse(anemia_tox == 1 & anemia_either_G1 == 1, "TP",
                                 ifelse(anemia_tox == 0 & anemia_either_G1 == 1, "FP", "FN"))))
b <- a %>% group_by(Vali_G3) %>% summarise(n = n())
c <- a %>% group_by(Vali_G2) %>% summarise(n = n())
d <- a %>% group_by(Vali_G1) %>% summarise(n = n())
## G1
f <- as.table(matrix(c(169,381,6,257), nrow = 2, byrow = TRUE))
epi.tests(f, conf.level = 0.95)
## G2
f <- as.table(matrix(c(147,147,28,491), nrow = 2, byrow = TRUE))
epi.tests(f, conf.level = 0.95)
## G3
f <- as.table(matrix(c(134,80,41,558), nrow = 2, byrow = TRUE))
epi.tests(f, conf.level = 0.95)

## neutropenia: either or
a <- adjuvant_tox[adjuvant_tox$missing_lab_neutropenia==0,] %>% # Ref = chart
  mutate(Vali_G3 = ifelse(neutropenia_tox == 0 & neutropenia_either_G3  == 0, "TN",
                          ifelse(neutropenia_tox == 1 & neutropenia_either_G3  == 1, "TP",
                                 ifelse(neutropenia_tox == 0 & neutropenia_either_G3  == 1, "FP", "FN"))),
         Vali_G2 = ifelse(neutropenia_tox == 0 & neutropenia_either_G2  == 0, "TN",
                          ifelse(neutropenia_tox == 1 & neutropenia_either_G2 == 1, "TP",
                                 ifelse(neutropenia_tox == 0 & neutropenia_either_G2 == 1, "FP", "FN"))),
         Vali_G1 = ifelse(neutropenia_tox == 0 & neutropenia_either_G1  == 0, "TN",
                          ifelse(neutropenia_tox == 1 & neutropenia_either_G1 == 1, "TP",
                                 ifelse(neutropenia_tox == 0 & neutropenia_either_G1 == 1, "FP", "FN"))))
b <- a %>% group_by(Vali_G3) %>% summarise(n = n())
c <- a %>% group_by(Vali_G2) %>% summarise(n = n())
d <- a %>% group_by(Vali_G1) %>% summarise(n = n())
## G1
f <- as.table(matrix(c(177,215,10,211), nrow = 2, byrow = TRUE))
epi.tests(f, conf.level = 0.95)
## G2
f <- as.table(matrix(c(173,158,14,268), nrow = 2, byrow = TRUE))
epi.tests(f, conf.level = 0.95)
## G3
f <- as.table(matrix(c(162,105,25,321), nrow = 2, byrow = TRUE))
epi.tests(f, conf.level = 0.95)

## neutropenia: either or after 2007
a <- adjuvant_tox[adjuvant_tox$dx_year_grp_more!="Jan 2004 - Dec 2007",] %>% # Ref = chart
  mutate(Vali_G3 = ifelse(neutropenia_tox == 0 & neutropenia_either_G3  == 0, "TN",
                          ifelse(neutropenia_tox == 1 & neutropenia_either_G3  == 1, "TP",
                                 ifelse(neutropenia_tox == 0 & neutropenia_either_G3  == 1, "FP", "FN"))),
         Vali_G2 = ifelse(neutropenia_tox == 0 & neutropenia_either_G2  == 0, "TN",
                          ifelse(neutropenia_tox == 1 & neutropenia_either_G2 == 1, "TP",
                                 ifelse(neutropenia_tox == 0 & neutropenia_either_G2 == 1, "FP", "FN"))),
         Vali_G1 = ifelse(neutropenia_tox == 0 & neutropenia_either_G1  == 0, "TN",
                          ifelse(neutropenia_tox == 1 & neutropenia_either_G1 == 1, "TP",
                                 ifelse(neutropenia_tox == 0 & neutropenia_either_G1 == 1, "FP", "FN"))))
b <- a %>% group_by(Vali_G3) %>% summarise(n = n())
c <- a %>% group_by(Vali_G2) %>% summarise(n = n())
d <- a %>% group_by(Vali_G1) %>% summarise(n = n())
## G1
f <- as.table(matrix(c(160,190,8,207), nrow = 2, byrow = TRUE))
epi.tests(f, conf.level = 0.95)
## G2
f <- as.table(matrix(c(157,144,11,253), nrow = 2, byrow = TRUE))
epi.tests(f, conf.level = 0.95)
## G3
f <- as.table(matrix(c(147,101,21,296), nrow = 2, byrow = TRUE))
epi.tests(f, conf.level = 0.95)

## thrombocytopenia: either or
a <- adjuvant_tox[adjuvant_tox$missing_lab_thrombocytopenia==0,] %>% # Ref = chart
  mutate(Vali_G3 = ifelse(thrombocytopenia_tox == 0 & thrombocytopenia_either_G3  == 0, "TN",
                          ifelse(thrombocytopenia_tox == 1 & thrombocytopenia_either_G3  == 1, "TP",
                                 ifelse(thrombocytopenia_tox == 0 & thrombocytopenia_either_G3  == 1, "FP", "FN"))),
         Vali_G2 = ifelse(thrombocytopenia_tox == 0 & thrombocytopenia_either_G2  == 0, "TN",
                          ifelse(thrombocytopenia_tox == 1 & thrombocytopenia_either_G2 == 1, "TP",
                                 ifelse(thrombocytopenia_tox == 0 & thrombocytopenia_either_G2 == 1, "FP", "FN"))),
         Vali_G1 = ifelse(thrombocytopenia_tox == 0 & thrombocytopenia_either_G1  == 0, "TN",
                          ifelse(thrombocytopenia_tox == 1 & thrombocytopenia_either_G1 == 1, "TP",
                                 ifelse(thrombocytopenia_tox == 0 & thrombocytopenia_either_G1 == 1, "FP", "FN"))))
b <- a %>% group_by(Vali_G3) %>% summarise(n = n())
c <- a %>% group_by(Vali_G2) %>% summarise(n = n())
d <- a %>% group_by(Vali_G1) %>% summarise(n = n())
## G1
f <- as.table(matrix(c(20,121,1,672), nrow = 2, byrow = TRUE))
epi.tests(f, conf.level = 0.95)
## G2
f <- as.table(matrix(c(16,16,5,777), nrow = 2, byrow = TRUE))
epi.tests(f, conf.level = 0.95)
## G3
f <- as.table(matrix(c(13,10,8,783), nrow = 2, byrow = TRUE))
epi.tests(f, conf.level = 0.95)

