##VACCINE EFFICACY ANALYSIS
library(tidyverse)
library(readxl)
library(WriteXLS)
library(writexl)
library(kableExtra)

if(R.version$os == "linux-gnu") setwd("~/Documents/Vaccine_Efficacy/Vaccine_Efficacy_LMICs")
## Read in data ----
vacc <- readxl::read_excel("Data/efficacy_unique.xlsx", sheet = 1)

vacc <- vacc %>% 
  janitor::clean_names()

## Clean duplicate labels ----
vacc <- vacc %>%
  rename(who_atc_lbl = who_atc_code) %>%
  mutate(who_atc_code = str_sub(who_atc_lbl, 1, 5),
         who_atc_lbl = if_else(who_atc_lbl == "J07BH Rota virus diarrhea vaccines", "J07BH Rotavirus diarrhea vaccines", who_atc_lbl))

sum(!duplicated(vacc$study_id)) #92 unique trials

## Drop trials where not at least one low income group in each; 63 trials ---- drop drug classes
vacc1 <- vacc %>%
  mutate(income_cat_fin = case_when(
    income_cat %in% c("Low Income", "Lower Middle Income") ~ "lm",
    income_cat == "Mixed Income" ~ "mx",
    TRUE ~ "hm"),
    income_cat_fin = factor(income_cat_fin, levels = c("hm", "mx", "lm"))) %>%
  group_by(who_atc_code) %>%
  mutate(gotli = if_else(any(income_cat_fin == "lm"), 1L, 0L)) %>%
  ungroup() %>%
  filter(gotli == 1) %>%
  select(-gotli)

## takes us down to 7 classes and 62 trials (1 trial has 4 estimates thus 66 observations)
vacc1 %>% 
  count(who_atc_lbl, income_cat_fin) %>% 
  spread(income_cat_fin, n, fill = "")

#drop unrelated trial (Safety, immunogenicity, and efficacy of live attenuated Vibrio cholerae 0139 vaccine prototype)
vacc1 <- vacc1 %>% 
  filter(!(study_id %in% c("4636")))

## takes us down to 7 classes and 61 trials (65 unique observations)
sum(!duplicated(vacc1$study_id)) #61 trials

## Recode to log-relative risk (i.e relative risk reduction)
## ce=log-relative risk (i.e relative risk reduction)
## Note sent 100 to 99, not actually possible to have 100% efficacy, would be infinite on original scale
vacc1 <- vacc1 %>% 
  mutate_at(vars(central_estimate_efficacy, upper_efficacy), function(x) if_else(x ==100, 99, x)) %>% 
  mutate(ce = log(1- central_estimate_efficacy/100),
         uci = log(1 - upper_efficacy/100),
         lci = log(1- lower_efficacy/100),
         se = (lci-uci)/(2*1.96))

#handling missing standard errors
vacc1 %>% filter(is.na(se)) %>% t
missingSE<-vacc1 %>% filter(is.na(se))

vacc1$p_value[vacc1$study_id %in% c(4699,8465,410)] =
  as.numeric(stringi::stri_replace_all_fixed(vacc1$p_value[vacc1$study_id %in% c(4699,8465,410)],"<",""))

list1=(vacc1$central_estimate_efficacy[vacc1$study_id %in% c(4699,8465,410)]/100) /((qnorm(as.numeric(vacc1$p_value[vacc1$study_id %in% c(4699,8465,410)]) ,lower.tail=FALSE)))

vacc1$se[vacc1$study_id %in% c(4699,8465,410)] = list1
missingse=vacc1 %>% filter(is.na(se))
missingse

#Harmonizing Covariates
vacc2=vacc1

#sample size
vacc1["participants_n_harmonised"]=vacc1$sample_size

#year
vacc1["year_harmonised"]=vacc1$study_period_from
table(vacc1$study_period_from)
table(vacc1$year_harmonised)

#Phase
vacc1["phase_harmonised"]=vacc1$phase
table(vacc1$phase)
vacc1$phase_harmonised[vacc1$phase %in% c("Phase II","Phase IIb")]="Phase II"
vacc1$phase_harmonised[vacc1$phase %in% c("Phase IIIb-IV")]="Phase IV"
vacc1$phase_harmonised[is.na(vacc1$phase)]="Phase III"
table(vacc1$phase_harmonised)

#Blinding
vacc1["blinding_harmonised"]=vacc1$blinding
table(vacc1$blinding)
vacc1$blinding_harmonised[vacc1$blinding %in% c("observer-blind","single (observer)-blind")]="double-blind"
vacc1$blinding_harmonised[vacc1$blinding %in% c("partially double-blind")]="single-blind"
vacc1$blinding_harmonised[vacc1$study_id %in% c("3675","332")]="single-blind"
vacc1$blinding_harmonised[is.na(vacc1$blinding_harmonised)]="open-label"
table(vacc1$blinding_harmonised)

#control
vacc1["controltype_harmonised"]=vacc1$control_type
table(vacc1$control_type)
vacc1$controltype_harmonised[vacc1$study_id %in% c("2901","2642","1405")]="placebo-controlled"
vacc1$controltype_harmonised[vacc1$study_id %in% c("10567","4308","3675","1863","75")]="usual care"
vacc1$controltype_harmonised[is.na(vacc1$controltype_harmonised)]="active-controlled trial"
table(vacc1$controltype_harmonised)

#method_of_administration
vacc1["administration_harmonised"]=vacc1$method_of_administration
table(vacc1$method_of_administration)
vacc1$administration_harmonised[vacc1$method_of_administration %in% c("Injection","Intramuscular Injection","Intramuscular or deep subcutaneous injection","Intranasally | Intramuscular")]="Intramuscular"
vacc1$administration_harmonised[vacc1$method_of_administration %in% c("Intranasally","Oral")]="Mucosal"
vacc1$administration_harmonised[is.na(vacc1$administration_harmonised)]="not specified"
table(vacc1$administration_harmonised)

#analysis
vacc1["analysis_harmonised"]=vacc1$analysis
table(vacc1$analysis)
vacc1$analysis_harmonised[is.na(vacc1$analysis)]="not specified"
table(vacc1$analysis_harmonised)

#vacc_type
table(vacc1$vaccine_type)

l1=c("adjuvanted","conjugate", "conjugate | saccharide", "polysaccharide", "polysaccharide | conjugate", 
     "recombinant","virus-like particle")
l2=c("attenuated","live human-bovine reassortant","Live-attenuated","lyophilized","lyophilized attenuated",
     "reassortant","reassortant | live-attenuated")
l3=c("Inactivated")

vacc1["vactype_harmonised"]="not specified"
vacc1$vactype_harmonised[vacc1$vaccine_type %in% l1 ]= "recombinant| polysaccharide | conjugate"
vacc1$vactype_harmonised[vacc1$vaccine_type %in% l2 ]= "Live-attenuated"
vacc1$vactype_harmonised[vacc1$vaccine_type %in% l3 ]="Inactivated Vaccines"
table(vacc1$vactype_harmonised)

#deal with spaces
vacc1$who_atc_lbl = gsub('[\r\n]', '', vacc1$who_atc_lbl) #new line spaces

vacc1$who_atc_lbl = stringi::stri_replace_all_fixed(vacc1$who_atc_lbl, "\u00a0" , " ") #replace it with space

#creating function
getVaccineTables = function(vacc1,vacname,env){
  t1=vacc1 %>% filter(who_atc_lbl==vacname)
  tt=aggregate(t1$participants_n_harmonised, list(t1$income_cat_fin), FUN=mean)
  colnames(tt)=c("grp","mean")
  m=matrix(0,nrow = 1,ncol = 3)
  colnames(m)=c("hm","mx","lm")
  tt$grp = as.character(tt$grp)
  for(i in 1:nrow(tt)){
    m[,tt$grp[i]]= tt$mean[i]
  }
  rownames(m) = "Sample Mean"
  
  k1=table(t1$phase_harmonised, t1$income_cat_fin)
  k2=table(t1$blinding_harmonised, t1$income_cat_fin)
  k3=table(t1$controltype_harmonised, t1$income_cat_fin)
  k4=table(t1$vactype_harmonised, t1$income_cat_fin)
  k5=table(t1$administration_harmonised, t1$income_cat_fin)
  k6=table(t1$analysis_harmonised, t1$income_cat_fin)
  
  k1=cbind(factors=c(""),levels=row.names(k1),k1)
  k2=cbind(factors=c(""),levels=row.names(k2),k2)
  k3=cbind(factors=c(""),levels=row.names(k3),k3)
  k4=cbind(factors=c(""),levels=row.names(k4),k4)
  k5=cbind(factors=c(""),levels=row.names(k5),k5)
  k6=cbind(factors=c(""),levels=row.names(k6),k6)
  m=cbind(factors=c(""),levels=row.names(m),m)
  
  k1[1,'factors']="Phase"
  k2[1,'factors']="Blinding"
  k3[1,'factors']="Control"
  k4[1,'factors']="Vacc_type"
  k5[1,'factors']="Method_of_administration"
  k6[1,'factors']="Type_of_Analysis"
  m[1,'factors']="Average"
  
  final = rbind(k1,k2,k3,k4,k5,k6,m)
  finald = as.data.frame(final,row.names = F)
  fname =  stringi::stri_replace_all_fixed(vacname," ",".")
  assign(fname,finald,envir = env)
  print(fname)
  
}

#creating druglist
druglist=c("J07AE Cholera vaccines","J07AL Pneumococcal vaccines","J07AP Typhoid vaccines",
           "J07BB Influenza vaccines","J07BC01 Hepatitis vaccines","J07BH Rotavirus diarrhea vaccines",
           "Hepatitis E")
listn = c() #declaring empty list

for(i in 1:length(druglist)){
  getVaccineTables(vacc1,druglist[i],environment())
  pname =  stringi::stri_replace_all_fixed(druglist[i]," ",".")
  listn = c(listn,pname)
}

#export to excel and setting sheet names
write_xlsx(setNames(lapply(listn,get),listn),
           path = "Outputs/Final.xlsx")

#export csv data file
write.csv(vacc1,"Data/Efficacy_harmonised.csv", row.names = F)

#Baseline characteristics table (proportions)
#percentages calculated in reference to numbers in income categories
#totals in each income category;hm=37, mx= 7, lm=21
t1=vacc1
tt=aggregate(t1$participants_n_harmonised, list(t1$income_cat_fin), FUN=mean)
colnames(tt)=c("grp","mean")
m=matrix(0,nrow = 1,ncol = 3)
colnames(m)=c("hm","mx","lm")
tt$grp = as.character(tt$grp)
for(i in 1:nrow(tt)){
  m[,tt$grp[i]]= tt$mean[i]
}
rownames(m) = "Sample Mean"
m=round(m,digits = 0)

k1=table(t1$phase_harmonised, t1$income_cat_fin)
k2=table(t1$blinding_harmonised, t1$income_cat_fin)
k3=table(t1$controltype_harmonised, t1$income_cat_fin)
k4=table(t1$vactype_harmonised, t1$income_cat_fin)
k5=table(t1$administration_harmonised, t1$income_cat_fin)
k6=table(t1$analysis_harmonised, t1$income_cat_fin)

p1=round(prop.table(table(t1$phase_harmonised, t1$income_cat_fin),2)*100,0)
p1=matrix(paste0(k1,"(",p1,"%)"),ncol = dim(k1)[2])
colnames(p1)=c("hm","mx","lm")
row.names(p1)=row.names(k1)

p2=round(prop.table(table(t1$blinding_harmonised, t1$income_cat_fin),2)*100,0)
p2=matrix(paste0(k2,"(",p2,"%)"),ncol = dim(k2)[2])
colnames(p2)=c("hm","mx","lm")
row.names(p2)=row.names(k2)

p3=round(prop.table(table(t1$controltype_harmonised, t1$income_cat_fin),2)*100,0)
p3=matrix(paste0(k3,"(",p3,"%)"),ncol = dim(k3)[2])
colnames(p3)=c("hm","mx","lm")
row.names(p3)=row.names(k3)

p4=round(prop.table(table(t1$vactype_harmonised, t1$income_cat_fin),2)*100,0)
p4=matrix(paste0(k4,"(",p4,"%)"),ncol = dim(k4)[2])
colnames(p4)=c("hm","mx","lm")
row.names(p4)=row.names(k4)

p5=round(prop.table(table(t1$administration_harmonised, t1$income_cat_fin),2)*100,0)
p5=matrix(paste0(k5,"(",p5,"%)"),ncol = dim(k5)[2])
colnames(p5)=c("hm","mx","lm")
row.names(p5)=row.names(k5)

p6=round(prop.table(table(t1$analysis_harmonised, t1$income_cat_fin),2)*100,0)
p6=matrix(paste0(k6,"(",p6,"%)"),ncol = dim(k6)[2])
colnames(p6)=c("hm","mx","lm")
row.names(p6)=row.names(k6)

p1=cbind(factors=c(""),levels=row.names(p1),p1)
p2=cbind(factors=c(""),levels=row.names(p2),p2)
p3=cbind(factors=c(""),levels=row.names(p3),p3)
p4=cbind(factors=c(""),levels=row.names(p4),p4)
p5=cbind(factors=c(""),levels=row.names(p5),p5)
p6=cbind(factors=c(""),levels=row.names(p6),p6)
m=cbind(factors=c(""),levels=row.names(m),m)

final = rbind(p1,p2,p3,p4,p5,p6,m)
finald = as.data.frame(final,row.names = F)

knitr::kable(finald) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed")) %>% 
  pack_rows(index=c("Phase"=3,"Blinding"=3,"Control"=3,"Vaccine Type"=4, "Method of Administration"=3, "Type of Analysis"=3, "Average"=1))

#introducing factor titles to save in csv
p1[1,'factors']="Phase"
p2[1,'factors']="Blinding"
p3[1,'factors']="Control"
p4[1,'factors']="Vacc_type"
p5[1,'factors']="Method_of_administration"
p6[1,'factors']="Type_of_Analysis"
m[1,'factors']="Average"  #participants' average

final = rbind(p1,p2,p3,p4,p5,p6,m)
finald = as.data.frame(final,row.names = F)

#export csv data file
write.csv(finald,"Outputs/base_characteristics.csv", row.names=F)


#MODELLING
#finding studies in mixed income (mx) group
mx_trials=vacc1[vacc1$income_cat_fin %in% c("mx"),]
mx.trials= as.data.frame(mx_trials,row.names = F)
write.csv(mx.trials,"Outputs/mx_trials.csv", row.names=F)

#drop the mx category
vacc1<-vacc1%>%
  filter(!(income_cat_fin %in% c("mx")))

## takes us down to 7 classes and 58 trials
sum(!duplicated(vacc1$study_id)) #58 unique
vacc1 %>% 
  count(who_atc_lbl, income_cat_fin) %>% 
  spread(income_cat_fin, n, fill = "")

#Drop high income trials with older age groups and no corresponding trials for low income countries)
# view(vacc1[vacc1$study_id %in% c("4542","4430","4428","4429","3675","3018","2165","974","977"),
#            c("study_id","condition","who_atc_lbl", "age_imputed_years")])
high_age_trials=vacc1[vacc1$study_id %in% c("4542","4430","4428","4429","3675","3018","2165","974","977"),]

# view(vacc1[!(vacc1$study_id %in% c("4542","4430","4428","4429","3675","3018","2165","974","977")),])

vacc1=vacc1[!(vacc1$study_id %in% c("4542","4430","4428","4429","3675","3018","2165","974","977")),]

## Fit Frequentist model as simple check, not valid estimate as lumping all together ----
library(metafor)
modf <- rma(ce ~ income_cat_fin , sei = se, data = vacc1)
summary(modf)
coef(modf)

## Run nested frequentist model (ie stratifying by WHO group) ----
vacc1 <- vacc1 %>% 
  mutate(year_harmonised = as.integer(year_harmonised),
         year_harmonised = if_else(is.na(year_harmonised), 
                                   as.integer(median(year_harmonised, na.rm = TRUE)),
                                   year_harmonised))
vacc1_nst <- vacc1 %>% 
  group_by(who_atc_lbl) %>% 
  nest() %>% 
  ungroup()
vacc1_nst$cnt <- map_int(vacc1_nst$data, nrow)
vacc1_nst <- vacc1_nst %>% 
  filter(cnt >2)
vacc1_nst$res <- map(vacc1_nst$data, ~ rma(ce ~ income_cat_fin , sei = se, data = .x))
vacc1_nst$coef <- map(vacc1_nst$res, broom::tidy)
vacc1_nst2 <- vacc1_nst %>% 
  select(who_atc_lbl, cnt, coef) %>% 
  unnest(coef)
write_csv(vacc1_nst2, "Outputs/freq_strat.csv")

#RDS
saveRDS(vacc1 %>% 
          select(ce, se, study_id, income_cat_fin, who_atc_lbl,
                 year_harmonised, phase_harmonised, blinding_harmonised, 
                 controltype_harmonised, administration_harmonised, analysis_harmonised, 
                 vactype_harmonised, age_imputed_years, sd_imputed_years), "Outputs/vacc_cleaned_data.Rds")
