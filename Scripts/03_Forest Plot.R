
################################################
##plot trial using ratio
##Forest plot of reported trial effect estimates by income category for each infection
vacc3 <- vacc1 %>% 
  mutate_at(vars(central_estimate_efficacy, upper_efficacy), function(x) if_else(x ==100, 99, x)) %>% 
  mutate(pe = 1- central_estimate_efficacy/100,
         u_ci = 1 - upper_efficacy/100,
         l_ci = 1- lower_efficacy/100,
         te = (l_ci-u_ci)/(2*1.96))

vacc3$p_value[vacc3$study_id %in% c(4699,8465,410)] =
  as.numeric(stringi::stri_replace_all_fixed(vacc3$p_value[vacc3$study_id %in% c(4699,8465,410)],"<",""))

list3=(vacc3$central_estimate_efficacy[vacc3$study_id %in% c(4699,8465,410)]/100) /((qnorm(as.numeric(vacc3$p_value[vacc3$study_id %in% c(4699,8465,410)]) ,lower.tail=FALSE)))

vacc3$te[vacc3$study_id %in% c(4699,8465,410)] = list3
# View(vacc3[vacc3$study_id %in% c(4699,8465,410),c("study_id","p_value","central_estimate_efficacy","pe","te")])
vacc3 %>% filter(is.na(te))

vacc3_plot <- vacc3 %>% 
  mutate(income_cat_fin = factor(income_cat_fin, levels = c("lm", "mx", "hm"), labels = c("Low", "Mixed", "High"))) %>% 
  arrange(who_atc_lbl, income_cat_fin) %>% 
  group_by(who_atc_lbl) %>% 
  mutate(study_rdr = letters[seq_along(who_atc_lbl)]) %>% 
  ungroup()

forestplot1 <- ggplot(vacc3_plot,
            aes(x = study_rdr,
                y = pe, ymin = pe-2*te, 
                ymax = pe+2*te, 
                colour = income_cat_fin)) + 
  geom_linerange() + 
  geom_point() + 
  facet_wrap(~who_atc_lbl, scales = "free_y") +
  scale_x_discrete("", breaks = NULL) +
  scale_y_continuous("rate ratio") +
  scale_colour_discrete("") +
  scale_linetype("") +
  geom_hline(yintercept = 1) +
  coord_flip() +
  theme_minimal()
forestplot1

###################################################################
#Reported Efficacy forest plot - ratio - infection titles
vacc3_plot1 <- vacc3 %>%
  mutate(income_cat_fin = factor(income_cat_fin, levels = c("lm", "mx", "hm"), labels = c("Low", "Mixed", "High"))) %>% 
  arrange(condition, income_cat_fin) %>% 
  group_by(condition) %>% 
  mutate(study_rdr = letters[seq_along(condition)]) %>% 
  ungroup()

forestplot2 <- ggplot(vacc3_plot1, 
            aes(x = study_rdr,
                y = pe, ymin = pe-2*te, 
                ymax = pe+2*te, 
                colour = income_cat_fin)) + 
  geom_linerange() + 
  geom_point() + 
  facet_wrap(~condition, scales = "free_y") +
  scale_x_discrete("", breaks = NULL) +
  scale_y_continuous("rate ratio") +
  scale_colour_discrete("") +
  scale_linetype("") +
  geom_hline(yintercept = 1) +
  coord_flip() +
  theme_minimal()
forestplot2

################################################################
##Forest plot - efficacy in high- and low-income settings and the difference between each, by infection and overall
data2 <- readxl::read_excel("Outputs/Covariate_Output_Model.xlsx", sheet = "Efficacy by infection")
data2_plot <- data2 %>% 
  arrange(Infection, title_cat) %>% 
  group_by(Infection) %>% 
  mutate(inf_rdr = letters[seq_along(Infection)]) %>% 
  ungroup()

forestplot3 <- ggplot(data2_plot %>% 
                        mutate(Infection=factor(Infection,
                                                levels = c("Hepatitis E","Cholera","Pneumococcal","Typhoid", "Influenza",
                                                           "Hepatitis B","Rotavirus","Overall" ))) %>%
                        mutate(title_cat = factor(title_cat, levels = c("Low Income", "Difference", "High Income"), labels = c("Low Income", "Difference", "High Income"))), 
                      aes(x = inf_rdr,
                          y = estimate, ymin = lower, 
                          ymax = upper, 
                          colour = title_cat)) + 
  geom_linerange() + 
  geom_point() +
  facet_wrap(~Infection, scales = "free_y") +
  scale_x_discrete("", breaks = NULL) +
  scale_color_discrete("category") +
  scale_y_continuous("odds ratio", limits = c(0,12), labels=seq(0,12,1),breaks = seq(0,12,1)) +
  #scale_y_continuous("odds ratio") +
  geom_hline(yintercept = 1) +
  labs(title="Forest Plot")+
  coord_flip() +
  theme_minimal()
forestplot3

########################################################################
##FIGURE 3: Forest Plots of the raw data and the modelled effects
#Raw Data / reported Efficacy Data
vacc_3=vacc3[,c("condition","who_atc_lbl","income_cat_fin","pe","te")]
vacc_3 <- vacc_3 %>% 
  mutate (lower = pe-2*te, upper = pe+2*te) %>%
  mutate (estimate_type="raw") %>%
  mutate(income_cat_fin = factor(income_cat_fin, levels = c("lm", "mx", "hm"), labels = c("Low Income", "Mixed Income", "High Income"))) %>%
  rename(c("Infection"="condition", "estimate"="pe", "title_cat"="income_cat_fin")) %>%
  select(Infection, estimate, lower, upper, title_cat, estimate_type)

#modelled effect Data
data2_plot <- data2 %>% 
  filter(!title_cat == "Difference") %>%
  arrange(Infection, title_cat) %>% 
  group_by(Infection) %>% 
  mutate (estimate_type="meta") %>%
  mutate(inf_rdr = letters[seq_along(Infection)]) %>% 
  ungroup() %>%
  select(Infection, estimate, lower, upper, title_cat, estimate_type)

#mergedData
combinedData <- rbind(vacc_3,data2_plot)

#forestplot 
vacc_data <- combinedData %>% 
  arrange(Infection, title_cat) %>% 
  group_by(Infection) %>% 
  mutate(inf_rdr = letters[seq_along(Infection)]) %>% 
  ungroup()

forestplot4 <- vacc_data %>% 
  mutate(Infection=factor(Infection,
                          levels = c("Hepatitis E","Cholera","Pneumococcal","Typhoid", "Influenza",
                                     "Hepatitis B","Rotavirus","Overall")), 
         estimate_type= factor(estimate_type, levels = c("raw", "meta"), ordered = T, labels = c("Trials", "Meta-analysis"))) %>% 
  ggplot(aes(x = inf_rdr,
             y = estimate, ymin = lower, 
             ymax = upper, 
             colour = title_cat)) + 
  geom_linerange(aes(linetype = estimate_type)) + 
  geom_point(aes(shape = estimate_type)) + 
  facet_wrap(~Infection, scales = "free_y") +
  scale_x_discrete("", breaks = NULL) +
  scale_y_continuous("rate ratio", limits = c(-1,2), labels=seq(-1,2,0.5),breaks = seq(-1,2,0.5)) +
  scale_colour_discrete("") +
  scale_linetype("") +
  scale_shape("")+
  geom_hline(yintercept = 1) +
  labs(title="Forest Plot")+
  coord_flip() +
  theme_minimal()
forestplot4   ##figure3

#######################################
# plot trials (using log odds ratio)
##Forest plot of reported trial effect estimates by income category for each infection
vacc1_plot <- vacc1 %>% 
  arrange(who_atc_lbl, income_cat_fin) %>% 
  group_by(who_atc_lbl) %>% 
  mutate(study_rdr = letters[seq_along(who_atc_lbl)]) %>% 
  ungroup()

forestplot5 <- ggplot(vacc1_plot %>% 
              mutate(income_cat_fin = factor(income_cat_fin, levels = c("lm", "mx", "hm"), labels = c("Low", "Mixed", "High"))), 
            aes(x = study_rdr,
                y = ce, ymin = ce-2*se, 
                ymax = ce+2*se, 
                colour = income_cat_fin)) + 
  geom_linerange() + 
  geom_point() + 
  facet_wrap(~who_atc_lbl, scales = "free_y") +
  scale_x_discrete("", breaks = NULL) +
  scale_color_discrete("Income category") +
  scale_y_continuous("rate ratio", breaks = log(c(0.1, 0.25, 1)), labels = c(0.1, 0.25, 1)) +
  geom_hline(yintercept = 0) +
  coord_flip() +
  theme_minimal()
forestplot5

########################################################
##Forest plot of model effect estimates
## Read in data ----
data1 <- readxl::read_excel("Outputs/Covariate_Output_Model.xlsx", sheet = "Effect Estimates")

forestplot6 <- ggplot(data1, 
            aes(x = Model,
                y = as.numeric(estimate), ymin = as.numeric(lower), 
                ymax = as.numeric(upper), 
                colour = Model)) + 
  geom_linerange() + 
  geom_point() +
  scale_x_discrete("", breaks = NULL ) +
  #scale_x_continuous(limits = c(0, 5)) +
  # xlim(0, 5)+
  scale_y_continuous("rate ratio", limits = c(0,5), labels=seq(0,5,0.5),breaks = seq(0,5,0.5)) +
  geom_hline(yintercept = 1) +
  labs(title="Forest Plot")+
  coord_flip() +
  theme_minimal()
forestplot6