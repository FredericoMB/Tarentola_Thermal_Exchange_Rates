##########################################################################################
###########        TARENTOLA MAURITANICA THERMAL EXCHANGE RATES ANALYSIS        ##########
###########                   Mochales-Riaño & Barroso et al.                   ##########
###########                    (built under R Version 4.3.2)                    ##########
##########################################################################################

rm(list = ls())

{ 
  library(plyr)             # data wrangling (must load before loading dplyr)
  library(dplyr)            # data wrangling
  library(tidyr)            # data wrangling
  library(ggplot2)          # plotting/ visualizing
    theme_set(theme_bw())   
  library(cowplot)          # plotting/ visualizing
  library(plotrix)          # plotting/ visualizing
  library(grid)             # plotting/ visualizing
  library(gridExtra)        # plotting/ visualizing
  library(lattice)          # plotting (for model diagnostics)
  library(lmodel2)          # fitting OLS models and getting CIs
  library(lme4)             # fitting mixed models
  library(lmerTest)         # extracting p-values from mixed models
  library(effectsize)       # calculating effect size for mixed models
  }

setwd("YOUR_WORKING_DIRECTORY")   # set your working directory


##########################################################
###########       1. FULL DATASET ANALYSIS       #########
##########################################################


#### 1.1. Data Importing and Wrangling ####

  cumul_t_change <- read.csv("cumulative_temperatures_data.csv", header = T)

  cumul_t_change <- cumul_t_change %>% filter(time!=0)  # Eliminating Time = 0 which would force intercept to be = 0

  { cumul_t_change$jump <- as.factor(cumul_t_change$jump)
    cumul_t_change$treat <- as.factor(cumul_t_change$treat)
    cumul_t_change$heat <- as.factor(cumul_t_change$heat)
    cumul_t_change$id <- as.factor(cumul_t_change$id)
    cumul_t_change$population <- as.factor(cumul_t_change$population)
    cumul_t_change$heat2 <- as.factor(cumul_t_change$heat2) 
    }

  str(cumul_t_change)
  
  cumul_t_change_tidy <- cumul_t_change %>% 
    pivot_longer(c("Snout","Eye","Head","Dorsum","Leg","Foot","Tail"), 
                 names_to = "position", values_to = "temp")
  
  cumul_t_change_tidy$position <- as.factor(cumul_t_change_tidy$position)
  
  str(cumul_t_change_tidy)

  
  # Separating the Cooling and Heating Datasets:
  
    cumul_t_change_tidy_H <- cumul_t_change_tidy[cumul_t_change_tidy$heat=="Heating",]
    cumul_t_change_tidy_H$heat <- droplevels(cumul_t_change_tidy_H$heat)

    cumul_t_change_tidy_C <- cumul_t_change_tidy[cumul_t_change_tidy$heat=="Cooling",]
    cumul_t_change_tidy_C$heat <- droplevels(cumul_t_change_tidy_C$heat)


#### 1.2. Sample Description ####

   # Mass
   
    data_long <- cumul_t_change[cumul_t_change$time==20 ,] %>% 
       pivot_longer(c("mass_initial","mass_final"), names_to = "i_vs_f", values_to = "mass")  
    
    summary_mass<- ddply(data_long, c("treat"), summarise,
                         N    = as.numeric(length(unique(data_long$id))),
                         Mean = mean(mass, na.rm=TRUE),
                         SD   = sd(mass, na.rm=TRUE),
                         SE   = SD / sqrt(N))
    summary_mass 
    
    t.test(data_long[data_long$treat=="Heliothermy",]$mass~data_long[data_long$treat=="Heliothermy",]$i_vs_f)
    t.test(data_long[data_long$treat=="Thigmothermy",]$mass~data_long[data_long$treat=="Thigmothermy",]$i_vs_f)
    
    mass_plot <- ggplot(data_long, aes(x=treat, y=mass)) +
       geom_boxplot(aes(fill=reorder(i_vs_f, -mass)), alpha=0.3) +
       stat_summary(aes(group=reorder(i_vs_f, -mass)), fun=mean, geom="point", shape=4, size=4, 
                    position = position_dodge2(width = 0.75, preserve = "single")) +
       xlab("") + ylab("\n Mass (grams) \n")  + 
       scale_x_discrete(labels=c("\nHeliothermy","\nThigmothermy")) +
       scale_fill_discrete(labels=c("Initial Mass   ", "Final Mass")) +
       geom_hline(aes(yintercept=-Inf), linewidth=0.9) + geom_vline(aes(xintercept=-Inf), linewidth=0.9) +
       theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
             axis.line = element_blank(), panel.border = element_blank()) +
       theme(axis.text.x = element_text(size = 11),
             axis.text.y = element_text(size = 11), axis.title = element_text(size=14),
             strip.text.x = element_text(size = 16), strip.text.y = element_text(size = 16)) +
       theme(legend.position = "bottom", legend.direction="horizontal", legend.title=element_text(size=12), 
             legend.text=element_text(size=11)) + guides(fill=guide_legend(title=" "))
    
     
   # Snout-to-Vent Length
   
     data_SVL <- cumul_t_change[cumul_t_change$time==20 & cumul_t_change$treat=="Heliothermy" & cumul_t_change$heat2==2,]    
     
     mean(data_SVL$svl)
     (sd(data_SVL$svl, na.rm=T))/(sqrt(as.numeric(length(unique(data_SVL$id)))))
     
     SVL_plot <- ggplot(data_SVL, aes(y=svl)) +
       geom_boxplot() + coord_flip() + geom_point(x=0, y=mean(data_SVL$svl), size=4, shape=4) + 
       xlab("\n") + ylab("\n SVL (mm) \n") +  ggtitle("\n") + 
       geom_vline(aes(xintercept=-Inf), linewidth=0.9) +
       theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
             axis.line = element_blank(), panel.border = element_blank()) +
       theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
             axis.text.x = element_text(size = 11), axis.title = element_text(size=14),
             strip.text.x = element_text(size = 16), strip.text.y = element_text(size = 16)) 
          
   # Visualizing the Mass and SVL Plots
   
     body_size_plots <- grid.arrange(mass_plot, SVL_plot, nrow = 2, layout_matrix = cbind(c(1,1,2), c(1,1,2)))
     
        
#### 1.3. Plotting the Data ####
        
  # Plotting the Time:Body_Part Interaction for Heliothermic and Thigmothermic Heating and Cooling 
    
    names1 <- c("Heliothermy"="HELIOTHERMY", "Thigmothermy"="THIGMOTHERMY", "2"="COOLING", "1"="HEATING")
    
    Plot_1 <- ggplot(cumul_t_change_tidy, aes(x=time, y=temp, group=factor(position))) +
      geom_jitter(col="grey80", size=0.2) +
      stat_smooth(aes(colour=position), method="lm", formula = y ~ poly(x, 2), se = T, alpha=0.25, linewidth=1.1) +
      facet_grid(rows = vars(heat2), cols = vars(treat), scales="free", margins=F, labeller = as_labeller(names1)) +
      scale_color_manual(values=c(Dorsum="cornflowerblue", Eye="forestgreen", Head="firebrick2", 
                                  Leg="burlywood4", Foot="orchid3", Snout="paleturquoise3", Tail="gold1")) +
      xlab("\n Time (s) \n") + ylab("\n Cumulative  ∆ Temperature (ºC) \n")  + 
      coord_cartesian(xlim = c(0,600)) + 
      geom_hline(aes(yintercept=-Inf), linewidth=0.9) + geom_vline(aes(xintercept=-Inf), linewidth=0.9) + 
      labs(colour="    Body Part:    ") + guides(colour=guide_legend(nrow=1), byrow = TRUE) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
            axis.text.y = element_text(size = 11), axis.title = element_text(size=14),
            strip.text.x = element_text(size = 16), strip.text.y = element_text(size = 16)) +
      theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
            axis.line = element_blank(), panel.border = element_blank(), panel.spacing = unit(1, "cm"),
            strip.background = element_rect(color = "white", linewidth = 1)) +
      theme(legend.position = "bottom", legend.direction="horizontal", legend.title=element_text(size=14), 
            legend.text=element_text(size=12, margin = margin(r = 0.5, unit = 'cm'))) 
    
    Plot_1
    
    
  # Plotting the Time:Treatment (Heliothermy vs Thigmothermy) Interaction for Heating and Cooling per Body Part 
    
    names2 <- c("Dorsum"="Dorsum", "Eye"="Eye", "Head"="Head", "Leg"="Leg", "Foot"="Foot", "Snout"="Snout", "Tail"="Tail",
                "2"="COOLING", "1"="HEATING")
    
    Plot_2 <- ggplot(cumul_t_change_tidy, aes(x=time,y=abs(temp), group=factor(treat))) +
      geom_jitter(aes(shape=treat), col="grey80", size=0.9) +
      stat_smooth(aes(colour=position, linetype=treat), method="lm", formula = y ~ poly(x, 2), se = T, alpha=0.25, linewidth=1.3) +
      facet_grid(rows = vars(position), cols = vars(heat2), scales="fixed", margins=F, labeller = as_labeller(names2)) +
      scale_color_manual(values=c(Dorsum="cornflowerblue", Eye="forestgreen", Head="firebrick2", 
                                  Leg="burlywood4", Foot="orchid3", Snout="paleturquoise3", Tail="gold1")) +
      scale_shape_manual(values=c(Heliothermy=19, Thigmothermy=1)) +
      xlab("\n Time (s)") + ylab("\n Absolute Cumulative ∆ Temperature (ºC) \n")  + 
      coord_cartesian(xlim = c(0,600)) + 
      geom_hline(aes(yintercept=-Inf), linewidth=0.9) + geom_vline(aes(xintercept=-Inf), linewidth=0.9) +
      labs(colour="Body Part: \n", linetype="\n\n Treatment: \n", shape="") + 
      guides(colour = guide_legend(order=1), byrow = T, 
             linetype = guide_legend(order=2),
             shape = guide_legend(order=3)) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11), 
            axis.text.y = element_text(size = 11), axis.title = element_text(size=15),
            strip.text.x = element_text(size = 15), strip.text.y = element_text(size = 14)) +
      theme(legend.position = "right", legend.box="vertical", legend.direction="vertical") +
      theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
            axis.line = element_blank(), panel.border = element_blank(), panel.spacing = unit(1, "line"),
            strip.background = element_rect(color = "white", linewidth = 1))
    
    Plot_2

        
#### 1.4. Fitting Mixed Effects Models and Model Selection ####
    
# For the Heating Dataset:
  
  # Full model 
   
    h_mod1 <- lmer(temp ~ poly(time,2) * treat * position + urine + posture + mass_initial*poly(time,2) + population + 
                     (1|id), data =  cumul_t_change_tidy_H, REML=T, na.action = na.omit)
    anova(h_mod1)
      
  # Removing ID Random Effect
    
    h_mod2 <- lm(temp ~ poly(time,2) * treat * position + urine + posture + mass_initial*poly(time,2) + population, 
                 data =  cumul_t_change_tidy_H,  na.action = na.omit)
    
    anova(h_mod1, h_mod2, test=T)  # SIGN. ChiSq = 843.43, p<2.2e-16
    
  # Removing MASS_INITIAL:POLY(TIME,2) Interaction 
    
    h_mod1A <- lmer(temp ~ poly(time,2) * treat * position + urine + posture + mass_initial + population + 
                      (1|id), data =  cumul_t_change_tidy_H, REML=T, na.action = na.omit)
    
    anova(h_mod1, h_mod1A, test=T)  # SIGN. ChiSq = 8.4281, p=0.01479
    
  # Removing POSITION from the POLY(TIME,2):TREAT:POSITION Interaction
    
    h_mod1B <- lmer(temp ~ poly(time,2) * treat + position + urine + posture + mass_initial*poly(time,2) + population + 
                      (1|id), data =  cumul_t_change_tidy_H, REML=T, na.action = na.omit)
    
    anova(h_mod1, h_mod1B, test=T)  # SIGN. ChiSq = 4735.5, p<2.2e-16
    
  # Removing TREAT from the POLY(TIME,2):TREAT:POSITION Interaction
    
    h_mod1C <- lmer(temp ~ poly(time,2) * position + treat + urine + posture + mass_initial*poly(time,2) + population + 
                      (1|id), data =  cumul_t_change_tidy_H, REML=T, na.action = na.omit)
    
    anova(h_mod1, h_mod1C, test=T)  # SIGN. ChiSq = 4594.8, p<2.2e-16
    
  # Removing POLY(TIME,2) from the POLY(TIME,2):TREAT:POSITION Interaction
    
    h_mod1D <- lmer(temp ~ poly(time,2) + treat * position + urine + posture + mass_initial*poly(time,2) + population + 
                      (1|id), data =  cumul_t_change_tidy_H, REML=T, na.action = na.omit)
    
    anova(h_mod1, h_mod1D, test=T)  # SIGN. ChiSq = 1317.6, p<2.2e-16
    
  # Removing POPULATION
    
    h_mod1E <- lmer(temp ~ poly(time,2) * treat * position + urine + posture + mass_initial*poly(time,2) + (1|id), 
                    data =  cumul_t_change_tidy_H, REML=T, na.action = na.omit)
    
    anova(h_mod1, h_mod1E, test=T)  # NOT SIGN. ChiSq = 2.4973, p=0.4758
       
  # Removing POSTURE
    
    h_mod1Ea <- lmer(temp ~ poly(time,2) * treat * position + urine + mass_initial*poly(time,2) + (1|id), 
                     data =  cumul_t_change_tidy_H, REML=T, na.action = na.omit)
    
    anova(h_mod1E, h_mod1Ea, test=T)  # NOT SIGN. ChiSq = 0.6666, p=0.4143
      
  # Removing URINE
    
    h_mod1Eaa <- lmer(temp ~ poly(time,2) * treat * position + mass_initial*poly(time,2) + (1|id), 
                      data =  cumul_t_change_tidy_H, REML=T, na.action = na.omit)
    
    anova(h_mod1Ea, h_mod1Eaa, test=T)  # SIGN. ChiSq = 14.159, p=0.000168
     
  # Removing TIME^2
    
    h_mod1Eab <- lmer(temp ~ poly(time,1) * treat * position + mass_initial*poly(time,1) + urine + (1|id), 
                      data =  cumul_t_change_tidy_H, REML=T, na.action = na.omit)
    
    anova(h_mod1Ea, h_mod1Eab, test=T)  # SIGN. ChiSq = 1944.4, p<2.2e-16
     
  # Final Model for Heating Dataset
    
    h_mod_final <- h_mod1Ea
    anova(h_mod_final)
        
  # Diagnostic plots to see if data follows model's assumptions
    
    plot(h_mod_final, resid(.,scaled=T)~fitted(.), abline=0, cex=0.5)
    hist(residuals(h_mod_final), main = "Histogram of Residuals")
    qqmath(resid(h_mod_final, scaled=T), abline=c(0,1), strip=T)
    qqmath(ranef(h_mod_final, condVar = T), strip = T, abline=c(0,1))$id
    

# For the Cooling Dataset: 

  # Full model
    
    c_mod1 <- lmer(temp ~ poly(time,2) * treat * position + urine + posture + mass_initial*poly(time,2) + population + (1|id), 
                   data =  cumul_t_change_tidy_C, REML=T, na.action = na.omit)
    anova(c_mod1)
    
  # Removing ID Random Effect
    
    c_mod2 <- lm(temp ~ poly(time,2) * treat * position + urine + posture + mass_initial*poly(time,2) + population, 
                 data =  cumul_t_change_tidy_C,  na.action = na.omit)
    
    anova(c_mod1, c_mod2, test=T)  # SIGN. ChiSq = 1116.8, p<2.2e-16
    
  # Removing MASS_INITIAL:POLY(TIME,2) Interaction
    
    c_mod1A <- lmer(temp ~ poly(time,2) * treat * position + urine + posture + mass_initial + population + (1|id), 
                    data =  cumul_t_change_tidy_C, REML=T, na.action = na.omit)
    
    anova(c_mod1, c_mod1A, test=T)  # SIGN. ChiSq = 25.798, p=2.5e-06
    
  # Removing POSITION from the POLY(TIME,2):TREAT:POSITION Interaction
    
    c_mod1B <- lmer(temp ~ poly(time,2) * treat + position + urine + posture + mass_initial*poly(time,2) + population + (1|id), 
                    data =  cumul_t_change_tidy_C, REML=T, na.action = na.omit)
    
    anova(c_mod1, c_mod1B, test=T)  # SIGN. ChiSq = 3673.6, p<2.2e-16
      
  # Removing TREAT from the POLY(TIME,2):TREAT:POSITION Interaction
    
    c_mod1C <- lmer(temp ~ poly(time,2) * position + treat + urine + posture + mass_initial*poly(time,2) + population + (1|id), 
                    data =  cumul_t_change_tidy_C, REML=T, na.action = na.omit)
    
    anova(c_mod1, c_mod1C, test=T)  # SIGN. ChiSq = 2904.8, p<2.2e-16
    
  # Removing POLY(TIME,2) from the POLY(TIME,2):TREAT:POSITION Interaction
    
    c_mod1D <- lmer(temp ~ poly(time,2) + treat * position + urine + posture + mass_initial*poly(time,2) + population + (1|id), 
                    data =  cumul_t_change_tidy_C, REML=T, na.action = na.omit)
    
    anova(c_mod1, c_mod1D, test=T)  # SIGN. ChiSq = 1020.3, p<2.2e-16
        
  # Removing POPULATION
    
    c_mod1E <- lmer(temp ~ poly(time,2) * treat * position + urine + posture + mass_initial*poly(time,2) + (1|id), 
                    data =  cumul_t_change_tidy_C, REML=T, na.action = na.omit)
    
    anova(c_mod1, c_mod1E, test=T)  # NOT SIGN. ChiSq = 6.2947, p=0.09812
       
  # Removing POSTURE
    
    c_mod1Ea <- lmer(temp ~ poly(time,2) * treat * position + urine + mass_initial*poly(time,2) + (1|id), 
                     data =  cumul_t_change_tidy_C, REML=T, na.action = na.omit)
    
    anova(c_mod1E, c_mod1Ea, test=T)  # SIGN. ChiSq = 28.283, p=1.048e-07
        
  # Removing URINE
    
    c_mod1Eb <- lmer(temp ~ poly(time,2) * treat * position + posture + mass_initial*poly(time,2) + (1|id), 
                     data =  cumul_t_change_tidy_C, REML=T, na.action = na.omit)
    
    anova(c_mod1E, c_mod1Eb, test=T)  # SIGN. ChiSq = 75.177, p<2.2e-16
       
  # Removing TIME^2
    
    c_mod1Ec <- lmer(temp ~ poly(time,1) * treat * position + urine + posture + mass_initial*poly(time,1) + (1|id), 
                     data =  cumul_t_change_tidy_C, REML=T, na.action = na.omit)
    
    anova(c_mod1E, c_mod1Ec, test=T)  # SIGN. ChiSq = 964.24, p<2.2e-16
    
  # Final Model 
    
    c_mod_final <- c_mod1E
    anova(c_mod_final)
        
  # Diagnostic plots to see if data follows model's assumptions
    
    plot(c_mod_final, resid(.,scaled=T)~fitted(.), abline=0, cex=0.5)
    hist(residuals(c_mod_final), main="Histogram of Residuals")
    qqmath(resid(c_mod_final, scaled=T), abline=c(0,1), strip=T)
    qqmath(ranef(c_mod_final, condVar = T), strip = T, abline=c(0,1))$id
    
    
#### 1.5. Plotting Effect Sizes for Final Models ####

# Effect Size from Heating Final Model
    
    size_effect_h = effectsize(h_mod_final)
    
    colnames(size_effect_h) = c('Parameters','Coefficient Values','CI','CI_low','CI_high')
    
    size_effect_h$Parameters[1:46]<-c("Intercept", "Time^1", "Time^2", "Thigmothermy", "Eye", "Head", "Leg", "Foot", "Snout", "Tail",
                                      "Excretion", "Mass", "Time^1 : Thigmothermy", "Time^2 : Thigmothermy",
                                      "Time^1 : Eye", "Time^2 : Eye", "Time^1 : Head", "Time^2 : Head", "Time^1 : Leg", 
                                      "Time^2 : Leg", "Time^1 : Foot", "Time^2 : Foot", "Time^1 : Snout", "Time^2 : Snout", 
                                      "Time^1 : Tail", "Time^2 : Tail", "Thigmothermy : Eye", "Thigmothermy : Head", 
                                      "Thigmothermy : Leg", "Thigmothermy : Foot", "Thigmothermy : Snout", "Thigmothermy : Tail",
                                      "Time^1 : Mass", "Time^2 : Mass", "Time^1 : Thigmothermy : Eye", "Time^2 : Thigmothermy : Eye", 
                                      "Time^1 : Thigmothermy : Head", "Time^2 : Thigmothermy : Head",
                                      "Time^1 : Thigmothermy : Leg", "Time^2 : Thigmothermy : Leg",
                                      "Time^1 : Thigmothermy : Foot", "Time^2 : Thigmothermy : Foot",
                                      "Time^1 : Thigmothermy : Snout", "Time^2 : Thigmothermy : Snout",
                                      "Time^1 : Thigmothermy : Tail", "Time^2 : Thigmothermy : Tail")
    
    Plot_EfS_h <- ggplot(size_effect_h, aes(`Coefficient Values`, reorder(Parameters, abs(desc(`Coefficient Values`))),
                                            xmin=size_effect_h$`Coefficient Values` - size_effect_h$CI_low,
                                            xmax = size_effect_h$`Coefficient Values` + size_effect_h$CI_high)) + 
                          geom_point() +  geom_errorbarh(height=.2) + coord_cartesian(xlim = c(-150,150)) +
                          xlab("\n Effect Size") + ylab("Variable(s)\n") + ggtitle("HEATING \n") +
                          theme_classic() + theme(plot.title = element_text(hjust = 0.5, size=20, face="bold")) +
                          geom_vline(xintercept = 0,linetype="dashed") +
                          geom_vline(xintercept = c(-2,-1,1,2,3,4),linetype="dotted") 
    
    
# Effect Size from Cooling Final Model
    
    size_effect_c = effectsize(c_mod_final)
    
    colnames(size_effect_c) = c('Parameters','Coefficient Values','CI','CI_low','CI_high')
    
    size_effect_c$Parameters[1:47]<-c("Intercept", "Time^1", "Time^2", "Thigmothermy", "Eye", "Head", "Leg", "Foot", "Snout", "Tail",
                                      "Excretion", "Posture", "Mass", "Time^1 : Thigmothermy", "Time^2 : Thigmothermy",
                                      "Time^1 : Eye", "Time^2 : Eye", "Time^1 : Head", "Time^2 : Head", "Time^1 : Leg", 
                                      "Time^2 : Leg", "Time^1 : Foot", "Time^2 : Foot", "Time^1 : Snout", "Time^2 : Snout", 
                                      "Time^1 : Tail", "Time^2 : Tail", "Thigmothermy : Eye", "Thigmothermy : Head", 
                                      "Thigmothermy : Leg", "Thigmothermy : Foot", "Thigmothermy : Snout", "Thigmothermy : Tail",
                                      "Time^1 : Mass", "Time^2 : Mass", "Time^1 : Thigmothermy : Eye", "Time^2 : Thigmothermy : Eye", 
                                      "Time^1 : Thigmothermy : Head", "Time^2 : Thigmothermy : Head",
                                      "Time^1 : Thigmothermy : Leg", "Time^2 : Thigmothermy : Leg",
                                      "Time^1 : Thigmothermy : Foot", "Time^2 : Thigmothermy : Foot",
                                      "Time^1 : Thigmothermy : Snout", "Time^2 : Thigmothermy : Snout",
                                      "Time^1 : Thigmothermy : Tail", "Time^2 : Thigmothermy : Tail")
    
    Plot_EfS_c <- ggplot(size_effect_c, aes(`Coefficient Values`, reorder(Parameters, abs(desc(`Coefficient Values`))),
                                            xmin=size_effect_c$`Coefficient Values` - size_effect_c$CI_low,
                                            xmax = size_effect_c$`Coefficient Values` + size_effect_c$CI_high)) + 
                        geom_point() + geom_errorbarh(height=.2) + coord_cartesian(xlim = c(-150,150)) +
                        xlab("\n Effect Size") + ylab("Variable(s)\n") + ggtitle("COOLING \n") +
                        theme_classic() + theme(plot.title = element_text(hjust = 0.5, size=20, face="bold")) +
                        geom_vline(xintercept = 0,linetype="dashed") +
                        geom_vline(xintercept = c(-2,-1,1,2,3,4),linetype="dotted")
    
    
# Visualizing the Effect Sizes for Both the Heating and Cooling Final Models:

    EfSs_plots <- grid.arrange(Plot_EfS_h, Plot_EfS_c, ncol=2)
    


    
###################################################################
###########       2. MAXIMUM RATES DATASET ANALYSIS       #########
###################################################################

    
#### 2.1. Data Importing and Wrangling ####
    
 #  This data set already has the data filtered to times between 60 and 120 seconds (re-named to time_2, 
 #  representing change in time from time 60, thus where time=60 is time_2=0) and temperature calculated 
 #  as cumulative temperature change between temperature at time x minus temperature at time = 60s, thus 
 #  making for equivalent/ comparable starting points for all lines 
    
  max_rt_data <- read.csv("max_rates_data.csv", header = T) 
    
  max_rt_data <- max_rt_data %>% filter(temp!="NA")  # getting rid of rows without Temperature values (i.e.NAs)
  
  max_rt_data <- max_rt_data %>% filter(time_2!=0)  # excludes time_2=0 which was forcing the y-intercept to zero
    
  { max_rt_data$treat <- as.factor(max_rt_data$treat)
    max_rt_data$id <- as.factor(max_rt_data$id)
    max_rt_data$population <- as.factor(max_rt_data$population)
    max_rt_data$heat_2 <- as.factor(max_rt_data$heat_2)
    max_rt_data$body_part <- as.factor(max_rt_data$body_part)
    max_rt_data$urine <- as.factor(max_rt_data$urine)
    max_rt_data$posture <- as.factor(max_rt_data$posture) 
    }
  str(max_rt_data)
  
  # Separating the Cooling and Heating Datasets:
  
    max_rt_data_H <- max_rt_data[max_rt_data$heat_2=="Heating",]
    max_rt_data_H$heat_2 <- droplevels(max_rt_data_H$heat_2)
  
    max_rt_data_C <- max_rt_data[max_rt_data$heat_2=="Cooling",]
    max_rt_data_C$heat_2 <- droplevels(max_rt_data_C$heat_2)
  
  
#### 2.2. Fitting OLS Models to Calculate Maximum Heat Exchange Rates ####
  
 # Ordinary Least Squares (OLS) models were fit for every Body Part, at every Treatment (Heliothermy vs Thigmothermy)
 # and for every Heat Exchange Process (Heating vs Cooling) to obtain the linear relationship between Time_2 (in 
 # seconds after 60s) and Cumulative Change in Temperature and hence calculate the maximum heat exchange rate (i.e. gradient)
  
 { l <- 2
   i <- 0
   j <- 0
   k <- 0
  
   bp <- c("Eye","Snout","Head","Dorsum","Leg","Foot","Tail")
   tt <- c("Heliothermy", "Thigmothermy")
   ht <- c("Cooling","Heating")
  
   mat <- matrix(ncol=10,nrow=29)
   mat[1,] <- c("Body_Part", "Treatment", "Heat", "Slope", "Slope_2.5%_CI", "Slope_97.5%_CI", 
                "y_Intercept", "Intercept_2.5%_CI", "Intercept_97.5%_CI", "R-squared")
  
   for(i in bp) {
     for(j in tt) {
       for(k in ht) {
        
         dt <- max_rt_data %>% filter(body_part==i) %>% filter(treat==j) %>% filter(heat_2==k)
        
         mod <- lmodel2(temp~time_2, data = dt, nperm=999, 
                        range.y="interval",range.x="interval")
        
         mat[l,1] <- i
         mat[l,2] <- j
         mat[l,3] <- k
         mat[l,4] <- mod$regression.results[1,3]
         mat[l,5] <- mod$confidence.intervals[1,4]
         mat[l,6] <- mod$confidence.intervals[1,5]
         mat[l,7] <- mod$regression.results[1,2]
         mat[l,8] <- mod$confidence.intervals[1,2]
         mat[l,9] <- mod$confidence.intervals[1,3]
         mat[l,10] <- mod$rsquare
        
         l <- l + 1
       }
     }
    }
   }
  
  # Cleaning up the models' outputs into a summary table 
  
    output_table <- as.data.frame(mat[2:29,])
  
    colnames(output_table) <- c("Body Part", "Treatment", "Heat", "Slope", "Slope 2.5% CI", "Slope 97.5% CI", 
                                "Y-Intercept", "Intercept 2.5% CI", "Intercept 97.5% CI", "R-squared")
    
    { output_table$`Body Part` <- as.factor(output_table$`Body Part`)
      output_table$Treatment <- as.factor(output_table$Treatment)
      output_table$Heat <- as.factor(output_table$Heat)
      output_table$Slope <- as.numeric(as.character(output_table$Slope))
      output_table$`Slope 2.5% CI` <- as.numeric(as.character(output_table$`Slope 2.5% CI`))
      output_table$`Slope 97.5% CI` <- as.numeric(as.character(output_table$`Slope 97.5% CI`))
      output_table$`Y-Intercept` <- as.numeric(as.character(output_table$`Y-Intercept`))
      output_table$`Intercept 2.5% CI` <- as.numeric(as.character(output_table$`Intercept 2.5% CI`))
      output_table$`Intercept 97.5% CI` <- as.numeric(as.character(output_table$`Intercept 97.5% CI`))
      output_table$`R-squared` <- as.numeric(as.character(output_table$`R-squared`))
      }
  
    Sup_Table1 <- output_table
    View(Sup_Table1)
      

#### 2.3. Plotting the Data ####
    
 #  Note that ABSOLUTE (i.e. the modulus of the calculated value) temperature changes were used in some of the plots,
 #  instead of the direct temperature change values (which would render negative slopes for the cooling
 #  trials) in order to facilitate the visual inspection and comparison between heating and cooling slopes


  # Plotting the Cumulative Absolute Change in Temperature for Each Body Part per Heat Treatment       #
  # (Heliothermy vs Thigmothermy) and per Heat Exchange Process (Heating vs Cooling)                   #
    
    names3 <- c("Heliothermy"="HELIOTHERMY", "Thigmothermy"="THIGMOTHERMY", "Cooling"="COOLING", "Heating"="HEATING")
    
    Mx_RT_Plot1 <- ggplot(max_rt_data, aes(x=time_2, y=abs(temp), group=body_part, colour=body_part)) +
      geom_jitter(size=0.5) + 
      stat_smooth(method="lm", formula = abs(y) ~ x, se = T, alpha = 0.1, linewidth = 1.1, fullrange=TRUE) +
      facet_grid(cols = vars(treat), rows = vars(factor(heat_2, levels=c('Heating','Cooling'))), scales="free", 
                 margins=F, labeller = as_labeller(names3)) +
      scale_color_manual(values=c(Dorsum="cornflowerblue", Eye="forestgreen", Head="firebrick2", 
                                  Leg="burlywood4", Foot="orchid3", Snout="paleturquoise3", Tail="gold1")) +
      xlab("\n Time After 60s (s) \n") + ylab("\n Cumulative Absolute  ∆ Temperature (ºC) \n")  + 
      labs(colour="    Body Part:  ") + 
      geom_hline(aes(yintercept=-Inf), linewidth=0.9) + geom_vline(aes(xintercept=-Inf), linewidth=0.9) + 
      scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
      coord_cartesian(xlim = c(0,70), ylim = c(0,3.5)) + expand_limits(x=c(0,68), y=c(0, 3))  +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
            axis.text.y = element_text(size = 11), axis.title = element_text(size=14),
            strip.text.x = element_text(size = 16), strip.text.y = element_text(size = 16)) +
      theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
            axis.line = element_blank(), panel.border = element_blank(), panel.spacing = unit(1, "cm"),
            strip.background = element_rect(color = "white", linewidth = 1)) +
      theme(legend.position = "bottom", legend.direction="horizontal", legend.title=element_text(size=14), 
            legend.text=element_text(size=12, margin = margin(r = 0.7, unit = 'cm'))) + guides(colour=guide_legend(nrow=1))         
    Mx_RT_Plot1
    
      
  # Plotting the Cumulative Absolute Change in Temperature for Each Heat Exchange Process (Heating vs        #
  # Cooling) per Heat Treatment (Heliothermy vs Thigmothermy) and per Body Part                              #
    
    names4 <- c("Heliothermy"="HELIOTHERMY", "Thigmothermy"="THIGMOTHERMY", "Dorsum"="DORSUM", "Eye"="EYE", 
                "Head"="HEAD", "Leg"="LEG", "Foot"="FOOT", "Snout"="SNOUT", "Tail"="TAIL")
    
    Mx_RT_Plot2 <- ggplot(max_rt_data, aes(x=time_2, y=abs(temp), group=heat_2, colour=heat_2)) +
      geom_jitter(size=0.5) +
      stat_smooth(method="lm", formula = abs(y) ~ x, se = T, alpha = 0.25, linewidth = 1.1, fullrange=TRUE) +
      scale_color_manual(labels = c("Cooling        ", "Heating"), values = c("steelblue2","tomato2")) + 
      facet_grid(cols = vars(treat), rows = vars(body_part), scales="free", margins=F, 
                 labeller = as_labeller(names4)) +
      xlab("\n Time After 60s (s)") + ylab("\n Cumulative Absolute  ∆ Temperature (ºC) \n")  + 
      labs(colour="Heat Exchange Process:   ") + 
      geom_hline(aes(yintercept=-Inf), linewidth=0.9) + geom_vline(aes(xintercept=-Inf), linewidth=0.9) + 
      scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
      coord_cartesian(xlim = c(0,70), ylim = c(0,3.5)) + expand_limits(x=c(0,68), y=c(0, 3))  +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
            axis.text.y = element_text(size = 11), axis.title = element_text(size=14),
            strip.text.x = element_text(size = 15), strip.text.y = element_text(size = 13)) +
      theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
            axis.line = element_blank(), panel.border = element_blank(), panel.spacing = unit(1, "line"),
            strip.background = element_rect(color = "white", linewidth = 1)) +
      theme(legend.position = "bottom", legend.direction="horizontal", legend.title=element_text(size=14), 
            legend.text=element_text(size=12)) + guides(colour=guide_legend(nrow=1))  
    Mx_RT_Plot2
    
       
  # Plotting the Absolute Slopes (with 2.5% and 97.5% CIs) for Each Body Part Under Each Treatment (Heliothermy     #
  # vs Thigmothermy) and for Each Heat Exchange Process (Heating vs Cooling)                                        #
    
    Mx_RT_Plot3 <- ggplot(output_table, aes(x=`Body Part`, y=abs(Slope), group=(Treatment:Heat), colour=Heat, shape=Treatment)) + 
      geom_point(position=position_dodge(0.5), size=2)+
      geom_errorbar(aes(ymin=abs(`Slope 2.5% CI`), ymax=abs(`Slope 97.5% CI`)), 
                    width=.2,position=position_dodge(0.5)) +
      scale_color_manual(labels = c("Cooling", "Heating"), values = c("steelblue2", "tomato2")) +
      scale_shape_manual(labels = c("Heliothermy", "Thigmothermy"), values = c(8,15)) +
      xlab("\n Body Part") + ylab("\n Absoulute Slope (∆ºC / second) \n") +
      labs(colour="Heat Exchange Process:   ", shape="Treatment:  ") +
      theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12), axis.title = element_text(size=14)) +
      theme(panel.border = element_blank(), axis.line = element_line(colour = "black"), 
            panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
      theme(legend.title=element_text(size=14), legend.text=element_text(size=12)) +
      theme(legend.position = "bottom", legend.box="horizontal", legend.direction="vertical",
            legend.title.align = 0, legend.justification = "right") + 
      guides(colour=guide_legend(nrow=2), shape=guide_legend(nrow=2), byrow = TRUE)
    Mx_RT_Plot3
        
    
#### 2.4. Fitting Mixed Effects Models and Model Selection ####
    
# For the Heating Dataset:
    
  # Full model
     
    h_max_mod1 <- lmer(temp ~ poly(time_2,2) * treat * body_part + urine + mass_initial*poly(time_2,2) + population + posture +
                         (1|id), data =  max_rt_data_H, REML=T, na.action = na.omit)
    
    anova(h_max_mod1)
        
  # Removing ID Random Effect
    
    h_max_mod2 <- lm(temp ~ poly(time_2,2) * treat * body_part + urine + mass_initial*poly(time_2,2) + population + posture, 
                     data =  max_rt_data_H,  na.action = na.omit)
    
    anova(h_max_mod1, h_max_mod2, test=T)  # SIGN. ChiSq = 4.2474, p=0.03931
      
  # Removing MASS_INITIAL:POLY(TIME,2) Interaction 
    
    h_max_mod1A <- lmer(temp ~ poly(time_2,2) * treat * body_part + urine + mass_initial + population + posture +
                          (1|id), data =  max_rt_data_H, REML=T, na.action = na.omit)
    
    anova(h_max_mod1, h_max_mod1A, test=T)  # NOT SIGN. ChiSq = 1.1876, p=0.5522
      
  # Removing TREAT from the POLY(TIME,2):TREAT:POSITION Interaction
    
    h_max_mod1Ab <- lmer(temp ~ poly(time_2,2) * body_part + treat + urine + mass_initial + population + posture +
                           (1|id), data =  max_rt_data_H, REML=T, na.action = na.omit)
    
    anova(h_max_mod1A, h_max_mod1Ab, test=T)  # SIGN. ChiSq = 136.4, p<2.2e-16
      
  # Removing POLY(TIME,2) from the POLY(TIME,2):TREAT:POSITION Interaction
    
    h_max_mod1Ac <- lmer(temp ~ poly(time_2,2) + body_part * treat + urine + mass_initial + population + posture +
                           (1|id), data =  max_rt_data_H, REML=T, na.action = na.omit)
    
    anova(h_max_mod1A, h_max_mod1Ac, test=T)  # NOT SIGN. ChiSq = 36.507, p=0.08273
     
  # Removing BODY PART from the BODY PART:TREAT Interaction
    
    h_max_mod1Aca <- lmer(temp ~ poly(time_2,2) + body_part + treat + urine + mass_initial + population + posture +
                            (1|id), data =  max_rt_data_H, REML=T, na.action = na.omit)
    
    anova(h_max_mod1Ac, h_max_mod1Aca, test=T)  # SIGN. ChiSq = 118.66, p<2.2e-16
       
  # Removing POPULATION
    
    h_max_mod1Acb <- lmer(temp ~ poly(time_2,2) + body_part * treat + urine + mass_initial + posture +
                            (1|id), data =  max_rt_data_H, REML=T, na.action = na.omit)
    
    anova(h_max_mod1Ac, h_max_mod1Acb, test=T)  # NOT SIGN. ChiSq = 2.2956, p=0.5134
      
  # Removing POSTURE
    
    h_max_mod1Acba <- lmer(temp ~ poly(time_2,2) + body_part * treat + urine + mass_initial + 
                             (1|id), data =  max_rt_data_H, REML=T, na.action = na.omit)
    
    anova(h_max_mod1Acb, h_max_mod1Acba, test=T)  # NOT SIGN. ChiSq = 5.3258, p=0.06975
      
  # Removing URINE
    
    h_max_mod1Acba1 <- lmer(temp ~ poly(time_2,2) + body_part * treat + mass_initial + 
                              (1|id), data =  max_rt_data_H, REML=T, na.action = na.omit)
    
    anova(h_max_mod1Acba, h_max_mod1Acba1, test=T)  # NOT SIGN. ChiSq = 0.0733, p=0.7867
      
  # Removing MASS_INITIAL
    
    h_max_mod1Acba1a <- lmer(temp ~ poly(time_2,2) + body_part * treat +  
                               (1|id), data =  max_rt_data_H, REML=T, na.action = na.omit)
    
    anova(h_max_mod1Acba1, h_max_mod1Acba1a, test=T)  # NOT SIGN. ChiSq = 1.7512, p=0.1857
      
  # Removing TIME^2
    
    h_max_mod1Acba1aA <- lmer(temp ~ poly(time_2,1) + body_part * treat +  
                                (1|id), data =  max_rt_data_H, REML=T, na.action = na.omit)
    
    anova(h_max_mod1Acba1a, h_max_mod1Acba1aA, test=T)  # NOT SIGN. ChiSq = 0.0508, p=0.8216
       
  # Removing TIME^1
    
    h_max_mod1Acba1aAa <- lmer(temp ~ body_part * treat +  
                                 (1|id), data =  max_rt_data_H, REML=T, na.action = na.omit)
    
    anova(h_max_mod1Acba1aA, h_max_mod1Acba1aAa, test=T)  # SIGN. ChiSq = 494.26, p<2.2e-16
       
  # Final Model for Heating Dataset
    
    h_max_mod_final <- h_max_mod1Acba1aA
    anova(h_max_mod_final)
      
  # Diagnostic plots to see if data follows model's assumptions
    
    plot(h_max_mod_final, resid(.,scaled=T)~fitted(.), abline=0, cex=0.5)
    hist(residuals(h_max_mod_final), main = "Histogram of Residuals")
    qqmath(resid(h_max_mod_final, scaled=T), abline=c(0,1), strip=T)
    qqmath(ranef(h_max_mod_final, condVar = T), strip = T, abline=c(0,1))$id


# For the Cooling Dataset:
    
  # Full model 
    
    c_max_mod1 <- lmer(temp ~ poly(time_2,2) * treat * body_part + urine + mass_initial*poly(time_2,2) + population + posture +
                         (1|id), data =  max_rt_data_C, REML=T, na.action = na.omit)  
    anova(c_max_mod1)
       
  # Removing ID Random Effect
    
    c_max_mod2 <- lm(temp ~ poly(time_2,2) * treat * body_part + urine + mass_initial*poly(time_2,2) + population + posture, 
                     data =  max_rt_data_C,  na.action = na.omit)
    
    anova(c_max_mod1, c_max_mod2, test=T)  # SIGN. ChiSq = 11.733, p=0.0006141
      
  # Removing MASS_INITIAL:POLY(TIME,2) Interaction 
    
    c_max_mod1A <- lmer(temp ~ poly(time_2,2) * treat * body_part + urine + mass_initial + population + posture +
                          (1|id), data =  max_rt_data_C, REML=T, na.action = na.omit)
    
    anova(c_max_mod1, c_max_mod1A, test=T)  # NOT SIGN. ChiSq = 4.5917, p=0.1007
      
  # Removing TREAT from the POLY(TIME,2):TREAT:POSITION Interaction
    
    c_max_mod1Ab <- lmer(temp ~ poly(time_2,2) * body_part + treat + urine + mass_initial + population + posture +
                           (1|id), data =  max_rt_data_C, REML=T, na.action = na.omit)
    
    anova(c_max_mod1A, c_max_mod1Ab, test=T)  # SIGN. ChiSq = 44.584, p=0.001256
       
  # Removing POLY(TIME,2) from the POLY(TIME,2):TREAT:POSITION Interaction
    
    c_max_mod1Ac <- lmer(temp ~ poly(time_2,2) + body_part * treat + urine + mass_initial + population + posture +
                           (1|id), data =  max_rt_data_C, REML=T, na.action = na.omit)
    
    anova(c_max_mod1A, c_max_mod1Ac, test=T)  # NOT SIGN. ChiSq = 36.303, p=0.08624
      
  # Removing BODY PART from the BODY PART:TREAT Interaction
    
    c_max_mod1Aca <- lmer(temp ~ poly(time_2,2) + body_part + treat + urine + mass_initial + population + posture +
                            (1|id), data =  max_rt_data_C, REML=T, na.action = na.omit)
    
    anova(c_max_mod1Ac, c_max_mod1Aca, test=T)  # SIGN. ChiSq = 27.931, p = 9.682e-05
       
  # Removing POPULATION
    
    c_max_mod1Acb <- lmer(temp ~ poly(time_2,2) + body_part * treat + urine + mass_initial + posture +
                            (1|id), data =  max_rt_data_C, REML=T, na.action = na.omit)
    
    anova(c_max_mod1Ac, c_max_mod1Acb, test=T)  # SIGN. ChiSq = 10.106, p=0.01768
      
  # Removing POSTURE
    
    c_max_mod1Acc <- lmer(temp ~ poly(time_2,2) + body_part * treat + urine + mass_initial + population +
                            (1|id), data =  max_rt_data_C, REML=T, na.action = na.omit)
    
    anova(c_max_mod1Ac, c_max_mod1Acc, test=T)  # NOT SIGN. ChiSq = 2.5185, p=0.2839
       
  # Removing URINE
    
    c_max_mod1Acca <- lmer(temp ~ poly(time_2,2) + body_part * treat + mass_initial + population +
                             (1|id), data =  max_rt_data_C, REML=T, na.action = na.omit)
    
    anova(c_max_mod1Acc, c_max_mod1Acca, test=T)  # NOT SIGN. ChiSq = 2.4843, p=0.115
      
  # Removing MASS_INITIAL
    
    c_max_mod1Acca1 <- lmer(temp ~ poly(time_2,2) + body_part * treat + population + 
                              (1|id), data =  max_rt_data_C, REML=T, na.action = na.omit)
    
    anova(c_max_mod1Acca, c_max_mod1Acca1, test=T)  # NOT SIGN. ChiSq = 1.127, p=0.2884
    
    
  # Removing TIME^2
    
    c_max_mod1Acca1a <- lmer(temp ~ poly(time_2,1) + body_part * treat + population + 
                               (1|id), data =  max_rt_data_C, REML=T, na.action = na.omit)
    
    anova(c_max_mod1Acca1, c_max_mod1Acca1a, test=T)  # NOT SIGN. ChiSq = 0.0159, p=0.8996
      
  # Removing TIME^1
    
    c_max_mod1Acca1ab <- lmer(temp ~ body_part * treat + population + 
                                (1|id), data =  max_rt_data_C, REML=T, na.action = na.omit)
    
    anova(c_max_mod1Acca1a, c_max_mod1Acca1ab, test=T)  # SIGN. ChiSq = 360.33, p<2.2e-16
      
  # Final Model for Heating Dataset
    
    c_max_mod_final <- c_max_mod1Acca1a
    anova(c_max_mod_final)
    
  # Diagnostic plots to see if data follows model's assumptions
    
    plot(c_max_mod_final, resid(.,scaled=T)~fitted(.), abline=0, cex=0.5)
    hist(residuals(c_max_mod_final), main = "Histogram of Residuals")
    qqmath(resid(c_max_mod_final, scaled=T), abline=c(0,1), strip=T)
    qqmath(ranef(c_max_mod_final, condVar = T), strip = T, abline=c(0,1))$id  
    

#### 2.5. Plotting Effect Sizes for Final Maximum Rates Models ####
    
# Effect Size from Maximum Heating Rate Final Model
    
    size_effect_h_max = effectsize(h_max_mod_final)

    colnames(size_effect_h_max) = c('Parameters','Coefficient Values','CI','CI_low','CI_high')
    
    size_effect_h_max$Parameters[1:15]<-c("Intercept", "Time^1", "Eye", "Head", "Leg", "Foot", "Snout", "Tail", "Thigmothermy",
                                          "Thigmothermy : Eye", "Thigmothermy : Head", "Thigmothermy : Leg", 
                                          "Thigmothermy : Foot", "Thigmothermy : Snout", "Thigmothermy : Tail")
                                          
    
    Plot_EfS_h_max <- ggplot(size_effect_h_max,aes(`Coefficient Values`, reorder(Parameters, abs(desc(`Coefficient Values`))),
                                                    xmin=size_effect_h_max$`Coefficient Values` - size_effect_h_max$CI_low,
                                                    xmax = size_effect_h_max$`Coefficient Values` + size_effect_h_max$CI_high)) + 
                              geom_point() + geom_errorbarh(height=.2) + coord_cartesian(xlim = c(-35,35)) +
                              xlab("\n Effect Size") + ylab("Variable(s)\n") + ggtitle("HEATING (Maximum Rates) \n") +
                              theme_classic() + theme(plot.title = element_text(hjust = 0.5, size=20, face="bold")) +
                              geom_vline(xintercept = 0,linetype="dashed") +
                              geom_vline(xintercept = c(-2,-1,1,2,3,4),linetype="dotted")
    
    
# Effect Size from Cooling Final Model
    
    size_effect_c_max = effectsize(c_max_mod_final)
    
    colnames(size_effect_c_max) = c('Parameters','Coefficient Values','CI','CI_low','CI_high')
    
    size_effect_c_max$Parameters[1:18]<-c("Intercept", "Time^1", "Eye", "Head", "Leg", "Foot", "Snout", "Tail", "Thigmothermy",
                                          "Population = Évora", "Population = Portimão", "Population = Torres Vedras",
                                          "Thigmothermy : Eye", "Thigmothermy : Head", "Thigmothermy : Leg", 
                                          "Thigmothermy : Foot", "Thigmothermy : Snout", "Thigmothermy : Tail")
    
    Plot_EfS_c_max <- ggplot(size_effect_c_max,aes(`Coefficient Values`, reorder(Parameters, abs(desc(`Coefficient Values`))),
                                                   xmin=size_effect_c_max$`Coefficient Values` - size_effect_c_max$CI_low,
                                                   xmax = size_effect_c_max$`Coefficient Values` + size_effect_c_max$CI_high)) + 
                              geom_point() + geom_errorbarh(height=.2) + coord_cartesian(xlim = c(-35,35)) +
                              xlab("\n Effect Size") + ylab("Variable(s)\n") + ggtitle("COOLING (Maximum Rates) \n") +
                              theme_classic() + theme(plot.title = element_text(hjust = 0.5, size=20, face="bold")) +
                              geom_vline(xintercept = 0,linetype="dashed") +
                              geom_vline(xintercept = c(-2,-1,1,2,3,4),linetype="dotted")
    
    
    # Visualizing the Effect Sizes for Both the Heating and Cooling Final Models:
    
    EfSs_max_plots <- grid.arrange(Plot_EfS_h_max, Plot_EfS_c_max, ncol=2)
   
    ### THE END! ###
