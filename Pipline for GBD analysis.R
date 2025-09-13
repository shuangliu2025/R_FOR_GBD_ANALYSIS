# Pipline for GBD analysis 
###########part 1########################
library(data.table)
library(tidyverse)
library(gtools)
analysedir <- "F:/workforSCI/252/dataanalyse/GBD/results/NAFLD/over45/test/global/"
dir.create(analysedir)
diseasedata <- "F:/workforSCI/252/dataanalyse/GBD/diseasedata_all/NAFLD/NAFLD45.csv"
GBD_Maps_match_path <- "F:/workforSCI/252/dataanalyse/GBD/gbd_rnaturalearth_mapping.csv"

digestdata<- fread(input =diseasedata,sep = "auto",header="auto")

###########select and modify colnames
if (length(colnames(digestdata))>15) {
  digestdata <- digestdata %>% select(c(2,4,6,8,10,12,13:16))
  colnames(digestdata) <- c("measure","location","sex","age","cause","metric","year","val","upper","lower")
}
gc()
#####choose regions (lable AIMDISEASELABEL)
colnames(digestdata)
unique(digestdata$cause)   #recheck data
mixedsort(unique(digestdata$year))   
mixedsort(unique(digestdata$age))
unique(digestdata$location) 
unique(digestdata$measure) 
unique(digestdata$metric)
location_a <- as.vector(unique(digestdata$location))
income_c <- c("World Bank High Income","World Bank Low Income", "World Bank Lower Middle Income","World Bank Upper Middle Income")
sum(location_a %in% income_c)
SDI_c <- c("Low SDI" ,"Low-middle SDI","Middle SDI" ,"High-middle SDI","High SDI")
sum(location_a %in% SDI_c)
main_c <- c('Global','High-income North America','Australasia','High-income Asia Pacific','Western Europe',
            'Southern Latin America','Eastern Europe','Central Europe','Central Asia',
            'Central Latin America','Andean Latin America','Caribbean',
            'Tropical Latin America','East Asia','Southeast Asia','Oceania',
            'North Africa and Middle East','South Asia','Southern Sub-Saharan Africa',
            'Western Sub-Saharan Africa','Eastern Sub-Saharan Africa',
            'Central Sub-Saharan Africa') ### 22 countries and regions
sum(location_a %in% main_c)
location_c <- location_a[!(location_a %in% c(income_c,SDI_c,main_c))]  ##select regions
agegroups <- mixedsort(unique(digestdata$age))[c(1:11)]
agegroups
aimdisease <- unique(digestdata$cause)
aimdisease
aimdiseaselabel <- c("NAFLD45")
measures <- mixedsort(unique(digestdata$measure))[c(1:4)]
measures
measureslabel<- c("DALYs","Deaths" ,"Incidence","Prevalence")
###############one##########################
library(ggmap)
library(rgdal) 
library(maps)
library(dplyr)
library(rnaturalearth)
library(geodata)
year1=1990
year2=2021
year3=2035
#four indexes_INCIDENCE#######################
sexindex <- unique(digestdata$sex)
sexindex

##match regions name
world_ne <- ne_countries(scale = "medium", returnclass = "sf")
tokelau_gadm <- gadm(country = "TKL", level = 0, path = tempdir())
tokelau_sf <- st_as_sf(tokelau_gadm)
# 
tokelau_sf <- tokelau_sf %>%
  mutate(
    name_long = "Tokelau",
    iso_a3 = "TKL",
    iso_a2 = "TK",
    economy = "6. Developing region",
    income_grp = "4. Lower middle income",
    # 添加其他必要字段，设置为NA或适当值
    geounit = "Tokelau",
    sovereignt = "New Zealand",
    continent = "Oceania",
    region_un = "Oceania",
    subregion = "Polynesia"
  ) %>%
  select(any_of(names(world_ne))) # 只保留与world_ne相同的列
world_ne <- world_ne %>% select(any_of(names(tokelau_sf)))
world_ne <- rbind(world_ne, tokelau_sf)

##reading matching table
GBD_Maps_match <- read.csv(file =  GBD_Maps_match_path)

for (d in c(1:length(unique(digestdata$cause)))) {
  for (s in c(1:length(sexindex))) {
    for (m in c(1:length(measures))) {
      year_ARS <- digestdata%>%
        filter(digestdata$year == year2,
               digestdata$age == 'Age-standardized',  
               digestdata$metric == 'Rate', #（per 10^5）
               digestdata$measure == measures[m],  
               digestdata$cause==aimdisease[d],
               digestdata$sex == sexindex[s],
               digestdata$location %in% location_c 
        )
      year_ARS2 <- year_ARS%>% select(cause,measure,sex,location,year,val,upper,lower)
      year_ARS2$val<-round(year_ARS2$val,1)#
      year_ARS2$lower<-round(year_ARS2$lower,1)#
      year_ARS2$upper<-round(year_ARS2$upper,1) #
      year_ARS2$ASR_2021 <-paste(year_ARS2$val,'(',year_ARS2$lower,'-',year_ARS2$upper,')',sep='')
      #########countries in 2021—ASR
      write.csv(x = year_ARS2[order(year_ARS2$val),],
                file = paste(analysedir,"ASRs_of_",measureslabel[m],"_",aimdiseaselabel[d],"_All_Countries_",sexindex[s],"2021.csv",sep = ''),row.names = F)

      sum(GBD_Maps_match$rnaturalearth_Country %in% unique(world_ne$name_long)) ##count the matching number with map
      sum(GBD_Maps_match$GBD_Country %in% unique(year_ARS2$location)) ##count the matching number with gbd data
      
      year_ARS2 <- left_join(x = year_ARS2,y = GBD_Maps_match,by = c("location"="GBD_Country"))
      total <- world_ne[world_ne$name_long %in% year_ARS2$rnaturalearth_Country,] 
      total <- full_join (total, year_ARS2, by =c('name_long'='rnaturalearth_Country')) 
      cat(paste("map data has",length(unique(total$name_long)),'regions names'))
      #split the range into seven parts
      range(total$val)
      label15 <- round((max(total$val)-min(total$val))*0.15+min(total$val),1)
      label30 <- round((max(total$val)-min(total$val))*0.30+min(total$val),1)
      label45 <- round((max(total$val)-min(total$val))*0.45+min(total$val),1)
      label60 <- round((max(total$val)-min(total$val))*0.60+min(total$val),1)
      label75 <- round((max(total$val)-min(total$val))*0.75+min(total$val),1)
      label90 <- round((max(total$val)-min(total$val))*0.90+min(total$val),1)
      
      rangelabels <- c(paste0("<",label15),paste0(label15,"<",label30),paste0(label30,"<",label45),
                       paste0(label45,"<",label60),paste0(label60,"<",label75),paste0(label75,"<",label90),paste0(">=",label90))
      #chunks <- split(mixedsort(unique(total$val)), cut(seq_along(mixedsort(unique(total$val))), 5 , labels= FALSE ))
      total2 <- total %>% mutate(val_label = cut(val, breaks = c(-Inf,label15,label30,label45,label60,label75,label90,Inf),
                                                 labels = rangelabels, 
                                                 include.lowest = T,right = T))
      
      COLORtake <- colorRampPalette(c("#99CCFF","#CCFFFF","#FFFFCC",'#FFCC99',"#CC3333"))
      #cc=RColorBrewer::brewer.pal(nrow(unique(total_case2["val_label"])),name = "Set1") ##分类变量取色
      cc <- COLORtake(length(unique(total2$val_label))) 
      scales::show_col(cc)
      
      p2 <- ggplot(data=total2,)+
       geom_sf (aes(fill=val_label), 
                             colour="black", 
                linewidth =.05) +
        #scale_fill_gradient( low = "#66CD00", high = "#FF4040")+ ###Coloring of continuous variables
        scale_fill_manual(values =cc,
                          na.value = "grey90")+ ##coloring of categorical variables 
        theme_void()+ 
        xlab(label = "")+
        ylab(label = "")+
        labs(title = paste0("ASRs_of_",measureslabel[m],"_",aimdiseaselabel[d],"_All_Countries_",sexindex[s]," 2021"),fixed = T,pattern = ".txt",replacement = "")+
        theme(plot.title = element_text(hjust = 0.5))+  
        guides(fill = guide_legend(title='Cases (/10^5)'))+ 
        #theme_bw(base_rect_size = 0.1,base_line_size = 0.1) +
        theme(
          legend.text = element_text(size = 16),   
          legend.title = element_text(size = 18)  
        )
      p2 
      ggsave(plot = p2,filename = paste(analysedir,"ASRs_of_",measureslabel[m],"_",aimdiseaselabel[d],"_All_Countries_",sexindex[s],"2021.png",sep = ''),
             width = 14,height =7)
      ggsave(plot = p2,filename = paste(analysedir,"ASRs_of_",measureslabel[m],"_",aimdiseaselabel[d],"_All_Countries_",sexindex[s],"2021.pdf",sep = ''),
             width = 14,height =7)
    }
  }
  
  ##############2、cases_to_agegroup_by_sex##################
  
  for (m in c(1:length(measures))) {
    
    cases_agegroup_both <- digestdata%>%
      filter(digestdata$year == year2,
             digestdata$age %in% agegroups, 
             digestdata$metric == 'Number', 
             digestdata$measure == measures[m],   
             digestdata$cause== aimdisease[d],
             digestdata$sex == "Both",
             digestdata$location == "Global" 
      ) %>% select(cause,measure,age,sex,
                   location,year,val,upper,lower)
    
    cases_agegroup_both[,7:9] <- cases_agegroup_both[,7:9]/100000
    cases_agegroup_both$val<-round(cases_agegroup_both$val,digits = 2)#choose zhe data precision 
    cases_agegroup_both$lower<-round(cases_agegroup_both$lower,2)#
    cases_agegroup_both$upper<-round(cases_agegroup_both$upper,2) #
    cases_agegroup_both$case_2021 <-paste(cases_agegroup_both$val,'(',cases_agegroup_both$lower,'-',cases_agegroup_both$upper,')',sep='')
    #male_and_female
    cases_agegroup_maf <- digestdata%>%
      filter(digestdata$year == year2,
             digestdata$age %in% agegroups,  #
             digestdata$metric == 'Number', #率（per 10^5）
             digestdata$measure == measures[m],  
             digestdata$cause== aimdisease[d],
             digestdata$sex %in% c("Male","Female"),
             digestdata$location == "Global" #selected region
      ) %>%
      select(cause,measure,age,sex,
             location,year,val,upper,lower)
    
    cases_agegroup_maf[,7:9] <- cases_agegroup_maf[,7:9]/100000
    cases_agegroup_maf$val   <- round(cases_agegroup_maf$val,digits = 2)
    cases_agegroup_maf$lower <- round(cases_agegroup_maf$lower,2)#
    cases_agegroup_maf$upper <- round(cases_agegroup_maf$upper,2) #
    cases_agegroup_maf$case_2021 <-paste(cases_agegroup_maf$val,'(',cases_agegroup_maf$lower,'-',cases_agegroup_maf$upper,')',sep='')
    
    #########combining data
    cases_agegroup_all <- merge(x = cases_agegroup_both,y = cases_agegroup_maf[sex=="Male",c(3,4,7:9)],by.x = "age",by.y = "age")
    cases_agegroup_all <- merge(x = cases_agegroup_all,y = cases_agegroup_maf[sex=="Female",c(3,4,7:9)],by.x = "age",by.y = "age")
    cases_agegroup_all <- cases_agegroup_all%>%select(-c(4,10,11,15))%>%
      rename("val_B"="val.x","upper_B"="upper.x","lower_B"="lower.x",
             "val_M"="val.y","upper_M"="upper.y","lower_M"="lower.y",
             "val_F"="val","upper_F"="upper","lower_F"="lower")
    cases_agegroup_all <- cases_agegroup_all[mixedsort(cases_agegroup_all$age),]
    
    write.csv(x = cases_agegroup_all, file = paste(analysedir,"Number of ",measureslabel[m]," of ",aimdiseaselabel[d],"for all by age_global_2021.csv",sep = ''),row.names = F)
    
    ####plot the graph
    COLORtake <- colorRampPalette(c("#99CCFF","#CC3333"))
    cc <- COLORtake(nrow(unique(cases_agegroup_maf$sex))) 
    scales::show_col(cc)
    
    agegroups <- mixedsort(unique(cases_agegroup_both$age))
    if (length(agegroups)>=21) {
      ageage <- agegroups[c(21,1:20)]
    }else{
      ageage <- agegroups
    }
    ageage
    
    #reorder the data####
    cases_agegroup_all <- cases_agegroup_all[match(ageage, cases_agegroup_all$age),]
    
    p <- ggplot(data = cases_agegroup_all)+
      geom_bar(aes(x=factor(cases_agegroup_all$age,levels = ageage,ordered = T),y= val_B),
               position="dodge",stat="identity",fill = "grey")+
      geom_errorbar(aes(x=factor(cases_agegroup_all$age,levels = ageage,ordered = T),y=val_B,
                        ymax = upper_B, ymin = lower_B),
                    position = position_dodge(0.9), width = 0.15)+
      geom_line (data = cases_agegroup_maf,
                 aes(x=factor(age,levels = ageage,ordered = T), y=val,
                     group=sex,color=sex), linewidth = 0.8) +
      geom_ribbon(data = cases_agegroup_maf,
                  aes(x=factor(age,levels = ageage,ordered = T), y=val,
                      ymin=lower, ymax=upper,group =sex, fill = sex), alpha=0.3, 
                  linetype = "blank")+
      theme_bw()+
      xlab(label = "Age groups")+
      ylab(label = "Cases per 100 000")+
      labs(title = paste0("Number of ",measureslabel[m]," of ",aimdiseaselabel[d],"for by age_global_2021"),fixed = T,pattern = ".txt",replacement = "")+
      theme(plot.title = element_text(hjust = 0.5))+  ##
      theme(legend.position ='right')
    p3 <- p+theme(axis.text.x = element_text(size = 8, family = "myFont", color = "black",
                                             face = "plain", vjust = 0.5, hjust = 0.5, angle = 30)) 
    p3  
    ggsave(plot = p3,filename = paste(analysedir,"Number of ",measureslabel[m]," of ",aimdiseaselabel[d],"for all by age_global_2021.png",sep = ''),
           width = 14,height =7)
    
  }
}

###############################3 calculate the EAPC values and plot the graph###################
sexindex <- unique(digestdata$sex)
for (d in c(1:length(unique(digestdata$cause)))) {
  for (m in c(1:length(unique(measures)))) {
    for (s in sexindex) {
      EAPC_data <- digestdata%>%
        filter(digestdata$year %in% c(year1:year2),
               digestdata$age == 'Age-standardized',  
               digestdata$metric == 'Rate', 
               digestdata$measure == 'Incidence',  
               digestdata$cause== aimdisease[d] ,
               digestdata$sex == s,
               location %in% location_c 
        )
      length(unique(EAPC_data$location))  
      EAPC_2021 <- data.table() #
      for (i in unique(EAPC_data$location)) {
        temp <- EAPC_data[location == i, ]
        temp$lograte <- log(temp$val)
        fit <- lm(lograte ~ year, temp)
        b <- coefficients(fit)
        CI <- confint(fit)
        EAPC0.5 <- format(round(100 * (exp(b[2]) - 1), digits = 2), nsmall = 2)
        EAPC0.975 <- format(round(100 * (exp(CI[4]) - 1), digits = 2), nsmall = 2)
        EAPC0.025 <- format(round(100 * (exp(CI[2]) - 1), digits = 2), nsmall = 2)
        #
        location <- EAPC_data$location[EAPC_data$location == i] %>% unique()
        EAPC_t <- cbind(location, EAPC0.025, EAPC0.5, EAPC0.975)
        EAPC_2021 <- rbind(EAPC_2021,EAPC_t)
      }
      EAPC_2021 <- EAPC_2021[,.(location, EAPC= EAPC_2021$EAPC0.5,EAPC_0.05_0.975 = paste(EAPC_2021$EAPC0.5,' (', 
                                                                                          EAPC_2021$EAPC0.025, ' to ', 
                                                                                          EAPC_2021$EAPC0.975, ')', sep = ''))]
      
      write.csv(x = EAPC_2021[order(EAPC_2021$EAPC),],
                file = paste(analysedir,"EAPC_of_Incidence_of_",aimdiseaselabel[d],"_",s,"2021.csv",sep = ''),row.names = F)
      
      EAPC_2021<- left_join(x = EAPC_2021,y = GBD_Maps_match,by = c("location"="GBD_Country"))
      total_EAPC <- world_ne[world_ne$name_long %in% EAPC_2021$rnaturalearth_Country,] 
      total_EAPC <- full_join (total_EAPC, EAPC_2021, by =c('name_long'='rnaturalearth_Country')) 
      cat(paste("map data has",length(unique(total_EAPC$name_long)),'regions'))
      total_EAPC$EAPC <- as.numeric(total_EAPC$EAPC)
      summary(total_EAPC$EAPC)#
      range(total_EAPC$EAPC)
      label15 <- round((max(total_EAPC$EAPC)-min(total_EAPC$EAPC))*0.15+min(total_EAPC$EAPC),2)
      label30 <- round((max(total_EAPC$EAPC)-min(total_EAPC$EAPC))*0.30+min(total_EAPC$EAPC),digits = 2)
      label45 <- round((max(total_EAPC$EAPC)-min(total_EAPC$EAPC))*0.45+min(total_EAPC$EAPC),digits = 2)
      label60 <- round((max(total_EAPC$EAPC)-min(total_EAPC$EAPC))*0.60+min(total_EAPC$EAPC),digits = 2)
      label75 <- round((max(total_EAPC$EAPC)-min(total_EAPC$EAPC))*0.75+min(total_EAPC$EAPC),digits = 2)
      label90 <- round((max(total_EAPC$EAPC)-min(total_EAPC$EAPC))*0.90+min(total_EAPC$EAPC),digits = 2)
      
      rangelabels <- c(paste0("<",label15),paste0(label15,"<",label30),paste0(label30,"<",label45),
                       paste0(label45,"<",label60),paste0(label60,"<",label75),paste0(label75,"<",label90),paste0(">=",label90))
      
      total2_EAPC <- total_EAPC %>% mutate(EAPC_label = cut(EAPC, breaks = c(-Inf,label15,label30,label45,label60,label75,label90,Inf),
                                                            labels = rangelabels, 
                                                            include.lowest = T,right = T))
      
      
      p2 <- ggplot(data=total_EAPC)+ 
      geom_sf (aes(fill=EAPC), 
                             colour="black", linewidth =.05) +
        scale_fill_gradient2( low = '#005792',mid='white', high = "#F95959")+ 
        theme_void()+ ##transparent background 
        xlab(label = "")+
        ylab(label = "")+
        labs(title = paste0('EAPC_in_',measureslabel[m]," of ",aimdiseaselabel[d],"_in_",s,'_global'),fixed = T,pattern = ".txt",replacement = "")+
        theme(plot.title = element_text(hjust = 0.5))+  #
        guides(fill = guide_legend(title='EAPC'))+ 
        #theme_bw(base_rect_size = 0.1,base_line_size = 0.1) +
        theme(
          legend.text = element_text(size = 16),   #
          legend.title = element_text(size = 18)   # 
        )
      p2 
      ggsave(plot = p2,filename = paste(analysedir,"EAPC_in_",measureslabel[m]," of ",aimdiseaselabel[d],"_in_",s,"_global.png",sep = ''),
             width = 14,height =7)
      ggsave(plot = p2,filename = paste(analysedir,"EAPC_in_",measureslabel[m]," of ",aimdiseaselabel[d],"_in_",s,"_global.pdf",sep = ''),
             width = 14,height =7)
    }
  }  
}
gc()

###########part 2###########################
library(dplyr)
library(data.table)
library(tidyverse)
library(gtools)
library(scales)

resultdir <- "F:/workforSCI/252/dataanalyse/GBD/results/NAFLD/over45/test/mainregions/"
dir.create(resultdir)
year1<- 1990
year2 <- 2021
year3 <- 2035
aimyears <- c(year1:year2)

##########################1、ASR############################################
index_data_both <- data.frame(location=c(SDI_c,main_c))

for (m in c(1:length(measures))) {
  index_data_both <- data.frame(location=c(SDI_c,main_c))
  for (d in c(1:length(unique(digestdata$cause)))) {
    #for (s in sex_index) {
    incidence_year1 <- digestdata%>%
      filter(digestdata$year == year1,
             digestdata$age == 'Age-standardized',  
             digestdata$metric == 'Rate', #
             digestdata$measure == measures[m], 
             digestdata$cause == aimdisease[d] ,
             digestdata$sex == "Both",
             location %in% c(main_c,SDI_c) #
      )
    
    incidence_year2 <- digestdata%>%
      filter(digestdata$year == year2,
             digestdata$age == 'Age-standardized',  #
             digestdata$metric == 'Rate', #
             digestdata$measure == measures[m],  #
             digestdata$cause== aimdisease[d] ,
             digestdata$sex == "Both",
             location %in% c(main_c,SDI_c) #
      )
  
    incidence_year1$val <- round(incidence_year1$val,2)
    incidence_year1$upper <- round(incidence_year1$upper,2)
    incidence_year1$lower <- round(incidence_year1$lower,2)
    incidence_year1$ASR_year1 <- paste(incidence_year1$val,'(',incidence_year1$lower,'-',incidence_year1$upper,')',sep='')
    
    incidence_year2$val <- round(incidence_year2$val,2)
    incidence_year2$upper <- round(incidence_year2$upper,2)
    incidence_year2$lower <- round(incidence_year2$lower,2)
    incidence_year2$ASR_year2 <- paste(incidence_year2$val,'(',incidence_year2$lower,'-',incidence_year2$upper,')',sep='')
    
    index_data_both <- merge(x = index_data_both,y = incidence_year1[,c("location","ASR_year1")],by.x = 'location',by.y="location")
    index_data_both <- merge(x = index_data_both,y = incidence_year2[,c("location","ASR_year2")],by.x = 'location',by.y="location")
    
    name1<-paste0("ASR_",aimdiseaselabel[d],"_",year1)
    name2<- paste0("ASR_",aimdiseaselabel[d],"_",year2)
    index_data_both <- plyr::rename(index_data_both,c("ASR_year1"= name1, "ASR_year2" = name2))
    
    ##EAPC_SDI
    
    EAPC_data <- digestdata%>%
      filter( digestdata$year %in% c(year1:year2),
              digestdata$age == 'Age-standardized',  
              digestdata$metric == 'Rate', 
              digestdata$measure == measures[m],  
              digestdata$cause == aimdisease[d] ,
              digestdata$sex == "Both",
              location %in% c(main_c,SDI_c)     
      )
    
    length(unique(EAPC_data$location))  
    EAPC_year2 <- data.table() 
    for (i in unique(EAPC_data$location)) {
      temp <- EAPC_data %>% filter(location == i)
      temp$lograte <- log(temp$val)
      fit <- lm(lograte ~ year, temp)
      b <- coefficients(fit)
      CI <- confint(fit)
      EAPC0.5 <- format(round(100 * (exp(b[2]) - 1), digits = 2), nsmall = 2)
      EAPC0.975 <- format(round(100 * (exp(CI[4]) - 1), digits = 2), nsmall = 2)
      EAPC0.025 <- format(round(100 * (exp(CI[2]) - 1), digits = 2), nsmall = 2)

      location <- EAPC_data$location[EAPC_data$location == i] %>% unique()
      EAPC_t <- cbind(location, EAPC0.025, EAPC0.5, EAPC0.975)
      EAPC_year2 <- rbind(EAPC_year2,EAPC_t)
    }
    ##bar plot of EAPC by different regions 
    EAPC_SDI <- EAPC_year2 %>% filter(location %in% c("Global",SDI_c))
    #
    EAPC_SDI[,c(2,3,4)] <- as.data.frame(apply(EAPC_SDI[,c(2,3,4)],2,as.numeric))
    EAPC_mainc <- EAPC_year2 %>% filter(location %in% main_c[-1])
    EAPC_mainc[,c(2,3,4)] <- as.data.frame(apply(EAPC_mainc[,c(2,3,4)],2,as.numeric))
    p <- ggplot(data = EAPC_SDI,aes(x=factor(location,levels = c("Global",SDI_c),ordered = T),
                                    y=EAPC0.5))+
      geom_bar(position="dodge",stat="identity",fill = "#99CCFF",colour = "black")+
      geom_errorbar(aes(ymax = EAPC0.975, ymin = EAPC0.025),
                    position = position_dodge(0.9), width = 0.15)+
      theme_bw()+
      xlab(label = "")+
      ylab(label = paste0("EAPC in age-standardised ",aimdiseaselabel[d]," rate of ",measureslabel[m]))+
      theme(plot.title = element_text(hjust = 0.5))+  ###
      theme(legend.position ='right')
    p2<- p+ scale_y_continuous(n.breaks = 6)  
    p3<- p2 + ylim(min(EAPC_SDI$EAPC0.025, 0)*1.2, max(EAPC_SDI$EAPC0.5,0)*1.2)
    p4 <- p3+theme(axis.text.x = element_text(size = 12, #family = "Serif", 
                                              color = "black",
                                              face = "plain", vjust = 0.5, hjust = 0.5, angle = 30)) 
    
    p4  
    ggsave(plot = p4,filename = paste0(resultdir,"EAPC in age-standardised ",aimdiseaselabel[d]," rate of ",measureslabel[m],' SDI_in_',year1,'_to_',year2,'.png'),width = 8,height = 8)
    ggsave(plot = p4,filename = paste0(resultdir,"EAPC in age-standardised ",aimdiseaselabel[d]," rate of ",measureslabel[m],' SDI_in_',year1,'_to_',year2,'.pdf'),width = 8,height = 8)
    
    #region EAPC
    EAPC_mainc <- EAPC_mainc[order(EAPC0.5),]
    p <- ggplot(data = EAPC_mainc,aes(x=factor(location,levels =location,ordered = T),
                                      y=EAPC0.5))+
      geom_bar(position="dodge",stat="identity",fill = "#99CCFF",colour = "black")+
      geom_errorbar(aes(ymax = EAPC0.975, ymin = EAPC0.025),
                    position = position_dodge(0.9), width = 0.15)+
      theme_bw()+
      xlab(label = "")+
      ylab(label = paste0("EAPC in age-standardised ",aimdiseaselabel[d]," rate of ",measureslabel[m]))+
      theme(plot.title = element_text(hjust = 0.5))+  ###
      theme(legend.position ='right')
    p2<- p+ scale_y_continuous(n.breaks = 6)  
    p3<- p2 + ylim(min(EAPC_SDI$EAPC0.025, 0)*1.2, max(EAPC_SDI$EAPC0.5,0)*1.2)
    p4 <- p3+theme(axis.text.x = element_text(size = 12, color = "black",
                                              face = "plain", vjust = 0.5, hjust = 0.5, angle = 30)) 
    
    p4  
    ggsave(plot = p4,filename = paste0(resultdir,"EAPC in age-standardised ",aimdiseaselabel[d]," rate of ",measureslabel[m],' Mainc_in_',year1,'_to_',year2,'.png'),width = 8,height = 8)
    ggsave(plot = p4,filename = paste0(resultdir,"EAPC in age-standardised ",aimdiseaselabel[d]," rate of ",measureslabel[m],' Mainc_in_',year1,'_to_',year2,'.pdf'),width = 8,height = 8)
    
    ##orgnize the data
    EAPC_year2 <- EAPC_year2[,.(location, EAPC_year2= EAPC_year2$EAPC0.5,EAPC = paste(EAPC_year2$EAPC0.5,' (', 
                                                                                      EAPC_year2$EAPC0.025,"-", 
                                                                                      EAPC_year2$EAPC0.975, ')', sep = ''))]
    
    index_data_both <- merge(x = index_data_both,y = EAPC_year2[,c("location","EAPC")],by.x = 'location',by.y="location")
    
    name3<-paste0("EAPC_",aimdiseaselabel[d])
    index_data_both <- plyr::rename(index_data_both,c("EAPC" = name3))
  }
  #
  index_data_both$location <- factor(index_data_both$location,levels = c(main_c[1],SDI_c,main_c[-c(1)]),ordered = T)
  index_data_both <- index_data_both %>% arrange(location)
  write.csv(x = index_data_both,file = paste(resultdir,"ASR and EAPC of_",measureslabel[m],'_four_all.csv',sep = ''),row.names = F)
  index_data_both <- data.frame(location=c(SDI_c,main_c))
}
gc()

###############2、ASR INDEXES of PAH by sex (year2)################
sex_index <- unique(digestdata$sex)
sex_index
index_data <- data.frame(location=main_c)

##bar plot by sex
for (d in c(1:length(unique(digestdata$cause)))) {
  
  for (m in c(1:length(measures))) {
    dataname_year1_year2 <- digestdata%>%
      filter(digestdata$year == year2,
             digestdata$age == 'Age-standardized', 
             digestdata$metric == 'Rate', 
             digestdata$measure == measures[m],  
             digestdata$cause==aimdisease[d] ,
             digestdata$sex %in% c("Male","Female"),
             location %in% main_c 
      )
    dataname_year1_year2 <- dataname_year1_year2%>% select(cause,measure,year,sex,location,val,upper,lower)##选择列
    dataname_year1_year2$val <- round(dataname_year1_year2$val,1)
    write.csv(x = dataname_year1_year2,file = paste(resultdir,measureslabel[m],'_of_',aimdiseaselabel[d],'_of_ASR_',year1,'_to_',year2,'_mainregions_.csv',sep = ''),row.names = F)

    COLORtake <- colorRampPalette(c("#E32D32","#005792"))
    cc <- COLORtake(length(unique(dataname_year1_year2$sex))) 
    scales::show_col(cc)
    #ploting 
    p2 <- ggplot(dataname_year1_year2, aes(
      x = factor(location,levels = main_c,ordered = T),            
      y = ifelse(sex == "Female", val, -val),  # 
      fill = sex)) +
      geom_bar(stat = 'identity')+    
      geom_errorbar(aes(ymax = ifelse(sex == "Female", upper, -upper), 
                        ymin = ifelse(sex == "Female", lower, -lower)),
                    position = position_dodge(0.1), width = 0.5)+
      theme_bw()+ #
      xlab(label = "")+
      ylab(label = "")+
      labs(title = paste("ASR ",measureslabel[m]," attributed to ",aimdiseaselabel[d]," by sex (",year2,") ",sep = ''),
           fixed = T,pattern = ".txt",replacement = "")+
      coord_flip()+                                               
      scale_fill_manual(values =cc )+
      scale_y_continuous(                                         
        labels = abs,                                             
        expand = expansion(mult = c(0.1, 0.1)                   
        )) 
    
    p3 <- p2 + theme(axis.text.x = element_text(size = 12,color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
      theme(axis.text.y = element_text(size = 14, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))
    p3   
    ggsave(plot = p3,filename = paste0(resultdir,"ASR ",measureslabel[m]," attributed to ",aimdiseaselabel[d]," by sex (",year2,").png"),
           width = 18,height =11)
    ggsave(plot = p3,filename = paste0(resultdir,"ASR ",measureslabel[m]," attributed to ",aimdiseaselabel[d]," by sex (",year2,").pdf"),
           width = 18,height =11)
  }
}
gc()

##########3、SDI ########################
#####SDI data from  https://ghdx.healthdata.org/record/global-burden-disease-study-2021-gbd-2021-socio-demographic-index-sdi-1950%E2%80%932021
library(dplyr)
library(data.table)
library(tidyverse)
library(gtools)
SDIdata <- fread(input = "F:/workforSCI/252/dataanalyse/GBD/IHME_GBD_SDI_2021_SDI_1950_2021_Y2024M05D16.csv")
colnames(SDIdata)
unique(SDIdata$sex)   # both only
unique(SDIdata$year)   # 1950-2021
unique(SDIdata$age_group_name) # all Ages only

sum(unique(SDIdata$location_name) %in% main_c )

main_SDI_year2 <- SDIdata %>% filter(location_name %in% main_c,
                                     year_id %in% c(year1:year2)
) %>% select(c(3,4,9:11)) %>% # location,year,and val
  rename("year" ="year_id","location"="location_name")
######remove duplicate rows
main_SDI_year2 <- main_SDI_year2[!duplicated(main_SDI_year2),]

for (m in c(1:length(measures))) {
  disease_mainc <- digestdata%>%
    filter(digestdata$year %in% c(year1:year2),
           digestdata$age == 'Age-standardized', 
           digestdata$metric == 'Rate',
           digestdata$measure == measures[m],  
           digestdata$sex == "Both",
           location %in% main_c    
    )%>% select(c("location","year","val","upper","lower"))
  
  sdi_asr <- merge(x = disease_mainc,y = main_SDI_year2,by = c("location","year"))
  sdi_asr$location <- factor(x =sdi_asr$location,
                             levels = main_c,ordered = T
  )
  #choose colors
  COLORtake <- colorRampPalette(c("#005792","#CCFFFF", "#E32D32"))
  cc <- COLORtake(length(main_c)) 
  scales::show_col(cc)
  p <- ggplot(data = sdi_asr,aes(x=mean_value, y=val))+
    geom_point ( aes( group=location,
                      color=location,shape = location), 
                 size=2 ) +
    scale_shape_manual(values = c(1:22))+  ###reflect shapes
    geom_smooth(colour='black',stat ="smooth",method='loess',se=T ,span=0.5)+
    scale_fill_manual(values =cc )+ ##for categorical variables
    xlab(label = "SDI")+
    ylab(label = paste0("ASR of ",measureslabel[m]," (/10^5)"))+
    labs(title = paste("SDI_",measureslabel[m],"_ASR_of_",aimdiseaselabel[d],'_',year1,'_to_',year2,sep = ''),
         fixed = T,pattern = ".txt",replacement = "")+
    theme(plot.title = element_text(hjust = 0.5))+  
    guides(fill = guide_legend(title='main_regions'))+ 
    theme_bw(base_rect_size = 0.1,base_line_size = 0.1) +
    theme(legend.position ='right')
  p + theme(axis.text.x = element_text(size = 12, color = "black",
                                       face = "plain", vjust = 0.5, hjust = 0.5, angle = 30))+
    theme(axis.text.y = element_text(size = 14, color = "black",
                                     face = "plain", vjust = 0.5, hjust = 0.5, angle = 30))
  p
  
  ggsave(plot = p,filename = paste(resultdir,"SDI_",measureslabel[m],"_ASR_of_",aimdiseaselabel[d],'_',year1,'_to_',year2,'.png',sep = ''),
         width = 15,height =8)
  ggsave(plot = p,filename = paste(resultdir,"SDI_",measureslabel[m],"_ASR_of_",aimdiseaselabel[d],'_',year1,'_to_',year2,'.pdf',sep = ''),
         width = 15,height =8)
  
}

######3.2 countries' SDI############
if (sum(unique(SDIdata$location_name) %in% location_c )<204) {
  SDIdata[location_name=="Türkiye",3] <- "Turkey"
}  ###“Türkiye”==Turkey 
location_c[which(!location_c %in% unique(SDIdata$location_name))]
sum(unique(SDIdata$location_name) %in% location_c)
#####year1 and year2
for (y in c(year1,year2)) {
  
  SDI_count <- SDIdata %>% filter(location_name %in% location_c,year_id == y
  ) %>% select(c(3,4,9:11)) %>% 
    rename("year" ="year_id","location"="location_name")
  ######remove the duplicate
  SDI_count <- SDI_count[!duplicated(SDI_count$location),]
  
  library(ggrepel)
  
  for (m in c(1:length(measures))) {

    disease_count <- digestdata%>%
      filter(digestdata$year == y,
             digestdata$age == 'Age-standardized',  
             digestdata$metric == 'Rate', 
             digestdata$measure == measures[m], 
             digestdata$sex == "Both",
             location %in% location_c   
      )%>% select(c("location","year","val","upper","lower"))
    
    sdi_asr <- merge(x = disease_count,y = SDI_count,by = c("location","year"))
    sdi_asr$location <- factor(x =sdi_asr$location,
                               levels = location_c,ordered = T
    )
    
    p <- ggplot(data = sdi_asr,aes(x=mean_value, y=val))+
      geom_point ( aes( group=location,
                        #color=location,
                        size = val),
                   color="#99CCFF"
      ) +
      geom_smooth(colour='#CC3333',stat ="smooth",method='loess',se=T ,span=0.5,size= 0.8)+
      theme_void()+ 
      xlab(label = "SDI")+
      ylab(label = paste0("ASR of ",measureslabel[m]," (/10^6)"))+
      labs(title = paste("SDI_",measureslabel[m],"_ASR_of_",aimdiseaselabel[d],'_',y,sep = ''),
           fixed = T,pattern = ".txt",replacement = "")+
      theme(plot.title = element_text(hjust = 0.5))+  ###
      guides(fill = guide_legend(title='countries and regions'))+ 
      theme_bw(base_rect_size = 0.1,base_line_size = 0.1) +
      theme(legend.position ='right')
    p + theme(axis.text.x = element_text(size = 12,  color = "black",
                                         face = "plain", vjust = 0.5, hjust = 0.5, angle = 30))+
      theme(axis.text.y = element_text(size = 14,  color = "black",
                                       face = "plain", vjust = 0.5, hjust = 0.5, angle = 30))
    p
    ggsave(plot = p,filename = paste(resultdir,"SDI_",measureslabel[m],"_ASR_of_",aimdiseaselabel[d],'_',y,'.png',sep = ''),
           width = 15,height =8)
    ggsave(plot = p,filename = paste(resultdir,"SDI_",measureslabel[m],"_ASR_of_",aimdiseaselabel[d],'_',y,'.pdf',sep = ''),
           width = 15,height =8)
    
  }
}

####3.3 90-21 SDI#############
library(ggrepel)
library(ggplot2)

if (sum(unique(SDIdata$location_name) %in% location_c )<204) {
  SDIdata[location=="Türkiye",3] <- "Turkey"
}  ###“Türkiye”==Turkey 
which(!location_c %in% unique(SDIdata$location_name))
sum(unique(SDIdata$location_name) %in% location_c )
#####year1
SDI_count_year1 <- SDIdata %>% filter(location_name %in% location_c,
                                      year_id == year1 ) %>% select(c(3,4,9:11))%>%
  rename("year" ="year_id","location"="location_name")
######remove the duplicate
SDI_count_year1 <- SDI_count_year1[!duplicated(SDI_count_year1$location),]
#####2021年
SDI_count_year2 <- SDIdata %>% filter(location_name %in% location_c,
                                      year_id == year2) %>% select(c(3,4,9:11)) %>%
  rename("year" ="year_id","location"="location_name")
######remove the duplicate
SDI_count_year2 <- SDI_count_year2[!duplicated(SDI_count_year2$location),]


COLORtake <- colorRampPalette(c("#99CCFF","#CC3333"))

cc <- COLORtake(c(1,2)) 
scales::show_col(cc)

for (m in c(1:length(measures))) {
  disease_countyear1 <- digestdata%>%
    filter(digestdata$year == year1,
           digestdata$age == 'Age-standardized',  
           digestdata$metric == 'Rate', 
           digestdata$measure == measures[m], 
           digestdata$sex == "Both",
           location %in% location_c    
    )%>% select(c("location","year","val","upper","lower"))
  
  disease_countyear2 <- digestdata%>%
    filter(digestdata$year == year2,
           digestdata$age == 'Age-standardized',  
           digestdata$metric == 'Rate', 
           digestdata$measure == measures[m],  
           digestdata$sex == "Both",
           location %in% location_c    
    )%>% select(c("location","year","val","upper","lower"))
  
  sdi_asryear1 <- merge(x = disease_countyear1,y = SDI_count_year1,by = c("location","year"))
  sdi_asryear2 <- merge(x = disease_countyear2,y = SDI_count_year2,by = c("location","year"))
  sdi_asr <- rbind(sdi_asryear1,sdi_asryear2)
  
  sdi_asr$location <- factor(x =sdi_asr$location,
                             levels = location_c,ordered = T)
  sdi_asr$year <- factor(x =sdi_asr$year,
                         levels = c(year1,year2),ordered = T)
  
  p <- ggplot(data = sdi_asr,aes(x=mean_value, y=val))+
    geom_point ( aes( group=year,
                      color= year,
                      size = val)
    ) +
    geom_smooth(aes(colour=year),stat ="smooth",method='loess',se=T ,span=0.5,size= 0.5)+
    scale_fill_manual(cc )+ #
    theme_void() + ##
    xlab(label = "SDI")+
    ylab(label = paste0("ASR of ",measureslabel[m]," (/10^6)"))+
    labs(title = paste("SDI_and_ASR_of_",measureslabel[m],"_of_",aimdiseaselabel[d],'_',year1,' and ',year2,sep = ''),
         fixed = T,pattern = ".txt",replacement = "")+
    theme(plot.title = element_text(hjust = 0.5))+  #
    guides(fill = guide_legend(title='countries and regions'))+ 
    theme_bw(base_rect_size = 0.1,base_line_size = 0.1) +
    theme(legend.position ='right')

  p + theme(axis.text.x = element_text(size = 12,  color = "black",
                                       face = "plain", vjust = 0.5, hjust = 0.5, angle = 30))+
    theme(axis.text.y = element_text(size = 14, color = "black",
                                     face = "plain", vjust = 0.5, hjust = 0.5, angle = 30))
  p
  
  ggsave(plot = p,filename = paste(resultdir,"SDI_and_ASR_of_",measureslabel[m],"_of_",aimdiseaselabel[d],'_',year1,' and ',year2,'.png',sep = ''),
         width = 12,height =9)
  ggsave(plot = p,filename = paste(resultdir,"SDI_and_ASR_of_",measureslabel[m],"_of_",aimdiseaselabel[d],'_',year1,' and ',year2,'.pdf',sep = ''),
         width = 12,height =9)
  
}


#######4、ASR and time ##########################
sex_index <- unique(digestdata$sex)
sex_index
aimdisease <- unique(digestdata$cause)
aimdisease
aimdiseaselabel

COLORtake <- colorRampPalette(c("#005792","#CC9966","#CCFFFF","#009966","#000000", "#E32D32"))
cc <- COLORtake(length(main_c)) 
scales::show_col(cc)
for (d in c(1:length(unique(digestdata$cause)))) {
  for (s in sex_index) {
    for (m in c(1:length(measures))) {
      dataname_year1_year2 <- digestdata%>%
        filter(digestdata$year %in%  aimyears,
               digestdata$age == 'Age-standardized',  #
               digestdata$metric == 'Rate', #
               digestdata$measure == measures[m],  #
               digestdata$cause==aimdisease[d] ,
               digestdata$sex == s,
               location %in% main_c #
        ) %>% #
        select(c("location","year","val","upper","lower"))
      
      dataname_year1_year2$val <- round(dataname_year1_year2$val,2)
      write.csv(x = dataname_year1_year2,file = paste(resultdir,measureslabel[m],'_of_',aimdiseaselabel[d],'_of_',s,'ASR_',year1,'_to_',year2,'_mainregions_.csv',sep = ''),row.names = F)
      
      dataname_year1_year2$location <- factor(dataname_year1_year2$location,levels = main_c,ordered = T)
      p <- ggplot(data = dataname_year1_year2)+
        geom_point ( aes(x=year, y=val, group=location,
                         color=location,shape = location), 
                     size=2 ) +
        scale_shape_manual(values = c(1:22))+  ##
        geom_line ( aes(x=year, y=val, group=location,color=location), ) +
        scale_fill_manual(values =cc )+ #
        theme_classic()+
        xlab(label = "")+
        ylab(label = "ASR(/10^6)")+
        labs(title = paste(measureslabel[m],'_of_',aimdiseaselabel[d],'_of_',s,'_ASR_',year1,'_to_',year2,'_mainregions',sep = ''),
             fixed = T,pattern = ".txt",replacement = "")+
        theme(plot.title = element_text(hjust = 0.5))+  #
        guides(fill = guide_legend(title='main_regions'))+ 
        theme_classic() +  # 
        theme(
          plot.title = element_text(hjust = 0.5),  # 
          legend.position = 'right',
          axis.text.x = element_text(size = 12, color = "black", angle = 30, hjust = 0.5),  # 
          axis.text.y = element_text(size = 14, color = "black")  #
        ) +
        scale_x_continuous(breaks = seq(year1, year2, 2)) +  # 
        scale_y_continuous(breaks = scales::breaks_extended(n = 9))  # 
      
      p 
      ggsave(plot = p,filename = paste(resultdir,measureslabel[m],'_of_',aimdiseaselabel[d],'_of_',s,'_ASR_',year1,'_to_',year2,'_mainregions.png',sep = ''),
             width = 15,height =11)
      ggsave(plot = p,filename = paste(resultdir,measureslabel[m],'_of_',aimdiseaselabel[d],'_of_',s,'_ASR_',year1,'_to_',year2,'_mainregions.pdf',sep = ''),
             width = 15,height =11)
    }
  }
}  
gc()

######part 3#############
library(dplyr)
library(data.table)
library(tidyverse)
library(car)# 
library(mgcv)# 
library(splines)# 
library(broom)
library(gtools)
library(ggpmisc)  #
library(healthequal)
resultdir <- "F:/workforSCI/252/dataanalyse/GBD/results/NAFLD/over45/test/inequality_analysis/"
dir.create(resultdir)

aimyears <- c(year1:year2)
unitage <- c(#"0-4" ,"5-9","10-14", "15-19", "20-24" ,"25-29", "30-34" ,"35-39" ,"40-44",
  "45-49", "50-54", "55-59", "60-64", "65-69","70-74", "75-79" ,"80-84", "85+")
agegroups <- mixedsort(unique(digestdata$age))[c(1:11)]
agegroups
sex_index <- c("Male","Female","Both")

#read SDI data ##########
SDIdata <- fread(input = "F:/workforSCI/252/dataanalyse/GBD/IHME_GBD_SDI_2021_SDI_1950_2021_Y2024M05D16.csv")
colnames(SDIdata)
unique(SDIdata$sex)   # both only
unique(SDIdata$year_id)   # 1950-2021
unique(SDIdata$age_group_name) # All Ages only
##select 204regions
if (sum(unique(SDIdata$location_name) %in% location_c )<204) {
  SDIdata[location_name=="Türkiye",3] <- "Turkey"
}  ###“Türkiye”==Turkey 
location_c[which(!location_c %in% unique(SDIdata$location_name))]
sum(unique(SDIdata$location_name) %in% location_c)

SDI_count <- SDIdata %>% filter(location_name %in% location_c,
                                year_id %in% c(year1,year2)
) %>% select(c("location_name","year_id", "mean_value")) %>% 
  rename(c(sdi=mean_value))
######remove the duplicated
SDI_count <- SDI_count[!duplicated(SDI_count[,c(1,2)]),]

for (m in c(1:length(unique(digestdata$measure)))) {
  for (d in c(1:length(unique(digestdata$cause)))) {
    colnames(digestdata)
    inequality_90_21 <- digestdata %>% 
      filter(year  %in% c(year1,year2),
             age == "All ages",  
             measure == measures[m],  
             cause == aimdisease[d],
             sex == "Both",
             location %in% location_c ) %>%
      select(c(location,metric,year,val))%>% 
      rename(c(year_id=year,location_name = location))
    
    
    # combine data
    data_combined <- left_join(x = inequality_90_21,y = SDI_count,by= c("location_name","year_id"))
    
    ##population data：https://population.un.org/wpp/Download/Standard/Population/
    ##2、reading population data############################
    populationdata <- fread(input = "F:/workforSCI/252/dataanalyse/GBD/population/WPP2024_Population1JanuaryByAge5GroupSex_Medium.csv.gz")
    colnames(populationdata)
    unique(populationdata$AgeGrp)
    locationmatch <- data.frame(GBDlocation=c("Taiwan (Province of China)","Palestine","Turkey","Micronesia (Federated States of)","Democratic People's Republic of Korea"),
                                poplocation=c("China, Taiwan Province of China","State of Palestine","Türkiye","Micronesia","Dem. People's Republic of Korea"))
    martch_location <-location_c
    if (sum(unique(populationdata$Location) %in% location_c)<204) {
      for (l in location_c[!location_c %in% unique(populationdata$Location)]) {
        populationdata[Location==locationmatch$poplocation[locationmatch$GBDlocation==l],10] <- l
      }
    }
    sum(unique(populationdata$Location) %in% martch_location)
    
   pop_90_21 <- populationdata %>% filter(Time%in% c(year1,year2),Location %in% martch_location)%>% ##指定预测的区域
      select(c(10,13,15,20))  ###*1000,PopMale,PopFemale,PopTotal
    pop_90_21[,4] <- pop_90_21[,4]*1000
    
    ##
    pop_90_21 <- aggregate(x = pop_90_21[,4],by=list(type=pop_90_21$Location,pop_90_21$Time),sum)
    colnames(pop_90_21) <- c("location_name","year_id","population ")
    #Convert the population into the unit of millions.
    pop_90_21$population_M <- round(pop_90_21$population/1000000,2) ##
    #Modify the regions name to align with the regions name in the GBD incidence data
    
    if (sum(unique(SDIdata$location_name) %in% location_c )<204) {
      SDIdata[location_name=="Türkiye",3] <- "Turkey"
    }  ###“Türkiye”==Turkey 地区
    location_c[which(!location_c %in% unique(pop_90_21$location))]
    sum(unique(pop_90_21$location) %in% location_c)

    data_combined <- left_join(x = data_combined,y = pop_90_21,by = c("location_name","year_id"))
    
    
    ## 1.plot preparation-----------------------------------------------------------------
    # calculate the total population
    a <- data_combined %>%
      filter(metric=="Number") %>%
      group_by(year_id) %>%
      summarise(sum=sum(population_M))
    pop2000 <- a$sum[1]
    pop2021 <- a$sum[2]
    # Calculate the weighted order
    rank <- data_combined %>%
      mutate(pop_global=ifelse(year_id==year1,pop2000,pop2021)) %>%
      group_by(year_id,metric) %>%
      arrange(sdi) %>%
      mutate(cummu=cumsum(population_M)) %>% # accumulated population
      mutate(half=population_M/2) %>% # 
      mutate(midpoint=cummu-half) %>% # midpoint of population number
      mutate(weighted_order=midpoint/pop_global) # 
    rank$year_id <- factor(rank$year_id)
    # 
    temp1 <- rank %>%
      filter(metric=="Rate") %>%
      filter(year_id==year1)
    temp2 <- rank %>%
      filter(metric=="Rate") %>%
      filter(year_id==year2)
    
    # 
    fit1 <- lm(data = temp1,val~weighted_order)
    fit2 <- lm(data = temp2,val~weighted_order)
    coef(fit1)
    coef(fit2)
    
    # Verify whether heteroscedasticity exists
    ncvTest(fit1)
    ncvTest(fit2)
    library(MASS)# 
    r.huber1 <- rlm(data = temp1,val~weighted_order)
    r.huber2 <- rlm(data = temp2,val~weighted_order)
    # 
    coef(r.huber1)
    coef(r.huber2)
    # Calculate the 95% confidence interval for robust regression
    confint.default(r.huber1)  
    confint.default(r.huber2)  
    detach("package:MASS", unload = TRUE)
    
    library(ggpubr)
    # 2.ploting  ----------------------------------------------------------------------
    color <- c("#6699FF","#990000")
    
    max(rank[rank$metric=="Rate",4])
    colnames(rank)
    p1 <- rank %>%
      filter(metric=="Rate") %>%
      ggplot(aes(x=weighted_order,y=val,fill=year_id,group=year_id,color=year_id))+
      geom_point(aes(color=year_id,size=population_M),alpha=0.8,shape=21)+
      scale_size_area("Population\n(million)",breaks=c(200,400,600,800,1000))+
      geom_smooth(method = "lm",size=0.6,alpha=0.1)+
      scale_fill_manual(values = color)+
      scale_color_manual(values = color)+
      annotate("text",label="Slope Index of Inequality",x=0.2,y=max(rank[rank$metric=="Rate",4])*0.9,size=8,angle=0)+
      annotate("text",label=round(coef(r.huber1)[[2]],2),x=0.15,y=max(rank[rank$metric=="Rate",4])*0.8,size=7,color="#6699FF")+ # coef(r.huber1) 权重的系数即斜率指数
      annotate("text",label=round(coef(r.huber2)[[2]],2),x=0.15,y=max(rank[rank$metric=="Rate",4])*0.75,size=7,color="#990000")+  ##coef(r.huber2)  权重的系数即斜率指数
      scale_x_continuous(limits = c(0,1),labels = c("0","0.25","0.50","0.75","1.00"))+
      xlab("Relative rank by SDI")+
      ylab(paste0("Crude ",measureslabel[m]," rate (per 100,000)"))+
      theme_bw()+
      theme(axis.text.x = element_text(size = 12,  color = "black",face = "plain", vjust = 0.5, hjust = 0.5))+
      theme(axis.text.y = element_text(size = 14, color = "black",face = "plain", vjust = 0.5, hjust = 0.5))+
      theme(axis.title.x  = element_text(size = 16,  color = "black",face = "plain", vjust = 0.5, hjust = 0.5))+
      theme(axis.title.y = element_text(size = 16, color = "black",face = "plain", vjust = 0.5, hjust = 0.5))
    
    p1
    ggsave(plot = p1,filename = paste0(resultdir," Health inequality regression curves of ",measureslabel[m], "for ",aimdiseaselabel[d],".png"),width = 12,height = 8)
    ggsave(plot = p1,filename = paste0(resultdir," Health inequality regression curves of ",measureslabel[m], "for ",aimdiseaselabel[d],".pdf"),width = 12,height = 8)
    
    # CI illustration ----------------------------------------------------------------
    
    # 1.data preparation ------------------------------------------------------------------
    a <- data_combined %>%
      filter(metric=="Number") %>%
      group_by(year_id) %>%
      summarise(sum=sum(val))
    daly2000 <- a$sum[1]
    daly2021 <- a$sum[2]
    
    ci <- rank %>%
      filter(metric=="Number") %>%
      mutate(total_daly=ifelse(year_id==year1,daly2000,daly2021)) %>%
      group_by(year_id) %>%
      arrange(sdi) %>%
      mutate(cummu_daly=cumsum(val)) %>% # 
      mutate(frac_daly=cummu_daly/total_daly) %>% # 
      mutate(frac_population=cummu/pop_global) # 
    #####calculate the ci
    # 
    temp3 <- ci %>%
      filter(metric=="Number") %>%
      filter(year_id == year1)
    temp4 <- ci %>%
      filter(metric=="Number") %>%
      filter(year_id == year2)
    ##
    CI_1990 <- 2 * (sum(temp3$frac_daly) / nrow(temp3)) - 1
    CI_1990
    CI_2021 <- 2 * (sum(temp4$frac_daly) / nrow(temp4)) - 1
    CI_2021
    # 2.ploting  --------------------------------------------------------------------
    p2 <- ci %>% 
      ggplot(aes(x=frac_population,y=frac_daly,fill=year_id,color=year_id,group=year_id))+
      geom_segment(x=0,xend=1,
                   y=0,yend=0,
                   linetype=1,size=1,color="gray")+
      geom_segment(x=1,xend=1,
                   y=0,yend=1,
                   linetype=1,size=1,color="gray")+
      # Diagonal
      geom_segment(x=0,xend=1,
                   y=0,yend=1,
                   color="#CD853F",linetype=1,size=0.7,alpha=1)+
      # dot
      geom_point(aes(fill=year_id,size=population_M),alpha=0.75,shape=21)+
      scale_fill_manual(values = color)+
      scale_size_area("Population\n(million)",breaks=c(200,400,600,800,1000))+
      geom_smooth(method = "gam", # 
                  formula = y ~ ns(x,
                                   knots = c(0.0000000001,0.25,0.5,0.75,0.9999999),# 设置节点为
                                   Boundary.knots = c(0,1)),
                  linetype=1,size=0.1,alpha=0.6,se=T)+
      scale_color_manual(values = color)+
      annotate("text",label="Concentration Index",x=0.2,y=0.8,size=8)+
      annotate("text",label=paste0(year1,"  ",round(CI_1990,2)),x=0.2,y=0.75,size=7,color="#6699FF")+
      annotate("text",label=paste0(year2,"  ",round(CI_2021,2)),x=0.2,y=0.7,size=7,color="#990000")+

      geom_text(aes(label=ifelse(location_name %in% a&year_id==year1,
                                 as.character(location_name),"")),
                hjust=-0.6,vjust=0.8,
                size=3)+
      geom_text(aes(label=ifelse(location_name%in%a&year_id==year2,
                                 as.character(location_name),"")),
                hjust=1.8,vjust=-0.0,
                size=3)+
      # xy 标签
      xlab("Cumulative fraction of population ranked by SDI")+
      ylab(paste0("Cumulative fraction of ",measureslabel[m]))+
      theme_bw()+
      theme(axis.text.x = element_text(size = 12,  color = "black",face = "plain", vjust = 0.5, hjust = 0.5))+
      theme(axis.text.y = element_text(size = 14, color = "black",face = "plain", vjust = 0.5, hjust = 0.5))+
      theme(axis.title.x  = element_text(size = 16,  color = "black",face = "plain", vjust = 0.5, hjust = 0.5))+
      theme(axis.title.y = element_text(size = 16, color = "black",face = "plain", vjust = 0.5, hjust = 0.5))
    p2
    
    ggsave(plot = p2,filename = paste0(resultdir,"CI of ",measureslabel[m], "for ",aimdiseaselabel[d],".png"),width = 15,height = 8)
    ggsave(plot = p2,filename = paste0(resultdir,"CI of ",measureslabel[m], "for ",aimdiseaselabel[d],".pdf"),width = 15,height = 8)
    
  }
}
gc()
##part 4 ###############
library(INLA)
library(BAPC)
library(reshape2)
library(epitools)
library(ggplot2)
library(dplyr)
library(gtools)   ###mixedsort排序
resultdir <- "F:/workforSCI/252/dataanalyse/GBD/results/NAFLD/over45/test/forecaste12/"
dir.create(resultdir)
aimyears <- c(year1:year2)
unitage <- c("45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85-89","90-94","95+")
gegroups <- mixedsort(unique(digestdata$age))[c(1:11)]
agegroups
##Time-segmented Function
function_year5_2021 <- function(table, start_year, end_year, current_year){
  
  remain <- current_year - round((current_year - start_year)/5) * 5 
  years <- NULL
  for (i in start_year:end_year) {
    if((i - current_year)/5 - round((i - current_year)/5) == 0){
      if(i == remain){
        temp <- paste(start_year, i, sep = '-')
        years <- append(years, temp)
      }
      else{
        temp <- paste(i-4, i, sep = '-')
        years <- append(years, temp)
      }
    }
  }
  
  table <- table %>% as.data.frame()
  new_years <- seq(start_year,end_year,1)
  new_table <- matrix(data = rep(0, length(years)*nrow(table)), ncol = length(years), nrow = nrow(table)) %>% as.data.frame()
  colnames(new_table) <- years
  
  j = 1
  for (i in 1:(end_year - start_year + 1)) {
    if((new_years[i] - 2021)/5 - round((new_years[i] - 2021)/5) != 0){
      new_table[, years[j]] <- new_table[,years[j]] + table[,as.character(new_years[i])]
    }
    else{
      if(j == 1){
        new_table[,years[j]] <- (new_table[,years[j]] + table[,as.character(new_years[i])]) / (remain - start_year + 1)
      }
      else{
        new_table[,years[j]] <- (new_table[,years[j]] + table[,as.character(new_years[i])]) / 5
      }
      j = j + 1
    }
  }
  return(new_table)
} ##注意修改时间和行数
function_year5_2050 <- function(table, start_year, end_year, current_year){
  years <- NULL
  year_seq <- seq(start_year,end_year,5)
  for (i in year_seq) {
    if (i < end_year-5) {
      
      temp <- paste(i,i+4, sep = '-')
      years <- append(years, temp)
    }else{
      temp <- paste(i, end_year, sep = '-')
      years <- append(years, temp)
    } 
  }
  
  
  table <- table %>% as.data.frame()
  new_years <- seq(start_year,end_year,1)
  new_table <- matrix(data = rep(0, length(years)*nrow(table)), ncol = length(years), nrow = nrow(table)) %>% as.data.frame()
  colnames(new_table) <- years
  rownames(new_table) <- rownames(table)
  
  for (i in c(1:length(year_seq))) {
    if (year_seq[i] < end_year-5) {
      new_table[,i]=rowMeans(table[,as.character(c(year_seq[i]:(year_seq[i]+4)))])
      
    }else{
      new_table[,i]=rowMeans(table[,as.character(c(year_seq[i]:end_year))])
    } 
  }
  return(new_table)
} ##注意 目前只有start endyears 参数，end year = current year

####2、data preparation########################################
populationdata <- fread(input = "F:/workforSCI/252/dataanalyse/GBD/population/WPP2024_Population1JanuaryByAge5GroupSex_Medium.csv.gz")
colnames(populationdata)
unique(populationdata$AgeGrp)
###specify the period 
for (g in c(18:20)) { ##18:20 for c("PopMale" ，"PopFemale"，"PopTotal")
  region_population <- populationdata%>% filter(Time%in% c(year1:year3),Location == 'World')%>% ##指定预测的区域
    select(c(13,15,g))  ##
  region_population[,3] <- region_population[,3]*1000
  region_population <- dcast(data = region_population, AgeGrp ~ Time,value.var = colnames(region_population)[3])
  #
  region_population <- as.data.frame(region_population)
  rownames(region_population) <- region_population$AgeGrp
  region_population <- region_population %>% select(-c("AgeGrp"))
  #
  region_population <- region_population[mixedsort(rownames(region_population)),]
  #######adjust the age segements
  region_population['95+', ] <- colSums(region_population[20:21, ])
  region_population <- region_population[-c(1:9,20:21),]
  
  region_population <- region_population[mixedsort(rownames(region_population)),]
  
  ######$assign#################
  ##get()
  assign(x = paste0(colnames(populationdata)[g],"_BAPC"),value = region_population[mixedsort(rownames(region_population)),])
  population_BAPC <- region_population[mixedsort(rownames(region_population)),]
  year1col <- which(as.numeric(colnames(population_BAPC))==year1)
  year2col <- which(as.numeric(colnames(population_BAPC))==year2)
  year3col <- which(as.numeric(colnames(population_BAPC))==year3)
  region_population_22_35 <- population_BAPC[,c(year2col:year3col)] ## 
  region_population_20_21 <- population_BAPC[,c(year1col:year2col)]  ##
  #
  population_nordpred21 <- function_year5_2021(table = population_BAPC,start_year = year1,end_year = year2,current_year = year2)
  population_nordpred35 <- function_year5_2050(table = region_population_22_35,start_year = year2+1,end_year = year3,current_year = year3)
  assign(x = paste0(colnames(populationdata)[g],"_nordpred"),value = cbind(population_nordpred21,population_nordpred35))
  
}
######Read in the weights of population age groups from WHO
age_stand <- fread(input = 'F:/workforSCI/252/dataanalyse/GBD/population/WHO_age_stand.csv')

wstand45 <- age_stand %>%
  mutate(age = str_remove(`Age Group`, " years")) %>%  # 
  filter(age %in% unitage) %>%
  rename(weight = `WHO world standard (&)`) %>%
  mutate(weight = weight / sum(weight))  # 

wstand45 <- wstand45$weight
sum(wstand45) #sum of wstand must be 1

######3、predictive analysis##################
aimregion = "Global"

for (d in c(1:length(unique(digestdata$cause)))) {
  for (m in c(1:length(measures))) {
    for (s in c(1:3)) {
      ages <- agegroups
      disease_BAPC <-digestdata%>% filter(digestdata$year %in%  aimyears,
                                          digestdata$age %in% agegroups,  #
                                          digestdata$metric == 'Number', #
                                          digestdata$measure == measures[m] ,  #
                                          digestdata$cause == aimdisease[d] ,
                                          digestdata$sex == sex_index[s],
                                          location == aimregion #
      ) %>% select(c(age, year, val))
      
      disease_BAPC <- dcast(data = disease_BAPC, age ~ year,value.var = 'val')

      disease_BAPC <- as.data.frame(disease_BAPC)
      rownames(disease_BAPC) <- disease_BAPC$age
      disease_BAPC <-disease_BAPC %>% select(-c(age))
      disease_BAPC <- disease_BAPC[mixedsort(rownames(disease_BAPC)),]

      rownames(disease_BAPC) <- unitage

      ################################4、BAPC plot################################
      #Population data must be integers (to satisfy the Poisson distribution).
      #The incidence rates for periods after the current time shall be replaced with NA.
      popbapcindex <- c("PopMale_BAPC","PopFemale_BAPC","PopTotal_BAPC")
      population_BAPC <- as.data.frame(t(get(popbapcindex[s])))
      population_BAPC <- round(population_BAPC,0)
      #
      rownames(disease_BAPC) <- unitage #unify row names
      disease_BAPC <- as.data.frame(t(disease_BAPC))

      disease_BAPC <- as.data.frame(apply(disease_BAPC, 2, FUN = as.integer))
      #
      disease_pro <- matrix(data = NA, nrow = c(year3-year2), ncol = length(unitage)) %>% as.data.frame() 
      rownames(disease_pro) <- seq(year2+1,2035,1)
      colnames(disease_pro) <-  colnames(disease_BAPC)
      disease_BAPC <- rbind(disease_BAPC, disease_pro)
      #
      colnames(disease_BAPC) <- colnames(population_BAPC)
      rownames(disease_BAPC) <- rownames(population_BAPC)
      ##BAPC ploting
      lc_disease <- APCList(disease_BAPC, population_BAPC, gf = 5,npred=11) #gf for age segments
      require(INLA)
      #
      bapc_result <- BAPC(lc_disease, predict = list(npredict = 5, retro = T),
                          secondDiff = F, stdweight = wstand45, verbose = T,
                          control.predictor()
      )
      pdf(file = paste0(resultdir,aimregion," ",sex_index[s]," ",measureslabel[m]," BAPC forecast.pdf"),width = 10,height = 8)
      plotBAPC(bapc_result, scale=10^5, type = 'ageStdRate', showdata = TRUE)
      dev.off()
      dev.new()
      pred_rates_matrix <- as.data.frame(bapc_result@agestd.rate)
      fwrite(x =pred_rates_matrix,file = paste0(resultdir,aimregion," ",sex_index[s]," ",measureslabel[m]," BAPC forecast.csv") ,row.names = T)
    }
  }
}
