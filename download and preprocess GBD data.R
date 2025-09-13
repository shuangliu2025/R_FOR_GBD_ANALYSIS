##R code for downloading data from GBD database and preprocess it.

library(gtools)
library(tidyverse)
####
desfile <- ("F:/pancreatitis/over45/") ###create a folder for storing data
disease_aim <- 'pancreatitis45'
dir.create(desfile)
urlpart1<-'https://dl.healthdata.org:443/gbd-api-2021-public/3d52787dd31eaf58209ec68084d0ddb5_files/' #excluding the name of compressed package 
for (i in c(1:200) ) {  ###usually less than 200
  filename1 <- paste('IHME-GBD_2021_DATA-3d52787d-',i,'.zip',sep = '') ####the name of the first compressed package
  download.file(url = paste(urlpart1,filename1,sep = ''),
                method = "wininet" , ##"internal", "libcurl", "wget", "curl" and "wininet" (Windows only)
                destfile = paste(desfile,filename1,sep = '')
  )
}
#######################decompress and combine GBD data####################
library(R.utils)##decompressing
library(gtools)### mixedsort
filenum <- length(list.files(path = desfile,pattern = ".zip"))
cat(paste("压缩包：",filenum,"个"))
for (o in c(1:filenum)) {
  filename12 <- mixedsort(list.files(path = desfile,pattern = ".zip"),decreasing = T)[o]
  filezip <- paste(desfile,filename12,sep = '')
  unzip(zipfile  = filezip,exdir = desfile, overwrite = F )
}

##combining
library(tidyverse)
library(data.table)
length(list.files(path = desfile,pattern = "*.csv"))
list.files(path = desfile,
           pattern = "*.csv", full.names=TRUE) %>% 
  lapply(fread) %>% bind_rows() %>% write_csv(.,file=paste(desfile,'/',disease_aim,".csv",sep = ''),quote="none")

##part 2 recalculate the  ASR of special ages
library(dplyr)
library(data.table)
library(tidyverse)
library(gtools)

resultdir <- "F:/pancreatitis/over45/"
dir.create(resultdir)
diseasedata <- "F:/F:/pancreatitis/over45/NAFLD.csv"
aimfile <- "NAFLD45.csv"
diseasename <- "Multidrug-resistant tuberculosis without extensive drug resistance"

target_ages <- c("45-49", "50-54", "55-59",
  "60-64",
  "65-69", "70-74","75-79", "80-84", "85-89", "90-94", "95+")

##
pancreatitis<- fread(input =diseasedata,sep = "auto",header="auto" )
if (length(colnames(pancreatitis))>15) {
  pancreatitis <- pancreatitis %>% select(c(2,4,6,8,10,12,13:16))
  colnames(pancreatitis) <- c("measure","location","sex","age","cause","metric","year","val","upper","lower")
}
unique(pancreatitis$cause)
pancreatitis <- pancreatitis %>% filter(cause==diseasename)

##

pancreatitis <- pancreatitis %>%  mutate(age = str_remove(age, " years"))%>%
                             filter(age %in% target_ages)
colnames(pancreatitis)
unique(pancreatitis$year)
unique(pancreatitis$age_name)
location_a <- as.vector(unique(pancreatitis$location))
income_c <- c("World Bank High Income","World Bank Low Income", "World Bank Lower Middle Income","World Bank Upper Middle Income")
SDI_c <- c("Low SDI" ,"Low-middle SDI","Middle SDI" ,"High-middle SDI","High SDI")
main_c <- c('Global','High-income North America','Australasia','High-income Asia Pacific','Western Europe',
            'Southern Latin America','Eastern Europe','Central Europe','Central Asia',
            'Central Latin America','Andean Latin America','Caribbean',
            'Tropical Latin America','East Asia','Southeast Asia','Oceania',
            'North Africa and Middle East','South Asia','Southern Sub-Saharan Africa',
            'Western Sub-Saharan Africa','Eastern Sub-Saharan Africa',
            'Central Sub-Saharan Africa') ###
sum(location_a %in% main_c)
location_c <- location_a[!(location_a %in% c(income_c,SDI_c,main_c))]  ##

##
std_weights <- fread(input = 'F:/workforSCI/252/dataanalyse/GBD/population/WHO_age_stand.csv')
#colnames(std_weights) <- c("age","weight")

##
population <- fread(input = "F:/workforSCI/252/dataanalyse/GBD/population/WPP2024_Population1JanuaryByAge5GroupSex_Medium.csv.gz")

#
locationmatch <- data.frame(GBDlocation=c("Taiwan (Province of China)","Palestine","Turkey","Micronesia (Federated States of)","Democratic People's Republic of Korea"),
                            poplocation=c("China, Taiwan Province of China","State of Palestine","Türkiye","Micronesia","Dem. People's Republic of Korea"))
martch_location <-location_c
if (sum(unique(population$Location) %in% location_c)<204) {
  for (l in location_c[!location_c %in% unique(population$Location)]) {
    population[Location==locationmatch$poplocation[locationmatch$GBDlocation==l],10] <- l
  }
}
sum(unique(population$Location) %in% martch_location)
###
population <- population %>% filter(Location %in% c(martch_location,"Global",main_c)) %>% select(c(10,13,15,18,19,20)) 
unique(population$Location)

# 
pancreatitis_clean <- pancreatitis %>%
  mutate(age = str_remove(age, " years"))  # 

pancreatitis_45plus <- pancreatitis_clean %>%
  filter(age %in% target_ages, metric == "Rate")

std_weights_45plus <- std_weights %>%
  mutate(age = str_remove(`Age Group`, " years")) %>%  # 
  filter(age %in% target_ages) %>%
  rename(weight = `WHO world standard (&)`) %>%
  mutate(weight = weight / sum(weight))  # 

# 
asr_45plus <- pancreatitis_45plus %>%
  left_join(std_weights_45plus, by = "age") %>%
  group_by(measure, location, sex, year, cause) %>%
  summarise(
    asr_45plus = sum(val * weight, na.rm = TRUE),
    asr_45up = sum(upper * weight, na.rm = TRUE),
    asr_45low = sum(lower * weight, na.rm = TRUE),
    .groups = "drop"
  )%>% mutate(age = "Age-standardized",metric="Rate") %>%
  rename("val"="asr_45plus","upper"="asr_45up","lower"="asr_45low")

#####
population <- population %>%
  rename(
    year = Time,
    age = AgeGrp,
    location =Location

  )

population_long <- population %>%
  pivot_longer(
    cols = c(PopMale, PopFemale, PopTotal),
    names_to = "sex",
    values_to = "Population",
    names_prefix = "Pop"
  ) %>%
  mutate(
    sex = case_when(      # 
      sex == "Male" ~ "Male",
      sex == "Female" ~ "Female",
      sex == "Total" ~ "Both"
    )
  )

# 
population_long <- population_long %>%
  mutate(
    age = case_when(
    age %in% c("95-99", "100+") ~ "95+",
    TRUE ~ as.character(age)
    )
  ) %>%
  group_by(location, year, sex, age) %>%
  summarise(Population = sum(Population), .groups = "drop")

# 
age_45_plus <- target_ages

# 
population_45plus <- population_long %>% filter(age %in% target_ages)

# 
gbd_events <- pancreatitis_clean %>%   # 
  filter(age %in% age_45_plus,metric=="Number")  

gbd_events <- gbd_events %>%
  group_by(measure,location, year, age, sex) %>%
  summarise(
    Events = sum(val, na.rm = TRUE),
    eventupper = sum(upper, na.rm = TRUE),
    eventlower = sum(lower, na.rm = TRUE),
    .groups = "drop"
  )

# 
merged_data <- gbd_events %>%
  left_join(
    population_45plus, 
    by = c("location", "year", "age", "sex")
  )

# 
crude_rate_45plus <- merged_data %>%
  group_by(measure,location, year, sex) %>%    #
  summarise(
    TotalEvents = sum(Events, na.rm = TRUE),
    totlupper= sum(eventupper, na.rm = TRUE),
    totllower= sum(eventlower, na.rm = TRUE),
    TotalPopulation = sum(Population, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    CrudeRate = (TotalEvents / TotalPopulation),  #
    Crudeupper = (totlupper / TotalPopulation),
    Crudelower = (totllower / TotalPopulation)
  )
##
all_agesdata1 <- crude_rate_45plus %>% 
  select(c(1:7)) %>% rename(val=TotalEvents,upper=totlupper,lower=totllower)%>%
  mutate("metric"="Number")%>% mutate(age = "All ages",cause = diseasename)
all_agesdata2<- crude_rate_45plus %>% 
  select(c(1:4,9:11)) %>% rename(val=CrudeRate,upper=Crudeupper,lower=Crudelower)%>%
  mutate("metric"="Rate")%>% mutate(age = "All ages",cause = diseasename)
##
all_agesdata <- rbind(all_agesdata1,all_agesdata2)
colnames(all_agesdata)
colnames(asr_45plus) 
###
pancreatitis_age <- pancreatitis_clean %>% filter(age %in% target_ages)
# 
unique(pancreatitis_age$age)
unique(asr_45plus$age)
unique(all_agesdata$age)

pancreatitisdata <- bind_rows(pancreatitis_age, asr_45plus, all_agesdata)
colnames(pancreatitis)
unique(pancreatitisdata$age)
unique(pancreatitisdata$measure)
unique(pancreatitisdata$cause)

##
fwrite(x = pancreatitisdata,file = paste0(resultdir,aimfile))
gc()