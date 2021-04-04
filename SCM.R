# Author: Gordon Burtch
# Date: March 13th, 2021
# Purpose: Replication of Mitze, T., Kosfeld, R., Rode, J., & Wälde, K. (2020). Face masks considerably reduce COVID-19 cases in Germany. Proceedings of the National Academy of Sciences, 117(51), 32293-32301.

# Note: datasets and replication material are available at this github repository - https://github.com/gburtch/scm_replication

# Load libraries
library(tidyr)
library(dplyr)
library(ggplot2)
library(glmnet)
library(janitor)
library(Synth)
library(ggthemes)
library(patchwork)

# Background: Jena germany was the first German district to institute mandatory face mask requirements.
# They did this on April 6, 2020, weeks before the federal government imposed the same rule across the country.
# We are going to collect the same data, and try to replicate the result. 

### Let's import Covid-19 case volume data from german districts, by day..
covid <- read.csv("./Data Files/Covid_DE_cases_by_district.csv")
dist_names <- read.csv("./Data Files/Covid_DE_district_IDs.csv")
covid <- covid %>% merge(dist_names, by.x = "dist_id", by.y = "id")

# Let's convert date column into date format. 
covid$date <- as.Date(covid$date)
covid <- covid[order(covid$dist_id,covid$date),]

# We can pull in other features, though these are not actually required. 
# E.g., we can pull in physicians per capita in each district. 
dist_physicians <- read.csv("./Data Files/district_features/dist_phys_pcap.csv",sep=";") %>% rename("dist_id" = "Kennziffer", "physicians" = "Ärzte.je..Einwohner") %>% filter(!is.na(dist_id)) %>% select(c(dist_id,physicians))
dist_physicians$physicians <- as.numeric(gsub(",",".",dist_physicians$physicians))
covid <- covid %>% merge(dist_physicians, by="dist_id", all.x=TRUE)

# We now pull in pharmacy information.
dist_pharma <- read.csv("./Data Files/district_features/dist_pharmas_2017.csv",sep=";") %>% filter(!is.na(Kennziffer)) %>% rename("dist_id" = "Kennziffer","pharmacies" = "Apotheken")
dist_pharma$pharmacies <- as.numeric(gsub(",",".",dist_pharma$pharmacies))
dist_pharma <- dist_pharma %>% select(c(dist_id,pharmacies))
covid <- covid %>% merge(dist_pharma, by="dist_id",all.x=TRUE)

# We can pull in population information by age, and so on now as well.
dist_pop <- read.csv("./Data Files/district_features/dist_pop_age.csv", sep=";") %>% select(-c(dist_name))
dist_pop$year <- substr(dist_pop$year,1,4) 
# Some missing values here.
dist_pop <- dist_pop %>% lapply(as.integer) %>% as.data.frame()
# We have yearly values for 3 years - let's just use the most recent set of values in 2019.
dist_pop <- dist_pop %>% filter(year==2019) %>% select(-c(year))
covid <- covid %>% merge(dist_pop, by="dist_id", all.x=TRUE)

# Jena is district id 16053. 
covid$treat <- covid$dist_id==16053

# Let's trim to a reasonable window around the event date. 
covid <- subset(covid,date > as.numeric(as.Date("2020-03-01")) & date <= as.numeric(as.Date("2020-05-21")))

# Hard to see much here in the descriptive plot of the time series.
# Jena is there at the bottom, though hard to see. 
ggplot(data=covid,aes(x=date,y=log(cum_cases+1),color=factor(treat),group=dist_id,alpha=treat)) +
  geom_line() +
  geom_vline(xintercept=as.numeric(as.Date("2020-04-06")),color="red") + 
  xlab(expression(bold(paste("Date (2020)")))) +  
  ylab(expression(bold(paste("Logarithm of Cases")))) + 
  scale_alpha_manual(guide=FALSE,values=c(0.25,1))+
  scale_color_manual(name="District",labels=c("Others", "Jena"),values=c("gray","blue"))+
  ggtitle("Cumulative COVID-19 Cases Over Time") +
  theme_economist() +
  #Comment the below line out if you don't have Economica fonts installed.
  theme(text = element_text(family = "Economica", size = 10), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))+
  NULL

### Synthetic Control using LASSO

# First we pivot the data from long to wide, to use other districts' time series as predictors.
covid.wide <- covid %>% pivot_wider(id_cols=c("date"),names_from=c("district"),values_from=c("cum_cases","under.3.years","yr.3.to.under.6.years"))
covid.wide.train <- subset(covid.wide,date<as.numeric(as.Date("2020-04-06")))

# We have many more predictors than time periods now (~400 vs. 35), so we use LASSO for feature selection.
covid.wide.train.lasso <- remove_empty(covid.wide.train, which=c("rows","cols"))
covid.wide.train_mm <- model.matrix(`cum_cases_ SK Jena`~.,covid.wide.train.lasso)
lasso <- cv.glmnet(covid.wide.train_mm, covid.wide.train$`cum_cases_ SK Jena`, standardize=TRUE,alpha=1,nfolds=5)
ests <- as.matrix(coef(lasso,lasso$lambda.1se))

# Here are the non-zero control panels that lasso selected.
names(ests[ests!=0,])

# We can use the resulting control districts to create our 'synthetic control'. 
fml.rhs <- paste(c(names(ests[ests!=0,]))[2:length(names(ests[ests!=0,]))],collapse="+")
fml <- as.formula(paste("`cum_cases_ SK Jena`~",fml.rhs))
synth <- lm(data=covid.wide.train,formula=fml)

# Last, we can synthesize the resulting control series into the post treatment period. 
covid.wide$synth <- predict(synth,newdata = covid.wide)

# And, finally, we plot the comparison between synthetic and actual.
OLS_plot <- ggplot(data=covid.wide,aes(y=synth,x=date,linetype="dashed")) + geom_line() + 
  geom_line(aes(y=`cum_cases_ SK Jena`,x=date,linetype="solid")) +
  geom_vline(xintercept=as.numeric(as.Date("2020-04-06")),color="red") + 
  xlab(expression(bold(paste("Date (2020)")))) +  
  ylab(expression(bold(paste("Cumulative COVID-19 Cases")))) + 
  scale_linetype_manual(name="Series",values=c("dashed","solid"),labels=c("Synth","Jena, DE"))+
  ggtitle("Effect of Masks on COVID-19 (SCUL)") +
  theme_economist() +
  #Comment the below line out if you don't have Economica fonts installed.
  theme(text = element_text(family = "Economica", size = 10), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))+
  xlim(as.Date("2020-03-27"),as.Date("2020-04-26")) +
  ylim(0,225)+
  NULL

# Let's see how we did.
OLS_plot

## Vanilla Synthetic Control Method using the Synth package.

# The Synth package will implement non-negativity constraints on the control weights, and optimize control selection
# by minimizing MSPE. This package will also make it easy to include additional covariates as predictors (distict features).

# Make sure our data is a data frame.
covid <- as.data.frame(covid)

# Remove missing values.
covid <- subset(covid,complete.cases(covid)==TRUE)

# Renaming the outcome to Y, because Synth was complaining (throwing errors) with the original variable name.
covid <- covid %>% rename(Y = cum_cases)

# Cast the panel and date variables as numeric.
covid$date2 <- as.numeric(covid$date)
covid$dist_id <- as.numeric(covid$dist_id)

# Synth requires a list of control ID's - let's pull them out here.
# Note: 16053 is Jena

# We will also omit a few of Jena's neighboring districts and districts that instituted their own mask policies soon after Jena. 
# Specifically, we will exclude 16071; Weimarer Land, 16062; Nordhausen, and 8325; Rottweil.
# Saale-Holzland can also be omitted, but we dropped it due to missing data on district features, earlier.
dist_ids <- unique(covid$dist_id)
control_ids <- dist_ids[dist_ids != 16053 & dist_ids != 16071 & dist_ids != 16062 & dist_ids != 8325]

# Now we can use Synth's data preparation package. 
dataprep.out=
  dataprep(foo = covid,
           dependent = "Y",
           unit.variable = "dist_id",
           time.variable = "date2",
           
           # The authors used a lot of seemingly irrelevant predictors
           # For example, average female age? Average male age? Why is gender important?
           # I am going to keep things simple here: pharmacies, physicians and elderly.
           predictors = c("pharmacies","physicians","yr.75.years.and.over"),
           predictors.op = "mean",
           
           # We can also predict using case volumes day before treatment and week before treatment.
           special.predictors = list(list("Y", 18356, "mean"),list("Y", 18350, "mean")),
           
           #which panel is treated?
           treatment.identifier = 16053,
           
           #which panels are we using to construct the synthetic control?
           # Controls here will be every other district.
           controls.identifier = control_ids,
           
           #what is the pre-treatment time period?
           #these numeric values correspond to 34 days before treatment.
           #the paper only uses the 14 days before treatment for some reason?
           time.predictors.prior = c(18323:18357),
           
           time.optimize.ssr = c(18323:18357),
           
           #name of panel units
           unit.names.variable = "district",
           
           #time period to generate the plot for.
           #paper only goes 20 days post treatment because other treatments started.
           #We will just see what this looks like, however. 
           time.plot = 18343:18403)

# And final, create the synthetic control.
# This will take a few minutes to run and identify the optimal weights.
synth.out = synth(dataprep.out)

# Synth's native plotting functions.
# Path.plot() plots the synthetic against the actual treated unit data. 
path.plot(dataprep.res = dataprep.out, synth.res = synth.out,Xlab="Date",Ylab="Cumulative COVID-19 Cases",Main="Comparison of Synth vs. Actual Cum. COVID-19 Cases in Jena, Germany")
abline(v=18358,lty=2,col="red")

# Let's now make a more attractive plot of the ATT?  
# One that lets us use proper date labels and things.

# Let's pull out the data from the result, to make our own nicer plots in ggplot of course
synth_data_out = data.frame(dataprep.out$Y0plot%*%synth.out$solution.w) 
date = as.numeric(row.names(synth_data_out))
plot.df = data.frame(y=covid$Y[covid$dist_id==16053 & covid$date2 %in% date])
plot.df$synth = synth_data_out$w.weight
plot.df$date <- covid$date[covid$dist_id==16053 & covid$date2 %in% date]
SCM_plot <- ggplot(plot.df,aes(y=y,x=date,linetype="solid")) + geom_line() + 
  geom_line(aes(y=synth,x=date,linetype="dashed")) +
  geom_vline(xintercept=18358,color="red") + 
  xlab(expression(bold(paste("Date (2020)")))) +  
  ylab(expression(bold(paste("Cumulative COVID-19 Cases")))) + 
  scale_linetype_manual(name="Series",values=c("dashed","solid"),labels=c("Synth","Jena, DE"))+
  ggtitle("Effect of Masks on COVID-19 (Proper Synth)") +
  theme_economist() +
  #Comment the below line out if you don't have Economica fonts installed.
  theme(text = element_text(family = "Economica", size = 10), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))+
  xlim(as.Date("2020-03-27"),as.Date("2020-04-26")) +
  ylim(0,225)+
  NULL

# Let's put the two plots side-by-side.
OLS_plot + SCM_plot
