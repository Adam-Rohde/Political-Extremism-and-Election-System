library(haven)
library(lmtest)
library(sensemakr)
library(plm)
library(tidyverse)
library(miceadds)
library(estimatr)
library(rms)



################################################################################
#McCrary Tests
################################################################################

#Read Data
data = read_dta("..\\dataset\\dual_ballot_replication.dta")

#clean data
popdata = data %>% 
  select(c(id_city_istat,year_election,pop_census1991,pop_census)) %>% 
  mutate(cenyear = ifelse(year_election<=2001,"pop1991","pop2001")) %>%
  group_by(id_city_istat,cenyear,pop_census1991) %>% 
  summarise(pop_census = mean(pop_census, na.rm=TRUE)) %>% 
  ungroup() %>%
  spread(cenyear, pop_census) %>%
  select(-c(pop_census1991)) %>%
  #drop_na() %>%
  mutate(dif = (pop2001 - pop1991)/pop1991) %>%
  mutate(bucket1991 = cut(pop1991, breaks=seq(10000,20000,250)),
         bucket2001 = cut(pop2001, breaks=seq(10000,20000,250)),
         popnorm = pop1991-15000)

#population density plot
popPlotdata = popdata %>% drop_na()
png(file="..\\output\\population density plot.png",width=600, height=350)
plot(density(popPlotdata$pop1991),main = "Density of City Populations from 1991 and 2001 Census",
     xlab = "Population Size", ylab="Density",col = "red",ylim=c(0,0.0002))
lines(density(popPlotdata$pop2001),col = "blue")
legend(18000, 0.0002, legend=c("1991 Census", "2001 Census"),col=c("red", "blue"), lty=1:1, cex=0.8)
dev.off()

#1991 density
popdata_collapse1991 = popdata %>%
  group_by(bucket1991) %>% 
  count(bucket1991) %>% 
  ungroup() %>%
  mutate(bucket1991_label = as.numeric(substr(bucket1991, 2, str_locate(bucket1991, ",")-1))) %>%
  select(-c(bucket1991)) %>% 
  drop_na()
total = popdata_collapse1991 %>% summarise(total = sum(n))
total = as.numeric(total[1,1])
popdata_collapse1991 = popdata_collapse1991 %>% mutate(w = n/total)

#2001 density
popdata_collapse2001 = popdata %>%
  group_by(bucket2001) %>% 
  count(bucket2001) %>% 
  ungroup() %>%
  mutate(bucket2001_label = as.numeric(substr(bucket2001, 2, str_locate(bucket2001, ",")-1))) %>%
  select(-c(bucket2001)) %>% 
  drop_na()
total = popdata_collapse2001 %>% summarise(total = sum(n))
total = as.numeric(total[1,1])
popdata_collapse2001 = popdata_collapse2001 %>% mutate(w = n/total)

#density difference and pop buckets
dif = popdata_collapse2001$w - popdata_collapse1991$w 
pop = popdata_collapse1991$bucket1991_label

#Plot Density Difference and Population Buckets
#Plot Splines
png(file="..\\output\\McCrary Test for Sorting.png",
    width=600, height=350)

plot(pop,dif,main = "Testing for sorting between 1991 and 2001 Census",
     xlab = "Population Size", ylab="Density Difference 2001−1991")
abline(v=15000, col="red")

spline3a <- smooth.spline(pop[pop<=15000], dif[pop<=15000], df = 3)
res_a <- (spline3a$yin - spline3a$y)/(1-spline3a$lev) 
sigma_a <- sqrt(var(res_a)) 
upper_a <- spline3a$y + 2.0*sigma_a*sqrt(spline3a$lev)  
lower_a <- spline3a$y - 2.0*sigma_a*sqrt(spline3a$lev)  

spline3b <- smooth.spline(pop[pop>=15000], dif[pop>=15000], df = 3)
res_b <- (spline3b$yin - spline3b$y)/(1-spline3b$lev) 
sigma_b <- sqrt(var(res_b)) 
upper_b <- spline3b$y + 2.0*sigma_b*sqrt(spline3b$lev)  
lower_b <- spline3b$y - 2.0*sigma_b*sqrt(spline3b$lev)  

lines(spline3a, lty = 1, col = "blue")
lines(cbind(spline3a$x,upper_a), lty = 1, col = "darkgray")
lines(cbind(spline3a$x,lower_a), lty = 1, col = "darkgray")

lines(spline3b, lty = 1, col = "blue")
lines(cbind(spline3b$x,upper_b), lty = 1, col = "darkgray")
lines(cbind(spline3b$x,lower_b), lty = 1, col = "darkgray")

dev.off()

#actual McCrary tests
library(rdd)
png(file="..\\output\\McCrary 1991.png",width=600, height=350)
McC1991 = DCdensity(popPlotdata$pop1991,15000,htest = TRUE,plot=TRUE)
title(paste0("McCrary Test for Density of City Populations from 1991 Census\np-Value: ",round(McC1991$p.value,2)),
      ylab = "Density",xlab="Population Size")
dev.off()

png(file="..\\output\\McCrary 2001.png",width=600, height=350)
McC2001 = DCdensity(popPlotdata$pop2001,15000,htest = TRUE,plot=TRUE)
title(paste0("McCrary Test for Density of City Populations from 2001 Census\np-Value: ",round(McC2001$p.value,2)),
      ylab = "Density",xlab="Population Size")
dev.off()
