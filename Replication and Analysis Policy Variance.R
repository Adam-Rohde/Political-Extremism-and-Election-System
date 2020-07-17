library(haven)
library(lmtest)
library(sensemakr)
library(plm)
library(tidyverse)
library(miceadds)
library(estimatr)
library(rdrobust)

################################################################################
# TABLE 4 - Runoff and Policy Volatility, RDD Estimates
################################################################################




##############################
# Intertemporal variation
##############################


################
#prep 

data = read_dta("..\\dataset\\dual_ballot_replication.dta")

data = data %>% 
  mutate("counter" = 1) %>% 
  mutate(id =  cumsum(counter)) %>%
  select(-c(counter,var_ord))

#calculate var_ord
auxiliary = data %>% 
  mutate("avg" = ifelse(is.na(avg_ordinaria),beg_ordinaria,avg_ordinaria)) %>% 
  group_by(id_city_istat) %>% #group by city
  summarise(var_ord = sd(avg, na.rm=TRUE)^2) %>% #take sd and square to get variance
  ungroup()
data2 = left_join(data, auxiliary, by="id_city_istat")

#prep data
data3 = data2 %>%
  arrange(id_city_istat,t15000) %>% #sort by city  
  mutate("counter" = 1) %>% #add counter
  group_by(id_city_istat,t15000) %>% #group by city and treatment value
  mutate("n" = cumsum(counter)) %>% #create cumulative counter
  filter(n==1) %>% #limit to one obs per city-treatment values
  select(-c(counter,n)) %>% #remove counter
  ungroup() %>% 
  select(-c(pop15000,pop15000_2,pop15000_3,pop15000_4,
            t15000,t15000_int1,t15000_int2,t15000_int3,t15000_int4))

#update population fields
data4 = data3 %>%
  mutate(pop15000 = pop_census1991-15000,
         pop15000_2 = pop15000^2,
         pop15000_3 = pop15000^3,
         pop15000_4 = pop15000^4,
         t15000 = ifelse(pop_census1991>15000,1,0),
         t15000_int1=t15000*pop15000,
         t15000_int2=t15000*pop15000_2,
         t15000_int3=t15000*pop15000_3,
         t15000_int4=t15000*pop15000_4) %>% 
  select(-c(area,alt_max,end_rev_transf_pc,income_pc,
            elderly_index,active_pop,family_size,duration,term_limit))

#update covariates
temp = data %>%
  group_by(id_city_istat) %>% 
  summarise(area = mean(area, na.rm=TRUE),
            alt_max = mean(alt_max, na.rm=TRUE),
            end_rev_transf_pc = mean(end_rev_transf_pc, na.rm=TRUE),
            income_pc = mean(income_pc, na.rm=TRUE),
            elderly_index = mean(elderly_index, na.rm=TRUE),
            active_pop = mean(active_pop, na.rm=TRUE),
            family_size = mean(family_size, na.rm=TRUE),
            duration = mean(duration, na.rm=TRUE),
            term_limit = mean(term_limit)) %>% 
  ungroup()
TimeVarData = left_join(data4, temp, by="id_city_istat")

#reg function
TimeVarReg = function(formula,dta){
  wgt__ = NULL
  model   = lm.cluster(data=dta,formula=formula,cluster="id_city_istat")
  results = summary(model)
  sense   = sensemakr(lm(data=dta, formula=formula), treatment = "t15000")
  output  = matrix(0,nrow=9,ncol=1)
  output[1,1] = results[2,1]
  output[2,1] = results[2,2]
  output[3,1] = results[2,3]
  output[4,1] = model$lm_res$df.residual+model$lm_res$rank
  output[5,1] = sense$sensitivity_stats$se
  output[6,1] = sense$sensitivity_stats$t_statistic
  output[7,1] = sense$sensitivity_stats$r2yd.x
  output[8,1] = sense$sensitivity_stats$rv_q
  output[9,1] = sense$sensitivity_stats$rv_qa
  return(output)
}




################
#no covariates regs

output1 = matrix(0,nrow=9,ncol=12)
optimalBW1 = rdrobust(TimeVarData$var_ord,TimeVarData$pop_census,c=15000)$bws[1]
colnames(output1) = c("Spline Third","Spline Fourth","Spline Second","LLR (2000)","LLR (1750)","LLR (1500)","LLR (1250)","LLR (1000)","LLR (750)","LLR (500)","LLR (250)",paste0("LLR (Optimal = ",round(optimalBW1),")"))
rownames(output1) = c("Est.","Robust Std. Error","t-Value","Observations",
                      "Std. Std. Error","Std. t-Value","Partial R2 of treatment with outcome","Robustness Value, q = 1","Robustness Value, q = 1, alpha = 0.05")

#Spline Third
limited = TimeVarData
formula = "var_ord ~ t15000 + pop15000 + pop15000_2 + pop15000_3 + t15000_int1 + t15000_int2 + t15000_int3"
lm.cluster(data=limited,formula=formula,cluster="id_city_istat")
output1[,1] = TimeVarReg(formula,limited)

#Spline Fourth
limited = TimeVarData
formula = "var_ord ~ t15000 + pop15000 + pop15000_2 + pop15000_3 + pop15000_4 + t15000_int1 + t15000_int2 + t15000_int3 + t15000_int4"
output1[,2] = TimeVarReg(formula,limited)

#Spline Second
limited = TimeVarData
formula = "var_ord ~ t15000 + pop15000 + pop15000_2 + t15000_int1 + t15000_int2"
output1[,3] = TimeVarReg(formula,limited)

#LLR
j=4
for(i in rev(seq(250,2000,250))){
  print(i)
  start = 15000-i
  end = 15000+i
  limited = TimeVarData %>% filter(pop_census1991>=start & pop_census1991<=end)
  formula = "var_ord ~ t15000 + pop15000 + t15000_int1"
  output1[,j] = TimeVarReg(formula,limited)
  j = j + 1
}

#Optimal BW
start = 15000-optimalBW1
end = 15000+optimalBW1
limited = TimeVarData %>% filter(pop_census1991>=start & pop_census1991<=end)
formula = "var_ord ~ t15000 + pop15000 + t15000_int1"
output1[,12] = TimeVarReg(formula,limited)




################
#covariates regs

covariates = "+ north + CE + area + alt_max + end_rev_transf_pc + income_pc + elderly_index + active_pop + family_size + duration + term_limit"

output2 = matrix(0,nrow=9,ncol=12)
colnames(output2) = c("Spline Third","Spline Fourth","Spline Second","LLR (2000)","LLR (1750)","LLR (1500)","LLR (1250)","LLR (1000)","LLR (750)","LLR (500)","LLR (250)",paste0("LLR (Optimal = ",round(optimalBW1),")"))
rownames(output2) = c("Est.","Robust Std. Error","t-Value","Observations",
                      "Std. Std. Error","Std. t-Value","Partial R2 of treatment with outcome","Robustness Value, q = 1","Robustness Value, q = 1, alpha = 0.05")

#Spline Third
limited = TimeVarData
formula = paste("var_ord ~ t15000 + pop15000 + pop15000_2 + pop15000_3 + t15000_int1 + t15000_int2 + t15000_int3",covariates)
output2[,1] = TimeVarReg(formula,limited)

#Spline Fourth
limited = TimeVarData
formula = paste("var_ord ~ t15000 + pop15000 + pop15000_2 + pop15000_3 + pop15000_4 + t15000_int1 + t15000_int2 + t15000_int3 + t15000_int4",covariates)
output2[,2] = TimeVarReg(formula,limited)

#Spline Second
limited = TimeVarData
formula = paste("var_ord ~ t15000 + pop15000 + pop15000_2 + t15000_int1 + t15000_int2",covariates)
output2[,3] = TimeVarReg(formula,limited)

#LLR
j=4
for(i in rev(seq(250,2000,250))){
  print(i)
  start = 15000-i
  end = 15000+i
  limited = TimeVarData %>% filter(pop_census1991>=start & pop_census1991<=end)
  formula = paste("var_ord ~ t15000 + pop15000 + t15000_int1",covariates)
  output2[,j] = TimeVarReg(formula,limited)
  j = j + 1
}

#Optimal BW
start = 15000-optimalBW1
end = 15000+optimalBW1
limited = TimeVarData %>% filter(pop_census1991>=start & pop_census1991<=end)
formula = paste("var_ord ~ t15000 + pop15000 + t15000_int1",covariates)
output2[,12] = TimeVarReg(formula,limited)








##############################
# Cross sectional variation
##############################


################
#prep 

data = read_dta("..\\dataset\\dual_ballot_replication.dta")

data = data %>% 
  mutate("counter" = 1) %>% 
  mutate(id =  cumsum(counter)) %>%
  select(-c(counter,var_ord))


#create bins
data$bin100 = NA
for(i in seq(-5000, 4900, 100)){
  data$bin100 = ifelse(data$pop15000>=i & data$pop15000<i+100,i,data$bin100)
}

#create size100
temp = data %>%
  group_by(bin100) %>% 
  count(bin100) %>% 
  ungroup() %>%
  rename(size100=n) 
data2 = left_join(data, temp, by="bin100")

#create var3_ordinaria
auxiliary = data2 %>% 
  group_by(bin100,year_election) %>% 
  summarise(var2_ordinaria = sd(ordinaria_cs, na.rm=TRUE)^2) %>% 
  ungroup()
data3 = left_join(data2, auxiliary, by=c("bin100","year_election"))

auxiliary2 = data3 %>% 
  group_by(bin100) %>% 
  summarise(var3_ordinaria = mean(var2_ordinaria, na.rm=TRUE)) %>% 
  ungroup() 
data3 = left_join(data3, auxiliary2, by="bin100")
data3 = data3 %>% 
  arrange(bin100,id) %>% 
  select(-c(pop15000,pop15000_2,pop15000_3,pop15000_4,
            t15000_int1,t15000_int2,t15000_int3,t15000_int4))

#update population fields
data4 = data3 %>%
  mutate(pop15000 = bin100,
         pop15000_2 = bin100^2,
         pop15000_3 = bin100^3,
         pop15000_4 = bin100^4,
         t15000_int1=t15000*pop15000,
         t15000_int2=t15000*pop15000_2,
         t15000_int3=t15000*pop15000_3,
         t15000_int4=t15000*pop15000_4) %>% 
  select(-c(area,alt_max,north,south,CE))

#update covariates
temp = data %>%
  group_by(bin100) %>% 
  summarise(area = mean(area, na.rm=TRUE),
            alt_max = mean(alt_max, na.rm=TRUE),
            north = mean(north, na.rm=TRUE),
            south = mean(south, na.rm=TRUE),
            CE = mean(CE, na.rm=TRUE)) %>% 
  ungroup()
data5 = left_join(data4, temp, by="bin100")


#limit to one observation per bin
#generate weights
CrossSecData = data5 %>%
  arrange(bin100) %>% 
  mutate("counter" = 1) %>% 
  group_by(bin100) %>% 
  mutate("n" = cumsum(counter)) %>% 
  filter(n==1) %>% 
  select(-c(counter,n)) %>% 
  ungroup() %>%
  mutate(w = 1/size100)


#reg function
CrossSecReg = function(formula,dta){
  model   = lm(data=dta,formula=formula,weights=w)
  results = coeftest(model, vcov=vcovHC(model, type = "HC1"))
  sense   = sensemakr(model, treatment = "t15000")
  output  = matrix(0,nrow=9,ncol=1)
  output[1,1] = results[2,1]
  output[2,1] = results[2,2]
  output[3,1] = results[2,3]
  output[4,1] = model$df.residual+model$rank
  output[5,1] = sense$sensitivity_stats$se
  output[6,1] = sense$sensitivity_stats$t_statistic
  output[7,1] = sense$sensitivity_stats$r2yd.x
  output[8,1] = sense$sensitivity_stats$rv_q
  output[9,1] = sense$sensitivity_stats$rv_qa
  return(output)
}


################
#no covariates regs

output3 = matrix(0,nrow=9,ncol=12)
optimalBW2 = rdrobust(CrossSecData$var3_ordinaria,CrossSecData$pop_census,c=15000)$bws[1]
colnames(output3) = c("Spline Third","Spline Fourth","Spline Second","LLR (2000)","LLR (1750)","LLR (1500)","LLR (1250)","LLR (1000)","LLR (750)","LLR (500)","LLR (300)",paste0("LLR (Optimal = ",round(optimalBW2),")"))
rownames(output3) = c("Est.","Robust Std. Error","t-Value","Observations",
                      "Std. Std. Error","Std. t-Value","Partial R2 of treatment with outcome","Robustness Value, q = 1","Robustness Value, q = 1, alpha = 0.05")

#Spline Third
limited = CrossSecData
formula = "var3_ordinaria ~ t15000 + pop15000 + pop15000_2 + pop15000_3 + t15000_int1 + t15000_int2 + t15000_int3"
output3[,1] = CrossSecReg(formula,limited)

#Spline Fourth
limited = CrossSecData
formula = "var3_ordinaria ~ t15000 + pop15000 + pop15000_2 + pop15000_3 + pop15000_4 + t15000_int1 + t15000_int2 + t15000_int3 + t15000_int4"
output3[,2] = CrossSecReg(formula,limited)

#Spline Second
limited = CrossSecData
formula = "var3_ordinaria ~ t15000 + pop15000 + pop15000_2 + t15000_int1 + t15000_int2"
output3[,3] = CrossSecReg(formula,limited)


#LLR
j=4
for(i in rev(seq(250,2000,250))){
  print(i)
  if(i==250){i=300}
  start = 15000-i
  end = 15000+i
  limited = CrossSecData %>% filter(pop_census>=start & pop_census<=end)
  formula = "var3_ordinaria ~ t15000 + pop15000 + t15000_int1"
  output3[,j] = CrossSecReg(formula,limited)
  j = j + 1
}


#Optimal BW
start = 15000-optimalBW2
end = 15000+optimalBW2
limited = CrossSecData %>% filter(pop_census>=start & pop_census<=end)
formula = "var3_ordinaria ~ t15000 + pop15000 + t15000_int1"
output3[,12] = CrossSecReg(formula,limited)







################
#covariates regs

covariates = "+ north + CE + area + alt_max"

output4 = matrix(0,nrow=9,ncol=12)
colnames(output4) = c("Spline Third","Spline Fourth","Spline Second","LLR (2000)","LLR (1750)","LLR (1500)","LLR (1250)","LLR (1000)","LLR (750)","LLR (500)","LLR (450)",paste0("LLR (Optimal = ",round(optimalBW2),")"))
rownames(output4) = c("Est.","Robust Std. Error","t-Value","Observations",
                      "Std. Std. Error","Std. t-Value","Partial R2 of treatment with outcome","Robustness Value, q = 1","Robustness Value, q = 1, alpha = 0.05")

#Spline Third
limited = CrossSecData
formula = paste("var3_ordinaria ~ t15000 + pop15000 + pop15000_2 + pop15000_3 + t15000_int1 + t15000_int2 + t15000_int3",covariates)
output4[,1] = CrossSecReg(formula,limited)

#Spline Fourth
limited = CrossSecData
formula = paste("var3_ordinaria ~ t15000 + pop15000 + pop15000_2 + pop15000_3 + pop15000_4 + t15000_int1 + t15000_int2 + t15000_int3 + t15000_int4",covariates)
output4[,2] = CrossSecReg(formula,limited)

#Spline Second
limited = CrossSecData
formula = paste("var3_ordinaria ~ t15000 + pop15000 + pop15000_2 + t15000_int1 + t15000_int2",covariates)
output4[,3] = CrossSecReg(formula,limited)


#LLR
j=4
for(i in rev(seq(250,2000,250))){
  print(i)
  if(i==250){i=450}
  start = 15000-i
  end = 15000+i
  limited = CrossSecData %>% filter(pop_census>=start & pop_census<=end)
  formula = paste("var3_ordinaria ~ t15000 + pop15000 + t15000_int1",covariates)
  output4[,j] = CrossSecReg(formula,limited)
  j = j + 1
}

#Optimal BW
start = 15000-optimalBW2
end = 15000+optimalBW2
limited = CrossSecData %>% filter(pop_census>=start & pop_census<=end)
formula = paste("var3_ordinaria ~ t15000 + pop15000 + t15000_int1",covariates)
output4[,12] = CrossSecReg(formula,limited)





##############################
# Combine Output
##############################

library(data.table)
library(xtable)

o1 = data.frame(t(output1))
o2 = data.frame(t(output2))
o3 = data.frame(t(output3))
o4 = data.frame(t(output4))

o1 = data.frame(setDT(o1, keep.rownames = TRUE)[])
o2 = data.frame(setDT(o2, keep.rownames = TRUE)[])
o3 = data.frame(setDT(o3, keep.rownames = TRUE)[])
o4 = data.frame(setDT(o4, keep.rownames = TRUE)[])

o1["Covariates"] = "Without Covariates"
o2["Covariates"] = "With Covariates"
o3["Covariates"] = "Without Covariates"
o4["Covariates"] = "With Covariates"

o1["Variance"] = "Time"
o2["Variance"] = "Time"
o3["Variance"] = "Cross-Sectional"
o4["Variance"] = "Cross-Sectional"

o1 = o1 %>% select(Variance,Covariates, everything())
o2 = o2 %>% select(Variance,Covariates, everything())
o3 = o3 %>% select(Variance,Covariates, everything())
o4 = o4 %>% select(Variance,Covariates, everything())


output = rbind(o1,o2,o3,o4)
output[,4] = paste0(round(output[,4], 3))
output[,5] = paste0(round(output[,5], 3))
output[,6] = paste0(round(output[,6], 3))
output[,7] = paste0(round(output[,7], 0))


output[,8] = paste0(round(output[,8], 3))
output[,9] = paste0(round(output[,9], 3))

output[,10] = paste(round(100*output[,10], 1), "%", sep="")
output[,11] = paste(round(100*output[,11], 1), "%", sep="")
output[,12] = paste(round(100*output[,12], 1), "%", sep="")


xtable(output)






##############################
# Scatter and Lowess
##############################




#################
#Time Var

dataset = TimeVarData[complete.cases(TimeVarData[,"var_ord"]),]
x = dataset$pop15000
c = 0
label = "Normalized Population"
y <- dataset$var_ord


### Bins: ###
breaks0 = seq(c-5000,c+5000,by=40)#250 bins
bins0 = cut(x,breaks=breaks0)
bins = tapply(x,bins0,mean)
bin.y = tapply(y,bins0,function(x){return(mean(x,na.rm = T))})


#all data together
fit0 = loess(y~x,span = 0.25)
xx = x[order(x)]
ypredict = predict(fit0,xx, se.fit = T)
yy = ypredict

png(file="..\\output\\Lowess - Pop vs Time Var - Together.png",width=600, height=350)
plot(bins,bin.y,pch = 20, col=ifelse(bins<=c,"green4","purple2"),
     xlim=c(-5000,5000),ylim=c(0,2),
     xlab =label, ylab = 'Intertemporal Business Tax Rate Variance',
     main = "Intertemporal Business Tax Rate Variance vs Population\nLowess Fit with Span 0.25 for Entire Population Range")
lines(xx[xx<c], yy[xx<c], type='l',
      lwd=2, col = 'red4')
lines(xx[xx>c],yy[xx>c],col="blue4",lwd=2)
abline(v=c,lwd=2,col=1,lty=2)
dev.off()

#data split at cutoff
x1 <- x[x<c]; x2 <- x[x>=c]
y1 <- y[x<c]; y2 <- y[x>=c]
mod1 <- loess(y1~x1,span = 0.5)
mod2 <- loess(y2~x2,span = 0.5)
xx1 <- x1[order(x1)]
yy1 <- mod1$fitted[order(x1)]
xx2 <- x2[order(x2)]
yy2 <- mod2$fitted[order(x2)]

png(file="..\\output\\Lowess - Pop vs Time Var - Seperate.png",width=600, height=350)
plot(bins,bin.y,pch = 20, col=ifelse(bins<=c,"green4","purple2"),
     xlim=c(-5000,5000),ylim=c(0,2),
     xlab =label, ylab = 'Intertemporal Business Tax Rate Variance',
     main = "Intertemporal Business Tax Rate Variance vs Population\nLowess Fit with Span 0.5 for Treated and Controls Seperately")
lines(xx1, yy1,lwd=2, col = 'red4')
lines(xx2,yy2,col="blue4",lwd=2)
abline(v=c,lwd=2,col=1,lty=2)
dev.off()





#################
#CS Var

dataset = CrossSecData[complete.cases(CrossSecData[,"var3_ordinaria"]),]
x = dataset$pop15000
c = 0
label = "Normalized Population"
y <- dataset$var3_ordinaria


### Bins: ###
breaks0 = seq(c-5000,c+5000,by=40)#250 bins
bins0 = cut(x,breaks=breaks0)
bins = tapply(x,bins0,mean)
bin.y = tapply(y,bins0,function(x){return(mean(x,na.rm = T))})


#all data together
fit0 = loess(y~x,span = 0.25)
xx = x[order(x)]
ypredict = predict(fit0,xx, se.fit = T)
yy = ypredict

png(file="..\\output\\Lowess - Pop vs CS Var - Together.png",width=600, height=350)
plot(bins,bin.y,pch = 20, col=ifelse(bins<=c,"green4","purple2"),
     xlim=c(-5000,5000),ylim=c(0,2),
     xlab =label, ylab = 'Cross-Sectional Business Tax Rate Variance',
     main = "Cross-Sectional Business Tax Rate Variance vs Population\nLowess Fit with Span 0.25 for Entire Population Range")
lines(xx[xx<=c], yy[xx<=c], type='l',
      lwd=2, col = 'red4')
lines(xx[xx>=c],yy[xx>=c],col="blue4",lwd=2)
abline(v=c,lwd=2,col=1,lty=2)
dev.off()

#data split at cutoff
x1 <- x[x<=c]; x2 <- x[x>=c]
y1 <- y[x<=c]; y2 <- y[x>=c]
mod1 <- loess(y1~x1,span = 0.5)
mod2 <- loess(y2~x2,span = 0.5)
xx1 <- x1[order(x1)]
yy1 <- mod1$fitted[order(x1)]
xx2 <- x2[order(x2)]
yy2 <- mod2$fitted[order(x2)]

png(file="..\\output\\Lowess - Pop vs CS Var - Seperate.png",width=600, height=350)
plot(bins,bin.y,pch = 20, col=ifelse(bins<=c,"green4","purple2"),
     xlim=c(-5000,5000),ylim=c(0,2),
     xlab =label, ylab = 'Cross-Sectional Business Tax Rate Variance',
     main = "Cross-Sectional Business Tax Rate Variance vs Population\nLowess Fit with Span 0.5 for Treated and Controls Seperately")
lines(xx1, yy1,lwd=2, col = 'red4')
lines(xx2,yy2,col="blue4",lwd=2)
abline(v=c,lwd=2,col=1,lty=2)
dev.off()







##############################
# RDD Covariate Balance
##############################

#################
#Time Var

limited = TimeVarData %>% filter(pop_census1991>=14000 & pop_census1991<=16000)
covariates = "+ north + CE + area + alt_max + end_rev_transf_pc + income_pc + elderly_index + active_pop + family_size + duration + term_limit"

cov_out1 = matrix(0,nrow=8,ncol=4)
colnames(cov_out1) = c("Est.","SE","t-Value","p-Value")
rownames(cov_out1) = c("area","alt_max","end_rev_transf_pc",
                      "income_pc","elderly_index","active_pop","family_size","duration")

formula = paste("area ~ t15000 + pop15000 + t15000_int1",covariates,"- area")
cov_out1[1,] = summary(lm.cluster(data=limited,formula=formula,cluster="id_city_istat"))[2,]
formula = paste("alt_max ~ t15000 + pop15000 + t15000_int1",covariates,"- alt_max")
cov_out1[2,] = summary(lm.cluster(data=limited,formula=formula,cluster="id_city_istat"))[2,]
formula = paste("end_rev_transf_pc ~ t15000 + pop15000 + t15000_int1",covariates,"- end_rev_transf_pc")
cov_out1[3,] = summary(lm.cluster(data=limited,formula=formula,cluster="id_city_istat"))[2,]
formula = paste("income_pc ~ t15000 + pop15000 + t15000_int1",covariates,"- income_pc")
cov_out1[4,] = summary(lm.cluster(data=limited,formula=formula,cluster="id_city_istat"))[2,]
formula = paste("elderly_index ~ t15000 + pop15000 + t15000_int1",covariates,"- elderly_index")
cov_out1[5,] = summary(lm.cluster(data=limited,formula=formula,cluster="id_city_istat"))[2,]
formula = paste("active_pop ~ t15000 + pop15000 + t15000_int1",covariates,"- active_pop")
cov_out1[6,] = summary(lm.cluster(data=limited,formula=formula,cluster="id_city_istat"))[2,]
formula = paste("family_size ~ t15000 + pop15000 + t15000_int1",covariates,"- family_size")
cov_out1[7,] = summary(lm.cluster(data=limited,formula=formula,cluster="id_city_istat"))[2,]
formula = paste("duration ~ t15000 + pop15000 + t15000_int1",covariates,"- duration")
cov_out1[8,] = summary(lm.cluster(data=limited,formula=formula,cluster="id_city_istat"))[2,]


cov_out1 = round(cov_out1,2)

xtable(cov_out1)



#################
#CS Var

limited = CrossSecData %>% filter(pop_census>=14000 & pop_census<=16000)
covariates = "+ north + CE + area + alt_max"

cov_out2 = matrix(0,nrow=2,ncol=4)
colnames(cov_out2) = c("Est.","SE","t-Value","p-Value")
rownames(cov_out2) = c("area","alt_max")

formula = paste("area ~ t15000 + pop15000 + t15000_int1",covariates,"- area")
model = lm(data=limited,formula=formula,weights=w)
cov_out2[1,] = coeftest(model, vcov=vcovHC(model, type = "HC1"))[2,]
formula = paste("alt_max ~ t15000 + pop15000 + t15000_int1",covariates,"- alt_max")
model = lm(data=limited,formula=formula,weights=w)
cov_out2[2,] = coeftest(model, vcov=vcovHC(model, type = "HC1"))[2,]


cov_out2 = round(cov_out2,2)

xtable(cov_out2)





##############################
# Sensemakr
##############################
library(sensemakr)


#################
#Time Var


limited = TimeVarData %>% filter(pop_census1991>=14000 & pop_census1991<=16000)
sense.model <- sensemakr(var_ord ~ t15000 + pop15000 + t15000_int1 + north + CE + area + alt_max + end_rev_transf_pc + income_pc + elderly_index + active_pop + family_size + duration + term_limit,
                         treatment = "t15000",
                         benchmark = c("area","alt_max","end_rev_transf_pc","income_pc","elderly_index","active_pop","family_size","duration"),
                         kd = 5, data = limited)

plt = plot(sense.model)
r2dz.x = plt$r2dz.x
r2yz.dx = plt$r2yz.dx
value = plt$value
bounds = plt$bounds
lim=1
lim   <- min(max(c(0.4, r2dz.x*1.2)), 1 - 1e-12)
lim.y <- min(max(c(0.4, r2yz.dx*1.2)), 1 - 1e-12)
grid_values.x = seq(0, lim, by = lim/400)
grid_values.y = seq(0, lim.y, by = lim.y/400)
z_axis = value

png(file="..\\output\\Covariate Bal Sensemakr - Time Var - Total Effect.png",width=600, height=350)
contour(grid_values.x, grid_values.y, z_axis,nlevels = 20,
        main = "Sensemakr Contour Plot for Total Effect\nIntertemporal Variance of Business Property Tax  - LLR (h=1000)",
        xlab = expression(paste("Partial ", R^2, " of Z with D")),
        ylab = expression(paste("Partial ", R^2, " of Z with Y")))
contour(grid_values.x, grid_values.y,z_axis,level = 0,label = 0, col="red", lwd=2,lty = 2,add=TRUE)     
points(0, 0, pch = 17, col = "black", cex = 1)
points(bounds$r2dz.x, bounds$r2yz.dx, pch = 24, col = "blue", cex = 1)
xlabelbump = c(0,0,0,0,0,0,0,0)
ylabelbump = c(0,0,0,0,0,0,0,0)
text(bounds$r2dz.x+0.03+xlabelbump,bounds$r2yz.dx+0.03+ylabelbump, 
     labels=bounds$bound_label, cex=0.9, font=1)
dev.off()






sense.model <- sensemakr(var_ord ~ t15000 + pop15000 + t15000_int1 + north + CE + area + alt_max + end_rev_transf_pc + income_pc + elderly_index + active_pop + family_size + duration + term_limit,
                         treatment = "t15000",
                         benchmark = c("area","alt_max","end_rev_transf_pc","income_pc","elderly_index","active_pop","family_size","duration"),
                         kd = 3, data = limited)

plt = plot(sense.model,sensitivity.of = "t-value")
r2dz.x = plt$r2dz.x
r2yz.dx = plt$r2yz.dx
value = plt$value
bounds = plt$bounds
lim=1
lim   <- min(max(c(0.4, r2dz.x*1.2)), 1 - 1e-12)
lim.y <- min(max(c(0.4, r2yz.dx*1.2)), 1 - 1e-12)
grid_values.x = seq(0, lim, by = lim/400)
grid_values.y = seq(0, lim.y, by = lim.y/400)
z_axis = value

estimate <- sense.model$sensitivity_stats$estimate
q <- sense.model$info$q
reduce <- sense.model$info$reduce
dof <- sense.model$sensitivity_stats$dof
alpha <- sense.model$info$alpha
t.thr <- abs(qt(alpha/2, df = dof - 1))*sign(sense.model$sensitivity_stats$t_statistic)


png(file="..\\output\\Covariate Bal Sensemakr - Time Var - Significance.png",width=600, height=350)
contour(grid_values.x, grid_values.y, z_axis,nlevels = 20,
        main = "Sensemakr Contour Plot for Significance\nIntertemporal Variance of Business Property Tax  - LLR (h=1000)",
        xlab = expression(paste("Partial ", R^2, " of Z with D")),
        ylab = expression(paste("Partial ", R^2, " of Z with Y")))
contour(grid_values.x, grid_values.y,z_axis,level = t.thr,label = round(t.thr,2), lwd=2,col="red",lty = 2,add=TRUE)     
points(0, 0, pch = 17, col = "black", cex = 1)
points(bounds$r2dz.x, bounds$r2yz.dx, pch = 24, col = "blue", cex = 1)
xlabelbump = c(0,0,0,0,0,0,0,0)
ylabelbump = c(0,0,0,0,0,0,0,0)
text(bounds$r2dz.x+0.025+xlabelbump,bounds$r2yz.dx+0.025+ylabelbump, 
     labels=bounds$bound_label, cex=0.9, font=1)
dev.off()










#################
#CS Var

limited = CrossSecData %>% filter(pop_census>=14000 & pop_census<=16000)
sense.model <- sensemakr(var3_ordinaria ~ t15000 + pop15000 + t15000_int1 + north + CE + area + alt_max,
                         treatment = "t15000",
                         benchmark = c("area","alt_max"),
                         kd = 7, data = limited)

plt = plot(sense.model)
r2dz.x = plt$r2dz.x
r2yz.dx = plt$r2yz.dx
value = plt$value
bounds = plt$bounds
lim=1
lim   <- min(max(c(0.4, r2dz.x*1.2)), 1 - 1e-12)
lim.y <- min(max(c(0.4, r2yz.dx*1.2)), 1 - 1e-12)
grid_values.x = seq(0, lim, by = lim/400)
grid_values.y = seq(0, lim.y, by = lim.y/400)
z_axis = value

png(file="..\\output\\Covariate Bal Sensemakr - CS Var - Total Effect.png",width=600, height=350)
contour(grid_values.x, grid_values.y, z_axis,nlevels = 10,
        main = "Sensemakr Contour Plot for Total Effect\nCross-Sectional Variance of Business Property Tax  - LLR (h=1000)",
        xlab = expression(paste("Partial ", R^2, " of Z with D")),
        ylab = expression(paste("Partial ", R^2, " of Z with Y")))
contour(grid_values.x, grid_values.y,z_axis,level = 0,label = 0, col="red", lwd=2,lty = 2,add=TRUE)     
points(0, 0, pch = 17, col = "black", cex = 1)
points(bounds$r2dz.x, bounds$r2yz.dx, pch = 24, col = "blue", cex = 1)
xlabelbump = c(0,0,0,0,0,0,0,0)
ylabelbump = c(0,0,0,0,0,0,0,0)
text(bounds$r2dz.x+0.03+xlabelbump,bounds$r2yz.dx-0.03+ylabelbump, 
     labels=bounds$bound_label, cex=0.9, font=1)
dev.off()




sense.model <- sensemakr(var3_ordinaria ~ t15000 + pop15000 + t15000_int1 + north + CE + area + alt_max,
                         treatment = "t15000",
                         benchmark = c("area","alt_max"),
                         kd = 2, data = limited)

plt = plot(sense.model,sensitivity.of = "t-value")
r2dz.x = plt$r2dz.x
r2yz.dx = plt$r2yz.dx
value = plt$value
bounds = plt$bounds
lim=1
lim   <- min(max(c(0.4, r2dz.x*1.2)), 1 - 1e-12)
lim.y <- min(max(c(0.4, r2yz.dx*1.2)), 1 - 1e-12)
grid_values.x = seq(0, lim, by = lim/400)
grid_values.y = seq(0, lim.y, by = lim.y/400)
z_axis = value

estimate <- sense.model$sensitivity_stats$estimate
q <- sense.model$info$q
reduce <- sense.model$info$reduce
dof <- sense.model$sensitivity_stats$dof
alpha <- sense.model$info$alpha
t.thr <- abs(qt(alpha/2, df = dof - 1))*sign(sense.model$sensitivity_stats$t_statistic)


png(file="..\\output\\Covariate Bal Sensemakr - CS Var - Significance.png",width=600, height=350)
contour(grid_values.x, grid_values.y, z_axis,nlevels = 20,
        main = "Sensemakr Contour Plot for Significance\nCross-Sectional Variance of Business Property Tax  - LLR (h=1000)",
        xlab = expression(paste("Partial ", R^2, " of Z with D")),
        ylab = expression(paste("Partial ", R^2, " of Z with Y")))
contour(grid_values.x, grid_values.y,z_axis,level = t.thr,label = round(t.thr,2), lwd=2,col="red",lty = 2,add=TRUE)     
points(0, 0, pch = 17, col = "black", cex = 1)
points(bounds$r2dz.x, bounds$r2yz.dx, pch = 24, col = "blue", cex = 1)
xlabelbump = c(0,0,0,0,0,0,0,0)
ylabelbump = c(0,0,0,0,0,0,0,0)
text(bounds$r2dz.x+0.025+xlabelbump,bounds$r2yz.dx+0.025+ylabelbump, 
     labels=bounds$bound_label, cex=0.9, font=1)
dev.off()








