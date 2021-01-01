install.packages("RDS")
library(RDS)
library(dplyr)

# import dataset into R
dat1 <- read.csv("..path../nji.csv") # change the path to your own
glimpse(dat1)
anyDuplicated(dat1$SURID)    # make sure important variables no duplications
summary(dat1$SURID)    # make sure important variables have reasonable range, for instance, no "-1"
anyNA(dat1$SURID)   # make sure important variables no missing values, for instance, ID

##datasets in the package "RDS"
#faux: simple simulated RDS data
#fauxmadrona: no seed dependency, 50% sampling fraction, 5 Homophily, 1.8 differential activity
#fauxsycamore: seed dependency, 70% sampling fraction, 5 Homo, 1.8 w
#fauxtime: has interview date
data("faux")
glimpse(faux)
data("fauxmadrona")
glimpse(fauxmadrona)
fauxmadrona.network
data("fauxsycamore")
glimpse(fauxsycamore)
fauxsycamore.network
data("fauxtime")
glimpse(fauxtime)   #for HCG because there is time in the dataset

#data_1 is not a rds.data.frame; it is a r.data.frame
data_1 <- as.data.frame(faux)
class(data_1)
glimpse(data_1)
#make data_1 as rds.data.frame
data_rds <- as.rds.data.frame(data_1, id= "id", recruiter.id = "recruiter.id", 
                              network.size = "network.size")
class(data_rds)
glimpse(data_rds)
#what if you do not know recruiter id, you only know coupon id
fpath <- system.file("extdata", "nyjazz.csv", package="RDS")
dat <- read.csv(fpath)
glimpse(dat)
dat$recruiter.id <- rid.from.coupons(dat,subject.coupon = "own.coupon",
                                     coupon.variables=paste0("coupon.",1:7),
                                     subject.id = "id")
summary(dat$recruiter.id)

#understand your data
#any missing? some variables might have NA, for instance, coupon id for seeds
anyNA(fauxmadrona)
glimpse(fauxmadrona)
plot(fauxmadrona, plot.type = "Recruitment tree")
plot(fauxmadrona, plot.type = "Network size by wave")
plot(fauxmadrona, plot.type = "Recruitment tree", stratify.by = "disease")
plot(fauxmadrona, plot.type = "Network size by wave", stratify.by = "disease")
glimpse(fauxsycamore)
plot(fauxsycamore, plot.type = "Recruitment tree", stratify.by = "disease")
plot(fauxsycamore, plot.type = "Network size by wave", stratify.by = "disease")

#use proper estimators to estimate population based on "fauxmadrona"
RDS.I.estimates(fauxmadrona, outcome.variable = "disease")
#population size as 1000 to calculate CI; for point estimate, population size doesn't matter
RDS.II.estimates(fauxmadrona, outcome.variable = "disease")
RDS.SS.estimates(fauxmadrona, outcome.variable = "disease", N=1000)
#real mean is 0.259; SS is better

#diagnosis
#differential activity
#RDS-I
differential.activity.estimates(fauxmadrona, outcome.variable = "disease", 
                                weight.type = "RDS-I")
#RDS-II
differential.activity.estimates(fauxmadrona, outcome.variable = "disease", 
                                weight.type = "RDS-II")
#Gile SS
differential.activity.estimates(fauxmadrona, outcome.variable = "disease", 
                                weight.type = "Gile's SS")
#the real differential activity is 1.8

#homophily: population + recruitment
#RDS-I, recruitment
homophily.estimates(fauxmadrona, outcome.variable = "disease", weight.type = "RDS-I",
                    recruitment = TRUE)
#population
homophily.estimates(fauxmadrona, outcome.variable = "disease", weight.type = "RDS-I")
#RDS-II, recruitment
homophily.estimates(fauxmadrona, outcome.variable = "disease", weight.type = "RDS-II",
                    recruitment = TRUE)
#population
homophily.estimates(fauxmadrona, outcome.variable = "disease", weight.type = "RDS-II")
#Gile's SS, recruitment
homophily.estimates(fauxmadrona, outcome.variable = "disease", weight.type = "Gile's SS",
                    recruitment = TRUE)
#population
homophily.estimates(fauxmadrona, outcome.variable = "disease", weight.type = "Gile's SS")
#the real homophily is 5; all three estimators give lower estimates

#bottleneck plot
#RDS-I
bottleneck.plot(fauxmadrona, outcome.variable = "disease", est.func = RDS.I.estimates)
#RDS-II
bottleneck.plot(fauxmadrona, outcome.variable = "disease", est.func = RDS.II.estimates)
#Gile's SS
bottleneck.plot(fauxmadrona, outcome.variable = "disease", est.func = RDS.SS.estimates)

#convergence plot
#RDS-I
convergence.plot(fauxmadrona, outcome.variable = "disease", est.func = RDS.I.estimates)
#RDS-II
convergence.plot(fauxmadrona, outcome.variable = "disease", est.func = RDS.II.estimates)
#Gile's SS
convergence.plot(fauxmadrona, outcome.variable = "disease", est.func = RDS.SS.estimates)


#fauxsycamore
RDS.I.estimates(fauxsycamore, outcome.variable = "disease")
#population size as 1000 to calculate CI; for point estimate, population size doesn't matter
RDS.II.estimates(fauxsycamore, outcome.variable = "disease")
RDS.SS.estimates(fauxsycamore, outcome.variable = "disease", N=1000)
#real mean is 0.24; SS is better

#diagnosis
#differential activity
#RDS-I
differential.activity.estimates(fauxsycamore, outcome.variable = "disease", 
                                weight.type = "RDS-I")
#RDS-II
differential.activity.estimates(fauxsycamore, outcome.variable = "disease", 
                                weight.type = "RDS-II")
#Gile SS
differential.activity.estimates(fauxsycamore, outcome.variable = "disease", 
                                weight.type = "Gile's SS")
#the real differential activity is 1.8; seeds dependence affects RDS-II

#homophily: population + recruitment
#RDS-I, recruitment
homophily.estimates(fauxsycamore, outcome.variable = "disease", weight.type = "RDS-I",
                    recruitment = TRUE)
#population
homophily.estimates(fauxsycamore, outcome.variable = "disease", weight.type = "RDS-I")
#RDS-II, recruitment
homophily.estimates(fauxsycamore, outcome.variable = "disease", weight.type = "RDS-II",
                    recruitment = TRUE)
#population
homophily.estimates(fauxsycamore, outcome.variable = "disease", weight.type = "RDS-II")
#Gile's SS, recruitment
homophily.estimates(fauxsycamore, outcome.variable = "disease", weight.type = "Gile's SS",
                    recruitment = TRUE)
#population
homophily.estimates(fauxsycamore, outcome.variable = "disease", weight.type = "Gile's SS")
#the real homophily is 5; all three estimators give lower estimates

#bottleneck plot
#RDS-I
bottleneck.plot(fauxsycamore, outcome.variable = "disease", est.func = RDS.I.estimates)
#RDS-II
bottleneck.plot(fauxsycamore, outcome.variable = "disease", est.func = RDS.II.estimates)
#Gile's SS
bottleneck.plot(fauxsycamore, outcome.variable = "disease", est.func = RDS.SS.estimates)

#convergence plot
#RDS-I
convergence.plot(fauxsycamore, outcome.variable = "disease", est.func = RDS.I.estimates)
#RDS-II
convergence.plot(fauxsycamore, outcome.variable = "disease", est.func = RDS.II.estimates)
#Gile's SS
convergence.plot(fauxsycamore, outcome.variable = "disease", est.func = RDS.SS.estimates)


#what about HCG?
#use "fauxtime" dataset
glimpse(fauxtime)

#plot
plot(fauxtime, plot.type = "Recruitment tree", stratify.by = "var1")
#no obvious seed bias; homophily seems low
plot(fauxtime, plot.type = "Network size by wave", stratify.by = "var1")
#lower differential activity

#estimates
RDS.I.estimates(fauxtime, outcome.variable = "var1", N=1000)
#population size as 1000 to calculate CI; for point estimate, population size doesn't matter
RDS.II.estimates(fauxtime, outcome.variable = "var1", N=1000)
RDS.SS.estimates(fauxtime, outcome.variable = "var1", N=1000)
RDS.HCG.estimates(fauxtime, outcome.variable = "var1", N=1000)

#diagnosis
#differential activity
#RDS-I
differential.activity.estimates(fauxtime, outcome.variable = "var1", 
                                weight.type = "RDS-I")
#RDS-II
differential.activity.estimates(fauxtime, outcome.variable = "var1", 
                                weight.type = "RDS-II")
#Gile SS
differential.activity.estimates(fauxtime, outcome.variable = "var1", 
                                weight.type = "Gile's SS")
#HCG
differential.activity.estimates(fauxtime, outcome.variable = "var1", 
                                weight.type = "HCG")

#homophily: population + recruitment
#RDS-I, recruitment
homophily.estimates(fauxtime, outcome.variable = "var1", weight.type = "RDS-I",
                    recruitment = TRUE)
#population
homophily.estimates(fauxtime, outcome.variable = "var1", weight.type = "RDS-I")
#RDS-II, recruitment
homophily.estimates(fauxtime, outcome.variable = "var1", weight.type = "RDS-II",
                    recruitment = TRUE)
#population
homophily.estimates(fauxtime, outcome.variable = "var1", weight.type = "RDS-II")
#Gile's SS, recruitment
homophily.estimates(fauxtime, outcome.variable = "var1", weight.type = "Gile's SS",
                    recruitment = TRUE)
#population
homophily.estimates(fauxtime, outcome.variable = "var1", weight.type = "Gile's SS")
#HCG
homophily.estimates(fauxtime, outcome.variable = "var1", weight.type = "HCG",
                    recruitment = TRUE)
#population
homophily.estimates(fauxtime, outcome.variable = "var1", weight.type = "HCG")

#bottleneck plot
#RDS-I
bottleneck.plot(fauxtime, outcome.variable = "var1", est.func = RDS.I.estimates)
#RDS-II
bottleneck.plot(fauxtime, outcome.variable = "var1", est.func = RDS.II.estimates)
#Gile's SS
bottleneck.plot(fauxtime, outcome.variable = "var1", est.func = RDS.SS.estimates)
#HCG
bottleneck.plot(fauxtime, outcome.variable = "var1", est.func = RDS.HCG.estimates)

#convergence plot
#RDS-I
convergence.plot(fauxtime, outcome.variable = "var1", est.func = RDS.I.estimates)
#RDS-II
convergence.plot(fauxtime, outcome.variable = "var1", est.func = RDS.II.estimates)
#Gile's SS
convergence.plot(fauxtime, outcome.variable = "var1", est.func = RDS.SS.estimates)
#HCG
convergence.plot(fauxtime, outcome.variable = "var1", est.func = RDS.HCG.estimates)

