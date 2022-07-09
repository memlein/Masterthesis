library(VineCopula)
library(tidyverse)
library(copula)
library(ggplot2)
library(ggplotify)
library(doParallel)
library(stargazer)
library(rugarch)
library(MASS)
library(segMGarch)
library(GAS)
library(rmgarch)
library(data.table)
library(matrixStats)
library(fGarch)
library(purrr)
library(tidyquant)
library(parallelly)
library(viridis)
library(gridExtra)
library(FinTS)
library(latex2exp)
library(tseries)
library(stringi)
set.seed(2021)


# 1. Get Data, calculate logreturns and portfolio returns --------------------------

# Set ticker symbols to download
tickers <- c("^STOXX","^IBEX","^FTSE","^GDAXI", "CL=F","^TNX")

# get data from yahoo finance
prices <- tidyquant::tq_get(tickers,
                            from = "2009-01-01",
                            to = "2021-12-31",
                            get = "stock.prices")



### backup data

#write.csv(prices, file="backup_prices.csv")

# load data
#prices <- read_csv("backup_prices.csv",col_types = cols(...1 = col_skip(), date = col_date(format = "%Y-%m-%d")))

# remove NAs
prices <- prices%>%
  drop_na()

# calculate log-returns
daily_returns_stocks <- prices %>%
  group_by(symbol) %>%
  tq_transmute(adjusted, periodReturn, period = "daily", type="log")

# drop NAs
daily_returns_stocks <- daily_returns_stocks%>%drop_na()


# combine logreturns with portfolio returns into one dataframe
daily_returns_stocks <- daily_returns_stocks%>%
  pivot_wider(names_from = symbol, values_from=daily.returns)

daily_returns_stocks <- daily_returns_stocks%>%
  drop_na()

daily_returns_stocks <- daily_returns_stocks%>%
  relocate(date, .after=last_col())

# remove first row as it is 0
daily_returns_stocks <- daily_returns_stocks[2:nrow(daily_returns_stocks),]



# 2. Implement unweighted tests -------------------------------------------

test.unweighted <- list()

test.unweighted$tau <- "tau"

test.unweighted$AIC <- "AIC"

# CvM test is the S_n test in the thesis
test.unweighted$CvM.test <- function(u1,u2,weights){
  fitted_copula <- suppressWarnings(VineCopula::BiCopSelect(u1=u1,u2=u2,familyset = c(1,2,3,4,5,6)))
  if (fitted_copula$family %in% c(13)) {#180 Clayton
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::claytonCopula(dim=2,param = fitted_copula$par)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 13
  }
  else if (fitted_copula$family %in% c(23)) {#90 Clayton
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::claytonCopula(dim=2,param = -fitted_copula$par),flip=c(TRUE,FALSE)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 23
  }
  else if (fitted_copula$family %in% c(33)) {#270 Clayton
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::claytonCopula(dim=2,param = -fitted_copula$par),flip=c(FALSE,TRUE)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 33
  }
  else if (fitted_copula$family %in% c(3)){#Clayton
    result <- copula::gofCopula(copula=copula::claytonCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 3
  }
  else if (fitted_copula$family %in% c(14)) {#180 Gumbel
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::gumbelCopula(dim=2,param = fitted_copula$par)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 14
  }
  else if (fitted_copula$family %in% c(24)) {#90 Gumbel
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::gumbelCopula(dim=2,param = -fitted_copula$par),flip=c(TRUE,FALSE)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 24
  }
  else if (fitted_copula$family %in% c(34)) {#270 Gumbel
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::gumbelCopula(dim=2,param = -fitted_copula$par),flip=c(FALSE,TRUE)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 34
  }
  else if (fitted_copula$family %in% c(4)){#Gumbel
    result <- copula::gofCopula(copula=copula::gumbelCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 4
  }
  else if (fitted_copula$family %in% c(1)){#Normal
    result <- copula::gofCopula(copula=copula::ellipCopula(family="normal",dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 1
  }
  else if (fitted_copula$family %in% c(2)){#Student
    result <- copula::gofCopula(copula=copula::ellipCopula(family="t",dim=2,param=fitted_copula$par,df=round(fitted_copula$par2),df.fixed=TRUE), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 2
  }
  else if (fitted_copula$family %in% c(5)){#Frank
    result <- copula::gofCopula(copula=copula::frankCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 5
  }
  else if (fitted_copula$family %in% c(16)) {#180 Joe
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::joeCopula(dim=2,param = fitted_copula$par)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 16
  }
  else if (fitted_copula$family %in% c(26)) {#90 Joe
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::joeCopula(dim=2,param = -fitted_copula$par),flip=c(TRUE,FALSE)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 26
  }
  else if (fitted_copula$family %in% c(36)) {#270 Joe
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::joeCopula(dim=2,param = -fitted_copula$par),flip=c(FALSE,TRUE)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 36
  }
  else if (fitted_copula$family %in% c(6)){#Joe
    result <- copula::gofCopula(copula=copula::joeCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="pb")$p.value
    fam <- 6
  }
  #print(fam)
  result
}

test.unweighted$Rn.test <- function(u1,u2,weights){
  fitted_copula <- suppressWarnings(VineCopula::BiCopSelect(u1=u1,u2=u2,familyset = c(1,2,3,4,5,6)))
  if (fitted_copula$family %in% c(13)) {#180 Clayton
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::claytonCopula(dim=2,param = fitted_copula$par)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    fam <- 13
  }
  else if (fitted_copula$family %in% c(23)) {#90 Clayton
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::claytonCopula(dim=2,param = -fitted_copula$par),flip=c(TRUE,FALSE)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    fam <- 23
  }
  else if (fitted_copula$family %in% c(33)) {#270 Clayton
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::claytonCopula(dim=2,param = -fitted_copula$par),flip=c(FALSE,TRUE)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    fam <- 33
  }
  else if (fitted_copula$family %in% c(3)){#Clayton
    result <- copula::gofCopula(copula=copula::claytonCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    fam <- 3
  }
  else if (fitted_copula$family %in% c(14)) {#180 Gumbel
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::gumbelCopula(dim=2,param = fitted_copula$par)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    fam <- 14
  }
  else if (fitted_copula$family %in% c(24)) {#90 Gumbel
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::gumbelCopula(dim=2,param = -fitted_copula$par),flip=c(TRUE,FALSE)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    fam <- 24
  }
  else if (fitted_copula$family %in% c(34)) {#270 Gumbel
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::gumbelCopula(dim=2,param = -fitted_copula$par),flip=c(FALSE,TRUE)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    fam <- 34
  }
  else if (fitted_copula$family %in% c(4)){#Gumbel
    result <- copula::gofCopula(copula=copula::gumbelCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    fam <- 4
  }
  else if (fitted_copula$family %in% c(1)){#Normal
    result <- copula::gofCopula(copula=copula::ellipCopula(family="normal",dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    fam <- 1
  }
  else if (fitted_copula$family %in% c(2)){#Student
    result <- copula::gofCopula(copula=copula::ellipCopula(family="t",dim=2,param=fitted_copula$par,df=round(fitted_copula$par2),df.fixed=TRUE), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    fam <- 2
  }
  else if (fitted_copula$family %in% c(5)){#Frank
    result <- copula::gofCopula(copula=copula::frankCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    fam <- 5
  }
  else if (fitted_copula$family %in% c(16)) {#180 Joe
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::joeCopula(dim=2,param = fitted_copula$par)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    fam <- 16
  }
  else if (fitted_copula$family %in% c(26)) {#90 Joe
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::joeCopula(dim=2,param = -fitted_copula$par),flip=c(TRUE,FALSE)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    fam <- 26
  }
  else if (fitted_copula$family %in% c(36)) {#270 Joe
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::joeCopula(dim=2,param = -fitted_copula$par),flip=c(FALSE,TRUE)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    fam <- 36
  }
  else if (fitted_copula$family %in% c(6)){#Joe
    result <- copula::gofCopula(copula=copula::joeCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    fam <- 6
  }
  #print(fam)
  result
}

test.unweighted$CvM.tau.mult <- function(u1,u2,weights){
  fitted_copula <- suppressWarnings(VineCopula::BiCopSelect(u1=u1,u2=u2,familyset = c(1,2,3,4,5,6)))
  if (fitted_copula$family %in% c(13)) {#180 Clayton
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::claytonCopula(dim=2,param = fitted_copula$par)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 13
  }
  else if (fitted_copula$family %in% c(23)) {#90 Clayton
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::claytonCopula(dim=2,param = -fitted_copula$par),flip=c(TRUE,FALSE)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 23
  }
  else if (fitted_copula$family %in% c(33)) {#270 Clayton
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::claytonCopula(dim=2,param = -fitted_copula$par),flip=c(FALSE,TRUE)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 33
  }
  else if (fitted_copula$family %in% c(3)){#Clayton
    result <- copula::gofCopula(copula=copula::claytonCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 3
  }
  else if (fitted_copula$family %in% c(14)) {#180 Gumbel
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::gumbelCopula(dim=2,param = fitted_copula$par)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 14
  }
  else if (fitted_copula$family %in% c(24)) {#90 Gumbel
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::gumbelCopula(dim=2,param = -fitted_copula$par),flip=c(TRUE,FALSE)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 24
  }
  else if (fitted_copula$family %in% c(34)) {#270 Gumbel
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::gumbelCopula(dim=2,param = -fitted_copula$par),flip=c(FALSE,TRUE)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 34
  }
  else if (fitted_copula$family %in% c(4)){#Gumbel
    result <- copula::gofCopula(copula=copula::gumbelCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 4
  }
  else if (fitted_copula$family %in% c(1)){#Normal
    result <- copula::gofCopula(copula=copula::ellipCopula(family="normal",dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 1
  }
  else if (fitted_copula$family %in% c(2)){#Student
    result <- copula::gofCopula(copula=copula::ellipCopula(family="t",dim=2,param=fitted_copula$par,df=round(fitted_copula$par2),df.fixed=TRUE), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 2
  }
  else if (fitted_copula$family %in% c(5)){#Frank
    result <- copula::gofCopula(copula=copula::frankCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 5
  }
  else if (fitted_copula$family %in% c(16)) {#180 Joe
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::joeCopula(dim=2,param = fitted_copula$par)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 16
  }
  else if (fitted_copula$family %in% c(26)) {#90 Joe
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::joeCopula(dim=2,param = -fitted_copula$par),flip=c(TRUE,FALSE)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 26
  }
  else if (fitted_copula$family %in% c(36)) {#270 Joe
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::joeCopula(dim=2,param = -fitted_copula$par),flip=c(FALSE,TRUE)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 36
  }
  else if (fitted_copula$family %in% c(6)){#Joe
    result <- copula::gofCopula(copula=copula::joeCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="pb")$p.value
    fam <- 6
  }
  #print(fam)
  result*abs(cor(u1, u2, method = "kendall", use = "complete.obs"))
}

test.unweighted$Rn.tau.mult <- function(u1,u2,weights){
  fitted_copula <- suppressWarnings(VineCopula::BiCopSelect(u1=u1,u2=u2,familyset = c(1,2,3,4,5,6)))
  if (fitted_copula$family %in% c(13)) {#180 Clayton
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::claytonCopula(dim=2,param = fitted_copula$par)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    fam <- 13
  }
  else if (fitted_copula$family %in% c(23)) {#90 Clayton
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::claytonCopula(dim=2,param = -fitted_copula$par),flip=c(TRUE,FALSE)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    fam <- 23
  }
  else if (fitted_copula$family %in% c(33)) {#270 Clayton
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::claytonCopula(dim=2,param = -fitted_copula$par),flip=c(FALSE,TRUE)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    fam <- 33
  }
  else if (fitted_copula$family %in% c(3)){#Clayton
    result <- copula::gofCopula(copula=copula::claytonCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    fam <- 3
  }
  else if (fitted_copula$family %in% c(14)) {#180 Gumbel
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::gumbelCopula(dim=2,param = fitted_copula$par)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    fam <- 14
  }
  else if (fitted_copula$family %in% c(24)) {#90 Gumbel
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::gumbelCopula(dim=2,param = -fitted_copula$par),flip=c(TRUE,FALSE)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    fam <- 24
  }
  else if (fitted_copula$family %in% c(34)) {#270 Gumbel
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::gumbelCopula(dim=2,param = -fitted_copula$par),flip=c(FALSE,TRUE)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    fam <- 34
  }
  else if (fitted_copula$family %in% c(4)){#Gumbel
    result <- copula::gofCopula(copula=copula::gumbelCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    fam <- 4
  }
  else if (fitted_copula$family %in% c(1)){#Normal
    result <- copula::gofCopula(copula=copula::ellipCopula(family="normal",dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    fam <- 1
  }
  else if (fitted_copula$family %in% c(2)){#Student
    result <- copula::gofCopula(copula=copula::ellipCopula(family="t",dim=2,param=fitted_copula$par,df=round(fitted_copula$par2),df.fixed=TRUE), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    fam <- 2
  }
  else if (fitted_copula$family %in% c(5)){#Frank
    result <- copula::gofCopula(copula=copula::frankCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    fam <- 5
  }
  else if (fitted_copula$family %in% c(16)) {#180 Joe
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::joeCopula(dim=2,param = fitted_copula$par)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    fam <- 16
  }
  else if (fitted_copula$family %in% c(26)) {#90 Joe
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::joeCopula(dim=2,param = -fitted_copula$par),flip=c(TRUE,FALSE)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    fam <- 26
  }
  else if (fitted_copula$family %in% c(36)) {#270 Joe
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::joeCopula(dim=2,param = -fitted_copula$par),flip=c(FALSE,TRUE)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    fam <- 36
  }
  else if (fitted_copula$family %in% c(6)){#Joe
    result <- copula::gofCopula(copula=copula::joeCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    fam <- 6
  }
  #print(fam)
  result*abs(cor(u1, u2, method = "kendall", use = "complete.obs"))
}

test.unweighted$Rn.CvM.mult <- function(u1,u2,weights){
  fitted_copula <- suppressWarnings(VineCopula::BiCopSelect(u1=u1,u2=u2,familyset = c(1,2,3,4,5,6)))
  if (fitted_copula$family %in% c(13)) {#180 Clayton
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::claytonCopula(dim=2,param = fitted_copula$par)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    result2 <- copula::gofCopula(copula=copula::rotCopula(copula = copula::claytonCopula(dim=2,param = fitted_copula$par)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 13
  }
  else if (fitted_copula$family %in% c(23)) {#90 Clayton
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::claytonCopula(dim=2,param = -fitted_copula$par),flip=c(TRUE,FALSE)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    result2 <- copula::gofCopula(copula=copula::rotCopula(copula = copula::claytonCopula(dim=2,param = -fitted_copula$par),flip=c(TRUE,FALSE)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 23
  }
  else if (fitted_copula$family %in% c(33)) {#270 Clayton
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::claytonCopula(dim=2,param = -fitted_copula$par),flip=c(FALSE,TRUE)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    result2 <- copula::gofCopula(copula=copula::rotCopula(copula = copula::claytonCopula(dim=2,param = -fitted_copula$par),flip=c(FALSE,TRUE)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 33
  }
  else if (fitted_copula$family %in% c(3)){#Clayton
    result <- copula::gofCopula(copula=copula::claytonCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    result2 <- copula::gofCopula(copula=copula::claytonCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 3
  }
  else if (fitted_copula$family %in% c(14)) {#180 Gumbel
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::gumbelCopula(dim=2,param = fitted_copula$par)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    result2 <- copula::gofCopula(copula=copula::rotCopula(copula = copula::gumbelCopula(dim=2,param = fitted_copula$par)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 14
  }
  else if (fitted_copula$family %in% c(24)) {#90 Gumbel
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::gumbelCopula(dim=2,param = -fitted_copula$par),flip=c(TRUE,FALSE)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    result2 <- copula::gofCopula(copula=copula::rotCopula(copula = copula::gumbelCopula(dim=2,param = -fitted_copula$par),flip=c(TRUE,FALSE)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 24
  }
  else if (fitted_copula$family %in% c(34)) {#270 Gumbel
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::gumbelCopula(dim=2,param = -fitted_copula$par),flip=c(FALSE,TRUE)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    result2 <- copula::gofCopula(copula=copula::rotCopula(copula = copula::gumbelCopula(dim=2,param = -fitted_copula$par),flip=c(FALSE,TRUE)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 34
  }
  else if (fitted_copula$family %in% c(4)){#Gumbel
    result <- copula::gofCopula(copula=copula::gumbelCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    result2 <- copula::gofCopula(copula=copula::gumbelCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 4
  }
  else if (fitted_copula$family %in% c(1)){#Normal
    result <- copula::gofCopula(copula=copula::ellipCopula(family="normal",dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    result2 <- copula::gofCopula(copula=copula::ellipCopula(family="normal",dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 1
  }
  else if (fitted_copula$family %in% c(2)){#Student
    result <- copula::gofCopula(copula=copula::ellipCopula(family="t",dim=2,param=fitted_copula$par,df=round(fitted_copula$par2),df.fixed=TRUE), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    result2 <- copula::gofCopula(copula=copula::ellipCopula(family="t",dim=2,param=fitted_copula$par,df=round(fitted_copula$par2),df.fixed=TRUE), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 2
  }
  else if (fitted_copula$family %in% c(5)){#Frank
    result <- copula::gofCopula(copula=copula::frankCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    result2 <- copula::gofCopula(copula=copula::frankCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 5
  }
  else if (fitted_copula$family %in% c(16)) {#180 Joe
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::joeCopula(dim=2,param = fitted_copula$par)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    result2 <- copula::gofCopula(copula=copula::rotCopula(copula = copula::joeCopula(dim=2,param = fitted_copula$par)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 16
  }
  else if (fitted_copula$family %in% c(26)) {#90 Joe
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::joeCopula(dim=2,param = -fitted_copula$par),flip=c(TRUE,FALSE)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    result2 <- copula::gofCopula(copula=copula::rotCopula(copula = copula::joeCopula(dim=2,param = -fitted_copula$par),flip=c(TRUE,FALSE)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 26
  }
  else if (fitted_copula$family %in% c(36)) {#270 Joe
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::joeCopula(dim=2,param = -fitted_copula$par),flip=c(FALSE,TRUE)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    result2 <- copula::gofCopula(copula=copula::rotCopula(copula = copula::joeCopula(dim=2,param = -fitted_copula$par),flip=c(FALSE,TRUE)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 36
  }
  else if (fitted_copula$family %in% c(6)){#Joe
    result <- copula::gofCopula(copula=copula::joeCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    result2 <- copula::gofCopula(copula=copula::joeCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 6
  }
  #print(fam)
  #result=Rn, result2=CvM
  result*result2
}


# 3. Setup ----------------------------------------------------------------

set.seed(2021)
cores = 1 # set cores for parallelization of RVineStructureSelect().
# However, more efficient to parallelize the estimation of several
# RVines in parallel rather than estimate one model in parallel sequentially

cores2 <- detectCores()-1 # leave one core for other stuff

# get logreturns from above
logreturns <- daily_returns_stocks
logreturns <- logreturns%>%drop_na()

# remove last 31 periods to get an even number of periods
logreturns <- logreturns[1:(nrow(logreturns)-31),]

# calculate portfolio returns
logreturns$PFReturns <- NA
for(i in 1:length(logreturns$PFReturns)){
  logreturns$PFReturns[i] <- mean(as.numeric(logreturns[i,1:6]))
}

ncol_assets <- ncol(logreturns)-2 # the first n-2 number of columns are where
                                  #the asset returns are stored, last two are
                                  # date and portfolio returns
                                  # use as indicator for models later
logreturns <- logreturns%>%
  relocate(c(PFReturns,date), .after=last_col())





# 4. ARMA-GARCH Modelling ----------------------------------------------------

# ARMA-GARCH Specs for ugarch

# three possible ARMA-GARCH specs: normal, t distribution, skewed t distr.
specnorm = ugarchspec(variance.model = list(model = "fGARCH", submodel="GARCH", garchOrder=c(1,1)), distribution.model = "norm",  mean.model=list(armaOrder=c(1,1),include.mean=TRUE, arfima=FALSE,archm = FALSE, archpow = 1, arfima = FALSE, external.regressors = NULL, archex = FALSE))
specstd = ugarchspec(variance.model = list(model = "fGARCH", submodel="GARCH", garchOrder=c(1,1)), distribution.model = "std",  mean.model=list(armaOrder=c(1,1),include.mean=TRUE, arfima=FALSE,archm = FALSE, archpow = 1, arfima = FALSE, external.regressors = NULL, archex = FALSE))
specsstd = ugarchspec(variance.model = list(model = "fGARCH", submodel="GARCH", garchOrder=c(1,1)), distribution.model = "sstd",  mean.model=list(armaOrder=c(1,1),include.mean=TRUE, arfima=FALSE,archm = FALSE, archpow = 1, arfima = FALSE, external.regressors = NULL, archex = FALSE))

WE <- 900 # Window of estimation
N <- 10000 # Number of simulations per period for VaR

TComplete <- dim(logreturns)[1] # number of total periods

# ARMA-GARCH for the complete dataset to decide on distribution model

# Plot all assets
logreturns.plots <- list()


logreturns.plots$stoxx <- logreturns%>%dplyr::select(c(1,8))%>%
  ggplot(aes(x = date, y = `^STOXX`))+
  geom_line(size = 0.1)+
  theme_bw()+
  labs(y = "logreturns", title = "^STOXX")+
  theme(plot.title = element_text(hjust = 0.5))

logreturns.plots$ibex <- logreturns%>%dplyr::select(c(2,8))%>%
  ggplot(aes(x = date, y = `^IBEX`))+
  geom_line(size = 0.1)+
  theme_bw()+
  labs(y = "logreturns", title = "^IBEX")+
  theme(plot.title = element_text(hjust = 0.5))

logreturns.plots$ftse <- logreturns%>%dplyr::select(c(3,8))%>%
  ggplot(aes(x = date, y = `^FTSE`))+
  geom_line(size = 0.1)+
  theme_bw()+
  labs(y = "logreturns", title = "^FTSE")+
  theme(plot.title = element_text(hjust = 0.5))

logreturns.plots$gdaxi <- logreturns%>%dplyr::select(c(4,8))%>%
  ggplot(aes(x = date, y = `^GDAXI`))+
  geom_line(size = 0.1)+
  theme_bw()+
  labs(y = "logreturns", title = "^GDAXI")+
  theme(plot.title = element_text(hjust = 0.5))

logreturns.plots$clf <- logreturns%>%dplyr::select(c(5,8))%>%
  ggplot(aes(x = date, y = `CL=F`))+
  geom_line(size = 0.1)+
  theme_bw()+
  labs(y = "logreturns", title = "CL=F")+
  theme(plot.title = element_text(hjust = 0.5))

logreturns.plots$tnx <- logreturns%>%dplyr::select(c(6,8))%>%
  ggplot(aes(x = date, y = `^TNX`))+
  geom_line(size = 0.1)+
  theme_bw()+
  labs(y = "logreturns", title = "^TNX")+
  theme(plot.title = element_text(hjust = 0.5))

g <- grid.arrange(logreturns.plots$stoxx,logreturns.plots$ibex,logreturns.plots$ftse,
                    logreturns.plots$gdaxi,logreturns.plots$clf,logreturns.plots$tnx,
                  ncol = 2)
ggsave("logreturns.pdf", plot = g, path = "../MA/Plots/VaR", units = "cm", width = 20, height = 20)

#summary statistics

summaryfunction <- function(x){
  out <- list(mean(x),sd(x),base::min(x),base::max(x),kurtosis(x),skewness(x),jarque.bera.test(x)$p.value,
              Box.test(x, type = "Ljung-Box")$p.value,Box.test(abs(x), type = "Ljung-Box")$p.value)
}

summarystatistics <- logreturns%>%dplyr::select(c(1:6))%>%map(.,summaryfunction)%>%
  map(~ set_names(., c("mean","sd","min","max","kurtosis","skewness", "JQ", "LB","$LB_{abs}$")))%>%map_dfr(.,`[`)


summarystatistics <- t(summarystatistics)
colnames(summarystatistics) <- tickers
summarystatistics <- rbind(summarystatistics, cor(logreturns[1:6], method = "kendall"))
stargazer(summarystatistics, summary = F, align = F, digits = 5)

# ACF plots of return series

plotacf <- function(x, title){
  x = x
  ggplot(mapping = aes(y = acf(x)$acf, x = acf(x)$lag))+
    theme_bw()+
    geom_hline(aes(yintercept = 0)) +
    geom_segment(mapping = aes(xend = acf(x)$lag, yend = 0)) +
    geom_hline(aes(yintercept = qnorm((1 + 0.95)/2)/sqrt(acf(x)$n.used)), linetype = 2, color = 'darkblue') +
    geom_hline(aes(yintercept = -qnorm((1 + 0.95)/2)/sqrt(acf(x)$n.used)), linetype = 2, color = 'darkblue') +
    labs(x = "Lag", y = "ACF", title = title)+
    theme(plot.title = element_text(hjust = 0.5))
}
logreturns.acf.plots <- list()
logreturns.acf.plots$stoxx <- plotacf(logreturns$`^STOXX`,"^STOXX")
logreturns.acf.plots$ibex <- plotacf(logreturns$`^IBEX`,"^IBEX")
logreturns.acf.plots$ftse <- plotacf(logreturns$`^FTSE`,"^FTSE")
logreturns.acf.plots$gdaxi <- plotacf(logreturns$`^GDAXI`,"^GDAXI")
logreturns.acf.plots$clf <- plotacf(logreturns$`CL=F`,"CL=F")
logreturns.acf.plots$tnx <- plotacf(logreturns$`^TNX`,"^TNX")

g <- grid.arrange(logreturns.acf.plots$stoxx,logreturns.acf.plots$ibex,
                  logreturns.acf.plots$ftse,logreturns.acf.plots$gdaxi,
                  logreturns.acf.plots$clf,logreturns.acf.plots$tnx, ncol = 2)

ggsave("logreturns.acf.pdf", plot = g, path = "../MA/Plots/VaR", units = "cm", width = 15.5)

plotacf.abs <- function(x, title){
  x <- abs(x)
  ggplot(mapping = aes(y = acf(x)$acf, x = acf(x)$lag))+
    theme_bw()+
    geom_hline(aes(yintercept = 0)) +
    geom_segment(mapping = aes(xend = acf(x)$lag, yend = 0)) +
    geom_hline(aes(yintercept = qnorm((1 + 0.95)/2)/sqrt(acf(x)$n.used)), linetype = 2, color = 'darkblue') +
    geom_hline(aes(yintercept = -qnorm((1 + 0.95)/2)/sqrt(acf(x)$n.used)), linetype = 2, color = 'darkblue') +
    labs(x = "Lag", y = "ACF", title = title)+
    theme(plot.title = element_text(hjust = 0.5))
}

logreturns.acf.abs.plots <- list()
for(i in 1:6){
  logreturns.acf.abs.plots[[names(logreturns)[[i]]]] <- plotacf.abs(logreturns[,i],names(logreturns)[[i]])
}

g <- grid.arrange(logreturns.acf.abs.plots$`^STOXX`,logreturns.acf.abs.plots$`^IBEX`,
                  logreturns.acf.abs.plots$`^FTSE`,logreturns.acf.abs.plots$`^GDAXI`,
                  logreturns.acf.abs.plots$`CL=F`,logreturns.acf.abs.plots$`^TNX`, ncol = 2)
ggsave("logreturns.acf.abs.pdf", plot = g, path = "../MA/Plots/VaR", units = "cm", width = 15.5)

# histograms
logreturns.hist <- list()
logreturns.hist$stoxx <- logreturns%>%dplyr::select(1)%>%
  ggplot(data = ., aes(`^STOXX`))+
  theme_bw()+
  geom_histogram(aes(y=..density..),bins = 200, color = "grey50", fill ="lightblue", size = 0.1)+
  labs(x = "Returns", y = "Density", title = "^STOXX")+
  theme(plot.title = element_text(hjust = 0.5))

logreturns.hist$ibex <- logreturns%>%dplyr::select(2)%>%
  ggplot(data = ., aes(`^IBEX`))+
  theme_bw()+
  geom_histogram(aes(y=..density..),bins = 200, color = "grey50", fill ="lightblue", size = 0.1)+
  labs(x = "Returns", y = "Density", title = "^IBEX")+
  theme(plot.title = element_text(hjust = 0.5))

logreturns.hist$ftse <- logreturns%>%dplyr::select(3)%>%
  ggplot(data = ., aes(`^FTSE`))+
  theme_bw()+
  geom_histogram(aes(y=..density..),bins = 200, color = "grey50", fill ="lightblue", size = 0.1)+
  labs(x = "Returns", y = "Density", title = "^FTSE")+
  theme(plot.title = element_text(hjust = 0.5))

logreturns.hist$gdaxi <- logreturns%>%dplyr::select(4)%>%
  ggplot(data = ., aes(`^GDAXI`))+
  theme_bw()+
  geom_histogram(aes(y=..density..),bins = 200, color = "grey50", fill ="lightblue", size = 0.1)+
  labs(x = "Returns", y = "Density", title = "^GDAXI")+
  theme(plot.title = element_text(hjust = 0.5))

logreturns.hist$clf <- logreturns%>%dplyr::select(5)%>%
  ggplot(data = ., aes(`CL=F`))+
  theme_bw()+
  geom_histogram(aes(y=..density..),bins = 200, color = "grey50", fill ="lightblue", size = 0.1)+
  labs(x = "Returns", y = "Density", title = "CL=F")+
  theme(plot.title = element_text(hjust = 0.5))

logreturns.hist$tnx <- logreturns%>%dplyr::select(6)%>%
  ggplot(data = ., aes(`^TNX`))+
  theme_bw()+
  geom_histogram(aes(y=..density..),bins = 200, color = "grey50", fill ="lightblue", size = 0.1)+
  labs(x = "Returns", y = "Density", title = "^TNX")+
  theme(plot.title = element_text(hjust = 0.5))

g <- grid.arrange(logreturns.hist$stoxx,logreturns.hist$ibex,logreturns.hist$ftse,
                  logreturns.hist$gdaxi,logreturns.hist$clf,logreturns.hist$tnx,
                  ncol=2)
ggsave("histograms.pdf", plot = g, path = "../MA/Plots/VaR", units = "cm", width = 20, height = 20)


# fit arma-garch models with three different  distributions for all assets

armagarch <- list()
cl <- parallel::makeCluster(cores2)
doParallel::registerDoParallel(cl)
armagarch$norm <- foreach(i = 1:ncol_assets,.packages = c("rugarch"), .errorhandling = "pass",.final = function(x) setNames(x, names(logreturns[1:6]))) %dopar%{
    set.seed(2021)
    rugarch::ugarchfit(spec=specnorm, data=logreturns[(1:WE),i],solver.control=list(solver=c(9,10)),solver = "hybrid")#,tol = 1e-10
}
parallel::stopCluster(cl)

cl <- parallel::makeCluster(cores2)
doParallel::registerDoParallel(cl)
armagarch$std <- foreach(i = 1:ncol_assets,.packages = c("rugarch"), .errorhandling = "pass",.final = function(x) setNames(x, names(logreturns[1:6]))) %dopar%{
  set.seed(2021)
  rugarch::ugarchfit(spec=specstd, data=logreturns[(1:WE),i],solver.control=list(solver=c(9,10)),solver = "hybrid")#,tol = 1e-10
}
parallel::stopCluster(cl)

cl <- parallel::makeCluster(cores2)
doParallel::registerDoParallel(cl)
armagarch$sstd <- foreach(i = 1:ncol_assets,.packages = c("rugarch"), .errorhandling = "pass",.final = function(x) setNames(x, names(logreturns[1:6]))) %dopar%{
  set.seed(2021)
  rugarch::ugarchfit(spec=specsstd, data=logreturns[(1:WE),i],solver.control=list(solver=c(9,10)),solver = "hybrid")#,tol = 1e-10
}
parallel::stopCluster(cl)

armagarch.report <- function(fit){

  Ljung.Box <- Box.test(fit@fit$z, type = "Ljung-Box")$p.value
  Ljung.Box.stat <- Box.test(fit@fit$z, type = "Ljung-Box")$statistic
  Ljung.Box.abs <- Box.test(abs(fit@fit$z), type = "Ljung-Box")$p.value
  Ljung.Box.abs.stat <- Box.test(abs(fit@fit$z), type = "Ljung-Box")$statistic
  Arch <- ArchTest(fit@fit$z)$p.value
  Arch.stat <- ArchTest(fit@fit$z)$statistic
  results <- t(fit@fit$matcoef)
  results <- results[-(2:3),]
  results <- as.data.frame(results)
  results["LB"] <- NA
  results["LB.abs"] <- NA
  results["Arch"] <- NA
  results$LB[1] <- Ljung.Box.stat
  results$LB[2] <- Ljung.Box
  results$LB.abs[1] <- Ljung.Box.stat
  results$LB.abs[2] <- Ljung.Box.abs
  results$Arch[1] <- Arch.stat
  results$Arch[2] <- Arch

  return(round(results,2))

}

#norm
armagarch.results <- map_dfr(armagarch$norm, armagarch.report)
armagarch.results[" "] <- rep(c("Estimate/Stat.", "p-value"))
armagarch.results["asset"] <- rep(c(tickers), each = 2)
armagarch.results$asset[seq(2,12,2)] <- NA
armagarch.results <- armagarch.results%>%
  relocate(" ")%>%
  relocate("asset", .before = 2)

stargazer(as.matrix(armagarch.results), digits = 2, summary = F,
          covariate.labels=c(" ","ticker","$\\hat{\\mu}$", "$\\hat{\\phi}$", "$\\hat{\\theta}$",
                             "$\\hat{\\omega}$","$\\hat{\\alpha}$","$\\hat{\\beta}$",
                              "$LB$", "$LB_{abs}$", "$LM$"),
          column.sep.width = "1pt", align = TRUE, rownames = F, float.env = "sidewaystable",
          title ="ARMA(1,1)-GARCH(1,1) models with standard normal distribution")

#std
armagarch.results <- map_dfr(armagarch$std, armagarch.report)
armagarch.results[" "] <- rep(c("Estimate/Stat.", "p-value"))
armagarch.results["asset"] <- rep(c(tickers), each = 2)
armagarch.results$asset[seq(2,12,2)] <- NA
armagarch.results <- armagarch.results%>%
  relocate(" ")%>%
  relocate("asset", .before = 2)

stargazer(as.matrix(armagarch.results), digits = 2, summary = F,
          covariate.labels=c(" ","ticker","$\\hat{\\mu}$", "$\\hat{\\phi}$", "$\\hat{\\theta}$",
                             "$\\hat{\\omega}$","$\\hat{\\alpha}$","$\\hat{\\beta}$", "$shape$", "$LB$", "$LB_{abs}$", "$LM$"),
          column.sep.width = "1pt", align = T, rownames = F,float.env = "sidewaystable",
          title = "ARMA(1,1)-GARCH(1,1) models with standard student t distribution")

#sstd
armagarch.results <- map_dfr(armagarch$sstd, armagarch.report)
armagarch.results[" "] <- rep(c("Estimate/Stat.", "p-value"))
armagarch.results["asset"] <- rep(c(tickers), each = 2)
armagarch.results$asset[seq(2,12,2)] <- NA
armagarch.results <- armagarch.results%>%
  relocate(" ")%>%
  relocate("asset", .before = 2)

stargazer(as.matrix(armagarch.results), digits = 2, summary = F,
          covariate.labels=c(" ","ticker","$\\hat{\\mu}$", "$\\hat{\\phi}$", "$\\hat{\\theta}$",
                             "$\\hat{\\omega}$","$\\hat{\\alpha}$","$\\hat{\\beta}$",
                             "$skew$", "$shape$", "$LB$", "$LB_{abs}$", "$LM$"),
          column.sep.width = "1pt", align = TRUE, rownames = F,
          title = "ARMA(1,1)-GARCH(1,1) models with skewed standard student t distribution", float.env = "sidewaystable")

# acf of std. residuals
armagarch.acf.norm <- list()

for(i in 1:6){
  armagarch.acf.norm[[names(logreturns)[[i]]]] <- plotacf(armagarch$norm[[i]]@fit$z,names(armagarch$norm)[[i]])
}
g <- grid.arrange(armagarch.acf.norm$`^STOXX`,armagarch.acf.norm$`^IBEX`,
                  armagarch.acf.norm$`^FTSE`,armagarch.acf.norm$`^GDAXI`,
                  armagarch.acf.norm$`CL=F`,armagarch.acf.norm$`^TNX`, ncol = 2)
ggsave("armagarch.acf.norm.pdf", plot = g, path = "../MA/Plots/VaR", units = "cm", width = 15.5)

armagarch.acf.std <- list()

for(i in 1:6){
  armagarch.acf.std[[names(logreturns)[[i]]]] <- plotacf(armagarch$std[[i]]@fit$z,names(armagarch$std)[[i]])
}
g <- grid.arrange(armagarch.acf.std$`^STOXX`,armagarch.acf.std$`^IBEX`,
                  armagarch.acf.std$`^FTSE`,armagarch.acf.std$`^GDAXI`,
                  armagarch.acf.std$`CL=F`,armagarch.acf.std$`^TNX`, ncol = 2)
ggsave("armagarch.acf.std.pdf", plot = g, path = "../MA/Plots/VaR", units = "cm", width = 15.5)

armagarch.acf.sstd <- list()

for(i in 1:6){
  armagarch.acf.sstd[[names(logreturns)[[i]]]] <- plotacf(armagarch$sstd[[i]]@fit$z,names(armagarch$sstd)[[i]])
}
g <- grid.arrange(armagarch.acf.sstd$`^STOXX`,armagarch.acf.sstd$`^IBEX`,
                  armagarch.acf.sstd$`^FTSE`,armagarch.acf.sstd$`^GDAXI`,
                  armagarch.acf.sstd$`CL=F`,armagarch.acf.sstd$`^TNX`, ncol = 2)
ggsave("armagarch.acf.sstd.pdf", plot = g, path = "../MA/Plots/VaR", units = "cm", width = 15.5)


# acf of abs std. residuals
armagarch.acf.abs.norm <- list()

for(i in 1:6){
  armagarch.acf.abs.norm[[names(logreturns)[[i]]]] <- plotacf.abs(armagarch$norm[[i]]@fit$z,names(armagarch$norm)[[i]])
}
g <- grid.arrange(armagarch.acf.abs.norm$`^STOXX`,armagarch.acf.abs.norm$`^IBEX`,
             armagarch.acf.abs.norm$`^FTSE`,armagarch.acf.abs.norm$`^GDAXI`,
             armagarch.acf.abs.norm$`CL=F`,armagarch.acf.abs.norm$`^TNX`, ncol = 2)
ggsave("armagarch.acf.abs.norm.pdf", plot = g, path = "../MA/Plots/VaR", units = "cm", width = 15.5)


armagarch.acf.abs.std <- list()
for(i in 1:6){
  armagarch.acf.abs.std[[names(logreturns)[[i]]]] <- plotacf.abs(armagarch$std[[i]]@fit$z,names(armagarch$std)[[i]])
}
g <- grid.arrange(armagarch.acf.abs.std$`^STOXX`,armagarch.acf.abs.std$`^IBEX`,
                  armagarch.acf.abs.std$`^FTSE`,armagarch.acf.abs.std$`^GDAXI`,
                  armagarch.acf.abs.std$`CL=F`,armagarch.acf.abs.std$`^TNX`, ncol = 2)
ggsave("armagarch.acf.abs.std.pdf", plot = g, path = "../MA/Plots/VaR", units = "cm", width = 15.5)

armagarch.acf.abs.sstd <- list()
for(i in 1:6){
  armagarch.acf.abs.sstd[[names(logreturns)[[i]]]] <- plotacf.abs(armagarch$sstd[[i]]@fit$z,names(armagarch$sstd)[[i]])
}
g <- grid.arrange(armagarch.acf.abs.sstd$`^STOXX`,armagarch.acf.abs.sstd$`^IBEX`,
                  armagarch.acf.abs.sstd$`^FTSE`,armagarch.acf.abs.sstd$`^GDAXI`,
                  armagarch.acf.abs.sstd$`CL=F`,armagarch.acf.abs.sstd$`^TNX`, ncol = 2)
ggsave("armagarch.acf.abs.sstd.pdf", plot = g, path = "../MA/Plots/VaR", units = "cm", width = 15.5)

#qq plots


qqplotnorm <- function(x, title){
  x = x
  ggplot(data = NULL, aes(sample = x@fit$z))+
    theme_bw()+
    geom_qq(distribution = stats::qnorm, color = "dodgerblue3", size = 1.2, fill = "transparent", shape = 21,
            alpha = 1)+
    geom_qq_line()+
    theme_bw()+
    labs(x = "Theoretical Quantiles", y = "Sample Quantiles", title = title)+
    theme(plot.title = element_text(hjust = 0.5))
}

qqplotstd <- function(x, title){
  x = x
  ggplot(data = NULL, aes(sample = x@fit$z))+
    geom_qq(distribution = qt, dparams = c(df=as.numeric(x@fit$coef[7])),
            color = "dodgerblue3", size = 1.2, fill = "transparent", shape = 21,
            alpha = 1)+
    geom_qq_line(distribution = qt, dparams = c(df=as.numeric(x@fit$coef[7])))+
    theme_bw()+
    labs(x = "Theoretical Quantiles", y = "Sample Quantiles", title = title)+
    theme(plot.title = element_text(hjust = 0.5))
}

qqplotsstd <- function(x, title){
  x = x
  ggplot(data = NULL, aes(sample = x@fit$z))+
    geom_qq(distribution = rugarch::qdist, dparams = c(distribution = "sstd", shape = as.numeric(x@fit$coef[8]),
                                                       skew = as.numeric(x@fit$coef[7])),
            color = "dodgerblue3", size = 1.2, fill = "transparent", shape = 21,
            alpha = 1)+
    geom_qq_line(distribution = rugarch::qdist, dparams = c(distribution = "sstd", shape = as.numeric(x@fit$coef[8]),
                                                            skew = as.numeric(x@fit$coef[7])))+
    theme_bw()+
    labs(x = "Theoretical Quantiles", y = "Sample Quantiles", title = title)+
    theme(plot.title = element_text(hjust = 0.5))
}


qqplots <- list()
qqplots$norm <- list()
for(i in 1:6){
  qqplots$norm[[names(logreturns)[[i]]]] <- qqplotnorm(armagarch$norm[[i]],names(armagarch$norm)[[i]])
}

g <- grid.arrange(qqplots$norm$`^STOXX`, qqplots$norm$`^IBEX`, qqplots$norm$`^FTSE`,
                  qqplots$norm$`^GDAXI`, qqplots$norm$`CL=F`, qqplots$norm$`^TNX`, ncol = 3)
ggsave("qqplots.norm.pdf", plot = g, path = "../MA/Plots/VaR", units = "cm", width = 30)

qqplots$std <- list()
for(i in 1:6){
  qqplots$std[[names(logreturns)[[i]]]] <- qqplotstd(armagarch$std[[i]],names(armagarch$std)[[i]])
}
g <- grid.arrange(qqplots$std$`^STOXX`, qqplots$std$`^IBEX`, qqplots$std$`^FTSE`,
                  qqplots$std$`^GDAXI`, qqplots$std$`CL=F`, qqplots$std$`^TNX`, ncol = 3)
ggsave("qqplots.std.pdf", plot = g, path = "../MA/Plots/VaR", units = "cm", width = 30)

qqplots$sstd <- list()
for(i in 1:6){
  qqplots$sstd[[names(logreturns)[[i]]]] <- qqplotsstd(armagarch$sstd[[i]],names(armagarch$sstd)[[i]])
}
g <- grid.arrange(qqplots$sstd$`^STOXX`, qqplots$sstd$`^IBEX`, qqplots$sstd$`^FTSE`,
                  qqplots$sstd$`^GDAXI`, qqplots$sstd$`CL=F`, qqplots$sstd$`^TNX`, ncol = 3)
ggsave("qqplots.sstd.pdf", plot = g, path = "../MA/Plots/VaR", units = "cm", width = 30)

# Settle on ARMA(1,1)-GARCH(1,1) with sstd distribution
# estimate for all time windows for all assets with moving window approach
cl <- parallel::makeCluster(cores2)
doParallel::registerDoParallel(cl)
GARCH.sstd <- foreach(i = 1:ncol_assets,.packages = c("rugarch"), .errorhandling = "pass") %:% #ncol_assets
  foreach(j = 1:(TComplete-WE)) %dopar%{# (T-WE)
    set.seed(2021)
    rugarch::ugarchfit(spec=specsstd, data=logreturns[(j):(j+WE-1),i],solver.control=list(solver=c(9,10)),solver = "hybrid")#,tol = 1e-10
  }
parallel::stopCluster(cl)


#check if convergence error in a model
convergence.sstd <- matrix(ncol=ncol_assets, nrow=(TComplete-WE))
for(i in 1:ncol_assets){
  for(j in 1:(TComplete-WE)){
    convergence.sstd[j,i] = GARCH.sstd[[i]][[j]]@fit$convergence
  }
}
any(convergence.sstd!=0) # If FALSE, then no convergence error

GARCH <- GARCH.sstd


# models for selection of weighing parameter (last 10 in-sample periods)
cl <- parallel::makeCluster(cores2)
doParallel::registerDoParallel(cl)
GARCH.for.weights <- foreach(i = 1:ncol_assets,.packages = c("rugarch"), .errorhandling = "pass") %:% #ncol_assets
  foreach(j = 1:(10)) %dopar%{# (T-WE)
    set.seed(2021)
    rugarch::ugarchfit(spec=specsstd, data=logreturns[(j):(j+WE-10),i],solver.control=list(solver=c(9,10)),solver = "hybrid")#,tol = 1e-10
  }
parallel::stopCluster(cl)

# 5. Fit Vines for first Period with unweighted Tests ---------------------

z_first <- matrix(nrow=WE, ncol = ncol_assets)
for(i in 1:ncol_assets){
  if("skew" %in% names(GARCH[[i]][[1]]@fit$coef)){
    z_first[,i] = psstd(GARCH[[i]][[1]]@fit$z, nu=GARCH[[i]][[1]]@fit$coef[8],
                        xi = GARCH[[i]][[1]]@fit$coef[7])
  }else{
    z_first[,i] = pstd(GARCH[[i]][[1]]@fit$z, nu=GARCH[[i]][[1]]@fit$coef[7])}
}


vine.first.unweighted <- list(length = length(test.unweighted))
cl <- parallel::makeCluster(cores2)
doParallel::registerDoParallel(cl)
vine.first.unweighted <- foreach(i = 1:length(test.unweighted),.final = function(x) setNames(x, names(test.unweighted))) %dopar% {
  set.seed(2021)
  VineCopula::RVineStructureSelect(data=z_first,familyset = c(1,2,3,4,5,6), treecrit = test.unweighted[[i]])
}
parallel::stopCluster(cl)

vine.first.unweighted.aics <- matrix(nrow = 1, ncol = length(test.unweighted))
colnames(vine.first.unweighted.aics) <- names(test.unweighted)
for(i in seq_along(vine.first.unweighted.aics)){
  vine.first.unweighted.aics[,names(test.unweighted)[i]] <- vine.first.unweighted[[i]]$AIC
}


#6 . Find best weighting Parameters for weighted tests --------------------

# weighing paramter
par <- seq(from=0,to=1,by=0.1)

# list of tests that have weighting parameter that needs to be determined
find.weights <- list()

find.weights$CvM.tau <- list(length=length(par))
for(i in 1:length(par)){
  find.weights$CvM.tau[[i]] <- local({
    mu <- par[i]
    function(u1,u2,weights){
      fitted_copula <- suppressWarnings(VineCopula::BiCopSelect(u1=u1,u2=u2,familyset = c(1,2,3,4,5,6)))
      if (fitted_copula$family %in% c(13)) {#180 Clayton
        result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::claytonCopula(dim=2,param = fitted_copula$par)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
        fam <- 13
      }
      else if (fitted_copula$family %in% c(23)) {#90 Clayton
        result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::claytonCopula(dim=2,param = -fitted_copula$par),flip=c(TRUE,FALSE)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
        fam <- 23
      }
      else if (fitted_copula$family %in% c(33)) {#270 Clayton
        result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::claytonCopula(dim=2,param = -fitted_copula$par),flip=c(FALSE,TRUE)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
        fam <- 33
      }
      else if (fitted_copula$family %in% c(3)){#Clayton
        result <- copula::gofCopula(copula=copula::claytonCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
        fam <- 3
      }
      else if (fitted_copula$family %in% c(14)) {#180 Gumbel
        result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::gumbelCopula(dim=2,param = fitted_copula$par)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
        fam <- 14
      }
      else if (fitted_copula$family %in% c(24)) {#90 Gumbel
        result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::gumbelCopula(dim=2,param = -fitted_copula$par),flip=c(TRUE,FALSE)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
        fam <- 24
      }
      else if (fitted_copula$family %in% c(34)) {#270 Gumbel
        result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::gumbelCopula(dim=2,param = -fitted_copula$par),flip=c(FALSE,TRUE)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
        fam <- 34
      }
      else if (fitted_copula$family %in% c(4)){#Gumbel
        result <- copula::gofCopula(copula=copula::gumbelCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
        fam <- 4
      }
      else if (fitted_copula$family %in% c(1)){#Normal
        result <- copula::gofCopula(copula=copula::ellipCopula(family="normal",dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
        fam <- 1
      }
      else if (fitted_copula$family %in% c(2)){#Student
        result <- copula::gofCopula(copula=copula::ellipCopula(family="t",dim=2,param=fitted_copula$par,df=round(fitted_copula$par2),df.fixed=TRUE), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
        fam <- 2
      }
      else if (fitted_copula$family %in% c(5)){#Frank
        result <- copula::gofCopula(copula=copula::frankCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
        fam <- 5
      }
      else if (fitted_copula$family %in% c(16)) {#180 Joe
        result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::joeCopula(dim=2,param = fitted_copula$par)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
        fam <- 16
      }
      else if (fitted_copula$family %in% c(26)) {#90 Joe
        result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::joeCopula(dim=2,param = -fitted_copula$par),flip=c(TRUE,FALSE)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
        fam <- 26
      }
      else if (fitted_copula$family %in% c(36)) {#270 Joe
        result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::joeCopula(dim=2,param = -fitted_copula$par),flip=c(FALSE,TRUE)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
        fam <- 36
      }
      else if (fitted_copula$family %in% c(6)){#Joe
        result <- copula::gofCopula(copula=copula::joeCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="pb")$p.value
        fam <- 6
      }
      #print(fam)
      mu*result+(1-mu)*abs(cor(u1, u2, method = "kendall", use = "complete.obs"))
    }
  })
}

find.weights$Rn.tau <- list(length=length(par))
for(i in 1:length(par)){
  find.weights$Rn.tau[[i]] <- local({
    mu <- par[i]
    function(u1,u2,weights){
      fitted_copula <- suppressWarnings(VineCopula::BiCopSelect(u1=u1,u2=u2,familyset = c(1,2,3,4,5,6)))
      if (fitted_copula$family %in% c(13)) {#180 Clayton
        result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::claytonCopula(dim=2,param = fitted_copula$par)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
        fam <- 13
      }
      else if (fitted_copula$family %in% c(23)) {#90 Clayton
        result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::claytonCopula(dim=2,param = -fitted_copula$par),flip=c(TRUE,FALSE)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
        fam <- 23
      }
      else if (fitted_copula$family %in% c(33)) {#270 Clayton
        result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::claytonCopula(dim=2,param = -fitted_copula$par),flip=c(FALSE,TRUE)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
        fam <- 33
      }
      else if (fitted_copula$family %in% c(3)){#Clayton
        result <- copula::gofCopula(copula=copula::claytonCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
        fam <- 3
      }
      else if (fitted_copula$family %in% c(14)) {#180 Gumbel
        result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::gumbelCopula(dim=2,param = fitted_copula$par)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
        fam <- 14
      }
      else if (fitted_copula$family %in% c(24)) {#90 Gumbel
        result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::gumbelCopula(dim=2,param = -fitted_copula$par),flip=c(TRUE,FALSE)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
        fam <- 24
      }
      else if (fitted_copula$family %in% c(34)) {#270 Gumbel
        result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::gumbelCopula(dim=2,param = -fitted_copula$par),flip=c(FALSE,TRUE)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
        fam <- 34
      }
      else if (fitted_copula$family %in% c(4)){#Gumbel
        result <- copula::gofCopula(copula=copula::gumbelCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
        fam <- 4
      }
      else if (fitted_copula$family %in% c(1)){#Normal
        result <- copula::gofCopula(copula=copula::ellipCopula(family="normal",dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
        fam <- 1
      }
      else if (fitted_copula$family %in% c(2)){#Student
        result <- copula::gofCopula(copula=copula::ellipCopula(family="t",dim=2,param=fitted_copula$par,df=round(fitted_copula$par2),df.fixed=TRUE), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
        fam <- 2
      }
      else if (fitted_copula$family %in% c(5)){#Frank
        result <- copula::gofCopula(copula=copula::frankCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
        fam <- 5
      }
      else if (fitted_copula$family %in% c(16)) {#180 Joe
        result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::joeCopula(dim=2,param = fitted_copula$par)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
        fam <- 16
      }
      else if (fitted_copula$family %in% c(26)) {#90 Joe
        result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::joeCopula(dim=2,param = -fitted_copula$par),flip=c(TRUE,FALSE)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
        fam <- 26
      }
      else if (fitted_copula$family %in% c(36)) {#270 Joe
        result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::joeCopula(dim=2,param = -fitted_copula$par),flip=c(FALSE,TRUE)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
        fam <- 36
      }
      else if (fitted_copula$family %in% c(6)){#Joe
        result <- copula::gofCopula(copula=copula::joeCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
        fam <- 6
      }
      #print(fam)
      mu*result+(1-mu)*abs(cor(u1, u2, method = "kendall", use = "complete.obs"))
    }
  })
}

find.weights$Rn.CvM <- list(length=length(par))
for(i in 1:length(par)){
  find.weights$Rn.CvM[[i]] <- local({
    mu <- par[i]
    function(u1,u2,weights){
      fitted_copula <- suppressWarnings(VineCopula::BiCopSelect(u1=u1,u2=u2,familyset = c(1,2,3,4,5,6)))
      if (fitted_copula$family %in% c(13)) {#180 Clayton
        result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::claytonCopula(dim=2,param = fitted_copula$par)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
        result2 <- copula::gofCopula(copula=copula::rotCopula(copula = copula::claytonCopula(dim=2,param = fitted_copula$par)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
        fam <- 13
      }
      else if (fitted_copula$family %in% c(23)) {#90 Clayton
        result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::claytonCopula(dim=2,param = -fitted_copula$par),flip=c(TRUE,FALSE)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
        result2 <- copula::gofCopula(copula=copula::rotCopula(copula = copula::claytonCopula(dim=2,param = -fitted_copula$par),flip=c(TRUE,FALSE)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
        fam <- 23
      }
      else if (fitted_copula$family %in% c(33)) {#270 Clayton
        result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::claytonCopula(dim=2,param = -fitted_copula$par),flip=c(FALSE,TRUE)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
        result2 <- copula::gofCopula(copula=copula::rotCopula(copula = copula::claytonCopula(dim=2,param = -fitted_copula$par),flip=c(FALSE,TRUE)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
        fam <- 33
      }
      else if (fitted_copula$family %in% c(3)){#Clayton
        result <- copula::gofCopula(copula=copula::claytonCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
        result2 <- copula::gofCopula(copula=copula::claytonCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
        fam <- 3
      }
      else if (fitted_copula$family %in% c(14)) {#180 Gumbel
        result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::gumbelCopula(dim=2,param = fitted_copula$par)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
        result2 <- copula::gofCopula(copula=copula::rotCopula(copula = copula::gumbelCopula(dim=2,param = fitted_copula$par)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
        fam <- 14
      }
      else if (fitted_copula$family %in% c(24)) {#90 Gumbel
        result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::gumbelCopula(dim=2,param = -fitted_copula$par),flip=c(TRUE,FALSE)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
        result2 <- copula::gofCopula(copula=copula::rotCopula(copula = copula::gumbelCopula(dim=2,param = -fitted_copula$par),flip=c(TRUE,FALSE)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
        fam <- 24
      }
      else if (fitted_copula$family %in% c(34)) {#270 Gumbel
        result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::gumbelCopula(dim=2,param = -fitted_copula$par),flip=c(FALSE,TRUE)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
        result2 <- copula::gofCopula(copula=copula::rotCopula(copula = copula::gumbelCopula(dim=2,param = -fitted_copula$par),flip=c(FALSE,TRUE)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
        fam <- 34
      }
      else if (fitted_copula$family %in% c(4)){#Gumbel
        result <- copula::gofCopula(copula=copula::gumbelCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
        result2 <- copula::gofCopula(copula=copula::gumbelCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
        fam <- 4
      }
      else if (fitted_copula$family %in% c(1)){#Normal
        result <- copula::gofCopula(copula=copula::ellipCopula(family="normal",dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
        result2 <- copula::gofCopula(copula=copula::ellipCopula(family="normal",dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
        fam <- 1
      }
      else if (fitted_copula$family %in% c(2)){#Student
        result <- copula::gofCopula(copula=copula::ellipCopula(family="t",dim=2,param=fitted_copula$par,df=round(fitted_copula$par2),df.fixed=TRUE), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
        result2 <- copula::gofCopula(copula=copula::ellipCopula(family="t",dim=2,param=fitted_copula$par,df=round(fitted_copula$par2),df.fixed=TRUE), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
        fam <- 2
      }
      else if (fitted_copula$family %in% c(5)){#Frank
        result <- copula::gofCopula(copula=copula::frankCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
        result2 <- copula::gofCopula(copula=copula::frankCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
        fam <- 5
      }
      else if (fitted_copula$family %in% c(16)) {#180 Joe
        result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::joeCopula(dim=2,param = fitted_copula$par)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
        result2 <- copula::gofCopula(copula=copula::rotCopula(copula = copula::joeCopula(dim=2,param = fitted_copula$par)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
        fam <- 16
      }
      else if (fitted_copula$family %in% c(26)) {#90 Joe
        result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::joeCopula(dim=2,param = -fitted_copula$par),flip=c(TRUE,FALSE)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
        result2 <- copula::gofCopula(copula=copula::rotCopula(copula = copula::joeCopula(dim=2,param = -fitted_copula$par),flip=c(TRUE,FALSE)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
        fam <- 26
      }
      else if (fitted_copula$family %in% c(36)) {#270 Joe
        result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::joeCopula(dim=2,param = -fitted_copula$par),flip=c(FALSE,TRUE)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
        result2 <- copula::gofCopula(copula=copula::rotCopula(copula = copula::joeCopula(dim=2,param = -fitted_copula$par),flip=c(FALSE,TRUE)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
        fam <- 36
      }
      else if (fitted_copula$family %in% c(6)){#Joe
        result <- copula::gofCopula(copula=copula::joeCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
        result2 <- copula::gofCopula(copula=copula::joeCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
        fam <- 6
      }
      #print(fam)
      #result=Rn, result2=CvM
      mu*result + (1-mu)*result2
    }
  })
}


# for each test, fit in parallel one vine for each possible weighting paramter


find.weights.results <- list()
cl <- parallel::makeCluster(cores2)
doParallel::registerDoParallel(cl)
find.weights.results <- foreach(i = 1:length(find.weights),.final = function(x) setNames(x, names(find.weights))) %:%
  foreach(k = 1:10)%:%
  foreach(j = 1:(length(par)))%dopar%{
    set.seed(2021)
    zvalues <- matrix(nrow=WE-9, ncol=ncol_assets)#nrow_assets
    for(l in 1:ncol_assets){#ncol_assets
      if("skew" %in% names(GARCH[[l]][[k]]@fit$coef)){
        zvalues[,l] = fGarch::psstd(GARCH.for.weights[[l]][[k]]@fit$z, nu=GARCH.for.weights[[l]][[k]]@fit$coef[8],
                                    xi = GARCH.for.weights[[l]][[k]]@fit$coef[7])
      }else{
         zvalues[,l] = fGarch::pstd(GARCH.for.weights[[j]][[i]]@fit$z, nu=GARCH.for.weights[[j]][[k]]@fit$coef[7])}
    }
    VineCopula::RVineStructureSelect(data=zvalues,familyset = c(1,2,3,4,5,6), treecrit = find.weights[[names(find.weights)[i]]][[j]])
  }
parallel::stopCluster(cl)



find.weights.aics <- list()
find.weights.plots <- list()
#CvM.tau
find.weights.aics$CvM.tau <- matrix(nrow = length(par), ncol = length(find.weights.results$CvM.tau))

for(i in 1:length(par)){
  for(j in 1:length(find.weights.results$CvM.tau)){
    find.weights.aics$CvM.tau[i,j] = find.weights.results$CvM.tau[[j]][[i]]$AIC
  }
}
colnames(find.weights.aics$CvM.tau) <- seq(ncol(find.weights.aics$CvM.tau))
find.weights.aics$CvM.tau <- as.data.frame(find.weights.aics$CvM.tau)
find.weights.aics$CvM.tau$parameter <- par


min <- find.weights.aics$CvM.tau%>%
  pivot_longer(cols = 1:10, names_to = "period")%>%
  arrange(as.numeric(period))%>%
  group_by(period)%>%
  dplyr::filter(value == min(value))

find.weights.plots$CvM.tau <- find.weights.aics$CvM.tau%>%
  pivot_longer(cols = 1:10, names_to = "period")%>%
  arrange(as.numeric(period))%>%
  ggplot(data = ., aes(x = parameter, y = value, colour = period))+
  scale_x_continuous(breaks = seq(0, 1, by = 0.1))+
  geom_line()+
  geom_point(data = min, aes(x=parameter, y = value, group = period), colour = "black", alpha = 1, shape = 2)+
  theme_bw()+
  theme(legend.position = "bottom")+
  labs(y = "AIC", title = TeX("$\\eta_{S_n,\\tau,w}p_{S_n}+(1-\\eta_{S_n,\\tau,w})\\tau$"), x = TeX("$\\eta_{S_n,\\tau,w}$"))+
  theme(plot.title = element_text(hjust = 0.5))+
  viridis::scale_color_viridis(discrete = T)+
  theme(legend.position = "none")

find.weights.plots$CvM.tau.legend <- find.weights.aics$CvM.tau%>%
  pivot_longer(cols = 1:10, names_to = "period")%>%
  arrange(as.numeric(period))%>%
  ggplot(data = ., aes(x = parameter, y = value, colour = period))+
  scale_x_continuous(breaks = seq(0, 1, by = 0.1))+
  geom_line()+
  geom_point(data = min, aes(x=parameter, y = value, group = period), colour = "black", alpha = 1, shape = 2)+
  theme_bw()+
  theme(legend.position = "bottom")+
  labs(y = "AIC", title = TeX("$\\eta_{S_n,\\tau,w}p_{S_n}+(1-\\eta_{S_n,\\tau,w})\\tau$"), x = TeX("$\\eta_{S_n,\\tau,w}$"))+
  theme(plot.title = element_text(hjust = 0.5))+
  viridis::scale_color_viridis(discrete = T)

#Rn.tau
find.weights.aics$Rn.tau <- matrix(nrow = length(par), ncol = length(find.weights.results$Rn.tau))

for(i in 1:length(par)){
  for(j in 1:length(find.weights.results$Rn.tau)){
    find.weights.aics$Rn.tau[i,j] = find.weights.results$Rn.tau[[j]][[i]]$AIC
  }
}
colnames(find.weights.aics$Rn.tau) <- seq(ncol(find.weights.aics$Rn.tau))
find.weights.aics$Rn.tau <- as.data.frame(find.weights.aics$Rn.tau)
find.weights.aics$Rn.tau$parameter <- par


min <- find.weights.aics$Rn.tau%>%
  pivot_longer(cols = 1:10, names_to = "period")%>%
  arrange(as.numeric(period))%>%
  group_by(period)%>%
  dplyr::filter(value == min(value))

find.weights.plots$Rn.tau <- find.weights.aics$Rn.tau%>%
  pivot_longer(cols = 1:10, names_to = "period")%>%
  arrange(as.numeric(period))%>%
  ggplot(data = ., aes(x = parameter, y = value, colour = period))+
  scale_x_continuous(breaks = seq(0, 1, by = 0.1))+
  geom_line()+
  geom_point(data = min, aes(x=parameter, y = value, group = period), colour = "black", alpha = 1, shape = 2)+
  theme_bw()+
  theme(legend.position = "bottom")+
  labs(y = "AIC", title = TeX("$\\eta_{R_n,\\tau,w}p_{R_n}+(1-\\eta_{R_n,\\tau,w})\\tau$"), x = TeX("$\\eta_{R_n,\\tau,w}$"))+
  theme(plot.title = element_text(hjust = 0.5))+
  viridis::scale_color_viridis(discrete = T)+
  theme(legend.position = "none")


#Rn.CvM
find.weights.aics$Rn.CvM <- matrix(nrow = length(par), ncol = length(find.weights.results$Rn.CvM))

for(i in 1:length(par)){
  for(j in 1:length(find.weights.results$Rn.CvM)){
    find.weights.aics$Rn.CvM[i,j] = find.weights.results$Rn.CvM[[j]][[i]]$AIC
  }
}
colnames(find.weights.aics$Rn.CvM) <- seq(ncol(find.weights.aics$Rn.CvM))
find.weights.aics$Rn.CvM <- as.data.frame(find.weights.aics$Rn.CvM)
find.weights.aics$Rn.CvM$parameter <- par


min <- find.weights.aics$Rn.CvM%>%
  pivot_longer(cols = 1:10, names_to = "period")%>%
  arrange(as.numeric(period))%>%
  group_by(period)%>%
  dplyr::filter(value == min(value))

find.weights.plots$Rn.CvM <- find.weights.aics$Rn.CvM%>%
  pivot_longer(cols = 1:10, names_to = "period")%>%
  arrange(as.numeric(period))%>%
  ggplot(data = ., aes(x = parameter, y = value, colour = period))+
  scale_x_continuous(breaks = seq(0, 1, by = 0.1))+
  geom_line()+
  geom_point(data = min, aes(x=parameter, y = value, group = period), colour = "black", alpha = 1, shape = 2)+
  theme_bw()+
  labs(y = "AIC", title = TeX("$\\eta_{R_n,S_n,w}p_{R_n}+(1-\\eta_{R_n,S_n,w})\\p_{S_n}$"), x = TeX("$\\eta_{R_n,\\S_n,w}$"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(aspect.ratio = 0.8)+
  theme(legend.position = "none")





get_only_legend <- function(plot) {

  # get tabular interpretation of plot
  plot_table <- ggplot_gtable(ggplot_build(plot))

  #  Mark only legend in plot
  legend_plot <- which(sapply(plot_table$grobs, function(x) x$name) == "guide-box")

  # extract legend
  legend <- plot_table$grobs[[legend_plot]]

  # return legend
  return(legend)
}
find.weights.plots$legend <- get_only_legend(find.weights.plots$CvM.tau.legend)
combined_plot <- grid.arrange(find.weights.plots$CvM.tau, find.weights.plots$Rn.tau, ncol = 2, nrow = 1)
g <- grid.arrange(combined_plot, find.weights.plots$legend, nrow = 2, heights = c(10,1))
ggsave("find.weights.pdf", plot = g, path = "../MA/Plots/VaR", units = "cm", width = 20, height = 11)

best.weights <- list()
best.weights$Rn.tau <- 0.1
best.weights$CvM.tau <- 0.1
best.weights$Rn.CvM <- 0.1

test.weighted.best <- list()
test.weighted.best$CvM.tau <- function(u1,u2,weights){
  fitted_copula <- suppressWarnings(VineCopula::BiCopSelect(u1=u1,u2=u2,familyset = c(1,2,3,4,5,6)))
  if (fitted_copula$family %in% c(13)) {#180 Clayton
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::claytonCopula(dim=2,param = fitted_copula$par)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 13
  }
  else if (fitted_copula$family %in% c(23)) {#90 Clayton
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::claytonCopula(dim=2,param = -fitted_copula$par),flip=c(TRUE,FALSE)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 23
  }
  else if (fitted_copula$family %in% c(33)) {#270 Clayton
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::claytonCopula(dim=2,param = -fitted_copula$par),flip=c(FALSE,TRUE)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 33
  }
  else if (fitted_copula$family %in% c(3)){#Clayton
    result <- copula::gofCopula(copula=copula::claytonCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 3
  }
  else if (fitted_copula$family %in% c(14)) {#180 Gumbel
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::gumbelCopula(dim=2,param = fitted_copula$par)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 14
  }
  else if (fitted_copula$family %in% c(24)) {#90 Gumbel
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::gumbelCopula(dim=2,param = -fitted_copula$par),flip=c(TRUE,FALSE)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 24
  }
  else if (fitted_copula$family %in% c(34)) {#270 Gumbel
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::gumbelCopula(dim=2,param = -fitted_copula$par),flip=c(FALSE,TRUE)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 34
  }
  else if (fitted_copula$family %in% c(4)){#Gumbel
    result <- copula::gofCopula(copula=copula::gumbelCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 4
  }
  else if (fitted_copula$family %in% c(1)){#Normal
    result <- copula::gofCopula(copula=copula::ellipCopula(family="normal",dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 1
  }
  else if (fitted_copula$family %in% c(2)){#Student
    result <- copula::gofCopula(copula=copula::ellipCopula(family="t",dim=2,param=fitted_copula$par,df=round(fitted_copula$par2),df.fixed=TRUE), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 2
  }
  else if (fitted_copula$family %in% c(5)){#Frank
    result <- copula::gofCopula(copula=copula::frankCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 5
  }
  else if (fitted_copula$family %in% c(16)) {#180 Joe
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::joeCopula(dim=2,param = fitted_copula$par)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 16
  }
  else if (fitted_copula$family %in% c(26)) {#90 Joe
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::joeCopula(dim=2,param = -fitted_copula$par),flip=c(TRUE,FALSE)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 26
  }
  else if (fitted_copula$family %in% c(36)) {#270 Joe
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::joeCopula(dim=2,param = -fitted_copula$par),flip=c(FALSE,TRUE)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 36
  }
  else if (fitted_copula$family %in% c(6)){#Joe
    result <- copula::gofCopula(copula=copula::joeCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="pb")$p.value
    fam <- 6
  }
  #print(fam)
  best.weights$CvM.tau*result+(1-best.weights$CvM.tau)*abs(cor(u1, u2, method = "kendall", use = "complete.obs"))
}
test.weighted.best$Rn.tau <- function(u1,u2,weights){
  fitted_copula <- suppressWarnings(VineCopula::BiCopSelect(u1=u1,u2=u2,familyset = c(1,2,3,4,5,6)))
  if (fitted_copula$family %in% c(13)) {#180 Clayton
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::claytonCopula(dim=2,param = fitted_copula$par)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    fam <- 13
  }
  else if (fitted_copula$family %in% c(23)) {#90 Clayton
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::claytonCopula(dim=2,param = -fitted_copula$par),flip=c(TRUE,FALSE)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    fam <- 23
  }
  else if (fitted_copula$family %in% c(33)) {#270 Clayton
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::claytonCopula(dim=2,param = -fitted_copula$par),flip=c(FALSE,TRUE)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    fam <- 33
  }
  else if (fitted_copula$family %in% c(3)){#Clayton
    result <- copula::gofCopula(copula=copula::claytonCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    fam <- 3
  }
  else if (fitted_copula$family %in% c(14)) {#180 Gumbel
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::gumbelCopula(dim=2,param = fitted_copula$par)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    fam <- 14
  }
  else if (fitted_copula$family %in% c(24)) {#90 Gumbel
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::gumbelCopula(dim=2,param = -fitted_copula$par),flip=c(TRUE,FALSE)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    fam <- 24
  }
  else if (fitted_copula$family %in% c(34)) {#270 Gumbel
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::gumbelCopula(dim=2,param = -fitted_copula$par),flip=c(FALSE,TRUE)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    fam <- 34
  }
  else if (fitted_copula$family %in% c(4)){#Gumbel
    result <- copula::gofCopula(copula=copula::gumbelCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    fam <- 4
  }
  else if (fitted_copula$family %in% c(1)){#Normal
    result <- copula::gofCopula(copula=copula::ellipCopula(family="normal",dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    fam <- 1
  }
  else if (fitted_copula$family %in% c(2)){#Student
    result <- copula::gofCopula(copula=copula::ellipCopula(family="t",dim=2,param=fitted_copula$par,df=round(fitted_copula$par2),df.fixed=TRUE), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    fam <- 2
  }
  else if (fitted_copula$family %in% c(5)){#Frank
    result <- copula::gofCopula(copula=copula::frankCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    fam <- 5
  }
  else if (fitted_copula$family %in% c(16)) {#180 Joe
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::joeCopula(dim=2,param = fitted_copula$par)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    fam <- 16
  }
  else if (fitted_copula$family %in% c(26)) {#90 Joe
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::joeCopula(dim=2,param = -fitted_copula$par),flip=c(TRUE,FALSE)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    fam <- 26
  }
  else if (fitted_copula$family %in% c(36)) {#270 Joe
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::joeCopula(dim=2,param = -fitted_copula$par),flip=c(FALSE,TRUE)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    fam <- 36
  }
  else if (fitted_copula$family %in% c(6)){#Joe
    result <- copula::gofCopula(copula=copula::joeCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    fam <- 6
  }
  #print(fam)
  best.weights$Rn.tau*result+(1-best.weights$Rn.tau)*abs(cor(u1, u2, method = "kendall", use = "complete.obs"))
}
test.weighted.best$Rn.CvM <- function(u1,u2,weights){
  fitted_copula <- suppressWarnings(VineCopula::BiCopSelect(u1=u1,u2=u2,familyset = c(1,2,3,4,5,6)))
  if (fitted_copula$family %in% c(13)) {#180 Clayton
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::claytonCopula(dim=2,param = fitted_copula$par)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    result2 <- copula::gofCopula(copula=copula::rotCopula(copula = copula::claytonCopula(dim=2,param = fitted_copula$par)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 13
  }
  else if (fitted_copula$family %in% c(23)) {#90 Clayton
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::claytonCopula(dim=2,param = -fitted_copula$par),flip=c(TRUE,FALSE)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    result2 <- copula::gofCopula(copula=copula::rotCopula(copula = copula::claytonCopula(dim=2,param = -fitted_copula$par),flip=c(TRUE,FALSE)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 23
  }
  else if (fitted_copula$family %in% c(33)) {#270 Clayton
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::claytonCopula(dim=2,param = -fitted_copula$par),flip=c(FALSE,TRUE)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    result2 <- copula::gofCopula(copula=copula::rotCopula(copula = copula::claytonCopula(dim=2,param = -fitted_copula$par),flip=c(FALSE,TRUE)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 33
  }
  else if (fitted_copula$family %in% c(3)){#Clayton
    result <- copula::gofCopula(copula=copula::claytonCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    result2 <- copula::gofCopula(copula=copula::claytonCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 3
  }
  else if (fitted_copula$family %in% c(14)) {#180 Gumbel
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::gumbelCopula(dim=2,param = fitted_copula$par)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    result2 <- copula::gofCopula(copula=copula::rotCopula(copula = copula::gumbelCopula(dim=2,param = fitted_copula$par)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 14
  }
  else if (fitted_copula$family %in% c(24)) {#90 Gumbel
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::gumbelCopula(dim=2,param = -fitted_copula$par),flip=c(TRUE,FALSE)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    result2 <- copula::gofCopula(copula=copula::rotCopula(copula = copula::gumbelCopula(dim=2,param = -fitted_copula$par),flip=c(TRUE,FALSE)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 24
  }
  else if (fitted_copula$family %in% c(34)) {#270 Gumbel
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::gumbelCopula(dim=2,param = -fitted_copula$par),flip=c(FALSE,TRUE)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    result2 <- copula::gofCopula(copula=copula::rotCopula(copula = copula::gumbelCopula(dim=2,param = -fitted_copula$par),flip=c(FALSE,TRUE)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 34
  }
  else if (fitted_copula$family %in% c(4)){#Gumbel
    result <- copula::gofCopula(copula=copula::gumbelCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    result2 <- copula::gofCopula(copula=copula::gumbelCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 4
  }
  else if (fitted_copula$family %in% c(1)){#Normal
    result <- copula::gofCopula(copula=copula::ellipCopula(family="normal",dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    result2 <- copula::gofCopula(copula=copula::ellipCopula(family="normal",dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 1
  }
  else if (fitted_copula$family %in% c(2)){#Student
    result <- copula::gofCopula(copula=copula::ellipCopula(family="t",dim=2,param=fitted_copula$par,df=round(fitted_copula$par2),df.fixed=TRUE), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    result2 <- copula::gofCopula(copula=copula::ellipCopula(family="t",dim=2,param=fitted_copula$par,df=round(fitted_copula$par2),df.fixed=TRUE), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 2
  }
  else if (fitted_copula$family %in% c(5)){#Frank
    result <- copula::gofCopula(copula=copula::frankCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    result2 <- copula::gofCopula(copula=copula::frankCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 5
  }
  else if (fitted_copula$family %in% c(16)) {#180 Joe
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::joeCopula(dim=2,param = fitted_copula$par)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    result2 <- copula::gofCopula(copula=copula::rotCopula(copula = copula::joeCopula(dim=2,param = fitted_copula$par)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 16
  }
  else if (fitted_copula$family %in% c(26)) {#90 Joe
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::joeCopula(dim=2,param = -fitted_copula$par),flip=c(TRUE,FALSE)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    result2 <- copula::gofCopula(copula=copula::rotCopula(copula = copula::joeCopula(dim=2,param = -fitted_copula$par),flip=c(TRUE,FALSE)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 26
  }
  else if (fitted_copula$family %in% c(36)) {#270 Joe
    result <- copula::gofCopula(copula=copula::rotCopula(copula = copula::joeCopula(dim=2,param = -fitted_copula$par),flip=c(FALSE,TRUE)), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    result2 <- copula::gofCopula(copula=copula::rotCopula(copula = copula::joeCopula(dim=2,param = -fitted_copula$par),flip=c(FALSE,TRUE)), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 36
  }
  else if (fitted_copula$family %in% c(6)){#Joe
    result <- copula::gofCopula(copula=copula::joeCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult",method="Rn")$p.value
    result2 <- copula::gofCopula(copula=copula::joeCopula(dim=2,param=fitted_copula$par), x=cbind(u1,u2),test.method="single", simulation="mult")$p.value
    fam <- 6
  }
  #print(fam)
  #result=Rn, result2=CvM
  best.weights$Rn.CvM*result + (1-best.weights$Rn.CvM)*result2
}


# 7. Find best weighting parameters for ranked tests ----------------------

# implement ranked tests

RVineStructureSelect_2ranked_alpha <- function (alpha=0.5, data, familyset = NA, type = 0, selectioncrit = "AIC",
                                                indeptest = FALSE, level = 0.05, trunclevel = NA, progress = FALSE,
                                                weights = NA, treecrit = "tau",treecrit2="tau", rotations = TRUE, se = FALSE,
                                                presel = TRUE, method = "mle", cores = 1)
{
  args <- preproc(c(as.list(environment()), call = match.call()),
                  check_data, check_nobs, check_if_01, prep_familyset,
                  check_twoparams)
  list2env(args, environment())
  d <- ncol(data)
  n <- nrow(data)
  if (d < 2)
    stop("Dimension has to be at least 2.")
  if (d == 2) {
    return(RVineCopSelect(data, familyset = familyset, Matrix = matrix(c(2,
                                                                         1, 0, 1), 2, 2), selectioncrit = selectioncrit,
                          indeptest = indeptest, level = level, trunclevel = trunclevel,
                          weights = weights, rotations = FALSE, se = se, presel = presel,
                          method = method, cores = cores))
  }
  if (!(selectioncrit %in% c("AIC", "BIC", "logLik")))
    stop("Selection criterion not implemented.")
  if (level < 0 & level > 1)
    stop("Significance level has to be between 0 and 1.")
  if (type == 0)
    type <- "RVine"
  if (type == 1)
    type <- "CVine"
  treecrit <- set_treecrit(treecrit, familyset)
  treecrit2 <- set_treecrit2(treecrit2, familyset)
  if (is.na(trunclevel))
    trunclevel <- d
  if (trunclevel == 0)
    familyset <- 0
  RVine <- list(Tree = NULL, Graph = NULL)
  warn <- NULL
  g <- initializeFirstGraph_2ranked_alpha(alpha, data, treecrit, treecrit2, weights)
  MST <- findMaxTree(g, mode = type)
  if (cores != 1 | is.na(cores)) {
    if (is.na(cores))
      cores <- max(1, detectCores() - 1)
    if (cores > 1) {
      cl <- makePSOCKcluster(cores)
      setDefaultCluster(cl)
      on.exit(try(stopCluster(cl), silent = TRUE))
    }
  }
  VineTree <- fit.FirstTreeCopulas(MST, data, familyset = familyset,
                                   selectioncrit = selectioncrit, indeptest = indeptest,
                                   level = level, weights = weights, se = se, presel = presel,
                                   method = method, cores = cores)
  if (!is.null(VineTree$warn))
    warn <- VineTree$warn
  RVine$Tree[[1]] <- VineTree
  RVine$Graph[[1]] <- g
  oldVineGraph <- VineTree
  for (tree in 2:(d - 1)) {
    RVine$Tree[[tree - 1]]$E$Copula.CondData.1 <- NULL
    RVine$Tree[[tree - 1]]$E$Copula.CondData.2 <- NULL
    if (trunclevel == tree - 1)
      familyset <- 0
    g <- buildNextGraph_2ranked_alpha(alpha,VineTree, weights, treecrit = treecrit,treecrit2=treecrit2,
                                      cores > 1, truncated = trunclevel < tree)
    MST <- findMaxTree(g, mode = type, truncated = trunclevel <
                         tree)
    VineTree <- fit.TreeCopulas(MST, VineTree, familyset,
                                selectioncrit, indeptest, level, progress, weights = weights,
                                se = se, presel = presel, method = method, cores = cores)
    if (!is.null(VineTree$warn))
      warn <- VineTree$warn
    RVine$Tree[[tree]] <- VineTree
    RVine$Graph[[tree]] <- g
  }
  if (any(is.na(data)))
    warning(" In ", args$call[1], ": ", "Some of the data are NA. ",
            "Only pairwise complete observations are used.",
            call. = FALSE)
  if (!is.null(warn))
    warning(" In ", args$call[1], ": ", warn, call. = FALSE)
  .RVine <- RVine
  .data <- data
  .callexp <- match.call()
  rm(list = ls())
  as.RVM2(.RVine, .data, .callexp)
}


environment(RVineStructureSelect_2ranked_alpha) <- asNamespace('VineCopula')


initializeFirstGraph_2ranked_alpha <- function (alpha, data, treecrit, treecrit2, weights)
{
  all.pairs <- combn(1:ncol(data), 2)
  edge.ws <- apply(all.pairs, 2, function(ind) treecrit(data[,
                                                             ind[1]], data[, ind[2]], weights))
  edge.ws2 <- apply(all.pairs, 2, function(ind) treecrit2(data[,
                                                               ind[1]], data[, ind[2]], weights))
  rel.nobs <- apply(all.pairs, 2, function(ind) mean(!is.na(data[,
                                                                 ind[1]] + data[, ind[2]])))
  weightmatrix <- as.data.frame(cbind(rank(edge.ws,na.last = "keep"),rank(edge.ws2,na.last = "keep")))
  for (i in 1:length(weightmatrix[,1])){
    weightmatrix[i,3] <- alpha*weightmatrix[i,1]+(1-alpha)*weightmatrix[i,2]
  }
  edge.ws <- weightmatrix[,3]
  W <- diag(ncol(data))
  W[lower.tri(W)] <- edge.ws
  W <- t(W)
  colnames(W) <- rownames(W) <- colnames(data)
  graphFromWeightMatrix(W)
}

environment(initializeFirstGraph_2ranked_alpha) <- asNamespace('VineCopula')

getEdgeInfo_2ranked_alpha <- function (alpha,i, g, oldVineGraph, treecrit, treecrit2, weights, truncated = FALSE)
{
  con <- g$E$nums[i, ]
  temp <- oldVineGraph$E$nums[con, ]
  ok <- FALSE
  if ((temp[1, 1] == temp[2, 1]) || (temp[1, 2] == temp[2,
                                                        1])) {
    ok <- TRUE
    same <- temp[2, 1]
  }
  else {
    if ((temp[1, 1] == temp[2, 2]) || (temp[1, 2] == temp[2,
                                                          2])) {
      ok <- TRUE
      same <- temp[2, 2]
    }
  }
  w <- nedSet <- ningSet <- name <- NA
  w2 <- nedSet <- ningSet <- name <- NA
  todel <- TRUE
  if (ok) {
    l1 <- c(g$V$conditionedSet[[con[1]]], g$V$conditioningSet[[con[1]]])
    l2 <- c(g$V$conditionedSet[[con[2]]], g$V$conditioningSet[[con[2]]])
    nedSet <- c(setdiff(l1, l2), setdiff(l2, l1))
    ningSet <- intersect(l1, l2)
    todel <- FALSE
    if (truncated == FALSE) {
      if (temp[1, 1] == same) {
        zr1 <- oldVineGraph$E$Copula.CondData.2[con[1]]
      }
      else {
        zr1 <- oldVineGraph$E$Copula.CondData.1[con[1]]
      }
      if (temp[2, 1] == same) {
        zr2 <- oldVineGraph$E$Copula.CondData.2[con[2]]
      }
      else {
        zr2 <- oldVineGraph$E$Copula.CondData.1[con[2]]
      }
      if (is.list(zr1)) {
        zr1a <- as.vector(zr1[[1]])
        zr2a <- as.vector(zr2[[1]])
      }
      else {
        zr1a <- zr1
        zr2a <- zr2
      }
      keine_nas <- !(is.na(zr1a) | is.na(zr2a))
      w <- treecrit(zr1a[keine_nas], zr2a[keine_nas],
                    weights)
      w2 <- treecrit2(zr1a[keine_nas], zr2a[keine_nas],
                      weights)
      name.node1 <- strsplit(g$V$names[con[1]], split = " *[,;] *")[[1]]
      name.node2 <- strsplit(g$V$names[con[2]], split = " *[,;] *")[[1]]
      nmdiff <- c(setdiff(name.node1, name.node2), setdiff(name.node2,
                                                           name.node1))
      nmsect <- intersect(name.node1, name.node2)
      name <- paste(paste(nmdiff, collapse = ","), paste(nmsect,
                                                         collapse = ","), sep = " ; ")
    }
    else {
      w <- 1
      w2 <- 1
    }
  }
  list(w = w, w2=w2, nedSet = nedSet, ningSet = ningSet, name = name,
       todel = todel)
}

environment(getEdgeInfo_2ranked_alpha) <- asNamespace('VineCopula')

buildNextGraph_2ranked_alpha <- function (alpha,oldVineGraph, treecrit, treecrit2, weights = NA, parallel, truncated = FALSE)
{
  d <- nrow(oldVineGraph$E$nums)
  g <- makeFullGraph(d)
  g$V$names <- oldVineGraph$E$names
  g$V$conditionedSet <- oldVineGraph$E$conditionedSet
  g$V$conditioningSet <- oldVineGraph$E$conditioningSet
  if (parallel)
    lapply <- function(...) parallel::parLapply(getDefaultCluster(),
                                                ...)
  out <- lapply(seq_len(nrow(g$E$nums)), getEdgeInfo_2ranked_alpha, g = g,alpha=alpha,
                oldVineGraph = oldVineGraph, treecrit = treecrit, treecrit2=treecrit2, weights = weights,
                truncated = truncated)
  weightmatrix <- as.data.frame(cbind(rank(sapply(out, function(x) x$w),na.last = "keep"),rank(sapply(out, function(x) x$w2),na.last = "keep")))
  for (i in 1:length(weightmatrix[,1])){
    weightmatrix[i,3] <- alpha*weightmatrix[i,1]+(1-alpha)*weightmatrix[i,2]
  }
  g$E$weights <- weightmatrix[,3]
  g$E$names <- sapply(out, function(x) x$name)
  g$E$conditionedSet <- lapply(out, function(x) x$nedSet)
  g$E$conditioningSet <- lapply(out, function(x) x$ningSet)
  g$E$todel <- sapply(out, function(x) x$todel)
  deleteEdges(g)
}

environment(buildNextGraph_2ranked_alpha) <- asNamespace('VineCopula')

set_treecrit2 <- function (treecrit2, famset)
{
  if (is.function(treecrit2)) {
    w <- try(treecrit2(u1 = runif(10), u2 = runif(10), weights = rep(1,
                                                                     10)), silent = TRUE)
    if (inherits(w, "try-error"))
      stop("treecrit2 must be of the form 'function(u1, u2, weights)'")
    if (!is.numeric(w) || length(w) > 1)
      stop("treecrit2 does not return a numeric scalar")
  }
  else if (all(treecrit2 == "tau")) {
    treecrit2 <- function(u1, u2, weights) {
      complete.i <- which(!is.na(u1 + u2))
      if (length(complete.i) < 10) {
        tau <- 0
      }
      else {
        complete.freq <- mean(!is.na(u1 + u2))
        tau <- abs(fasttau(u1[complete.i], u2[complete.i],
                           weights))
        tau * sqrt(complete.freq)
      }
    }
  }
  else if (all(treecrit2 == "rho")) {
    treecrit2 <- function(u1, u2, weights) {
      complete.i <- which(!is.na(u1 + u2))
      if (length(complete.i) < 10) {
        tau <- 0
      }
      else {
        complete.freq <- mean(!is.na(u1 + u2))
        rho <- abs(cor(u1, u2, method = "spearman",
                       use = "complete.obs"))
        rho * sqrt(complete.freq)
      }
    }
  }
  else if (all(treecrit2 == "AIC")) {
    treecrit2 <- function(u1, u2, weights) {
      complete.i <- which(!is.na(u1 + u2))
      if (length(complete.i) < 2) {
        AIC <- 0
      }
      else {
        AIC <- -suppressWarnings(BiCopSelect(u1[complete.i],
                                             u2[complete.i], familyset = famset, weights = weights)$AIC)
      }
      AIC
    }
  }
  else if (all(treecrit2 == "BIC")) {
    treecrit2 <- function(u1, u2, weights) {
      complete.i <- which(!is.na(u1 + u2))
      if (length(complete.i) < 2) {
        BIC <- 0
      }
      else {
        BIC <- -suppressWarnings(BiCopSelect(u1[complete.i],
                                             u2[complete.i], familyset = famset, weights = weights)$BIC)
      }
      BIC
    }
  }
  else if (all(treecrit2 == "cAIC")) {
    treecrit2 <- function(u1, u2, weights) {
      complete.i <- which(!is.na(u1 + u2))
      if (length(complete.i) < 2) {
        cAIC <- 0
      }
      else {
        fit <- suppressWarnings(BiCopSelect(u1[complete.i],
                                            u2[complete.i], familyset = famset, weights = weights))
        n <- fit$nobs
        p <- fit$npars
        cAIC <- -(fit$AIC + 2 * p * (p + 1)/(n - p -
                                               1))
      }
      cAIC
    }
  }
  else {
    txt1 <- "treecrit2 must be one of \"tau\", \"rho\", \"AIC\", \"BIC\", \"cAIC\""
    txt2 <- "or a function like function(u1, u2, weights) ... returning a numeric value."
    stop(paste(txt1, txt2))
  }
  treecrit2
}
environment(set_treecrit2) <- asNamespace('VineCopula')


#Find best weights
find.weights.ranked.results <- list()
cl <- parallel::makeCluster(cores2)
doParallel::registerDoParallel(cl)
clusterExport(cl=cl, c('RVineStructureSelect_2ranked_alpha', 'initializeFirstGraph_2ranked_alpha','getEdgeInfo_2ranked_alpha',
                       'buildNextGraph_2ranked_alpha','set_treecrit2'))
find.weights.ranked.results$Rn.tau <- foreach(k = 1:10,.packages = "VineCopula") %:%
  foreach(i=1:length(par)) %dopar% {
  set.seed(2021)
    zvalues <- matrix(nrow=WE-9, ncol=ncol_assets)#nrow_assets
    for(l in 1:ncol_assets){#ncol_assets
      if("skew" %in% names(GARCH[[l]][[k]]@fit$coef)){
        zvalues[,l] = fGarch::psstd(GARCH.for.weights[[l]][[k]]@fit$z, nu=GARCH.for.weights[[l]][[k]]@fit$coef[8],
                                    xi = GARCH.for.weights[[l]][[k]]@fit$coef[7])
      }else{
        zvalues[,l] = fGarch::pstd(GARCH.for.weights[[j]][[i]]@fit$z, nu=GARCH.for.weights[[j]][[k]]@fit$coef[7])}
    }

  RVineStructureSelect_2ranked_alpha(alpha=par[i],data=zvalues, familyset = c(1,2,3,4,5,6), treecrit = test.unweighted$Rn.test, treecrit2="tau")
}
parallel::stopCluster(cl)


cl <- parallel::makeCluster(cores2)
doParallel::registerDoParallel(cl)
clusterExport(cl=cl, c('RVineStructureSelect_2ranked_alpha', 'initializeFirstGraph_2ranked_alpha','getEdgeInfo_2ranked_alpha',
                       'buildNextGraph_2ranked_alpha','set_treecrit2'))
find.weights.ranked.results$CvM.tau <- foreach(k = 1:10,.packages = "VineCopula") %:%
  foreach(i=1:length(par)) %dopar% {
    set.seed(2021)
    zvalues <- matrix(nrow=WE-9, ncol=ncol_assets)#nrow_assets
    for(l in 1:ncol_assets){#ncol_assets
      if("skew" %in% names(GARCH[[l]][[k]]@fit$coef)){
        zvalues[,l] = fGarch::psstd(GARCH.for.weights[[l]][[k]]@fit$z, nu=GARCH.for.weights[[l]][[k]]@fit$coef[8],
                                    xi = GARCH.for.weights[[l]][[k]]@fit$coef[7])
      }else{
        zvalues[,l] = fGarch::pstd(GARCH.for.weights[[j]][[i]]@fit$z, nu=GARCH.for.weights[[j]][[k]]@fit$coef[7])}
    }

    RVineStructureSelect_2ranked_alpha(alpha=par[i],data=zvalues, familyset = c(1,2,3,4,5,6), treecrit = test.unweighted$CvM.test, treecrit2="tau")
  }
parallel::stopCluster(cl)

cl <- parallel::makeCluster(cores2)
doParallel::registerDoParallel(cl)
clusterExport(cl=cl, c('RVineStructureSelect_2ranked_alpha', 'initializeFirstGraph_2ranked_alpha','getEdgeInfo_2ranked_alpha',
                       'buildNextGraph_2ranked_alpha','set_treecrit2'))
find.weights.ranked.results$Rn.AIC <- foreach(k = 1:10,.packages = "VineCopula") %:%
  foreach(i=1:length(par)) %dopar% {
    set.seed(2021)
    zvalues <- matrix(nrow=WE-9, ncol=ncol_assets)#nrow_assets
    for(l in 1:ncol_assets){#ncol_assets
      if("skew" %in% names(GARCH[[l]][[k]]@fit$coef)){
        zvalues[,l] = fGarch::psstd(GARCH.for.weights[[l]][[k]]@fit$z, nu=GARCH.for.weights[[l]][[k]]@fit$coef[8],
                                    xi = GARCH.for.weights[[l]][[k]]@fit$coef[7])
      }else{
        zvalues[,l] = fGarch::pstd(GARCH.for.weights[[j]][[i]]@fit$z, nu=GARCH.for.weights[[j]][[k]]@fit$coef[7])}
    }

    RVineStructureSelect_2ranked_alpha(alpha=par[i],data=zvalues, familyset = c(1,2,3,4,5,6), treecrit = test.unweighted$Rn.test, treecrit2="AIC")
  }
parallel::stopCluster(cl)


cl <- parallel::makeCluster(cores2)
doParallel::registerDoParallel(cl)
clusterExport(cl=cl, c('RVineStructureSelect_2ranked_alpha', 'initializeFirstGraph_2ranked_alpha','getEdgeInfo_2ranked_alpha',
                       'buildNextGraph_2ranked_alpha','set_treecrit2'))
find.weights.ranked.results$CvM.AIC <- foreach(k = 1:10,.packages = "VineCopula") %:%
  foreach(i=1:length(par)) %dopar% {
    set.seed(2021)
    zvalues <- matrix(nrow=WE-9, ncol=ncol_assets)#nrow_assets
    for(l in 1:ncol_assets){#ncol_assets
      if("skew" %in% names(GARCH[[l]][[k]]@fit$coef)){
        zvalues[,l] = fGarch::psstd(GARCH.for.weights[[l]][[k]]@fit$z, nu=GARCH.for.weights[[l]][[k]]@fit$coef[8],
                                    xi = GARCH.for.weights[[l]][[k]]@fit$coef[7])
      }else{
        zvalues[,l] = fGarch::pstd(GARCH.for.weights[[j]][[i]]@fit$z, nu=GARCH.for.weights[[j]][[k]]@fit$coef[7])}
    }

    RVineStructureSelect_2ranked_alpha(alpha=par[i],data=zvalues, familyset = c(1,2,3,4,5,6), treecrit = test.unweighted$CvM.test, treecrit2="AIC")
  }
parallel::stopCluster(cl)



cl <- parallel::makeCluster(cores2)
doParallel::registerDoParallel(cl)
clusterExport(cl=cl, c('RVineStructureSelect_2ranked_alpha', 'initializeFirstGraph_2ranked_alpha','getEdgeInfo_2ranked_alpha',
                       'buildNextGraph_2ranked_alpha','set_treecrit2'))
find.weights.ranked.results$Rn.CvM <- foreach(k = 1:10,.packages = "VineCopula") %:%
  foreach(i=1:length(par)) %dopar% {
    set.seed(2021)
    zvalues <- matrix(nrow=WE-9, ncol=ncol_assets)#nrow_assets
    for(l in 1:ncol_assets){#ncol_assets
      if("skew" %in% names(GARCH[[l]][[k]]@fit$coef)){
        zvalues[,l] = fGarch::psstd(GARCH.for.weights[[l]][[k]]@fit$z, nu=GARCH.for.weights[[l]][[k]]@fit$coef[8],
                                    xi = GARCH.for.weights[[l]][[k]]@fit$coef[7])
      }else{
        zvalues[,l] = fGarch::pstd(GARCH.for.weights[[j]][[i]]@fit$z, nu=GARCH.for.weights[[j]][[k]]@fit$coef[7])}
    }

    RVineStructureSelect_2ranked_alpha(alpha=par[i],data=zvalues, familyset = c(1,2,3,4,5,6), treecrit = test.unweighted$Rn.test, treecrit2=test.unweighted$CvM.test)
  }
parallel::stopCluster(cl)

## Plots

find.weights.ranked.aics <- list()
find.weights.ranked.plots <- list()

#Rn.tau
find.weights.ranked.aics$Rn.tau <- matrix(nrow = length(par), ncol = length(find.weights.ranked.results$Rn.tau))

for(i in 1:length(par)){
  for(j in 1:length(find.weights.ranked.results$Rn.tau)){
    find.weights.ranked.aics$Rn.tau[i,j] = find.weights.ranked.results$Rn.tau[[j]][[i]]$AIC
  }
}
colnames(find.weights.ranked.aics$Rn.tau) <- seq(ncol(find.weights.ranked.aics$Rn.tau))
find.weights.ranked.aics$Rn.tau <- as.data.frame(find.weights.ranked.aics$Rn.tau)
find.weights.ranked.aics$Rn.tau$parameter <- par


min <- find.weights.ranked.aics$Rn.tau%>%
  pivot_longer(cols = 1:10, names_to = "period")%>%
  arrange(as.numeric(period))%>%
  group_by(period)%>%
  dplyr::filter(value == min(value))

find.weights.ranked.plots$Rn.tau <- find.weights.ranked.aics$Rn.tau%>%
  pivot_longer(cols = 1:10, names_to = "period")%>%
  arrange(as.numeric(period))%>%
  ggplot(data = ., aes(x = parameter, y = value, colour = period))+
  scale_x_continuous(breaks = seq(0, 1, by = 0.1))+
  geom_line()+
  geom_point(data = min, aes(x=parameter, y = value, group = period), colour = "black", alpha = 1, shape = 2)+
  theme_bw()+
  theme(legend.position = "bottom")+
  labs(y = "AIC", title = TeX("$\\eta_{R_n,\\tau,r}f(p_{R_n})+(1-\\eta_{R_n,\\tau,r})f(\\tau)$"), x = TeX("$\\eta_{R_n,\\tau,r}$"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(aspect.ratio = 0.8)+
  viridis::scale_color_viridis(discrete = T)+
  theme(legend.position = "none")

find.weights.ranked.plots$Rn.tau.legend <- find.weights.ranked.aics$Rn.tau%>%
  pivot_longer(cols = 1:10, names_to = "period")%>%
  arrange(as.numeric(period))%>%
  ggplot(data = ., aes(x = parameter, y = value, colour = period))+
  scale_x_continuous(breaks = seq(0, 1, by = 0.1))+
  geom_line()+
  geom_point(data = min, aes(x=parameter, y = value, group = period), colour = "black", alpha = 1, shape = 2)+
  theme_bw()+
  theme(legend.position = "bottom")+
  labs(y = "AIC", title = TeX("$\\eta_{R_n,\\tau,r}f(p_{R_n})+(1-\\eta_{R_n,\\tau,r})f(\\tau)$"), x = TeX("$\\eta_{R_n,\\tau,r}$"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(aspect.ratio = 0.8)+
  viridis::scale_color_viridis(discrete = T)


#CvM.tau
find.weights.ranked.aics$CvM.tau <- matrix(nrow = length(par), ncol = length(find.weights.ranked.results$CvM.tau))

for(i in 1:length(par)){
  for(j in 1:length(find.weights.ranked.results$CvM.tau)){
    find.weights.ranked.aics$CvM.tau[i,j] = find.weights.ranked.results$CvM.tau[[j]][[i]]$AIC
  }
}
colnames(find.weights.ranked.aics$CvM.tau) <- seq(ncol(find.weights.ranked.aics$CvM.tau))
find.weights.ranked.aics$CvM.tau <- as.data.frame(find.weights.ranked.aics$CvM.tau)
find.weights.ranked.aics$CvM.tau$parameter <- par


min <- find.weights.ranked.aics$CvM.tau%>%
  pivot_longer(cols = 1:10, names_to = "period")%>%
  arrange(as.numeric(period))%>%
  group_by(period)%>%
  dplyr::filter(value == min(value))

find.weights.ranked.plots$CvM.tau <- find.weights.ranked.aics$CvM.tau%>%
  pivot_longer(cols = 1:10, names_to = "period")%>%
  arrange(as.numeric(period))%>%
  ggplot(data = ., aes(x = parameter, y = value, colour = period))+
  scale_x_continuous(breaks = seq(0, 1, by = 0.1))+
  geom_line()+
  geom_point(data = min, aes(x=parameter, y = value, group = period), colour = "black", alpha = 1, shape = 2)+
  theme_bw()+
  theme(legend.position = "bottom")+
  labs(y = "AIC", title = TeX("$\\eta_{S_n,\\tau,r}f(p_{S_n})+(1-\\eta_{S_n,\\tau,r})f(\\tau)$"), x = TeX("$\\eta_{S_n,\\tau,r}$"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(aspect.ratio = 0.8)+
  viridis::scale_color_viridis(discrete = T)+
  theme(legend.position = "none")


#Rn.AIC
find.weights.ranked.aics$Rn.AIC <- matrix(nrow = length(par), ncol = length(find.weights.ranked.results$Rn.AIC))

for(i in 1:length(par)){
  for(j in 1:length(find.weights.ranked.results$Rn.AIC)){
    find.weights.ranked.aics$Rn.AIC[i,j] = find.weights.ranked.results$Rn.AIC[[j]][[i]]$AIC
  }
}
colnames(find.weights.ranked.aics$Rn.AIC) <- seq(ncol(find.weights.ranked.aics$Rn.AIC))
find.weights.ranked.aics$Rn.AIC <- as.data.frame(find.weights.ranked.aics$Rn.AIC)
find.weights.ranked.aics$Rn.AIC$parameter <- par


min <- find.weights.ranked.aics$Rn.AIC%>%
  pivot_longer(cols = 1:10, names_to = "period")%>%
  arrange(as.numeric(period))%>%
  group_by(period)%>%
  dplyr::filter(value == min(value))

find.weights.ranked.plots$Rn.AIC <- find.weights.ranked.aics$Rn.AIC%>%
  pivot_longer(cols = 1:10, names_to = "period")%>%
  arrange(as.numeric(period))%>%
  ggplot(data = ., aes(x = parameter, y = value, colour = period))+
  scale_x_continuous(breaks = seq(0, 1, by = 0.1))+
  geom_line()+
  geom_point(data = min, aes(x=parameter, y = value, group = period), colour = "black", alpha = 1, shape = 2)+
  theme_bw()+
  theme(legend.position = "bottom")+
  labs(y = "AIC", title = TeX("$\\eta_{R_n,AIC,r}f(p_{R_n})+(1-\\eta_{R_n,AIC,r})f(\\AIC)$"), x = TeX("$\\eta_{R_n,AIC,r}$"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(aspect.ratio = 0.8)+
  viridis::scale_color_viridis(discrete = T)+
  theme(legend.position = "none")


#CvM.AIC
find.weights.ranked.aics$CvM.AIC <- matrix(nrow = length(par), ncol = length(find.weights.ranked.results$CvM.AIC))

for(i in 1:length(par)){
  for(j in 1:length(find.weights.ranked.results$CvM.AIC)){
    find.weights.ranked.aics$CvM.AIC[i,j] = find.weights.ranked.results$CvM.AIC[[j]][[i]]$AIC
  }
}
colnames(find.weights.ranked.aics$CvM.AIC) <- seq(ncol(find.weights.ranked.aics$CvM.AIC))
find.weights.ranked.aics$CvM.AIC <- as.data.frame(find.weights.ranked.aics$CvM.AIC)
find.weights.ranked.aics$CvM.AIC$parameter <- par


min <- find.weights.ranked.aics$CvM.AIC%>%
  pivot_longer(cols = 1:10, names_to = "period")%>%
  arrange(as.numeric(period))%>%
  group_by(period)%>%
  dplyr::filter(value == min(value))

find.weights.ranked.plots$CvM.AIC <- find.weights.ranked.aics$CvM.AIC%>%
  pivot_longer(cols = 1:10, names_to = "period")%>%
  arrange(as.numeric(period))%>%
  ggplot(data = ., aes(x = parameter, y = value, colour = period))+
  scale_x_continuous(breaks = seq(0, 1, by = 0.1))+
  geom_line()+
  geom_point(data = min, aes(x=parameter, y = value, group = period), colour = "black", alpha = 1, shape = 2)+
  theme_bw()+
  theme(legend.position = "bottom")+
  labs(y = "AIC", title = TeX("$\\eta_{S_n,AIC,r}f(p_{S_n})+(1-\\eta_{S_n,AIC,r})f(\\AIC)$"), x = TeX("$\\eta_{S_n,AIC,r}$"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(aspect.ratio = 0.8)+
  viridis::scale_color_viridis(discrete = T)+
  theme(legend.position = "none")

#Rn.CvM
find.weights.ranked.aics$Rn.CvM <- matrix(nrow = length(par), ncol = length(find.weights.ranked.results$Rn.CvM))

for(i in 1:length(par)){
  for(j in 1:length(find.weights.ranked.results$Rn.CvM)){
    find.weights.ranked.aics$Rn.CvM[i,j] = find.weights.ranked.results$Rn.CvM[[j]][[i]]$AIC
  }
}
colnames(find.weights.ranked.aics$Rn.CvM) <- seq(ncol(find.weights.ranked.aics$Rn.CvM))
find.weights.ranked.aics$Rn.CvM <- as.data.frame(find.weights.ranked.aics$Rn.CvM)
find.weights.ranked.aics$Rn.CvM$parameter <- par


min <- find.weights.ranked.aics$Rn.CvM%>%
  pivot_longer(cols = 1:10, names_to = "period")%>%
  arrange(as.numeric(period))%>%
  group_by(period)%>%
  dplyr::filter(value == min(value))

g <- find.weights.ranked.aics$Rn.CvM%>%
  pivot_longer(cols = 1:10, names_to = "period")%>%
  arrange(as.numeric(period))%>%
  ggplot(data = ., aes(x = parameter, y = value, colour = period))+
  scale_x_continuous(breaks = seq(0, 1, by = 0.1))+
  geom_line()+
  geom_point(data = min, aes(x=parameter, y = value, group = period), colour = "black", alpha = 0.4)+
  theme_bw()+
  theme(legend.position = "bottom")+
  labs(title = "Rn.CvM ranked")+
  theme(legend.position = "none")



find.weights.ranked.plots$legend <- get_only_legend(find.weights.ranked.plots$Rn.tau.legend)
combined_plot <- grid.arrange(find.weights.ranked.plots$CvM.tau,find.weights.ranked.plots$Rn.tau,
                              find.weights.ranked.plots$CvM.AIC, find.weights.ranked.plots$Rn.AIC, ncol = 2, nrow = 2)
g <- grid.arrange(combined_plot, find.weights.ranked.plots$legend, nrow = 2, heights = c(10,1))
ggsave("find.weights.ranked.pdf", plot = g, path = "../MA/Plots/VaR", units = "cm", width = 20)


# save best weights for ranked tests
best.weights.ranked <- list()
best.weights.ranked$Rn.tau <- 0.2
best.weights.ranked$CvM.tau <- 0.1
best.weights.ranked$Rn.AIC <- 0.2
best.weights.ranked$CvM.AIC <- 0.1
best.weights.ranked$Rn.CvM <- 0.1


# 8. Vines ----------------------------------------------------------------

# Fit vines with all tree selection criteria

Vines <- list()


Vines$tau <- vector("list", length = TComplete-WE)
cl <- parallel::makeCluster(cores2)
doParallel::registerDoParallel(cl)
Vines$tau <- foreach(i = 1:(TComplete-WE)) %dopar% { #(T-WE)
  set.seed(2021)
  zvalues <- matrix(nrow=WE, ncol=ncol_assets)#nrow_assets
  for(j in 1:ncol_assets){#ncol_assets
    if("skew" %in% names(GARCH[[j]][[i]]@fit$coef)){
      zvalues[,j] = fGarch::psstd(GARCH[[j]][[i]]@fit$z, nu=GARCH[[j]][[i]]@fit$coef[8],
                                  xi = GARCH[[j]][[i]]@fit$coef[7])
    }else{
      zvalues[,j] = fGarch::pstd(GARCH[[j]][[i]]@fit$z, nu=GARCH[[j]][[i]]@fit$coef[7])}
  }
  VineCopula::RVineStructureSelect(data=zvalues,familyset = c(1,2,3,4,5,6), treecrit = test.unweighted$tau, indeptest = FALSE)
}
parallel::stopCluster(cl)

Vines$AIC <- vector("list", length = TComplete-WE)
cl <- parallel::makeCluster(cores2)
doParallel::registerDoParallel(cl)
Vines$AIC <- foreach(i = 1:(TComplete-WE)) %dopar% { #(T-WE)
  set.seed(2021)
  zvalues <- matrix(nrow=WE, ncol=ncol_assets)#nrow_assets
  for(j in 1:ncol_assets){#ncol_assets
    if("skew" %in% names(GARCH[[j]][[i]]@fit$coef)){
      zvalues[,j] = fGarch::psstd(GARCH[[j]][[i]]@fit$z, nu=GARCH[[j]][[i]]@fit$coef[8],
                                  xi = GARCH[[j]][[i]]@fit$coef[7])
    }else{
      zvalues[,j] = fGarch::pstd(GARCH[[j]][[i]]@fit$z, nu=GARCH[[j]][[i]]@fit$coef[7])}
  }
  VineCopula::RVineStructureSelect(data=zvalues,familyset = c(1,2,3,4,5,6), treecrit = test.unweighted$AIC, indeptest = FALSE)
}
parallel::stopCluster(cl)

Vines$Rn <- vector("list", length = TComplete-WE)
cl <- parallel::makeCluster(cores2)
doParallel::registerDoParallel(cl)
Vines$Rn <- foreach(i = 1:(TComplete-WE)) %dopar% { #(T-WE)
  set.seed(2021)
  zvalues <- matrix(nrow=WE, ncol=ncol_assets)#nrow_assets
  for(j in 1:ncol_assets){#ncol_assets
    if("skew" %in% names(GARCH[[j]][[i]]@fit$coef)){
      zvalues[,j] = fGarch::psstd(GARCH[[j]][[i]]@fit$z, nu=GARCH[[j]][[i]]@fit$coef[8],
                                  xi = GARCH[[j]][[i]]@fit$coef[7])
    }else{
      zvalues[,j] = fGarch::pstd(GARCH[[j]][[i]]@fit$z, nu=GARCH[[j]][[i]]@fit$coef[7])}
  }
  VineCopula::RVineStructureSelect(data=zvalues,familyset = c(1,2,3,4,5,6), treecrit = test.unweighted$Rn.test, indeptest = FALSE)
}
parallel::stopCluster(cl)

Vines$CvM <- vector("list", length = TComplete-WE)
cl <- parallel::makeCluster(cores2)
doParallel::registerDoParallel(cl)
Vines$CvM <- foreach(i = 1:(TComplete-WE)) %dopar% { #(T-WE)
  set.seed(2021)
  zvalues <- matrix(nrow=WE, ncol=ncol_assets)#nrow_assets
  for(j in 1:ncol_assets){#ncol_assets
    if("skew" %in% names(GARCH[[j]][[i]]@fit$coef)){
      zvalues[,j] = fGarch::psstd(GARCH[[j]][[i]]@fit$z, nu=GARCH[[j]][[i]]@fit$coef[8],
                                  xi = GARCH[[j]][[i]]@fit$coef[7])
    }else{
      zvalues[,j] = fGarch::pstd(GARCH[[j]][[i]]@fit$z, nu=GARCH[[j]][[i]]@fit$coef[7])}
  }
  VineCopula::RVineStructureSelect(data=zvalues,familyset = c(1,2,3,4,5,6), treecrit = test.unweighted$CvM.test, indeptest = FALSE)
}
parallel::stopCluster(cl)

Vines$Rn.tau.mult <- vector("list", length = TComplete-WE)
cl <- parallel::makeCluster(cores2)
doParallel::registerDoParallel(cl)
Vines$Rn.tau.mult <- foreach(i = 1:(TComplete-WE)) %dopar% { #(T-WE)
  set.seed(2021)
  zvalues <- matrix(nrow=WE, ncol=ncol_assets)#nrow_assets
  for(j in 1:ncol_assets){#ncol_assets
    if("skew" %in% names(GARCH[[j]][[i]]@fit$coef)){
      zvalues[,j] = fGarch::psstd(GARCH[[j]][[i]]@fit$z, nu=GARCH[[j]][[i]]@fit$coef[8],
                                  xi = GARCH[[j]][[i]]@fit$coef[7])
    }else{
      zvalues[,j] = fGarch::pstd(GARCH[[j]][[i]]@fit$z, nu=GARCH[[j]][[i]]@fit$coef[7])}
  }
  VineCopula::RVineStructureSelect(data=zvalues,familyset = c(1,2,3,4,5,6), treecrit = test.unweighted$Rn.tau.mult, indeptest = FALSE)
}
parallel::stopCluster(cl)

Vines$CvM.tau.mult <- vector("list", length = TComplete-WE)
cl <- parallel::makeCluster(cores2)
doParallel::registerDoParallel(cl)
Vines$CvM.tau.mult <- foreach(i = 1:(TComplete-WE)) %dopar% { #(T-WE)
  set.seed(2021)
  zvalues <- matrix(nrow=WE, ncol=ncol_assets)#nrow_assets
  for(j in 1:ncol_assets){#ncol_assets
    if("skew" %in% names(GARCH[[j]][[i]]@fit$coef)){
      zvalues[,j] = fGarch::psstd(GARCH[[j]][[i]]@fit$z, nu=GARCH[[j]][[i]]@fit$coef[8],
                                  xi = GARCH[[j]][[i]]@fit$coef[7])
    }else{
      zvalues[,j] = fGarch::pstd(GARCH[[j]][[i]]@fit$z, nu=GARCH[[j]][[i]]@fit$coef[7])}
  }
  VineCopula::RVineStructureSelect(data=zvalues,familyset = c(1,2,3,4,5,6), treecrit = test.unweighted$CvM.tau.mult, indeptest = FALSE)
}
parallel::stopCluster(cl)

Vines$CvM.Rn.mult <- vector("list", length = TComplete-WE)
cl <- parallel::makeCluster(cores2)
doParallel::registerDoParallel(cl)
Vines$CvM.Rn.mult <- foreach(i = 1:(TComplete-WE)) %dopar% { #(T-WE)
  set.seed(2021)
  zvalues <- matrix(nrow=WE, ncol=ncol_assets)#nrow_assets
  for(j in 1:ncol_assets){#ncol_assets
    if("skew" %in% names(GARCH[[j]][[i]]@fit$coef)){
      zvalues[,j] = fGarch::psstd(GARCH[[j]][[i]]@fit$z, nu=GARCH[[j]][[i]]@fit$coef[8],
                                  xi = GARCH[[j]][[i]]@fit$coef[7])
    }else{
      zvalues[,j] = fGarch::pstd(GARCH[[j]][[i]]@fit$z, nu=GARCH[[j]][[i]]@fit$coef[7])}
  }
  VineCopula::RVineStructureSelect(data=zvalues,familyset = c(1,2,3,4,5,6), treecrit = test.unweighted$CvM.Rn.mult, indeptest = FALSE)
}
parallel::stopCluster(cl)


test.weighted.best.CvM.tau <- test.weighted.best[["CvM.tau"]]

Vines$CvM.tau.weighted <- vector("list", length = TComplete-WE)
cl <- parallel::makeCluster(cores2)
doParallel::registerDoParallel(cl)
Vines$CvM.tau.weighted <- foreach(i = 1:(TComplete-WE)) %dopar% { #(T-WE)
  set.seed(2021)
  zvalues <- matrix(nrow=WE, ncol=ncol_assets)#nrow_assets
  for(j in 1:ncol_assets){#ncol_assets
    if("skew" %in% names(GARCH[[j]][[i]]@fit$coef)){
      zvalues[,j] = fGarch::psstd(GARCH[[j]][[i]]@fit$z, nu=GARCH[[j]][[i]]@fit$coef[8],
                                  xi = GARCH[[j]][[i]]@fit$coef[7])
    }else{
      zvalues[,j] = fGarch::pstd(GARCH[[j]][[i]]@fit$z, nu=GARCH[[j]][[i]]@fit$coef[7])}
  }
  VineCopula::RVineStructureSelect(data=zvalues,familyset = c(1,2,3,4,5,6), treecrit = test.weighted.best.CvM.tau, indeptest = FALSE)
}
parallel::stopCluster(cl)


test.weighted.best.Rn.tau <- test.weighted.best[["Rn.tau"]]
Vines$Rn.tau.weighted <- vector("list", length = TComplete-WE)
cl <- parallel::makeCluster(cores2)
doParallel::registerDoParallel(cl)
Vines$Rn.tau.weighted <- foreach(i = 1:(TComplete-WE)) %dopar% { #(T-WE)
  set.seed(2021)
  zvalues <- matrix(nrow=WE, ncol=ncol_assets)#nrow_assets
  for(j in 1:ncol_assets){#ncol_assets
    if("skew" %in% names(GARCH[[j]][[i]]@fit$coef)){
      zvalues[,j] = fGarch::psstd(GARCH[[j]][[i]]@fit$z, nu=GARCH[[j]][[i]]@fit$coef[8],
                                  xi = GARCH[[j]][[i]]@fit$coef[7])
    }else{
      zvalues[,j] = fGarch::pstd(GARCH[[j]][[i]]@fit$z, nu=GARCH[[j]][[i]]@fit$coef[7])}
  }
  VineCopula::RVineStructureSelect(data=zvalues,familyset = c(1,2,3,4,5,6), treecrit = test.weighted.best.Rn.tau, indeptest = FALSE)
}
parallel::stopCluster(cl)

test.weighted.best.Rn.CvM <- test.weighted.best[["Rn.CvM"]]
Vines$Rn.CvM.weighted <- vector("list", length = TComplete-WE)
cl <- parallel::makeCluster(cores2)
doParallel::registerDoParallel(cl)
Vines$Rn.CvM.weighted <- foreach(i = 1:(TComplete-WE)) %dopar% { #(T-WE)
  set.seed(2021)
  zvalues <- matrix(nrow=WE, ncol=ncol_assets)#nrow_assets
  for(j in 1:ncol_assets){#ncol_assets
    if("skew" %in% names(GARCH[[j]][[i]]@fit$coef)){
      zvalues[,j] = fGarch::psstd(GARCH[[j]][[i]]@fit$z, nu=GARCH[[j]][[i]]@fit$coef[8],
                                  xi = GARCH[[j]][[i]]@fit$coef[7])
    }else{
      zvalues[,j] = fGarch::pstd(GARCH[[j]][[i]]@fit$z, nu=GARCH[[j]][[i]]@fit$coef[7])}
  }
  VineCopula::RVineStructureSelect(data=zvalues,familyset = c(1,2,3,4,5,6), treecrit = test.weighted.best.Rn.CvM, indeptest = FALSE)
}
parallel::stopCluster(cl)


Vines$Rn.tau.ranked <- vector("list", length = TComplete-WE)
cl <- parallel::makeCluster(cores2)
doParallel::registerDoParallel(cl)
clusterExport(cl=cl, c('RVineStructureSelect_2ranked_alpha', 'initializeFirstGraph_2ranked_alpha','getEdgeInfo_2ranked_alpha',
                       'buildNextGraph_2ranked_alpha','set_treecrit2'))
Vines$Rn.tau.ranked <- foreach(i = 1:(TComplete-WE), .packages = "VineCopula") %dopar% { #(T-WE)
  set.seed(2021)
  zvalues <- matrix(nrow=WE, ncol=ncol_assets)#nrow_assets
  for(j in 1:ncol_assets){#ncol_assets
    if("skew" %in% names(GARCH[[j]][[i]]@fit$coef)){
      zvalues[,j] = fGarch::psstd(GARCH[[j]][[i]]@fit$z, nu=GARCH[[j]][[i]]@fit$coef[8],
                                  xi = GARCH[[j]][[i]]@fit$coef[7])
    }else{
      zvalues[,j] = fGarch::pstd(GARCH[[j]][[i]]@fit$z, nu=GARCH[[j]][[i]]@fit$coef[7])}
  }
  RVineStructureSelect_2ranked_alpha(alpha = best.weights.ranked$Rn.tau, data=zvalues,familyset = c(1,2,3,4,5,6),
                                     treecrit = test.unweighted$Rn.test, treecrit2 = "tau", indeptest = FALSE)
}
parallel::stopCluster(cl)

Vines$CvM.tau.ranked <- vector("list", length = TComplete-WE)
cl <- parallel::makeCluster(cores2)
doParallel::registerDoParallel(cl)
clusterExport(cl=cl, c('RVineStructureSelect_2ranked_alpha', 'initializeFirstGraph_2ranked_alpha','getEdgeInfo_2ranked_alpha',
                       'buildNextGraph_2ranked_alpha','set_treecrit2'))
Vines$CvM.tau.ranked <- foreach(i = 1:(TComplete-WE),.packages = "VineCopula") %dopar% { #(T-WE)
  set.seed(2021)
  zvalues <- matrix(nrow=WE, ncol=ncol_assets)#nrow_assets
  for(j in 1:ncol_assets){#ncol_assets
    if("skew" %in% names(GARCH[[j]][[i]]@fit$coef)){
      zvalues[,j] = fGarch::psstd(GARCH[[j]][[i]]@fit$z, nu=GARCH[[j]][[i]]@fit$coef[8],
                                  xi = GARCH[[j]][[i]]@fit$coef[7])
    }else{
      zvalues[,j] = fGarch::pstd(GARCH[[j]][[i]]@fit$z, nu=GARCH[[j]][[i]]@fit$coef[7])}
  }
  RVineStructureSelect_2ranked_alpha(alpha = best.weights.ranked$CvM.tau, data=zvalues,familyset = c(1,2,3,4,5,6),
                                     treecrit = test.unweighted$CvM.test, treecrit2 = "tau", indeptest = FALSE)
}
parallel::stopCluster(cl)

Vines$Rn.AIC.ranked <- vector("list", length = TComplete-WE)
cl <- parallel::makeCluster(cores2)
doParallel::registerDoParallel(cl)
clusterExport(cl=cl, c('RVineStructureSelect_2ranked_alpha', 'initializeFirstGraph_2ranked_alpha','getEdgeInfo_2ranked_alpha',
                       'buildNextGraph_2ranked_alpha','set_treecrit2'))
Vines$Rn.AIC.ranked <- foreach(i = 1:(TComplete-WE), .packages = "VineCopula") %dopar% { #(T-WE)
  set.seed(2021)
  zvalues <- matrix(nrow=WE, ncol=ncol_assets)#nrow_assets
  for(j in 1:ncol_assets){#ncol_assets
    if("skew" %in% names(GARCH[[j]][[i]]@fit$coef)){
      zvalues[,j] = fGarch::psstd(GARCH[[j]][[i]]@fit$z, nu=GARCH[[j]][[i]]@fit$coef[8],
                                  xi = GARCH[[j]][[i]]@fit$coef[7])
    }else{
      zvalues[,j] = fGarch::pstd(GARCH[[j]][[i]]@fit$z, nu=GARCH[[j]][[i]]@fit$coef[7])}
  }
  RVineStructureSelect_2ranked_alpha(alpha = best.weights.ranked$Rn.AIC, data=zvalues,familyset = c(1,2,3,4,5,6),
                                     treecrit = test.unweighted$Rn.test, treecrit2 = "AIC", indeptest = FALSE)
}
parallel::stopCluster(cl)

Vines$CvM.AIC.ranked <- vector("list", length = TComplete-WE)
cl <- parallel::makeCluster(cores2)
doParallel::registerDoParallel(cl)
clusterExport(cl=cl, c('RVineStructureSelect_2ranked_alpha', 'initializeFirstGraph_2ranked_alpha','getEdgeInfo_2ranked_alpha',
                       'buildNextGraph_2ranked_alpha','set_treecrit2'))
Vines$CvM.AIC.ranked <- foreach(i = 1:(TComplete-WE), .packages = "VineCopula") %dopar% { #(T-WE)
  set.seed(2021)
  zvalues <- matrix(nrow=WE, ncol=ncol_assets)#nrow_assets
  for(j in 1:ncol_assets){#ncol_assets
    if("skew" %in% names(GARCH[[j]][[i]]@fit$coef)){
      zvalues[,j] = fGarch::psstd(GARCH[[j]][[i]]@fit$z, nu=GARCH[[j]][[i]]@fit$coef[8],
                                  xi = GARCH[[j]][[i]]@fit$coef[7])
    }else{
      zvalues[,j] = fGarch::pstd(GARCH[[j]][[i]]@fit$z, nu=GARCH[[j]][[i]]@fit$coef[7])}
  }
  RVineStructureSelect_2ranked_alpha(alpha = best.weights.ranked$CvM.AIC, data=zvalues,familyset = c(1,2,3,4,5,6),
                                     treecrit = test.unweighted$CvM.test, treecrit2 = "AIC", indeptest = FALSE)
}
parallel::stopCluster(cl)

Vines$Rn.CvM.ranked <- vector("list", length = TComplete-WE)
cl <- parallel::makeCluster(cores2)
doParallel::registerDoParallel(cl)
clusterExport(cl=cl, c('RVineStructureSelect_2ranked_alpha', 'initializeFirstGraph_2ranked_alpha','getEdgeInfo_2ranked_alpha',
                       'buildNextGraph_2ranked_alpha','set_treecrit2'))
Vines$Rn.CvM.ranked <- foreach(i = 1:(TComplete-WE),.packages = "VineCopula") %dopar% { #(T-WE)
  set.seed(2021)
  zvalues <- matrix(nrow=WE, ncol=ncol_assets)#nrow_assets
  for(j in 1:ncol_assets){#ncol_assets
    if("skew" %in% names(GARCH[[j]][[i]]@fit$coef)){
      zvalues[,j] = fGarch::psstd(GARCH[[j]][[i]]@fit$z, nu=GARCH[[j]][[i]]@fit$coef[8],
                                  xi = GARCH[[j]][[i]]@fit$coef[7])
    }else{
      zvalues[,j] = fGarch::pstd(GARCH[[j]][[i]]@fit$z, nu=GARCH[[j]][[i]]@fit$coef[7])}
  }
  RVineStructureSelect_2ranked_alpha(alpha = best.weights.ranked$Rn.CvM, data=zvalues,familyset = c(1,2,3,4,5,6),
                                     treecrit = test.unweighted$Rn.test, treecrit2 = test.unweighted$Rn.test, indeptest = FALSE)
}
parallel::stopCluster(cl)


# 9. VaR ------------------------------------------------------------------

# calculate VaR for all tree selection criteria

VaR <- list()

VaR$tau <- vector("list", length = TComplete-WE)
cl <- parallel::makeCluster(cores2)
doParallel::registerDoParallel(cl)
VaR$tau <- foreach(i = 1:(TComplete-WE),.combine='rbind') %dopar% { #(T-WE)
  set.seed(2021)
  sim.rv <- VineCopula::RVineSim(N=N, RVM = Vines$tau[[i]])

  simZ <- matrix(nrow=N, ncol=ncol_assets)
  for(j in 1:ncol_assets){
    if("skew" %in% names(GARCH[[j]][[i]]@fit$coef)){
      simZ[,j] = rugarch::qdist(distribution = "sstd", p = sim.rv[,j], mu = 0, sigma= 1, shape=GARCH[[j]][[i]]@fit$coef[8], skew = GARCH[[j]][[i]]@fit$coef[7])
    }else{
      simZ[,j] = rugarch::qdist(distribution = "std", p = sim.rv[,j], mu = 0, sigma= 1, shape=GARCH[[j]][[i]]@fit$coef[7])
    }
  }

  sigma <- vector(length = ncol_assets)
  for(j in 1:ncol_assets){
    sigma[j] = as.numeric(sqrt(pmax(0,(as.numeric(GARCH[[j]][[i]]@fit$coef[4]+GARCH[[j]][[i]]@fit$coef[2]*GARCH[[j]][[i]]@fit$sigma[WE]^2*GARCH[[j]][[i]]@fit$z[WE]^2+GARCH[[j]][[i]]@fit$coef[6]*GARCH[[j]][[i]]@fit$sigma[WE]^2)))))
  }

  returns <- matrix(nrow = N, ncol = ncol_assets)
  for(j in 1:ncol_assets){
    returns[,j] = as.numeric(GARCH[[j]][[i]]@fit$coef[1] + GARCH[[j]][[i]]@fit$coef[2]*logreturns[(WE-1+i),j] +simZ[,j]*sigma[j] + GARCH[[j]][[i]]@fit$coef[3]*GARCH[[j]][[i]]@fit$sigma[WE]*GARCH[[j]][[i]]@fit$z[WE])
  }

  returnssorted <- apply(returns, 2, sort)

  VaR01 <- mean(returnssorted[N*0.01,])
  VaR025 <- mean(returnssorted[N*0.025,])
  VaR05 <- mean(returnssorted[N*0.05,])

  c(VaR01, VaR025, VaR05)


}
parallel::stopCluster(cl)

VaR$AIC <- vector("list", length = TComplete-WE)
cl <- parallel::makeCluster(cores2)
doParallel::registerDoParallel(cl)
VaR$AIC <- foreach(i = 1:(TComplete-WE),.combine='rbind') %dopar% { #(T-WE)
  set.seed(2021)
  sim.rv <- Vi <- eCopula::RVineSim(N=N, RVM = Vines$AIC[[i]])

  simZ <- matrix(nrow=N, ncol=ncol_assets)
  for(j in 1:ncol_assets){
    if("skew" %in% names(GARCH[[j]][[i]]@fit$coef)){
      simZ[,j] = rugarch::qdist(distribution = "sstd", p = sim.rv[,j], mu = 0, sigma= 1, shape=GARCH[[j]][[i]]@fit$coef[8], skew = GARCH[[j]][[i]]@fit$coef[7])
    }else{
      simZ[,j] = rugarch::qdist(distribution = "std", p = sim.rv[,j], mu = 0, sigma= 1, shape=GARCH[[j]][[i]]@fit$coef[7])
    }
  }

  sigma <- vector(length = ncol_assets)
  for(j in 1:ncol_assets){
    sigma[j] = as.numeric(sqrt(pmax(0,(as.numeric(GARCH[[j]][[i]]@fit$coef[4]+GARCH[[j]][[i]]@fit$coef[2]*GARCH[[j]][[i]]@fit$sigma[WE]^2*GARCH[[j]][[i]]@fit$z[WE]^2+GARCH[[j]][[i]]@fit$coef[6]*GARCH[[j]][[i]]@fit$sigma[WE]^2)))))
  }

  returns <- matrix(nrow = N, ncol = ncol_assets)
  for(j in 1:ncol_assets){
    returns[,j] = as.numeric(GARCH[[j]][[i]]@fit$coef[1] + GARCH[[j]][[i]]@fit$coef[2]*logreturns[(WE-1+i),j] +simZ[,j]*sigma[j] + GARCH[[j]][[i]]@fit$coef[3]*GARCH[[j]][[i]]@fit$sigma[WE]*GARCH[[j]][[i]]@fit$z[WE])
  }

  returnssorted <- apply(returns, 2, sort)

  VaR01 <- mean(returnssorted[N*0.01,])
  VaR025 <- mean(returnssorted[N*0.025,])
  VaR05 <- mean(returnssorted[N*0.05,])

  c(VaR01, VaR025, VaR05)


}
parallel::stopCluster(cl)

VaR$Rn <- vector("list", length = TComplete-WE)
cl <- parallel::makeCluster(cores2)
doParallel::registerDoParallel(cl)
VaR[["Rn"]] <- foreach(i = 1:(TComplete-WE),.combine='rbind') %dopar% { #(T-WE)
  set.seed(2021)
  sim.rv <- VineCopula::RVineSim(N=N, RVM = Vines$Rn[[i]])

  simZ <- matrix(nrow=N, ncol=ncol_assets)
  for(j in 1:ncol_assets){
    if("skew" %in% names(GARCH[[j]][[i]]@fit$coef)){
      simZ[,j] = rugarch::qdist(distribution = "sstd", p = sim.rv[,j], mu = 0, sigma= 1, shape=GARCH[[j]][[i]]@fit$coef[8], skew = GARCH[[j]][[i]]@fit$coef[7])
    }else{
      simZ[,j] = rugarch::qdist(distribution = "std", p = sim.rv[,j], mu = 0, sigma= 1, shape=GARCH[[j]][[i]]@fit$coef[7])
    }
  }

  sigma <- vector(length = ncol_assets)
  for(j in 1:ncol_assets){
    sigma[j] = as.numeric(sqrt(pmax(0,(as.numeric(GARCH[[j]][[i]]@fit$coef[4]+GARCH[[j]][[i]]@fit$coef[2]*GARCH[[j]][[i]]@fit$sigma[WE]^2*GARCH[[j]][[i]]@fit$z[WE]^2+GARCH[[j]][[i]]@fit$coef[6]*GARCH[[j]][[i]]@fit$sigma[WE]^2)))))
  }

  returns <- matrix(nrow = N, ncol = ncol_assets)
  for(j in 1:ncol_assets){
    returns[,j] = as.numeric(GARCH[[j]][[i]]@fit$coef[1] + GARCH[[j]][[i]]@fit$coef[2]*logreturns[(WE-1+i),j] +simZ[,j]*sigma[j] + GARCH[[j]][[i]]@fit$coef[3]*GARCH[[j]][[i]]@fit$sigma[WE]*GARCH[[j]][[i]]@fit$z[WE])
  }

  returnssorted <- apply(returns, 2, sort)

  VaR01 <- mean(returnssorted[N*0.01,])
  VaR025 <- mean(returnssorted[N*0.025,])
  VaR05 <- mean(returnssorted[N*0.05,])

  c(VaR01, VaR025, VaR05)


}
parallel::stopCluster(cl)


VaR$CvM <- vector("list", length = TComplete-WE)
cl <- parallel::makeCluster(cores2)
doParallel::registerDoParallel(cl)
VaR$CvM <- foreach(i = 1:(TComplete-WE),.combine='rbind') %dopar% { #(T-WE)
  set.seed(2021)
  sim.rv <- VineCopula::RVineSim(N=N, RVM = Vines$CvM[[i]])

  simZ <- matrix(nrow=N, ncol=ncol_assets)
  for(j in 1:ncol_assets){
    if("skew" %in% names(GARCH[[j]][[i]]@fit$coef)){
      simZ[,j] = rugarch::qdist(distribution = "sstd", p = sim.rv[,j], mu = 0, sigma= 1, shape=GARCH[[j]][[i]]@fit$coef[8], skew = GARCH[[j]][[i]]@fit$coef[7])
    }else{
      simZ[,j] = rugarch::qdist(distribution = "std", p = sim.rv[,j], mu = 0, sigma= 1, shape=GARCH[[j]][[i]]@fit$coef[7])
    }
  }

  sigma <- vector(length = ncol_assets)
  for(j in 1:ncol_assets){
    sigma[j] = as.numeric(sqrt(pmax(0,(as.numeric(GARCH[[j]][[i]]@fit$coef[4]+GARCH[[j]][[i]]@fit$coef[2]*GARCH[[j]][[i]]@fit$sigma[WE]^2*GARCH[[j]][[i]]@fit$z[WE]^2+GARCH[[j]][[i]]@fit$coef[6]*GARCH[[j]][[i]]@fit$sigma[WE]^2)))))
  }

  returns <- matrix(nrow = N, ncol = ncol_assets)
  for(j in 1:ncol_assets){
    returns[,j] = as.numeric(GARCH[[j]][[i]]@fit$coef[1] + GARCH[[j]][[i]]@fit$coef[2]*logreturns[(WE-1+i),j] +simZ[,j]*sigma[j] + GARCH[[j]][[i]]@fit$coef[3]*GARCH[[j]][[i]]@fit$sigma[WE]*GARCH[[j]][[i]]@fit$z[WE])
  }

  returnssorted <- apply(returns, 2, sort)

  VaR01 <- mean(returnssorted[N*0.01,])
  VaR025 <- mean(returnssorted[N*0.025,])
  VaR05 <- mean(returnssorted[N*0.05,])

  c(VaR01, VaR025, VaR05)


}
parallel::stopCluster(cl)

VaR$Rn.tau.mult <- vector("list", length = TComplete-WE)
cl <- parallel::makeCluster(cores2)
doParallel::registerDoParallel(cl)
VaR$Rn.tau.mult <- foreach(i = 1:(TComplete-WE),.combine='rbind') %dopar% { #(T-WE)
  set.seed(2021)
  sim.rv <- VineCopula::RVineSim(N=N, RVM = Vines$Rn.tau.mult[[i]])

  simZ <- matrix(nrow=N, ncol=ncol_assets)
  for(j in 1:ncol_assets){
    if("skew" %in% names(GARCH[[j]][[i]]@fit$coef)){
      simZ[,j] = rugarch::qdist(distribution = "sstd", p = sim.rv[,j], mu = 0, sigma= 1, shape=GARCH[[j]][[i]]@fit$coef[8], skew = GARCH[[j]][[i]]@fit$coef[7])
    }else{
      simZ[,j] = rugarch::qdist(distribution = "std", p = sim.rv[,j], mu = 0, sigma= 1, shape=GARCH[[j]][[i]]@fit$coef[7])
    }
  }

  sigma <- vector(length = ncol_assets)
  for(j in 1:ncol_assets){
    sigma[j] = as.numeric(sqrt(pmax(0,(as.numeric(GARCH[[j]][[i]]@fit$coef[4]+GARCH[[j]][[i]]@fit$coef[2]*GARCH[[j]][[i]]@fit$sigma[WE]^2*GARCH[[j]][[i]]@fit$z[WE]^2+GARCH[[j]][[i]]@fit$coef[6]*GARCH[[j]][[i]]@fit$sigma[WE]^2)))))
  }

  returns <- matrix(nrow = N, ncol = ncol_assets)
  for(j in 1:ncol_assets){
    returns[,j] = as.numeric(GARCH[[j]][[i]]@fit$coef[1] + GARCH[[j]][[i]]@fit$coef[2]*logreturns[(WE-1+i),j] +simZ[,j]*sigma[j] + GARCH[[j]][[i]]@fit$coef[3]*GARCH[[j]][[i]]@fit$sigma[WE]*GARCH[[j]][[i]]@fit$z[WE])
  }

  returnssorted <- apply(returns, 2, sort)

  VaR01 <- mean(returnssorted[N*0.01,])
  VaR025 <- mean(returnssorted[N*0.025,])
  VaR05 <- mean(returnssorted[N*0.05,])

  c(VaR01, VaR025, VaR05)


}
parallel::stopCluster(cl)

VaR$CvM.tau.mult <- vector("list", length = TComplete-WE)
cl <- parallel::makeCluster(cores2)
doParallel::registerDoParallel(cl)
VaR$CvM.tau.mult <- foreach(i = 1:(TComplete-WE),.combine='rbind') %dopar% { #(T-WE)
  set.seed(2021)
  sim.rv <- VineCopula::RVineSim(N=N, RVM = Vines$CvM.tau.mult[[i]])

  simZ <- matrix(nrow=N, ncol=ncol_assets)
  for(j in 1:ncol_assets){
    if("skew" %in% names(GARCH[[j]][[i]]@fit$coef)){
      simZ[,j] = rugarch::qdist(distribution = "sstd", p = sim.rv[,j], mu = 0, sigma= 1, shape=GARCH[[j]][[i]]@fit$coef[8], skew = GARCH[[j]][[i]]@fit$coef[7])
    }else{
      simZ[,j] = rugarch::qdist(distribution = "std", p = sim.rv[,j], mu = 0, sigma= 1, shape=GARCH[[j]][[i]]@fit$coef[7])
    }
  }

  sigma <- vector(length = ncol_assets)
  for(j in 1:ncol_assets){
    sigma[j] = as.numeric(sqrt(pmax(0,(as.numeric(GARCH[[j]][[i]]@fit$coef[4]+GARCH[[j]][[i]]@fit$coef[2]*GARCH[[j]][[i]]@fit$sigma[WE]^2*GARCH[[j]][[i]]@fit$z[WE]^2+GARCH[[j]][[i]]@fit$coef[6]*GARCH[[j]][[i]]@fit$sigma[WE]^2)))))
  }

  returns <- matrix(nrow = N, ncol = ncol_assets)
  for(j in 1:ncol_assets){
    returns[,j] = as.numeric(GARCH[[j]][[i]]@fit$coef[1] + GARCH[[j]][[i]]@fit$coef[2]*logreturns[(WE-1+i),j] +simZ[,j]*sigma[j] + GARCH[[j]][[i]]@fit$coef[3]*GARCH[[j]][[i]]@fit$sigma[WE]*GARCH[[j]][[i]]@fit$z[WE])
  }

  returnssorted <- apply(returns, 2, sort)

  VaR01 <- mean(returnssorted[N*0.01,])
  VaR025 <- mean(returnssorted[N*0.025,])
  VaR05 <- mean(returnssorted[N*0.05,])

  c(VaR01, VaR025, VaR05)


}
parallel::stopCluster(cl)

VaR$CvM.Rn.mult <- vector("list", length = TComplete-WE)
cl <- parallel::makeCluster(cores2)
doParallel::registerDoParallel(cl)
VaR$CvM.Rn.mult <- foreach(i = 1:(TComplete-WE),.combine='rbind') %dopar% { #(T-WE)
  set.seed(2021)
  sim.rv <- VineCopula::RVineSim(N=N, RVM = Vines$CvM.Rn.mult[[i]])

  simZ <- matrix(nrow=N, ncol=ncol_assets)
  for(j in 1:ncol_assets){
    if("skew" %in% names(GARCH[[j]][[i]]@fit$coef)){
      simZ[,j] = rugarch::qdist(distribution = "sstd", p = sim.rv[,j], mu = 0, sigma= 1, shape=GARCH[[j]][[i]]@fit$coef[8], skew = GARCH[[j]][[i]]@fit$coef[7])
    }else{
      simZ[,j] = rugarch::qdist(distribution = "std", p = sim.rv[,j], mu = 0, sigma= 1, shape=GARCH[[j]][[i]]@fit$coef[7])
    }
  }

  sigma <- vector(length = ncol_assets)
  for(j in 1:ncol_assets){
    sigma[j] = as.numeric(sqrt(pmax(0,(as.numeric(GARCH[[j]][[i]]@fit$coef[4]+GARCH[[j]][[i]]@fit$coef[2]*GARCH[[j]][[i]]@fit$sigma[WE]^2*GARCH[[j]][[i]]@fit$z[WE]^2+GARCH[[j]][[i]]@fit$coef[6]*GARCH[[j]][[i]]@fit$sigma[WE]^2)))))
  }

  returns <- matrix(nrow = N, ncol = ncol_assets)
  for(j in 1:ncol_assets){
    returns[,j] = as.numeric(GARCH[[j]][[i]]@fit$coef[1] + GARCH[[j]][[i]]@fit$coef[2]*logreturns[(WE-1+i),j] +simZ[,j]*sigma[j] + GARCH[[j]][[i]]@fit$coef[3]*GARCH[[j]][[i]]@fit$sigma[WE]*GARCH[[j]][[i]]@fit$z[WE])
  }

  returnssorted <- apply(returns, 2, sort)

  VaR01 <- mean(returnssorted[N*0.01,])
  VaR025 <- mean(returnssorted[N*0.025,])
  VaR05 <- mean(returnssorted[N*0.05,])

  c(VaR01, VaR025, VaR05)


}
parallel::stopCluster(cl)

VaR$CvM.tau.weighted <- vector("list", length = TComplete-WE)
cl <- parallel::makeCluster(cores2)
doParallel::registerDoParallel(cl)
VaR$CvM.tau.weighted <- foreach(i = 1:(TComplete-WE),.combine='rbind') %dopar% { #(T-WE)
  set.seed(2021)
  sim.rv <- VineCopula::RVineSim(N=N, RVM = Vines$CvM.tau.weighted[[i]])

  simZ <- matrix(nrow=N, ncol=ncol_assets)
  for(j in 1:ncol_assets){
    if("skew" %in% names(GARCH[[j]][[i]]@fit$coef)){
      simZ[,j] = rugarch::qdist(distribution = "sstd", p = sim.rv[,j], mu = 0, sigma= 1, shape=GARCH[[j]][[i]]@fit$coef[8], skew = GARCH[[j]][[i]]@fit$coef[7])
    }else{
      simZ[,j] = rugarch::qdist(distribution = "std", p = sim.rv[,j], mu = 0, sigma= 1, shape=GARCH[[j]][[i]]@fit$coef[7])
    }
  }

  sigma <- vector(length = ncol_assets)
  for(j in 1:ncol_assets){
    sigma[j] = as.numeric(sqrt(pmax(0,(as.numeric(GARCH[[j]][[i]]@fit$coef[4]+GARCH[[j]][[i]]@fit$coef[2]*GARCH[[j]][[i]]@fit$sigma[WE]^2*GARCH[[j]][[i]]@fit$z[WE]^2+GARCH[[j]][[i]]@fit$coef[6]*GARCH[[j]][[i]]@fit$sigma[WE]^2)))))
  }

  returns <- matrix(nrow = N, ncol = ncol_assets)
  for(j in 1:ncol_assets){
    returns[,j] = as.numeric(GARCH[[j]][[i]]@fit$coef[1] + GARCH[[j]][[i]]@fit$coef[2]*logreturns[(WE-1+i),j] +simZ[,j]*sigma[j] + GARCH[[j]][[i]]@fit$coef[3]*GARCH[[j]][[i]]@fit$sigma[WE]*GARCH[[j]][[i]]@fit$z[WE])
  }

  returnssorted <- apply(returns, 2, sort)

  VaR01 <- mean(returnssorted[N*0.01,])
  VaR025 <- mean(returnssorted[N*0.025,])
  VaR05 <- mean(returnssorted[N*0.05,])

  c(VaR01, VaR025, VaR05)


}
parallel::stopCluster(cl)

VaR$Rn.tau.weighted <- vector("list", length = TComplete-WE)
cl <- parallel::makeCluster(cores2)
doParallel::registerDoParallel(cl)
VaR$Rn.tau.weighted <- foreach(i = 1:(TComplete-WE),.combine='rbind') %dopar% { #(T-WE)
  set.seed(2021)
  sim.rv <- VineCopula::RVineSim(N=N, RVM = Vines$Rn.tau.weighted[[i]])

  simZ <- matrix(nrow=N, ncol=ncol_assets)
  for(j in 1:ncol_assets){
    if("skew" %in% names(GARCH[[j]][[i]]@fit$coef)){
      simZ[,j] = rugarch::qdist(distribution = "sstd", p = sim.rv[,j], mu = 0, sigma= 1, shape=GARCH[[j]][[i]]@fit$coef[8], skew = GARCH[[j]][[i]]@fit$coef[7])
    }else{
      simZ[,j] = rugarch::qdist(distribution = "std", p = sim.rv[,j], mu = 0, sigma= 1, shape=GARCH[[j]][[i]]@fit$coef[7])
    }
  }

  sigma <- vector(length = ncol_assets)
  for(j in 1:ncol_assets){
    sigma[j] = as.numeric(sqrt(pmax(0,(as.numeric(GARCH[[j]][[i]]@fit$coef[4]+GARCH[[j]][[i]]@fit$coef[2]*GARCH[[j]][[i]]@fit$sigma[WE]^2*GARCH[[j]][[i]]@fit$z[WE]^2+GARCH[[j]][[i]]@fit$coef[6]*GARCH[[j]][[i]]@fit$sigma[WE]^2)))))
  }

  returns <- matrix(nrow = N, ncol = ncol_assets)
  for(j in 1:ncol_assets){
    returns[,j] = as.numeric(GARCH[[j]][[i]]@fit$coef[1] + GARCH[[j]][[i]]@fit$coef[2]*logreturns[(WE-1+i),j] +simZ[,j]*sigma[j] + GARCH[[j]][[i]]@fit$coef[3]*GARCH[[j]][[i]]@fit$sigma[WE]*GARCH[[j]][[i]]@fit$z[WE])
  }

  returnssorted <- apply(returns, 2, sort)

  VaR01 <- mean(returnssorted[N*0.01,])
  VaR025 <- mean(returnssorted[N*0.025,])
  VaR05 <- mean(returnssorted[N*0.05,])

  c(VaR01, VaR025, VaR05)


}
parallel::stopCluster(cl)

VaR$Rn.CvM.weighted <- vector("list", length = TComplete-WE)
cl <- parallel::makeCluster(cores2)
doParallel::registerDoParallel(cl)
VaR$Rn.CvM.weighted <- foreach(i = 1:(TComplete-WE),.combine='rbind') %dopar% { #(T-WE)
  set.seed(2021)
  sim.rv <- VineCopula::RVineSim(N=N, RVM = Vines$Rn.CvM.weighted[[i]])

  simZ <- matrix(nrow=N, ncol=ncol_assets)
  for(j in 1:ncol_assets){
    if("skew" %in% names(GARCH[[j]][[i]]@fit$coef)){
      simZ[,j] = rugarch::qdist(distribution = "sstd", p = sim.rv[,j], mu = 0, sigma= 1, shape=GARCH[[j]][[i]]@fit$coef[8], skew = GARCH[[j]][[i]]@fit$coef[7])
    }else{
      simZ[,j] = rugarch::qdist(distribution = "std", p = sim.rv[,j], mu = 0, sigma= 1, shape=GARCH[[j]][[i]]@fit$coef[7])
    }
  }

  sigma <- vector(length = ncol_assets)
  for(j in 1:ncol_assets){
    sigma[j] = as.numeric(sqrt(pmax(0,(as.numeric(GARCH[[j]][[i]]@fit$coef[4]+GARCH[[j]][[i]]@fit$coef[2]*GARCH[[j]][[i]]@fit$sigma[WE]^2*GARCH[[j]][[i]]@fit$z[WE]^2+GARCH[[j]][[i]]@fit$coef[6]*GARCH[[j]][[i]]@fit$sigma[WE]^2)))))
  }

  returns <- matrix(nrow = N, ncol = ncol_assets)
  for(j in 1:ncol_assets){
    returns[,j] = as.numeric(GARCH[[j]][[i]]@fit$coef[1] + GARCH[[j]][[i]]@fit$coef[2]*logreturns[(WE-1+i),j] +simZ[,j]*sigma[j] + GARCH[[j]][[i]]@fit$coef[3]*GARCH[[j]][[i]]@fit$sigma[WE]*GARCH[[j]][[i]]@fit$z[WE])
  }

  returnssorted <- apply(returns, 2, sort)

  VaR01 <- mean(returnssorted[N*0.01,])
  VaR025 <- mean(returnssorted[N*0.025,])
  VaR05 <- mean(returnssorted[N*0.05,])

  c(VaR01, VaR025, VaR05)


}
parallel::stopCluster(cl)

VaR$Rn.tau.ranked <- vector("list", length = TComplete-WE)
cl <- parallel::makeCluster(cores2)
doParallel::registerDoParallel(cl)
VaR$Rn.tau.ranked <- foreach(i = 1:(TComplete-WE),.combine='rbind') %dopar% { #(T-WE)
  set.seed(2021)
  sim.rv <- VineCopula::RVineSim(N=N, RVM = Vines$Rn.tau.ranked[[i]])

  simZ <- matrix(nrow=N, ncol=ncol_assets)
  for(j in 1:ncol_assets){
    if("skew" %in% names(GARCH[[j]][[i]]@fit$coef)){
      simZ[,j] = rugarch::qdist(distribution = "sstd", p = sim.rv[,j], mu = 0, sigma= 1, shape=GARCH[[j]][[i]]@fit$coef[8], skew = GARCH[[j]][[i]]@fit$coef[7])
    }else{
      simZ[,j] = rugarch::qdist(distribution = "std", p = sim.rv[,j], mu = 0, sigma= 1, shape=GARCH[[j]][[i]]@fit$coef[7])
    }
  }

  sigma <- vector(length = ncol_assets)
  for(j in 1:ncol_assets){
    sigma[j] = as.numeric(sqrt(pmax(0,(as.numeric(GARCH[[j]][[i]]@fit$coef[4]+GARCH[[j]][[i]]@fit$coef[2]*GARCH[[j]][[i]]@fit$sigma[WE]^2*GARCH[[j]][[i]]@fit$z[WE]^2+GARCH[[j]][[i]]@fit$coef[6]*GARCH[[j]][[i]]@fit$sigma[WE]^2)))))
  }

  returns <- matrix(nrow = N, ncol = ncol_assets)
  for(j in 1:ncol_assets){
    returns[,j] = as.numeric(GARCH[[j]][[i]]@fit$coef[1] + GARCH[[j]][[i]]@fit$coef[2]*logreturns[(WE-1+i),j] +simZ[,j]*sigma[j] + GARCH[[j]][[i]]@fit$coef[3]*GARCH[[j]][[i]]@fit$sigma[WE]*GARCH[[j]][[i]]@fit$z[WE])
  }

  returnssorted <- apply(returns, 2, sort)

  VaR01 <- mean(returnssorted[N*0.01,])
  VaR025 <- mean(returnssorted[N*0.025,])
  VaR05 <- mean(returnssorted[N*0.05,])

  c(VaR01, VaR025, VaR05)

}
parallel::stopCluster(cl)

VaR$CvM.tau.ranked <- vector("list", length = TComplete-WE)
cl <- parallel::makeCluster(cores2)
doParallel::registerDoParallel(cl)
VaR$CvM.tau.ranked <- foreach(i = 1:(TComplete-WE),.combine='rbind') %dopar% { #(T-WE)
  set.seed(2021)
  sim.rv <- VineCopula::RVineSim(N=N, RVM = Vines$CvM.tau.ranked[[i]])

  simZ <- matrix(nrow=N, ncol=ncol_assets)
  for(j in 1:ncol_assets){
    if("skew" %in% names(GARCH[[j]][[i]]@fit$coef)){
      simZ[,j] = rugarch::qdist(distribution = "sstd", p = sim.rv[,j], mu = 0, sigma= 1, shape=GARCH[[j]][[i]]@fit$coef[8], skew = GARCH[[j]][[i]]@fit$coef[7])
    }else{
      simZ[,j] = rugarch::qdist(distribution = "std", p = sim.rv[,j], mu = 0, sigma= 1, shape=GARCH[[j]][[i]]@fit$coef[7])
    }
  }

  sigma <- vector(length = ncol_assets)
  for(j in 1:ncol_assets){
    sigma[j] = as.numeric(sqrt(pmax(0,(as.numeric(GARCH[[j]][[i]]@fit$coef[4]+GARCH[[j]][[i]]@fit$coef[2]*GARCH[[j]][[i]]@fit$sigma[WE]^2*GARCH[[j]][[i]]@fit$z[WE]^2+GARCH[[j]][[i]]@fit$coef[6]*GARCH[[j]][[i]]@fit$sigma[WE]^2)))))
  }

  returns <- matrix(nrow = N, ncol = ncol_assets)
  for(j in 1:ncol_assets){
    returns[,j] = as.numeric(GARCH[[j]][[i]]@fit$coef[1] + GARCH[[j]][[i]]@fit$coef[2]*logreturns[(WE-1+i),j] +simZ[,j]*sigma[j] + GARCH[[j]][[i]]@fit$coef[3]*GARCH[[j]][[i]]@fit$sigma[WE]*GARCH[[j]][[i]]@fit$z[WE])
  }

  returnssorted <- apply(returns, 2, sort)

  VaR01 <- mean(returnssorted[N*0.01,])
  VaR025 <- mean(returnssorted[N*0.025,])
  VaR05 <- mean(returnssorted[N*0.05,])

  c(VaR01, VaR025, VaR05)


}
parallel::stopCluster(cl)

VaR$Rn.AIC.ranked <- vector("list", length = TComplete-WE)
cl <- parallel::makeCluster(cores2)
doParallel::registerDoParallel(cl)
VaR$Rn.AIC.ranked <- foreach(i = 1:(TComplete-WE),.combine='rbind') %dopar% { #(T-WE)
  set.seed(2021)
  sim.rv <- VineCopula::RVineSim(N=N, RVM = Vines$Rn.AIC.ranked[[i]])

  simZ <- matrix(nrow=N, ncol=ncol_assets)
  for(j in 1:ncol_assets){
    if("skew" %in% names(GARCH[[j]][[i]]@fit$coef)){
      simZ[,j] = rugarch::qdist(distribution = "sstd", p = sim.rv[,j], mu = 0, sigma= 1, shape=GARCH[[j]][[i]]@fit$coef[8], skew = GARCH[[j]][[i]]@fit$coef[7])
    }else{
      simZ[,j] = rugarch::qdist(distribution = "std", p = sim.rv[,j], mu = 0, sigma= 1, shape=GARCH[[j]][[i]]@fit$coef[7])
    }
  }

  sigma <- vector(length = ncol_assets)
  for(j in 1:ncol_assets){
    sigma[j] = as.numeric(sqrt(pmax(0,(as.numeric(GARCH[[j]][[i]]@fit$coef[4]+GARCH[[j]][[i]]@fit$coef[2]*GARCH[[j]][[i]]@fit$sigma[WE]^2*GARCH[[j]][[i]]@fit$z[WE]^2+GARCH[[j]][[i]]@fit$coef[6]*GARCH[[j]][[i]]@fit$sigma[WE]^2)))))
  }

  returns <- matrix(nrow = N, ncol = ncol_assets)
  for(j in 1:ncol_assets){
    returns[,j] = as.numeric(GARCH[[j]][[i]]@fit$coef[1] + GARCH[[j]][[i]]@fit$coef[2]*logreturns[(WE-1+i),j] +simZ[,j]*sigma[j] + GARCH[[j]][[i]]@fit$coef[3]*GARCH[[j]][[i]]@fit$sigma[WE]*GARCH[[j]][[i]]@fit$z[WE])
  }

  returnssorted <- apply(returns, 2, sort)

  VaR01 <- mean(returnssorted[N*0.01,])
  VaR025 <- mean(returnssorted[N*0.025,])
  VaR05 <- mean(returnssorted[N*0.05,])

  c(VaR01, VaR025, VaR05)

}
parallel::stopCluster(cl)

VaR$CvM.AIC.ranked <- vector("list", length = TComplete-WE)
cl <- parallel::makeCluster(cores2)
doParallel::registerDoParallel(cl)
VaR$CvM.AIC.ranked <- foreach(i = 1:(TComplete-WE),.combine='rbind') %dopar% { #(T-WE)
  set.seed(2021)
  sim.rv <- VineCopula::RVineSim(N=N, RVM = Vines$CvM.AIC.ranked[[i]])

  simZ <- matrix(nrow=N, ncol=ncol_assets)
  for(j in 1:ncol_assets){
    if("skew" %in% names(GARCH[[j]][[i]]@fit$coef)){
      simZ[,j] = rugarch::qdist(distribution = "sstd", p = sim.rv[,j], mu = 0, sigma= 1, shape=GARCH[[j]][[i]]@fit$coef[8], skew = GARCH[[j]][[i]]@fit$coef[7])
    }else{
      simZ[,j] = rugarch::qdist(distribution = "std", p = sim.rv[,j], mu = 0, sigma= 1, shape=GARCH[[j]][[i]]@fit$coef[7])
    }
  }

  sigma <- vector(length = ncol_assets)
  for(j in 1:ncol_assets){
    sigma[j] = as.numeric(sqrt(pmax(0,(as.numeric(GARCH[[j]][[i]]@fit$coef[4]+GARCH[[j]][[i]]@fit$coef[2]*GARCH[[j]][[i]]@fit$sigma[WE]^2*GARCH[[j]][[i]]@fit$z[WE]^2+GARCH[[j]][[i]]@fit$coef[6]*GARCH[[j]][[i]]@fit$sigma[WE]^2)))))
  }

  returns <- matrix(nrow = N, ncol = ncol_assets)
  for(j in 1:ncol_assets){
    returns[,j] = as.numeric(GARCH[[j]][[i]]@fit$coef[1] + GARCH[[j]][[i]]@fit$coef[2]*logreturns[(WE-1+i),j] +simZ[,j]*sigma[j] + GARCH[[j]][[i]]@fit$coef[3]*GARCH[[j]][[i]]@fit$sigma[WE]*GARCH[[j]][[i]]@fit$z[WE])
  }

  returnssorted <- apply(returns, 2, sort)

  VaR01 <- mean(returnssorted[N*0.01,])
  VaR025 <- mean(returnssorted[N*0.025,])
  VaR05 <- mean(returnssorted[N*0.05,])

  c(VaR01, VaR025, VaR05)

}
parallel::stopCluster(cl)


VaR$Rn.CvM.ranked <- vector("list", length = TComplete-WE)
cl <- parallel::makeCluster(cores2)
doParallel::registerDoParallel(cl)
VaR$Rn.CvM.ranked <- foreach(i = 1:(TComplete-WE),.combine='rbind') %dopar% { #(T-WE)
  set.seed(2021)
  sim.rv <- VineCopula::RVineSim(N=N, RVM = Vines$Rn.CvM.ranked[[i]])

  simZ <- matrix(nrow=N, ncol=ncol_assets)
  for(j in 1:ncol_assets){
    if("skew" %in% names(GARCH[[j]][[i]]@fit$coef)){
      simZ[,j] = rugarch::qdist(distribution = "sstd", p = sim.rv[,j], mu = 0, sigma= 1, shape=GARCH[[j]][[i]]@fit$coef[8], skew = GARCH[[j]][[i]]@fit$coef[7])
    }else{
      simZ[,j] = rugarch::qdist(distribution = "std", p = sim.rv[,j], mu = 0, sigma= 1, shape=GARCH[[j]][[i]]@fit$coef[7])
    }
  }

  sigma <- vector(length = ncol_assets)
  for(j in 1:ncol_assets){
    sigma[j] = as.numeric(sqrt(pmax(0,(as.numeric(GARCH[[j]][[i]]@fit$coef[4]+GARCH[[j]][[i]]@fit$coef[2]*GARCH[[j]][[i]]@fit$sigma[WE]^2*GARCH[[j]][[i]]@fit$z[WE]^2+GARCH[[j]][[i]]@fit$coef[6]*GARCH[[j]][[i]]@fit$sigma[WE]^2)))))
  }

  returns <- matrix(nrow = N, ncol = ncol_assets)
  for(j in 1:ncol_assets){
    returns[,j] = as.numeric(GARCH[[j]][[i]]@fit$coef[1] + GARCH[[j]][[i]]@fit$coef[2]*logreturns[(WE-1+i),j] +simZ[,j]*sigma[j] + GARCH[[j]][[i]]@fit$coef[3]*GARCH[[j]][[i]]@fit$sigma[WE]*GARCH[[j]][[i]]@fit$z[WE])
  }

  returnssorted <- apply(returns, 2, sort)

  VaR01 <- mean(returnssorted[N*0.01,])
  VaR025 <- mean(returnssorted[N*0.025,])
  VaR05 <- mean(returnssorted[N*0.05,])

  c(VaR01, VaR025, VaR05)


}
parallel::stopCluster(cl)



# 10. Results -------------------------------------------------------------

#de-select models that combine both tests with each other (not needed, just fitted as a test)
VaRcomplete <- VaR
VaR <- VaRcomplete[c(1:6,8:11,13:14)]


report.testresult <- function(VaRs_){ # collect results
  VaR.01 <- VaRs_[,1]
  VaR.025 <- VaRs_[,2]
  VaR.05 <- VaRs_[,3]
  expected.exceed01 <- rugarch::VaRTest(alpha=0.01,actual=logreturns$PFReturns[(WE+1):(TComplete)], VaR=VaRs_[,1], conf.level=0.95)$expected.exceed
  expected.exceed025 <- rugarch::VaRTest(alpha=0.025,actual=logreturns$PFReturns[(WE+1):(TComplete)], VaR=VaRs_[,2], conf.level=0.95)$expected.exceed
  expected.exceed05 <- rugarch::VaRTest(alpha=0.05,actual=logreturns$PFReturns[(WE+1):(TComplete)], VaR=VaRs_[,3], conf.level=0.95)$expected.exceed

  actual.exceed01 <- rugarch::VaRTest(alpha=0.01,actual=logreturns$PFReturns[(WE+1):(TComplete)], VaR=VaRs_[,1], conf.level=0.95)$actual.exceed
  actual.exceed025 <- rugarch::VaRTest(alpha=0.025,actual=logreturns$PFReturns[(WE+1):(TComplete)], VaR=VaRs_[,2], conf.level=0.95)$actual.exceed
  actual.exceed05 <- rugarch::VaRTest(alpha=0.05,actual=logreturns$PFReturns[(WE+1):(TComplete)], VaR=VaRs_[,3], conf.level=0.95)$actual.exceed

  uc.LRp01 <- rugarch::VaRTest(alpha=0.01,actual=logreturns$PFReturns[(WE+1):(TComplete)], VaR=VaRs_[,1], conf.level=0.95)$uc.LRp
  uc.LRp025 <- rugarch::VaRTest(alpha=0.025,actual=logreturns$PFReturns[(WE+1):(TComplete)], VaR=VaRs_[,2], conf.level=0.95)$uc.LRp
  uc.LRp05 <- rugarch::VaRTest(alpha=0.05,actual=logreturns$PFReturns[(WE+1):(TComplete)], VaR=VaRs_[,3], conf.level=0.95)$uc.LRp
  cc.LRp01 <- rugarch::VaRTest(alpha=0.01,actual=logreturns$PFReturns[(WE+1):(TComplete)], VaR=VaRs_[,1], conf.level=0.95)$cc.LRp
  cc.LRp025 <- rugarch::VaRTest(alpha=0.025,actual=logreturns$PFReturns[(WE+1):(TComplete)], VaR=VaRs_[,2], conf.level=0.95)$cc.LRp
  cc.LRp05 <- rugarch::VaRTest(alpha=0.05,actual=logreturns$PFReturns[(WE+1):(TComplete)], VaR=VaRs_[,3], conf.level=0.95)$cc.LRp
  vd.LRp01 <- rugarch::VaRDurTest(alpha=0.01,actual=logreturns$PFReturns[(WE+1):(TComplete)], VaR=VaRs_[,1], conf.level=0.95)$LRp
  vd.LRp025 <- rugarch::VaRDurTest(alpha=0.025,actual=logreturns$PFReturns[(WE+1):(TComplete)], VaR=VaRs_[,2], conf.level=0.95)$LRp
  vd.LRp05 <- rugarch::VaRDurTest(alpha=0.05,actual=logreturns$PFReturns[(WE+1):(TComplete)], VaR=VaRs_[,3], conf.level=0.95)$LRp

  resultlist <- list()
  resultlist$expected <- c(expected.exceed01,expected.exceed025,expected.exceed05)
  resultlist$actual <- c(actual.exceed01,actual.exceed025,actual.exceed05)
  resultlist$UC <- c(uc.LRp01,uc.LRp025,uc.LRp05)
  resultlist$CC <- c(cc.LRp01,cc.LRp025,cc.LRp05)
  resultlist$WB <- c(vd.LRp01, vd.LRp025, vd.LRp05)



  return(resultlist)
}

testresults <- list()

for(i in 1:length(VaR)){
  testresults[[names(VaR)[i]]] <- report.testresult(VaR[[i]])
}


results <- map_dfr(testresults, ~as.data.frame(.x))
results["weight criterion"] <- rep(names(VaR), each=3)
results["VaR-level (%)"] <- rep(c("1","2.5","5"),length(VaR))

ordernames <- c("tau", "AIC", "CvM", "Rn", "CvM.tau.mult", "Rn.tau.mult",
                "CvM.tau.weighted", "Rn.tau.weighted", "CvM.tau.ranked",
                "Rn.tau.ranked", "CvM.AIC.ranked", "Rn.AIC.ranked")

results <- results%>%arrange(factor(`weight criterion`, levels = ordernames))

results <- results%>%relocate(`weight criterion`, "VaR-level (%)")
results <- results%>%arrange(`VaR-level (%)`)

results.1 <- results%>%dplyr::filter(`VaR-level (%)` == 1)
results.25 <- results%>%dplyr::filter(`VaR-level (%)` == 2.5)
results.5 <- results%>%dplyr::filter(`VaR-level (%)` == 5)

results.overall <- cbind(results.1,results.25,results.5)
colnames(results.overall) <- c("tree crit.01", "VaR-level (%).01", "expected.01",
                               "actual.01", "UC.01", "CC.01", "WB.01",
                               "tree crit.025", "VaR-level (%).025", "expected.025",
                               "actual.025", "UC.025", "CC.025", "WB.025",
                               "tree crit.05", "VaR-level (%).05", "expected.05",
                               "actual.05", "UC.05", "CC.05", "WB.05")

results.overall <- results.overall%>%dplyr::select(-c(2,3,9,10,16,17,8,15))

colnames(results.overall) <- c("tree crit", "actual", "UC", "CC", "WB", "actual", "UC", "CC", "WB", "actual", "UC", "CC", "WB")



results.stargazer <- stargazer(results.overall, summary = FALSE, digits = 3, model.numbers = F, model.names = F,
                               rownames = F, covariate.labels = c("tree crit", "actual",
                                  "$LR_{uc}$", "$LR_{cc}$", "$WB$","actual",
                                  "$LR_{uc}$", "$LR_{cc}$", "$WB$","actual",
                                  "$LR_{uc}$", "$LR_{cc}$", "$WB$"))



latextable <- function(x){
  x <- gsub("tau ", "$w_\\\\tau$ ", x)
  x <- gsub("AIC ", "$w_{AIC}$ ", x)
  x <- gsub("Rn ", "$w_{R_n}$ ", x)
  x <- gsub("CvM ", "$w_{S_n}$ ", x)
  x <- gsub("Rn.tau.mult ", "$w_{R_n \\\\times \\\\tau}$ ", x)
  x <- gsub("CvM.tau.mult ", "$w_{S_n \\\\times \\\\tau}$ ", x)
  x <- gsub("Rn.tau.ranked ", "$w_{R_n, \\\\tau, r}$ ", x)
  x <- gsub("CvM.tau.ranked ", "$w_{S_n, \\\\tau, r}$ ", x)
  x <- gsub("CvM.tau.weighted ", "$w_{S_n, \\\\tau, w}$ ", x)
  x <- gsub("Rn.tau.weighted ", "$w_{R_n, \\\\tau, w}$ ", x)
  x <- gsub("Rn.AIC.ranked ", "$w_{R_n, AIC, r}$ ", x)
  x <- gsub("CvM.AIC.ranked ", "$w_{S_n, AIC, r}$ ", x)

  cat(x)
}

latextable(results.stargazer)





VaR.df <- map(VaR, as.data.frame)

for(i in 1:length(VaR.df)){
  VaR.df[[names(VaR.df[i])]]$name <- names(VaR.df[i])
}

for(i in 1:length(VaR.df)){
  VaR.df[[names(VaR.df[i])]]["date"] <- logreturns$date[(WE+1):TComplete]
}

for(i in 1:length(VaR.df)){
  VaR.df[[names(VaR.df[i])]]["PFReturns"] <- logreturns$PFReturns[(WE+1):TComplete]
  VaR.df[[names(VaR.df[i])]]["PFReturns2"] <- logreturns$PFReturns[(WE+1):TComplete]
}


for(i in 1:length(VaR.df)){
  colnames(VaR.df[[i]]) <- c("VaR-1%", "VaR-2.5%", "VaR-5%","treecrit","date","PFReturns", "Returns")
}


VaR.df <- map_dfr(VaR.df, `[`)


# 11. Plots ---------------------------------------------------------------



# arrange data for plots
VaR.df2 <- pivot_longer(VaR.df, cols = c("VaR-1%", "VaR-2.5%", "VaR-5%", "Returns"), names_to = "VaR")

# add hit indicator
VaR.df2 <- VaR.df2%>%
  mutate(hit = if_else(VaR=="VaR-1%" & value>PFReturns,-0.5,
                        if_else(VaR=="VaR-2.5%" & value > PFReturns,-0.45,
                                if_else(VaR=="VaR-5%" & value > PFReturns,-0.4,
                                        NA_real_))))

alphaplot=0.8
linesize=0.1



VaR.plot <- function(x, title){
  x = x
  VaR.df2%>%
    dplyr::filter(treecrit==x)%>%
    ggplot(data = ., mapping = aes(x=date, y = value, linetype = VaR, colour = VaR))+
    geom_line(alpha = alphaplot, size = linesize)+
    geom_point(mapping = aes(y=hit), shape = 2, show.legend = F)+
    scale_fill_manual(values = c("Returns" = "grey50", "VaR-1%" = "#ece51b","VaR-2.5%" = "#21918c","VaR-5%" = "#440154"), aesthetics = "colour")+ #
    theme(plot.title = element_text(hjust = 0.5))+
    theme_bw()+
    theme(legend.title=element_blank())+
    theme(plot.title = element_text(hjust = 0.5))+
    labs(title = title, y = "Log return \n / VaR-forecast", size = 0.5)+
    guides(linetype = guide_legend(override.aes = list(size = 0.5)))+
    annotate(geom = "text", x = as.Date("2021-11-11", format = "%Y-%m-%d"), y = -0.4, hjust = -0.5, label = VaR.df2%>%dplyr::filter(treecrit==x & hit==-0.4)%>%summarise(n = n()), size = 2.5)+
    annotate(geom = "text", x = as.Date("2021-11-11", format = "%Y-%m-%d"), y = -0.45, hjust = -0.9, label = VaR.df2%>%dplyr::filter(treecrit==x & hit==-0.45)%>%summarise(n = n()), size = 2.5)+
    annotate(geom = "text", x = as.Date("2021-11-11", format = "%Y-%m-%d"), y = -0.5, hjust = -0.9, label = VaR.df2%>%dplyr::filter(treecrit==x & hit==-0.5)%>%summarise(n = n()), size = 2.5)
}

alphaplot=0.8
linesize=0.2

VaR.plots <- list()
VaR.plots$tau <- VaR.plot("tau", TeX("$w_{\\tau}$"))
VaR.plots$AIC <- VaR.plot("AIC", TeX("$w_{AIC}$"))
VaR.plots$Sn <- VaR.plot("CvM", TeX("$w_{S_n}$"))
VaR.plots$Rn <- VaR.plot("Rn", TeX("$w_{R_n}$"))
VaR.plots$Sn.tau.mult <- VaR.plot("CvM.tau.mult", TeX("$w_{S_n \\times \\tau}$"))
VaR.plots$Rn.tau.mult <- VaR.plot("Rn.tau.mult", TeX("$w_{R_n \\times \\tau}$"))
VaR.plots$Sn.tau.weighted <- VaR.plot("CvM.tau.weighted", TeX("$w_{S_n,\\tau,w}$"))
VaR.plots$Rn.tau.weighted <- VaR.plot("Rn.tau.weighted", TeX("$w_{R_n,\\tau,w}$"))
VaR.plots$Sn.tau.ranked <- VaR.plot("CvM.tau.ranked", TeX("$w_{S_n,\\tau,r}$"))
VaR.plots$Rn.tau.ranked <- VaR.plot("Rn.tau.ranked", TeX("$w_{R_n,\\tau,r}$"))
VaR.plots$Sn.AIC.ranked <- VaR.plot("CvM.AIC.ranked", TeX("$w_{S_n,AIC,r}$"))
VaR.plots$Rn.AIC.ranked <- VaR.plot("Rn.AIC.ranked", TeX("$w_{R_n,AIC,r}$"))




g <- grid.arrange(VaR.plots$tau, VaR.plots$AIC, VaR.plots$Sn, VaR.plots$Rn,
                  VaR.plots$Sn.tau.mult, VaR.plots$Rn.tau.mult
                  , ncol=1)
ggsave("VaR.pdf", plot = g, path = "../MA/Plots/VaR", units = "cm", width = 20,
       height = 30)

g <- grid.arrange(VaR.plots$Sn.tau.weighted,VaR.plots$Rn.tau.weighted,
                  VaR.plots$Sn.tau.ranked,VaR.plots$Rn.tau.ranked,
                  VaR.plots$Sn.AIC.ranked, VaR.plots$Rn.AIC.ranked
                  , ncol=1)
ggsave("VaR2.pdf", plot = g, path = "../MA/Plots/VaR", units = "cm", width = 20,
       height = 30)




# exemplary R-Vine
z_first2 <- z_first
colnames(z_first2) <- tickers
tauvine <- RVineStructureSelect(data = z_first2, familyset = c(1,2,3,4,5,6))



pdf(file = "../MA/Plots/VaR/tautree1.pdf", width = 11, height = 10)
plot(tauvine, tree = 1, type = 1, edge.labels = "family-par")
dev.off()

pdf(file = "../MA/Plots/VaR/tautree2.pdf", width = 11, height = 10)
plot(tauvine, tree = 2, type = 1, edge.labels = "family-par")
dev.off()

pdf(file = "../MA/Plots/VaR/tautree3.pdf", width = 11, height = 10)
plot(tauvine, tree = 3, type = 1, edge.labels = "family-par")
dev.off()

pdf(file = "../MA/Plots/VaR/tautree4.pdf", width = 11, height = 10)
plot(tauvine, tree = 4, type = 1, edge.labels = "family-par")
dev.off()

pdf(file = "../MA/Plots/VaR/tautree5.pdf", width = 11, height = 10)
plot(tauvine, tree = 5, type = 1, edge.labels = "family-par")
dev.off()

# vuong test
vuong_self <- function (data, RVM1, RVM2)
{
  args <- preproc(c(as.list(environment()), call = match.call()),
                  check_data, remove_nas, check_if_01, check_nobs, check_RVMs,
                  na.txt = " Only complete observations are used.")
  list2env(args, environment())
  N <- args$n
  Model1.ll <- RVineLogLik(data, RVM1, separate = TRUE, calculate.V = FALSE)$loglik
  Model2.ll <- RVineLogLik(data, RVM2, separate = TRUE, calculate.V = FALSE)$loglik
  anz.1 <- sum(RVM1$family >= 1, na.rm = TRUE) + sum(RVM1$family %in%
                                                       c(2, 7, 8, 9, 10, 17, 18, 19, 20, 27, 28, 29, 30, 37,
                                                         38, 39, 40, 104, 114, 124, 134, 204, 214, 224, 234),
                                                     na.rm = TRUE)
  anz.2 <- sum(RVM2$family >= 1, na.rm = TRUE) + sum(RVM2$family %in%
                                                       c(2, 7, 8, 9, 10, 17, 18, 19, 20, 27, 28, 29, 30, 37,
                                                         38, 39, 40, 104, 114, 124, 134, 204, 214, 224, 234),
                                                     na.rm = TRUE)
  if (all(Model1.ll - Model2.ll == 0)) {
    V <- 0
    V.Schwarz <- 0
    V.Akaike <- 0
    p <- 1
    p.Schwarz <- 1
    p.Akaike <- 1
    model.decision <- "no decision"
    model.decision.Schwarz <- "no decision"
    model.decision.Akaike <- "no decision"
  }
  else {
    LR <- sum(Model1.ll) - sum(Model2.ll)
    LR.Schwarz <- LR - ((anz.1/2 * log(N) - anz.2/2 * log(N)))
    LR.Akaike <- LR - (anz.1 - anz.2)
    w <- sd(Model1.ll - Model2.ll)
    V <- LR/(sqrt(N) * w)
    V.Schwarz <- LR.Schwarz/(sqrt(N) * w)
    V.Akaike <- LR.Akaike/(sqrt(N) * w)
    p <- 2 * pnorm(-abs(V))
    p.Schwarz <- 2 * pnorm(-abs(V.Schwarz))
    p.Akaike <- 2 * pnorm(-abs(V.Akaike))
    if(V > qnorm(1-0.05/2)){
      model.decision <- "choose M1"}
    if (V < -qnorm(1-0.05/2)){
      model.decision <- "choose M2"
    }
    if (abs(V) <= qnorm(1-0.05/2)){
      model.decision <- "no decision"
    }

    if(V.Schwarz > qnorm(1-0.05/2)){
      model.decision.Schwarz <- "choose M1"}
    if (V.Schwarz < -qnorm(1-0.05/2)){
      model.decision.Schwarz <- "choose M2"
    }
    if (abs(V.Schwarz) <= qnorm(1-0.05/2)){
      model.decision.Schwarz <- "no decision"
    }

    if(V.Akaike > qnorm(1-0.05/2)){
      model.decision.Akaike <- "choose M1"}
    if (V.Akaike < -qnorm(1-0.05/2)){
      model.decision.Akaike <- "choose M2"
    }
    if (abs(V.Akaike) <= qnorm(1-0.05/2)){
      model.decision.Akaike <- "no decision"
    }


  }
  return(list(statistic = V, statistic.Akaike = V.Akaike,
              statistic.Schwarz = V.Schwarz, p.value = p, p.value.Akaike = p.Akaike,
              p.value.Schwarz = p.Schwarz, model.decision = model.decision,
              model.decision.Schwarz = model.decision.Schwarz,
              model.decision.Akaike = model.decision.Akaike))
}

environment(vuong_self) <- asNamespace('VineCopula')


vuong.results <- list()

vuong.results$AIC <- matrix(nrow = (TComplete-WE), ncol = 6)
vuong.results$AIC <- as.data.frame(vuong.results$AIC)
colnames(vuong.results$AIC) <- c("p.value", "model.decision", "p.value.Schwarz",
                                "model.decision.Schwarz", "p.value.Akaike",
                                "model.decision.Akaike")

for(i in 1:(TComplete-WE)){#,0
  zvalues <- matrix(nrow=WE, ncol=ncol_assets)#nrow_assets
  for(j in 1:ncol_assets){#ncol_assets
    if("skew" %in% names(GARCH[[j]][[i]]@fit$coef)){
      zvalues[,j] = fGarch::psstd(GARCH[[j]][[i]]@fit$z, nu=GARCH[[j]][[i]]@fit$coef[8],
                                  xi = GARCH[[j]][[i]]@fit$coef[7])
    }else{
      zvalues[,j] = fGarch::pstd(GARCH[[j]][[i]]@fit$z, nu=GARCH[[j]][[i]]@fit$coef[7])}
  }
  vuongtest <-  vuong_self(zvalues, Vines$tau[[i]], Vines$AIC[[i]])
  vuong.results$AIC[i,1] = vuongtest$p.value
  vuong.results$AIC[i,2] = vuongtest$model.decision
  vuong.results$AIC[i,3] = vuongtest$p.value.Schwarz
  vuong.results$AIC[i,4] = vuongtest$model.decision.Schwarz
  vuong.results$AIC[i,5] = vuongtest$p.value.Akaike
  vuong.results$AIC[i,6] = vuongtest$model.decision.Akaike
}



vuong.results$Rn <- matrix(nrow = (TComplete-WE), ncol = 6)
vuong.results$Rn <- as.data.frame(vuong.results$Rn)
colnames(vuong.results$Rn) <- c("p.value", "model.decision", "p.value.Schwarz",
                                "model.decision.Schwarz", "p.value.Akaike",
                                "model.decision.Akaike")

for(i in 1:(TComplete-WE)){#
  zvalues <- matrix(nrow=WE, ncol=ncol_assets)#nrow_assets
  for(j in 1:ncol_assets){#ncol_assets
    if("skew" %in% names(GARCH[[j]][[i]]@fit$coef)){
      zvalues[,j] = fGarch::psstd(GARCH[[j]][[i]]@fit$z, nu=GARCH[[j]][[i]]@fit$coef[8],
                                  xi = GARCH[[j]][[i]]@fit$coef[7])
    }else{
      zvalues[,j] = fGarch::pstd(GARCH[[j]][[i]]@fit$z, nu=GARCH[[j]][[i]]@fit$coef[7])}
  }
   vuongtest <-  vuong_self(zvalues, Vines$tau[[i]], Vines$Rn[[i]])
   vuong.results$Rn[i,1] = vuongtest$p.value
   vuong.results$Rn[i,2] = vuongtest$model.decision
   vuong.results$Rn[i,3] = vuongtest$p.value.Schwarz
   vuong.results$Rn[i,4] = vuongtest$model.decision.Schwarz
   vuong.results$Rn[i,5] = vuongtest$p.value.Akaike
   vuong.results$Rn[i,6] = vuongtest$model.decision.Akaike
}


vuong.results$CvM <- matrix(nrow = (TComplete-WE), ncol = 6)
vuong.results$CvM <- as.data.frame(vuong.results$CvM)
colnames(vuong.results$CvM) <- c("p.value", "model.decision", "p.value.Schwarz",
                                "model.decision.Schwarz", "p.value.Akaike",
                                "model.decision.Akaike")

for(i in 1:(TComplete-WE)){#
  zvalues <- matrix(nrow=WE, ncol=ncol_assets)#nrow_assets
  for(j in 1:ncol_assets){#ncol_assets
    if("skew" %in% names(GARCH[[j]][[i]]@fit$coef)){
      zvalues[,j] = fGarch::psstd(GARCH[[j]][[i]]@fit$z, nu=GARCH[[j]][[i]]@fit$coef[8],
                                  xi = GARCH[[j]][[i]]@fit$coef[7])
    }else{
      zvalues[,j] = fGarch::pstd(GARCH[[j]][[i]]@fit$z, nu=GARCH[[j]][[i]]@fit$coef[7])}
  }
  vuongtest <-  vuong_self(zvalues, Vines$tau[[i]], Vines$CvM[[i]])
  vuong.results$CvM[i,1] = vuongtest$p.value
  vuong.results$CvM[i,2] = vuongtest$model.decision
  vuong.results$CvM[i,3] = vuongtest$p.value.Schwarz
  vuong.results$CvM[i,4] = vuongtest$model.decision.Schwarz
  vuong.results$CvM[i,5] = vuongtest$p.value.Akaike
  vuong.results$CvM[i,6] = vuongtest$model.decision.Akaike
}


vuong.results$Rn.tau.mult <- matrix(nrow = (TComplete-WE), ncol = 6)
vuong.results$Rn.tau.mult <- as.data.frame(vuong.results$Rn.tau.mult)
colnames(vuong.results$Rn.tau.mult) <- c("p.value", "model.decision", "p.value.Schwarz",
                                 "model.decision.Schwarz", "p.value.Akaike",
                                 "model.decision.Akaike")

for(i in 1:(TComplete-WE)){#
  zvalues <- matrix(nrow=WE, ncol=ncol_assets)#nrow_assets
  for(j in 1:ncol_assets){#ncol_assets
    if("skew" %in% names(GARCH[[j]][[i]]@fit$coef)){
      zvalues[,j] = fGarch::psstd(GARCH[[j]][[i]]@fit$z, nu=GARCH[[j]][[i]]@fit$coef[8],
                                  xi = GARCH[[j]][[i]]@fit$coef[7])
    }else{
      zvalues[,j] = fGarch::pstd(GARCH[[j]][[i]]@fit$z, nu=GARCH[[j]][[i]]@fit$coef[7])}
  }
  vuongtest <-  vuong_self(zvalues, Vines$tau[[i]], Vines$Rn.tau.mult[[i]])
  vuong.results$Rn.tau.mult[i,1] = vuongtest$p.value
  vuong.results$Rn.tau.mult[i,2] = vuongtest$model.decision
  vuong.results$Rn.tau.mult[i,3] = vuongtest$p.value.Schwarz
  vuong.results$Rn.tau.mult[i,4] = vuongtest$model.decision.Schwarz
  vuong.results$Rn.tau.mult[i,5] = vuongtest$p.value.Akaike
  vuong.results$Rn.tau.mult[i,6] = vuongtest$model.decision.Akaike
}


vuong.results$CvM.tau.mult <- matrix(nrow = (TComplete-WE), ncol = 6)
vuong.results$CvM.tau.mult <- as.data.frame(vuong.results$CvM.tau.mult)
colnames(vuong.results$CvM.tau.mult) <- c("p.value", "model.decision", "p.value.Schwarz",
                                         "model.decision.Schwarz", "p.value.Akaike",
                                         "model.decision.Akaike")

for(i in 1:(TComplete-WE)){#
  zvalues <- matrix(nrow=WE, ncol=ncol_assets)#nrow_assets
  for(j in 1:ncol_assets){#ncol_assets
    if("skew" %in% names(GARCH[[j]][[i]]@fit$coef)){
      zvalues[,j] = fGarch::psstd(GARCH[[j]][[i]]@fit$z, nu=GARCH[[j]][[i]]@fit$coef[8],
                                  xi = GARCH[[j]][[i]]@fit$coef[7])
    }else{
      zvalues[,j] = fGarch::pstd(GARCH[[j]][[i]]@fit$z, nu=GARCH[[j]][[i]]@fit$coef[7])}
  }
  vuongtest <-  vuong_self(zvalues, Vines$tau[[i]], Vines$CvM.tau.mult[[i]])
  vuong.results$CvM.tau.mult[i,1] = vuongtest$p.value
  vuong.results$CvM.tau.mult[i,2] = vuongtest$model.decision
  vuong.results$CvM.tau.mult[i,3] = vuongtest$p.value.Schwarz
  vuong.results$CvM.tau.mult[i,4] = vuongtest$model.decision.Schwarz
  vuong.results$CvM.tau.mult[i,5] = vuongtest$p.value.Akaike
  vuong.results$CvM.tau.mult[i,6] = vuongtest$model.decision.Akaike
}

vuong.results$Rn.tau.ranked <- matrix(nrow = (TComplete-WE), ncol = 6)
vuong.results$Rn.tau.ranked <- as.data.frame(vuong.results$Rn.tau.ranked)
colnames(vuong.results$Rn.tau.ranked) <- c("p.value", "model.decision", "p.value.Schwarz",
                                          "model.decision.Schwarz", "p.value.Akaike",
                                          "model.decision.Akaike")

for(i in 1:(TComplete-WE)){#
  zvalues <- matrix(nrow=WE, ncol=ncol_assets)#nrow_assets
  for(j in 1:ncol_assets){#ncol_assets
    if("skew" %in% names(GARCH[[j]][[i]]@fit$coef)){
      zvalues[,j] = fGarch::psstd(GARCH[[j]][[i]]@fit$z, nu=GARCH[[j]][[i]]@fit$coef[8],
                                  xi = GARCH[[j]][[i]]@fit$coef[7])
    }else{
      zvalues[,j] = fGarch::pstd(GARCH[[j]][[i]]@fit$z, nu=GARCH[[j]][[i]]@fit$coef[7])}
  }
  vuongtest <-  vuong_self(zvalues, Vines$tau[[i]], Vines$Rn.tau.ranked[[i]])
  vuong.results$Rn.tau.ranked[i,1] = vuongtest$p.value
  vuong.results$Rn.tau.ranked[i,2] = vuongtest$model.decision
  vuong.results$Rn.tau.ranked[i,3] = vuongtest$p.value.Schwarz
  vuong.results$Rn.tau.ranked[i,4] = vuongtest$model.decision.Schwarz
  vuong.results$Rn.tau.ranked[i,5] = vuongtest$p.value.Akaike
  vuong.results$Rn.tau.ranked[i,6] = vuongtest$model.decision.Akaike
}


vuong.results$CvM.tau.ranked <- matrix(nrow = (TComplete-WE), ncol = 6)
vuong.results$CvM.tau.ranked <- as.data.frame(vuong.results$CvM.tau.ranked)
colnames(vuong.results$CvM.tau.ranked) <- c("p.value", "model.decision", "p.value.Schwarz",
                                           "model.decision.Schwarz", "p.value.Akaike",
                                           "model.decision.Akaike")

for(i in 1:(TComplete-WE)){#
  zvalues <- matrix(nrow=WE, ncol=ncol_assets)#nrow_assets
  for(j in 1:ncol_assets){#ncol_assets
    if("skew" %in% names(GARCH[[j]][[i]]@fit$coef)){
      zvalues[,j] = fGarch::psstd(GARCH[[j]][[i]]@fit$z, nu=GARCH[[j]][[i]]@fit$coef[8],
                                  xi = GARCH[[j]][[i]]@fit$coef[7])
    }else{
      zvalues[,j] = fGarch::pstd(GARCH[[j]][[i]]@fit$z, nu=GARCH[[j]][[i]]@fit$coef[7])}
  }
  vuongtest <-  vuong_self(zvalues, Vines$tau[[i]], Vines$CvM.tau.ranked[[i]])
  vuong.results$CvM.tau.ranked[i,1] = vuongtest$p.value
  vuong.results$CvM.tau.ranked[i,2] = vuongtest$model.decision
  vuong.results$CvM.tau.ranked[i,3] = vuongtest$p.value.Schwarz
  vuong.results$CvM.tau.ranked[i,4] = vuongtest$model.decision.Schwarz
  vuong.results$CvM.tau.ranked[i,5] = vuongtest$p.value.Akaike
  vuong.results$CvM.tau.ranked[i,6] = vuongtest$model.decision.Akaike
}


vuong.results$CvM.tau.weighted <- matrix(nrow = (TComplete-WE), ncol = 6)
vuong.results$CvM.tau.weighted <- as.data.frame(vuong.results$CvM.tau.weighted)
colnames(vuong.results$CvM.tau.weighted) <- c("p.value", "model.decision", "p.value.Schwarz",
                                            "model.decision.Schwarz", "p.value.Akaike",
                                            "model.decision.Akaike")

for(i in 1:(TComplete-WE)){#
  zvalues <- matrix(nrow=WE, ncol=ncol_assets)#nrow_assets
  for(j in 1:ncol_assets){#ncol_assets
    if("skew" %in% names(GARCH[[j]][[i]]@fit$coef)){
      zvalues[,j] = fGarch::psstd(GARCH[[j]][[i]]@fit$z, nu=GARCH[[j]][[i]]@fit$coef[8],
                                  xi = GARCH[[j]][[i]]@fit$coef[7])
    }else{
      zvalues[,j] = fGarch::pstd(GARCH[[j]][[i]]@fit$z, nu=GARCH[[j]][[i]]@fit$coef[7])}
  }
  vuongtest <-  vuong_self(zvalues, Vines$tau[[i]], Vines$CvM.tau.weighted[[i]])
  vuong.results$CvM.tau.weighted[i,1] = vuongtest$p.value
  vuong.results$CvM.tau.weighted[i,2] = vuongtest$model.decision
  vuong.results$CvM.tau.weighted[i,3] = vuongtest$p.value.Schwarz
  vuong.results$CvM.tau.weighted[i,4] = vuongtest$model.decision.Schwarz
  vuong.results$CvM.tau.weighted[i,5] = vuongtest$p.value.Akaike
  vuong.results$CvM.tau.weighted[i,6] = vuongtest$model.decision.Akaike
}


vuong.results$Rn.tau.weighted <- matrix(nrow = (TComplete-WE), ncol = 6)
vuong.results$Rn.tau.weighted <- as.data.frame(vuong.results$Rn.tau.weighted)
colnames(vuong.results$Rn.tau.weighted) <- c("p.value", "model.decision", "p.value.Schwarz",
                                              "model.decision.Schwarz", "p.value.Akaike",
                                              "model.decision.Akaike")

for(i in 1:(TComplete-WE)){#
  zvalues <- matrix(nrow=WE, ncol=ncol_assets)#nrow_assets
  for(j in 1:ncol_assets){#ncol_assets
    if("skew" %in% names(GARCH[[j]][[i]]@fit$coef)){
      zvalues[,j] = fGarch::psstd(GARCH[[j]][[i]]@fit$z, nu=GARCH[[j]][[i]]@fit$coef[8],
                                  xi = GARCH[[j]][[i]]@fit$coef[7])
    }else{
      zvalues[,j] = fGarch::pstd(GARCH[[j]][[i]]@fit$z, nu=GARCH[[j]][[i]]@fit$coef[7])}
  }
  vuongtest <-  vuong_self(zvalues, Vines$tau[[i]], Vines$Rn.tau.weighted[[i]])
  vuong.results$Rn.tau.weighted[i,1] = vuongtest$p.value
  vuong.results$Rn.tau.weighted[i,2] = vuongtest$model.decision
  vuong.results$Rn.tau.weighted[i,3] = vuongtest$p.value.Schwarz
  vuong.results$Rn.tau.weighted[i,4] = vuongtest$model.decision.Schwarz
  vuong.results$Rn.tau.weighted[i,5] = vuongtest$p.value.Akaike
  vuong.results$Rn.tau.weighted[i,6] = vuongtest$model.decision.Akaike
}


vuong.results$Rn.AIC.ranked <- matrix(nrow = (TComplete-WE), ncol = 6)
vuong.results$Rn.AIC.ranked <- as.data.frame(vuong.results$Rn.AIC.ranked)
colnames(vuong.results$Rn.AIC.ranked) <- c("p.value", "model.decision", "p.value.Schwarz",
                                             "model.decision.Schwarz", "p.value.Akaike",
                                             "model.decision.Akaike")

for(i in 1:(TComplete-WE)){#
  zvalues <- matrix(nrow=WE, ncol=ncol_assets)#nrow_assets
  for(j in 1:ncol_assets){#ncol_assets
    if("skew" %in% names(GARCH[[j]][[i]]@fit$coef)){
      zvalues[,j] = fGarch::psstd(GARCH[[j]][[i]]@fit$z, nu=GARCH[[j]][[i]]@fit$coef[8],
                                  xi = GARCH[[j]][[i]]@fit$coef[7])
    }else{
      zvalues[,j] = fGarch::pstd(GARCH[[j]][[i]]@fit$z, nu=GARCH[[j]][[i]]@fit$coef[7])}
  }
  vuongtest <-  vuong_self(zvalues, Vines$tau[[i]], Vines$Rn.AIC.ranked[[i]])
  vuong.results$Rn.AIC.ranked[i,1] = vuongtest$p.value
  vuong.results$Rn.AIC.ranked[i,2] = vuongtest$model.decision
  vuong.results$Rn.AIC.ranked[i,3] = vuongtest$p.value.Schwarz
  vuong.results$Rn.AIC.ranked[i,4] = vuongtest$model.decision.Schwarz
  vuong.results$Rn.AIC.ranked[i,5] = vuongtest$p.value.Akaike
  vuong.results$Rn.AIC.ranked[i,6] = vuongtest$model.decision.Akaike
}



vuong.results$CvM.AIC.ranked <- matrix(nrow = (TComplete-WE), ncol = 6)
vuong.results$CvM.AIC.ranked <- as.data.frame(vuong.results$CvM.AIC.ranked)
colnames(vuong.results$CvM.AIC.ranked) <- c("p.value", "model.decision", "p.value.Schwarz",
                                           "model.decision.Schwarz", "p.value.Akaike",
                                           "model.decision.Akaike")

for(i in 1:(TComplete-WE)){#
  zvalues <- matrix(nrow=WE, ncol=ncol_assets)#nrow_assets
  for(j in 1:ncol_assets){#ncol_assets
    if("skew" %in% names(GARCH[[j]][[i]]@fit$coef)){
      zvalues[,j] = fGarch::psstd(GARCH[[j]][[i]]@fit$z, nu=GARCH[[j]][[i]]@fit$coef[8],
                                  xi = GARCH[[j]][[i]]@fit$coef[7])
    }else{
      zvalues[,j] = fGarch::pstd(GARCH[[j]][[i]]@fit$z, nu=GARCH[[j]][[i]]@fit$coef[7])}
  }
  vuongtest <-  vuong_self(zvalues, Vines$tau[[i]], Vines$CvM.AIC.ranked[[i]])
  vuong.results$CvM.AIC.ranked[i,1] = vuongtest$p.value
  vuong.results$CvM.AIC.ranked[i,2] = vuongtest$model.decision
  vuong.results$CvM.AIC.ranked[i,3] = vuongtest$p.value.Schwarz
  vuong.results$CvM.AIC.ranked[i,4] = vuongtest$model.decision.Schwarz
  vuong.results$CvM.AIC.ranked[i,5] = vuongtest$p.value.Akaike
  vuong.results$CvM.AIC.ranked[i,6] = vuongtest$model.decision.Akaike
}




for(i in 1:length(vuong.results)){
  vuong.results[[names(vuong.results)[[i]]]]$name <- names(vuong.results)[[i]]
}



vuong.results.overall <- bind_rows(vuong.results$tau, vuong.results$AIC, vuong.results$CvM,
                                   vuong.results$Rn, vuong.results$CvM.tau.mult,
                                   vuong.results$Rn.tau.mult, vuong.results$CvM.tau.weighted,
                                   vuong.results$Rn.tau.weighted, vuong.results$CvM.tau.ranked,
                                   vuong.results$Rn.tau.ranked, vuong.results$CvM.AIC.ranked,
                                   vuong.results$Rn.AIC.ranked)



vuong.results.df4 <- vuong.results.overall%>%group_by(name)%>%
  dplyr::count(model.decision, .drop = F)%>%pivot_wider(names_from = model.decision, values_from = n, values_fill = 0)%>%
  mutate(`no decision %` = paste0("(",round(`no decision` / 2200 * 100,2), ")"))%>%
  mutate(`choose M1 %` = paste0("(",round(`choose M1` / 2200 * 100,2), ")"))%>%
  mutate(`choose M2 %` = paste0("(",round(`choose M2` / 2200 * 100,2), ")"))
vuong.results.nodecision <- vuong.results.df4%>%mutate(across(where(is.numeric), as.character))%>%
  pivot_longer(cols = c(`no decision`, `no decision %`), names_to = "test", values_to = "no decision")%>%
  dplyr::select(c(`name`, `no decision`))
vuong.results.nodecision

vuong.results.m1 <- vuong.results.df4%>%mutate(across(where(is.numeric), as.character))%>%
  pivot_longer(cols = c(`choose M1`, `choose M1 %`), names_to = "test2", values_to = "choose tau")%>%
  dplyr::select(c(`name`, `choose tau`))

vuong.results.m2 <- vuong.results.df4%>%mutate(across(where(is.numeric), as.character))%>%
  pivot_longer(cols = c(`choose M2`, `choose M2 %`), names_to = "test3", values_to = "choose alternative")%>%
  dplyr::select(c(`name`, `choose alternative`))

vuong.results.nodecision["choose tau"] <- vuong.results.m1$`choose tau`
vuong.results.nodecision["choose alternative"] <- vuong.results.m2$`choose alternative`

ordernames <- c("AIC", "CvM", "Rn", "CvM.tau.mult", "Rn.tau.mult",
                "CvM.tau.weighted", "Rn.tau.weighted", "CvM.tau.ranked",
                "Rn.tau.ranked", "CvM.AIC.ranked", "Rn.AIC.ranked")

vuong.results.nodecision <- vuong.results.nodecision%>%arrange(factor(name, levels = ordernames))

vuong.results.nodecision$name[seq(2,22,2)] <- NA

vuong.results.nodecision.stargazer <- stargazer(vuong.results.nodecision, summary = F, rownames = F, title = "Vuong test")

latextable(vuong.results.nodecision.stargazer)


aicsfirstperiod <- matrix(nrow = 15, ncol = 4)
aicsfirstperiod <- as.data.frame(aicsfirstperiod)
colnames(aicsfirstperiod) <- c("treecrit", "AIC", "BIC", "logLik")
for(i in 1:nrow(aicsfirstperiod)){
  aicsfirstperiod[i,1] <- names(Vines)[[i]]
  aicsfirstperiod[i,2] <- Vines[[names(Vines)[[i]]]][[1]]$AIC
  aicsfirstperiod[i,3] <- Vines[[names(Vines)[[i]]]][[1]]$BIC
  aicsfirstperiod[i,4] <- Vines[[names(Vines)[[i]]]][[1]]$logLik
}

aicsfirstperiod <- aicsfirstperiod%>%dplyr::filter(treecrit != "CvM.Rn.mult"
                                                    & treecrit != "Rn.CvM.weighted"
                                                    & treecrit != "Rn.CvM.ranked")

ordernames <- c("tau", "AIC", "CvM", "Rn", "CvM.tau.mult", "Rn.tau.mult",
                "CvM.tau.weighted", "Rn.tau.weighted", "CvM.tau.ranked",
                "Rn.tau.ranked", "CvM.AIC.ranked", "Rn.AIC.ranked")
aicsfirstperiod <- aicsfirstperiod%>%arrange(factor(treecrit, levels = ordernames))
aicsfirstperiod.stargazer <- stargazer(aicsfirstperiod, summary = F,digits = 2, model.numbers = F, model.names = F,
                                       rownames = F)
latextable(aicsfirstperiod.stargazer)



# packages and versions

usedpackages <- matrix(ncol = 2, nrow = 23)
usedpackages <- as.data.frame(usedpackages)
colnames(usedpackages) <- c("package", "version")

usedpackages$package <- c("VineCopula", "tidyverse", "copula", "ggplot2",
                          "ggplotify", "doParallel", "stargazer",
                          "rugarch", "MASS", "segMGarch", "GAS", "rmgarch",
                          "data.table", "matrixStats", "fGarch", "purrr",
                          "tidyquant", "parallelly", "viridis", "gridExtra",
                          "latex2exp", "tseries", "stringi")
for(i in 1:nrow(usedpackages)){
  usedpackages[i,2] <- as.character(packageVersion(usedpackages[i,1]))
}

usedpackages <- dplyr::arrange(usedpackages, usedpackages$package)
usedpackages <- as.matrix(usedpackages)
stargazer(usedpackages)





