#Setting the working directory
mywd <- 'D:/UNI/Math3624'
setwd(mywd)
getwd()

## ---- DCCest
load('USequity.R')

#Installing packages
#as ccgarch is no longer available 
#we had to use an archive copy
install.packages("ccgarch")
install.packages("quantmod")
install.packages("tseries")
install.packages("ggplot2")
install.packages("tidyverse")
install.packages("GGally")
install.packages("urca")
install.packages("vars")
install.packages('pracma')
install.packages('PerformanceAnalytics')

library(ccgarch)
library(quantmod)
library(tseries)
library(MASS)
library(ggplot2)
library(tidyverse)
library(GGally)
library(urca)
library(pracma)
library(TTR)
library(PerformanceAnalytics)
library(vars)

###Plotting the starting data###
us.sect_plot<-function(){
  #Making us.sect into a data frame
  us_sect_df <- fortify(us.sect)
  
  # sorts the data into a 3 columns table; 
  #date,sector and price
  us_sect_long_df <- us_sect_df %>%
    pivot_longer(cols = 2:5,
                 names_to = "Sector",
                 values_to = "Price")
 
   #Improving names of the sectors by creating a new 
  #column(Sector_new) and renaming the sectors
  us_sect_long_df <- us_sect_long_df %>%
    mutate(Sector_new = factor(Sector, 
                               levels = c("ConsumerServ",
                                          "Financials",
                                          "Healtcare",
                                          "RealEstate"),
                               labels = c("Consumer Services",
                                          "Financial Services",
                                          "Healthcare",
                                          "Real Estate")))
  # Plotting the graph using ggplot
  us_sect_plot<-ggplot(us_sect_long_df,
                       aes(x = Index,
                           y = Price,
                           colour = Sector_new)) +
    geom_line(size=1.25) +
    theme(text = element_text(size=20))+
    scale_color_manual(values=c("Blue","Red","black","orange"))+
    labs(x = "Date",
         y = "Closing Price($)",
         colour = "Sector Type")+
    #Adding the start of Covid line 
    geom_vline(xintercept = as.numeric(
                            us_sect_long_df$Index[3205]), 
               linetype="dotted",color = "violetred", size=1.5)+
    annotate("text",x=c(us_sect_long_df$Index[3250]),y=250,
             label=c("March 11th 2020"),size=6,hjust=0)+
    #Adding First and last closing price for each Sector
    annotate("text",x=c(us_sect_long_df$Index[1]),y=140,
             label=c("145.70"),size=6,hjust=0)+
    annotate("text",x=c(us_sect_long_df$Index[2]),y=83,
             label=c("77.25"),size=6,hjust=0)+
    annotate("text",x=c(us_sect_long_df$Index[3]),y=58,
             label=c("51.115"),size=6,hjust=0)+
    annotate("text",x=c(us_sect_long_df$Index[4]),y=34,
             label=c("37.8825"),size=6,hjust=0)+
    
    annotate("text",x=c(us_sect_long_df$Index[4129]),y=250,
             label=c("255.17"),size=6,hjust=0)+
    annotate("text",x=c(us_sect_long_df$Index[4130]),y=90,
             label=c("88.87"),size=6,hjust=0)+
    annotate("text",x=c(us_sect_long_df$Index[4131]),y=67,
             label=c("69.930"),size=6,hjust=0)+
    annotate("text",x=c(us_sect_long_df$Index[4132]),y=75,
             label=c("73.77"),size=6,hjust=0)
  return(us_sect_plot)
}
us.sect_plot()

###DCC###
#A function that calculates the DCC of the data set
DCC<-function(){
  #Log return vector step 1
  
  muhat<-function(xts){
    logts<-log(xts)
    logtsdiff<-diff(logts)*100
    logtsdiff[is.na(logtsdiff)]<-0  
    mu_hat <- apply(logtsdiff, 2, mean)
    return(mu_hat)
  }
  
  #Function 1 complete
  
  mhut<-muhat(us.sect)
  
  #Conditional Covariance Matrix equation step 2
  
  h_return<-function(xts){
    a<-c(0.004,0.005,0.006,0.007)
    A<-diag(0.2,4)
    B<-diag(0.35,4)
    inip<-c(0.15,0.8)
    rtnval<-dcc.estimation(inia=a,iniA=A,iniB=B,ini.dcc=inip,
                           dvar=xts,model='diagonal')
    hmtx<-rtnval$h
    DCCmtx<-rtnval$DCC
    htail<-sqrt(tail(hmtx,n=1))
    tDCCmtx<-tail(DCCmtx,1)
    R=matrix(tDCCmtx,nrow=4,byrow=T)
    DCCM <- (t(htail)%*%htail)*R
    return(solve(DCCM))
    
  }
  Health<-  us.sect$Healtcare
  #Lowest 158.63
  #peak 221.16
  #Function 2 complete
  Cosumer<-us.sect$ConsumerServ
  #peak59.6725
  #lowest39.7950
  h_return(us.sect)
  
  
  
  #Weights equation Step 3
  
  weights<-function(xts,q){
    mut<-muhat(xts)
    Hinv <-h_return(xts)
    onea <- c(1,1,1,1)
    a<- t(onea) %*% Hinv %*% onea
    b<- t(onea) %*% Hinv %*% mut
    (Hinv %*% onea)/(c(a))+q*(Hinv%*%(mut-c((b/a))))
  }
  
  #function 3 complete
  
  weights(us.sect,1)
  
  #portfolio Returns Function Step 4
  
  Portret<- function(xts,q){
    wi<- weights(xts,q)
    ei<- muhat(xts)
    ret<- c(ei%*%wi)
    return(ret)
    
  }
  Portret(us.sect,1)
  
  
  
  ################
  #===============# Section for long evenly weighted returns. 
  #take the first 300 values to zero
  ################## Take values of max and min returns,
  #against the weighted
  
  logret<-diff(log(us.sect))
  logret$Healtcare[1:300]<-0
  logret$RealEstate[1:300]<-0
  logret$Financials[1:300]<-0
  logret$ConsumerServ[1:300]<-0
  w<-c(0.25,0.25,0.25,0.25)
  ret.long<-cumsum(logret%*%w)
  plot.ts(ret.long)
  
  ###############################
  #=============================# Section for q value 0
  ###############################
  
  portrets.q0<- rep(0,1033)
  
  window=300
  
  
  for(t in(window+1):1033){
    winxts<- us.sect[(t-window):t]
    
    portrets.q0[t]<- Portret(winxts,0)
  } 
  portrets.q0
  port.retq0<-cumsum(portrets.q0)
  plot.ts(port.retq0)
  weight.zero<-weights(us.sect,0)
  ###############################
  #===============================# Section for q value 1
  ###############################
  
  portrets<- rep(0,1033)
  
  window=300
  
  
  for(t in(window+1):1033){
    winxts<- us.sect[(t-window):t]
    
    portrets[t]<- Portret(winxts,1)
  }  
  
  weight.one<-weights(us.sect,1)
  
  port.ret<-cumsum(portrets)
  
  
  
  
  ###############################
  #=========================#    section for q value 2
  #############################
  portrets.q2<- rep(0,1033)
  window=300
  for(t in(window+1):1033){
    winxts<- us.sect[(t-window):t]
    
    portrets.q2[t]<- Portret(winxts,2)
  } 
  weight.two<-weights(us.sect,2)
  port.retq2<-cumsum(portrets.q2)
  
  
  ###############################
  #=============================# q value of 3.
  ###############################
  
  portrets.q10<- rep(0,1033)
  window=300
  for(t in(window+1):1033){
    winxts<- us.sect[(t-window):t]
    
    portrets.q10[t]<- Portret(winxts,10)
  } 
  
  weight.three<-weights(us.sect,10)
  
  port.retq10<-cumsum(portrets.q10)
  
  
  ##############################
  #===============================#  weights given for various 
  #risk appetite parameters 
  ################################ q=0
  
  us.sect.weights<-us.sect
  
  window=300
  for (t in 1:1033){
    
    us.sect.weights$Healtcare[t]<-0
    us.sect.weights$RealEstate[t]<-0
    us.sect.weights$Financials[t]<-0
    us.sect.weights$ConsumerServ[t]<-0
    
    
  }
  us.sect.weights.q0<-us.sect.weights
  
  for (t in (window+1):1033){
    
    winxts<-us.sect[(t-window):t]
    weights.test<-weights(winxts,0)
    us.sect.weights.q0$Healtcare[t]<-weights.test[1]
    us.sect.weights.q0$RealEstate[t]<-weights.test[2]
    us.sect.weights.q0$Financials[t]<-weights.test[3]
    us.sect.weights.q0$ConsumerServ[t]<-weights.test[4]
    
    
  }
  
  ##############################
  #===============================#  weights given for various 
  #risk appetite parameters
  ################################ q=1
  us.sect.weights.q1<-us.sect.weights
  
  for (t in (window+1):1033){
    
    winxts<-us.sect[(t-window):t]
    weights.test<-weights(winxts,1)
    us.sect.weights.q1$Healtcare[t]<-weights.test[1]
    us.sect.weights.q1$RealEstate[t]<-weights.test[2]
    us.sect.weights.q1$Financials[t]<-weights.test[3]
    us.sect.weights.q1$ConsumerServ[t]<-weights.test[4]
    
    
  }

  ###Plotting portfolio returns for various risk appetites###
  #Putting data into data frames
  port.retq0_df<-as.data.frame(port.retq0)
  port.ret_df<-as.data.frame(port.ret)
  port.retq2_df<-as.data.frame(port.retq2)
  port.retq10_df<-as.data.frame(port.retq10)
  us_sect_df <- fortify(us.sect)
  #Deleting the healthcare,RealEstate,Financials and 
  #ConsumerServ columns from us.sect
  # We did this as the dates we had was not everyday 
  #between 03-01-2017 and 09-02-2021 
  # as there is no trading on the weekends 
  date<-select(us_sect_df,-Healtcare,-RealEstate,-Financials,
               -ConsumerServ)
  
  #Putting all of the portfolio returns into a data frame
  q_df<-data.frame(date,port.retq0_df,port.ret_df,
                   port.retq2_df,port.retq10_df)
  
  #Creates 3 columns; date,portfolio_returns 
  #and cumulative_log_returns
  q_long_df <- q_df %>%
    pivot_longer(cols = 2:5,
                 names_to = "portfolio_returns",
                 values_to = "cumulative_log_returns")
  
  #Improving names of the portfolio returns by creating
  #a new column(portfolio_returns_2) and renaming the sectors 
  q_long_df_2<- q_long_df %>%
    mutate(portfolio_returns_2 = factor(portfolio_returns, 
                      levels = c("port.retq0",
                                 "port.ret",
                                 "port.retq2",
                                 "port.retq10"),
                      labels = c("portfolio returns for q=0",
                                 "portfolio returns for q=1",
                                 "portfolio returns for q=2",
                                 "portfolio returns for q=10")))
  
  #Creating the graph
  q_plot<-ggplot(q_long_df_2,
                 aes(x = Index,
                     y = cumulative_log_returns,
                     colour = portfolio_returns_2)) +
    geom_line(size=1.25) +
    theme(text = element_text(size=20))+
    scale_color_manual(values=c("Blue","Red","black","orange"))+
    labs(x = "Date",
         y = "Cumulative log returns",
         colour = "Cumulative log returns")+
    #Adding the start of Covid line
    geom_vline(xintercept = as.numeric(q_long_df_2$Index[3205]),
               linetype="dotted",color = "violetred", size=1.5)+
    annotate("text",x=c(q_long_df_2$Index[3220]),y=70,
             label=c("March 11th 2020"),size=6,hjust=0)+
    #Adding last value for each Sector
    annotate("text",x=c(q_long_df_2$Index[4129]),y=28.33528,
             label=c("28.335"),size=6,hjust=0)+
    annotate("text",x=c(q_long_df_2$Index[4130]),y=33.67388,
             label=c("33.674"),size=6,hjust=0)+
    annotate("text",x=c(q_long_df_2$Index[4131]),y=39.01248,
             label=c("39.012"),size=6,hjust=0)+
    annotate("text",x=c(q_long_df_2$Index[4132]),y=81.72132,
             label=c("81.721"),size=6,hjust=0)

  ###Plotting portfolio returns for various risk appetites
  ###and the long position###
  
  #Putting long returns into a data frame
  #we have data frames already for the 
  #different risk appetites cause of the previous plot 
  
  port.ret_long<-as.data.frame(ret.long)
  #Deleting the healthcare,RealEstate,Financials 
  #and ConsumerServ columns from us.sect
  # We did this as the dates we had was not everyday
  #between 03-01-2017 and 09-02-2021 
  #as there is no trading on the weekends 
  date<-select(us_sect_df,-Healtcare,-RealEstate,-Financials,
               -ConsumerServ) 
  
  #Putting all of the portfolio returns into one data frame
  q_all_df<-data.frame(date,port.retq0,port.ret,
                       port.retq2,port.retq10,port.ret_long)
  
  #Creates 3 columns; date,portfolio_returns
  #and cumulative_log_returns
  q_all_long_df <- q_all_df %>%
    pivot_longer(cols = 2:6,
                 names_to = "portfolio_returns",
                 values_to = "cumulative_log_returns")
  
  #Improving names of the portfolio returns by
  #creating a new column(portfolio_returns_2) 
  #and renaming the sectors 
  q_all_long_df_2<- q_all_long_df %>%
    mutate(portfolio_returns_2 = factor(portfolio_returns, 
             levels = c("port.retq0",
                        "port.ret",
                        "port.retq2",
                        "port.retq10",
                        "ret.long"),
             labels = c("portfolio returns for q=0",
                        "portfolio returns for q=1",
                        "portfolio returns for q=2",
                        "portfolio returns for q=10",
                        "portfolio returns for long return")))
  
  #Creating the graph
  q_all_plot<-ggplot(q_all_long_df_2,
                     aes(x = Index,
                         y = cumulative_log_returns,
                         colour = portfolio_returns_2)) +
    geom_line(size=1.25) +
    theme(text = element_text(size=20))+
    scale_color_manual(values=c("Blue","Red","black",
                                "orange","green"))+
    labs(x = "Date",
         y = "Cumulative log returns",
         colour = "Cumulative log returns")+
    #Adding the start of Covid line
    geom_vline(xintercept = as.numeric(q_all_long_df_2$
                                         Index[4009]),
               linetype="dotted",color = "violetred", size=1)+
    annotate("text",x=c(q_all_long_df_2$Index[4009]),y=70,
             label=c("March 11th 2020"),size=6,hjust=0)+
    #Adding last value for each Sector
    annotate("text",x=c(q_all_long_df_2$Index[5151]),y=28.33528,
             label=c("28.335"),size=6,hjust=0)+
    annotate("text",x=c(q_all_long_df_2$Index[5152]),y=33.67388,
             label=c("33.674"),size=6,hjust=0)+
    annotate("text",x=c(q_all_long_df_2$Index[5153]),y=39.01248,
             label=c("39.012"),size=6,hjust=0)+
    annotate("text",x=c(q_all_long_df_2$Index[5154]),y=81.721,
             label=c("81.721"),size=6,hjust=0)+
    annotate("text",x=c(q_all_long_df_2$Index[5154]),y=0.263376,
             label=c("0.263"),size=6,hjust=0)
  q_all_plot
  ###Plotting portfolio weights for Q=0###
  
  #Making us.sect.weights.q0_df into a data frame
  us.sect.weights.q0_df <- fortify(us.sect.weights.q0)
  
  # sorts the data into a 3 columns table; date,sector
  #and Proportion
  us.sect.q0_long_df <- us.sect.weights.q0_df %>%
    pivot_longer(cols = 2:5,
                 names_to = "Sector",
                 values_to = "Proportion")
  #Improving names of the sectors by creating
  #a new column(Sector_new) and renaming the sectors
  us.sect.q0_long_df <- us.sect.q0_long_df %>%
    mutate(Sector_new = factor(Sector, 
                               levels = c("ConsumerServ",
                                          "Financials",
                                          "Healtcare",
                                          "RealEstate"),
                               labels = c("Consumer Services",
                                          "Financial Services",
                                          "Health Care",
                                          "Real Estate")))
  #Deleting the first 300 values of each sector
  #as these are 0 due to backtesting 
  us.sect.q0_long_df2<-us.sect.q0_long_df[-c(1:1200),]
  
  #Plotting the graph using ggplot
  us.sect.q0_long_plot<-ggplot(us.sect.q0_long_df2,
                                    aes(x = Index,
                                        y = Proportion,
                                        colour = Sector_new)) +
    geom_line(size=1.25) +
    theme(text = element_text(size=20))+
    scale_color_manual(values=c("Blue",
                                "Red",
                                "black",
                                "orange"))+
    labs(x = "Date",
         y = "proportion of investment",
         colour = "Sector Type")+
    #Adding the start of Covid line 
    geom_vline(xintercept = as.numeric(
                            us.sect.q0_long_df$Index[3205]),
               linetype="dotted",color = "violetred", size=1.5)+
    annotate("text",x=c(us.sect.q0_long_df$Index[3220]),y=3,
             label=c("March 11th 2020"),size=6,hjust=0)
  
  ###Plotting portfolio weights for Q=1###
  
  #Making us.sect.weights.q1_df into a data frame
  us.sect.q1_df <- fortify(us.sect.weights.q1)
  
  # sorts the data into a 3 columns table; date,sector 
  #and Proportion
  us.sect.q1_long_df <- us.sect.q1_df %>%
    pivot_longer(cols = 2:5,
                 names_to = "Sector",
                 values_to = "Proportion")
  
  #Improving names of the sectors by creating
  #a new column(Sector_new) and renaming the sectors
  us.sect.q1_long_df <- us.sect.q1_long_df %>%
    mutate(Sector_new = factor(Sector, 
                               levels = c("ConsumerServ",
                                          "Financials",
                                          "Healtcare",
                                          "RealEstate"),
                               labels = c("Consumer Services",
                                          "Financial Services",
                                          "Health Care",
                                          "Real Estate")))
  #Deleting the first 300 values of each sector as these are 0
  #due to backtesting
  us.sect.q1_long_df2<-us.sect.q1_long_df[-c(1:1200),]
  
  
  #Plotting the graph using ggplot
  us.sect.q1_long_plot<-ggplot(us.sect.q1_long_df2,
                                  aes(x = Index,
                                      y = Proportion,
                                      colour = Sector_new)) +
    geom_line(size=1.25) +
    theme(text = element_text(size=20))+
    scale_color_manual(values=c("Blue",
                                "Red",
                                "black",
                                "orange"))+
    labs(x = "Date",
         y = "proportion of investment",
         colour = "Sector Type")+
    geom_vline(xintercept = as.numeric(
                            us.sect.q1_long_df$Index[3205]),
               linetype="dotted",color = "violetred", size=1.5)+
    annotate("text",x=c(us.sect.q1_long_df$Index[3220]),y=3,
             label=c("March 11th 2020"),size=6,hjust=0)
  
  
  ###Correlation plots###
  #Correlation plot for 31/01/2020 to 09/02/2021
  pairsdata <- data.frame(us.sect$Healtcare[775:1033],
                          us.sect$RealEstate[775:1033],
                          us.sect$Financials[775:1033],
                          us.sect$ConsumerServ[775:1033])
  #Correlation plot for 03/01/2017 to 09/02/2021
  pairsdata1 <- data.frame(us.sect$Healtcare,
                           us.sect$RealEstate,
                           us.sect$Financials,
                           us.sect$ConsumerServ)
  #plotting pairsdata1
  pairsdata_plot1<-ggpairs(pairsdata1,
                           lower = list(continuous = 
                                          wrap("smooth",
                                               alpha = 0.3,
                                               size=0.3)),
                          columnLabels = c("Healthcare",
                                           "Real Estate",
                                           "Financials",
                                           "Consumer Services"),
                        title = 'Pairs Correlation for US.Sect')
  
  #plotting pairsdata
  pairsdata_plot<-ggpairs(pairsdata,
                          lower = list(continuous = 
                                         wrap("smooth",
                                              alpha = 0.3,
                                              size=0.3)),
                          columnLabels = c("Healthcare",
                                           "Real Estate",
                                           "Financials",
                                           "Consumer Services"),
  title = 'Pairs Correlation for US.Sect During COVID Pandemic')
  
  plot_list<-list(q_plot,
                  q_all_plot,
                  us.sect.q0_long_plot,
                  us.sect.q1_long_plot,
                  pairsdata_plot,
                  pairsdata_plot1)
  return(plot_list)
}
DCC()

##---- EWMADCC
#A function that calculates 
#Exponential Weighted Moving Average of the data set
EWMA<-function(){
  window<-300
  ret<-rep(0,1033)
  EWMA.DCC<-function(xts,q){
    for(t in(window+1):1033){
      winxts<- xts[(t-window):t]
      #Calculating the inverse conditional covariance matrix
      a<-c(0.004,0.005,0.006,0.007)
      A<-diag(0.2,4)
      B<-diag(0.35,4)
      inip<-c(0.15,0.8)
      rtnval<-dcc.estimation(inia=a,iniA=A,iniB=B,ini.dcc=inip,
                             dvar=winxts,model='diagonal')
      hmtx<-rtnval$h
      DCCmtx<-rtnval$DCC
      htail<-sqrt(tail(hmtx,n=1))
      tDCCmtx<-tail(DCCmtx,1)
      R=matrix(tDCCmtx,nrow=4,byrow=T)
      DCCM <- (t(htail)%*%htail)*R
      Hinv<-solve(DCCM)
      #Calculating the forecasts for the data set
      Health.Hma<-HMA(winxts$Healtcare,7)
      Real.Hma<-HMA(winxts$RealEstate,7)
      Finance.Hma<-HMA(winxts$Financials,7)
      Consumer.Hma<-HMA(winxts$ConsumerServ,7)
      us.HMA<-cbind(Health.Hma,Real.Hma,Finance.Hma,
                    Consumer.Hma)
      #Locating the most recent forecast in this case
      forecast<-as.numeric(tail(us.HMA,1))
      
      #Applying the weights equation to the data
      
      onea <- c(1,1,1,1)
      a<- t(onea) %*% Hinv %*% onea
      b<- t(onea) %*% Hinv %*% forecast 
      wi<-(Hinv %*% onea)/(c(a))+q*(Hinv%*%(forecast-c((b/a))))
      
      #Calculating average returns and
      #applying the weights to them
      
      logts<-log(xts)
      logtsdiff<-diff(logts)*100
      logtsdiff[is.na(logtsdiff)]<-0  
      ei<- apply(logtsdiff, 2, mean)
      ei
      ei%*%wi
      ret[t]<- c(ei%*%wi)
      
    }
    return(ret)
  }
  ewma.ret.q0<-EWMA.DCC(us.sect,0)
  ewma.ret.q1<-EWMA.DCC(us.sect,1)
  ewma.ret.q2<-EWMA.DCC(us.sect,2)
  ewma.ret.q10<-EWMA.DCC(us.sect,10)
  
  
  #Creating data frames 
  ewma.ret.q0_df<-as.data.frame(ewma.ret.q0)
  ewma.ret.q1_df<-as.data.frame(ewma.ret.q1)
  ewma.ret.q2_df<-as.data.frame(ewma.ret.q2)
  ewma.ret.q10_df<-as.data.frame(ewma.ret.q10)
  us_sect_df <- fortify(us.sect)
  
  #Deleting the healthcare,RealEstate,Financials 
  #and ConsumerServ columns from us.sect
  # We did this as the dates we had was not
  #everyday between 03-01-2017 and 09-02-2021 
  # as there is no trading on the weekends 
  date<-select(us_sect_df,
               -Healtcare,
               -RealEstate,
               -Financials,
               -ConsumerServ)
  
  #Putting all of the portfolio returns into a data frame
  #along with date
  ewma_df<-data.frame(date,
                      ewma.ret.q0,
                      ewma.ret.q1,
                      ewma.ret.q2,
                      ewma.ret.q10_df)
  
  #Creates 3 columns; date,portfolio_returns 
  #and cumulative_log_return
  ewma_long_df <- ewma_df %>%
    pivot_longer(cols = 2:5,
                 names_to = "portfolio_returns",
                 values_to = "log_returns")
  
  #Improving names of the portfolio returns by
  #creating a new column(portfolio_returns_2) 
  #and renaming the returns 
  ewma_long_df_2<- ewma_long_df %>%
    mutate(portfolio_returns_2 = factor(portfolio_returns,
              levels = c("ewma.ret.q0",
                         "ewma.ret.q1",
                         "ewma.ret.q2",
                         "ewma.ret.q10"),
              labels = c("EWMA portfolio returns for q=0",
                         "EWMA portfolio returns for q=1",
                         "EWMA portfolio returns for q=2",
                         "EWMA portfolio returns for q=10")))
  
  #Creating the graph
  ewma_plot<-ggplot(ewma_long_df_2,
                    aes(x = Index,
                        y = log_returns,
                        colour = portfolio_returns_2)) +
    geom_line(size=1.25) +
    theme(text = element_text(size=20))+
    scale_color_manual(values=c("black",
                                "blue",
                                "green",
                                "orange"))+
    labs(x = "Date",
         y = "log returns",
         colour = "Risk appetite")+
    
    #Adding the start of COVID-19 line
    geom_vline(xintercept = as.numeric(
                 ewma_long_df_2$Index[3205]),
               linetype="dotted",color = "violetred", size=1.5)+
    annotate("text",x=c(ewma_long_df_2$Index[3205]),y=2,
             label=c("March 11th 2020"),size=6,hjust=0)

  
  ###Calculating the cumulative sum
  ewma.ret.q0_cumsum_df<-as.data.frame(cumsum(ewma.ret.q0))
  ewma.ret.q1_cumsum_df<-as.data.frame(cumsum(ewma.ret.q1))
  ewma.ret.q2_cumsum_df<-as.data.frame(cumsum(ewma.ret.q2))
  ewma.ret.q10_cumsum_df<-as.data.frame(cumsum(ewma.ret.q10))
  
  ewma_cumsum_df<-data.frame(date,
                             ewma.ret.q0_cumsum_df,
                             ewma.ret.q1_cumsum_df,
                             ewma.ret.q2_cumsum_df,
                             ewma.ret.q10_cumsum_df)
  
  #Creates 3 columns; date,portfolio_returns 
  #and cumulative_log_return
  ewma_cumsum_long_df <- ewma_cumsum_df %>%
    pivot_longer(cols = 2:5,
                 names_to = "portfolio_returns",
                 values_to = "cumulative_log_returns")
  
  #Improving names of the portfolio returns by creating
  #a new column(portfolio_returns_2) 
  #and renaming the returns 
  ewma_cumsum_long_df_2<- ewma_cumsum_long_df %>%
    mutate(portfolio_returns_2 = factor(portfolio_returns, 
                 levels = c("cumsum.ewma.ret.q0.",
                            "cumsum.ewma.ret.q1.",
                            "cumsum.ewma.ret.q2.",
                            "cumsum.ewma.ret.q10."),
                 labels = c("EWMA portfolio returns for q=0",
                            "EWMA portfolio returns for q=1",
                            "EWMA portfolio returns for q=2",
                            "EWMA portfolio returns for q=10")))
  
  #Creating the graph
  ewma_cumsum_plot<-ggplot(ewma_cumsum_long_df_2,
                           aes(x = Index,
                               y = cumulative_log_returns,
                               colour = portfolio_returns_2)) +
    geom_line(size=1.25) +
    theme(text = element_text(size=20))+
    scale_color_manual(values=c("black",
                                "blue",
                                "green",
                                "orange"))+
    labs(x = "Date",
         y = "cumulative log returns",
         colour = "Risk appetite")+
    #Adding the start of COVID-19 line
    geom_vline(xintercept = as.numeric(
          ewma_cumsum_long_df_2$Index[3205]),
          linetype="dotted",color = "violetred", size=1.5)+
    annotate("text",x=c(ewma_cumsum_long_df_2$Index[3205]),
             y=47,label=c("March 11th 2020"),size=6,hjust=0)+
    #Adding last value for each Sector
    annotate("text",x=c(ewma_cumsum_long_df_2$Index[4129]),
             y=43,label=c("44.18484"),size=6,hjust=0)+
    annotate("text",x=c(ewma_cumsum_long_df_2$Index[4130]),
             y=44.85993,label=c("44.85993"),size=6,hjust=0)+
    annotate("text",x=c(ewma_cumsum_long_df_2$Index[4131]),
             y=46.5,label=c("45.53501"),size=6,hjust=0)+
    annotate("text",x=c(ewma_cumsum_long_df_2$Index[4132]),
             y=50.93566,label=c("50.93566"),size=6,hjust=0)
  EWMA_list<-list(ewma_plot,ewma_cumsum_plot)
  return(EWMA_list)
}
EWMA()

## ---- CaJo
#A function that calculates the co-integration
#of  the  data  set 
cajo<-function(){
  lsect<-log(us.sect) #Take log of data
  #Locates the lag length for the VAR
  var.fit<-VAR(lsect,lag.max=6,ic='AIC')
  var.fit$p
  #Lag length is 5
  us.coint<-ca.jo(lsect, type='eigen',ecdet="trend",K=5) 
  # The data can be considered to follow a restricted trend. 
  #This is due to the falls and steady recovery 
  #due to the pandemic.
  summary(us.coint)
  #This shows there are no co-integrations. In order to 
  #continue with this line of forecasting, we will use r=1,
  #as it has the smallest value between the test statistic and
  #the critical values.
  us.var<-vec2var(us.coint, r=1)
  # We use one step ahead for these forecasts for accuracy.
  pred<-predict(us.var, n.ahead=1) 
  tai.sect<-tail(lsect,n=1)
  pred.mat<-c(pred$fcst$Healtcare[1],pred$fcst$RealEstate[1], 
          pred$fcst$Financials[1],pred$fcst$ConsumerServ[1])
  ans<-pred.mat-tai.sect
  #Algorithmic window test for our forecast
  
  
  fore.ret<-NULL
  window<-300
  us.ret<-us.sect
  
  
  for(t in(window+1):1033){
    winxts<-lsect[(t-window):t]
    us.coint<-ca.jo(winxts, type='eigen',ecdet="trend",K=5)
    us.var<-vec2var(us.coint, r=1)
    pred<-predict(us.var, n.ahead=1)
    tai.sect<-tail(winxts,n=1)
    pred.mat<-c(pred$fcst$Healtcare[1],
                pred$fcst$RealEstate[1], 
                pred$fcst$Financials[1],
                pred$fcst$ConsumerServ[1])
    val<-pred.mat-tai.sect
    
    fore.ret<-rbind(fore.ret,val)
  }
  
  
  #Here we see from our data that we aren't able to correctly
  #predict and circumvent the impactsof the pandemic. 
  #This means that our original idea of GARCH modelling 
  #with an adaptive portfolio gets better results.
  
  ###Plotting Cajo graph###
  
  #Making the cumulative sum of fore.ret into a data frame
  fore.ret_df<-fortify(cumsum(fore.ret))
  # sorts the data into a 3 columns table; date,sector 
  #and Cumulative_log_return
  fore.ret_long_df <- fore.ret_df %>%
    pivot_longer(cols = 2:5,
                 names_to = "Sector",
                 values_to = "Cumulative_log_return")
  #Improving names of the sectors by creating 
  #a new column(Sector_new) and renaming the sectors
  fore.ret_long_df <- fore.ret_long_df %>%
    mutate(Sector_new = factor(Sector, 
                               levels = c("ConsumerServ",
                                          "Financials",
                                          "Healtcare",
                                          "RealEstate"),
                               labels = c("Consumer Services",
                                          "Financial Services",
                                          "Healthcare",
                                          "Real Estate")))
  # Plotting the graph using ggplot
  fore.ret_plot<-ggplot(fore.ret_long_df,
                        aes(x = Index,
                            y = Cumulative_log_return,
                            colour = Sector_new)) +
    geom_line(size=1.25) +
    theme(text = element_text(size=20))+
    scale_color_manual(values=c("Blue","Red","black","orange"))+
    labs(x = "Date",
         y = "Cumulative log return",
         colour = "Sector Type")+
    #Adding the start of Covid line
    geom_vline(xintercept = as.numeric(
                fore.ret_long_df$Index[2005]),
               linetype="dotted",color = "violetred", size=1.5)+
    annotate("text",x=c(fore.ret_long_df$Index[2020]),y=0.8,
             label=c("March 11th 2020"),size=6,hjust=0)+
    #Adding First and last closing price for each Sector
    annotate("text",x=c(fore.ret_long_df$Index[2929]),y=0.92,
             label=c("0.9291"),size=6,hjust=0)+
    annotate("text",x=c(fore.ret_long_df$Index[2930]),y=0.95,
             label=c("0.9335"),size=6,hjust=0)+
    annotate("text",x=c(fore.ret_long_df$Index[2931]),y=0.6348,
             label=c("0.6348"),size=6,hjust=0)+
    annotate("text",x=c(fore.ret_long_df$Index[2932]),y=0.8416,
             label=c("0.8416"),size=6,hjust=0)
  return(fore.ret_plot)
}
#Runs the cajo function
cajo()

