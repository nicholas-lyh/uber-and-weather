### load all packages
library(dplyr)
library(tidyr)
library(tseries)
library(forecast)
library(fma)
library(xts)
library(ggplot2)
library(lubridate)
library(calendar)
library(timeDate)
library(zoo)
library(corrplot)
library(stats)
library(tibble)
library(summarytools)
library(lmtest)
library(xtable)
library(stargazer)
library(sarima)
library(EnvStats)
### load the data and create a xts object for the data
data <- read.csv("uberLondon.csv", header = TRUE, sep = " ", colClasses = c(Date = "Date"))
dataTS <- xts(data[,c(2:7)], data$Date)

### create a training set and test dataset for the time series
trainingset <- window(dataTS, end = as.Date("2017-11-30"))
testingset <- window(dataTS, start = as.Date("2017-12-01")) 

### summary statistics for the training set
summary <- summary(trainingset)
write.table(summary, file = "summary for training set.txt", sep = ",", row.names = F, quote = F)
dfSummary(trainingset)

### correlation for all variables
traindroppedna <- na.omit(trainingset) 
correlation <- cor(traindroppedna)
corrplot(correlation, method = "number", type = "lower", 
         title = "", diag = F, outline = F,
         tl.cex = 0.7, tl.col = "BLACK", tl.srt = 360)
correlation

### qq plot and histogram for Mean Travel Time 
# Q-Q Plot
qqnorm(trainingset$MeanTravelTimeSeconds,
       main = "Mean Travel Time (Seconds) - Q-Q Plot",
       ylab = "Mean Travel Time (Seconds)")
qqline(trainingset$MeanTravelTimeSeconds, col = "RED", lwd = 1)

# Histogram
ggplot(trainingset$MeanTravelTimeSeconds, aes(x = MeanTravelTimeSeconds)) + 
  geom_freqpoly(binwidth = 60, color = "red", lwd = 1) + theme_minimal() + 
  ggtitle("Mean Travel Time (Seconds) - Histogram") + 
  xlab("Mean Travel Time (Seconds)") +  ylab("Count") +
  geom_histogram(binwidth = 60, fill = "white", color = "black", alpha = 0.5) +
  geom_vline(aes(xintercept = mean(trainingset$MeanTravelTimeSeconds)), 
             color = "blue", linetype = "dashed", size = 1)


### find out outliers of Mean Travel Time and their respective date
outlierdf <- data.frame(date=as.Date(index(trainingset)),y = trainingset$MeanTravelTimeSeconds)
outlier <- boxplot.stats(outlierdf$MeanTravelTimeSeconds)$out
outlier_ind <- which(outlierdf$MeanTravelTimeSecond %in% c(outlier))
outlier
outlierdf[outlier_ind,]
# Rosner's Test of the outliers
rosneroutlier <- rosnerTest(trainingset$MeanTravelTimeSeconds, k = 5)
rosneroutlier$all.stats

### histogram of tavg
par(mfrow=c(1,2))
ggplot(trainingset$tavg, aes(x = trainingset$tavg)) + 
  theme_minimal() + ggtitle("Temperature - Histogram") + 
  xlab("Degree Celsius") +  ylab("Count") +
  geom_histogram(bins = 30, fill = "grey", color = "black", alpha = 0.7) +
  geom_vline(aes(xintercept = mean(trainingset$tavg)), 
             color = "blue", linetype = "dashed", size = 1)

### Frequency Plot of Precipitation
ggplot(trainingset$prcp, aes(x = trainingset$prcp)) + theme_minimal() +
  geom_density() + xlab("Precipitation (mm)") + xlim(0,20) 
  ggtitle("Precipitation - Density Plot") + 
par(mfrow=c(1,1))

### Find out all NA values and their respective dates
trainingdf <- as.data.frame(trainingset)
allnadf <- apply(is.na(trainingdf),2,which)
allnadf

### --------------------------------------------------
### Q3
### Auto ARIMA model 
autoarima1 <- auto.arima(coredata(trainingset$MeanTravelTimeSeconds), seasonal = T)
autoarima1
autoarima1_forecast <- as.data.frame(forecast(autoarima1, h = 31))
autoarima1_forecast_value <- autoarima1_forecast[,1]
RMSEautoarima1 <- sqrt(sum(autoarima1_forecast_value - 
                             coredata(testingset$MeanTravelTimeSeconds))^2/31)
MAPEautoarima1 <- sum(abs(autoarima1_forecast_value - 
                            coredata(testingset$MeanTravelTimeSeconds))/
                        coredata(testingset$MeanTravelTimeSeconds))/31
MAEautoarima1 <- sum(abs(autoarima1_forecast_value - 
                           coredata(testingset$MeanTravelTimeSeconds)))/31

### auto arima with approximation = F, stepwise = F
autoarima2 <- auto.arima(coredata(trainingset$MeanTravelTimeSeconds), 
                         seasonal = T, approximation = F, stepwise = F)
autoarima2_forecast <- as.data.frame(forecast(autoarima2, h = 31))
autoarima2_forecast_value <- autoarima2_forecast[,1]
RMSEautoarima2 <- sqrt(sum(autoarima2_forecast_value - 
                             coredata(testingset$MeanTravelTimeSeconds))^2/31)
MAPEautoarima2 <- sum(abs(autoarima2_forecast_value - 
          coredata(testingset$MeanTravelTimeSeconds))/
            coredata(testingset$MeanTravelTimeSeconds))/31
MAEautoarima2 <- sum(abs(autoarima2_forecast_value - 
          coredata(testingset$MeanTravelTimeSeconds)))/31

### Dynamic Harmonic regression 
trainingts <- ts(trainingset$MeanTravelTimeSeconds, frequency = 365, c(2016,2))
dhr1 <- auto.arima(trainingts, xreg = fourier(trainingts, K = 1), 
                   seasonal = FALSE, lambda = 0)

dhr2 <- auto.arima(trainingts, xreg = fourier(trainingts, K = 2), 
                           seasonal = FALSE, lambda = 0)

dhr3 <- auto.arima(trainingts, xreg = fourier(trainingts, K = 3), 
                   seasonal = FALSE, lambda = 0)

dhr4 <- auto.arima(trainingts, xreg = fourier(trainingts, K = 4), 
                   seasonal = FALSE, lambda = 0)

dhr5 <- auto.arima(trainingts, xreg = fourier(trainingts, K = 5), 
                   seasonal = FALSE, lambda = 0)

# calculating RMSE, MAE, MAPE for ARIMA(5,1,1) with K = 2
dhr2_forecast_df <- as.data.frame(forecast(
  dhr2, xreg = fourier(trainingts, K = 2), h = 31))
dhr2_forecast_value <- dhr2_forecast_df[1:31,1]
dhr2_forecast <- as.data.frame(dhr2_forecast_value)
RMSEdhr2 <- sqrt(sum(dhr2_forecast_value - 
                       coredata(testingset$MeanTravelTimeSeconds))^2/31)
MAPEdhr2 <- sum(abs(dhr2_forecast_value - coredata(testingset$MeanTravelTimeSeconds))/
      coredata(testingset$MeanTravelTimeSeconds))/31
MAEdhr2 <- sum(abs(dhr2_forecast_value - coredata(testingset$MeanTravelTimeSeconds)))/31

### Holt-Winters' Method
trainingzoo <- ts(trainingset$MeanTravelTimeSeconds, frequency = 12)
hw_a <- hw(trainingzoo,seasonal="additive")
plot(as.vector(hw_a$fitted), 
     as.vector(hw_a$residuals), xlab = "Fitted Values", ylab = "Residuals")
hist(hw_a$residuals)
plot(hw_a$residuals)
acf(hw_a$residuals, main = "HW Model (Additive)")

hw_m <- hw(trainingzoo,seasonal="multiplicative")
plot(as.vector(hw_m$fitted), 
     as.vector(hw_m$residuals), xlab = "Fitted Values", ylab = "Residuals")
hist(hw_m$residuals, main = "")
plot(hw_m$residuals)
acf(hw_m$residuals, main = "HW Model (Multiplicative)")
plot(hw_m$fitted)

par(mfrow = c(1,2))
acf(hw_a$residuals, main = "HW Model (Additive)")
acf(hw_m$residuals, main = "HW Model (Multiplicative)")
par(mfrow = c(1,1))


### Forecasting via Neural Network
set.seed(101)
neural <- nnetar(trainingset$MeanTravelTimeSeconds)
neural_forecast <- as.data.frame(forecast(neural, h = 31))
RMSEneural <- sqrt(sum(neural_forecast - 
                         coredata(testingset$MeanTravelTimeSeconds))^2/31)

MAPEneural <- sum(abs(neural_forecast - 
                        coredata(testingset$MeanTravelTimeSeconds))/
                    coredata(testingset$MeanTravelTimeSeconds))/31

MAEneural <- sum(abs(neural_forecast - 
                       coredata(testingset$MeanTravelTimeSeconds)))/31


### Comparing all models
neural_forecast_value <- lapply(neural_forecast, as.numeric)
testing_value <- lapply(testingset$MeanTravelTimeSeconds, as.numeric)
comparison <- data.frame(testing_value, 
                         autoarima1_forecast_value, autoarima2_forecast_value,
                         dhr2_forecast_value, neural_forecast_value)   
colnames(comparison) <- c("Actual", "Autoarima1", 
                          "Autoarima2", 
                          "DHR",
                          "Neural")
plot(comparison$Actual, type = "l", 
     main = "Comparison of all Forecsats", lwd = 2, 
     ylab = "Mean Travel Time (Seconds)", xlab = "")
lines(comparison$Autoarima1, col = "RED", lwd = 2)
lines(comparison$Autoarima2, col = "DARKGREEN", lwd = 2)
lines(comparison$DHR, col = "BLUE", lwd = 2)
lines(comparison$Neural, col = "PURPLE", lwd = 2)
legend("bottomleft", legend=c("Actual Values", "Autoarima1", "Autoarima2", "DHR", "Neural"), col = c("BLACK", "RED", "DARKGREEN", "BLUE", "PURPLE"), lty=1:1, cex=0.5)
