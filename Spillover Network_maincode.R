# clear all variables
rm(list = ls(all = TRUE))
graphics.off()
# set the working directory
setwd("~/")

libraries <-  c("vars","urca","frequencyConnectedness","MASS", 'parallel',
                "strucchange","zoo","sandwich","lmtest","BigVAR","RcppEigen",
                "evaluate","highr","yaml","xfun","knitr","pbapply") 
lapply(libraries, function(x) if (!(x %in% installed.packages())) {
  install.packages(x)
})
lapply(libraries, library, quietly = TRUE, character.only = TRUE)


## 0. load data
Date <- as.data.frame(read.csv('GBN_Price_Daily.csv', header = TRUE)[-1,1]); colnames(date) <- 'Date'
Data <- as.data.frame(read.csv('GBN_Price_Daily.csv', header = TRUE)[,-1]) 
Bankname <- as.data.frame(read.csv('Bank name.csv', header = FALSE))
ins <- 215   # 215 sample banks
obs <- 2902  # 2902 observations
Return <- matrix(0, nrow = obs, ncol = ins); colnames(Return) <- t(bankname)
Volatility <- matrix(0, nrow = obs, ncol = ins); colnames(Volatility) <- t(bankname)

## 1. Calculate the return and volatility (Garman and Klass, 1980)
for (i in 1:ins) {
  tempins <- 5*i-1
  for (t in 1:obs) {
    Return[t,i] <- log(data[t+1,tempins]/data[t, tempins])
    PO <- log(data[t, 1+5*(i-1)])  
    PH <- log(data[t, 2+5*(i-1)])  
    PL <- log(data[t, 3+5*(i-1)]) 
    PC <- log(data[t, 4+5*(i-1)])  
    Volatility[t,i] <- 0.511*(PH-PL)^2-0.019*((PC-PO)*(PH+PL-2*PO)-2*(PH-PO)*(PL-PO))-0.383*(PC-PO)^2
  }
}

# Reorder the matrix by region
RF_Index      <- read.table('regionlist_Index.txt', header = T)
RF_Return     <- Return[, RF_Index$Index]
RF_Volatility <- Volatility[, RF_Index$Index]

## 2. construct the connectendess network by using the DY spillover method
#  2.1. construct BigVAR model with the Elastic network estimator
big_var_est <- function(data) {
  library(BigVAR)
  Model1 = constructModel(as.matrix(data), p = 4, struct = "BasicEN", gran = c(50, 50), VARX = list(), verbose = F)
  Model1Results = cv.BigVAR(Model1)
}


#  2.2. Perform the Full-sample estimation
#  2.2.1. Full sample return spillover network
ReturnVAR     <- big_var_est(RF_Return)
ReturnSpilloverFullSample <- spilloverDY12(ReturnVAR, n.ahead = 100, no.corr = F)          # 全样本的收益溢出网络

#  2.2.2. Full sample volatility spillover network
VolatilityVAR <- big_var_est(RF_Volatility)
VolatilitySpilloverFullSample <- spilloverDY12(VolatilityVAR, n.ahead = 100, no.corr = F)  # 全样本的波动溢出网络

#  2.2.3. save the full-sample result
ReturnFS_Table   <- ReturnSpilloverFullSample$tables
VolatilityFS_Table <- VolatilitySpilloverFullSample$tables
# write.csv(ReturnFS_Table, file = '~/ReturnFS_Table.csv')
# write.csv(VolatilityFS_Table, file = '~/VolatilityFS_Table.csv')


## 3. Perform the Rolling windows estimation
Hor = 10  # predictive horizon H
Win = 200 # window size w

# 3.1. parallel computation for return spillover network
params_est1 <- list()
clnum <- detectCores(logical = FALSE)
cl1 <- makeCluster(clnum)
clusterExport(cl1, c('big_var_est', 'params_est', 'RF_Return', 'Hor', 'Win'))
clusterEvalQ(cl1, c("vars","urca","frequencyConnectedness","MASS", 'parallel',
                   "strucchange","zoo","sandwich","lmtest","BigVAR","RcppEigen",
                   "evaluate","highr","yaml","xfun","knitr","pbapply"))
ReturnSpilloversRoiling <- spilloverRollingDY12(RF_Return, n.ahead = Hor, no.corr = F, func_est = "big_var_est", params_est = params_est1, window = Win, cluster = cl1)
stopCluster(cl1)

# 3.2. parallel computation for volatility spillover network
params_est2 <- list()
clnum <- detectCores(logical = FALSE)    #查看CPU可使用核的数量
cl2 <- makeCluster(clnum)
clusterExport(cl2, c('big_var_est', 'params_est', 'RF_Volatility', 'Hor', 'Win'))
clusterEvalQ(cl2, c("vars","urca","frequencyConnectedness","MASS", 'parallel',
                   "strucchange","zoo","sandwich","lmtest","BigVAR","RcppEigen",
                   "evaluate","highr","yaml","xfun","knitr","pbapply"))
VolatilitySpilloversRoiling <- spilloverRollingDY12(RF_Volatility, n.ahead = Hor, no.corr = F, func_est = "big_var_est", params_est = params_est2, window = Win, cluster = cl2)stopCluster(cl2)


# 3.3. save the roll-sample result
ReturnRS_Table   <- ReturnSpilloversRoiling$list_of_tables
VolatilityRS_Table <- VolatilitySpilloversRoiling$list_of_tables
Winnum <- length(ReturnRS_Table)
# set the roll-sample path
path1 <- '~/'  # path for the result of dynamic return spillover network
path2 <- '~/'  # path for the result of dynamic volatility spillover network

for (i in 1:Winnum) {
  ReturnRS_Table_temp <- ReturnRS_Table[[i]][["tables"]]
  VolatilityRS_Table_temp <- VolatilityRS_Table[[i]][["tables"]]
  write.csv(ReturnRS_Table_temp, file = paste(path1, 'ReturnRS_Table_W200_', i, '.csv', sep = ''), row.names = FALSE)
  write.csv(VolatilityRS_Table_temp, file = paste(path2, 'VolatilityRS_Table_W200_', i, '.csv', sep = ''), row.names = FALSE)
}

# Total, From, To Connectedness
RSRS_Total <- overall(ReturnSpilloversRoiling)
RSRS_From  <- from(ReturnSpilloversRoiling)
RSRS_To    <- to(ReturnSpilloversRoiling)

VSRS_Total <- overall(VolatilitySpilloversRoiling)
VSRS_From  <- from(VolatilitySpilloversRoiling)
VSRS_To    <- to(VolatilitySpilloversRoiling)

# write.csv(RSRS_Total, file = paste(path1, 'RSRS_W200_Total.csv', sep = ''), row.names = FALSE)
# write.csv(RSRS_From,  file = paste(path1, 'RSRS_W200_From.csv', sep = ''), row.names = FALSE)
# write.csv(RSRS_To,    file = paste(path1, 'RSRS_W200_To.csv', sep = ''), row.names = FALSE)
# 
# write.csv(VSRS_Total, file = paste(path2, 'VSRS_W200_Total.csv', sep = ''), row.names = FALSE)
# write.csv(VSRS_From,  file = paste(path2, 'VSRS_W200_From.csv', sep = ''), row.names = FALSE)
# write.csv(VSRS_To,    file = paste(path2, 'VSRS_W200_To.csv', sep = ''), row.names = FALSE)


## 4. Block aggregation for full-sample
banknum    <- 215
countrynum <- 46
regionnum  <- 6
country_name  <- read.table("country_name.txt", header = T, sep = '\t')
region_name <- read.table("region_name.txt", header = T, sep = '\t')

sub_region  <- c(0, 6, 104, 52, 38, 8, 7)
subIndex_region <- cumsum(sub_region)
sub_country <- c(0, 4, 2, 29, 20, 12, 4, 4, 3, 6, 4, 5, 3, 3, 1, 
                 5, 3, 1, 1, 2, 5, 3, 8, 2, 5, 1, 4, 3, 3, 1, 2, 
                 1, 2, 2, 1, 1, 1, 4, 1, 31, 6, 1, 6, 2, 2, 3, 2)
subIndex_country <- cumsum(sub_country)

# 4.1. Full-sample block aggregation
ReturnspilloversFS     <- ReturnFS_Table
VolatilityspilloversFS <- VolatilityFS_Table
# ReturnspilloversFS     <- read.csv('ReturnFS_Table.csv', header = TRUE)[,-1]
# VolatilityspilloversFS <- read.csv('VolatilityFS_Table.csv', header = TRUE)[,-1]

# 4.1.1. re-normalisation
ReturnspilloversFST <- ReturnspilloversFS/banknum
VolatilityspilloversFST <- VolatilityspilloversFS/banknum

Returnblock_country     <- matrix(0, nrow = countrynum, ncol = countrynum)
Returnblock_region      <- matrix(0, nrow = regionnum,  ncol = regionnum) 
Volatilityblock_country <- matrix(0, nrow = countrynum, ncol = countrynum)
Volatilityblock_region  <- matrix(0, nrow = regionnum,  ncol = regionnum) 

colnames(Returnblock_country) <- country_name$name
colnames(Returnblock_region) <- region_name$name
colnames(Volatilityblock_country) <- country_name$name
colnames(Volatilityblock_region) <- region_name$name

rownames(Returnblock_country) <- country_name$name
rownames(Returnblock_region) <- region_name$name
rownames(Volatilityblock_country) <- country_name$name
rownames(Volatilityblock_region) <- region_name$name

# 4.1.2. Block Aggregation for the region and country level
# (1) Region Level Block Aggregation
for (i in 1:regionnum) {
  for (j in 1:regionnum) {
    Returnblock_region[i, j] <- sum(ReturnspilloversFST[(subIndex_region[[i]]+1):subIndex_region[[i+1]], (subIndex_region[[j]]+1):subIndex_region[[j+1]]])
    Volatilityblock_region[i, j] <- sum(VolatilityspilloversFST[(subIndex_region[[i]]+1):subIndex_region[[i+1]], (subIndex_region[[j]]+1):subIndex_region[[j+1]]])
  }
}

# (2) Country Level Block Aggregation
for (i in 1:countrynum) {
  for (j in 1:countrynum) {
    Returnblock_country[i, j] <- sum(ReturnspilloversFST[(subIndex_country[[i]]+1):subIndex_country[[i+1]], (subIndex_country[[j]]+1):subIndex_country[[j+1]]])
    Volatilityblock_country[i, j] <- sum(VolatilityspilloversFST[(subIndex_country[[i]]+1):subIndex_country[[i+1]], (subIndex_country[[j]]+1):subIndex_country[[j+1]]])
  }
}

# (3) save the full-sample block aggregation results in the region and country level
write.csv(Returnblock_region, file = 'Returnblock_region.csv')
write.csv(Returnblock_country, file = 'Returnblock_country.csv')
write.csv(Volatilityblock_region, file = 'Volatilityblock_region.csv')
write.csv(Volatilityblock_country, file = 'Volatilityblock_country.csv')

# 4.1.3. Connectedness measures for different aggregation level
# (1) Region level (Within, From, To, Net, Dep, Infl)
RSB_region_con <- matrix(0, nrow = regionnum, ncol = 6)
VSB_region_con <- matrix(0, nrow = regionnum, ncol = 6)
RSBC_consta_name <- c('RSBR_Within', 'RSBR_From', 'RSBR_To', 'RSBR_Net', 'RSBR_Dep', 'RSBR_Infl')
VSBC_consta_name <- c('VSBR_Within', 'VSBR_From', 'VSBR_To', 'VSBR_Net', 'VSBR_Dep', 'VSBR_Infl')
colnames(RSB_region_con) <- RSBC_consta_name
colnames(VSB_region_con) <- VSBC_consta_name
rownames(RSB_region_con) <- region_name$name
rownames(VSB_region_con) <- region_name$name

for (i in 1:regionnum) {
  RSB_region_con[i, 1] <- Returnblock_region[i, i]      # Within
  VSB_region_con[i, 1] <- Volatilityblock_region[i, i]  # Within
  RSB_region_con[i, 2] <- sum(Returnblock_region[i, ]) - Returnblock_region[i, i]         # From
  VSB_region_con[i, 2] <- sum(Volatilityblock_region[i, ]) - Volatilityblock_region[i, i] # From
  RSB_region_con[i, 3] <- sum(Returnblock_region[, i]) - Returnblock_region[i, i]         # To
  VSB_region_con[i, 3] <- sum(Volatilityblock_region[, i]) - Volatilityblock_region[i, i] # To
  RSB_region_con[i, 4] <- RSB_region_con[i, 3] - RSB_region_con[i, 2]             # Net
  VSB_region_con[i, 4] <- VSB_region_con[i, 3] - VSB_region_con[i, 2]             # Net
  RSB_region_con[i, 5] <- RSB_region_con[i, 2] / sum(Returnblock_region[i, ])     # Dep
  VSB_region_con[i, 5] <- VSB_region_con[i, 2] / sum(Volatilityblock_region[i, ]) # Dep
  RSB_region_con[i, 6] <- RSB_region_con[i, 4] / (RSB_region_con[i, 3] + RSB_region_con[i, 2])     # Infl
  VSB_region_con[i, 6] <- VSB_region_con[i, 4] / (VSB_region_con[i, 3] + VSB_region_con[i, 2])     # Infl
}
# write.csv(RSB_region_con, file = '~\RSBR_ConnectednessTable.csv')
# write.csv(VSB_region_con, file = '~\VSBR_ConnectednessTable.csv')

# (2) Country level (Within, From, To, Net, Dep, Infl)
RSB_country_con <- matrix(0, nrow = countrynum, ncol = 6)
VSB_country_con <- matrix(0, nrow = countrynum, ncol = 6)
RSBC_consta_name <- c('RSBC_Within', 'RSBC_From', 'RSBC_To', 'RSBC_Net', 'RSBC_Dep', 'RSBC_Infl')
VSBC_consta_name <- c('VSBC_Within', 'VSBC_From', 'VSBC_To', 'VSBC_Net', 'VSBC_Dep', 'VSBC_Infl')
colnames(RSB_country_con) <- RSBC_consta_name
colnames(VSB_country_con) <- VSBC_consta_name
rownames(RSB_country_con) <- country_name$name
rownames(VSB_country_con) <- country_name$name
for (i in 1:countrynum) {
  RSB_country_con[i, 1] <- Returnblock_country[i, i]       # Wintin
  VSB_country_con[i, 1] <- Volatilityblock_country[i, i]   # Within
  RSB_country_con[i, 2] <- sum(Returnblock_country[i, ]) - Returnblock_country[i, i]         # From
  VSB_country_con[i, 2] <- sum(Volatilityblock_country[i, ]) - Volatilityblock_country[i, i] # From
  RSB_country_con[i, 3] <- sum(Returnblock_country[, i]) - Returnblock_country[i, i]         # To
  VSB_country_con[i, 3] <- sum(Volatilityblock_country[, i]) - Volatilityblock_country[i, i] # To
  RSB_country_con[i, 4] <- RSB_country_con[i, 3] - RSB_country_con[i, 2]                  # Net
  VSB_country_con[i, 4] <- VSB_country_con[i, 3] - VSB_country_con[i, 2]                  # Net
  RSB_country_con[i, 5] <- RSB_country_con[i, 2] / sum(Returnblock_country[i, ])          # Dep
  VSB_country_con[i, 5] <- VSB_country_con[i, 2] / sum(Volatilityblock_country[i,])       # Dep
  RSB_country_con[i, 6] <- RSB_country_con[i, 4] / (RSB_country_con[i, 3] + RSB_country_con[i, 2]) # Infi
  VSB_country_con[i, 6] <- VSB_country_con[i, 4] / (VSB_country_con[i, 3] + VSB_country_con[i, 2]) # Infi
}

aver_RSB_con <- matrix(0, nrow = 1, ncol = 6); rownames(aver_RSB_con) <- 'Average'
aver_VSB_con <- matrix(0, nrow = 1, ncol = 6); rownames(aver_VSB_con) <- 'Average'
for (i in 1:6) {
  aver_RSB_con[1, i] = sum(RSB_country_con[, i])/countrynum
  aver_VSB_con[1, i] = sum(VSB_country_con[, i])/countrynum
}

RSB_country_consta <- rbind(RSB_country_con, aver_RSB_con)
VSB_country_consta <- rbind(VSB_country_con, aver_VSB_con)

write.csv(RSB_country_consta, file = 'RSBC_ConnectednessTable.csv')
write.csv(VSB_country_consta, file = 'VSBC_ConnectednessTable.csv')

# (3) Bank level (Within, From, To, Net, Dep, Infl)
bank_info <- read.csv('Bank info.csv', header = TRUE)
RSB_inst_con <- matrix(0, nrow = banknum, ncol = 4)
VSB_inst_con <- matrix(0, nrow = banknum, ncol = 4)
colnames(RSB_inst_con) <- c('RSBI_Within', 'RSBI_From', 'RSBI_To', 'RSBI_Net')
colnames(VSB_inst_con) <- c('RSBI_Within', 'RSBI_From', 'RSBI_To', 'RSBI_Net')

for (i in 1:banknum) {
  RSB_inst_con[i, 1] <- ReturnspilloversFST[i, i]                                         # Within
  VSB_inst_con[i, 1] <- VolatilityspilloversFST[i, i]                                     # Within
  RSB_inst_con[i, 2] <- sum(ReturnspilloversFST[i, ]) - ReturnspilloversFST[i, i]         # From
  VSB_inst_con[i, 2] <- sum(VolatilityspilloversFST[i, ]) - VolatilityspilloversFST[i, i] # From
  RSB_inst_con[i, 3] <- sum(ReturnspilloversFST[, i]) - ReturnspilloversFST[i, i]         # To
  VSB_inst_con[i, 3] <- sum(VolatilityspilloversFST[, i]) - VolatilityspilloversFST[i, i] # To  
  RSB_inst_con[i, 4] <- RSB_inst_con[i, 3] - RSB_inst_con[i, 2]                           # Net
  VSB_inst_con[i, 4] <- VSB_inst_con[i, 3] - VSB_inst_con[i, 2]                           # Net 
}
RSB_inst_consta <- data.frame(bank_info$Region, bank_info$Country, bank_info$Bank.symbol,RSB_inst_con)
VSB_inst_consta <- data.frame(bank_info$Region, bank_info$Country, bank_info$Bank.symbol,VSB_inst_con)

# write.csv(RSB_inst_consta, file = '~\RSBI_ConnectednessTable.csv')
# write.csv(VSB_inst_consta, file = '~\VSBI_ConnectednessTable.csv')

# 4.2. Rolling-sample block aggregation
TimeLength <- 2902
windowsize <- 200
windowstep <- 1
windownum  <- TimeLength - windowsize + 1

country_name <- read.table("country_name.txt", header = T, sep = '\t')
region_name  <- read.table("region_name.txt", header = T, sep = '\t')


# 4.2.1. re-normalisation
RSRS_Institution <- array(0, dim = c(Institutionnum, Institutionnum, windownum))
RSRS_Country     <- array(0, dim = c(Countrynum, Countrynum, windownum))
RSRS_Region      <- array(0, dim = c(Regionnum, Regionnum, windownum))
VSRS_Institution <- array(0, dim = c(Institutionnum, Institutionnum, windownum))
VSRS_Country     <- array(0, dim = c(Countrynum, Countrynum, windownum))
VSRS_Region      <- array(0, dim = c(Regionnum, Regionnum, windownum))

for (i in 1:windownum) {
  ReturnRS_Table_temp     <- as.matrix(read.csv(file = paste('~/ReturnRS_Table_W200_', i, '.csv', sep = '')))
  VolatilityRS_Table_temp <- as.matrix(read.csv(file = paste('~/VolatilityRS_Table_W200_', i, '.csv', sep = '')))
  ReturnRS_Table_temp     <- ReturnRS_Table_temp / Institutionnum
  VolatilityRS_Table_temp <- VolatilityRS_Table_temp / Institutionnum
  RSRS_Institution[, , i] <- ReturnRS_Table_temp
  VSRS_Institution[, , i] <- VolatilityRS_Table_temp
}

# 4.2.2. Block Aggregation for the region and country level
# (1) Region Level Block Aggregation
for (k in 1:windownum) {
  ReturnRS_temp   <- RSRS_Institution[, , k]
  VolatilityRS_temp <- VSRS_Institution[, , k]
  for (i in 1:Regionnum) {
    for (j in 1:Regionnum) {
      RSRS_Region[i, j, k] <- sum(ReturnRS_temp[(subIndex_region[[i]]+1):subIndex_region[[i+1]], (subIndex_region[[j]]+1):subIndex_region[[j+1]]])
      VSRS_Region[i, j, k] <- sum(VolatilityRS_temp[(subIndex_region[[i]]+1):subIndex_region[[i+1]], (subIndex_region[[j]]+1):subIndex_region[[j+1]]])
    }
  }
}

# (2) Country Level Block Aggregation
for (k in 1:windownum) {
  ReturnRS_temp   <- RSRS_Institution[, , k]
  VolatilityRS_temp <- VSRS_Institution[, , k]
  for (i in 1:Countrynum) {
    for (j in 1:Countrynum) {
      RSRS_Country[i, j, k] <- sum(ReturnRS_temp[(subIndex_country[[i]]+1):subIndex_country[[i+1]], (subIndex_country[[j]]+1):subIndex_country[[j+1]]])
      VSRS_Country[i, j, k] <- sum(VolatilityRS_temp[(subIndex_country[[i]]+1):subIndex_country[[i+1]], (subIndex_country[[j]]+1):subIndex_country[[j+1]]])
    }
  }
}

# (3) save the rolling-sample block aggregation results in the region and country level
for (k in 1:windownum) {
  RSRS_Region_temp  <- RSRS_Region[, , k];  colnames(RSRS_Region_temp)  <- region_name$name
  VSRS_Region_temp  <- VSRS_Region[, , k];  colnames(VSRS_Region_temp)  <- region_name$name
  RSRS_Country_temp <- RSRS_Country[, , k]; colnames(RSRS_Country_temp) <- country_name$name
  VSRS_Country_temp <- VSRS_Country[, , k]; colnames(VSRS_Country_temp) <- country_name$name
  
  write.csv(RSRS_Region_temp,  file = paste('~/RSRS_Region_', k, '.csv', sep = ''), row.names = FALSE)
  write.csv(VSRS_Region_temp,  file = paste('~/VSRS_Region_', k, '.csv', sep = ''), row.names = FALSE)
  write.csv(RSRS_Country_temp, file = paste('~/RSRS_Country_', k, '.csv', sep = ''), row.names = FALSE)
  write.csv(VSRS_Country_temp, file = paste('~/VSRS_Country_', k, '.csv', sep = ''), row.names = FALSE)
}

# 4.2.3. Connectedness measures for different aggregation level

RSRS_con_Table <- matrix(NA, nrow = windownum, ncol = 9)  
VSRS_con_Table <- matrix(NA, nrow = windownum, ncol = 9)
con_Tablename <- c('Ins_spillover', 'Ins_spillself', 'Ins_total', 
                   'Country_spillover', 'Country_spillself', 'Country_total',
                   'Region_spillover', 'Region_spillself', 'Region_total')
colnames(RSRS_con_Table) <- con_Tablename
colnames(VSRS_con_Table) <- con_Tablename

for (k in 1:windownum) {
  # return panel
  RSRS_Institution_temp <- RSRS_Institution[, , k]
  RSRS_Country_temp <- RSRS_Country[, , k]
  RSRS_Region_temp  <- RSRS_Region[, , k]
  RSRS_con_Table[k,1] <- sum(RSRS_Institution_temp[!row(RSRS_Institution_temp)==col(RSRS_Institution_temp)])
  RSRS_con_Table[k,2] <- sum(RSRS_Institution_temp[row(RSRS_Institution_temp)==col(RSRS_Institution_temp)])
  RSRS_con_Table[k,3] <- RSRS_con_Table[k,1] + RSRS_con_Table[k,2]
  RSRS_con_Table[k,4] <- sum(RSRS_Country_temp[!row(RSRS_Country_temp)==col(RSRS_Country_temp)])
  RSRS_con_Table[k,5] <- sum(RSRS_Country_temp[row(RSRS_Country_temp)==col(RSRS_Country_temp)]) - RSRS_con_Table[k,2]
  RSRS_con_Table[k,6] <- RSRS_con_Table[k,4] + RSRS_con_Table[k,5]
  RSRS_con_Table[k,7] <- sum(RSRS_Region_temp[!row(RSRS_Region_temp)==col(RSRS_Region_temp)])
  RSRS_con_Table[k,8] <- sum(RSRS_Region_temp[row(RSRS_Region_temp)==col(RSRS_Region_temp)]) - RSRS_con_Table[k,2]
  RSRS_con_Table[k,9] <- RSRS_con_Table[k,7] + RSRS_con_Table[k,8]
  
  # volatility panel
  VSRS_Institution_temp <- VSRS_Institution[, , k]
  VSRS_Country_temp <- VSRS_Country[, , k]
  VSRS_Region_temp  <- VSRS_Region[, , k]
  VSRS_con_Table[k,1] <- sum(VSRS_Institution_temp[!row(VSRS_Institution_temp)==col(VSRS_Institution_temp)])
  VSRS_con_Table[k,2] <- sum(VSRS_Institution_temp[row(VSRS_Institution_temp)==col(VSRS_Institution_temp)])
  VSRS_con_Table[k,3] <- VSRS_con_Table[k,1] + VSRS_con_Table[k,2]
  VSRS_con_Table[k,4] <- sum(VSRS_Country_temp[!row(VSRS_Country_temp)==col(VSRS_Country_temp)])
  VSRS_con_Table[k,5] <- sum(VSRS_Country_temp[row(VSRS_Country_temp)==col(VSRS_Country_temp)]) - VSRS_con_Table[k,2]
  VSRS_con_Table[k,6] <- VSRS_con_Table[k,4] + VSRS_con_Table[k,5]
  VSRS_con_Table[k,7] <- sum(VSRS_Region_temp[!row(VSRS_Region_temp)==col(VSRS_Region_temp)])
  VSRS_con_Table[k,8] <- sum(VSRS_Region_temp[row(VSRS_Region_temp)==col(VSRS_Region_temp)]) - VSRS_con_Table[k,2]
  VSRS_con_Table[k,9] <- VSRS_con_Table[k,7] + VSRS_con_Table[k,8]
}

write.csv(RSRS_con_Table,  file = '~/RSRS_con_Table.csv', row.names = FALSE)
write.csv(VSRS_con_Table,  file = '~/VSRS_con_Table.csv', row.names = FALSE)


# (1) Country level (From、To、Net、Within、Dep，Infl)
RSRS_Country_From   <- matrix(0, nrow = windownum, ncol = Countrynum)
RSRS_Country_To     <- matrix(0, nrow = windownum, ncol = Countrynum)
RSRS_Country_Net    <- matrix(0, nrow = windownum, ncol = Countrynum)
RSRS_Country_Within <- matrix(0, nrow = windownum, ncol = Countrynum)
RSRS_Country_Dep    <- matrix(0, nrow = windownum, ncol = Countrynum)
RSRS_Country_Infl   <- matrix(0, nrow = windownum, ncol = Countrynum)

VSRS_Country_From   <- matrix(0, nrow = windownum, ncol = Countrynum)
VSRS_Country_To     <- matrix(0, nrow = windownum, ncol = Countrynum)
VSRS_Country_Net    <- matrix(0, nrow = windownum, ncol = Countrynum)
VSRS_Country_Within <- matrix(0, nrow = windownum, ncol = Countrynum)
VSRS_Country_Dep    <- matrix(0, nrow = windownum, ncol = Countrynum)
VSRS_Country_Infl   <- matrix(0, nrow = windownum, ncol = Countrynum)

colnames(RSRS_Country_From)   <- country_name$name
colnames(RSRS_Country_To)     <- country_name$name
colnames(RSRS_Country_Net)    <- country_name$name
colnames(RSRS_Country_Within) <- country_name$name
colnames(RSRS_Country_Dep)    <- country_name$name
colnames(RSRS_Country_Infl)   <- country_name$name

colnames(VSRS_Country_From)   <- country_name$name
colnames(VSRS_Country_To)     <- country_name$name
colnames(VSRS_Country_Net)    <- country_name$name
colnames(VSRS_Country_Within) <- country_name$name
colnames(VSRS_Country_Dep)    <- country_name$name
colnames(VSRS_Country_Infl)   <- country_name$name

for (k in 1:windownum) {
  for (i in 1:Countrynum) {
    RSRS_Country_From[k,i]   <- sum(RSRS_Country[i,,k]) - RSRS_Country[i,i,k]
    RSRS_Country_To[k,i]     <- sum(RSRS_Country[,i,k]) - RSRS_Country[i,i,k]
    RSRS_Country_Net[k,i]    <- RSRS_Country_To[k,i] - RSRS_Country_From[k,i]
    RSRS_Country_Within[k,i] <- RSRS_Country[i,i,k]
    RSRS_Country_Dep[k,i]    <- RSRS_Country_From[k,i] / sum(RSRS_Country[i,,k])
    RSRS_Country_Infl[k,i]   <- RSRS_Country_Net[k,i] / (RSRS_Country_To[k,i] + RSRS_Country_From[k,i])
    
    VSRS_Country_From[k,i]   <- sum(VSRS_Country[i,,k]) - VSRS_Country[i,i,k]
    VSRS_Country_To[k,i]     <- sum(VSRS_Country[,i,k]) - VSRS_Country[i,i,k]
    VSRS_Country_Net[k,i]    <- VSRS_Country_To[k,i] - VSRS_Country_From[k,i]
    VSRS_Country_Within[k,i] <- VSRS_Country[i,i,k]
    VSRS_Country_Dep[k,i]    <- VSRS_Country_From[k,i] / sum(VSRS_Country[i,,k])
    VSRS_Country_Infl[k,i]   <- VSRS_Country_Net[k,i] / (VSRS_Country_To[k,i] + VSRS_Country_From[k,i])
  }
}

# save the country level connectedness (From、To、Net、Within、Dep、Infl)
write.csv(RSRS_Country_From,   file = '~/RSRS_Country_From.csv',   row.names = FALSE)
write.csv(RSRS_Country_To,     file = '~/RSRS_Country_To.csv',     row.names = FALSE)
write.csv(RSRS_Country_Net,    file = '~/RSRS_Country_Net.csv',    row.names = FALSE)
write.csv(RSRS_Country_Within, file = '~/RSRS_Country_Within.csv', row.names = FALSE)
write.csv(RSRS_Country_Dep,    file = '~/RSRS_Country_Dep.csv',    row.names = FALSE)
write.csv(RSRS_Country_Infl,   file = '~/RSRS_Country_Infl.csv',   row.names = FALSE)

write.csv(VSRS_Country_From,   file = '~/VSRS_Country_From.csv',   row.names = FALSE)
write.csv(VSRS_Country_To,     file = '~/VSRS_Country_To.csv',     row.names = FALSE)
write.csv(VSRS_Country_Net,    file = '~/VSRS_Country_Net.csv',    row.names = FALSE)
write.csv(VSRS_Country_Within, file = '~/VSRS_Country_Within.csv', row.names = FALSE)
write.csv(VSRS_Country_Dep,    file = '~/VSRS_Country_Dep.csv',    row.names = FALSE)
write.csv(VSRS_Country_Infl,   file = '~/VSRS_Country_Infl.csv',   row.names = FALSE)


# (2) GDP Top 20 countries and other countries
country_info <- read.csv('country info.csv', header = T)
GDP_TOP20 <- as.data.frame(country_info[1:20, 1]); colnames(GDP_TOP20) = 'GDP_TOP20_Country_Index'
GDP_Other <- as.data.frame(country_info[21:46, 1]); colnames(GDP_Other) = 'GDP_Other_Country_Index'
GDP_Other_banknum <- country_info[21:46,9]

# return panel
RSRS_TOP20Country_From   <- as.matrix(RSRS_Country_From[,GDP_TOP20$GDP_TOP20_Country_Index])
RSRS_TOP20Country_To     <- as.matrix(RSRS_Country_To[,  GDP_TOP20$GDP_TOP20_Country_Index])
RSRS_TOP20Country_Net    <- as.matrix(RSRS_Country_Net[, GDP_TOP20$GDP_TOP20_Country_Index])
RSRS_TOP20Country_Within <- as.matrix(RSRS_Country_Within[,GDP_TOP20$GDP_TOP20_Country_Index])
RSRS_TOP20Country_Dep    <- as.matrix(RSRS_Country_Dep[, GDP_TOP20$GDP_TOP20_Country_Index])
RSRS_TOP20Country_Infl   <- as.matrix(RSRS_Country_Infl[,GDP_TOP20$GDP_TOP20_Country_Index])    

RSRS_OtherCountry_From   <- as.matrix(RSRS_Country_From[,GDP_Other$GDP_Other_Country_Index])
RSRS_OtherCountry_To     <- as.matrix(RSRS_Country_To[,  GDP_Other$GDP_Other_Country_Index])
RSRS_OtherCountry_Net    <- as.matrix(RSRS_Country_Net[, GDP_Other$GDP_Other_Country_Index])
RSRS_OtherCountry_Within <- as.matrix(RSRS_Country_Within[,GDP_Other$GDP_Other_Country_Index])
RSRS_OtherCountry_Dep    <- as.matrix(RSRS_Country_Dep[, GDP_Other$GDP_Other_Country_Index])
RSRS_OtherCountry_Infl   <- as.matrix(RSRS_Country_Infl[,GDP_Other$GDP_Other_Country_Index])

RSRS_Others_From   <- as.matrix(apply(RSRS_OtherCountry_From  , 1, sum)); colnames(RSRS_Others_From)   = 'OtherFrom'
RSRS_Others_To     <- as.matrix(apply(RSRS_OtherCountry_To    , 1, sum)); colnames(RSRS_Others_To)     = 'OtherTo'
RSRS_Others_Net    <- as.matrix(apply(RSRS_OtherCountry_Net   , 1, sum)); colnames(RSRS_Others_Net)    = 'OtherNet'
RSRS_Others_Within <- as.matrix(apply(RSRS_OtherCountry_Within, 1, sum)); colnames(RSRS_Others_Within) = 'OtherWithin'

RSRS_Others_Dep    <- matrix(0,nrow = windownum, ncol = 1); colnames(RSRS_Others_Dep)  = 'OtherDep'
RSRS_Others_Infl   <- matrix(0,nrow = windownum, ncol = 1); colnames(RSRS_Others_Infl) = 'OtherInfl'

for (k in 1:windownum) {
  RSRS_Others_Dep[k,1]  <- weighted.mean(RSRS_OtherCountry_Dep[k,], GDP_Other_banknum)
  RSRS_Others_Infl[k,1] <- weighted.mean(RSRS_OtherCountry_Infl[k,], GDP_Other_banknum)
}

# volatility panel
VSRS_TOP20Country_From   <- as.matrix(VSRS_Country_From[,GDP_TOP20$GDP_TOP20_Country_Index])
VSRS_TOP20Country_To     <- as.matrix(VSRS_Country_To[, GDP_TOP20$GDP_TOP20_Country_Index])
VSRS_TOP20Country_Net    <- as.matrix(VSRS_Country_Net[, GDP_TOP20$GDP_TOP20_Country_Index])
VSRS_TOP20Country_Within <- as.matrix(VSRS_Country_Within[, GDP_TOP20$GDP_TOP20_Country_Index])
VSRS_TOP20Country_Dep    <- as.matrix(VSRS_Country_Dep[, GDP_TOP20$GDP_TOP20_Country_Index])
VSRS_TOP20Country_Infl   <- as.matrix(VSRS_Country_Infl[, GDP_TOP20$GDP_TOP20_Country_Index])

VSRS_OtherCountry_From   <- as.matrix(VSRS_Country_From[,GDP_Other$GDP_Other_Country_Index])
VSRS_OtherCountry_To     <- as.matrix(VSRS_Country_To[,  GDP_Other$GDP_Other_Country_Index])
VSRS_OtherCountry_Net    <- as.matrix(VSRS_Country_Net[, GDP_Other$GDP_Other_Country_Index])
VSRS_OtherCountry_Within <- as.matrix(VSRS_Country_Within[, GDP_Other$GDP_Other_Country_Index])
VSRS_OtherCountry_Dep    <- as.matrix(VSRS_Country_Dep[, GDP_Other$GDP_Other_Country_Index])
VSRS_OtherCountry_Infl   <- as.matrix(VSRS_Country_Infl[, GDP_Other$GDP_Other_Country_Index])

VSRS_Others_From   <- as.matrix(apply(VSRS_OtherCountry_From  , 1, sum)); colnames(VSRS_Others_From)   = 'OtherFrom'
VSRS_Others_To     <- as.matrix(apply(VSRS_OtherCountry_To    , 1, sum)); colnames(VSRS_Others_To)     = 'OtherTo'
VSRS_Others_Net    <- as.matrix(apply(VSRS_OtherCountry_Net   , 1, sum)); colnames(VSRS_Others_Net)    = 'OtherNet'
VSRS_Others_Within <- as.matrix(apply(VSRS_OtherCountry_Within, 1, sum)); colnames(VSRS_Others_Within) = 'OtherWithin'

VSRS_Others_Dep    <- matrix(0,nrow = windownum, ncol = 1); colnames(VSRS_Others_Dep)  <- 'OtherDep'
VSRS_Others_Infl   <- matrix(0,nrow = windownum, ncol = 1); colnames(VSRS_Others_Infl) <- 'OtherInfl'

for (k in 1:windownum) {
  VSRS_Others_Dep[k,1]  <- weighted.mean(VSRS_OtherCountry_Dep[k,], GDP_Other_banknum)
  VSRS_Others_Infl[k,1] <- weighted.mean(VSRS_OtherCountry_Infl[k,], GDP_Other_banknum)
}

# save GDP Top20 countries and other countries' connectedness (From、To、Net、Within、Dep、Infl)
write.csv(RSRS_TOP20Country_From,   file = '~/RSRS_TOP20Country_From.csv',   row.names = FALSE)
write.csv(RSRS_TOP20Country_To,     file = '~/RSRS_TOP20Country_To.csv',     row.names = FALSE)
write.csv(RSRS_TOP20Country_Net,    file = '~/RSRS_TOP20Country_Net.csv',    row.names = FALSE)
write.csv(RSRS_TOP20Country_Within, file = '~/RSRS_TOP20Country_Within.csv', row.names = FALSE)
write.csv(RSRS_TOP20Country_Dep,    file = '~/RSRS_TOP20Country_Dep.csv',    row.names = FALSE)
write.csv(RSRS_TOP20Country_Infl,   file = '~/RSRS_TOP20Country_Infl.csv',   row.names = FALSE)

write.csv(cbind(RSRS_OtherCountry_From, RSRS_Others_From),     file = '~/RSRS_OtherCountry_From.csv', row.names = FALSE)
write.csv(cbind(RSRS_OtherCountry_To, RSRS_Others_To),         file = '~/RSRS_OtherCountry_To.csv',   row.names = FALSE)
write.csv(cbind(RSRS_OtherCountry_Net, RSRS_Others_Net),       file = '~/RSRS_OtherCountry_Net.csv',  row.names = FALSE)
write.csv(cbind(RSRS_OtherCountry_Within, RSRS_Others_Within), file = '~/RSRS_OtherCountry_Within.csv',  row.names = FALSE)
write.csv(cbind(RSRS_OtherCountry_Dep, RSRS_Others_Dep),       file = '~/RSRS_OtherCountry_Dep.csv',  row.names = FALSE)
write.csv(cbind(RSRS_OtherCountry_Infl, RSRS_Others_Infl),     file = '~/RSRS_OtherCountry_Infl.csv',  row.names = FALSE)

write.csv(VSRS_TOP20Country_From,   file = '~/VSRS_TOP20Country_From.csv',   row.names = FALSE)
write.csv(VSRS_TOP20Country_To,     file = '~/VSRS_TOP20Country_To.csv',     row.names = FALSE)
write.csv(VSRS_TOP20Country_Net,    file = '~/VSRS_TOP20Country_Net.csv',    row.names = FALSE)
write.csv(VSRS_TOP20Country_Within, file = '~/VSRS_TOP20Country_Within.csv', row.names = FALSE)
write.csv(VSRS_TOP20Country_Dep,    file = '~/VSRS_TOP20Country_Dep.csv',    row.names = FALSE)
write.csv(VSRS_TOP20Country_Infl,   file = '~/VSRS_TOP20Country_Infl.csv',   row.names = FALSE)

write.csv(cbind(VSRS_OtherCountry_From,   VSRS_Others_From),   file = '~/VSRS_OtherCountry_From.csv',   row.names = FALSE)
write.csv(cbind(VSRS_OtherCountry_To,     VSRS_Others_To),     file = '~/VSRS_OtherCountry_To.csv',     row.names = FALSE)
write.csv(cbind(VSRS_OtherCountry_Net,    VSRS_Others_Net),    file = '~/VSRS_OtherCountry_Net.csv',    row.names = FALSE)
write.csv(cbind(VSRS_OtherCountry_Within, VSRS_Others_Within), file = '~/VSRS_OtherCountry_Within.csv', row.names = FALSE)
write.csv(cbind(VSRS_OtherCountry_Dep,    VSRS_Others_Dep),    file = '~/VSRS_OtherCountry_Dep.csv',    row.names = FALSE)
write.csv(cbind(VSRS_OtherCountry_Infl,   VSRS_Others_Infl),   file = '~/VSRS_OtherCountry_Infl.csv',   row.names = FALSE)


## (3) G-SIBs' connectedness (From、To、Net)
# load G-SIBS Index and name
GSIBs_name <- read.table('G-SIBs name.txt', header = T)
GSIBs_Index  <- read.table('G-SIBs order.txt', header = T)

RSRS_GSIB_From   <- matrix(0, nrow = windownum, ncol = 29); colnames(RSRS_GSIB_From)   = GSIBs_name$G.SIBs_name
RSRS_GSIB_To     <- matrix(0, nrow = windownum, ncol = 29); colnames(RSRS_GSIB_To)     = GSIBs_name$G.SIBs_name
RSRS_GSIB_Net    <- matrix(0, nrow = windownum, ncol = 29); colnames(RSRS_GSIB_Net)    = GSIBs_name$G.SIBs_name

VSRS_GSIB_From   <- matrix(0, nrow = windownum, ncol = 29); colnames(VSRS_GSIB_From)   = GSIBs_name$G.SIBs_name
VSRS_GSIB_To     <- matrix(0, nrow = windownum, ncol = 29); colnames(VSRS_GSIB_To)     = GSIBs_name$G.SIBs_name
VSRS_GSIB_Net    <- matrix(0, nrow = windownum, ncol = 29); colnames(VSRS_GSIB_Net)    = GSIBs_name$G.SIBs_name

for (k in 1:windownum) {
  for (i in 1:29) {
    j <- GSIBs_Index$G.SIBs_Index[i]
    # return spillover panel
    RSRS_GSIB_From[k,i]   <- sum(RSRS_Institution[j,,k]) - RSRS_Institution[j,j,k]
    RSRS_GSIB_To[k,i]     <- sum(RSRS_Institution[,j,k]) - RSRS_Institution[j,j,k]
    RSRS_GSIB_Net[k,i]    <- RSRS_GSIB_To[k,i] - RSRS_GSIB_From[k,i]
    
    # volatility spillover panel
    VSRS_GSIB_From[k,i] <- sum(VSRS_Institution[j,,k]) - VSRS_Institution[j,j,k]
    VSRS_GSIB_To[k,i]   <- sum(VSRS_Institution[,j,k]) - VSRS_Institution[j,j,k]
    VSRS_GSIB_Net[k,i]  <- VSRS_GSIB_To[k,i] - VSRS_GSIB_From[k,i]
  }
}

write.csv(RSRS_GSIB_From, file = '~/RSRS_GSIB_From.csv', row.names = FALSE)
write.csv(RSRS_GSIB_To,   file = '~/RSRS_GSIB_To.csv',   row.names = FALSE)
write.csv(RSRS_GSIB_Net,  file = '~/RSRS_GSIB_Net.csv',  row.names = FALSE)

write.csv(VSRS_GSIB_From, file = '~/VSRS_GSIB_From.csv', row.names = FALSE)
write.csv(VSRS_GSIB_To,   file = '~/VSRS_GSIB_To.csv', row.names = FALSE)
write.csv(VSRS_GSIB_Net,  file = '~/VSRS_GSIB_Net.csv', row.names = FALSE)

## 5. Recursive Estimate
# Recursive analysis selects a recursive sample according to the desired sample interval, and then performs a Recursive-sample estimate
# In this case, the code do not provide all different recursive-sample estimations; but it can be provided on demand.
# Some examples of recursive analysis is provided below

# 5.1. construct the recursive samples
# T1，20130318-20140318，[315:575,]
# T2，20140318-20150318，[575:835,]
# T3，20181217-20191216，[1811:2070,]
# T4，20191216-20221230，[2070:2860,]
# T5，20210224-20220224，[2379:2639,]
# T6，20220224-20230228，[2639:2902,]

RF_ReturnT1 <- RF_Return[315:575,]
RF_ReturnT2 <- RF_Return[575:835,]
RF_ReturnT3 <- RF_Return[1811:2070,]
RF_ReturnT4 <- RF_Return[2070:2860,]
RF_ReturnT5 <- RF_Return[2379:2639,]
RF_ReturnT6 <- RF_Return[2639:2902,]

RF_VolatilityT1 <- RF_Volatility[315:575,]
RF_VolatilityT2 <- RF_Volatility[575:835,]
RF_VolatilityT3 <- RF_Volatility[1811:2070,]
RF_VolatilityT4 <- RF_Volatility[2070:2860,]
RF_VolatilityT5 <- RF_Volatility[2379:2639,]
RF_VolatilityT6 <- RF_Volatility[2639:2902,]

# 5.2. performs Recursive-sample estimation 
big_var_est <- function(data) {
  library(BigVAR)
  Model1 = constructModel(as.matrix(data), p = 4, struct = "BasicEN", gran = c(50, 50), VARX = list(), verbose = F)
  Model1Results = cv.BigVAR(Model1)
}

# 5.2.1. Return spillover panel 
ReturnVAR1  <- big_var_est(RF_ReturnT1)
ReturnVAR2  <- big_var_est(RF_ReturnT2)
ReturnVAR3  <- big_var_est(RF_ReturnT3)
ReturnVAR4  <- big_var_est(RF_ReturnT4)
ReturnVAR5  <- big_var_est(RF_ReturnT5)
ReturnVAR6  <- big_var_est(RF_ReturnT6)

ReturnSpilloverRecursive1  <- spilloverDY12(ReturnVAR1,  n.ahead = 10, no.corr = F)
ReturnSpilloverRecursive2  <- spilloverDY12(ReturnVAR2,  n.ahead = 10, no.corr = F)
ReturnSpilloverRecursive3  <- spilloverDY12(ReturnVAR3,  n.ahead = 10, no.corr = F)
ReturnSpilloverRecursive4  <- spilloverDY12(ReturnVAR4,  n.ahead = 10, no.corr = F)
ReturnSpilloverRecursive5  <- spilloverDY12(ReturnVAR5,  n.ahead = 10, no.corr = F)
ReturnSpilloverRecursive6  <- spilloverDY12(ReturnVAR6,  n.ahead = 10, no.corr = F)

ReturnRA_Table1  <- ReturnSpilloverRecursive1$tables
ReturnRA_Table2  <- ReturnSpilloverRecursive2$tables
ReturnRA_Table3  <- ReturnSpilloverRecursive3$tables
ReturnRA_Table4  <- ReturnSpilloverRecursive4$tables
ReturnRA_Table5  <- ReturnSpilloverRecursive5$tables
ReturnRA_Table6  <- ReturnSpilloverRecursive6$tables

# 5.2.2. Volatility spillover panel 
VolatilityVAR1  <- big_var_est(RF_VolatilityT1)
VolatilityVAR2  <- big_var_est(RF_VolatilityT2)
VolatilityVAR3  <- big_var_est(RF_VolatilityT3)
VolatilityVAR4  <- big_var_est(RF_VolatilityT4)
VolatilityVAR5  <- big_var_est(RF_VolatilityT5)
VolatilityVAR6  <- big_var_est(RF_VolatilityT6)

VolatilitySpilloverRecursive1  <- spilloverDY12(VolatilityVAR1, n.ahead = 10, no.corr = F)
VolatilitySpilloverRecursive2  <- spilloverDY12(VolatilityVAR2, n.ahead = 10, no.corr = F)
VolatilitySpilloverRecursive3  <- spilloverDY12(VolatilityVAR3, n.ahead = 10, no.corr = F)
VolatilitySpilloverRecursive4  <- spilloverDY12(VolatilityVAR4, n.ahead = 10, no.corr = F)
VolatilitySpilloverRecursive5  <- spilloverDY12(VolatilityVAR5, n.ahead = 10, no.corr = F)
VolatilitySpilloverRecursive6  <- spilloverDY12(VolatilityVAR6, n.ahead = 10, no.corr = F)

VolatilityRA_Table1  <- VolatilitySpilloverRecursive1$tables
VolatilityRA_Table2  <- VolatilitySpilloverRecursive2$tables
VolatilityRA_Table3  <- VolatilitySpilloverRecursive3$tables
VolatilityRA_Table4  <- VolatilitySpilloverRecursive4$tables
VolatilityRA_Table5  <- VolatilitySpilloverRecursive5$tables
VolatilityRA_Table6  <- VolatilitySpilloverRecursive6$tables

# 5.2.3. save the results of recursive-sample estimation 
# write.csv(ReturnRA_Table1,  file = '~/ReturnRA_Table_1.csv')
# write.csv(ReturnRA_Table2,  file = '~/ReturnRA_Table_2.csv')
# write.csv(ReturnRA_Table3,  file = '~/ReturnRA_Table_3.csv')
# write.csv(ReturnRA_Table4,  file = '~/ReturnRA_Table_4.csv')
# write.csv(ReturnRA_Table5,  file = '~/ReturnRA_Table_5.csv')
# write.csv(ReturnRA_Table6,  file = '~/ReturnRA_Table_6.csv')

# write.csv(VolatilityRA_Table1,  file = '~/VolatilityRA_Table_1.csv')
# write.csv(VolatilityRA_Table2,  file = '~/VolatilityRA_Table_2.csv')
# write.csv(VolatilityRA_Table3,  file = '~/VolatilityRA_Table_3.csv')
# write.csv(VolatilityRA_Table4,  file = '~/VolatilityRA_Table_4.csv')
# write.csv(VolatilityRA_Table5,  file = '~/VolatilityRA_Table_5.csv')
# write.csv(VolatilityRA_Table6,  file = '~/VolatilityRA_Table_6.csv')

# 5.3. Block aggregation for the recursive-sample estimation 
# The block aggregation of recursive estimation is consistent with full sample analysis... 
# and rolling sample analysis, which will not be detailed here.



