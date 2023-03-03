#绘制时间-turnoverate图
#获取时间序列
sample <- sample_data(ADBiofilmIG)
time <- sample[,5]
time <- data.frame(time)
time$Date <- as.Date(time$Date, format="%Y-%m-%d")
time_matrix <- matrix(0, nrow=nrow(sample), ncol=nrow(sample))
rownames(time_matrix) <- colnames(time_matrix) <- rownames(sample)
for(i in 1:nrow(time)) {
  for(j in 1:nrow(time)) {
    time_matrix[i,j] <- as.numeric(difftime(time[i,1], time[j,1], units="days"))
  }
}
logtime_matrix <- log(abs(time_matrix))
#获取数据temporal turnover nestedness 并生成列表
ADBiofilmIG01 <- ifelse(t(otu_table(ADBiofilmIG)) > 1, 1, 0) #出现次数小于1都是0
row <- 1
ADBiofilmIGresult <- data.frame()#创建新dataframe
for(i in 1:(nrow(time)-1)) {
  for(j in (i+1):nrow(time)) {
    temp1 <- ADBiofilmIG01[i, , drop = FALSE]
    rownames(temp1) <- NULL
    temp2 <- ADBiofilmIG01[j, , drop = FALSE]
    rownames(temp2) <- NULL
    tempdata <- data.frame(beta.temp(temp1,temp2,index.family="sor"))#得到结果 1sne是nestedness 2sim是turnover
    new_row <- data.frame(nestedness = log(tempdata[1,2]),turnover = log(tempdata[1,1]),timediff = logtime_matrix[i,j]) # 将结果添加到新的数据框
    Rowname = paste(rownames(time)[i],rownames(time)[j], sep = "_")
    ADBiofilmIGresult <- rbind(ADBiofilmIGresult, new_row)
    rownames(ADBiofilmIGresult)[row] <- Rowname
    row <- row + 1
  }
}
#绘制散点图
p = ggplot(ADBiofilmIGresult, aes(x=timediff, y=turnover, color=type)) +
  geom_point(alpha=1, size=0.7) + geom_jitter()+
  labs(x="logtime", y="logturnover") +
  theme(axis.text=element_text(angle=45,vjust=1, hjust=1)) +
  # 拟合只有loess方法初始为曲线，lm方法默认为直线，但可以添加公式
  main_theme + geom_smooth(method = "lm", formula = y ~ poly(x,3)) # , formula = y ~ ns(x,3)