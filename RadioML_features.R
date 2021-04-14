library("survival")
library("survminer")
library(ggpubr)
library(reshape2)
library(dplyr)
library(purrr)
library(optparse)
library(data.table)

option_list <- list(
  make_option(c("-p", "--prefix"), type="character", default=FALSE,
              help="make the prefix of the result, default is prefix"),
  make_option(c("-t", "--partial"), type = "integer", default=FALSE,
              help="the ratio of split the samples into train and test (integer: 1-9, default 5)"),
  make_option(c("-n", "--numperm"), type="integer", default=FALSE,
              help="the repeat times  (integer, default 100)"),
  make_option(c("-R", "--RData"), type="character", default=FALSE,
              help="the previous main outout RData")
			  
)

opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

prefix <- opt$prefix
partial1 <- opt$partial
num_perm <- opt$numperm
RData <- opt$RData

if(!num_perm) num_perm <- 100
if(prefix == FALSE ) prefix <- "prefix"

#argv <- commandArgs(TRUE)
#prefix <- argv[1]
#partial1 <- argv[2]
#partial1 <- as.numeric(partial1)
#files <- argv[3]
#prefix <- "test"


prefix1 <- paste0("perm", partial1)
# extract the feature and feature coe


#load( paste0(prefix, "_", num_perm, "_", partial1, "_training.RData") )
load(RData)

dim(matrix)
row_count <- nrow(matrix)
total_table_test <- NULL
total_table_p <- data.frame(0, nrow = row_count, ncol = 1)
total_table_coe <- data.frame(0, nrow = row_count, ncol = 1)
total_table_coe_test <- data.frame(0, nrow = row_count, ncol = 1)
total_table_HR <- data.frame(0, nrow = row_count, ncol = 1)


total_table_test <- cbind(total_table_test, table_test)
total_table_p <- cbind(total_table_p, table_p)
total_table_coe <- cbind(total_table_coe, table_coe)
total_table_coe_test <- cbind(total_table_coe_test, table_coe_test)
total_table_HR <- cbind(total_table_HR, table_HR)
rm(table_test)
rm(table_p)
rm(table_coe_test)
rm(table_coe)
rm(table_HR)

table_coe_test <- total_table_coe_test[,-c(1:3)]
table_coe <- total_table_coe[,-c(1:3)]
table_p <- total_table_p[, -c(1:3)]
table_HR <- total_table_HR[, -c(1:3)]
table_test <- total_table_test
table_name <- rownames(table_p)

sink("chose_threshlod.txt")
table(rowSums(table_coe_test))
sink()


print(table(rowSums(table_coe_test)))
max <- names(table(rowSums(table_coe_test)))
max <- as.numeric(max[length(max)])
max <- floor(max*0.01)*100
max <- 800
cat("choose the high frequency threshold, and meanwhile about several features(5 ~ 20 ), input the threshold(integer):\n")
max <- scan("stdin", n = 1, what = integer(0) )


 
feature_final <- rownames(table_coe_test[rowSums(table_coe_test)>=max,])
feature_final
while(length(feature_final) < 10) {
  max <- max - 10
  feature_final <- rownames(table_coe_test[rowSums(table_coe_test)>=max,])
  feature_final
}

feature_coe_final <- apply(table_coe[feature_final,which(table_test[1,] <= 0.05)],1,mean)
write.table(data.frame(feature_coe_final), paste0(prefix1, ".", max, ".feature_coe.txt"), sep = "\t", quote = F)

data <- cbind(table_name, table_p, table_HR)
colnames(data) <- c("feature", "pvalue", "HR")
data$chose <- 1
data$chose[which(data$feature %in% feature_final)] <- 2

freq <- rowSums(table_coe_test)                    
freq <- data.frame(freq)                                                                                                                                                                                                                                                                  
rownames(freq) <- rownames(table_coe_test)
freq$name <- rownames(freq)
data <- cbind(data, freq)                                                               

pdf(paste0(prefix1, ".", max, ".HR_pvalue.pdf"), width = 6, height = 6)

ggplot(data, aes(x = pvalue, y = as.numeric(as.character(freq)), color = as.factor(data$chose)) ) + geom_point(size = data$chose*1) + theme_classic() + scale_color_manual(name = "features", values = c("black", "red"), labels = c("not chosed", "chosed")) + geom_hline(yintercept=1, linetype="dashed", color = "grey") + labs(fill = "features") + geom_hline(yintercept=max, linetype="dashed", color = "red")
ggplot(data, aes(x = pvalue, y = HR, color = as.factor(data$chose)) ) + geom_point(size = data$chose*1) + theme_classic() + scale_color_manual(name = "features", values = c("black", "red"), labels = c("not chosed", "chosed")) + ylim(c(0,10)) + geom_hline(yintercept=1, linetype="dashed", color = "grey") + labs(fill = "features")
ggplot(data, aes(x = pvalue, y = HR, color = as.factor(data$chose)) ) + geom_point(size = data$chose*1) + theme_classic() + scale_color_manual(name = "features", values = c("black", "red"), labels = c("not chosed", "chosed")) + ylim(c(0,10)) + geom_hline(yintercept=1, linetype="dashed", color = "grey") + labs(fill = "features") + geom_label( label="Hazard Ratio = 1",  x=1, y=1, label.padding = unit(0.55, "lines"), label.size = 0.35, color = "black" )
dev.off()


library(gplots)
table_p1 <- data.frame(table_p)
#table_p1 <- apply(table_p, 2, as.numeric)
table_p1[is.na(table_p1)] <- 1
threshold_count <- apply(table_p1, 1, function(x) unlist(table(x < 0.05)[2]))
threshold_count <- cbind(rownames(table_p1),threshold_count)
threshold_count <- threshold_count[complete.cases(threshold_count),]
#plot(density(as.numeric(threshold_count[,2])))
#threshold_count <- data.frame(threshold_count)
#colnames(threshold_count) <- c("feature", "freq")
pdf( paste0(prefix1, ".", max, ".threshold_count.pdf"), height=4,width=6)
plot(density(as.numeric(threshold_count[,2])))
ggplot(, aes(x=as.numeric(threshold_count[,2]))) + geom_density() + theme_classic() + xlim(c(0,1000)) + geom_vline(xintercept=max, linetype="dashed", color = "red") + geom_vline(xintercept=max, linetype="dashed", color = "red")
dev.off()

rm(list = ls())
