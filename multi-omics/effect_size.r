#!/usr/bin/env Rscript

# -------------- NOTE: ---------------
# The last column should be the class.
# ------------------------------------

suppressPackageStartupMessages(library(vegan))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(psych))
parser <- ArgumentParser()

parser$add_argument("-i", help = "input affect file")
parser$add_argument("-d", help = "input affected file")
parser$add_argument("-m", default = "bray",help = "input adonis method default = bray")
parser$add_argument("-p", type="integer",default = 999,help = "input adonis permutations default = 999")
parser$add_argument("-o", help = "output file")
parser$add_argument("-a", help = "each variable adonis result")
parser$add_argument("-p_v",type = "double",default = 0.05,help = "cutoff filter adonis p_value")
parser$add_argument("-c_v",type = "double",default = 0.5,help = "cutoff pearson corr")
args <- parser$parse_args()

file1 <- file.path(args$i)
file2 <- file.path(args$d)
med <- file.path(args$m)
per <- file.path(args$p)
out_file <- file.path(args$o)
each_file <- file.path(args$a)
pvalue <- file.path(args$p_v)
cutoff_p <- file.path(args$c_v)


cat("\nUsing the following setting:\n")
cat("------------------------------------\n")
cat("input profile_test: ", file1, "\n")
cat("input file with the classified: ", file2, "\n")
cat("input file with the classified: ", out_file, "\n")
cat("input permutations: ", per, "\n")
cat("input method: ", med, "\n")
cat("------------------------------------\n")

#pvalue = 0.05
#cutoff_p = 0.5

library(vegan)
library(psych)
###########Get rid of the correlation > 0.5#########
corr_self <- function(x,y){
    name = x[,1]
    MGS1 = y[,name]
    corr_MGS_FM = corr.test(MGS1)
    a = corr_MGS_FM$r
    for(i in 1:nrow(a)){
        for(j in i:nrow(a)){
            a[j,i] = 0
        }
    }
    for(k in 1:nrow(x)){
        a = a[!a[k,]>cutoff_p,!a[k,]>cutoff_p]
        if(nrow(a)==k){
            break
        }
    }
    return(row.names(a))
}
##################################################################

############The selected variables are calculated adonis##########
calculta_R <- function(x,y){
    x <- x[order(x[,2],decreasing = T),]
    x <- x[x$P<=pvalue,]
    y_value <- y[,corr_self(x,y)]
    y_BM_adonis <- adonis(file2_data~.,data = y_value)
    y_BM_adonis
}

deal_NA <- function(x){
    for(i in 1:ncol(x)){
        x[,i][which(is.na(x[,i]))] = mean(x[,i],na.rm = T)
    }
    return(x)
}
################file###############################
#file1 <- "BM_CON.txt"
#file2 <- "FM_CON.txt"
#each_file <- "aaa" 

file1_data = read.table(file1,sep = "\t",row.names = 1,header = T,check.names = F,
                        stringsAsFactors = F)
file2_data = read.table(file2,sep = "\t",row.names = 1,header = T,check.names = F,
                        stringsAsFactors = F)

#############################
##The same row name#########
#############################
file2_data = file2_data[row.names(file2_data)%in%row.names(file1_data),]
file1_data = file1_data[row.names(file2_data),]

############################
##remove zero colname######
###########################
file1_data = file1_data[,apply(file1_data,2,sum)!=0]
file2_data = file2_data[,apply(file2_data,2,sum)!=0]

############################
###deal NA -> mean##########
###########################
file1_data = deal_NA(file1_data)
file2_data = deal_NA(file2_data)



each_adonis = read.table(each_file,sep = '\t',check.names = F,
                          stringsAsFactors = F,header = T)#index_all_CKD_index.txt


adonis_value = calculta_R(each_adonis,file1_data)

adj_adonis_value = RsquareAdj(1-adonis_value[[1]][5][,1][length(adonis_value[[1]][5][,1])-1],
                                 nrow(file1_data),length(adonis_value[[1]][5][,1])-2)
adonis_all_result = adonis_value[[1]]

value_table = data.frame(matrix(rep(0,ncol(adonis_all_result)*2),nrow=2))
row.names(value_table) = c("effect_size","adj_effect_size")
colnames(value_table) = colnames(adonis_all_result)
value_table[1,5] = 1-adonis_value[[1]][5][,1][length(adonis_value[[1]][5][,1])-1]
value_table[2,5] = adj_adonis_value
adonis_all_result = rbind(adonis_all_result,value_table)

write.table(adonis_all_result,out_file,quote = F,sep = "\t")
