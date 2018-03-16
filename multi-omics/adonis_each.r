#!/usr/bin/env Rscript

# -------------- NOTE: ---------------
# The last column should be the class.
# ------------------------------------

suppressPackageStartupMessages(library(vegan))
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()

parser$add_argument("-i", help = "input affect file")
parser$add_argument("-d", help = "input affected file")
parser$add_argument("-m", default = "bray",help = "input adonis method default = bray")
parser$add_argument("-p", type="integer",default = 999,help = "input adonis permutations default = 999")
parser$add_argument("-o", help = "output file")

args <- parser$parse_args()

file1 <- file.path(args$i)
file2 <- file.path(args$d)
med <- file.path(args$m)
per <- file.path(args$p)
out_file <- file.path(args$o)


cat("\nUsing the following setting:\n")
cat("------------------------------------\n")
cat("input profile_test: ", file1, "\n")
cat("input file with the classified: ", file2, "\n")
cat("input file with the classified: ", out_file, "\n")
cat("input permutations: ", per, "\n")
cat("input method: ", med, "\n")
cat("------------------------------------\n")


adonis_self <- function(x,y){
    data1 = matrix(rep(0,ncol(x)*3),ncol = 3)
    data1[,1] = colnames(x)
    colnames(data1) = c("name","R2","P")
    for(i in 1:ncol(x)){
        a = adonis(y~x[,i],permutations = as.numeric(per),method = med)
        data1[i,2] = a$aov.tab$R2[1]
        data1[i,3] = a$aov.tab$`Pr(>F)`[1]
    }
    return(data1)
}

deal_NA <- function(x){
    for(i in 1:ncol(x)){
        x[,i][which(is.na(x[,i]))] = mean(x[,i],na.rm = T)
    }
    return(x)
}

#file1 <- "BM_CON.txt"
#file2 <- "FM_CON.txt"

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

###############################################
###Calculate the impact size of each variable.#
###############################################
adonis_data = adonis_self(file1_data,file2_data)


###########save##############################
write.table(adonis_data,out_file,quote = F,sep = "\t")

