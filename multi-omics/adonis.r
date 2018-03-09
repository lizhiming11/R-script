
#!/usr/bin/env Rscript

# -------------- NOTE: ---------------
# The last column should be the class.
# ------------------------------------

suppressPackageStartupMessages(library(vegan))
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()

parser$add_argument("-i", help = "input profile_test")
parser$add_argument("-t", help = "input file with the classified")
parser$add_argument("-f", help = "input fold_change citoff,default:2",default=2)
parser$add_argument("-p", help = "input p value method(number),default:11",default=11)
parser$add_argument("-r", help = "input Whether to flip",default=FALSE)


args <- parser$parse_args()

infile <- file.path(args$i)
inclass <- file.path(args$t)
inFC <- file.path(args$f)
inp <- file.path(args$p)
inr <- file.path(args$r)

cat("\nUsing the following setting:\n")
cat("------------------------------------\n")
cat("input profile_test: ", infile, "\n")
cat("input file with the classified: ", inclass, "\n")
cat("------------------------------------\n")


adonis_self <- function(x,y){
  data1 = matrix(rep(0,ncol(x)*3),ncol = 3)
  data1[,1] = colnames(x)
  colnames(data1) = c("name","R2","P")
  for(i in 1:ncol(x)){
    a = adonis(y~x[,i])
    data1[i,2] = a$aov.tab$R2[1]
    data1[i,3] = a$aov.tab$`Pr(>F)`[1]
  }
  return(data1)
}

merge_self <- function(a,b,b1,c,c1){
  data1 = adonis_self(a,b)
  data2 = adonis_self(a,c)
  data5 = cbind(data1,data2[,c(2,3)])
  colnames(data5)[c(2,4)] = c(paste(b1,"R2",sep = "_"),
                                  paste(c1,"R2",sep = "_")
                                  )
  data5
}
clinical_procedure <- function(x){
    x <- x[,!apply(x,2,function(y){sum(is.na(y))/nrow(x)>0.5})]
    for(i in 1:ncol(x)){
        x[,i][which(is.na(x[,i]))] = mean(x[,i],na.rm = T)
    }
    return(x)
}




BM_file <- "BM_CON.txt"
FM_file <- "FM_CON.txt"
index_file <- "mgs_con.txt"#"index_ckd.txt"#
host_file <- "module_con.txt"#host_prop_ckd.txt
file1_name <- "mgs_adonis_CON.txt"#index_all_CKD_index.txt
file2_name <- "module_adonis_CON.txt"#host_prop_all_CKD_index.txt


BM = read.table(BM_file,sep = "\t",row.names = 1,header = T,check.names = F)
FM = read.table(FM_file,sep = "\t",row.names = 1,header = T,check.names = F)
index = read.table(index_file,sep = "\t",row.names = 1,header = T,check.names = F,
                   stringsAsFactors = F)
#index$`AST/ALT` = as.double(index$`AST/ALT`)
host = read.table(host_file,sep = "\t",row.names = 1,header = T,check.names = F,
                  stringsAsFactors = F)
host = host[row.names(BM),]
host = host[,apply(host,2,sum)!=0]

#index = clinical_procedure(index)
#host = clinical_procedure(host)

index_all = merge_self(index,FM,"FM",BM,"BM")
host_all = merge_self(host,FM,"FM",BM,"BM")

write.table(host_all,file2_name,quote = F,sep = "\t")
write.table(index_all,file1_name,quote = F,sep = "\t")
