#This script calculates clonality on n number of files in a given directory
    # Each file has TCR sequencing data for one sample
           # Input is the directory with the files
           # Output is a list of dataframes with clonality calculations
                    # TCR files are .txt files 

#These are the functions used to calculate Clonality
rm(list=ls())
options(stringsAsFactors=FALSE)

          normalize <- function(data) {
  nc = ncol(data)
  for (i in 1:nc) {
    data[,i] = data[,i] / sum(data[,i])
}
  return(data)
}
          
          shannon.entropy <- function(p){
  if (min(p) < 0 || sum(p) <= 0)
    return(NA)
  p.norm <- p[p>0]/sum(p)
  -sum(log2(p.norm)*p.norm)
}
 
         Clonality <- function(p) {
  x = p[p>0] / sum(p)
  l = length(x)
  entropy = shannon.entropy(p)
  maxentropy = -log2(1/l)
  return(signif(1 - entropy / maxentropy, 3)) 
}
          
#Simpsons Index
         
         calcSI<-function(vals){
                   vals=vals[vals>0]
                   fq=vals/sum(vals)
                   si=sum(fq^2)
                   return(si)
         }
#R20
         calcr20 = function(X){
                   X=sort(X,decreasing=T)
                   X=X[X>0]
                   CX=cumsum(X)
                   num=length(which(CX/sum(X)<=0.2))
                   den=length(X)
                   return(num/den)
         }
#R50
         calcr50 = function(X){
                   X=sort(X,decreasing=T)
                   X=X[X>0]
                   CX=cumsum(X)
                   num=length(which(CX/sum(X)<=0.5))
                   den=length(X)
                   return(num/den)
         }
         
         
# Input your path to files here:
           File_path <- "/Users/michellemiron/Desktop/TCR data/All TCR data/All Reps pooled/"
    
# Get a list of all files in the directory    
           files <- list.files(path=File_path, pattern="*.txt")
           file <- files[[1]]
           file <- "D229_BM_CD4+CD69-rep1_2.txt"
# Function to calculate clonality on a given file
          # output is dataframe with clonality and file name
                     outputclonality_data <- function(file) {
                              file_location <- paste(File_path, file, sep = "")
                              all_data<-read.table(file_location, header=T )
                              counts <- all_data[,2]
                              countsdf <-as.data.frame(counts)
                              normalizedcounts <- normalize(countsdf)
                              entropy <- shannon.entropy(normalizedcounts)
                              clonalitycalc <- Clonality(normalizedcounts)
                              SI <-calcSI(counts)
                              R20<- calcr20(counts)
                              R50<- calcr50(counts)
                              Output <- data.frame(file,clonalitycalc,SI,R20,R50)
                              Output
                     }

                     
                     
                     
# Apply function to all files in a given directory
          data_compiled_list <- lapply(files, outputclonality_data)
          data_compiled_table <- my.matrix<-do.call("rbind", data_compiled_list)
          Path_save = "/Users/michellemiron/Desktop/TCR data/All TCR data/All Reps pooled/results/"
          file_output <- paste(Path_save,"clonality_R20_R50_SI.txt", sep = "")
          write.csv(data_compiled_table, file=file_output)
          
