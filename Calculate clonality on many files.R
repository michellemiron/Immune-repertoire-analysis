#This script calculates clonality on n number of files in a given directory
    # Each file has TCR sequencing data for one sample
           # Input is the directory with the files
           # Output is a list of dataframes with clonality calculations
                    # TCR files are .txt files 

#These are the functions used to calculate Clonality
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
          
# Input your path to files here:
           File_path <- "/Users/michellemiron/Dropbox/For Michelle/D233/Pooled All Reps/Raw Data/"
    

# Function to calculate clonality on a given file
          # output is dataframe with clonality and file name
                     outputclonality_data <- function(file) {
                              file_location <- paste(File_path, file, sep = "")
                              all_data<-read.table(file_location, header=T )
                              counts <- all_data[,2]
                              counts<-as.data.frame(counts)
                              normalizedcounts <- normalize(counts)
                              entropy <- shannon.entropy(normalizedcounts)
                              clonalitycalc <- Clonality(normalizedcounts)
                              Output <- data.frame(file,clonalitycalc)
                              Output
                     }
# Get a list of all files in the directory    
          files <- list.files(path=File_path, pattern="*.txt")
          file <- files[[1]]
          
# Apply function to all files in a given directory
          data_compiled_list <- lapply (files, outputclonality_data)
          data_compiled_table <- my.matrix<-do.call("rbind", data_compiled_list)
          file_output <- paste(File_path,"clonality.txt", sep = "")
          write.csv(data_compiled_table, file=file_output)
