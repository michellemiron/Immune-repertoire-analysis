
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
  return(signif(1 - entropy / maxentropy, 3)) }

    File_path <- "/Users/michellemiron/Desktop/TCR test data/Pooled 1 + 2/"
    
    files <- list.files(path=File_path, pattern="*.txt")
    file <- files[[1]]
    
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
      #outputclonality_data("D287_LN_CD8.txt")
    
    files <- list.files(path=File_path, pattern="*.txt")
    data_compiled_list <- lapply (files, outputclonality_data)
    data_compiled_table <- my.matrix<-do.call("rbind", data_compiled_list)
   file_output <- "/Users/mm4556/Desktop/Master/MasterClonalityD299.csv"
    write.csv(data_compiled_table, file=file_output)
  