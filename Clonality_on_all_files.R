#Purpose: Measure clonality from TCR sequencing data and output csv file with all clonality calculations
#Step 1: Process files into correct format:
  # a csv with 6 columns: VJgene, Count, VGeneName, JGeneName, Amino Acid, and Nucleotide sequence
    #Instructions for PC computer
      #Python script: process.py
      #Open powershell
        #$env:path="$env:Path;C:\Python27"
        #python C:\Users\mm4556\Desktop\process.py C:\Users\mm4556\Desktop\D299\BM-4-69neg.txt C:\Users\mm4556\Desktop\D299_processed\D299_BM_CD4_TEM.txt
    #Instructions for Mac computer
      #open terminal
      #python process.py Desktop/LN_CD8+CD69+.txt Desktop/processed.txt 
  # columns of output file are (1) vJ gene (2) count (3) vGene (4) jGene (5) Aminoacid (6) Nucleotide
#Step 2: Calculate metrics from count column of data
  #Functions
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
  #enter in data from saved location on P drive

    #do this on the desktop
      #move TCR files to desktop with folder for each donor with all files, reps and combined
        
    #edit file location for donor
  outputclonality_data <- function(file_name) {
      file_location <- "/Users/mm4556/Desktop/Master/D299/"
          #file_name <-"BM_CD4+CD69-.txt"
      all_data<-read.table(file = paste(file_location, file_name, sep = ""), header=T )
      counts <- all_data[,2]
      counts<-as.data.frame(counts)
      normalizedcounts <- normalize(counts)
      entropy <- shannon.entropy(normalizedcounts)
      clonalitycalc <- Clonality(normalizedcounts)
      clonality_output_data <- list(c(a=file_name, b=clonalitycalc))
      write.table(t(unlist(clonality_output_data)),"output_of_data_table")
      return(read.table("output_of_data_table"))
    }
      #outputclonality_data("D287_LN_CD8.txt")
    
    files <- list.files(path="/Users/mm4556/Desktop/Master/D299", pattern="*.txt")
    data_compiled_list <- lapply (files, outputclonality_data)
    data_compiled_table <- my.matrix<-do.call("rbind", data_compiled_list)
   file_output <- "/Users/mm4556/Desktop/Master/MasterClonalityD299.csv"
    write.csv(data_compiled_table, file=file_output)
  