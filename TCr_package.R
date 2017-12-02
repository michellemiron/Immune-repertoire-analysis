#Use tcR package for analysis of TCR sequencing data

library("tcR")

filepath<-
        
Location <- "P:\\CCTI_USERS/Michelle  Miron/TCR paper/Analysis and Data/For Michelle/D255/Pooled All Reps/"
List_of_files<- list.files(Location)
listofdf <- list()

for (i in 1:length(List_of_files)){
        filepath <- paste(Location,List_of_files[i], sep="")
        data <- parse.cloneset(filepath,"cdr3nt","cdr3aa","count","count","v","j",
                       "d", "VEnd", "JStart",c("DStart","DEnd"),NA,NA
                       ,NA)
        name <- List_of_files[i]
        tmp <- list(data)
        listofdf[[name]] <- tmp
}


vis.top.proportions(listofdf)



#D229

List_of_files<- list.files("P:\\CCTI_USERS/Michelle  Miron/Lab experiments/TCR sequencing Penn/TCR data/d229/Pooled/")


file1 <- paste("P:\\CCTI_USERS/Michelle  Miron/Lab experiments/TCR sequencing Penn/TCR data/d229/Pooled/",List_of_files[1], sep="")
file2 <- paste("P:\\CCTI_USERS/Michelle  Miron/Lab experiments/TCR sequencing Penn/TCR data/d229/Pooled/",List_of_files[2], sep="")
file3 <- paste("P:\\CCTI_USERS/Michelle  Miron/Lab experiments/TCR sequencing Penn/TCR data/d229/Pooled/",List_of_files[3], sep="")
file4 <- paste("P:\\CCTI_USERS/Michelle  Miron/Lab experiments/TCR sequencing Penn/TCR data/d229/Pooled/",List_of_files[4], sep="")
file5 <- paste("P:\\CCTI_USERS/Michelle  Miron/Lab experiments/TCR sequencing Penn/TCR data/d229/Pooled/",List_of_files[5], sep="")
file6 <- paste("P:\\CCTI_USERS/Michelle  Miron/Lab experiments/TCR sequencing Penn/TCR data/d229/Pooled/",List_of_files[6], sep="")
file7 <- paste("P:\\CCTI_USERS/Michelle  Miron/Lab experiments/TCR sequencing Penn/TCR data/d229/Pooled/",List_of_files[7], sep="")
file8 <- paste("P:\\CCTI_USERS/Michelle  Miron/Lab experiments/TCR sequencing Penn/TCR data/d229/Pooled/",List_of_files[8], sep="")


data<-read.table("P:\\CCTI_USERS/Michelle  Miron/Lab experiments/TCR sequencing Penn/TCR data/d229/Pooled/BM_CD4+CD69-.txt", header=T )

inputdata1 <- parse.cloneset(file1,"cdr3nt","cdr3aa","count","count","v","j",
                             "d", "VEnd", "JStart",c("DStart","DEnd"),NA,NA
                             ,NA)
inputdata2 <- parse.cloneset(file2,"cdr3nt","cdr3aa","count","count","v","j",
                             "d", "VEnd", "JStart",c("DStart","DEnd"),NA,NA
                             ,NA)
inputdata3 <- parse.cloneset(file3,"cdr3nt","cdr3aa","count","count","v","j",
                             "d", "VEnd", "JStart",c("DStart","DEnd"),NA,NA
                             ,NA)
inputdata4 <- parse.cloneset(file4,"cdr3nt","cdr3aa","count","count","v","j",
                             "d", "VEnd", "JStart",c("DStart","DEnd"),NA,NA
                             ,NA)
inputdata5 <- parse.cloneset(file5,"cdr3nt","cdr3aa","count","count","v","j",
                             "d", "VEnd", "JStart",c("DStart","DEnd"),NA,NA
                             ,NA)
inputdata6 <- parse.cloneset(file6,"cdr3nt","cdr3aa","count","count","v","j",
                             "d", "VEnd", "JStart",c("DStart","DEnd"),NA,NA
                             ,NA)
inputdata7 <- parse.cloneset(file7,"cdr3nt","cdr3aa","count","count","v","j",
                             "d", "VEnd", "JStart",c("DStart","DEnd"),NA,NA
                             ,NA)
inputdata8 <- parse.cloneset(file8,"cdr3nt","cdr3aa","count","count","v","j",
                             "d", "VEnd", "JStart",c("DStart","DEnd"),NA,NA
                             ,NA)

mylist <- list(inputdata1,inputdata2,inputdata3,inputdata4,inputdata5,inputdata6,inputdata7,inputdata8)
vis.top.proportions(mylist)



#D280

List_of_files<- list.files("P:\\CCTI_USERS/Michelle  Miron/Lab experiments/TCR sequencing Penn/TCR data/D280/Pooled")


file1 <- paste("P:\\CCTI_USERS/Michelle  Miron/Lab experiments/TCR sequencing Penn/TCR data/D280/Pooled/",List_of_files[1], sep="")
file2 <- paste("P:\\CCTI_USERS/Michelle  Miron/Lab experiments/TCR sequencing Penn/TCR data/D280/Pooled/",List_of_files[2], sep="")
file3 <- paste("P:\\CCTI_USERS/Michelle  Miron/Lab experiments/TCR sequencing Penn/TCR data/D280/Pooled/",List_of_files[3], sep="")
file4 <- paste("P:\\CCTI_USERS/Michelle  Miron/Lab experiments/TCR sequencing Penn/TCR data/D280/Pooled/",List_of_files[4], sep="")
file5 <- paste("P:\\CCTI_USERS/Michelle  Miron/Lab experiments/TCR sequencing Penn/TCR data/D280/Pooled/",List_of_files[5], sep="")
file6 <- paste("P:\\CCTI_USERS/Michelle  Miron/Lab experiments/TCR sequencing Penn/TCR data/D280/Pooled/",List_of_files[6], sep="")
file7 <- paste("P:\\CCTI_USERS/Michelle  Miron/Lab experiments/TCR sequencing Penn/TCR data/D280/Pooled/",List_of_files[7], sep="")
file8 <- paste("P:\\CCTI_USERS/Michelle  Miron/Lab experiments/TCR sequencing Penn/TCR data/D280/Pooled/",List_of_files[8], sep="")


data<-read.table("P:\\CCTI_USERS/Michelle  Miron/Lab experiments/TCR sequencing Penn/TCR data/D280/Pooled/BM_CD4+CD69-.txt", header=T )

inputdata1 <- parse.cloneset(file1,"cdr3nt","cdr3aa","count","count","v","j",
                             "d", "VEnd", "JStart",c("DStart","DEnd"),NA,NA
                             ,NA)
inputdata2 <- parse.cloneset(file2,"cdr3nt","cdr3aa","count","count","v","j",
                             "d", "VEnd", "JStart",c("DStart","DEnd"),NA,NA
                             ,NA)
inputdata3 <- parse.cloneset(file3,"cdr3nt","cdr3aa","count","count","v","j",
                             "d", "VEnd", "JStart",c("DStart","DEnd"),NA,NA
                             ,NA)
inputdata4 <- parse.cloneset(file4,"cdr3nt","cdr3aa","count","count","v","j",
                             "d", "VEnd", "JStart",c("DStart","DEnd"),NA,NA
                             ,NA)
inputdata5 <- parse.cloneset(file5,"cdr3nt","cdr3aa","count","count","v","j",
                             "d", "VEnd", "JStart",c("DStart","DEnd"),NA,NA
                             ,NA)
inputdata6 <- parse.cloneset(file6,"cdr3nt","cdr3aa","count","count","v","j",
                             "d", "VEnd", "JStart",c("DStart","DEnd"),NA,NA
                             ,NA)
inputdata7 <- parse.cloneset(file7,"cdr3nt","cdr3aa","count","count","v","j",
                             "d", "VEnd", "JStart",c("DStart","DEnd"),NA,NA
                             ,NA)
inputdata8 <- parse.cloneset(file8,"cdr3nt","cdr3aa","count","count","v","j",
                             "d", "VEnd", "JStart",c("DStart","DEnd"),NA,NA
                             ,NA)

mylist <- list(inputdata1,inputdata2,inputdata3,inputdata4,inputdata5,inputdata6,inputdata7,inputdata8)
vis.top.proportions(mylist)



#D287

List_of_files<- list.files("P:\\CCTI_USERS/Michelle  Miron/Lab experiments/TCR sequencing Penn/TCR data/D287/Pooled/")

file1 <- paste("P:\\CCTI_USERS/Michelle  Miron/Lab experiments/TCR sequencing Penn/TCR data/D287/Pooled/",List_of_files[1], sep="")
file2 <- paste("P:\\CCTI_USERS/Michelle  Miron/Lab experiments/TCR sequencing Penn/TCR data/D287/Pooled/",List_of_files[2], sep="")
file3 <- paste("P:\\CCTI_USERS/Michelle  Miron/Lab experiments/TCR sequencing Penn/TCR data/D287/Pooled/",List_of_files[3], sep="")
file4 <- paste("P:\\CCTI_USERS/Michelle  Miron/Lab experiments/TCR sequencing Penn/TCR data/D287/Pooled/",List_of_files[4], sep="")
file5 <- paste("P:\\CCTI_USERS/Michelle  Miron/Lab experiments/TCR sequencing Penn/TCR data/D287/Pooled/",List_of_files[5], sep="")
file6 <- paste("P:\\CCTI_USERS/Michelle  Miron/Lab experiments/TCR sequencing Penn/TCR data/D287/Pooled/",List_of_files[6], sep="")
file7 <- paste("P:\\CCTI_USERS/Michelle  Miron/Lab experiments/TCR sequencing Penn/TCR data/D287/Pooled/",List_of_files[7], sep="")
file8 <- paste("P:\\CCTI_USERS/Michelle  Miron/Lab experiments/TCR sequencing Penn/TCR data/D287/Pooled/",List_of_files[8], sep="")

inputdata1 <- parse.cloneset(file1,"cdr3nt","cdr3aa","count","count","v","j",
                                  "d", "VEnd", "JStart",c("DStart","DEnd"),NA,NA
                                  ,NA)
inputdata2 <- parse.cloneset(file2,"cdr3nt","cdr3aa","count","count","v","j",
                            "d", "VEnd", "JStart",c("DStart","DEnd"),NA,NA
                            ,NA)
inputdata3 <- parse.cloneset(file3,"cdr3nt","cdr3aa","count","count","v","j",
                            "d", "VEnd", "JStart",c("DStart","DEnd"),NA,NA
                            ,NA)
inputdata4 <- parse.cloneset(file4,"cdr3nt","cdr3aa","count","count","v","j",
                            "d", "VEnd", "JStart",c("DStart","DEnd"),NA,NA
                            ,NA)
inputdata5 <- parse.cloneset(file5,"cdr3nt","cdr3aa","count","count","v","j",
                            "d", "VEnd", "JStart",c("DStart","DEnd"),NA,NA
                            ,NA)
inputdata6 <- parse.cloneset(file6,"cdr3nt","cdr3aa","count","count","v","j",
                            "d", "VEnd", "JStart",c("DStart","DEnd"),NA,NA
                            ,NA)
inputdata7 <- parse.cloneset(file7,"cdr3nt","cdr3aa","count","count","v","j",
                             "d", "VEnd", "JStart",c("DStart","DEnd"),NA,NA
                             ,NA)
inputdata8 <- parse.cloneset(file8,"cdr3nt","cdr3aa","count","count","v","j",
                             "d", "VEnd", "JStart",c("DStart","DEnd"),NA,NA
                             ,NA)

mylist <- list(inputdata1,inputdata2,inputdata3,inputdata4,inputdata5,inputdata6,inputdata7,inputdata8)
vis.top.proportions(mylist)

#enter in data from saved location on P drive
data<-read.table(file ="P:\\CCTI_USERS/Michelle  Miron/Lab experiments/TCR sequencing Penn/TCR data/D287/Pooled/16p040-Donna-53-287BM4-69.txt", header=T )
counts <- data[,2]
#keep counts as a dataframe
counts<-as.data.frame(counts)
normalizedcounts <- normalize(counts)
entropy <- shannon.entropy(normalizedcounts)
clonalitycalc <- Clonality(normalizedcounts)
clonalitycalc
entropy