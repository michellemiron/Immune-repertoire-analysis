#Use tcR package for analysis of TCR sequencing data

library("tcR")
        
#function to generate a list of df
Gen_listofdf <- function(List_of_files) {
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
          return(listofdf)
}

#input these paramters:
Location <- "/Users/michellemiron/Desktop/testing MiXCR/TCR data/"
##previous locations
#"/Users/michellemiron/Desktop/TCR data/All TCR data/All Reps pooled/"
List_of_files<- list.files(Location,pattern="txt")
CD8Tcells <- Gen_listofdf(List_of_files)
vis.top.proportions(CD8Tcells)

##second plot
Location <- "/Users/michellemiron/Desktop/testing MiXCR/TCR data/cd4 data/"
List_of_files<- list.files(Location,pattern="txt")
CD4Tcells <- Gen_listofdf(List_of_files)
vis.top.proportions(CD4Tcells)

clonal.proportion(listofdf, 25)
twb.space <- clonal.space.homeostasis(listofdf[[3]])
vis.clonal.space(twb.space)

