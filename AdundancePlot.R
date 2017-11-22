##troubleshooting abundance plot
#test <- read.csv("adunance_testfile.csv", header=T)
file <- "adunance_testfile.csv"

##This function inputs column of counts
          #each row being a unique clone sequence
          #output has two columns
          #first column in "value" which is a number of reads
          # second column is "count" which is numbe of unqiue clones
                   #that have that read number
          #v<- test[,3]
countUniq <- function(v) {
          b = aggregate(data.frame(count = v), list(value = v), length)
          
          return(b)
          
          
}



# args<-commandArgs(TRUE)
library(getopt)

spec <- matrix(c(
          'input'     , 'i', 1, "character", "input csv or tsv file (required)",
          'input2'     , 'j', 1, "character", "input csv or tsv file (optional)",
          'type'     , 't', 2, "character", "c or t (for csv or tsv) (optional); default: tsv",
          "freq1" , "m", 2, "double", "min freq to be considered as top clones",
          "freq2", "d", 2, "double", "min freq to be considered as detectable", 
          "fold", "f", 2, "integer", "min fold change to be considered as expanded", 
          "vj", 'v', 2, "integer", "[1/0] do analysis on V/J (default 1 means yes)",
          'help'   , 'h', 0, "logical",   "this help"
),ncol=5,byrow=T)

opt = getopt(spec);

#if input to help is not null or if no input
          #then show options and quit R session
if (!is.null(opt$help) || is.null(opt$input)) {
          cat(paste(getopt(spec, usage=T),"\n"));
          q();
}

#if no vj input, then set to 0
if (is.null(opt$vj)) {
          opt$vj = 0
}

file = opt$input
file2 = opt$input2

freq1 = 1e-4  # freq1 = 1e-4,  "top"
freq2 = 1e-5  # freq2 = 1e-5, "present"
fold = 5  # fold change default threshold

#if not set, type of file is tsv
itype = "t"
if (!is.null(opt$type)) {
          itype = opt$type
          
}

#if not set, freq1 is set as above
if(!is.null(opt$freq1)) {
          freq1 = opt$freq1
}

#same for freq 2 and fold
if(!is.null(opt$freq2)) {
          freq2 = opt$freq2
}

if(!is.null(opt$fold)) {
          fold = opt$fold
          
}

#if file type is comma separated do below
if( itype == "c")  {
          tcr = read.table(file, header=T, sep=",")
} else if (itype == "t") {
          tcr = read.table(file, header=T, sep="\t")                      
} else {
          tcr = read.table(file, header=T)                        
}

if(!is.null(opt$input2))  {
          if( itype == "c") {
                    tcr2 = read.table(file2, header=T, sep=",")
          } else if (itype == "t") {
                    tcr2 = read.table(file2, header=T, sep="\t")                      
          } else {
                    tcr2 = read.table(file2, header=T)                        
          }
          
          #merge TCR and TCR2
          #tcr <- read.table("D233/Pooled 1+2/D233-Spl-CD8-TEM-67-6ng.txt", header = TRUE)
          #tcr2 <- read.table("D233/Pooled 1+2/D233-Spl-CD8-TRM-67-6ng.txt", header = TRUE)

          file_1 <- tcr[,c(3,1)]
          print(head(file_1))
          file_2 <- tcr2[,c(3,1)]
          print(head(file_2))
          total <- merge(file_1,file_2,by="cdr3nt", all.x = TRUE, all.y = TRUE)
          total[is.na(total)] <- 0
          print(head(total))
          tcr <- total
}

print(tail(tcr))
n = ncol(tcr)
print(n)
mc = max(apply(tcr[,2:n],2, max))
for (i in 2:n) {
          mc[i-2] = mc[i-2] / sum(tcr[,i])
}

mf = max(mc, na.rm=TRUE)
fname = paste("abd", paste(file,file2), "pdf", sep=".")
pdf(fname, width=6, height=8)

flag = 0


for (i in 2:n) {
          cname = colnames(tcr)[i]
          s = tcr[,i]
          counts = countUniq(s)
          if (flag == 0) {
                    # plot
                    plot(counts[,1]/sum(tcr[,i]), counts[,2], log="xy", col=i, xlab="clone frequency", ylab="# of clones", main=paste(file,file2) , cex.lab=1.2, cex.axis=1.2, xlim=range(c(1e-7, 0.1)))
                    flag = 1
          } else {
                    points(counts[,1]/sum(tcr[,i]), counts[,2], col=i,)
          }
}
grid()
legend("topright", c(file,file2), col=c(2:n), pch=20)
dev.off()
