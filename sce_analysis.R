##########################################################
######         Workflow for SCE analysis      #####
##########################################################
args= commandArgs(trailingOnly = T)
library(dplyr)
library(GenomicRanges)
library(tidyr)
library(ggplot2)
library(breakpointR)

#######################################################
####       Part 1: Collect breakpoints from BPR    ###
#######################################################

## Function for collect breakpoints from RData files
collectBreaksAllFiles <- function(datapath="path/to/BreakpointR_output/data/"){
  
  # List of full names of files
  files <- list.files(datapath,full.names=T)
  # Initialize empty dataframe
  breaks.all.files <- data.frame()
  
  n=1 # Set a counter
  
  for (file in files) {
    
    # Print counter statement and add one
    message("Reading ... " , basename(file), " ... ",round(  (n/length(files))*100  ,  digits = 1  ) , "%"  )
    n=n+1
    
    # Load just the breakpoints from RData files
    data <- get(load(file))[c('breaks','ID')]
    breakpoints <- as.data.frame(data$breaks)
    
    # Insert breakpoints with filename into initialized dataframe
    if (nrow(breakpoints)) {
      breakpoints$Library = data$ID
      breaks.all.files <- rbind(breakpoints,breaks.all.files)
    }  
  }
  
  # Return all breakpoints
  return(breaks.all.files)
}


# Running function using user-defined data path
breakpoints = collectBreaksAllFiles(datapath=args[1])

#######################################################
####       Part 2: Filtering breakpoints           ###
#######################################################

# 1) Remove homozygous strand-state transition events (SCEs can't occur in the same position on two homologs in the same cell)
breakpoints = filter(breakpoints,genoT!= "cc-ww" | genoT!= "cc-ww")



# 2) Filter out events that are too close to each other (2Mb), more likely to be background than real SCE
breakpoints$Library<-as.factor(breakpoints$Library) # Convert filenames to factor
breakpoints_2 = data.frame() # Initialize empty dataframe
n=1 # Start another counter for progress report
# Iterate through list of library names
for (level in levels(breakpoints$Library)){
  
  # Print progress counter
  message(  round(     (      (n/length(levels(breakpoints$Library)))   *100        )     ,digits = 2)   ,"% ... complete"   )
  # Create temporary dataframe with individual cell's breakpoints
  tmp = filter(breakpoints, Library==level)
  
  # Drop unused levels of chromosomes
  tmp$seqnames <- droplevels(tmp$seqnames)
  
  # Iterate through chromosomes that have an SCE
  for (level2 in levels(tmp$seqnames)){
    
    # Refilter dataframe for only SCEs in cell at hand and chromosome at hand
    tmp2 = filter(tmp, seqnames==level2)
    tmp3 <- GRanges(tmp2) # Convert to GRange
    
    overlaps = countOverlaps(tmp3, tmp3, type="any",maxgap = 2500000) # Identify events closer than 2.5Mb
    overlap_rows = which(overlaps %in% c(1))
    
    # Remove both events if closer than 2.5 Mb on the same chromosome
    tmp4 = tmp2[overlap_rows,]

    # Add all that remains to empty dataframe
    breakpoints_2 <- rbind(breakpoints_2,tmp4)
    
  }
  n=n+1
}



# 3) Filter out events near centromeres
centromeres <- read.table("centromeres2.txt",header=F,fill=T)
centromeres <-centromeres %>% dplyr::rename("seqnames"=V1,"start"=V2,"end"=V3)
# Convert centromere coordinates to Genomic Range object
centroGRange <- GRanges(centromeres) 
# Convert filtered breakpoints to Genomic Range object
summaryBreaks.df <- GRanges(breakpoints_2)

breakpoints_3 <-as.data.frame(summaryBreaks.df[-queryHits(findOverlaps(summaryBreaks.df, centroGRange, type="any")),])


#######################################################
####       Part 3: Expoprting putative SCEs        ###
#######################################################

# Write to file
write.table(breakpoints_3, "breakpoint_sces.bed", col.names = T, row.names = F, quote = F, sep="\t")



