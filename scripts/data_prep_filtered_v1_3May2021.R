#####################################################
###by: C. Mallon
###purpose: load breseq data into R and tidy data
#####################################################

#load packages
library(tidyverse)
library(here)

#get folder where data is located
folder <- here("summary_txt_files")

#upload all breseq files. list.files() produces a character vector
file_list <- list.files(path=folder,full.names=TRUE, pattern="*.txt")

#read.delim() reads all files into a single list. This takes a long time--don't unnessecarily execute
#read.delim() also looks for those files in the working directory, all we need to do is supply the file_list variable we made.
data0 <- lapply(file_list,read.delim, header=TRUE, sep="\t")

#clean up file_list
file_list1 <- sub("_breseq_table_v3.txt","",file_list) 
file_list2 <- sub("D:/Google_Drive/Chemostat/data/breseq_filtered_run_v1/mutational_analysis/summary_txt_files/output_","",file_list1)


#change names in file list to account for contaminated reps
#file_list is a vector without names, so you have to change 
#the value of the character vector
#Only do this after you've already uploaded all files!
#(otherwise they can't upload)
#Also, below is what the changes should be and were in the origional data prep folder

#[17] "AB_r1_D14_breseq_table_v2.txt" <- "ABD_r4_D14_breseq_table_v2.txt"
#[18] "AB_r1_D21_breseq_table_v2.txt" <- "ABD_r4_D21_breseq_table_v2.txt"
#[19] "AB_r1_D31_breseq_table_v2.txt" <- "ABD_r4_D31_breseq_table_v2.txt"
#[20] "AB_r1_D7_breseq_table_v2.txt"  <- "ABD_r4_D7_breseq_table_v2.txt"

#[26] "AB_R3_D14_breseq_table_v2.txt" <- "ABD_R5_D14_breseq_table_v2.txt"
#[27] "AB_R3_D21_breseq_table_v2.txt" <- "ABD_R5_D21_breseq_table_v2.txt"
#[28] "AB_R3_D28_breseq_table_v2.txt" <- "ABD_R5_D28_breseq_table_v2.txt"
#[29] "AB_R3_D31_breseq_table_v2.txt" <- "ABD_R5_D31_breseq_table_v2.txt"
#[30] "AB_R3_D7_breseq_table_v2.txt"  <- "ABD_R5_D7_breseq_table_v2.txt"

#[36] "ABC_r2_D14_breseq_table_v2.txt" <- "ABD_r6_D14_breseq_table_v2.txt"
#[37] "ABC_r2_D21_breseq_table_v2.txt" <- "ABD_r6_D21_breseq_table_v2.txt"
#[38] "ABC_r2_D28_breseq_table_v2.txt" <- "ABD_r6_D28_breseq_table_v2.txt"
#[39] "ABC_r2_D31_breseq_table_v2.txt" <- "ABD_r6_D31_breseq_table_v2.txt"
#[40] "ABC_r2_D7_breseq_table_v2.txt"  <- "ABD_r6_D7_breseq_table_v2.txt"

#[41] "ABC_R3_D14_breseq_table_v2.txt" <- "ABD_R7_D14_breseq_table_v2.txt"
#[42] "ABC_R3_D21_breseq_table_v2.txt" <- "ABD_R7_D21_breseq_table_v2.txt"
#[43] "ABC_R3_D28_breseq_table_v2.txt" <- "ABD_R7_D28_breseq_table_v2.txt"
#[44] "ABC_R3_D31_breseq_table_v2.txt" <- "ABD_R7_D31_breseq_table_v2.txt"
#[45] "ABC_R3_D7_breseq_table_v2.txt"  <- "ABD_R7_D7_breseq_table_v2.txt"


file_list2[17] <- "ABD_r4_D14"
file_list2[18] <- "ABD_r4_D21"
file_list2[19] <- "ABD_r4_D31"
file_list2[20] <- "ABD_r4_D7" #note there is no AB_R1_D28

file_list2[27] <- "ABD_R5_D14"
file_list2[28] <- "ABD_R5_D21"
file_list2[29] <- "ABD_R5_D28"
file_list2[30] <- "ABD_R5_D31"
file_list2[31] <- "ABD_R5_D7"

file_list2[37] <- "ABD_r6_D14"
file_list2[38] <- "ABD_r6_D21"
file_list2[39] <- "ABD_r6_D28"
file_list2[40] <- "ABD_r6_D31"
file_list2[41] <- "ABD_r6_D7"

file_list2[42] <- "ABD_R7_D14"
file_list2[43] <- "ABD_R7_D21"
file_list2[44] <- "ABD_R7_D28"
file_list2[45] <- "ABD_R7_D31"
file_list2[45] <- "ABD_R7_D7"


#file_list without ancestral species, needed for downstream analysis
file_list3 <- file_list2[-c(1,122,178,202)]

#set the names for each dataframe in the list
data0 <- setNames(data0,file_list2)

#check data names
names(data0)

#add a column in each dataframe with the SampleID 
#the sub() function is just a substitue function,
#where I use part of the name in file_list to create a SampleID 
data1 <- mapply(cbind,data0,"SampleID"=file_list2,SIMPLIFY=F)

#make new colums for IDs
data2 <- data1
data2<- mapply(cbind,data2,"replicate"=NA, "time"=NA, "community" = NA, "SpRich"=NA, "mut_sp_origin"=NA, "position_scaffold"=NA, SIMPLIFY = FALSE)

#rename each df in list
data3 <- data2

#add rep number
for(i in seq_along(data3)){
  data3[[i]]$replicate[grepl("R1",data3[[i]]$SampleID,ignore.case=T)] <- "Rep_1"
  data3[[i]]$replicate[grepl("R2",data3[[i]]$SampleID,ignore.case=T)] <- "Rep_2"
  data3[[i]]$replicate[grepl("R3",data3[[i]]$SampleID,ignore.case=T)] <- "Rep_3"
  data3[[i]]$replicate[grepl("R4",data3[[i]]$SampleID,ignore.case=T)] <- "Rep_4"
  data3[[i]]$replicate[grepl("R5",data3[[i]]$SampleID,ignore.case=T)] <- "Rep_5"
  data3[[i]]$replicate[grepl("R6",data3[[i]]$SampleID,ignore.case=T)] <- "Rep_6"
  data3[[i]]$replicate[grepl("R7",data3[[i]]$SampleID,ignore.case=T)] <- "Rep_7"
}

#add day
for(i in seq_along(data3)){
  data3[[i]]$time[grepl("D7",data3[[i]]$SampleID,ignore.case=T)] <- 7
  data3[[i]]$time[grepl("D13",data3[[i]]$SampleID,ignore.case=T)] <- 14
  data3[[i]]$time[grepl("D14",data3[[i]]$SampleID,ignore.case=T)] <- 14
  data3[[i]]$time[grepl("D21",data3[[i]]$SampleID,ignore.case=T)] <- 21
  data3[[i]]$time[grepl("D28",data3[[i]]$SampleID,ignore.case=T)] <- 28
  data3[[i]]$time[grepl("D31",data3[[i]]$SampleID,ignore.case=T)] <- 31
}

#add community name
for(i in seq_along(data3)){
  data3[[i]]$community[grepl("ABCD_",data3[[i]]$SampleID, fixed = TRUE)] <- "ABCD"
  data3[[i]]$community[grepl("A_",data3[[i]]$SampleID, fixed=TRUE)] <- "A"
  data3[[i]]$community[grepl("B_",data3[[i]]$SampleID, fixed=TRUE)] <- "B"
  data3[[i]]$community[grepl("C_",data3[[i]]$SampleID, fixed=TRUE)] <- "C"
  data3[[i]]$community[grepl("D_",data3[[i]]$SampleID, fixed=TRUE)] <- "D"
  data3[[i]]$community[grepl("AB_",data3[[i]]$SampleID, fixed=TRUE)] <- "AB"
  data3[[i]]$community[grepl("AC_",data3[[i]]$SampleID, fixed=TRUE)] <- "AC"
  data3[[i]]$community[grepl("AD_",data3[[i]]$SampleID, fixed=TRUE)] <- "AD"
  data3[[i]]$community[grepl("BC_",data3[[i]]$SampleID, fixed=TRUE)] <- "BC"
  data3[[i]]$community[grepl("BD_",data3[[i]]$SampleID, fixed=TRUE)] <- "BD"
  data3[[i]]$community[grepl("CD_",data3[[i]]$SampleID, fixed=TRUE)] <- "CD"
  data3[[i]]$community[grepl("ABC_",data3[[i]]$SampleID, fixed = TRUE)] <- "ABC"
  data3[[i]]$community[grepl("ACD_",data3[[i]]$SampleID, fixed = TRUE)] <- "ACD"
  data3[[i]]$community[grepl("BCD_",data3[[i]]$SampleID, fixed = TRUE)] <- "BCD"
  data3[[i]]$community[grepl("ABD_",data3[[i]]$SampleID, fixed = TRUE)] <- "ABD"
  data3[[i]]$community[grepl("ABCD_",data3[[i]]$SampleID, fixed = TRUE)] <- "ABCD"
}

#add richness value
for(i in seq_along(data3)){
  data3[[i]]$SpRich <- ifelse(nchar(data3[[i]]$community)==1, 1,
                              ifelse(nchar(data3[[i]]$community)==2,2,
                                     ifelse(nchar(data3[[i]]$community)==3,3,
                                            ifelse(nchar(data3[[i]]$community)==4,4,NA
                                            )
                                     )
                              )
  )
}

#add species origin of mutation
for(i in seq_along(data3)){
  data3[[i]]$mut_sp_origin[grepl("A",data3[[i]]$seq_id, fixed = TRUE,ignore.case = FALSE)] <- "A"
  data3[[i]]$mut_sp_origin[grepl("B",data3[[i]]$seq_id, fixed = TRUE,ignore.case = FALSE)] <- "B"
  data3[[i]]$mut_sp_origin[grepl("C",data3[[i]]$seq_id, fixed = TRUE,ignore.case = FALSE)] <- "C"
  data3[[i]]$mut_sp_origin[grepl("D",data3[[i]]$seq_id, fixed = TRUE,ignore.case = FALSE)] <- "D"
}

#add position_scaffold. This column allows unique labels for positions
#of mutations since somtimes mutations of different scaffolds have the same position,
#making unnesseary repeats when plotting data
for(i in seq_along(data3)){
  data3[[i]]$position_scaffold <- paste(data3[[i]]$seq_id,data3[[i]]$position, sep=",")
}

#remove rows with NA in coverage because they are artifacts
#of species not in the community but included in the reference scaffold
data4 <- list()
for(i in seq_along(data3)){
  x <- data3[[i]][complete.cases(data3[[i]][,6]),]
  data4[[i]]<- x
}


#remove rows that are due to variation from the ancestor
#first make dataframe with all mutations present in ancestors
#then use anti_join to filter them out
ansc_df <- rbind(data4[[1]],data4[[122]],data4[[178]],data4[[202]])
data4a <- list()
for (i in seq_along(data4)){
  x <-dplyr::anti_join(data4[[i]],ansc_df,by="position_scaffold")
  data4a[[i]] <- x
} 

#rename df
data4b <- data4a

data4b <- setNames(data4b,file_list2)

#remove ancestral dataframes
data4b$"A-Ancestral_w_filters_v1"<-NULL
data4b$"B-Ancestral_w_filters_v1"<-NULL
data4b$"C-Ancestral_w_filters_v1"<-NULL
data4b$"D-Ancestral_w_filters_v1"<-NULL

#rename data
data4c <- data4b

#remove contaminated samples
data4c <- data4c[-c(16:19,26:30,36:44)]

#########################################################
#tidy data where ancestral mutations are not filtered
#remove ancestral dataframes

#set names
data4 <- setNames(data4,file_list2)

#remove ancestral data
data4$"A-Ancestral_w_filters_v1"<-NULL
data4$"B-Ancestral_w_filters_v1"<-NULL
data4$"C-Ancestral_w_filters_v1"<-NULL
data4$"D-Ancestral_w_filters_v1"<-NULL

save(data4,file = "data_unfiltered_mutations.RData")
save(data4c,file = "data_filtered_mutations.RData")
