library(tidyverse)
library(here)

#loading unfiltered mutatations (data4)
load(here("tidied_data/data_unfiltered_mutations.RData"))

#loading filtered mutatations (data4c)
load(here("tidied_data/data_filtered_mutations.RData"))

#create unfiltered dataset
data_unfilt <- bind_rows(data4)

#create filtered dataset
data_filt <- bind_rows(data4c)

x <- data4[[34]]
head(x)

#################################################
#################FUNCTIONS TO USE###############

#############################################
#function to count number of fixed mutations
count_fixed <- function(x){
  nrow(x[x$frequency > 99,])
}

########################################################
#function to count number of fixed mutations per species
count_fixed_per_species <- function(x,species){
  nrow(x[x$frequency > 99 & x$mut_sp_origin == species,])
}

#########################################################
#function to count mutation categories

#find all types of mutations
mutation_categories <- unique(data_unfilt$mutation_category)
#drop the <NA>, the last value
mutation_categories[-length(mutation_categories)]
#function to find types of mutation categories
count_mut_cat <- function(x,mutation_category){
  nrow(x[x$mutation_category == as.character(mutation_category),])
}

#make species vector
species <- c("A", "B", "C", "D")

#function to calculate the origin species of mutations
count_mut_sp_origin <- function(x, species) {
  nrow(x[x$mut_sp_origin == species,])
}





data_unfilt %>%
  group_by(replicate,time,community) %>%
  nest() %>%
  mutate(
    time = gsub("(.*)","D\\1",time),
    replicate = gsub("ep_","",replicate),
    all_mutations = map(data,nrow),
    fixed_mutations = map(data,count_fixed),
    snp_nonsynon = map(.x = data,mutation_category = mutation_categories[1], .f = count_mut_cat),
    snp_synon = map(.x = data,mutation_category = mutation_categories[2], .f = count_mut_cat),
    larger_insertion = map(.x = data,mutation_category = mutation_categories[3], .f = count_mut_cat),
    small_indel = map(.x = data,mutation_category = mutation_categories[4], .f = count_mut_cat),
    snp_intergenic = map(.x = data,mutation_category = mutation_categories[5], .f = count_mut_cat),
    snp_nonsense = map(.x = data,mutation_category = mutation_categories[6], .f = count_mut_cat),
    snp_noncoding = map(.x = data,mutation_category = mutation_categories[7], .f = count_mut_cat),
    large_amplification = map(.x = data,mutation_category = mutation_categories[8], .f = count_mut_cat),
    large_substitution = map(.x = data,mutation_category = mutation_categories[9], .f = count_mut_cat),
    large_deletion = map(.x = data,mutation_category = mutation_categories[10], .f = count_mut_cat),
    snp_nonsynonymous_nonsynonymous = map(.x = data,mutation_category = mutation_categories[11], .f = count_mut_cat),
    snp_nonsynonymous_synonymous = map(.x = data,mutation_category = mutation_categories[12], .f = count_mut_cat),
    snp_synonymous_synonymous = map(.x = data,mutation_category = mutation_categories[13], .f = count_mut_cat),
    snp_synonymous_nonsynonymous  = map(.x = data,mutation_category = mutation_categories[14], .f = count_mut_cat),
    snp_synonymous_nonsense   = map(.x = data,mutation_category = mutation_categories[15], .f = count_mut_cat),
    mutations_from_A = map(.x = data, species = species[1], .f = count_mut_sp_origin),
    mutations_from_B = map(.x = data, species = species[2], .f = count_mut_sp_origin),
    mutations_from_C = map(.x = data, species = species[3], .f = count_mut_sp_origin),
    mutations_from_D = map(.x = data, species = species[4], .f = count_mut_sp_origin),
    fixed_mutations_from_A = map(.x = data, species = species[1], .f = count_fixed_per_species),
    fixed_mutations_from_B = map(.x = data, species = species[2], .f = count_fixed_per_species),
    fixed_mutations_from_C = map(.x = data, species = species[3], .f = count_fixed_per_species),
    fixed_mutations_from_D = map(.x = data, species = species[4], .f = count_fixed_per_species),
        ) %>%
unnest(cols = c(5:ncol(.))) %>%
select(-data) %>%
unite("Label",c("community", "replicate", "time"), sep = "_", remove = FALSE)

data_filt %>%
  group_by(replicate,time,community) %>%
  nest() %>%
  mutate(
    time = gsub("(.*)","D\\1",time),
    replicate = gsub("ep_","",replicate),
    all_mutations = map(data,nrow),
    fixed_mutations = map(data,count_fixed),
    snp_nonsynon = map(.x = data,mutation_category = mutation_categories[1], .f = count_mut_cat),
    snp_synon = map(.x = data,mutation_category = mutation_categories[2], .f = count_mut_cat),
    larger_insertion = map(.x = data,mutation_category = mutation_categories[3], .f = count_mut_cat),
    small_indel = map(.x = data,mutation_category = mutation_categories[4], .f = count_mut_cat),
    snp_intergenic = map(.x = data,mutation_category = mutation_categories[5], .f = count_mut_cat),
    snp_nonsense = map(.x = data,mutation_category = mutation_categories[6], .f = count_mut_cat),
    snp_noncoding = map(.x = data,mutation_category = mutation_categories[7], .f = count_mut_cat),
    large_amplification = map(.x = data,mutation_category = mutation_categories[8], .f = count_mut_cat),
    large_substitution = map(.x = data,mutation_category = mutation_categories[9], .f = count_mut_cat),
    large_deletion = map(.x = data,mutation_category = mutation_categories[10], .f = count_mut_cat),
    snp_nonsynonymous_nonsynonymous = map(.x = data,mutation_category = mutation_categories[11], .f = count_mut_cat),
    snp_nonsynonymous_synonymous = map(.x = data,mutation_category = mutation_categories[12], .f = count_mut_cat),
    snp_synonymous_synonymous = map(.x = data,mutation_category = mutation_categories[13], .f = count_mut_cat),
    snp_synonymous_nonsynonymous  = map(.x = data,mutation_category = mutation_categories[14], .f = count_mut_cat),
    snp_synonymous_nonsense   = map(.x = data,mutation_category = mutation_categories[15], .f = count_mut_cat),
    mutations_from_A = map(.x = data, species = species[1], .f = count_mut_sp_origin),
    mutations_from_B = map(.x = data, species = species[2], .f = count_mut_sp_origin),
    mutations_from_C = map(.x = data, species = species[3], .f = count_mut_sp_origin),
    mutations_from_D = map(.x = data, species = species[4], .f = count_mut_sp_origin),
    fixed_mutations_from_A = map(.x = data, species = species[1], .f = count_fixed_per_species),
    fixed_mutations_from_B = map(.x = data, species = species[2], .f = count_fixed_per_species),
    fixed_mutations_from_C = map(.x = data, species = species[3], .f = count_fixed_per_species),
    fixed_mutations_from_D = map(.x = data, species = species[4], .f = count_fixed_per_species),
  ) %>%
  unnest(cols = c(5:ncol(.))) %>%
  select(-data) %>%
  unite("Label",c("community", "replicate", "time"), sep = "_", remove = FALSE)

x %>%
  filter(community == "ABCD")

gsub("ep_","",data_filt$replicate)

#merge mutational counts with coverage data
#load DNA sequencing values
load(here("coverage_data/data_cov.RData"))

head(data_cov)
     
as_tibble(data_cov) %>%
  filter(community == "ABCD")
