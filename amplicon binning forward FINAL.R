############# FINAL
############# FORWARD
############# AMPLICON
############# BINNING
############# Sydney Grant

library(dplyr)
library(plyr)
library(hash)
library(rlist)
library(R.utils)

# load mutation dataset
# dataset must include chromosome, column name "chr"
# and base pair location, column name "pos"

data <- read.csv("E:\\102621\\Fowler_NSES_Dataset.csv") 



panel_length = 100000 #set panel and amplicon length
amp_length = 125


pos <- data$pos
pos_freq <- data.frame(ftable(pos)) #make frequency table for each position

data <- merge(data, pos_freq, by = "pos") #merge frequency table with original df
data <- unique(data) #keep only unique positions


## amplicon finder function:
## 1) Defines hotspot regions
## 2) Generates list of amplicons which can cover each hotspot regions

amplicon_finder <- function(data, chr){
  pos <- data$pos
  
  keys <- make.keys(data$pos)
  dict <- hash(keys = keys, values = data$Freq) #create dictionary for mutation locations and mutation frequency
  
  make_bins <- sort(unique(pos))
  bins <- data.frame()
  while (length(make_bins) > 1){   ## groups all mutations by relative position
    temp_bins <- make_bins         ## starts with first mutation, checks for next closes mutation
    i <- 1                         ## mutations are added to the group if they are less than one amplicon length away from neighboring mutation
    start <- make_bins[1]          ## if mutations are further than one amplicon away, new group is started
    p_dif <- 1
    go = "yes"
    while (go == "yes"){
      p1 <- temp_bins[i]
      p_next <- temp_bins[i + 1]
      i <- i + 1
      if (i <= length(make_bins)){
        p_dif <- p_next - p1
      }
      if (i > length(make_bins)){
        p_dif = amp_length + 100
        go = "no"
      }
      if (p_dif > amp_length){
        p_next <- p1
        go = "no"
      }
    }
    end <- p_next
    vec <- as.list(start:end)
    keep <- setdiff(make_bins, vec)
    make_bins <- keep
    chr = chr
    row <- data.frame("lowerbound" = start, "upperbound" = end, "chromosome" = chr)
    bins <- rbind(bins, row)
  }
  if (length(make_bins) == 1){
    row <- data.frame("lowerbound" = make_bins, "upperbound" = make_bins, "chromosome" = chr)  # all groups are added to dataframe of bins
    bins <- rbind(bins, row)
  }
  
  mut_locations <- data.frame()
  possible_bins <- data.frame()
  
  for (i in 1:nrow(bins)){             ## for each bin, finds amplicons which may optimally cover that region
    up_mut <- bins$upperbound[[i]]
    low_mut <- bins$lowerbound[[i]]   # find length of bin
    vec <- up_mut:low_mut
    dif <- length(vec) - amp_length
    if (dif <= 0){                    # if the bin is <= amplicon length, we only need to generate a single amplicon
      up <- bins$upperbound[[i]]
      low <- bins$lowerbound[[i]]
      vec <- up:low
      mutations <- intersect(vec, pos)    #find all mutations in total bin region
      mutations_keys <- make.keys(mutations)
      count <- sum(as.vector(values(dict, keys = mutations_keys))) #get total mutations for total bin region
      mut_all <- c()
      for (m in 1:length(mutations)){
        mut_all <- c(mut_all, rep(mutations[[m]], as.vector(values(dict, keys = mutations_keys))[[m]]))
      }
      weighted_mid <- round(mean(mut_all))
      weighted_mid_max <- (ceiling(mean(vec)) + (amp_length/2))        # the single amplicon we choose to cover this bin is weighted based on mutation distribution
      weighted_mid_min <- (floor(mean(vec)) - (amp_length/2))
      if (weighted_mid < weighted_mid_min){final_mid <- weighted_mid_min}
      if (weighted_mid > weighted_mid_max){final_mid <- weighted_mid_max}
      if (weighted_mid > weighted_mid_min & weighted_mid < weighted_mid_max){final_mid <- weighted_mid}
      up <- final_mid + ceiling((amp_length - 1) / 2)
      low <- final_mid - floor((amp_length - 1) / 2)
      vec <- up:low
      mutations <- intersect(vec, pos) #find all mutations in total bin region
      mutations_keys <- make.keys(mutations)
      count <- sum(as.vector(values(dict, keys = mutations_keys))) #get total mutations for total bin region
      row <- data.frame("lowerbound" = low, "upperbound" = up, 
                        "chromosome" = bins$chromosome[[i]], "count" = count,"id" = "x")
      possible_bins <- rbind(possible_bins, row)
    }
    else{                                    # for bins greater than one amplicon length, we generate 5 possible amplicons for each mutation within the bin
      up_mut <- bins$upperbound[[i]]
      low_mut <- bins$lowerbound[[i]]
      vec <- up_mut:low_mut
      mut_list <- intersect(vec, pos)
      id = paste(as.character(bins$chromosome[[i]]), as.character(i), sep = "-") # every time a new bin is created, it is assigned a unique id
      b_list <- list()
      for (mut in mut_list){                      # creating 5 amplicons per mutation
        bin1 <- (mut - amp_length):(mut - 1)
        b_list[[1]] <- bin1
        bin2 <- (mut - (amp_length - 1)):(mut)
        b_list[[2]] <- bin2
        bin3 <- (mut + amp_length):(mut + 1)
        b_list[[3]] <- bin3
        bin4 <- (mut + (amp_length - 1)):(mut)
        b_list[[4]] <- bin4
        bin5 <- (mut - ceiling(amp_length / 2)):(mut - floor(amp_length / 2))
        b_list[[5]] <- bin5
        
        for (b in b_list){
          low <- min(b)
          up <- max(b)
          mutations <- intersect(b, pos)
          if (length(mutations > 0)){
            mutations_keys <- make.keys(mutations)
            count <- sum(as.vector(values(dict, keys = mutations_keys))) #get total mutations for total bin region
            row <- data.frame("lowerbound" = low, "upperbound" = up, 
                              "chromosome" = bins$chromosome[[i]], "count" = count,"id" = id)
            possible_bins <- rbind(possible_bins, row)
          }
        }
      }
      
    }
  }
  return(possible_bins)       # all amplicons are added to dataframe
}




##################################################################################


chrom_list <- unique(as.list(data$chr)) # list of all chromosomes found in this dataset
all_pos_bins <- data.frame()

for (chrom in chrom_list) {
  chromosome <- subset(data, chr == chrom)
  chromsort <- chromosome[order(chromosome$pos), ] #order df from lowest to highest position
  chrom_bins <- amplicon_finder(chromsort, chrom) #run through above function
  all_pos_bins <- rbind(all_pos_bins, chrom_bins) #add all bins together from all chromosomes
}

###################################################################################
ordered_bin <- all_pos_bins[order(-all_pos_bins$count), ]

big_spots <- ordered_bin[!ordered_bin$id == "x", ]
unique_spots <- unique(as.list(big_spots$id))      # we then need to find the optimal combination of amplicons to cover each large hotspot

keys <- make.keys(data$pos)
dict <- hash(keys = keys, values = data$Freq)      # create a new dictionary of mutation positions and frequency of mutations
possible_bins <- data.frame()
for (spot in unique_spots){
  spot_sub <- subset(big_spots, id == spot)
  big_bin_all <- list()
  while (nrow(spot_sub) > 0){
    row <- data.frame()                            # for each hotspot, subset a dataframe with only amplicons which cover that hotspot
    new_sub <- data.frame()
    big_spot <- spot_sub[1,]
    big_bin <- big_spot$upperbound:big_spot$lowerbound
    big_bin_all <- c(big_bin_all, big_bin)
    possible_bins <- rbind(possible_bins, big_spot)   # first find the amplicon which contains the most mutations and add to final list
    if (nrow(spot_sub) > 1){
      for (k in 2:nrow(spot_sub)){
        row <- data.frame()
        test_bin <- spot_sub$upperbound[k]:spot_sub$lowerbound[k]
        overlap <- intersect(big_bin_all, test_bin)    # then for each remaining bin, adjust mutation count so they are only showing unique mutations
        if (length(overlap > 0)){                      # continue process until there are no amplicons left which would capture additional unique mutations
          unique_bin <- setdiff(test_bin, overlap)
          mutations <- intersect(unique_bin, pos) 
          if (length(mutations) > 0){
            mutations_keys <- make.keys(mutations)
            new_count <- sum(as.vector(values(dict, keys = mutations_keys))) 
            row <- data.frame("lowerbound" = min(test_bin), "upperbound" = max(test_bin), 
                              "chromosome" = spot_sub$chromosome[1], "count" = new_count,"id" = spot_sub$id[1])
          }
        }
        if (length(overlap) <= 0) {row <- spot_sub[k,]}
        new_sub <- rbind(new_sub, row)
      }
    }
    spot_sub <- new_sub
    if (nrow(spot_sub > 0)){spot_sub <- spot_sub[order(-spot_sub$count), ]}
  }
}


final_bin <- rbind(subset(ordered_bin, id == "x"), possible_bins)
final_bin <- final_bin[order(-final_bin$count), ]

# reorder to favor amplicons with mutations in center

cutoff_point <- final_bin$count[[round(panel_length / amp_length)]]

subset_final_bin <- subset(final_bin, count == cutoff_point)
dif_list <- list()
for (i in nrow(subset_final_bin)){
  whole_region <- subset_final_bin$lowerbound[[i]]:subset_final_bin$upperbound[[i]]
  bin_mutations <- intersect(pos, whole_region)
  mid_diff <- abs((max(whole_region) - max(bin_mutations)) - (min(bin_mutations) - min(whole_region)))
  dif_list <- c(dif_list, mid_diff)
}
subset_final_bin$dif <- unlist(dif_list)
subset_final_bin <- subset_final_bin[order(-subset_final_bin$dif), ]
subset_final_bin <- subset_final_bin[,-ncol(subset_final_bin)]

df_range_min <- min(which(final_bin$count == cutoff_point))
df_range_max <- max(which(final_bin$count == cutoff_point))

final_bin[df_range_min:df_range_max,] <- subset_final_bin[1:nrow(subset_final_bin),]

##############################################################


bin_len_list <- list()
for (y in 1:nrow(final_bin)) {      #calculating bin length using upper and lower bound positions
  upper = final_bin$upperbound[[y]]
  lower = final_bin$lowerbound[[y]]
  bin_len <- upper - lower + 1# add 1 to include both upper and lowerbound plus all bp inbetween
  bin_len_list <- c(bin_len_list, bin_len) #make list of all bin lengths
}
final_bin$bin_length <- unlist(bin_len_list) # add to dataframe


bin_len_list <- as.numeric(final_bin$bin_length) #new lists since dataframe has been reordered
count_list <- as.numeric(final_bin$count)


bin_len_cum <- list()
cum_mut_list <- list()
cum_mut_list[[1]] <- count_list[[1]]  #start list with first bin length and first mutation count
bin_len_cum[[1]] <- bin_len_list[[1]]

for (u in 2:length(count_list)) {     #repeat for all other bins
  cum_mut_list[[u]] <- cum_mut_list[[u-1]] + count_list[[u]]  #add all calculated bin lengths and mutation counts to value one row above
  bin_len_cum[[u]] <- bin_len_cum[[u-1]] + bin_len_list[[u]]
}

final_bin$Cummulative_Bin_Length <- unlist(bin_len_cum)
final_bin$Cummulative_Mutations <- unlist(cum_mut_list)

#######################################
## Apply Cutoff of BP Length if Needed
#######################################

final_bin <- final_bin[final_bin$Cummulative_Bin_Length <= panel_length,] ##apply bp length cutoff if needed

final_bin <- final_bin[,c(1:4,7:8)]

write.csv(final_bin, "E:\\F_FW_Bins_011122.csv")