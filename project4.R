#

# read table, there IS a header, use tab as separator
data <- read.table("/work2/07475/vagheesh/stampede2/forOthers/forBIO321G/humanGeno.txt", header = TRUE, sep = "\t")

# generate all combinations of 8 strings in a vector taken 4 at a time
combinations <- combn(c("Chimp","Neanderthal", "Denisovan", "Dai", "British", "Yoruba", "Mbuti", "Papuan"), 4)
# transpose it to make a dataframe
combinations <- data.frame(t(combinations))

# get all of the "anatomically modern humans"
amh_table <- combinations[6:25,]
amh_table <- amh_table
colnames(amh_table) <- c("Pop1","Pop2", "Pop3", "Pop4")
# add up the remaining columns
amh_table <- cbind(amh_table, ABBA=NA, BABA=NA, JackMean=NA, JackError=NA, Z=NA)

# get the counts of ABBA and BABA given the row and the dataset
get_counts <- function(row, dataset) {
  # put the columns from the dataset in a dataframe
  temp_df <- data.frame(dataset[,row[,1]], dataset[,row[,2]], dataset[,row[,3]], dataset[,row[,4]])
  # turn the rows from the dataframe into strings for counting
  temp_paste <- paste(temp_df[,1], temp_df[,2], temp_df[,3], temp_df[,4])
  # make it a table so we can get the frequencies
  all_counts <- table(temp_paste)
  # add up the counts
  abba_counts <- all_counts[["0 2 2 0"]] + all_counts[["2 0 0 2"]]
  baba_counts <- all_counts[["2 0 2 0"]] + all_counts[["0 2 0 2"]]
  # return the counts
  return (c(abba_counts, baba_counts))
}

# fill in the table row by row
for (row in 1:nrow(amh_table)) {
  # get total ABBA and BABA counts and put them in the table
  totals <- get_counts(amh_table[row,], data)
  amh_table[row, 5:6] <- totals
  
  # store the D_stats here, we'll need them for both the mean and error
  D_stats <- c()
  
  # calculate the D_stats
  for (block in 1:5) {
    # quick mafs
    start <- ifelse(block == 1, 1, (block * 2 - 2) * 10000 + 1)
    end <- min(nrow(data), start + 19999)
    
    # this is the entire table excluding the current block
    curr_block <- data[-(start:end),]
    
    # get the counts and calculate the d_stat
    vals <- get_counts(amh_table[row,], curr_block)
    d_stat <- (vals[2] - vals[1]) / (vals[1] + vals[2])
    
    # append to the vector and move on
    D_stats <- c(D_stats, d_stat)
  }
  
  # calculate the JackMean
  first_leave_one_out <- c()
  for (stat in 1:5) {
    first_leave_one_out <- c(first_leave_one_out, mean(D_stats[-stat]))
  }
  jack_mean <- mean(first_leave_one_out)
  amh_table[row, 7] <- jack_mean
  
  # the JackError is then simple
  jack_error <- sqrt(sum((first_leave_one_out - jack_mean)^2) / 4)
  amh_table[row, 8] <- jack_error
  
  # the Z score is the mean divided by the error
  Z <- jack_mean / jack_error
  amh_table[row, 9] <- Z
}
write.table(amh_table, "/work2/09269/saip/stampede2/projects/project4/humanSeqD.txt", sep="\t", col.names=TRUE,row.names=FALSE)


