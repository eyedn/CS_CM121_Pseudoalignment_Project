#=============================================================================#
#       CS CM121
#       Winter 2023
#       Project 2
#       Aydin Karayas
#=============================================================================#
set.seed(108)
args <- commandArgs(trailingOnly = TRUE)
library(ggplot2)

#=============================================================================#
#       DEFINE GLOBAL VARIABLE(S)
#=============================================================================#
# define all reads data
all_reads <- read.csv(args[1], header = FALSE)
all_ec <- as.data.frame(table(all_reads$V3))
all_ec$Var1 <- as.character(all_ec$Var1)

# define correct reads data
cor_reads <- all_reads[all_reads[, 4] == "True", ]
cor_ec <- as.data.frame(table(cor_reads$V3))
cor_ec$Var1 <- as.character(cor_ec$Var1)

# define false reads data
err_reads <- all_reads[all_reads[, 4] == "False", ]
err_ec <- as.data.frame(table(err_reads$V3))
err_ec$Var1 <- as.character(err_ec$Var1)

#=============================================================================#
#       FUNCTIONS
#=============================================================================#
# format equiv counts table with no. occurrences, no. transcripts, equiv class
create_equiv_class_table <- function(equiv_class_counts) {
  # counts (no. occurrences) is the no. times that equiv class is mapped to
  counts <- as.numeric(equiv_class_counts[, 2])
  # no_tps (no. transcripts) is the size of the equiv class
  no_tps <- sapply(equiv_class_counts[, 1], count_transcripts)
  # equiv_class (equiv class) contains all transcripts for that equiv class
  equiv_class <- as.character(equiv_class_counts[, 1])
  # format the new equiv class occurances data frame
  equiv_class_occ <- data.frame(counts, no_tps, equiv_class, 
                                stringsAsFactors = FALSE, row.names = NULL)
  colnames(equiv_class_occ) <- c("counts", 
                          "number of items in equivalence class", 
                          "isoforms in equivalence class")
  equiv_class_occ[equiv_class_occ == ""] <- NA
  return(equiv_class_occ)
}

# given a comma separated list of transcripts, return the count of transcripts
count_transcripts <- function(equiv_class_str) {
  transcript_list <- strsplit(equiv_class_str, ",")[[1]]
  return(length(transcript_list))
}

# plot of number of items (noi) in an ec vs number of mapping to that noi
equiv_class_size_summary <- function(equiv_class_occ, subset) {
  # noi stores equiv class sizes vs. mapped reads or equiv classes of that size
  # this will be plotted by ggplot
  noi <- find_number_of_items(equiv_class_occ)
  ggplot(noi, aes(x = noi_index, y = noi_count)) +
    # generate graph visuals
    geom_segment( aes(x=noi_index ,xend=noi_index, y=0, yend=noi_count)) + 
    geom_point(size = 3, pch = 21, fill = "lightsalmon", color = "black") +
    # create specific label for noi = 0
    geom_point(data = noi[noi[,1] == 0,],
               size = 3, pch = 10, color = "black") +
    geom_text(aes(label=ifelse(noi_index==0,
                               paste("NULL SET",noi_count,sep = ": "),'')),
              hjust=0.1,vjust=1.8, size = 3, color = "black") +
    # create specific label for noi = 1
    geom_point(data = noi[noi[,1] == 1,], 
               size = 3, pch = 10, color = "black") +
    geom_text(aes(label=ifelse(noi_index==1,
                               paste(noi_index,noi_count,sep = ": "),'')),
              hjust=0.1,vjust=-1, size = 3, color = "black") +
    # create specific label for noi = 2
    geom_point(data = noi[noi[,1] == 2,],
               size = 3, pch = 10, color = "black") +
    geom_text(aes(label=ifelse(noi_index==2,
                               paste(noi_index,noi_count,sep = ": "),'')),
              hjust=0.1,vjust=-1, size = 3, color = "black") +
    # create specific label for noi = 3
    geom_point(data = noi[noi[,1] == 3,], 
               size = 3, pch = 10, color = "black") +
    geom_text(aes(label=ifelse(noi_index==3,
                               paste(noi_index,noi_count,sep = ": "),'')),
              hjust=0.1,vjust=-1, size = 3, color = "black") +
    # create specific label for noi = 10
    geom_point(data = noi[noi[,1] == 10,], 
               size = 3, pch = 10, color = "black") +
    geom_text(aes(label=ifelse(noi_index==10,
                               paste(noi_index,noi_count,sep = ": "),'')),
              hjust=0.2,vjust=-1.2, size = 3, color = "black") +
    # create specific label for noi = 20
    geom_point(data = noi[noi[,1] == 20,],
               size = 3, pch = 10, color = "black") +
    geom_text(aes(label=ifelse(noi_index==20,
                               paste(noi_index,noi_count,sep = ": "),'')),
              hjust=0.2,vjust=-1.2, size = 3, color = "black") +
    # create specific label for noi = 59
    geom_point(data = noi[noi[,1] == max(noi[,1]), ],
               size = 3, pch = 10, color = "black") +
    geom_text(aes(label=ifelse(noi_index==max(noi[,1]),
                               paste(noi_index,noi_count,sep = ": "),'')),
              hjust=0.2,vjust=-1.2, size = 3) +
    # adjust aesthetics
    ggtitle(paste("Relationship of Equivalance Class Size vs. Mapped Reads of", 
                   subset, "Data", sep = " ")) +
    xlab("Number of Items (transcripts) in a given Equivalence Class") +
    ylab("Number of Mapped Reads") +
    theme(plot.margin = margin(1,1,1,1, "cm"))
}

# print out summary statistics for equiv class size and occurrences
print_stats <- function(read_counts, equiv_class_occ, subset){
  cat(paste0(subset, " Reads: ", nrow(read_counts), "\n"))
  cat(paste0("Equivalence classes mapped to: ", nrow(equiv_class_occ), "\n"))
  cat("Equivalence class sizes:\n")
  print(summary(equiv_class_occ[,2]))
  cat("Equivalence class occurances:\n")
  print(summary(equiv_class_occ[,1]))
  cat("----------------------------------------------------------------\n")
}

# find the number of times a read map to equiv class of a specific size
find_number_of_items <- function(ec_occ_df) {
  # extract all possible equiv class sizes
  noi_index <- unique(ec_occ_df[,2])
  # find number of reads mapping to each equiv class size
  noi_count <- c()
  for (i in seq_len(length(noi_index))) {
    # extract rows of equiv class occurrences that have equiv class of size i
    to_be_summed <- as.matrix(ec_occ_df[ec_occ_df[,2] == noi_index[i], c(1,2)])
    # sum mapped reads of equiv class size i
    noi_count[i] <- sum(to_be_summed[,1])
  }
  # format noi data frame
  organized_noi <- data.frame(noi_index, noi_count, 
                              stringsAsFactors = FALSE, row.names = NULL)
  return(organized_noi)
}

#=============================================================================#
#       MAIN CODE
#=============================================================================#
# format equiv class occurrences in this table
all_ec_occ <- create_equiv_class_table(all_ec)
write.csv(all_ec_occ, file = "equiv_classes.csv", row.names = FALSE)
cor_ec_occ <- create_equiv_class_table(cor_ec)
err_ec_occ <- create_equiv_class_table(err_ec)

# graph noi in an ec vs. number of reads mapping to that noi in an ec
equiv_class_size_summary(all_ec_occ, subset = "All")
ggsave("all_noi.png", height = 6.14, width = 8.9, units = "in")
equiv_class_size_summary(cor_ec_occ, subset = "Correct")
ggsave("cor_noi.png", height = 6.14, width = 8.9, units = "in")
equiv_class_size_summary(err_ec_occ, subset = "Incorrect")
ggsave("err_noi.png", height = 6.14, width = 8.9, units = "in")

# print summary statistics
print_stats(all_reads, all_ec_occ, "All")
print_stats(cor_reads, cor_ec_occ, "Correct")
print_stats(err_reads, err_ec_occ, "Incorrect")