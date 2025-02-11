#start_progress function
#Version: 3.0.0
#Last_update:2024/10/09
#Programmer: David_P_Quevedo

start_progress <- function() {
  
  set_working_directory <<- function(Folder, Subfolder, Group, Type) {
    home_dir <- path.expand("~")
    home_dir <- gsub("\\\\", "/", home_dir)
    parts <- strsplit(home_dir, "/")[[1]]
    Local_disk <- paste0(parts[1], "/")
    Users <- parts[2]
    Users_value <- parts[3]
    sep <- "/"
    Environment <- paste0(Local_disk, sep, Users, sep, Users_value, sep, Folder, sep, 
                          Subfolder, sep, Group, sep, Type)
    setwd(Environment)
  }
  
  set_working_directory("Dropbox", "david_quevedo_SexNet", "thesis", "text")
  
  netdata <- read.table("manual_correction96.txt", header = TRUE)
  recordingdata <- read.table("recorddataALLCG.txt", header = TRUE)
  alldata <- merge.data.frame(recordingdata, netdata, by = "video")
  alldata <- subset(alldata, !(video %in% c(2, 15, 17, 21, 23, 24, 27, 31, 57, 58)))
  
  bcd <- alldata %>%
    dplyr::select(-starts_with("bcf")) %>%
    rename_all(~ gsub("bcd.", "", .))
  
  bcd_gath <- bcd %>%
    gather("partner_colour", "bcd", red, violet, green, darkblue, orange, blue, yellow, pink) %>%
    filter(colour != partner_colour, bcd != "-") %>%
    mutate(sexselF = ifelse(sexsel == "1", "POL", "MON"),
           popsubF = ifelse(popsub == "1", "ST", "NST"))
  
  bcf <- alldata %>%
    dplyr::select(-starts_with("bcd")) %>%
    rename_all(~ gsub("bcf.", "", .))
  
  bcf_gath <- bcf %>%
    gather("partner_colour", "bcf", red, violet, green, darkblue, orange, blue, yellow, pink) %>%
    filter(colour != partner_colour, bcf != "-") %>%
    mutate(sexselF = ifelse(sexsel == "1", "POL", "MON"),
           popsubF = ifelse(popsub == "1", "ST", "NST"))
  
  # Function to compare pairs for bcd and bcf
  compare_pairs <- function(data, value_column) {
    same_count <- 0
    different_count <- 0
    
    videos <- unique(data$video)
    
    for (video in videos) {
      video_data <- subset(data, video == !!video)
      colours <- unique(video_data$colour)
      
      for (i in 1:length(colours)) {
        for (j in (i + 1):length(colours)) {
          col1 <- colours[i]
          col2 <- colours[j]
          
          value1 <- video_data %>% filter(colour == col1 & partner_colour == col2) %>% pull(!!sym(value_column))
          value2 <- video_data %>% filter(colour == col2 & partner_colour == col1) %>% pull(!!sym(value_column))
          
          if (length(value1) > 0 & length(value2) > 0) {
            if (value1 == value2) {
              same_count <- same_count + 1
            } else {
              different_count <- different_count + 1
            }
          }
        }
      }
    }
    
    return(list(same = same_count, different = different_count))
  }
  
  # Pre-matrix comparison
  bcd_comparison_pre <- compare_pairs(bcd_gath, "bcd")
  bcf_comparison_pre <- compare_pairs(bcf_gath, "bcf")
  
  cat("Pre-matrix comparison\n")
  cat("Number of pairs with the same data for bcd: ", bcd_comparison_pre$same, "\n")
  cat("Number of pairs with different data for bcd: ", bcd_comparison_pre$different, "\n")
  cat("Number of pairs with the same data for bcf: ", bcf_comparison_pre$same, "\n")
  cat("Number of pairs with different data for bcf: ", bcf_comparison_pre$different, "\n")
  
  # Filter and process data
  
  common_columns <- intersect(names(bcd_gath), names(bcf_gath))
  both_weights <- left_join(bcd_gath, bcf_gath, by = common_columns)
  
  both_weights$bcd[both_weights$bcd< 11] <- 0
  both_weights$bcf[both_weights$bcd == 0] <- 0
  both_gaths_filtered<-both_weights
  # condition <- both_weights$bcd >= 11
  # both_gaths_filtered <- both_weights[condition, ]
  both_gaths_filtered <- as.data.frame(both_gaths_filtered)
  bcd_gath_2 <- both_gaths_filtered %>% dplyr::select(-bcf)
  bcf_gath_2 <- both_gaths_filtered %>% dplyr::select(-bcd)
  bcd_data <- as.data.frame(bcd_gath_2)
  bcf_data <- as.data.frame(bcf_gath_2)
  
  min_bcd_value <- min(both_gaths_filtered$bcd)
  
  spread_purged_data_bcd <- spread(bcd_data, key = partner_colour, value = bcd) %>%
    dplyr::select("colour", "video", red, violet, green, darkblue, orange, blue, yellow, pink) %>%
    mutate_all(~replace(., is.na(.), 0))
  
  target <- c("red", "violet", "green", "darkblue", "orange", "blue", "yellow", "pink")
  spread_purged_data_bcd <- spread_purged_data_bcd %>% arrange(factor(colour, levels = target))
  colnames(spread_purged_data_bcd)[3:10] <- paste0("bcd.", colnames(spread_purged_data_bcd)[3:10])
  
  spread_purged_data_bcf <- spread(bcf_data, key = partner_colour, value = bcf) %>%
    dplyr::select("colour", "video", red, violet, green, darkblue, orange, blue, yellow, pink) %>%
    mutate_all(~replace(., is.na(.), 0))
  
  spread_purged_data_bcf <- spread_purged_data_bcf %>% arrange(factor(colour, levels = target))
  colnames(spread_purged_data_bcf)[3:10] <- paste0("bcf.", colnames(spread_purged_data_bcf)[3:10])
  
  spread_purged_data <- left_join(spread_purged_data_bcd, spread_purged_data_bcf, by = c("colour", "video"))
  spread_purged_data <- mutate_all(spread_purged_data, ~replace(., is.na(.), 0))
  
  dataset_details <- alldata[, c("video", "line", "treatment", "starthour", "endhour", "generation", 
                                 "daterecord", "avgage", "trial", "sexsel", "popsub", "line.unique", 
                                 "colour", "id", "sex", "interpolation", "sbnf", "totaldistmov.mm", 
                                 "velocity", "movingtime")]
  
  fulldata_purged_data <- merge(dataset_details, spread_purged_data, by = c("colour", "video"))
  fulldata_purged_data <- as.data.frame(fulldata_purged_data)
  
  # Post-matrix comparison
  bcd_post_matrix <- fulldata_purged_data %>%
    dplyr::select(video, colour, starts_with("bcd.")) %>%
    gather("partner_colour", "bcd", starts_with("bcd.")) %>%
    mutate(partner_colour = gsub("bcd.", "", partner_colour)) %>%
    filter(colour != partner_colour)
  
  bcf_post_matrix <- fulldata_purged_data %>%
    dplyr::select(video, colour, starts_with("bcf.")) %>%
    gather("partner_colour", "bcf", starts_with("bcf.")) %>%
    mutate(partner_colour = gsub("bcf.", "", partner_colour)) %>%
    filter(colour != partner_colour)
  
  bcd_comparison_post <- compare_pairs(bcd_post_matrix, "bcd")
  bcf_comparison_post <- compare_pairs(bcf_post_matrix, "bcf")
  
  cat("Post-matrix comparison\n")
  cat("Number of pairs with the same data for bcd: ", bcd_comparison_post$same, "\n")
  cat("Number of pairs with different data for bcd: ", bcd_comparison_post$different, "\n")
  cat("Number of pairs with the same data for bcf: ", bcf_comparison_post$same, "\n")
  cat("Number of pairs with different data for bcf: ", bcf_comparison_post$different, "\n")
  
  return(list(min_bcd_value = min_bcd_value,
              bcd_data = bcd_data,
              bcf_data = bcf_data,
              both_gaths_filtered = both_gaths_filtered,
              data_to_calculate_net = fulldata_purged_data,
              alldata = alldata))
}

