#swap_values function
#Version: 1.0
#Last_update:2024/10/11
#Programmer: David_P_Quevedo
swap_values <- function(data) {
  library(dplyr)
  library(tidyr)
  library(stringr)
  
  bcd <- data %>%
    dplyr::select(-starts_with("bcf")) %>%
    rename_all(~ gsub("bcd.", "", .))
  
  bcd_gath <- bcd %>%
    gather("partner_colour", "bcd", red, violet, green, darkblue, orange, blue, yellow, pink) %>%
    filter(colour != partner_colour, bcd != "-") %>%
    mutate(sexselF = ifelse(sexsel == "1", "POL", "MON"),
           popsubF = ifelse(popsub == "1", "ST", "NST"))
  
  bcf <- data %>%
    dplyr::select(-starts_with("bcd")) %>%
    rename_all(~ gsub("bcf.", "", .))
  
  bcf_gath <- bcf %>%
    gather("partner_colour", "bcf", red, violet, green, darkblue, orange, blue, yellow, pink) %>%
    filter(colour != partner_colour, bcf != "-") %>%
    mutate(sexselF = ifelse(sexsel == "1", "POL", "MON"),
           popsubF = ifelse(popsub == "1", "ST", "NST"))
		   
  
  # Remove repeated contacts per video
  clean_repeated_contacts <- function(data, value_column) {
    data %>%
      group_by(video) %>%
      mutate(pair = if_else(colour < partner_colour,
                            paste(colour, partner_colour, sep = "-"),
                            paste(partner_colour, colour, sep = "-"))) %>%
      group_by(video, pair) %>%
      summarize(across(all_of(value_column), mean, na.rm = TRUE), .groups = "drop") %>%
      separate(pair, into = c("colour", "partner_colour"), sep = "-") %>%
      ungroup()
  }
  
  bcd_data <- clean_repeated_contacts(bcd_gath, "bcd")
  bcf_data <- clean_repeated_contacts(bcf_gath, "bcf")
  
  
  # Shuffle 'bcd' column while keeping the structure
  
  bcd_data <- bcd_data %>%
  mutate(bcd = sample(bcd))  
  
  bcf_data <- bcf_data %>%
  mutate(bcf = sample(bcf))  

  
  # Recreate the undirected network contacts per video
  recreate_undirected_contacts <- function(data) {
    reversed_data <- data %>%
      rename(colour = partner_colour, partner_colour = colour)
    
    bind_rows(data, reversed_data) %>%
      distinct()
  }
  
  bcd_data <- recreate_undirected_contacts(bcd_data)
  bcf_data <- recreate_undirected_contacts(bcf_data)
  
  # Spread and rename for 'bcd'
  spread_purged_data_bcd <- bcd_data %>%
    spread(key = partner_colour, value = bcd) %>%
    dplyr::select("colour", "video", red, violet, green, darkblue, orange, blue, yellow, pink) %>%
    mutate_all(~replace(., is.na(.), 0))
  
  target <- c("red", "violet", "green", "darkblue", "orange", "blue", "yellow", "pink")
  spread_purged_data_bcd <- spread_purged_data_bcd %>%
    arrange(factor(colour, levels = target))
  colnames(spread_purged_data_bcd)[3:10] <- paste0("bcd.", colnames(spread_purged_data_bcd)[3:10])
  
  # Spread and rename for 'bcf'
  spread_purged_data_bcf <- bcf_data %>%
    spread(key = partner_colour, value = bcf) %>%
    dplyr::select("colour", "video", red, violet, green, darkblue, orange, blue, yellow, pink) %>%
    mutate_all(~replace(., is.na(.), 0))
  
  spread_purged_data_bcf <- spread_purged_data_bcf %>%
    arrange(factor(colour, levels = target))
  colnames(spread_purged_data_bcf)[3:10] <- paste0("bcf.", colnames(spread_purged_data_bcf)[3:10])
  
  # Join both spreads
  spread_purged_data <- left_join(spread_purged_data_bcd, spread_purged_data_bcf, by = c("colour", "video"))
  spread_purged_data <- mutate_all(spread_purged_data, ~replace(., is.na(.), 0))
  
  # Join with recording data details
  dataset_details <- alldata[, c("video", "line", "treatment", "starthour", "endhour", "generation", 
                                 "daterecord", "avgage", "trial", "sexsel", "popsub", "line.unique", 
                                 "colour", "id", "sex", "interpolation", "sbnf", "totaldistmov.mm", 
                                 "velocity", "movingtime")]
  
  fulldata_purged_data <- merge(dataset_details, spread_purged_data, by = c("colour", "video"))
  fulldata_purged_data <- as.data.frame(fulldata_purged_data)
  fulldata_purged_data <- fulldata_purged_data %>% 
							arrange(video, line.unique)
  # Clean up
  rm(spread_purged_data_bcd, spread_purged_data_bcf, spread_purged_data, dataset_details, bcd_data, bcf_data)
  
  return(fulldata_purged_data)
}

