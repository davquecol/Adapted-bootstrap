#Net_data_calculation function
#Version: 2.1.1
#Last_update:2024/10/09
#Programmer: David_P_Quevedo

net_data_calculation <- function(x) {       
  

  #GETTING THE PROPER STRUCTURE OF THE DATA FOR BCD
  
  bcd <- x %>%
    dplyr::select(!(names(x)[ substr(names(x), 13, 15) %in% "bcd"]))  
  
  names(bcd) <- gsub("bcd.", "", names(bcd)) 
  
  # head(bcd)
  
  bcd_gath <- bcd %>% 
    gather("partner_colour", "bcd", red, violet, green, darkblue, orange, blue, yellow, pink) %>%
    filter(colour != partner_colour) %>% 
    filter(!(bcd %in% "-"))
  
  # print(unique(bcd_gath$colour)  %in% unique(bcd_gath$partner_colour))
  
  bcd_gath<-bcd_gath%>% dplyr::select (-c(bcf.red, bcf.violet, bcf.green, bcf.darkblue,bcf.orange, bcf.blue,bcf.yellow,bcf.pink))
  
  
  #GETTING THE PROPER STRUCTURE OF THE DATA FOR BCF
  
  
  bcf <- x %>%
    dplyr::select(!(names(x)[substr(names(x), 13, 15) %in% "bcf"]))  
  
  names(bcf) <- gsub("bcf.", "", names(bcf)) 
  
  # head(bcf)
  
  bcf_gath <- bcf %>% 
    gather("partner_colour", "bcf", red, violet, green, darkblue, orange, blue, yellow, pink) %>%
    filter(colour != partner_colour) %>% 
    filter(!(bcf %in% "-"))
  
  # print(unique(bcf_gath$colour)  %in% unique(bcf_gath$partner_colour))
  
  bcf_gath<-bcf_gath%>% dplyr::select (-c(bcd.red, bcd.violet, bcd.green, bcd.darkblue,bcd.orange, bcd.blue,bcd.yellow,bcd.pink))
  
  
  full_gath <- full_join(bcd_gath,bcf_gath, by=c("video","line","treatment","starthour","endhour","generation", 
                                                 "daterecord","avgage","trial","sexsel","popsub","line.unique", 
                                                 "colour","id", "sex", "interpolation","sbnf","totaldistmov.mm", 
                                                 "velocity", "movingtime","partner_colour"))
  threshold_data<-full_gath
  
  # # BUILDING NETWORKS 
  
  # bothways <- threshold_data
  # oneways <- bothways[!duplicated(t(apply(bothways[c("colour", "partner_colour")], 1, sort))),]

  
  # GETTING THE STRUCTURE FOR TNET ------- BCD
  
  from.bcd <-threshold_data$colour
  to.bcd <-  threshold_data$partner_colour
  weight.bcd <- as.numeric(threshold_data$bcd)
  tnatdata.bcd <- cbind.data.frame(from.bcd, to.bcd, weight.bcd) 
    
  # unique(from.bcd)
  tnatdata.bcd <- replace(tnatdata.bcd,tnatdata.bcd == "red", "1") 
  tnatdata.bcd <- replace(tnatdata.bcd,tnatdata.bcd == "green", "2")
  tnatdata.bcd <- replace(tnatdata.bcd,tnatdata.bcd == "violet", "3")
  tnatdata.bcd <- replace(tnatdata.bcd,tnatdata.bcd == "darkblue", "4")
  tnatdata.bcd <- replace(tnatdata.bcd,tnatdata.bcd == "orange", "5")
  tnatdata.bcd <- replace(tnatdata.bcd,tnatdata.bcd == "blue", "6")
  tnatdata.bcd <- replace(tnatdata.bcd,tnatdata.bcd == "yellow", "7")
  tnatdata.bcd <- replace(tnatdata.bcd,tnatdata.bcd == "pink", "8")
  
  colnames(tnatdata.bcd)[1] <- "from.bcd"
  colnames(tnatdata.bcd)[2] <- "to.bcd"
  colnames(tnatdata.bcd)[3] <- "weight.bcd"
  
  
  tnatnet.bcd <- as.tnet(tnatdata.bcd, type="weighted one-mode tnet")
  # str(tnatnet.bcd)
  tnatnet.bcd$i <- as.integer(tnatnet.bcd$i)
  tnatnet.bcd$j <- as.integer(tnatnet.bcd$j)
  
  
  # GETTING THE STRUCTURE FOR TNET ------- BCF
  
  from.bcf <-threshold_data$colour
  to.bcf <-  threshold_data$partner_colour
  weight.bcf <- as.numeric(threshold_data$bcf)
  
  tnatdata.bcf <- cbind.data.frame(from.bcf, to.bcf, weight.bcf) 
  
  # unique(from.bcf)
  tnatdata.bcf <- replace(tnatdata.bcf,tnatdata.bcf == "red", "1") 
  tnatdata.bcf <- replace(tnatdata.bcf,tnatdata.bcf == "green", "2")
  tnatdata.bcf <- replace(tnatdata.bcf,tnatdata.bcf == "violet", "3")
  tnatdata.bcf <- replace(tnatdata.bcf,tnatdata.bcf == "darkblue", "4")
  tnatdata.bcf <- replace(tnatdata.bcf,tnatdata.bcf == "orange", "5")
  tnatdata.bcf <- replace(tnatdata.bcf,tnatdata.bcf == "blue", "6")
  tnatdata.bcf <- replace(tnatdata.bcf,tnatdata.bcf == "yellow", "7")
  tnatdata.bcf <- replace(tnatdata.bcf,tnatdata.bcf == "pink", "8")
  
  colnames(tnatdata.bcf)[1] <- "from.bcf"
  colnames(tnatdata.bcf)[2] <- "to.bcf"
  colnames(tnatdata.bcf)[3] <- "weight.bcf"
  
  
  tnatnet.bcf <- as.tnet(tnatdata.bcf, type="weighted one-mode tnet")
  # str(tnatnet.bcf)
  tnatnet.bcf$i <- as.integer(tnatnet.bcf$i)
  tnatnet.bcf$j <- as.integer(tnatnet.bcf$j)
  
  
  
  # GETTING CENTRALITY MEASURES FROM BCD
  
  
  clustering_bcd <- as.data.frame(clustering_w(tnatnet.bcd, measure=c("bi", "gm")))
  clustering_bcf <- as.data.frame(clustering_w(tnatnet.bcf, measure=c("bi", "gm")))
  
  betweenness_bcd <-  as.data.frame(betweenness_w(tnatnet.bcd, directed=FALSE, alpha=1))
  betweenness_bcf <-  as.data.frame(betweenness_w(tnatnet.bcf, directed=FALSE, alpha=1))
  
  degree_bcd <-  as.data.frame(degree_w(tnatnet.bcd,measure=c("degree","output"), type="all", alpha=1))
  # degree.bcd <- subset(degree.bcd, degree.bcd$output>0)
  
  degree_bcf <-  as.data.frame(degree_w(tnatnet.bcf,measure=c("degree","output"), type="all", alpha=1))
  # degree.bcf <- subset(degree.bcf, degree.bcf$output>0)
  
  node<-unique(from.bcd)
  node<- 1:length(node)
  
  num_edges_bcd <- nrow(tnatnet.bcd)
  num_nodes_bcd <- length(unique(c(tnatnet.bcd$i, tnatnet.bcd$j)))
  possible_edges_bcd <- num_edges_bcd * (num_edges_bcd - 1)
  network_density_bcd <- num_edges_bcd / possible_edges_bcd
  density_bcd<- rep(network_density_bcd, max(node))
  density_bcd<- cbind.data.frame(node, density_bcd)
  
  num_edges_bcf <- nrow(tnatnet.bcf)
  num_nodes_bcf <- length(unique(c(tnatnet.bcf$i, tnatnet.bcf$j)))
  possible_edges_bcf <- num_edges_bcf * (num_edges_bcf - 1)
  network_density_bcf <- num_edges_bcf / possible_edges_bcf
  density_bcf<- rep(network_density_bcf,max(node))
  density_bcf<- cbind.data.frame(node, density_bcf)

  colnames(clustering_bcd)[1]<- "values"
  clustering_bcd <- tibble::rownames_to_column(clustering_bcd, "rn") 
  clustering_bcd <- spread(clustering_bcd, key = rn, value = values)
  clustering_bcd <- cbind.data.frame(node,clustering_bcd)
  
  colnames(clustering_bcf)[1]<- "values"
  clustering_bcf <- tibble::rownames_to_column(clustering_bcf, "rn") 
  clustering_bcf <- spread(clustering_bcf, key = rn, value = values)
  clustering_bcf <- cbind.data.frame(node,clustering_bcf)
  
  
  gm_divby_bi_bcd <- as.data.frame(clustering_bcd$gm/clustering_bcd$bi) 
  gm_divby_bi_bcd <- cbind.data.frame(node, gm_divby_bi_bcd)
  colnames(gm_divby_bi_bcd)[2]<- "gm_divby_bi_bcd"
  
  gm_divby_bi_bcf <- as.data.frame(clustering_bcf$gm/clustering_bcf$bi) 
  gm_divby_bi_bcf <- cbind.data.frame(node, gm_divby_bi_bcf)
  colnames(gm_divby_bi_bcf)[2]<- "gm_divby_bi_bcf"
  
  ####part for diversity####
	  	  
  transform_tnatdata_bcf<-as.data.frame(tnatnet.bcf)
  colnames(transform_tnatdata_bcf)[1]<-"from"
  colnames(transform_tnatdata_bcf)[2]<-"to"
  colnames(transform_tnatdata_bcf)[3]<-"weight"

  transform_tnatdata_bcd<-tnatnet.bcd
  colnames(transform_tnatdata_bcd)[1]<-"from"
  colnames(transform_tnatdata_bcd)[2]<-"to"
  colnames(transform_tnatdata_bcd)[3]<-"weight"

  # Function to calculate entropy for a node in an undirected graph
  calculate_entropy <- function(data, node) {
	  # Filter edges connected to the node
	  edges <- which(data$from == node | data$to == node)
	  
	  if (length(edges) == 0) return(0)
	  
	  # Retrieve weights of these edges
	  weights <- data$weight[edges]
	  
	  # Compute total weight
	  total_weight <- sum(weights)
	  
	  # Compute entropy
	  proportions <- weights / total_weight
	  entropy <- -sum(proportions * log(proportions))
	  
	  return(entropy)
	}

  # Calculate entropy for each node
  nodes.bcf <- unique(c(transform_tnatdata_bcf$from, transform_tnatdata_bcf$to))
  nodes.bcd<- unique(c(transform_tnatdata_bcd$from, transform_tnatdata_bcd$to))

  entropy_values_bcf <- sapply(nodes.bcf, function(node) calculate_entropy(transform_tnatdata_bcf, node))
  entropy_values_bcd <- sapply(nodes.bcd, function(node) calculate_entropy(transform_tnatdata_bcd, node))

  # Create a data frame to store the results
  entropy_bcf <-cbind.data.frame(
  node = nodes.bcf,
  entropy_bcf = entropy_values_bcf
  )

  entropy_bcd<- cbind.data.frame(
  node = nodes.bcd,
  entropy_bcd = entropy_values_bcd
  )

  diversity<-full_join(entropy_bcd,entropy_bcf, by="node")
  diversity$node<-as.character(diversity$node)
 
  measures1 <- full_join(clustering_bcd,betweenness_bcd, by = "node")
  measures2 <- full_join(degree_bcd, density_bcd, by = "node")
  
  measures11 <- full_join(clustering_bcf,betweenness_bcf, by = "node")
  measures22 <- full_join(degree_bcf, density_bcf, by = "node")
  
  allmeasures.bcd <-  full_join(measures1,measures2, by = "node")
  allmeasures.bcd <-  full_join(allmeasures.bcd,gm_divby_bi_bcd, by = "node")
  colnames(allmeasures.bcd)[2]<-"bi_bcd"
  colnames(allmeasures.bcd)[3]<-"gm_bcd"
  colnames(allmeasures.bcd)[4]<-"betweenness_bcd"
  colnames(allmeasures.bcd)[5]<-"degree_bcd"
  colnames(allmeasures.bcd)[6]<-"strength_bcd"
  
  allmeasures.bcf <-  full_join(measures11,measures22, by = "node")
  allmeasures.bcf <-  full_join(allmeasures.bcf,gm_divby_bi_bcf, by = "node")
  colnames(allmeasures.bcf)[2]<-"bi_bcf"
  colnames(allmeasures.bcf)[3]<-"gm_bcf"
  colnames(allmeasures.bcf)[4]<-"betweenness_bcf"
  colnames(allmeasures.bcf)[5]<-"degree_bcf"
  colnames(allmeasures.bcf)[6]<-"strength_bcf"
  
  allmeasures <-  full_join(allmeasures.bcd,allmeasures.bcf,by = "node")
  allmeasures$node<-as.character(allmeasures$node)
  allmeasures <-  full_join(allmeasures, diversity, by="node")

  allmeasures$node[allmeasures$node == "1"] <- "red"
  allmeasures$node[allmeasures$node == "2"] <- "green"
  allmeasures$node[allmeasures$node == "3"] <- "violet"
  allmeasures$node[allmeasures$node == "4"] <- "darkblue"        
  allmeasures$node[allmeasures$node == "5"] <- "orange"
  allmeasures$node[allmeasures$node == "6"] <- "blue"   
  allmeasures$node[allmeasures$node == "7"] <- "yellow"
  allmeasures$node[allmeasures$node == "8"] <- "pink"
                  
  video <- unique(x$video)
  video <- rep(video,length(allmeasures$node))
  allmeasures <- cbind.data.frame(video, allmeasures)
  colnames(allmeasures)[2]<-"colour"
        
  x$sexselF <-ifelse(x$sexsel=="1", "POL", "MON")
  x$popsubF <-ifelse(x$popsub=="1", "ST","NST")
  x$sexselF <-as.factor(x$sexselF)
  x$popsubF <-as.factor(x$popsubF)
  x<-x%>% dplyr::select (-c(bcf.red, bcf.violet, bcf.green, bcf.darkblue,bcf.orange, bcf.blue,bcf.yellow,bcf.pink,
  bcd.red, bcd.violet, bcd.green, bcd.darkblue,bcd.orange, bcd.blue,bcd.yellow,bcd.pink))
                
  
 ## Check the dimensions of the data frames before joining
  # cat("Rows in allmeasures for video", unique(x$video), ":", nrow(allmeasures), "\n")
  # cat("Rows in x for video", unique(x$video), ":", nrow(x), "\n")  
  original_and_centrality <-  full_join(allmeasures,x, by=c("video","colour"))
  
}

