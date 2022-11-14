set.seed(1)
options(tz = "Asia/Jerusalem")
Sys.setenv(TZ=getOption("tz"))
Sys.getenv("TZ")
library(lubridate)
library(sp)
library(openxlsx)
library(readr)
library(questionr)
library(dplyr)
library(plyr)
library(psych)
library(gdata)
library(DescTools)
library(data.table)
library(rlist)
library(exactRankTests)
library(asnipe)
library(sna)
library(igraph)
library(tidygraph)
library(devtools)
library(transformr)
library(ggplot2)
library(ggthemes)
library(ggridges)
library(ggpubr)
library(gridExtra)
library(RColorBrewer)
library(wesanderson)
theme_set(theme_base())



#############################################################################
######################## DATA CALL & CLEANING ###############################
#############################################################################
# Experiment dates
exp_start   <- as.POSIXct("14-07-2017 06:00:00", format = "%d-%m-%Y %H:%M:%S", tz = 'Asia/Jerusalem')
exp_end     <- as.POSIXct("11-08-2017 06:00:00", format = "%d-%m-%Y %H:%M:%S", tz = 'Asia/Jerusalem')


# Build the list of individuals
LoggersBrut <- read.csv(file = paste0("~/01_DeployedLoggers.csv"), header = TRUE, 
                        sep = ";", na = c("","NA"))
Loggers          <- LoggersBrut[!is.na(LoggersBrut$Date_on),]
Loggers$Date_on  <- as.Date(strptime(Loggers$Date_on,  "%d/%m/%Y %H:%M", tz = 'Asia/Jerusalem'))
Loggers$Date_off <- as.Date(strptime(Loggers$Date_off, "%d/%m/%Y %H:%M", tz = 'Asia/Jerusalem'))
Trapped          <- Loggers[!is.na(Loggers$Date_off),]
List             <- Loggers[is.na(Loggers$Date_off),]
List$Date_off    <- List$Date_on + months(8)
list_ind         <- rbind(Trapped[data.table::between(exp_start, Trapped$Date_on, Trapped$Date_off),], 
                          List[data.table::between(exp_start, List$Date_on,    List$Date_off),])
list_ind         <- list_ind[!is.na(list_ind$Chip),]


# Metadata
sunlight    <- read.csv("~/02_sunlightTimes.csv", header = T)
dyads_list  <- read.csv("~/03_Dyads.csv", header = T, sep = ",")
dyads_list  <- dyads_list[dyads_list$deployment==year(exp_start),]
dyads_list  <- dyads_list[dyads_list$retreived=="2_logger",]


# List of proximity events
Contactlist <- read.csv("~/04_ContactList.csv", header = T, sep = ",")
Contactlist$interval <- as.POSIXct(strptime(as.character(Contactlist$interval), format = "%Y-%m-%d %H:%M:%S", tz = 'Asia/Jerusalem'))
Contactlist$start <- as.POSIXct(strptime(as.character(Contactlist$start), format = "%Y-%m-%d %H:%M:%S", tz = 'Asia/Jerusalem'))
Contactlist$date  <- as.POSIXct(strptime(as.character(Contactlist$date), format = "%Y-%m-%d", tz = 'Asia/Jerusalem'))


# Adjust for missing collars
# Randomly removes one collar per day (for dyad where both collars were retreived)
li_f              <- c()
doubles <- Contactlist %>% group_by(date, dyad) %>% 
  dplyr::summarise(n = length(unique(individual1))) %>%
  mutate(individual1 = sapply(strsplit(dyad, ":"), "[", 1),
         individual2 = sapply(strsplit(dyad, ":"), "[", 2),
         to_remove   = sample(c(4,5), size=n(), replace=TRUE)) %>%
  subset(n>1)
doubles$ind_to_remove <- as.data.frame(doubles)[cbind(seq_len(nrow(doubles)), doubles$to_remove)]
li_f <- c()
for (i in 1:nrow(doubles)) {
  dy <- doubles$dyad[i]
  da <- doubles$date[i]
  id <- doubles$ind_to_remove[i]
  j <- which(Contactlist$dyad==dy & Contactlist$date==da & Contactlist$individual1==id)
  li_f <- c(li_f, j)
}
Contactlist <- Contactlist[-li_f,]
Contactlist$chip     <- Contactlist$dyad
Contactlist$dyad     <- strsplit(as.character(Contactlist$dyad), ":")



# Final subset
TimeEdgelist         <- Contactlist[between(Contactlist$interval, exp_start, exp_end),]
TimeEdgelist         <- TimeEdgelist[TimeEdgelist$individual1%in%list_ind$Chip & TimeEdgelist$individual2%in%list_ind$Chip,]
TimeEdgelist$Mome    <- (TimeEdgelist$date + dhours(TimeEdgelist$hour)) %>% as.character() %>% as.factor()
TimeEdgelist$site    <- list_ind$Canyon[match(TimeEdgelist$individual1, list_ind$Chip)]
TimeEdgelist$phase   <- "other"
TimeEdgelist$phase[TimeEdgelist$day_phase %in% c("morning", "afternoon")] <- "day"
TimeEdgelist$phase[TimeEdgelist$day_phase %in% c("early_night", "late_night")] <- "night"
TimeEdgelist         <- TimeEdgelist[TimeEdgelist$duration > 10,] 
for (ip in 1:nrow(TimeEdgelist)) {
  TimeEdgelist$n5[ip] <- sunlight$n5[TimeEdgelist$phase[ip]==sunlight$phase2 & TimeEdgelist$date[ip]==sunlight$date]
}




#############################################################################
###################### FULL NETWORK & COMMUNITY #############################
#############################################################################
# Full network
gbi_full <- get_group_by_individual(TimeEdgelist$dyad, data_format = "groups")
gbi_full <- gbi_full[,order(colnames(gbi_full))]
invisible(
  capture.output(
    net_full <- get_network(association_data  =  gbi_full, 
                            data_format       = "GBI", 
                            association_index = 'SRI', 
                            identities        = colnames(gbi_full))
  )
)
gra_full <- graph_from_adjacency_matrix(net_full, "undirected", weighted = T)



# Overlapping community detection + network plot
Community    <- linkcomm::getOCG.clusters(get.edgelist(gra_full))
adds         <- setdiff(names(V(gra_full)), Community$nodeclusters$node)
max          <- as.numeric(max(Community$nodeclusters$cluster))
if (length(adds > 0)) {
  Community$nodeclusters    <- rbind(Community$nodeclusters,
                                     data.frame(node = adds,
                                                cluster = max+c(1:length(adds))))
}
Community$nodeclusters$ID <- 1:nrow(Community$nodeclusters)  # unique variable
lst <- Community$nodeclusters %>% tidyr::spread(cluster, node) %>% select(-ID) %>% as.list()
lapply(lst, function(x) x[!is.na(x)])
linkcomm::plotOCGraph(Community, pal=brewer.pal(10,"Paired"), vertex.radius=0.05, vlabel=FALSE)



#############################################################################
########################## DATA EXPLORATION #################################
#############################################################################
### Length of proximity contacts
TimeEdgelist$g1 <- "between"
for (dd in unique(TimeEdgelist$chip)) {
  
  inds <- strsplit(dd, split = ":")
  ind1 <- inds[[1]][1]
  ind2 <- inds[[1]][2]
  p1   <- 1 * sapply(lst, `%in%`, x = ind1)
  p2   <- 1 * sapply(lst, `%in%`, x = ind2)
  if (any(p1&p2)) {TimeEdgelist$g1[which(TimeEdgelist$chip==dd)] <- "within"}
  
}
TimeEdgelist$category <- cut(TimeEdgelist$duration, 
                             breaks = c(10, 25*60, 50*60, 75*60, 100*60, 125*60, 150*60, Inf), 
                             labels = c("<25", "25-50", "50-75", "75-100", "100-125", "125-150", ">150"))
table(TimeEdgelist$category)/nrow(TimeEdgelist)
table(TimeEdgelist$category[TimeEdgelist$phase=="night"])/nrow(TimeEdgelist[TimeEdgelist$phase=="night",])
table(TimeEdgelist$category[TimeEdgelist$phase=="day"])/nrow(TimeEdgelist[TimeEdgelist$phase=="day",])


### Temporal distribution of social contacts
# Day-night differences in number and length of social contacts
ggggg <- c()
TimeEdgelist <- TimeEdgelist[!duplicated(TimeEdgelist[,c("chip", "start", "duration")]),]
for (i in unique(as.character(TimeEdgelist$date))) {
  i  <- as.POSIXct(i, tz = "Asia/Jerusalem")
  d  <- which(TimeEdgelist$date==i & TimeEdgelist$phase=="day")
  n1 <- intersect(
    intersect(
      which(TimeEdgelist$date==i), 
      which(TimeEdgelist$phase=="night")),
    which(hour(TimeEdgelist$interval)>12)
  )
  n2 <- intersect(
    intersect(
      which(TimeEdgelist$date==(i+ddays(1))), 
      which(TimeEdgelist$phase=="night")),
    which(hour(TimeEdgelist$interval)<12)
  )
  mu_d <- mean(TimeEdgelist$duration[d])
  mu_n <- mean(TimeEdgelist$duration[c(n1,n2)])
  nb_d <- length(TimeEdgelist$duration[d])
  nb_n <- length(TimeEdgelist$duration[c(n1,n2)])
  ggggg <- rbind.data.frame(ggggg,
                            data.frame(date = i,
                                       day_duration = mu_d,
                                       nig_duration = mu_n,
                                       day_interact = nb_d,
                                       nig_interact = nb_n))
}
ggggg1 <- ggggg[complete.cases(ggggg),]
acf(ggggg1$nig_interact)
acf(ggggg1$nig_duration)
acf(ggggg1$day_interact)
acf(ggggg1$day_duration)

# observed test statistics
# Wilcoxon test on contact length
nnnn <- c(ggggg1$nig_duration)
dddd <- c(ggggg1$day_duration)
staow <- wilcox.test(nnnn, dddd, paired = TRUE)
mumuow <- mean(nnnn - dddd)
# Student t-test on number of contacts
nnnn <- c(ggggg1$nig_interact)
dddd <- c(ggggg1$day_interact)
staot <- t.test(nnnn, dddd, paired = TRUE)
mumuot <- mean(nnnn - dddd)

# permutation test for paired samples
perm.test(abs(ggggg$nig_duration-ggggg$day_duration), alternative="greater", exact=TRUE)
perm.test(abs(ggggg$nig_interact-ggggg$day_interact), alternative="greater", exact=TRUE)





