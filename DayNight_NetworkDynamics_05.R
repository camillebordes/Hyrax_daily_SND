################################################################################
###################### SUPPORTING INFORMATION PLOTS ############################
################################################################################
# Neighbors stability plots
data.day.night    <- openxlsx::loadWorkbook("~/08_networkTraits_dayNight_active.xlsx")
dat <- data.provi <- openxlsx::readWorkbook(data.day.night, sheet = 1, detectDates = TRUE)

plot.data <- c()
for (index in 1:501) {
  data.provi <- openxlsx::readWorkbook(data.day.night, 
                                       sheet = index,
                                       detectDates = TRUE)
  data.provi <- data.provi[order(data.provi$date, data.provi$phase, data.provi$name),]
  
  plot.data <- rbind.data.frame(
    plot.data,
    cbind.data.frame(
      data.provi[,c("assoc_stab_w", "assoc_stab_b", "name", "date", "phase")],
      data.frame(type = ifelse(index==1, "obs", "rand"))
    )
  )
}
plot.data$date     <- as.POSIXct(strptime(plot.data$date, format = "%Y-%m-%d"))
plot.data$phase    <- factor(plot.data$phase, c("day", "night"))
plot.data$hours    <- ifelse(plot.data$phase == "day", 11, 23)
plot.data$time     <- as.POSIXct(plot.data$date + hours(plot.data$hours), origin = origin)

a <- ggplot(plot.data,
            aes(x = time,
                y = assoc_stab_w,
                group = interaction(time, type),
                fill = type,
                color = type)) +
  geom_boxplot() +
  ylab("Individual selectivity") +
  theme_pubr() +
  rremove("legend") +
  theme(axis.title.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x  = element_blank(),
        axis.text.x  = element_text(size = 15))
b <- ggplot(plot.data,
            aes(x = time,
                y = assoc_stab_b,
                group = interaction(time, type),
                fill = type,
                color = type)) +
  scale_colour_manual(labels = c("observed network traits", "predicted under the null hypothesis"),
                      values = c("#E69F00", "#56B4E9")) +
  scale_fill_manual(labels = c("observed network traits", "predicted under the null hypothesis"),
                    values = c("#E69F00", "#56B4E9")) + 
  geom_boxplot() +
  ylab("Neighbors' stability") +
  theme_pubr() +
  theme(legend.title = element_blank())


library("cowplot")
plot_grid(a, b, NULL, NULL,
          ncol = 1,
          nrow = 4,
          align = "hv")


# Anomalies in daily correlations of "active" networks
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

TimeEdgelist$g1 <- "between"
for (dd in unique(TimeEdgelist$chip)) {
  
  inds <- strsplit(dd, split = ":")
  ind1 <- inds[[1]][1]
  ind2 <- inds[[1]][2]
  p1   <- 1 * sapply(lst, `%in%`, x = ind1)
  p2   <- 1 * sapply(lst, `%in%`, x = ind2)
  if (any(p1&p2)) {TimeEdgelist$g1[which(TimeEdgelist$chip==dd)] <- "within"}
  
}
plo.te <- TimeEdgelist %>% group_by(date, phase) %>% dplyr::summarise(n = n())
plo.te$time    <- ifelse(plo.te$phase=="day", 11, 23)
plo.te$timestp <- as.POSIXct(plo.te$date) + dhours(plo.te$time)
plot.te        <- filter(plo.te, phase != "other")
plot(plo.te$timestp, plo.te$n, type = "l")


layout <- layout_nicely(gra_full)
edges  <- data.frame()
for (i in as.character(unique(TimeEdgelist$date))) {
  exe <- TimeEdgelist[TimeEdgelist$date==i & TimeEdgelist$phase=="night",]
  gbx <- get_group_by_individual(exe$dyad, data_format = "groups")
  
  names <- colnames(gbx)
  diffs <- setdiff(colnames(gbi_full), names)
  gbx <- cbind(gbx, matrix(data=0, nrow=nrow(gbx), ncol=length(diffs)))
  colnames(gbx) <- c(names, diffs)
  gbx <- gbx[,order(colnames(gbx))]
  
  adx <- get_network(gbx, data_format = "GBI", association_index = "SRI")
  grx <- graph_from_adjacency_matrix(adx, mode="undirected", weighted=TRUE)
  temp_plot <- plot(grx, vertex.label=NA, layout = layout, main = paste0("date:", as.character(i)))
  
  r_full <- sum(unname(unmatrix(adx)!=0) & unname(unmatrix(net_full)!=0))
  
  edges <- rbind.data.frame(edges,
                            data.frame(date = as.POSIXct(i),
                                       n_edges = length(unique(exe$chip)),
                                       triangles = length(triangles(grx))/3,
                                       r_full = r_full))
  
  ggsave(temp_plot, file = paste0("~/Desktop/plot_", i,".png"), width = 14, height = 10, units = "cm")
  
}
ggplot(edges, aes(x = date)) +
  geom_line(aes(y = n_edges), color="blue", size = 1.5) +
  geom_line(aes(y = triangles), color = "red", size = 1.5) +
  geom_line(aes(y = r_full), color = "orange", size = 1.5) +
  ylab("values") + xlab("") +
  scale_color_manual(values = c("Number of edges" = "blue", 
                                "Number of triads" = "red", 
                                "Edges' placement" = "orange"))

