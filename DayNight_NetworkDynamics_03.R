###############################################################################
###################### NIGHT vs DAY ## ACTIVE vs PASSIVE ######################
###############################################################################
# Building networks per day/night, for passive("long")/active("short") contacts
dates  <- as.character(unique(TimeEdgelist$date))
phases <- c("day", "night")
n      <- length(unique(list_ind$Chip))
ind    <- unique(list_ind$Chip)

# for (i in dates) {

for (j in phases) {
  
  for (k in c("short", "long")) {
    
    #if (j=="day") {
    # d <- which(TimeEdgelist$date == i)
    #p <- which(TimeEdgelist$phase == j)
    #l <- intersect(d, p)
    #} else {
    # d1 <- which(TimeEdgelist$date ==i)
    #p1 <- which(TimeEdgelist$phase==j & TimeEdgelist$hour>16)
    #l1 <- intersect(d1, p1)
    
    # d2 <- which(TimeEdgelist$date == as.character(as.POSIXct(i, format = "%Y-%m-%d") + ddays(1)))
    #p2 <- which(TimeEdgelist$phase ==j & TimeEdgelist$hour<9)
    #l2 <- intersect(d2, p2)
    #l <- union(l1, l2)
    #}
    
    l <- which(TimeEdgelist$phase==j)
    if (k == "short") {su <- which(TimeEdgelist$duration < 25*60)}
    if (k == "long")  {su <- which(TimeEdgelist$duration >= 25*60)}
    ll <- intersect(l, su)
    list <- TimeEdgelist[ll,]
    
    if (length(ll) > 0) {
      
      Adj <- randomize_within_list(list, perm = perm, list_ind = list_ind, community = Community)
      
    }
    
    if (j==phases[1] & k=="short") {
      Arr <- Adj
      nam <- rep(paste0(j, "_", k), (perm+1))
    } else {
      Arr <- DescTools::Abind(Arr, Adj, along = 1)
      nam <- c(nam, rep(paste0(j, "_", k), (perm+1)))
    }
    
  }
  
}

#print(paste("Done: ", i))

#}
dimnames(Arr)[[1]] <- nam



#for (i in unique(dimnames(Arr)[[1]])) {

#  mats <- which(dimnames(Arr)[[1]] == i)

for (j in seq(1, (dim(Arr)[1]), (perm+1))) {
  
  k <- j
  while (k<dim(Arr)[1]) {
    
    invisible(capture.output(A <- lowerMat(Arr[j,,])))
    invisible(capture.output(B <- lowerMat(Arr[k,,])))
    r <- lsa::cosine(A,B)
    t <- 0
    for (p in 1:perm) {
      invisible(capture.output(A <- psych::lowerMat(Arr[j+p,,])))
      invisible(capture.output(B <- psych::lowerMat(Arr[k+p,,])))
      v <- lsa::cosine(A, B)
      if (v >= r) {t <- t+1}
    }
    
    if (j==1 & k==1) {
      mosaic.data <- data.frame(X = dimnames(Arr)[[1]][j],
                                Y = dimnames(Arr)[[1]][k],
                                R.cosine = r,
                                P.cosine = t/perm)
    } else {
      mosaic.data <- rbind.data.frame(mosaic.data,
                                      data.frame(X = dimnames(Arr)[[1]][j],
                                                 Y = dimnames(Arr)[[1]][k],
                                                 R.cosine = r,
                                                 P.cosine = t/perm))
    }
    
    k <- k + 101
    
  }
  
  # print(paste0("Type = ", i, ", Date = ", j))
  
}

#}
data.file2 <- mosaic.data
data.file2$R.cosine[data.file2$X==data.file2$Y] <- NA
data.file2$P.cosine[data.file2$X==data.file2$Y] <- NA
data.file2$class <- "white"
data.file2$class[data.file2$P.cosine<0.05] <- "tomato"
data.file2$class[data.file2$P.cosine>0.95] <- "cornflowerblue"
data.file2 <- rbind(data.file2,
                    data.frame(X = data.file2$Y,
                               Y = data.file2$X,
                               R.cosine = data.file2$R.cosine,
                               P.cosine = data.file2$P.cosine,
                               class = data.file2$class))
colnames(data.file2) <- c("phase1", "phase2", "r_mean", "p_value", "class")

#mosaic.data$P.transform <- mosaic.data$P.cosine
#mosaic.data$P.transform[mosaic.data$P.transform > 0.5] <- 1 - mosaic.data$P.transform[mosaic.data$P.transform > 0.5]

#data.file2 <- data.frame()
#for (i in unique(mosaic.data$X)) {

#  for (j in unique(mosaic.data$Y)) {
#    
#    r <- mosaic.data$R.cosine[mosaic.data$X==i & mosaic.data$Y==j]
#    s <- mosaic.data$P.transform[mosaic.data$X==i & mosaic.data$Y==j]
#    q <- mosaic.data$P.cosine[mosaic.data$X==i & mosaic.data$Y==j]
#    object <- MeanCI(r, na.rm = T)
#    a <- round(object[1], digits = 4) # mean coefficient of correlation
#    b <- object[2]                    # lower boundary of the CI
#    c <- object[3]                    # upper boundary of the CI
#    d <- MeanAD(r, na.rm = T)         # mean absolute deviation around the mean
#    e <- poolr::bonferroni(s)$p
#    g <- ifelse(e<0.05 & mean(q)<0.5, "cornflowerblue", ifelse(e<0.05 & mean(q)>0.5, "tomato", "white"))
#    
#    data.file2 <- rbind.data.frame(data.file2,
#                                   data.frame(phase1 = i,
#                                              phase2 = j,
#                                              r_mean = a,
#                                              ci_low = b,
#                                              ci_up  = c,
#                                              r_absd = d,
#                                              p_val  = e,
#                                              class  = g))
#    
#  }
#  
#}


# Plotting the results
# Pre-made files
data.file2 <- read.csv("~/06_pairwiseCorrelation_dayNight_activePassive.csv", header = T)
for (i in 1:nrow(data.file2)) {
  data.file2$group.label.y[i] <- ifelse(length(grep(pattern = "day", x = data.file2$phase1[i]))>0, "day", "night")
  data.file2$group.label.x[i] <- ifelse(length(grep(pattern = "day", x = data.file2$phase2[i]))>0, "day", "night")
  data.file2$tick.label.y[i] <- ifelse(length(grep(pattern = "long", x = data.file2$phase1[i]))>0, ">= 25'", "< 25'")
  data.file2$tick.label.x[i] <- ifelse(length(grep(pattern = "long", x = data.file2$phase2[i]))>0, ">= 25'", "< 25'")
}
pha <- c(paste0(phases, "_short"), paste0(phases, "_long"))
data.file2$phase1 <- factor(data.file2$phase1, pha[order(pha)])
data.file2$phase2 <- factor(data.file2$phase2, pha[order(pha)])

# plotting
a <- ggplot(data.file2, aes(x  = phase2, 
                            y    = phase1, 
                            fill = r_mean)) + 
  geom_tile(width = 0.9, 
            height = 0.9) +
  scale_fill_gradient2(low      = "cornflowerblue", 
                       high     = "tomato", 
                       mid      = "lightyellow", 
                       na.value = "white",
                       midpoint = mean(data.file2$r_mean, na.rm = T), 
                       limit    = c(min(data.file2$r_mean, na.rm = T), max(data.file2$r_mean, na.rm = T)), 
                       breaks   = c(min(data.file2$r_mean, na.rm = T), mean(data.file2$r_mean, na.rm = T), max(data.file2$r_mean, na.rm = T)), 
                       labels   = c(as.character(round(min(data.file2$r_mean, na.rm = T)*100)/100), 
                                    as.character(round(mean(data.file2$r_mean, na.rm = T)*100)/100), 
                                    as.character(round(max(data.file2$r_mean, na.rm = T)*100)/100)), 
                       space    = "Lab",
                       name     = "Cosine similarity index") +
  theme_few() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal")
b <- ggplot(data.file2, 
            aes(x    = phase2, 
                y    = phase1, 
                fill = class)) + 
  geom_tile(width = 0.9, 
            height = 0.9) +
  scale_fill_manual(breaks = c("cornflowerblue", "tomato", "white"),
                    values = c("cornflowerblue", "tomato", "white"),
                    labels = c("Lower than random", 
                               "Higher than random",
                               "Non significant"),
                    name = "Observed value is:") +
  geom_text(aes(x     = phase2, 
                y     = phase1, 
                label = ""), 
            color = "black", 
            size  = 1) +
  coord_fixed() +
  theme_few() +
  theme(
    plot.title          = element_text(hjust = 0.5,
                                       size  = 18),
    axis.text.x         = element_text(angle = 45, 
                                       vjust = 1, 
                                       size  = 15, 
                                       hjust = 1,
                                       color = "white"),
    axis.text.y         = element_text(angle = 0, 
                                       vjust = 1, 
                                       size  = 15, 
                                       hjust = 1,
                                       color = "white"),
    axis.title.x         = element_blank(),
    axis.title.y         = element_blank(),
    panel.grid.major     = element_blank(),
    panel.border         = element_blank(),
    panel.background     = element_blank(),
    axis.ticks           = element_blank(),
    legend.direction     = "vertical",
    legend.position      = "none") +
  ggtitle("Associated p-values")

df.t <- data.file2
a <- ggplot(data = df.t, 
            mapping = aes(x = interaction(tick.label.x, group.label.x), 
                          y = interaction(tick.label.y, group.label.y))) +
  geom_tile(mapping = aes(fill  = r_mean, 
                          width = 0.9, 
                          height = 0.9)) +
  geom_text(data = df.t, aes(x = interaction(tick.label.x, group.label.x), 
                             y = interaction(tick.label.y, group.label.y), 
                             label = round(r_mean*100)/100), size = 9) +
  scale_fill_gradient2(low      = "cornflowerblue", 
                       high     = "tomato", 
                       mid      = "lightyellow", 
                       na.value = "white",
                       midpoint = mean(data.file2$r_mean, na.rm = T), 
                       limit    = c(min(data.file2$r_mean, na.rm = T), max(data.file2$r_mean, na.rm = T)), 
                       breaks   = c(min(data.file2$r_mean, na.rm = T), mean(data.file2$r_mean, na.rm = T), max(data.file2$r_mean, na.rm = T)), 
                       labels   = c(as.character(round(min(data.file2$r_mean, na.rm = T)*100)/100), 
                                    as.character(round(mean(data.file2$r_mean, na.rm = T)*100)/100), 
                                    as.character(round(max(data.file2$r_mean, na.rm = T)*100)/100)), 
                       name     = "Cosine similarity index") +
  scale_x_discrete(labels = c("< 25'", "> 25'", "< 25'", "> 25'")) +
  scale_y_discrete(labels = c("< 25'", "> 25'", "< 25'", "> 25'")) +
  ylab("day                     night") +
  xlab("day                                night") +
  theme(axis.title.x  = element_text(vjust = -5, size = 28),
        axis.title.y  = element_text(vjust = 5, size = 28),
        axis.text.x   = element_text(vjust = -1, size = 25),
        axis.text.y   = element_text(hjust = -1, size = 25),
        plot.margin   = unit(c(1.5, 1.5, 1.5, 1.5), "cm"),
        legend.position = "none")




###############################################################################
##################### CORRELATIONS BETWEEN DAILY NETSWORKS ####################
###############################################################################
# Building daily "active" networks (day/night)
sho.list <- subset(TimeEdgelist, duration <= 25*60)
data.day <- subset(sho.list, day_phase %in% c("morning","afternoon"))
data.nig <- subset(sho.list, day_phase %in% c("early_night", "late_night"))
dates    <- levels(as.factor(c(data.day$date, data.nig$date)))[1:28]
perm     <- 1000

for (i in dates) {
  
  # Select the lines for day network and night network that day
  i  <- as.POSIXct(i, format = "%Y-%m-%d")
  d  <- which(data.day$date == i)
  n1 <- intersect(which(data.nig$date==i),         which(data.nig$hour>12))
  n2 <- intersect(which(data.nig$date==i+days(1)), which(data.nig$hour<12))
  
  # Observed lists of interactions
  list_day   <- data.day[d,]
  list_night <- data.nig[c(n1,n2),]
  
  # Creation of the permuted networks
  repeat {
    net_day    <- randomize_within_list(list_day, perm, list_ind, Community)
    test       <- apply(net_day, 1, sum)
    if ( !(0%in%test) ) {break;} 
  }
  repeat {
    net_night  <- randomize_within_list(list_night, perm, list_ind, Community)
    test       <- apply(net_night, 1, sum)
    if ( !(0%in%test) ) {break;}
  }
  net_day_b   <- net_day
  net_night_b <- net_night
  net_day_b[net_day_b!=0] <- 1
  net_night_b[net_night_b!=0] <- 1
  
  if (i == dates[1]) {
    Nets.mosaic.w <- Abind(net_day, net_night, along = 1)
    Nets.mosaic.b <- Abind(net_day_b, net_night_b, along = 1)
    names.dim     <- c(as.POSIXct(i + dhours(11)), as.POSIXct(i + dhours(23)))
  } else {
    Nets.mosaic.w <- Abind(Abind(Nets.mosaic.w, net_day,   along = 1), net_night,   along = 1)
    Nets.mosaic.b <- Abind(Abind(Nets.mosaic.b, net_day_b, along = 1), net_night_b, along = 1)
    names.dim     <- c(names.dim, c(as.POSIXct(i + dhours(11)), as.POSIXct(i + dhours(23))))
  }
  
}

# Calculates network similarity between all potential pairs of networks
m <- 1
for (k in seq(1, dim(Nets.mosaic.w)[1]-((2*perm)+1), (perm+1))) {
  
  n <- m + 1
  for (l in seq(k+(perm+1), dim(Nets.mosaic.w)[1]-perm, (perm+1))) {
    
    invisible(capture.output(A <- lowerMat(Nets.mosaic.w[k,,])))
    invisible(capture.output(B <- lowerMat(Nets.mosaic.w[l,,])))
    invisible(capture.output(C <- lowerMat(Nets.mosaic.b[k,,])))
    invisible(capture.output(D <- lowerMat(Nets.mosaic.b[l,,])))
    r <- lsa::cosine(A, B)
    s <- lsa::cosine(C, D)
    t <- 0
    u <- 0
    for (p in 1:perm) {
      
      invisible(capture.output(A <- lowerMat(Nets.mosaic.w[k+p,,])))
      invisible(capture.output(B <- lowerMat(Nets.mosaic.w[l+p,,])))
      invisible(capture.output(C <- lowerMat(Nets.mosaic.b[k+p,,])))
      invisible(capture.output(D <- lowerMat(Nets.mosaic.b[l+p,,])))
      v <- lsa::cosine(A, B)
      w <- lsa::cosine(C, D)
      if (v < r) {t <- t+1}
      if (w < s) {u <- u+1}
    }
    
    if (m == 1 & n == 2) {
      mosaic.data <- data.frame(X = names.dim[m],
                                Y = names.dim[n],
                                R.cosine.w = r,
                                P.cosine.w = t/perm,
                                R.cosine.b = s,
                                P.cosine.b = u/perm)
    } else {
      mosaic.data <- rbind.data.frame(mosaic.data,
                                      data.frame(X = names.dim[m],
                                                 Y = names.dim[n],
                                                 R.cosine.w = r,
                                                 P.cosine.w = t/perm,
                                                 R.cosine.b = s,
                                                 P.cosine.b = u/perm))
    }
    
    n <- n + 1
  }
  
  print(paste0("You have completed ", m, " date(s) out of 54."))
  m <- m + 1
}

mosaic.data5 <- mosaic.data
mosaic.data5$P.transform <- mosaic.data5$P.cosine.w
mosaic.data5$P.transform[mosaic.data5$P.transform > 0.5] <- 1 - mosaic.data5$P.transform[mosaic.data5$P.transform > 0.5]
mosaic.data5$P.adjust <- stats::p.adjust(mosaic.data5$P.transform, method = "fdr")
mosaic.data5$class.w <- "white"
mosaic.data5$class.w[mosaic.data5$P.transform<0.05 & mosaic.data5$P.cosine.w<0.5] <- "cornflowerblue"
mosaic.data5$class.w[mosaic.data5$P.transform<0.05 & mosaic.data5$P.cosine.w>0.5] <- "tomato"

mosaic.data5$P.transform <- mosaic.data5$P.cosine.b
mosaic.data5$P.transform[mosaic.data5$P.transform > 0.5] <- 1 - mosaic.data5$P.transform[mosaic.data5$P.transform > 0.5]
mosaic.data5$P.adjust <- stats::p.adjust(mosaic.data5$P.transform, method = "fdr")
mosaic.data5$class.b <- "white"
mosaic.data5$class.b[mosaic.data5$P.transform<0.05 & mosaic.data5$P.cosine.b<0.5] <- "cornflowerblue"
mosaic.data5$class.b[mosaic.data5$P.transform<0.05 & mosaic.data5$P.cosine.b>0.5] <- "tomato"

mosaic.data2 <- rbind.data.frame(mosaic.data5,
                                 data.frame(X = mosaic.data5[,2],
                                            Y = mosaic.data5[,1],
                                            R.cosine.w = mosaic.data5[,3],
                                            P.cosine.w = mosaic.data5[,4],
                                            R.cosine.b = mosaic.data5[,5],
                                            P.cosine.b = mosaic.data5[,6],
                                            P.transform = mosaic.data5[,7],
                                            P.adjust   = mosaic.data5[,8],
                                            class.w    = mosaic.data5[,9],
                                            class.b    = mosaic.data5[,10]))


# Summarise correlation coefficients
e.b <-c()
for (k in seq(1, dim(Nets.mosaic.w)[1]-((2*perm)+1), (perm+1))) {
  e.b <- c(e.b, sum(Nets.mosaic.b[k,,]>0)/2)
}
data.frame(mean = mean(e.b,na.rm=T), sd = sd(e.b,na.rm=T), tot = (18*17/2)+(10*9/2))

# average R2 night to day
lls <- intersect(which(as.Date(mosaic.data2$X) == as.Date(mosaic.data2$Y)-lubridate::ddays(1)),
                 which(hour(mosaic.data2$X)==23))
mean(mosaic.data2$R.cosine.w[lls])

lls <- intersect(which(as.Date(mosaic.data2$X) == as.Date(mosaic.data2$Y)+lubridate::ddays(1)),
                 which(hour(mosaic.data2$X)==11))
mean(mosaic.data2$R.cosine.w[lls])


# Plot the results
# Pre-made files
mosaic.data2 <- read.csv("~/07_pairwiseCorrelation_1month_active.csv", header = T, stringsAsFactors = F)
mosaic.data2$X <- as.POSIXct(mosaic.data2$X, format = "%Y-%m-%d %H:%M:%S")
mosaic.data2$Y <- as.POSIXct(mosaic.data2$Y, format = "%Y-%m-%d %H:%M:%S")

# plotting
a <- ggplot(mosaic.data2, 
            aes(x    = Y, 
                y    = X, 
                fill = R.cosine.w)) + 
  geom_tile() +
  scale_fill_gradient2(
    low      = "cornflowerblue", 
    high     = "tomato", 
    mid      = "lightyellow", 
    na.value = "white",
    midpoint = mean(c(mosaic.data2$R.cosine.w, mosaic.data2$R.cosine.b) , na.rm = T), 
    limit    = c(min(c(mosaic.data2$R.cosine.w, mosaic.data2$R.cosine.b), na.rm = T), 
                 max(c(mosaic.data2$R.cosine.w, mosaic.data2$R.cosine.b), na.rm = T)),
    breaks   = c(min(c(mosaic.data2$R.cosine.w, mosaic.data2$R.cosine.b), na.rm = T),
                 mean(c(mosaic.data2$R.cosine.w, mosaic.data2$R.cosine.b) , na.rm = T),
                 max(c(mosaic.data2$R.cosine.w, mosaic.data2$R.cosine.b), na.rm = T)),
    labels   = c(as.character(round(min(c(mosaic.data2$R.cosine.w, mosaic.data2$R.cosine.b), na.rm = T)*100)/100),
                 as.character(round(mean(c(mosaic.data2$R.cosine.w, mosaic.data2$R.cosine.b) , na.rm = T)*100)/100),
                 as.character(round(max(c(mosaic.data2$R.cosine.w, mosaic.data2$R.cosine.b), na.rm = T)*100)/100)),
    name     = "Similarity index") +
  geom_text(aes(x     = Y, 
                y     = X, 
                label = ""), 
            color = "black", 
            size  = 1) +
  coord_fixed() +
  theme_few() +  
  theme(
    plot.title          = element_text(hjust = 0.5,
                                       size = 18),
    axis.text.x         = element_blank(),
    axis.text.y         = element_text(angle = 0, 
                                       vjust = 1, 
                                       size  = 15, 
                                       hjust = 1,
                                       color = "black"),
    axis.title.x         = element_blank(),
    axis.title.y         = element_text(size = 20, margin = margin(r=20)),
    plot.margin          = unit(c(0,0,0,0), "cm"),
    panel.spacing        = unit(c(0,0,0,0), "cm"),
    panel.grid.major     = element_blank(),
    panel.border         = element_blank(),
    panel.background     = element_blank(),
    axis.ticks           = element_blank(),
    legend.justification = c(1, 0),
    legend.direction     = "horizontal",
    legend.text          = element_text(size = 12),
    legend.title         = element_text(size  = 15)) +
  guides(fill = guide_colorbar(barwidth       = 7, 
                               barheight      = 1,
                               title.position = "top",
                               title.hjust    = 0.5)) +
  ylab("Weighted networks") +
  ggtitle("Cosine similarity indexes")


b <- ggplot(mosaic.data2, 
            aes(x    = Y, 
                y    = X, 
                fill = class.w)) + 
  geom_tile() +
  scale_fill_manual(breaks = c("cornflowerblue", "tomato", "white"),
                    values = c("cornflowerblue", "tomato", "white"),
                    labels = c("Lower than random", 
                               "Higher than random",
                               "Non significant"),
                    name = "Observed value is:") +
  geom_text(aes(x     = Y, 
                y     = X, 
                label = ""), 
            color = "black", 
            size  = 1) +
  coord_fixed() +
  theme_few() +
  theme(
    plot.title          = element_text(hjust = 0.5,
                                       size  = 18),
    axis.text.x         = element_blank(),
    axis.text.y         = element_blank(),
    axis.title.x         = element_blank(),
    axis.title.y         = element_blank(),
    plot.margin          = unit(c(0,0,0,0), "cm"),
    panel.spacing        = unit(c(0,0,0,0), "cm"),
    panel.grid.major     = element_blank(),
    panel.border         = element_blank(),
    panel.background     = element_blank(),
    axis.ticks           = element_blank(),
    legend.direction     = "vertical",
    legend.key           = element_rect(colour = "black", size = 1)) +
  ggtitle("Associated p-values")

c <- ggplot(mosaic.data2, 
            aes(x    = Y, 
                y    = X, 
                fill = R.cosine.b)) + 
  geom_tile() +
  scale_fill_gradient2(
    low      = "cornflowerblue", 
    high     = "tomato", 
    mid      = "lightyellow", 
    na.value = "white",
    midpoint = mean(c(mosaic.data2$R.cosine.w, mosaic.data2$R.cosine.b) , na.rm = T), 
    limit    = c(min(c(mosaic.data2$R.cosine.w, mosaic.data2$R.cosine.b), na.rm = T), 
                 max(c(mosaic.data2$R.cosine.w, mosaic.data2$R.cosine.b), na.rm = T)),
    breaks   = c(min(c(mosaic.data2$R.cosine.w, mosaic.data2$R.cosine.b), na.rm = T),
                 mean(c(mosaic.data2$R.cosine.w, mosaic.data2$R.cosine.b) , na.rm = T),
                 max(c(mosaic.data2$R.cosine.w, mosaic.data2$R.cosine.b), na.rm = T)),
    labels   = c(as.character(round(min(c(mosaic.data2$R.cosine.w, mosaic.data2$R.cosine.b), na.rm = T)*100)/100),
                 as.character(round(mean(c(mosaic.data2$R.cosine.w, mosaic.data2$R.cosine.b) , na.rm = T)*100)/100),
                 as.character(round(max(c(mosaic.data2$R.cosine.w, mosaic.data2$R.cosine.b), na.rm = T)*100)/100)),
    name     = "Similarity index") +
  geom_text(aes(x     = Y, 
                y     = X, 
                label = ""), 
            color = "black", 
            size  = 1) +
  coord_fixed() +
  theme_few() +  
  theme(
    plot.title          = element_blank(),
    axis.text.x         = element_text(angle = 45, 
                                       vjust = 1, 
                                       size  = 15, 
                                       hjust = 1,
                                       color = "black"),
    axis.text.y         = element_text(angle = 0, 
                                       vjust = 1, 
                                       size  = 15, 
                                       hjust = 1,
                                       color = "black"),
    axis.title.x         = element_blank(),
    axis.title.y         = element_text(size = 20, margin = margin(r=20)),
    plot.margin          = unit(c(0,0,0,0), "cm"),
    panel.spacing        = unit(c(0,0,0,0), "cm"),
    panel.grid.major     = element_blank(),
    panel.border         = element_blank(),
    panel.background     = element_blank(),
    axis.ticks           = element_blank(),
    legend.justification = c(1, 0),
    legend.direction     = "horizontal",
    legend.position      = "none", 
    legend.text          = element_text(size = 14),
    legend.title         = element_text(size  = 18)) +
  guides(fill = guide_colorbar(barwidth       = 7, 
                               barheight      = 1,
                               title.position = "top",
                               title.hjust    = 0.5)) +
  ylab("Binary networks")


d <- ggplot(mosaic.data2, 
            aes(x    = Y, 
                y    = X, 
                fill = class.b)) + 
  geom_tile() +
  scale_fill_manual(breaks = c("cornflowerblue", "tomato", "white"),
                    values = c("cornflowerblue", "tomato", "white"),
                    labels = c("Lower than random", 
                               "Higher than random",
                               "Non significant"),
                    name = "Observed value is:") +
  geom_text(aes(x     = Y, 
                y     = X, 
                label = ""), 
            color = "black", 
            size  = 1) +
  coord_fixed() +
  theme_few() +
  theme(
    plot.title          = element_blank(),
    axis.text.x         = element_text(angle = 45, 
                                       vjust = 1, 
                                       size  = 15, 
                                       hjust = 1,
                                       color = "black"),
    axis.text.y          = element_blank(),
    axis.title.x         = element_blank(),
    axis.title.y         = element_blank(),
    plot.margin          = unit(c(0,0,0,0), "cm"),
    panel.spacing        = unit(c(-1,-1,-1,-1), "cm"),
    panel.grid.major     = element_blank(),
    panel.border         = element_blank(),
    panel.background     = element_blank(),
    axis.ticks           = element_blank(),
    legend.justification = c(1, 0),
    legend.direction     = "vertical",
    legend.position      = "none",
    legend.key           = element_rect(colour = "black", size = 1) ,
    legend.text          = element_text(size = 12),
    legend.title         = element_text(size  = 15)) +
  guides(fill = guide_legend(title.position = "top",
                             title.hjust    = 0.5))

library("patchwork")
a + b + c + d +
  plot_layout(design = "AB
              CD",
              widths = c(1, 1),
              guides = 'collect') +
  plot_annotation()




