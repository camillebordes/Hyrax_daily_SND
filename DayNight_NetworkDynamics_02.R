################################################################################
################################ CUSTOM FUNCTIONS ##############################
################################################################################
randomize_within_list <- function(interaction_list, perm, list_ind, community) {
  
  interact <- interaction_list
  GBI      <- get_group_by_individual(interact$dyad, data_format = "groups")
  names    <- colnames(GBI)
  missing  <- c(setdiff(as.character(list_ind$Chip), names), 
                setdiff(names, as.character(list_ind$Chip)))
  if (length(missing) > 0) {
    GBI           <- cbind(GBI, matrix(0, 
                                       ncol = length(missing), 
                                       nrow = nrow(GBI)))
    colnames(GBI) <- c(names, as.character(missing))
  }
  invisible(
    capture.output(
      Adjacency <- get_network(
        association_data  =  GBI, 
        data_format       = "GBI", 
        association_index = 'SRI', 
        identities        = colnames(GBI)
      )
    )
  )
  Adjacency   <- Adjacency[order(colnames(Adjacency)),] # sort the rows by names
  Adjacency   <- Adjacency[,order(colnames(Adjacency))] # sort the cols by names
  nline       <- nrow(Adjacency)
  A           <- array(0, c(perm+1, nline, nline))
  A[1,,]      <- Adjacency
  
  int         <- interaction_list
  for (iteration in 1:perm) {
    
    # Iterations are limited within each community
    for (com in unique(community$nodeclusters$cluster)) {
      
      swap_id <- as.numeric(community$nodeclusters$node[which(community$nodeclusters$cluster==com)])
      lin_co1 <- which(int$individual1%in%swap_id)
      lin_co2 <- which(int$individual2%in%swap_id)
      lin_com <- union(lin_co1, lin_co2)
      total_c <- length(unique(c(int$individual1[lin_com], int$individual2[lin_com])))
      first_c <- length(unique(int$individual1[lin_com]))
      stop_if <- 0
      
      if (total_c > 3) {
        if (first_c > 1) {
          repeat {
            focals <- sample(as.character(unique(int$individual1[lin_com])), 2, replace = F)
            ind1   <- as.character(focals[1])
            ind2   <- as.character(focals[2])
            if (ind1!=ind2) {
              cand3  <- unique(int$individual2[int$individual1==ind1 & int$individual2!=ind2])
              lines3 <- intersect(which(int$individual1==ind1), which(int$individual2%in%cand3))
              if (length(lines3)>0) {
                l1     <- sample(lines3, 1, replace = F)
                ind3   <- as.character(int$individual2[l1])
                cand4  <- unique(int$individual2[int$individual1==ind2 & int$individual2!=ind1])
                lines4 <- intersect(which(int$individual1==ind2), which(int$individual2%in%cand4))
                lines4 <- intersect(which(int$individual2!=ind3), lines4)
                if (length(lines4)>0) {
                  l2     <- sample(lines4, 1, replace = F)
                  ind4   <- as.character(int$individual2[l2])
                  if (!(ind3%in%c(ind1, ind2, ind4)) & !(ind4%in%c(ind1, ind2, ind3)) & int$site[l1]==int$site[l2]) {
                    break;
                  }
                }
              }
            }
            stop_if <- stop_if + 1
            if (stop_if==1000) {break}
          }
        } else {
          repeat {
            focals <- rep(unique(int$individual1[lin_com]), 2)
            ind1   <- as.character(focals[1])
            ind2   <- as.character(focals[2])
            cand3  <- unique(int$individual2[int$individual1==ind1 & int$individual2!=ind2])
            lines3 <- intersect(which(int$individual1==ind1), which(int$individual2%in%cand3))
            if (length(lines3)>0) {
              l1     <- sample(lines3, 1, replace = F)
              ind3   <- as.character(int$individual2[l1])
              cand4  <- unique(int$individual2[int$individual1==ind2 & int$individual2!=ind1])
              lines4 <- intersect(which(int$individual1==ind2), which(int$individual2%in%cand4))
              lines4 <- intersect(which(int$individual2!=ind3), lines4)
              if(length(lines4)>0) {
                l2     <- sample(lines4, 1, replace = F)
                ind4   <- as.character(int$individual2[l2])
                if (ind1!=ind4 & ind2!=ind3 & ind3!=ind4 & int$site[l1]==int$site[l2]) {
                  break;
                }
              }
            }
            stop_if <- stop_if + 1
            if (stop_if==1000) {break}
          }
        }
        if (stop_if==1000 & length(l2)==0) {
          l2   <- l1
          ind4 <- ind3
        }
        sl1 <- int$start[l1]
        i1l1<- int$individual1[l1]
        i2l1<- int$individual2[l1]
        all1<- intersect(which(int$start==sl1), 
                         intersect(which(int$individual1==i1l1), 
                                   which(int$individual2==i2l1)
                         )
        )
        sl2 <- int$start[l2]
        i1l2<- int$individual1[l2]
        i2l2<- int$individual2[l2]
        all2<- intersect(which(int$start==sl2), 
                         intersect(which(int$individual1==i1l2), 
                                   which(int$individual2==i2l2)
                         )
        )
        int$individual2[all1] <- ind4
        int$individual2[all2] <- ind3
      }
      
    }
    int$dyad  <- strsplit(paste0(int$individual1, ":", int$individual2), ":")
    
    interact <- int
    GB       <- get_group_by_individual(interact$dyad, data_format = "groups")
    names    <- colnames(GB)
    missing  <- c(setdiff(as.character(list_ind$Chip), names), 
                  setdiff(names, as.character(list_ind$Chip)))
    if (length(missing) > 0) {
      GB           <- cbind(GB, matrix(0, ncol = length(missing), nrow = nrow(GB)))
      colnames(GB) <- c(names, as.character(missing))
    }
    invisible(
      capture.output(
        Adjacence <- get_network(
          association_data  = GB, 
          data_format       = "GBI", 
          association_index = 'SRI', 
          identities        = colnames(GB)
        )
      )
    )
    Adjacence <- Adjacence[order(colnames(Adjacence)),] # sort the rows by names
    Adjacence <- Adjacence[,order(colnames(Adjacence))] # sort the cols by names
    A[(iteration+1),,] <- Adjacence
    
  }
  
  dimnames(A) <- list(c(1:dim(A)[1]), colnames(Adjacence), colnames(Adjacence))
  return(A)
  
}

unique_swap <- function(interaction_list, list_ind, community) {
  
  # Iterative data-stream permutations
  int         <- interaction_list
  
  # Iterations are limited within each canyon
  for (com in unique(community$nodeclusters$cluster)) {
    
    swap_id <- as.numeric(community$nodeclusters$node[which(community$nodeclusters$cluster==com)])
    lin_co1 <- which(int$individual1%in%swap_id)
    lin_co2 <- which(int$individual2%in%swap_id)
    lin_com <- union(lin_co1, lin_co2)
    total_c <- length(unique(c(int$individual1[lin_com], int$individual2[lin_com])))
    first_c <- length(unique(int$individual1[lin_com]))
    stop_if <- 0
    
    if (total_c > 3) {
      if (first_c > 1) {
        repeat {
          focals <- sample(as.character(unique(int$individual1[lin_com])), 2, replace = F)
          ind1   <- as.character(focals[1])
          ind2   <- as.character(focals[2])
          if (ind1!=ind2) {
            cand3  <- unique(int$individual2[int$individual1==ind1 & int$individual2!=ind2])
            lines3 <- intersect(which(int$individual1==ind1), which(int$individual2%in%cand3))
            if (length(lines3)>0) {
              l1     <- sample(lines3, 1, replace = F)
              ind3   <- as.character(int$individual2[l1])
              cand4  <- unique(int$individual2[int$individual1==ind2 & int$individual2!=ind1])
              lines4 <- intersect(which(int$individual1==ind2), which(int$individual2%in%cand4))
              lines4 <- intersect(which(int$individual2!=ind3), lines4)
              if (length(lines4)>0) {
                l2     <- sample(lines4, 1, replace = F)
                ind4   <- as.character(int$individual2[l2])
                if (!(ind3%in%c(ind1, ind2, ind4)) & !(ind4%in%c(ind1, ind2, ind3)) & int$site[l1]==int$site[l2]) {
                  break;
                }
              }
            }
          }
          stop_if <- stop_if + 1
          if (stop_if==1000) {break}
        }
      } else {
        repeat {
          focals <- rep(unique(int$individual1[lin_com]), 2)
          ind1   <- as.character(focals[1])
          ind2   <- as.character(focals[2])
          cand3  <- unique(int$individual2[int$individual1==ind1 & int$individual2!=ind2])
          lines3 <- intersect(which(int$individual1==ind1), which(int$individual2%in%cand3))
          if (length(lines3)>0) {
            l1     <- sample(lines3, 1, replace = F)
            ind3   <- as.character(int$individual2[l1])
            cand4  <- unique(int$individual2[int$individual1==ind2 & int$individual2!=ind1])
            lines4 <- intersect(which(int$individual1==ind2), which(int$individual2%in%cand4))
            lines4 <- intersect(which(int$individual2!=ind3), lines4)
            if(length(lines4)>0) {
              l2     <- sample(lines4, 1, replace = F)
              ind4   <- as.character(int$individual2[l2])
              if (ind1!=ind4 & ind2!=ind3 & ind3!=ind4 & int$site[l1]==int$site[l2]) {
                break;
              }
            }
          }
          stop_if <- stop_if + 1
          if (stop_if==1000) {break}
        }
      }
      if (stop_if==1000 & length(l2)==0) {
        l2   <- l1
        ind4 <- ind3
      }
      sl1 <- int$start[l1]
      i1l1<- int$individual1[l1]
      i2l1<- int$individual2[l1]
      all1<- intersect(which(int$start==sl1), 
                       intersect(which(int$individual1==i1l1), 
                                 which(int$individual2==i2l1)
                       )
      )
      sl2 <- int$start[l2]
      i1l2<- int$individual1[l2]
      i2l2<- int$individual2[l2]
      all2<- intersect(which(int$start==sl2), 
                       intersect(which(int$individual1==i1l2), 
                                 which(int$individual2==i2l2)
                       )
      )
      int$individual2[all1] <- ind4
      int$individual2[all2] <- ind3
      int$dyad[all1] <- paste0(min(int$individual1[l1], ind4), ":", max(int$individual1[l1], ind4))
      int$dyad[all2] <- paste0(min(int$individual1[l2], ind3), ":", max(int$individual1[l2], ind3))
      
    }
    
  }
  
  int$dyad <- strsplit(paste0(int$individual1, ":", int$individual2), ":")
  return(int)
  
}

random_net <- function(Array, perm) {
  
  t                <- 0
  diag(Array[1,,]) <- NA
  sd_obs           <-   sd(Array[1,,], na.rm = T)
  mean_obs         <- mean(Array[1,,], na.rm = T)
  CV_observed      <- 100*(sd_obs/mean_obs)
  for (index in 1:perm) {
    diag(Array[(index+1),,]) <- NA
    sd_rand        <-   sd(Array[(index+1),,], na.rm = T)
    mean_rand      <- mean(Array[(index+1),,], na.rm = T)
    CV_random      <- 100*(sd_rand/mean_rand)
    if (sum(CV_random, na.rm = T)>0 & sum(CV_observed, na.rm = T)>0) {
      if (CV_random >= CV_observed) {t <- t+1}
    }
  }
  
  return(t/perm)
  
}


################################################################################
######################### ACTIVE vs PASSIVE CONTACTS ###########################
################################################################################
# Network filtered by encounter length
perm     <- 1000
for (m in seq(0, 100, 5)) {
  
  encounters <- TimeEdgelist[data.table::between(TimeEdgelist$duration,(m*60) ,(max(TimeEdgelist$duration)*60)) ,]
  encounters <- encounters[encounters$phase=="night" ,]
  
  if (dim(encounters)[1] != 0) {
    
    net <- randomize_within_list(encounters, perm = perm, list_ind = list_ind, community = Community)
    gri <- graph_from_adjacency_matrix(net[1,,], mode = "undirected", weighted = T)
    
    if (m==0) {
      Nets.mosaic <- net
      names.dim <- c(m)
      average.s <- c(mean(unname(igraph::strength(gri)), na.rm = T))
      std.err.s <- c(MeanSE(unname(igraph::strength(gri)), na.rm = T))
      std.dev.s <- c(unname(igraph::strength(gri)))
    } else {
      Nets.mosaic <- DescTools::Abind(Nets.mosaic, net, along = 1)
      names.dim <- c(names.dim, m)
      average.s <- c(average.s, mean(unname(igraph::strength(gri)), na.rm = T))
      std.err.s <- c(std.err.s, MeanSE(unname(igraph::strength(gri)), na.rm = T))
      std.dev.s <- c(std.dev.s, unname(igraph::strength(gri)))
    }
    
    print(paste("Done: ", m))
    
  }
  
}


# Pairwise cosine similarity index between networks
m <- 1
for (k in seq(1, dim(Nets.mosaic)[1]-((2*perm)+1), (perm+1))) {
  
  n <- m + 1
  for (l in seq(k+(perm+1), dim(Nets.mosaic)[1]-perm, (perm+1))) {
    
    invisible(capture.output(A <- psych::lowerMat(Nets.mosaic[k,,])))
    invisible(capture.output(B <- psych::lowerMat(Nets.mosaic[l,,])))
    r <- ecodist::mantel(A ~ B)[1]
    s <- lsa::cosine(A, B)
    
    t1 <- 0
    t2 <- 0
    for (p in 1:perm) {
      invisible(capture.output(A <- psych::lowerMat(Nets.mosaic[k+p,,])))
      invisible(capture.output(B <- psych::lowerMat(Nets.mosaic[l+p,,])))
      u <- ecodist::mantel(A ~ B)[1]
      v <- lsa::cosine(A, B)
      if (u < r) {t1 <- t1+1}
      if (v < s) {t2 <- t2+1}
    }
    
    if (m == 1 & n == 2) {
      mosaic.data <- data.frame(X = names.dim[m],
                                Y = names.dim[n],
                                R.mantel = r,
                                P.mantel = t1/perm,
                                R.cosine = s,
                                P.cosine = t2/perm)
    } else {
      mosaic.data <- rbind.data.frame(mosaic.data,
                                      data.frame(X = names.dim[m],
                                                 Y = names.dim[n],
                                                 R.mantel = r,
                                                 P.mantel = t1/perm,
                                                 R.cosine = s,
                                                 P.cosine = t2/perm))
    }
    
    n <- n + 1
    
  }
  
  print(paste0("You have completed ", m, "/50 window scans."))
  m <- m + 1
  
}
mosaic.data5 <- mosaic.data
mosaic.data5$P.transform <- mosaic.data5$P.cosine
mosaic.data5$P.transform[mosaic.data5$P.transform > 0.5] <- 1 - mosaic.data5$P.transform[mosaic.data5$P.transform > 0.5]
mosaic.data5$P.adjust <- stats::p.adjust(mosaic.data5$P.transform, method = "fdr")
mosaic.data5$class.w <- "white"
mosaic.data5$class.w[mosaic.data5$P.transform<0.05 & mosaic.data5$P.cosine<0.5] <- "cornflowerblue"
mosaic.data5$class.w[mosaic.data5$P.transform<0.05 & mosaic.data5$P.cosine>0.5] <- "tomato"
mosaic.data1 <- rbind.data.frame(mosaic.data5,
                                 data.frame(X = mosaic.data5[,2],
                                            Y = mosaic.data5[,1],
                                            R.mantel = mosaic.data5[,3],
                                            P.mantel = mosaic.data5[,4],
                                            R.cosine = mosaic.data5[,5],
                                            P.cosine = mosaic.data5[,6],
                                            P.transform = mosaic.data5[,7],
                                            P.adjust   = mosaic.data5[,8],
                                            class.w    = mosaic.data5[,9]))


# Average strength per network filtered on encounter lengths
filter.strength <- data.frame(window = seq(0, 100, 5), 
                              mean_s = average.s, 
                              vari_s = std.err.s,
                              std_s  = std.dev.s)


################################################################################
################################ PLOTTING RESULTS ##############################
################################################################################
# Pre-made files
mosaic.data1    <- read.csv("~/05_pairwiseCorrelation_activePassive.csv", header = T)
filter.strength <- read.csv("~/05_averageStrength_activePassive.csv", header = T)

# Plots
a <- ggplot(mosaic.data1, 
            aes(x    = Y, 
                y    = X, 
                fill = R.cosine)) + 
  geom_tile(width = 4, 
            height = 4) +
  geom_rect(aes(xmin = 0-2.5, xmax = 25-2.5,
                ymin = 0-2.5, ymax = 25-2.5),
            color = "red",
            fill = NA) +
  geom_rect(aes(xmin = 25-2.5, xmax = 100+2.5,
                ymin = 25-2.5, ymax = 100+2.5),
            color = "red",
            fill = NA) +
  scale_fill_gradient2(low      = "cornflowerblue", 
                       high     = "tomato", 
                       mid      = "lightyellow", 
                       na.value = "white",
                       midpoint = mean( mosaic.data1$R.cosine, na.rm = T), 
                       limit    = c(min(mosaic.data1$R.cosine, na.rm = T), 
                                    max(mosaic.data1$R.cosine, na.rm = T)),
                       breaks   = c(min(mosaic.data1$R.cosine, na.rm = T),
                                    mean( mosaic.data1$R.cosine, na.rm = T),
                                    max(mosaic.data1$R.cosine, na.rm = T)),
                       labels   = c(as.character(round(min(mosaic.data1$R.cosine, na.rm = T)*100)/100),
                                    as.character(round(mean( mosaic.data1$R.cosine, na.rm = T)*100)/100),
                                    as.character(round(max(mosaic.data1$R.cosine, na.rm = T)*100)/100)),
                       name     = "") +
  geom_text(aes(x     = Y, 
                y     = X, 
                label = ""), 
            color = "black", 
            size  = 1) +
  coord_fixed() +
  ylab("cut-off (minutes)") +
  xlab("cut-off (minutes)") +
  theme_few() + 
  theme(
    plot.title          = element_blank(),
    axis.text.x         = element_text(angle = 0, 
                                       size  = 15,
                                       color = "black"),
    axis.text.y         = element_text(angle = 0,
                                       size  = 15,
                                       color = "black"),
    axis.title.x         = element_text(vjust = -3, size = 18),
    axis.title.y         = element_text(angle = 90,
                                        color = "black",
                                        vjust = 5,
                                        size = 18),
    panel.grid.major     = element_blank(),
    panel.border         = element_blank(),
    panel.background     = element_blank(),
    plot.margin          = unit(c(0, 1, 0.5, 0), "cm"),
    axis.ticks           = element_blank(),
    legend.justification = c(1, 0),
    legend.direction     = "horizontal",
    legend.position      = "bottom",
    legend.text          = element_text(size = 12),
    legend.title         = element_text(size  = 15)) +
  guides(fill = guide_colorbar(barwidth       = 7, 
                               barheight      = 1,
                               title.position = "left",
                               title.hjust    = 0.5))


b <- ggplot(mosaic.data1, 
            aes(x    = Y, 
                y    = X, 
                fill = class.w)) + 
  geom_tile() +
  geom_vline(xintercept = 25-2.5,
             col = "darkred") +
  geom_hline(yintercept = 25-2.5,
             color = "darkred") +
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

c <- ggplot(filter.strength, aes(x = window, 
                                 y = mean_s)) +
  geom_errorbar(aes(x = window,
                    ymin = mean_s - vari_s,
                    ymax = mean_s + vari_s),
                width = 0) +
  geom_point() +
  ylab("strength centrality") +
  theme_pubr() +
  theme(axis.title.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.line.x  = element_blank(), 
        axis.title.y = element_text(vjust = 5, size = 18),
        plot.margin  = unit(c(2,1,-0.5,1), "cm"))

ggarrange(c, a, ncol = 1, nrow = 2, align = "v", heights = c(1, 3.4), legend = "bottom")




