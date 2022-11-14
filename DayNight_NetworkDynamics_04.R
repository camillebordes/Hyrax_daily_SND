#################################################################################
#################### DAYTIME vs NIGHTTIME NETWORK TRAITS ########################
#################################################################################
# Measure daily network traits on "active" networks
library("cowplot")
sho.list <- subset(TimeEdgelist, duration <= 25*60)
dates    <- levels(as.factor(sho.list$date))
perm     <- 10000
periods  <- expand.grid(date = dates, phase = c("day", "night"))
periods  <- periods[order(periods$date),]
periods  <- periods[-c(57,58),]
ind_in_l <- unique(c(TimeEdgelist$individual1, TimeEdgelist$individual2))

data.day.night <- list()
for (p in 1:(perm+1)) {
  
  d <- data.frame()
  for (i in 1:nrow(periods)) {
    
    # Calculate the social network from interaction_list
    da       <- as.POSIXct(periods$date[i], format = "%Y-%m-%d")
    ph       <- as.character(periods$phase[i])
    
    if (ph == "day") {
      ll     <- which(sho.list$date==da & sho.list$phase==ph)
    } else {
      n1 <- intersect(which(sho.list$date==da), intersect(which(sho.list$phase==ph), which(sho.list$hour>16)))
      n2 <- intersect(which(sho.list$date==(da+ddays(1))), intersect(which(sho.list$phase==ph), which(sho.list$hour<10)))
      ll <- union(n1, n2)
    }
    interact <- sho.list[ll,]
    
    gbi     <- get_group_by_individual(interact$dyad, data_format = "groups")
    invisible(
      capture.output(
        net <- get_network(association_data  =  gbi, 
                           data_format       = "GBI", 
                           association_index = 'SRI', 
                           identities        = colnames(gbi))
      )
    )
    gra <- graph_from_adjacency_matrix(net, "undirected", weighted = T)
    
    
    names    <- colnames(gbi)
    missing  <- c(setdiff(as.character(ind_in_l), names), 
                  setdiff(names, as.character(ind_in_l)))
    if (length(missing) > 0) {
      gbi1    <- cbind(gbi, 
                       matrix(0, 
                              ncol = length(missing), 
                              nrow = nrow(gbi)))
      colnames(gbi1) <- c(names, as.character(missing))
    }
    invisible(
      capture.output(
        net1 <- get_network(association_data =  gbi1, 
                            data_format       = "GBI", 
                            association_index = 'SRI', 
                            identities        = colnames(gbi1))
      )
    )
    net1 <- net1[order(colnames(net1)),] # sort the rows by names
    net1 <- net1[,order(colnames(net1))] # sort the cols by names
    gra1 <- graph_from_adjacency_matrix(net1, "undirected", weighted = T)
    
    lines         <- match(names(V(gra1)), list_ind$Chip)
    V(gra1)$date   <- as.character(as.POSIXct(periods$date[i], format = "%Y-%m-%d"))
    V(gra1)$site   <- as.character(list_ind$Canyon[lines])
    V(gra1)$phase  <- as.character(periods$phase[i])
    V(gra1)$group  <- as.character(list_ind$Group[lines])
    V(gra1)$sex    <- as.character(list_ind$Sex[lines])
    V(gra1)$weight <- as.numeric(as.character(list_ind$Weight[lines]))
    V(gra1)$smi    <- as.numeric(as.character(list_ind$SMI[lines]))
    V(gra1)$status <- as.factor(list_ind$Social_status[lines])
    
    for (v in names(V(gra1))) {
      
      j      <- as.numeric(Community$nodeclusters$cluster[which(Community$nodeclusters$node==v)])
      soc_en <- unique(Community$nodeclusters$node[which(Community$nodeclusters$cluster %in% j)])
      # node.e <- unique(names(unlist(ego(gra, order = 1, nodes = v, mode = "all"))))
      # soc_en <- unique(c(node.e, node.j))
      
      site <- V(gra1)$site[match(v, names(V(gra1)))]
      invisible(
        capture.output(
          sub_glo <- igraph::subgraph(gra1, v = names(V(gra1))[V(gra1)$site == site])
        )
      )
      sub_glo <- as.undirected(sub_glo, mode = c("collapse"), edge.attr.comb = "ignore")
      sub_adg <- as_adjacency_matrix(sub_glo, attr = "weight")
      
      invisible(
        capture.output(
          sub_loc <- igraph::subgraph(gra1, v = names(V(gra1))[names(V(gra1)) %in% soc_en])
        )
      )
      sub_loc <- as.undirected(sub_loc, mode = c("collapse"), edge.attr.comb = "ignore")
      sub_adl <- as_adjacency_matrix(sub_loc, attr = "weight")
      a       <- as_tbl_graph(sub_loc) %>% activate(edges) %>% activate(nodes)
      a       <- as_tbl_graph(sub_loc) %>% activate(edges) %>% activate(nodes) %>% 
        mutate(id_num_loc           = paste0(j, collapse = "-"),
               diff_bonds_glo       = sd(sub_adg)/mean(sub_adg),
               diff_bonds_loc       = sd(sub_adl)/mean(sub_adl),
               diff_bonds_ind       = sd(sub_adl[match(v, colnames(sub_adl)),])/mean(sub_adl[match(v, colnames(sub_adl)),]),
               size_tot             = length(names),
               size_glo             = nrow(sub_adg[rowSums(sub_adg)!=0,]),
               size_loc             = nrow(sub_adl[rowSums(sub_adl)!=0,]),
               transitivity_glo     = igraph::transitivity(sub_glo, 
                                                           type    = "globalundirected", 
                                                           weights = E(sub_glo)$weight, 
                                                           vids    = V(sub_glo)),
               transitivity_loc     = igraph::transitivity(sub_loc,   
                                                           type    = "globalundirected", 
                                                           weights = E(sub_loc)$weight,   
                                                           vids    = V(sub_loc)),
               transitivity_ind     = local_transitivity(weights = E(sub_loc)$weight),
               density_glo          = graph.density(sub_glo),
               density_loc          = graph.density(sub_loc),
               degree               = igraph::degree(sub_glo)[match(names(V(sub_loc)), names(igraph::degree(sub_glo)))],
               strength             = centrality_degree(weights = E(sub_loc)$weight, 
                                                        mode    = "all", 
                                                        loops   = FALSE,
                                                        normalized = TRUE),
               eigen                = centrality_eigen(weights  = E(sub_loc)$weight, 
                                                       directed = FALSE, 
                                                       scale    = TRUE))
      to_keep <- which(v == as.data.frame(a)$name)                                  
      if (v==names(V(gra1))[1]) {
        c <- as.data.frame(as.data.frame(a)[to_keep,])
      } else {
        c <- rbind.data.frame(c, as.data.frame(a)[to_keep,])
      }
      
    }
    
    # Calculates social differentiation
    if (i > 1) {
      for (k in 1:nrow(c)) {
        line              <- match(c$name[k], colnames(net1))
        c$assoc_stab_w[k] <- cor(x = net1[line,], y = net_less[line,], method = "spearman")
        net2 <- net1[line,]
        net3 <- net_less[line,]
        net2[net2>0] <- 1
        net3[net3>0] <- 1
        both   <- length(which(net2+net3==2))
        either <- length(which(net2+net3==2)) + length(which(net2+net3==1))
        c$assoc_stab_b[k] <- both/either
      }
    } else {
      for (k in 1:nrow(c)) {
        line              <- match(c$name[k], colnames(net1))
        c$assoc_stab_w[k] <- NA
        c$assoc_stab_b[k] <- NA
      }
    }
    d <- rbind.data.frame(d, as.data.frame(c))
    net_less <- net1
    
    # Before leaving this loop, perform one swap in the interaction list subset.
    sho.list[ll,]     <- unique_swap(interaction_list = sho.list[ll,], 
                                     list_ind         = list_ind, 
                                     community        = Community)
    print(paste0("period: ", i))
    
  }
  data.day.night <- list.append(data.day.night, d)
  print(paste0("permutation: ", p-1))
  
}
names(data.day.night) <- c("observed", paste0("perm", c(1:(p-2))))




# Overall statistical significance (individual-level network traits)
# Pre-made files are too BIG and exceed our repository size limit per file.
# Files can be retrieved by email the author: bordescamille93[at]gmail[dot]com
data.day.night <- openxlsx::loadWorkbook("~/08_networkTraits_dayNight_active.xlsx")
centralities <- c("density_loc", "diff_bonds_loc", "diff_bonds_ind",
                  "degree", "strength", "eigen")
Omnibus <- data.frame()
Graph_diffs <- data.frame()
perm <- 1000
for (central in centralities) {
  
  r.cen.d <- c()
  r.cen.n <- c()
  r.cen.diff <- c()
  r.dat.diff <- c()
  for (index in 1:(perm+1)) {
    data.provi <- openxlsx::readWorkbook(data.day.night, 
                                         sheet = index,
                                         detectDates = TRUE)
    data.provi <- data.provi[order(data.provi$name, data.provi$date, data.provi$phase),]
    data.provi$date  <- as.POSIXct(data.provi$date, "%Y-%m-%d", tz = "Asia/Jerusalem")
    data.provi$hours <- ifelse(data.provi$phase == "day", 7, 14)
    data.provi$time  <- as.POSIXct(data.provi$date + hours(data.provi$hours), origin = origin)
    col   <- match(central, colnames(data.provi))
    
    if (length(grep(pattern = "assoc_stab", x = central))>0) {
      cen.d <- data.provi[data.provi$phase=="day", col]
      cen.n <- data.provi[data.provi$phase=="night", col]
      cen.diff <- data.provi[, col]
      dat.diff <- data.provi$time
    } else {
      if (length(grep(pattern = "glo", x = central))>0) {
        data.provi <- data.provi[!duplicated(data.provi[,c("date", "phase", "site")]) ,]
        cen.d <- data.provi[data.provi$phase=="day", col]
        cen.n <- data.provi[data.provi$phase=="night", col]
      } else {
        if (length(grep(pattern = "loc", x = central))>0) {
          data.provi <- data.provi[!duplicated(data.provi[,c("date", "phase", "site", "id_num_loc")]) ,]
          cen.d <- data.provi[data.provi$phase=="day", col]
          cen.n <- data.provi[data.provi$phase=="night", col]
        } else {
          cen.d <- data.provi[data.provi$phase=="day", col]
          cen.n <- data.provi[data.provi$phase=="night", col]
        }
      }
      cen.diff <- c(lag(cen.d,1)-cen.n, cen.d-cen.n)
      dat.diff <- c(data.provi$time[data.provi$phase=="night"], data.provi$time[data.provi$phase=="day"])
    }
    
    print(paste(central, "| permutation", index))
    if (index == 1) {
      o.cen.d <- cen.d
      o.cen.n <- cen.n
      o.cen.diff <- cen.diff
      o.dat.diff <- dat.diff
      if (length(grep(pattern = "assoc_stab", x = central))>0) {
        o.cen.diff.dat <- data.frame(observed.diff = cen.diff,
                                     permuted.diff = NA,
                                     p = 0,
                                     d = data.provi$degree[data.provi$phase=="day"],
                                     n = data.provi$degree[data.provi$phase=="night"])
      } else {
        o.cen.diff.dat <- data.frame(observed.diff = abs(cen.diff), 
                                     permuted.diff = NA,
                                     p = 0,
                                     d = data.provi$degree[data.provi$phase=="day"],
                                     n = data.provi$degree[data.provi$phase=="night"])
      }
    } else {
      r.cen.d <- c(r.cen.d, cen.d)
      r.cen.n <- c(r.cen.n, cen.n)
      r.cen.diff <- c(r.cen.diff, abs(cen.diff))
      r.dat.diff <- dat.diff
      if (length(grep(pattern = "assoc_stab", x = central))>0) {
        o.cen.diff.dat$permuted.diff <- cen.diff
      } else {
        o.cen.diff.dat$permuted.diff <- abs(cen.diff)
      }
      for (ll in 1:nrow(o.cen.diff.dat)) {
        if (!is.na(o.cen.diff.dat$permuted.diff[ll]) & !is.na(o.cen.diff.dat$observed.diff[ll])) {
          if (o.cen.diff.dat$permuted.diff[ll] >= o.cen.diff.dat$observed.diff[ll]) {
            o.cen.diff.dat$p[ll] <- o.cen.diff.dat$p[ll] + 1
          }
        }
      }
    }
  }
  o.cen.diff.dat$p <- o.cen.diff.dat$p/perm
  o.cen.diff.dat$p[o.cen.diff.dat$p==0] <- 2e-16
  o.cen.diff.dat <- o.cen.diff.dat[complete.cases(o.cen.diff.dat),]
  o.cen.diff.dat <- o.cen.diff.dat[(o.cen.diff.dat$d>0 & o.cen.diff.dat$n>0),]
  hist(o.cen.diff.dat$p, main = central)
  
  #n.r  <- length(r.cen.diff)
  #r.r  <- rep(1:ceiling(n.r/perm), each = perm)[1:n.r]
  #d.r  <- split(r.cen.diff, r.r)
  #r.cen.diff  <- unname(unlist(lapply(d.r, FUN=mean, na.rm=TRUE)))
  r.dat.diff  <- as.POSIXct(r.dat.diff, origin = origin)
  
  pooled.p <- poolr::bonferroni(na.omit(o.cen.diff.dat$p), adjust = "none")
  compet.p <- CombinePValue::competitive.test(na.omit(o.cen.diff.dat$p), Weight = NA)
  Omnibus <- rbind.data.frame(Omnibus, 
                              data.frame(
                                centrality = central,
                                ave_obs_diff = mean(abs(o.cen.diff), na.rm = TRUE),
                                std_obs_diff = sd(abs(o.cen.diff), na.rm = TRUE),
                                ave_rand_diff = mean(abs(r.cen.diff), na.rm = TRUE),
                                std_rand_diff = sd(abs(r.cen.diff), na.rm = TRUE),
                                p_value_1 = unlist(unname(compet.p)),
                                p_value_2 = pooled.p$p,
                                p_value_n = pooled.p$k,
                                p_value_s = pooled.p$statistic
                              ))
  
  o.diff  <- data.frame(date = o.dat.diff, centrality = central, value = o.cen.diff, type = "obs")
  r.diff  <- data.frame(date = r.dat.diff, centrality = central, value = r.cen.diff, type = "rand")
  Graph_diffs <- rbind.data.frame(Graph_diffs, o.diff, r.diff)
  
}



# Plotting the results
# Pre-made files are too BIG and exceed our repository size limit per file.
# Files can be retreived by email the author: bordescamille93[at]gmail[dot]com
data.day.night    <- openxlsx::loadWorkbook("~/08_networkTraits_dayNight_active.xlsx")
dat <- data.provi <- openxlsx::readWorkbook(data.day.night, sheet = 1, detectDates = TRUE)

Graph_diffs       <- read.csv("~/08_networkTraits_dayNight_dailyVariations.csv", header = TRUE)[,-1]
Graph_diffs$date  <- as.POSIXct(Graph_diffs$date, "%Y-%m-%d %H:%M:%S")

centralities <- c("density_loc", "diff_bonds_loc", "diff_bonds_ind", "degree", "strength", "eigen")
significance <- c("p = 0.712", "p = 0.619", "p = 0.038", "p = 1", "p < 0.001", "p = 0.999")
y.labs.text  <- c("c) network density (p = 0.712)", "b) group differentiation (p = 0.619)", 
                  "c) individual selectivity (p = 0.038)", "a) degree centrality (p = 1)", 
                  "a) strength centrality (p < 0.001)", "b) eigenvector centrality (p = 0.999)")
brackets.h   <- c(0.4, 2.5, 1, 3.2, 0.69, 1)

for (ce in centralities) {
  
  dat$interest <- as.numeric(as.character(dat[, match(ce, colnames(dat))]))
  dat$date     <- as.POSIXct(strptime(dat$date, format = "%Y-%m-%d"))
  dat$phase    <- factor(dat$phase, c("day", "night"))
  dat$hours    <- ifelse(dat$phase == "day", 7, 14)
  dat$time     <- as.POSIXct(dat$date + hours(dat$hours), origin = origin)
  blibli <- dat %>% group_by(time, phase) %>% dplyr::summarize(mean = mean(interest, na.rm = T),
                                                               sd   = sd(interest, na.rm = T),
                                                               se   = MeanSE(interest, na.rm = T),
                                                               CI.u = MeanCI(interest, conf.level = .95, na.rm = T)[3],
                                                               CI.l = MeanCI(interest, conf.level = .95, na.rm = T)[2])
  
  dat1 <- Graph_diffs[Graph_diffs$centrality==ce,]
  pp  <- ggplot(dat1, aes(x = date,
                          y = abs(value),
                          color = type)) +
    stat_summary(fun.data = mean_sdl, 
                 fun.args = list(mult = 1),
                 alpha = 0.5) +
    geom_smooth(se = FALSE) + 
    scale_color_discrete(labels = c("observed", "random")) +
    ggtitle(significance[match(ce, centralities)]) +
    theme_pubr()
  l.pp <- get_legend(pp + theme(legend.position= "bottom",
                                legend.title = element_blank()))
  
  sp <- ggplot(blibli, aes(x = time, 
                           y = mean, 
                           color = phase)) +
    geom_line(color = "grey80") +
    geom_errorbar(aes(ymin = mean - se,
                      ymax = mean + se), 
                  width = 0) +
    geom_point() +
    scale_color_manual(labels = c("day", "night"),
                       values = c("#E69F00", "#56B4E9")) + 
    ylab("") +
    theme_pubr()
  l.sp <- get_legend(sp + theme(legend.position = "bottom",
                                legend.title = element_blank()))
  
  yplot <- ggplot(blibli, aes(y = mean, fill = phase)) +
    geom_density(alpha = 0.3) + 
    scale_fill_manual(labels = c("day", "night"),
                      values = c("#E69F00", "#56B4E9")) +
    theme_pubr()
  
  
  if (ce %in% c("diff_bonds_ind", "density_loc")) {
    
    sp    <- sp + 
      theme(plot.margin = unit(c(0.4, 0.4, 0.4, 0.4), "cm"),
            title = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y  = element_text(size = 11),
            axis.title.x = element_blank(),
            axis.text.x  = element_text(size = 11),
            legend.position = "none",
            legend.title = element_blank()) + panel_border(remove = TRUE)
    
    yplot <- yplot + 
      theme(plot.margin = unit(c(0.4, 0.4, 0.4, 0.4), "cm"),
            title = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y  = element_blank(),
            axis.line.y  = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.line.x  = element_blank(),
            axis.text.x  = element_blank(),
            legend.position = "none",
            legend.title = element_blank()) + panel_border(remove = TRUE)
    
    pp <- pp + 
      theme(plot.margin = unit(c(0.4, 0.4, 0.4, 0.4), "cm"),
            title = element_text(size = 10, hjust = -1),
            axis.title.y = element_blank(),
            axis.text.y  = element_text(size = 11),
            axis.title.x = element_blank(),
            axis.text.x  = element_text(size = 11),
            legend.position = "none",
            legend.title = element_blank()) + panel_border(remove = TRUE)
    
  } else {
    
    if (ce %in% c("strength", "degree")) {
      
      sp    <- sp + 
        theme(plot.margin = unit(c(0.4, 0.4, 0.4, 0.4), "cm"),
              title = element_blank(),
              axis.title.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.x  = element_blank(),
              axis.line.x  = element_blank(),
              axis.text.y  = element_text(size = 11),
              axis.title.y = element_blank(),
              legend.position = "none",
              legend.title = element_blank()) + panel_border(remove = TRUE)
      
      yplot <- yplot + 
        theme(plot.margin = unit(c(0.4, 0.4, 0.4, 0.4), "cm"),
              title = element_blank(),
              axis.title.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.y  = element_blank(),
              axis.line.y  = element_blank(),
              axis.title.x = element_blank(),
              axis.text.x  = element_blank(),
              axis.ticks.x = element_blank(),
              axis.line.x  = element_blank(),
              legend.position = "none",
              legend.title = element_blank()) + panel_border(remove = TRUE)
      
      pp <- pp +
        theme(plot.margin = unit(c(0.4, 0.4, 0.4, 0.4), "cm"),
              title = element_text(size = 10, hjust = -1),
              axis.title.x = element_blank(),
              axis.text.x  = element_blank(),
              axis.ticks.x = element_blank(),
              axis.line.x  = element_blank(),
              axis.text.y  = element_text(size = 11),
              axis.title.y = element_blank(),
              legend.position = "none",
              legend.title = element_blank()) + panel_border(remove = TRUE)
      
    } else {
      
      sp    <- sp + 
        theme(plot.margin = unit(c(0.4, 0.4, 0.4, 0.4), "cm"),
              title = element_blank(),
              axis.title.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.x  = element_blank(),
              axis.line.x  = element_blank(),
              axis.text.y  = element_text(size = 11),
              axis.title.y = element_blank(),
              legend.position = "none",
              legend.title = element_blank()) + panel_border(remove = TRUE)
      
      yplot <- yplot + 
        theme(plot.margin = unit(c(0.4, 0.4, 0.4, 0.4), "cm"),
              title = element_blank(),
              axis.title.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.y  = element_blank(),
              axis.line.y  = element_blank(),
              axis.title.x = element_blank(),
              axis.text.x  = element_blank(),
              axis.ticks.x = element_blank(),
              axis.line.x  = element_blank(),
              legend.position = "none",
              legend.title = element_blank()) + panel_border(remove = TRUE)
      
      pp <- pp +
        theme(plot.margin = unit(c(0.4, 0.4, 0.4, 0.4), "cm"),
              title = element_text(size = 10, hjust = -1),
              axis.title.x = element_blank(),
              axis.text.x  = element_blank(),
              axis.ticks.x = element_blank(),
              axis.line.x  = element_blank(),
              axis.text.y  = element_text(size = 11),
              axis.title.y = element_blank(),
              legend.position = "none",
              legend.title = element_blank()) + panel_border(remove = TRUE)
      
    }
    
  }
  
  assign(x = as.character(paste0("sp.", ce)),    value = sp)
  assign(x = as.character(paste0("yplot.", ce)), value = yplot)
  assign(x = as.character(paste0(ce, ".plot")),  value = pp)
  
}

# annotated composite plots
rects <- data.frame(x = c(1:2), text = c('Absolute Day-Night variations', 'Observed daily network traits'))
p.t1  <- ggplot(rects, aes(x = 0, y = 0)) + 
  geom_tile(fill = "white", height = .1,) + 
  geom_text(aes(label = "Absolute Day-Night variations"), hjust = 0.2, size = 5) + 
  scale_fill_identity(guide = "none") + 
  theme_void()
p.t2  <- ggplot(rects, aes(x = 0, y = 0)) + 
  geom_tile(fill = "white", height = .1,) + 
  geom_text(aes(label = "Observed daily network traits"), hjust = 0.2, size = 5) + 
  scale_fill_identity(guide = "none") + 
  theme_void()
p.t1                + p.t2              + patchwork::plot_spacer() +
  strength.plot       + sp.strength       + yplot.strength + 
  eigen.plot          + sp.eigen          + yplot.eigen + 
  diff_bonds_ind.plot + sp.diff_bonds_ind + yplot.diff_bonds_ind + 
  l.pp                + l.sp              + patchwork::plot_spacer() +
  patchwork::plot_layout(ncol = 3,
                         nrow = 5,
                         widths = c(3, 3, 0.5),
                         heights = c(2, 5, 5, 5, 1))
p.t1                + p.t2              + patchwork::plot_spacer() +
  degree.plot         + sp.degree         + yplot.degree + 
  diff_bonds_loc.plot + sp.diff_bonds_loc + yplot.diff_bonds_loc + 
  density_loc.plot    + sp.density_loc    + yplot.density_loc + 
  l.pp                + l.sp              + patchwork::plot_spacer() +
  patchwork::plot_layout(ncol = 3,
                         nrow = 5,
                         widths = c(3, 3, 0.5),
                         heights = c(2, 5, 5, 5, 1))


# non-annotated composite plots
strength.plot       + sp.strength       + yplot.strength +
  eigen.plot          + sp.eigen          + yplot.eigen + 
  diff_bonds_ind.plot + sp.diff_bonds_ind + yplot.diff_bonds_ind + 
  patchwork::plot_layout(ncol = 3,
                         nrow = 3,
                         widths  = c(3, 3, 0.5),
                         heights = c(5, 5, 5),
                         guides = "collect")
degree.plot         + sp.degree         + yplot.degree + 
  diff_bonds_loc.plot + sp.diff_bonds_loc + yplot.diff_bonds_loc + 
  density_loc.plot    + sp.density_loc    + yplot.density_loc + 
  patchwork::plot_layout(ncol = 3,
                         nrow = 3,
                         widths  = c(3, 3, 0.5),
                         heights = c(5, 5, 5),
                         guides = "collect")



# Overall statistical significance (group-level standard deviation in network traits)
centralities <- c("sd.eigen", "sd.degree", "sd.strength")
Omnibus <- data.frame()
Ds <- c()
Dats <- data.frame()
Difs <- data.frame()
Graph_diffs <- data.frame()
perm <- 1000

for (central in centralities) {
  
  r.cen.d <- c()
  r.cen.n <- c()
  r.cen.diff <- c()
  for (index in 1:(perm+1)) {
    data.provi <- openxlsx::readWorkbook(data.day.night, 
                                         sheet = index,
                                         detectDates = TRUE)
    data.provi <- data.provi[order(data.provi$name, data.provi$date, data.provi$phase),]
    if (central == "sd.eigen") {
      data.prov <- data.provi %>% 
        group_by(date, phase, id_num_loc) %>%
        dplyr::summarize(sd = sd(eigen, na.rm = T))
    }
    if (central == "sd.degree") {
      data.prov <- data.provi %>% 
        group_by(date, phase, id_num_loc) %>%
        dplyr::summarize(sd = sd(degree, na.rm = T))
    }
    if (central == "sd.strength") {
      data.prov <- data.provi %>% 
        group_by(date, phase, id_num_loc) %>%
        dplyr::summarize(sd = sd(strength, na.rm = T))
    }
    if (central == "sd.diff_bonds_ind") {
      data.prov <- data.provi %>% 
        group_by(date, phase, id_num_loc) %>%
        dplyr::summarize(sd = sd(diff_bonds_ind, na.rm = T))
    }
    #data.prov <- data.prov %>%
    #             group_by(date, phase) %>%
    #             dplyr::summarise(sd.m = mean(sd, na.rm = TRUE))
    cen.d <- data.prov$sd[data.prov$phase=="day"]
    cen.n <- data.prov$sd[data.prov$phase=="night"]
    cen.diff <- c(lag(cen.d,1)-cen.n, cen.d-cen.n)
    
    print(paste(central, "| permutation", index))
    if (index == 1) {
      o.cen.diff <- cen.diff
      o.cen.diff.dat <- data.frame(observed.diff = abs(cen.diff), 
                                   permuted.diff = NA,
                                   p = 0)
    } else {
      r.cen.d <- c(r.cen.d, cen.d)
      r.cen.n <- c(r.cen.n, cen.n)
      r.cen.diff <- c(r.cen.diff, cen.diff)
      o.cen.diff.dat$permuted.diff <- abs(cen.diff)
      for (ll in 1:nrow(o.cen.diff.dat)) {
        if (!is.na(o.cen.diff.dat$permuted.diff[ll]) & !is.na(o.cen.diff.dat$observed.diff[ll])) {
          if (o.cen.diff.dat$permuted.diff[ll] >= o.cen.diff.dat$observed.diff[ll]) {
            o.cen.diff.dat$p[ll] <- o.cen.diff.dat$p[ll] + 1
          }
        }
      }
    }
  }
  o.cen.diff.dat$p <- o.cen.diff.dat$p/perm
  #n.r  <- length(r.cen.diff)
  #r.r  <- rep(1:ceiling(n.r/perm), each = perm)[1:n.r]
  #d.r  <- split(r.cen.diff, r.r)
  #r.cen.diff  <- lapply(d.r, 1, mean, na.rm=TRUE) 
  
  o.diff  <- data.frame(centrality = central, value = o.cen.diff, type = "obs")
  r.diff  <- data.frame(centrality = central, value = r.cen.diff, type = "rand")
  Difs    <- rbind.data.frame(Difs, o.diff, r.diff)
  
  o.cen.diff.dat$p[o.cen.diff.dat$p==0] <- 2e-16
  o.cen.diff.dat <- o.cen.diff.dat[complete.cases(o.cen.diff.dat),]
  
  pooled.p <- poolr::bonferroni(na.omit(o.cen.diff.dat$p), adjust = "none")
  compet.p <- CombinePValue::competitive.test(na.omit(o.cen.diff.dat$p), Weight = NA)
  Omnibus <- rbind.data.frame(Omnibus, 
                              data.frame(
                                centrality = central,
                                ave_obs_diff = mean(abs(o.cen.diff), na.rm = TRUE),
                                std_obs_diff = sd(abs(o.cen.diff), na.rm = TRUE),
                                ave_rand_diff = mean(abs(r.cen.diff), na.rm = TRUE),
                                std_rand_diff = sd(abs(r.cen.diff), na.rm = TRUE),
                                p_value_1 = unlist(unname(compet.p)),
                                p_value_2 = pooled.p$p,
                                p_value_n = pooled.p$k,
                                p_value_s = pooled.p$statistic
                              ))
  Graph_diffs <- rbind.data.frame(Graph_diffs, Difs)
  
}

# Plotting the results
dat.sd <- dat %>% 
  group_by(date, phase, id_num_loc) %>%
  dplyr::summarize(sd.eigen = sd(eigen, na.rm = T),
                   sd.degree = sd(degree, na.rm = T),
                   sd.strength = sd(strength, na.rm = T)) %>%
  as.data.frame()
centralities <- c("sd.degree", "sd.strength", "sd.eigen")
significance <- c("p = 1", "p < 0.001", "p < 0.001")
y.labs.text  <- c("degree centrality", "strength centrality", "eigenvector centrality")
brackets.h   <- c(2.4, 0.9, 1)
for (ce in centralities) {
  
  dat.sd$interest <- as.numeric(as.character(dat.sd[, match(ce, colnames(dat.sd))]))
  dat.sd$date     <- as.POSIXct(strptime(dat.sd$date, format = "%Y-%m-%d"))
  dat.sd$phase    <- factor(dat.sd$phase, c("day", "night"))
  dat.sd$hours    <- ifelse(dat.sd$phase == "day", 7, 14)
  dat.sd$time     <- as.POSIXct(dat.sd$date + hours(dat.sd$hours), origin = origin)
  
  sp <- ggplot(dat.sd, aes(x = time, 
                           y = interest, 
                           fill = phase,
                           color = phase,
                           group = interaction(time, phase))) +
    geom_violin() +
    ylab(y.labs.text[match(ce, centralities)]) +
    theme_pubr()
  
  yplot <- ggpubr::ggboxplot(dat.sd, 
                             y = "interest",
                             x = "phase",
                             fill = "phase") +
    scale_fill_manual(labels = c("day", "night"),
                      values = c("#E69F00", "#56B4E9")) + 
    geom_bracket(xmin = 0.95, 
                 xmax = 2.05, 
                 y.position = brackets.h[match(ce, centralities)],
                 label = as.character(significance[match(ce, centralities)]),
                 label.size = 4,
                 tip.length = c(0, 0),
                 size = 0.9) +
    ylab(as.character(y.labs.text[match(ce, centralities)]))
  
  if (ce %in% c("sd.select")) {
    sp    <- sp + 
      rremove("legend") + 
      theme(plot.margin = unit(c(0.2,0.5,0.2,0.5), "cm"),
            axis.title.y = element_text(size = 18,
                                        vjust = 3),
            axis.text.y  = element_text(size = 15),
            axis.title.x = element_text(size = 18,
                                        vjust = -1),
            axis.text.x  = element_text(size = 15))
    yplot <- yplot + 
      rremove("legend") + 
      theme(plot.margin = unit(c(0.2,0.5,0.5,0.5), "cm"),
            axis.title.y = element_text(size = 15, color = "black"),
            axis.title.x = element_blank(),
            axis.text.x  = element_text(size = 15, color = "black"))
  } else {
    sp    <- sp + 
      rremove("legend") + 
      theme(plot.margin = unit(c(0.2,0.5,0.2,0.5), "cm"),
            axis.title.x  = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x  = element_blank(),
            axis.line.x  = element_blank(),
            axis.text.y  = element_text(size = 15),
            axis.title.y = element_text(size = 18, vjust = 3))
    yplot <- yplot + 
      rremove("legend") + 
      theme(plot.margin = unit(c(0.2,0.5,0.5,0.5), "cm"),
            axis.title.y = element_text(size = 15, color = "black"),
            axis.title.x = element_blank(),
            axis.text.x  = element_text(size = 15, color = "black"))
  }
  assign(x = as.character(paste0("sp.", ce)), value = sp)
  assign(x = as.character(paste0("yplot.", ce)), value = yplot)
  
}

plot_grid(sp.sd.degree, yplot.sd.degree,
          sp.sd.strength, yplot.sd.strength,
          sp.sd.eigen, yplot.sd.eigen,
          ncol = 2,
          nrow = 3,
          align = "h")
plot_grid(yplot.sd.degree, 
          yplot.sd.strength, 
          yplot.sd.eigen,
          ncol = 3, 
          nrow = 1, 
          align = "h")

dat.dat <- dat %>% 
  group_by(phase) %>%
  dplyr::summarise(eigen    = mean(eigen, na.rm = TRUE),
                   sd.eig = sd(eigen, na.rm = TRUE),
                   degree   = mean(degree, na.rm = TRUE),
                   sd.deg = sd(degree, na.rm = TRUE),
                   strength = mean(strength, na.rm = TRUE),
                   sd.stren = sd(strength, na.rm = TRUE),
                   diff_ind = mean(diff_bonds_ind, na.rm = TRUE),
                   sd.diffi = sd(diff_bonds_ind, na.rm = TRUE),
                   diff_loc = mean(diff_bonds_loc, na.rm = TRUE),
                   sd_diffl = sd(diff_bonds_loc, na.rm = TRUE),
                   dens_loc = mean(density_loc, na.rm = TRUE),
                   sd.densl = sd(density_loc, na.rm = TRUE)) %>%
  as.data.frame()

