
# Consts ----------------------------------------------------------------
paper_output_path <- "~/Dropbox (BGU)/Apps/Overleaf/Rumen microbiome coocurrence/"


# HPC ---------------------------------------------------------------------

write_to_log <- function(x, exp_id, JOB_ID, level, level_name, append=T){
  write_lines(x, paste(exp_id,JOB_ID,level,level_name,'log.txt',sep='_'), append = append)
}

parse_networks_HPC <- function(){
  # Create a node list from the original matrix
  nodes <- tibble(node_id=1:nrow(ASV_occurrence_mat_cow), node_name=sort(rownames(ASV_occurrence_mat_cow)))
  
  # The positive network. All ASVs are included but interactions are only the POSITIVE ones
  print('Building positive coocurrence matrix')
  # Represent network as matrix to account for singletons
  network_pos <- matrix(0, nrow(nodes), nrow(nodes), dimnames = list(nodes$node_name, nodes$node_name))
  realized_interactions <- ASV_co %>% filter(edge_type=='pos') %>% select(sp1_name, sp2_name, weight)
  g <- graph.data.frame(realized_interactions, directed = F)
  print(graph.density(g, loops = F))
  realized_interactions <- as_adjacency_matrix(g, names = T, sparse = F, attr = 'weight')
  network_pos[rownames(realized_interactions),colnames(realized_interactions)] <- realized_interactions
  # There could be singletons
  singletons_pos <- names(which(colSums(network_pos)==0))
  # Represent network as edge list
  edge_list_pos <- as_tibble(igraph::as_data_frame(g, what = 'edges'))
  
  # The negative network. All ASVs are included but interactions are only the NEGATIVE ones
  print('Building negative coocurrence matrix')
  # Represent network as matrix to account for possible singletons
  network_neg <- matrix(0, nrow(nodes), nrow(nodes), dimnames = list(nodes$node_name, nodes$node_name))
  realized_interactions <- ASV_co %>% filter(edge_type=='neg') %>% select(sp1_name, sp2_name, weight)
  g <- graph.data.frame(realized_interactions, directed = F)
  print(graph.density(g, loops = F))
  realized_interactions <- as_adjacency_matrix(g, names = T, sparse = F, attr = 'weight')
  network_neg[rownames(realized_interactions),colnames(realized_interactions)] <- realized_interactions
  # There could be singletons
  singletons_neg <- names(which(colSums(network_neg)==0))
  # Represent network as edge list
  edge_list_neg <- as_tibble(igraph::as_data_frame(g, what = 'edges'))
  
  singletons <- tibble(ASV_ID=c(singletons_pos,singletons_neg), edge_type=c(rep('pos',length(singletons_pos)),rep('neg',length(singletons_neg))))
  edge_list <- bind_rows(edge_list_pos %>% mutate(edge_type='pos'), edge_list_neg %>% mutate(edge_type='neg'))
  out <- list(nodes=nodes, ASV_co=ASV_co, mat_pos=network_pos, mat_neg=network_neg, singletons=singletons, edge_list=edge_list)
  return(out)
}

parse_networks <- function(e_id, Level, Level_name){
  JOB_ID <- subset(as.data.frame(run_summary), exp_id==e_id & level==Level & level_name==Level_name)$JOB_ID
  # Get nodes
  nodes <- suppressMessages(read_csv(paste(e_id,JOB_ID,Level,Level_name,"nodes.csv",sep = "_")))
  print(paste('Nodes (microbes):',nrow(nodes)))
  # Get singletons
  singletons <- suppressMessages(read_csv(paste(e_id,JOB_ID,Level,Level_name,"singletons.csv",sep = "_")))
  # Get edge list
  edge_list <- suppressMessages(read_csv(paste(e_id,JOB_ID,Level,Level_name,"edge_list.csv",sep = "_")))
  # Get positive matrix
  mat_pos <- read.csv(paste(e_id,JOB_ID,Level,Level_name,"mat_pos.csv",sep = "_"))
  mat_pos <- data.matrix(mat_pos)
  rownames(mat_pos) <- colnames(mat_pos)
  # Get negative matrix
  mat_neg <- read.csv(paste(e_id,JOB_ID,Level,Level_name,"mat_neg.csv",sep = "_"))
  mat_neg <- data.matrix(mat_neg)
  rownames(mat_neg) <- colnames(mat_neg)
  print(paste('Density of positive network:', round(density_uni(mat_pos),2)))
  print(paste('Density of negative network:', round(density_uni(mat_neg),2)))
  
  # Get edge lists to igraph
  g_pos <- graph.data.frame(edge_list %>% filter(edge_type=='pos'), directed = F)
  print(paste('Density of positive network (without singletons):', round(graph.density(g_pos, loops = F),2)))
  g_neg <- graph.data.frame(edge_list %>% filter(edge_type=='neg'), directed = F)
  print(paste('Density of negative network (without singletons):', round(graph.density(g_neg, loops = F),2)))
  
  print('Singletons: ')
  print(table(singletons$edge_type))
  
  print('# of edges: ')
  print(table(edge_list$edge_type))
  
  out <- list(nodes=nodes, edge_list=edge_list, g_pos=g_pos, g_neg=g_neg)
  return(out)
}

parse_networks_from_cooc <- function(e_id, Level="Farm", Level_name="UK1", pos_only=TRUE){
  JOB_ID <- subset(as.data.frame(run_summary), exp_id==e_id & level==Level & level_name==Level_name)$JOB_ID
  # Get nodes
  nodes <- suppressMessages(read_csv(paste(e_id,JOB_ID,Level,Level_name,"nodes.csv",sep = "_")))
  print(paste('Nodes (microbes):',nrow(nodes)))
  # Get singletons
  singletons <- suppressMessages(read_csv(paste(e_id,JOB_ID,Level,Level_name,"singletons.csv",sep = "_")))
  # Get edge list
  cooc_pos_list <- suppressMessages(read_csv(paste(e_id,JOB_ID,Level,Level_name,"COOC.csv",sep = "_")))
  edge_list <- cooc_pos_list %>% 
        select(from=sp1_name, to=sp2_name, weight, p_val=p_gt, edge_type, level, level_name, everything())
  if (pos_only){
    edge_list <- edge_list %>% filter(edge_type=="pos")
  } 
  
  print('Singletons: ')
  print(table(singletons$edge_type))
  
  print('# of positive edges: ')
  print(table(edge_list$edge_type))
  nodes <- NULL
  out <- list(nodes=nodes, edge_list=edge_list)
  return(out)
}


# Density -----------------------------------------------------------------

density_bip <- function(x){
  sum(x)/(ncol(x)*nrow(x))
}

density_uni <- function(x){
  x <- 1*(x>0)
  2*sum(x)/(nrow(x)*(nrow(x)-1))
}


# Plotting ----------------------------------------------------------------

paper_figs_theme <- 
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black",fill = NA,size = 1),
        panel.spacing = unit(0.5, "cm", data = NULL),
        axis.text = element_text(size=14, color='black'),
        axis.title = element_text(size=14, color='black'),
        axis.line = element_blank())
paper_figs_theme_no_legend <- 
  paper_figs_theme +
  theme(legend.position = 'none')


html_figs_theme <- 
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black",fill = NA,size = 1),
        panel.spacing = unit(0.5, "cm", data = NULL),
        axis.text = element_text(size=10, color='black'),
        axis.title = element_text(size=10, color='black'),
        axis.line = element_blank())
html_figs_theme_no_legend <- 
  html_figs_theme +
  theme(legend.position = 'none')

# shuffling by scale ---------------------------

shuffle_farm_microbe <- function(data_file_name, nsim=100, output_folder="shuff_farm_005") {
  asv_data <- read_csv(data_file_name)
  
  # This permutes the microbes between cows WITHIN a farm. Microbes are not
  # shuffled between farms.
  for (i in 1:nsim){
    print(i)
    asv_data %>% 
      group_by(Farm) %>% # Within a farm
      mutate(ASV_ID=sample(ASV_ID, replace = F)) %>% 
      ungroup() %>% 
      mutate(shuff_id=i) %>% 
      write_csv(paste(output_folder,'/shuff_farm_',str_pad(i, 3, '0',side='left'),'.csv', sep=''))
  }
  write_csv(
    tibble(e_id=1:nsim, data_file=paste('shuff_farm_',str_pad(1:nsim, 3, '0',side='left'),'.csv', sep='')),
    paste(output_folder,'/experiments.csv', sep=''))
}

shuffle_region_microbe <- function(data_file_name, nsim=100, output_folder="shuff_farm_005") {
  asv_data <- read_csv(data_file_name)
  
  # This permutes the microbes between cows between farms. Microbes are 
  # shuffled between farms.
  for (i in 1:nsim){
    print(i)
    asv_data %>% 
      mutate(ASV_ID=sample(ASV_ID, replace = F)) %>% 
      ungroup() %>% 
      mutate(shuff_id=i) %>% 
      write_csv(paste(output_folder,'/shuff_farm_',str_pad(i, 3, '0',side='left'),'.csv', sep=''))
  }
  write_csv(
    tibble(e_id=1:nsim, data_file=paste('shuff_farm_',str_pad(1:nsim, 3, '0',side='left'),'.csv', sep='')),
    paste(output_folder,'/experiments.csv', sep=''))
}

suffle_cow_microb_vegan <- function(data_file_name, nsim=500, output_folder="shuff_vegan_30", shuff_method='curveball'){
  asv_data <- read_csv(data_file_name)
  
  # get layers names
  farm_names <- asv_data %>% select(Farm) %>% distinct() %>% pull(Farm)

  layers <- list()
  # make shuffled matrices -
  # have to separate to farms in order to not shuffle between them
  for (frm in farm_names) {
    lyr <- asv_data %>% filter(Farm==frm)

    # convert edgelist to matrix
    lyr_mat <- lyr %>%
      select(Cow_Code, ASV_ID) %>% add_column(weight=1) %>%
      dcast(Cow_Code ~ ASV_ID, value.var = "weight", fill = 0) %>%
      column_to_rownames(var="Cow_Code")

    # run curveball shuffling
    null <- vegan::nullmodel(lyr_mat, method = shuff_method)
    suff <- simulate(null, nsim = nsim, burnin = 5000, seed = 1234)

    layers[[frm]] <- suff
  }

  # make files from shuffled
  for (i in 1:nsim) {
    all_farms <- tibble(Country = character(),
                        Farm = character(),
                        Cow_Code = character(),
                        ASV_ID = character(),
                        Abundance = numeric(),
                        shuff_id = numeric())
    # convert multi-mats to "edge lists"
    for (frm in farm_names) {
       matt <- layers[[frm]][,,i]
       # convert
       edge <- melt(matt) %>% filter(value == 1) %>%
               select(Cow_Code=Var1, ASV_ID=Var2, Abundance=value) %>%
               add_column(shuff_id=i) %>%
               add_column(Farm=frm, .before = 1) %>%
               add_column(Country=substr(frm, 1, 2), .before = 1)
       all_farms <- rbind(all_farms, edge)
    }
    # save to file
    write_csv(all_farms, paste(output_folder,'/shuff_farm_',str_pad(i, 3, '0',side='left'),'.csv', sep=''))
  }

  write_csv(tibble(e_id=1:nsim, 
                   data_file=paste('shuff_farm_',str_pad(1:nsim, 3, '0',side='left'),'.csv', sep=''),
                   Abundance_file=basename(data_file_name)),
            paste(output_folder,'/experiments.csv', sep=''))
}

# intra layer building functions -------
cooccur <- function (mat, type = "spp_site", thresh = TRUE, spp_names = FALSE, 
         true_rand_classifier = 0.1, prob = "hyper", site_mask = NULL, 
         only_effects = FALSE, eff_standard = TRUE, eff_matrix = FALSE) 
{
  if (type == "spp_site") {
    spp_site_mat <- mat
  }
  if (type == "site_spp") {
    spp_site_mat <- t(mat)
  }
  if (spp_names == TRUE) {
    spp_key <- data.frame(num = 1:nrow(spp_site_mat), spp = row.names(spp_site_mat))
  }
  if (!is.null(site_mask)) {
    if (nrow(site_mask) == nrow(spp_site_mat) & ncol(site_mask) == 
        ncol(spp_site_mat)) {
      N_matrix <- create.N.matrix(site_mask)
    }
    else {
      stop("Incorrect dimensions for site_mask, aborting.")
    }
  }
  else {
    site_mask <- matrix(data = 1, nrow = nrow(spp_site_mat), 
                        ncol = ncol(spp_site_mat))
    N_matrix <- matrix(data = ncol(spp_site_mat), nrow = nrow(spp_site_mat), 
                       ncol = nrow(spp_site_mat))
  }
  spp_site_mat[spp_site_mat > 0] <- 1
  tsites <- ncol(spp_site_mat)
  nspp <- nrow(spp_site_mat)
  spp_pairs <- choose(nspp, 2)
  incidence <- prob_occur <- obs_cooccur <- prob_cooccur <- exp_cooccur <- matrix(nrow = spp_pairs, 
                                                                                  ncol = 3)
  incidence <- prob_occur <- matrix(nrow = nrow(N_matrix), 
                                    ncol = ncol(N_matrix))
  for (spp in 1:nspp) {
    if (spp < nspp) {
      for (spp_next in (spp + 1):nspp) {
        incidence[spp, spp_next] <- sum(site_mask[spp, 
        ] * site_mask[spp_next, ] * mat[spp, ])
        incidence[spp_next, spp] <- sum(site_mask[spp, 
        ] * site_mask[spp_next, ] * mat[spp_next, ])
      }
    }
  }
  prob_occur <- incidence/N_matrix
  pb <- txtProgressBar(min = 0, max = (nspp + nrow(obs_cooccur)), 
                       style = 3)
  row <- 0
  for (spp in 1:nspp) {
    if (spp < nspp) {
      for (spp_next in (spp + 1):nspp) {
        pairs <- sum(as.numeric(mat[spp, site_mask[spp, 
        ] * site_mask[spp_next, ] == 1] == 1 & mat[spp_next, 
                                                   site_mask[spp, ] * site_mask[spp_next, ] == 
                                                     1] == 1))
        row <- row + 1
        obs_cooccur[row, 1] <- spp
        obs_cooccur[row, 2] <- spp_next
        obs_cooccur[row, 3] <- pairs
        prob_cooccur[row, 1] <- spp
        prob_cooccur[row, 2] <- spp_next
        prob_cooccur[row, 3] <- prob_occur[spp, spp_next] * 
          prob_occur[spp_next, spp]
        exp_cooccur[row, 1] <- spp
        exp_cooccur[row, 2] <- spp_next
        exp_cooccur[row, 3] <- prob_cooccur[row, 3] * 
          N_matrix[spp, spp_next]
      }
    }
    setTxtProgressBar(pb, spp)
  }
  if (thresh == TRUE) {
    n_pairs <- nrow(prob_cooccur)
    prob_cooccur <- prob_cooccur[exp_cooccur[, 3] >= 1, ]
    obs_cooccur <- obs_cooccur[exp_cooccur[, 3] >= 1, ]
    exp_cooccur <- exp_cooccur[exp_cooccur[, 3] >= 1, ]
    n_omitted <- n_pairs - nrow(prob_cooccur)
    pb <- txtProgressBar(min = 0, max = (nspp + nrow(obs_cooccur)), 
                         style = 3)
  }
  output <- data.frame(matrix(nrow = 0, ncol = 9))
  colnames(output) <- c("sp1", "sp2", "sp1_inc", "sp2_inc", 
                        "obs_cooccur", "prob_cooccur", "exp_cooccur", "p_lt", 
                        "p_gt")
  for (row in 1:nrow(obs_cooccur)) {
    sp1 <- obs_cooccur[row, 1]
    sp2 <- obs_cooccur[row, 2]
    sp1_inc <- incidence[sp1, sp2]
    sp2_inc <- incidence[sp2, sp1]
    max_inc <- max(sp1_inc, sp2_inc)
    min_inc <- min(sp1_inc, sp2_inc)
    nsite <- N_matrix[sp1, sp2]
    psite <- as.numeric(nsite + 1)
    prob_share_site <- rep(x = 0, times = psite)
    if (prob == "hyper") {
      if (only_effects == FALSE) {
        all.probs <- phyper(0:min_inc, min_inc, nsite - 
                              min_inc, max_inc)
        prob_share_site[1] <- all.probs[1]
        for (j in 2:length(all.probs)) {
          prob_share_site[j] <- all.probs[j] - all.probs[j - 
                                                           1]
        }
      }
      else {
        for (j in 0:nsite) {
          if ((sp1_inc + sp2_inc) <= (nsite + j)) {
            if (j <= min_inc) {
              prob_share_site[(j + 1)] <- 1
            }
          }
        }
      }
    }
    if (prob == "comb") {
      if (only_effects == FALSE) {
        for (j in 0:nsite) {
          if ((sp1_inc + sp2_inc) <= (nsite + j)) {
            if (j <= min_inc) {
              prob_share_site[(j + 1)] <- coprob(max_inc = max_inc, 
                                                 j = j, min_inc = min_inc, nsite = nsite)
            }
          }
        }
      }
      else {
        for (j in 0:nsite) {
          if ((sp1_inc + sp2_inc) <= (nsite + j)) {
            if (j <= min_inc) {
              prob_share_site[(j + 1)] <- 1
            }
          }
        }
      }
    }
    p_lt <- 0
    p_gt <- 0
    for (j in 0:nsite) {
      if (j <= obs_cooccur[row, 3]) {
        p_lt <- prob_share_site[(j + 1)] + p_lt
      }
      if (j >= obs_cooccur[row, 3]) {
        p_gt <- prob_share_site[(j + 1)] + p_gt
      }
      if (j == obs_cooccur[row, 3]) {
        p_exactly_obs <- prob_share_site[(j + 1)]
      }
    }
    p_exactly_obs <- round(p_exactly_obs, 5)
    prob_cooccur[row, 3] <- round(prob_cooccur[row, 3], 3)
    exp_cooccur[row, 3] <- round(exp_cooccur[row, 3], 1)
    output[row, ] <- c(sp1, sp2, sp1_inc, sp2_inc, obs_cooccur[row, 
                                                               3], prob_cooccur[row, 3], exp_cooccur[row, 3], p_lt, 
                       p_gt)
    setTxtProgressBar(pb, nspp + row)
  }
  close(pb)
  if (spp_names == TRUE) {
    sp1_name <- merge(x = data.frame(order = 1:length(output$sp1), 
                                     sp1 = output$sp1), y = spp_key, by.x = "sp1", by.y = "num", 
                      all.x = T, sort = FALSE)
    sp2_name <- merge(x = data.frame(order = 1:length(output$sp2), 
                                     sp2 = output$sp2), y = spp_key, by.x = "sp2", by.y = "num", 
                      all.x = T, sort = FALSE)
    output$sp1_name <- sp1_name[with(sp1_name, order(order)), 
                                "spp"]
    output$sp2_name <- sp2_name[with(sp2_name, order(order)), 
                                "spp"]
  }
  true_rand <- (nrow(output[(output$p_gt >= 0.05 & output$p_lt >= 
                               0.05) & (abs(output$obs_cooccur - output$exp_cooccur) <= 
                                          (tsites * true_rand_classifier)), ]))
  output_list <- list(call = match.call(), results = output, 
                      positive = nrow(output[output$p_gt < 0.05, ]), negative = nrow(output[output$p_lt < 
                                                                                              0.05, ]), co_occurrences = (nrow(output[output$p_gt < 
                                                                                                                                        0.05 | output$p_lt < 0.05, ])), pairs = nrow(output), 
                      random = true_rand, unclassifiable = nrow(output) - (true_rand + 
                                                                             nrow(output[output$p_gt < 0.05, ]) + nrow(output[output$p_lt < 
                                                                                                                                0.05, ])), sites = N_matrix, species = nspp, percent_sig = (((nrow(output[output$p_gt < 
                                                                                                                                                                                                            0.05 | output$p_lt < 0.05, ])))/(nrow(output))) * 
                        100, true_rand_classifier = true_rand_classifier)
  if (spp_names == TRUE) {
    output_list$spp_key <- spp_key
    output_list$spp.names = row.names(spp_site_mat)
  }
  else {
    output_list$spp.names = c(1:nrow(spp_site_mat))
  }
  if (thresh == TRUE) {
    output_list$omitted <- n_omitted
    output_list$pot_pairs <- n_pairs
  }
  class(output_list) <- "cooccur"
  if (only_effects == F) {
    output_list
  }
  else {
    effect.sizes(mod = output_list, standardized = eff_standard, 
                 matrix = eff_matrix)
  }
}


# interlayer building functions ---------------------------

# Interlayer edges value: proportion difference in <abundance> between farms for each ASV
interlayer_definition <- function(x){1-ifelse(x[1]/x[2]<1,
                                              x[1]/x[2],
                                              x[2]/x[1])}

# calculate the interlayer edges and return a tibble with an edge list.
calc_abund_inter_ASV <- function(y){
  node_interlayer_edges <- combn(y$rel_abund,2, FUN=interlayer_definition)
  # make an output matrix
  # out <- matrix(0, nrow(y), nrow(y), dimnames = list(y$Farm, y$Farm))
  # out[lower.tri(out)] <- interlayer_edges
  # Or make an output list
  out <- tibble(
    layer_from=combn(y$Farm, 2)[1,], 
    layer_to=combn(y$Farm, 2)[2,], 
    weight=node_interlayer_edges)
  return(out)
}

# calculate the interlayer edges for nodes degree and return a tibble with an edge list
calc_deg_inter_ASV <- function(z){
  node_interlayer_edges <- combn(z$node_degree,2, FUN=interlayer_definition)
  # make an output matrix
  # out <- matrix(0, nrow(z), nrow(z), dimnames = list(z$Farm, z$Farm))
  # out[lower.tri(out)] <- interlayer_edges
  # Or make an output list
  out <- tibble(
    layer_from=combn(z$layer_id, 2)[1,], 
    layer_to=combn(z$layer_id, 2)[2,], 
    weight=node_interlayer_edges)
  return(out)
}

# analysis on nodes ---------------------------

# Function to calculate clustering coefficient
calc_CC_local <- function(x){
  g <- graph_from_data_frame(x, directed = FALSE, vertices = NULL)
  CC=transitivity(g, type = 'local')
  tibble(ASV_ID=V(g)$name, CC=CC, k=degree(g))
}

# # This function calculates SIMILARITY
calculate_PF_J <- function(x) {
  mat_ASV=
    x %>%
    group_by(to) %>%
    select(c(to,level_name)) %>%
    mutate(present=1) %>%
    spread(to, present, fill = 0) %>%
    column_to_rownames("level_name")
  beta_ASV <- 1-vegdist(mat_ASV, "jaccard") # Convert to similarity
  PF_J <- mean(beta_ASV)
  PF_J_sd <- sd(beta_ASV)
  num_layers <- nrow(as.matrix(beta_ASV))
  out <- data.frame(PF_J, PF_J_sd, num_layers=num_layers)
  return(out)
}

# This function calculates DISSIMILARITY
calculate_PF_U <-  function(x, tree) {
  # prune the tree
  included_asvs <- unique(x$to)
  unincluded <- tree$tip.label[!tree$tip.label %in% included_asvs]
  pruned <- dendextend::prune(tree, unincluded)
  
  
  mat_format <- x %>%
    group_by(to) %>%
    select(c(to,level_name)) %>%
    mutate(present=1) %>%
    spread(to, present, fill = 0) %>%
    column_to_rownames("level_name")
  
  # run unifrec
  unifracs <- GUniFrac(mat_format, pruned, alpha=c(0, 0.5, 1))$unifracs
  du <- unifracs[, , "d_UW"]
  n_habitats <- nrow(unifracs)
  
  unifrac_summ <- NULL
  for (d in c("d_1", "d_UW", "d_0", "d_0.5")){
    d_vals <- unifracs[, , d]
    du_low <- d_vals[lower.tri(d_vals)]
    
    d_Mean <- mean(du_low)
    d_SD <- sd(du_low)
    num_habitats <- nrow(unifracs)
    unifrac_summ <- rbind(unifrac_summ, data.frame(PF_U=d_Mean, 
                                                   PF_U_sd=d_SD, 
                                                   num_layers=n_habitats, 
                                                   unifrec_type=d))
  }
  
  return(unifrac_summ)
}


# running Infomap -------------------------------------------------
run_infomap_multilayer_multilevel <- function (M, infomap_executable = "Infomap", flow_model = NULL, 
                                               two_level=T,
                                               silent = T, trials = 100, seed = NULL, relax = F, multilayer_relax_rate = 0.1, 
                                               multilayer_relax_limit = NULL, multilayer_relax_limit_up = NULL, 
                                               multilayer_relax_limit_down = NULL, temporal_network = F, 
                                               run_standalone = T, remove_auxilary_files = T, ...) {
  if (check_infomap(infomap_executable) == F) {
    stop("Error in Infomap stand-alone file.")
  }
  if (class(M) != "multilayer") {
    stop("M must be of class multilayer")
  }
  arguments <- paste(" --tree -N ", trials, sep = "")
  arguments <- ifelse(!is.null(seed), paste(arguments, "--seed", seed), arguments)
  arguments <- ifelse(!is.null(flow_model), paste(arguments, "-f", flow_model), arguments)
  arguments <- ifelse(silent, paste(arguments, "--silent"), arguments)
  arguments <- ifelse(two_level, paste(arguments, "-2"), arguments)
  arguments <- paste(arguments, ...)
  if (relax == F) {
    print("Using interlayer edge values to determine flow between layers.")
    write_lines("*Multilayer", "infomap_multilayer.txt")
    write_delim(M$intra, "infomap_multilayer.txt", delim = " ", 
                append = T)
    write_delim(M$inter, "infomap_multilayer.txt", delim = " ", 
                append = T)
  }  else {
    if (ncol(M$intra) == 5) {
      stop("Cannot use relax rates with extended format of intralayer edges. See function create_multilayer_object.")
    }
    print("Using global relax to determine flow between layers.")
    write_lines("*Intra", "infomap_multilayer.txt")
    write_delim(M$intra, "infomap_multilayer.txt", delim = " ", append = T)
    if (!is.null(M$inter)) {
      if (ncol(M$inter) == 5) {
        stop("Cannot use relax rates with extended format of interlayer edges. See function create_multilayer_object.")
      }
      print("Global relax will be constrained by interlayer edges.")
      write_lines("*Inter", "infomap_multilayer.txt", 
                  append = T)
      write_delim(M$inter, "infomap_multilayer.txt", delim = " ", 
                  append = T)
    }
    arguments <- ifelse(!is.null(multilayer_relax_rate), 
                        paste(arguments, "--multilayer-relax-rate", multilayer_relax_rate), 
                        arguments)
    arguments <- ifelse(!is.null(multilayer_relax_limit), 
                        paste(arguments, "--multilayer-relax-limit", multilayer_relax_limit), 
                        arguments)
    arguments <- ifelse(!is.null(multilayer_relax_limit_up), 
                        paste(arguments, "--multilayer-relax-limit-up", 
                              multilayer_relax_limit_up), arguments)
    arguments <- ifelse(!is.null(multilayer_relax_limit_down), 
                        paste(arguments, "--multilayer-relax-limit-down", 
                              multilayer_relax_limit_down), arguments)
  }
  
  call <- paste("./", infomap_executable, " infomap_multilayer.txt . ", 
                arguments, sep = "")
  if (run_standalone == T) {
    print(call)
    system(call)
  } else {
    print("Please run Infomap online at https://www.mapequation.org/infomap/ using the following arguments (copy-paste):")
    print(arguments)
    invisible(readline(prompt = "After running, download statenodes results and press [ENTER] when done"))
    if (!file.exists("network_states.tree")) {
      stop("Result file network_states.tree was not found. Did you download results?")
    }
    file.rename(from = "network_states.tree", to = "infomap_multilayer_states.tree")
  }
  L_output <- parse_number(read_lines("infomap_multilayer_states.tree")[6])
  modules <- suppressMessages(read_delim("infomap_multilayer_states.tree", 
                                         delim = " ", 
                                         skip = 11, 
                                         col_names = c("path", "flow", "name", "state_id", "node_id", "layer_id")))
  
  modules %<>% filter(flow > 0) 
  num_levels <- max(str_count(modules$path, pattern = ':'))
  level_columns <- c(paste("level",1:num_levels,sep=''), 'leaf_id')
  modules %<>%
    select(path, node_id, layer_id, flow) %>% 
    separate(path, into = level_columns, sep = ":") %>% 
    mutate(node_id=as.integer(node_id)) %>%
    mutate(leaf_id=as.integer(leaf_id)) %>%
    mutate(layer_id=as.integer(layer_id)) %>%
    mutate_if(is.character, as.integer) %>%
    full_join(M$nodes, "node_id") %>%
    select(node_id, starts_with("level"), everything()) %>% arrange(node_id, layer_id)
  
  if (temporal_network) {
    print("Reorganizing modules...")
    renamed_moduels <- modules %>% distinct(module, layer_id) %>% 
      arrange(module, layer_id)
    x <- c(1, table(renamed_moduels$module))
    module_birth_layers <- renamed_moduels %>% slice(cumsum(x)) %>% 
      arrange(layer_id, module)
    module_renaming <- data.frame(module = module_birth_layers$module, 
                                  module_renamed = 1:max(module_birth_layers$module))
    modules %<>% left_join(module_renaming, "module") %>% 
      select(-module) %>% rename(module = module_renamed)
  }
  if (remove_auxilary_files) {
    print("Removing auxilary files...")
    file.remove("infomap_multilayer_states.tree")
    file.remove("infomap_multilayer.txt")
    file.remove("infomap_multilayer.tree")
  }
  print(paste("Partitioned into ", max(modules$level1), " top modules.", 
              sep = ""))
  out <- list(call = call, L = L_output, m = max(modules$level1), seed=seed,
              modules = modules)
  class(out) <- "infomap_multilayer"
  return(out)
}
