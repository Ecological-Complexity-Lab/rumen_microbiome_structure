library(attempt)

# -- temporary run infomap from functions instread of package
create_multilayer_object <- function(extended=NULL, intra=NULL, inter=NULL, nodes=NULL, intra_output_extended=T, layers=NULL, write_to_file=F, filename_prefix=NULL){
  if(!is.null(layers)){
    if(names(layers)[1]!='layer_id') {stop('First column in layers must be named layer_id')}
  }
  if(names(nodes)[1]!='node_id') {stop('First column in nodes must be named node_id')}
  
  if (!is.null(extended)){
    if(names(extended)[5]!='weight') {stop('5th column should be "weight"')}
    
    intra <- extended %>% filter(layer_from==layer_to)
    inter <- extended %>% filter(layer_from!=layer_to)
    # If there are no interlayer edges
    if (nrow(inter)==0){inter <- NULL}
    
    # Set the output formats
    if (intra_output_extended==F){intra %<>% select(layer=layer_from, node_from, node_to, weight)}
    # Write edge list to a txt file
    if (write_to_file==T){write_delim(extended, paste(filename_prefix,'.txt',sep=''), delim = ' ')}
  } else {
    if (!is.null(intra)){if(names(intra)[4]!='weight') {stop('4th column in Intralayer edge list should be weight')}}
    if (!is.null(inter)){if(names(inter)[5]!='weight') {stop('5th column in Interlayer edge list should be weight')}}
    
    # for non-extended output of intralayer edges
    if (intra_output_extended) {intra %<>% select(layer_from=layer, node_from, layer_to=layer, node_to, weight)}
    # Write edge list to a txt file
    write_delim(intra, paste(filename_prefix,'.txt',sep=''), delim = ' ')
    if (!is.null(inter)){write_delim(inter, paste(filename_prefix,'.txt',sep=''), delim = ' ', append = T)}
  }
  
  out <- list(intra=intra,
              inter=inter,
              nodes=nodes,
              layers=layers)
  class(out) <- 'multilayer'
  return(out)
}


check_infomap <- function(x='Infomap'){
  out <- attempt(system(paste('./',x,' -V',sep='')), msg = 'Infomap not installed correctly. See www.mapequation.org for instructions on how to install.')
  if (out==0) {
    return(T)
  } else {
    return(F)
  }
}


run_infomap_multilayer <- function(M,
                                   infomap_executable='Infomap',
                                   flow_model=NULL,
                                   silent=T,
                                   trials=100,
                                   seed=NULL,
                                   relax=F,
                                   multilayer_relax_rate=0.1,
                                   multilayer_relax_limit=NULL,
                                   multilayer_relax_limit_up=NULL,
                                   multilayer_relax_limit_down=NULL,
                                   temporal_network=F,
                                   run_standalone=T,
                                   remove_auxilary_files=T,
                                   ...){
  if(check_infomap(infomap_executable)==F){stop('Error in Infomap stand-alone file.')}
  if(class(M)!='multilayer'){stop('M must be of class multilayer')}
  
  # Infomap arguments
  arguments <- paste('--tree -2 -N ',trials, sep='')
  arguments <- ifelse(!is.null(seed), paste(arguments, '--seed',seed), arguments)
  arguments <- ifelse(!is.null(flow_model), paste(arguments, '-f',flow_model), arguments)
  arguments <- ifelse(silent, paste(arguments, '--silent'), arguments)
  arguments <- paste(arguments,...)
  
  # If using interlayer edges to determine flow
  if (relax==F){
    print('Using interlayer edge values to determine flow between layers.')
    # Write file for Infomap
    write_lines('*Multilayer', 'infomap_multilayer.txt')
    write_delim(M$intra, 'infomap_multilayer.txt', delim = ' ', append = T)
    write_delim(M$inter, 'infomap_multilayer.txt', delim = ' ', append = T)
  } else { # If using relax rates
    if (ncol(M$intra)==5){stop('Cannot use relax rates with extended format of intralayer edges. See function create_multilayer_object.')}
    print('Using global relax to determine flow between layers.')
    # Write file for Infomap
    write_lines('*Intra', 'infomap_multilayer.txt')
    write_delim(M$intra, 'infomap_multilayer.txt', delim = ' ', append = T)
    if(!is.null(M$inter)){
      if (ncol(M$inter)==5){stop('Cannot use relax rates with extended format of interlayer edges. See function create_multilayer_object.')}
      print('Global relax will be constrained by interlayer edges.')
      write_lines('*Inter', 'infomap_multilayer.txt', append = T)
      write_delim(M$inter, 'infomap_multilayer.txt', delim = ' ', append = T)
    }
    # Add arguments for relax rates and limits
    arguments <- ifelse(!is.null(multilayer_relax_rate), paste(arguments, '--multilayer-relax-rate',multilayer_relax_rate), arguments)
    arguments <- ifelse(!is.null(multilayer_relax_limit), paste(arguments, '--multilayer-relax-limit',multilayer_relax_limit), arguments)
    arguments <- ifelse(!is.null(multilayer_relax_limit_up), paste(arguments, '--multilayer-relax-limit-up',multilayer_relax_limit_up), arguments)
    arguments <- ifelse(!is.null(multilayer_relax_limit_down), paste(arguments, '--multilayer-relax-limit-down',multilayer_relax_limit_down), arguments)
  }
  
  # Run Infomap
  call <- paste('./',infomap_executable,' infomap_multilayer.txt . ', arguments, sep='')
  
  # If running within R
  if (run_standalone==T){
    print(call)
    system(call)
  } else {
    print('Please run Infomap online at https://www.mapequation.org/infomap/ using the following arguments (copy-paste):')
    print(arguments)
    invisible(readline(prompt="After running, download statenodes results and press [ENTER] when done"))
    if (!file.exists('network_states.tree')){stop('Result file network_states.tree was not found. Did you download results?')}
    file.rename(from = 'network_states.tree', to = 'infomap_multilayer_states.tree')
  }
  # Get L
  L_output <- parse_number(read_lines('infomap_multilayer_states.tree')[6])
  #Read infomap's output file
  modules <- suppressMessages(read_delim('infomap_multilayer_states.tree', delim = ' ', skip = 8, col_names = c('path', 'flow', 'name', 'state_id', 'node_id', 'layer_id')))
  # Parse modules
  modules %<>% 
    filter(flow > 0) %>%
    select(path, node_id, layer_id, flow) %>%
    separate(path, into = c("module", "leaf_id"), sep = ":") %>% 
    mutate_at(.vars = 1:4, as.integer) %>% 
    full_join(M$nodes, "node_id") %>% 
    select(node_id, starts_with("module"),  everything(), -leaf_id) %>% 
    arrange(node_id, layer_id)
  
  # For temporal networks, need to rename modules to be in a temporal order
  # because Infomap gives names by flow and not by order of appearence.
  if (temporal_network){
    print('Reorganizing modules...')
    renamed_moduels <- modules %>%
      distinct(module,layer_id) %>%
      arrange(module,layer_id)
    x <- c(1,table(renamed_moduels$module))
    module_birth_layers <- renamed_moduels %>% slice(cumsum(x)) %>% arrange(layer_id,module)
    module_renaming <- data.frame(module=module_birth_layers$module, module_renamed = 1:max(module_birth_layers$module))
    modules %<>%
      left_join(module_renaming, 'module') %>%
      select(-module) %>%
      rename(module=module_renamed)
  }
  
  if (remove_auxilary_files){
    print('Removing auxilary files...')
    file.remove('infomap_multilayer_states.tree')
    file.remove('infomap_multilayer.txt')
    file.remove('infomap_multilayer.tree')
  }
  # Output
  print(paste('Partitioned into ', max(modules$module),' modules.', sep=''))
  out <- list(call=call, L=L_output, m=max(modules$module), modules=modules)
  class(out) <- 'infomap_multilayer'
  return(out)
}