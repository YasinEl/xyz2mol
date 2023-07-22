library(data.table)


fill_empty_rows = function(li){
  
  max_rows = max(unlist(lapply(li, nrow)))
  
  op_li = 
  lapply(li, function(x){
    while(nrow(x) < max_rows) {
      x <- rbind(x, data.table::data.table(matrix(NA, nrow = 1, ncol = ncol(x), dimnames = list(NULL, names(x)))))
    }
    return(x)
  })
  
  return(op_li)
}


remove_equilibrium_reactions = function(dt){
  
  #remove the N2
  dt[, SMILES := sapply(SMILES, function(x) {
    # Split each string in SMILES by '.'
    split_strs = unlist(strsplit(x, split = "\\.", fixed = FALSE))
    
    # Remove 'N#N' from the list
    filtered_strs = split_strs[split_strs != "N#N"]
    
    # Collaps back the remaining strings by '.' and return
    return(paste(filtered_strs, collapse = "."))
  })]
  
  dt[, idx := .I]
  
  # Create an intermediate index
  dt[, indx := .GRP, by = SMILES]
  
  dt = dt[, .(idx_min = min(idx),
              idx_max = max(idx),
              collisions_min = min(collisions),
              collisions_max = max(collisions)), by =.(SMILES)]
  
  dt = dt[order(idx_min)]
  dt[, id := .I]
  
  dt[, flag := FALSE]
  smiles_entry = 1
  while(smiles_entry < nrow(dt)){
    dt[idx_max < dt[id == smiles_entry, idx_max] &
         id > smiles_entry, flag := TRUE]
    dt = dt[flag == FALSE]
    dt[, id := .I]
    smiles_entry = smiles_entry + 1
  }
  dt = dt[, .(collisions = paste0(c(collisions_min, collisions_max), collapse = '_')), by =.(id, SMILES)]

  return(dt)
  
}




root = 'C:/PostDoc/Ming_time/example_files/csvs'

folders = list.files(root, full.names = TRUE, pattern = '2_protonated_mol_3')
reduce_to_relevant = TRUE

li_allTrj = list()
length(li_allTrj) = length(folders)

#handle all folders/trajectories
for (folder_idx in seq(length(folders))){
  
  cid_files = list.files(folders[folder_idx], full.names = TRUE, pattern = 'CID')
  MDtrj_files = list.files(folders[folder_idx], full.names = TRUE, pattern = 'MDtrj')
  
  li_singleTrj = list()
  length(li_singleTrj) = length(cid_files) + length(MDtrj_files)
  li_idx = 1
  #simulation startup
  dt = fread(MDtrj_files[1])
  dt[, file := basename(MDtrj_files[1])]
  dt[, collisions := 0]
  li_singleTrj[[li_idx]] = dt
  li_idx = 2
  #handle 
  for (cid_idx in seq(length(cid_files))){
    for (MD_type in c('CID', 'MDtrj')){
      if(MD_type == 'CID'){
        current_file = cid_files[cid_idx]
      } else if (MD_type == 'MDtrj'){
        current_file = MDtrj_files[cid_idx + 1]
      }
      
      if(!is.na(current_file)){
        dt = fread(current_file)
        dt[, file := basename(current_file)]
        dt[, collisions := cid_idx]
        li_singleTrj[[li_idx]] = dt
        li_idx = li_idx + 1
      }
    }
  }
  
  dt_trj = rbindlist(li_singleTrj)
  number = sub('.*TMP\\.(\\d+)$', '\\1', folders[folder_idx])
  if (reduce_to_relevant == TRUE){
    setnames(dt_trj, 'V2', 'SMILES')
    dt_trj = remove_equilibrium_reactions(dt_trj)
    colnames(dt_trj) = paste0(c('id', 'SMILES', 'collisions'), '_', number)
  } else if (reduce_to_relevant == FALSE){
    colnames(dt_trj) = paste0(c('id', 'SMILES', 'file', 'collisions'), '_', number)
  }

  
  li_allTrj[[folder_idx]] = dt_trj
  
  print(paste0(folder_idx, '/', length(folders)))
} 
li_allTrj = fill_empty_rows(li_allTrj)

dt_combined <- do.call(cbind, li_allTrj)

fwrite(dt_combined, paste0('C:/PostDoc/Ming_time/example_files/csvs_summary/summary_,', '2_protonated_mol_3_reduced', '.csv'))

