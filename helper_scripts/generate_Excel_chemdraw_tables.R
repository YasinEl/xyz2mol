library(data.table)

root = 'C:/PostDoc/Ming_time/example_files/csvs'

collision = 'CID3'

files = list.files(root, pattern = collision, full.names = TRUE, recursive = TRUE)
files

li_dt = list()
length(li_dt) = length(files)
max_rows = 0
for( i in seq(length(files))){
  dt = fread(files[i])
  max_rows = max(max_rows, nrow(dt))
}


for( i in seq(length(files))){
  number = sub('.*TMP\\.(\\d+)/.*', '\\1', files[i])
  
  
  dt = fread(files[i])
  colnames(dt) = paste0(rep(paste0(c(collision, '_'), collapse = ''), 2), '_', number)
  while(nrow(dt) < max_rows) {
    dt <- rbind(dt, data.table::data.table(matrix(NA, nrow = 1, ncol = ncol(dt), dimnames = list(NULL, names(dt)))))
  }
  li_dt[[i]]  = dt
  
  print(i)
  
}

dt_combined <- do.call(cbind, li_dt)


fwrite(dt_combined, paste0('C:/PostDoc/Ming_time/example_files/csvs_summary/summary_,', collision, '.csv'))
