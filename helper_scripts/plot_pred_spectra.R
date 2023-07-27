library(data.table)
library(ggplot2)
library(ChemmineOB)
library(ChemmineR)
library(rcdk)        
library(cowplot)


plot_molecule <- function(molecule, name = NULL, sma = NULL, ...){
  #' molecule an object as returned by rcdk::load.molecules or rcdk::parse.smiles()
  #' name a character for the name of the molecule, 
  #' sma a character witht the smarts string as passed onto get.depictor()
  #' ... other arguments for get.depictor()
  
  tryCatch({
    # Image aesthetics 
    dep <- get.depictor(
      width = 1000, height = 1000,
      zoom = 7, sma = sma, ...
    )
    molecule_sdf <- view.image.2d(molecule[[1]], depictor = dep)
    
    ## Remove extra margins around the molecule
    par(mar=c(0,0,0,0))
    plot(NA, 
         xlim=c(1, 10), ylim=c(1, 10), 
         # Remove the black bounding boxes around the molecule
         axes = F)
    rasterImage(molecule_sdf, 1,1, 10,10)
    # Annotate the molecule
    text(x = 5.5, y = 1.1,  deparse(substitute(molecule)))
  },
  error=function(cond) {
    message("An error occurred: ", cond)
    plot(NA, xlim=c(1, 10), ylim=c(1, 10), axes = F)  # Empty plot
  },
  finally={
    message("Finished plotting.")
  })
}


fl = list.files('C:/PostDoc/Ming_time/example_files/tim_mol1', full.names = TRUE, pattern = '\\.csv')
mols = fread('C:/PostDoc/Ming_time/example_files/tim_mol1/mols/mols.csv')
li = list()
length(li) = length(fl)

li_smiles = list()
length(li_smiles) = length(fl)

for(i in seq(length(li))){
  
  dt = fread(fl[i])
  colnames(dt) = c('mz', 'int')
  dt[, smiles := mols[i]$smiles]
  li[[i]] = dt
  
  #create smiles depiction
  png(filename="p.png")
  smiles_mol <- parse.smiles(mols[i]$smiles)
  plot_molecule(molecule = smiles_mol ,
                #abbr = "reagents", 
                #annotate = "number",
                suppressh = F
  )
  dev.off()
  plot_mol <- ggdraw() + draw_image("p.png")
  li_smiles[[i]] = plot_mol
  
  
}


dt = rbindlist(li)






p = ggplot(dt[int > 0.001], aes(x = mz, ymin = 0, ymax = int)) +
  geom_linerange()  +
  geom_text(data = dt[int > 0.2], aes(x = mz, y = int, label = round(mz, 3)), 
            vjust = -0.5, size = 3) +
  theme_classic() +
  facet_wrap(~smiles, ncol = 1)


p_mols = eval(parse(text = paste0('plot_grid(', paste0(paste0('li_smiles[[', seq(length(li_smiles)), ']]'), collapse = ', '), ', ncol = 1)')))


final_plot <- plot_grid(p_mols, p, ncol = 2, rel_widths = c(1, 3))


final_plot



























