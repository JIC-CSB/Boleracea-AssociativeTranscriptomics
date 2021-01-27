## This script creates Manhattan Plots for TASSEL output

rm(list = ls())

require('tidyr')
require('dplyr')
require('stringr')

VERSION = "1.0.0"

MODEL.MLM = 'MLM'
MODEL.GLM = 'GLM'

require('rstudioapi')
  
  script_path = rstudioapi::getActiveDocumentContext()$path
  setwd(dirname(script_path))

getDataFiles <- function( model_type ) {
    tassel_stats = tcltk::tk_choose.files(caption = paste("Select the `stats` file from", model_type, 'analysis'), multi = FALSE)
    tassel_effects = tcltk::tk_choose.files(caption = paste("Select the `effects` file from", model_type, 'analysis'), multi = FALSE)
    annotation_file = tcltk::tk_choose.files(caption = "Select Arabidopsis annotation file", multi = FALSE)
    direction_file = tcltk::tk_choose.files(caption = "Select graph direction file", multi = FALSE)
    
    tassel_stats = paste( tassel_stats, collapse = " " )
    tassel_effects = paste( tassel_effects, collapse = " " )
    annotation_file = paste( annotation_file, collapse = " " )
    direction_file = paste( direction_file, collapse = " " )
    
    l = list( "tassel_stats" = tassel_stats,
              "tassel_effects" = tassel_effects,
              "annotation_file" = annotation_file,
              "direction_file" = direction_file )
    
    return(l)
  }


getGrapherDirections <- function( direction_file ) {
  ## Read GrapherDirections file with graphing orders
  graph_df = read.delim( direction_file, stringsAsFactors = FALSE )
  
    graph_df = graph_df[, c('Unigene.Marker', 'Chr', 'sort.SNP')]
    colnames(graph_df)[ colnames(graph_df) == 'Unigene.Marker' ] <- 'Marker'
    colnames(graph_df)[ colnames(graph_df) == 'sort.SNP' ] <- 'sort'

  graph_df = graph_df %>% drop_na(sort)
  graph_df = subset( graph_df, Chr == "C1" | Chr == "C2" | Chr == "C3" | Chr == "C4" | Chr == "C5" | Chr == "C6" | Chr == "C7" | Chr == "C8" | Chr == "C9" )

  
  return (graph_df)
}

load_data <- function( model_type, tassel_stats, tassel_effects ) {
  ## Load the effects data - size of effect for each trait/marker/allele
  ##
  print( 'Loading the effects file...' )
  
  effects_df = read.delim( file = tassel_effects )
  if ( model_type == MODEL.GLM ) {
    # use the MLM column names - adapt the GLM column differences to the MLM format
    colnames(effects_df)[ colnames(effects_df) == 'Estimate' ] <- 'Effect'
  }
  effects_df = effects_df[, c( 'Trait', 'Marker', 'Allele', 'Effect', 'Obs' ) ]
  
  ## Load the stats data from the GLM/MLM analysis
  ## 
  print( 'Loading the statistics file...' )
  
  stats_df = read.delim( file = tassel_stats )
  stats_df = stats_df[, c( 'Trait', 'Marker', 'p' ) ]
  colnames(stats_df)[ colnames(stats_df) == 'p' ] <- 'P.value'
  stats_df = stats_df %>% drop_na()
  
  ## Merge and post-process...
  print( 'Merging data...' )
  
  result_df = merge( effects_df, stats_df, by = c( 'Trait', 'Marker' ) )
  
  # set NonEffAllele -> the Allele at the end of Marker
  result_df$NonEffAllele = str_sub( result_df$Marker,-1,-1)
  
  # cat the trait/marker to make a unique marker
  result_df$UniqueMarker = paste0( result_df$Marker, ',', result_df$Allele )
  
  return(result_df)
}

getTRAIT <- function ( result_df, trait_select ) {
  trait_df = result_df[  result_df$Trait == trait_select , ]
}

getTasselData <- function( trait_df ) {
  ## Get the output file from TASSEL and augment with "-log10(p)"

  tassel_df = trait_df[ !is.nan(trait_df$P.value), ]
  
  tassel_df$log10P = -log10(tassel_df$P.value)
  
  tassel_df = tassel_df %>% rename( Allele1 = Allele,
                                    Allele2 = NonEffAllele)
  
  tassel_df = separate(
    tassel_df,
    col = "Marker",
    into = c("Marker"),
    sep = ",",
    extra = 'drop'
  )
  
  tassel_df = tassel_df[, c("Marker", "Allele1", "Allele2", "log10P")]
  
  return(tassel_df)
}

getChromoInfo <- function( graph_df ) {
  min_vals = setNames( aggregate( sort ~ Chr, data = graph_df, FUN = min ), nm = c( 'Chr', 'Start' ) )
  max_vals = setNames( aggregate( sort ~ Chr, data = graph_df, FUN = max ), nm = c( 'Chr', 'End' ) )
  
  chromo_df = merge( min_vals, max_vals, by = 'Chr' )
  chromo_df = separate( data = chromo_df,
                        col = 'Chr',
                        into = c( 'Genome', 'ChromoNum' ),
                        sep = c( 1 ),
                        remove = FALSE, 
                        convert = TRUE )
  
  chromo_df = chromo_df[ order(chromo_df$Genome, chromo_df$ChromoNum), ]
  
  return(chromo_df)
}

getAnnotationData <- function( annotation_file, mode ) {
  ## Read in annotations file
  ## Format:  Unigene         Unigene.Marker          Chr     AGI
  ##          BnaA01g00650D   BnaA01g00650D:144:A     A1      AT4G37280.1
  annos_df = read.delim( file = annotation_file )
  
  annos_df = annos_df[, c('Unigene.Marker', 'Chr', 'AGI')]
  colnames(annos_df)[ colnames(annos_df) == 'Unigene.Marker' ] <- 'Marker'

    annos_df = subset( annos_df, Chr == "C1" | Chr == "C2" | Chr == "C3" | Chr == "C4" | Chr == "C5" | Chr == "C6" | Chr == "C7" | Chr == "C8" | Chr == "C9" )
  
  return(annos_df)
}

createGenomeManhattan <- function(filename, df, ylim) {
  ## Saves output as jpg
  ##
  jpeg(
    filename = filename,
    quality = 2000,
    width = 1500,
    height = 800,
    units = "px",
    bg = "white",
    res = NA
  )
  
  ## Plot graphs
  
  
  genomes = c("C")
  
  for (genome in  genomes) {
    genome_df = df[df$Genome == genome, ]
    xlim = c(1, max(genome_df$End))
    
    plot (
      genome_df$sort,
      genome_df$log10P,
      type = "p",
      pch = 19,
      xlim = xlim,
      ylim = ylim,
      main = genome,
      cex.main = 2,
      ylab = "-log10(p)",
      xlab = " ",
      cex.lab = 1.5,
      col = ifelse(genome_df$ChromoNum %% 2 == 0, "red", "black"),
      xaxt = 'n'
    )
  }
  
  close.screen(all = TRUE)
  dev.off()
  
  print( paste0( "Saved Manhattan plot to ", filename ) )
  
  return()
}

createChromosomeManhattan <- function(df, ylim, choices, annos_df) {
  ## plot chosen graph and allow point identification
  print( 'Generating Manhattan plot for chromosomes. Please choose one from the list.')
  
  the_terminator = "No thanks"
  choices = c(choices, the_terminator)
  
  choice = "Yes please!"
  while (choice != the_terminator) {
    choice = select.list(choices = choices,
                         title = "Chromosome",
                         graphics = TRUE)
    
    if (choice == the_terminator) {            
      close.screen(all = TRUE)
      dev.off()
      
      break
    }
    
    chrom_df = df[df$Chr == choice,]
    
    x = chrom_df$sort
    y = chrom_df$log10P
    labels = chrom_df$Marker
    
    # 'Start' and 'End' are all the same
    xlim = c(chrom_df$Start[1], chrom_df$End[1])
    
    plot(
      x = x,
      y = y,
      type = "p",
      pch = 20,
      xlim = xlim,
      ylim = ylim,
      main = choice,
      ylab = "-log10(p)",
      xlab = "",
      col = "black",
      xaxt = 'n'
    )
    
    print( "Select points on plot - press `escape' when finished")
    
    idx = identify(
      x = x,
      y = y,
      labels = labels,
      pos = TRUE,
      n = 10,
      plot = TRUE,
      atpen = TRUE,
      tolerance = 0.1,
      col = "red",
      cex = 0.6
    )
    
    if (select.list( choices = c("Yes", "No"), title = "Save as jpeg?", graphics = TRUE ) == "Yes") {
      filename = paste0 (outputName, "_", "Edited_", choice, ".jpeg")
      dev.copy(jpeg, filename)
      dev.off()
      
      print( paste("Plot saved to", filename ))
    }
    
    selectedMarkers = chrom_df[idx$ind,]
    results = getHitsAndEffects( selectedMarkers = selectedMarkers, 
                                 annos_df = annos_df )
    print(results)
    
    csv_file = paste0 (outputName, "_", "Edited_", choice, ".tsv")
    write.csv( x = results, file = csv_file, row.names = FALSE )
    
    print( paste( "Information on selected points saved to", csv_file ) )
  }
}

getHitsAndEffects <- function(selectedMarkers,
                              annos_df) {
  ## Find markers matching Arabidopsis
  ath_hits = annos_df[ annos_df$Marker %in% selectedMarkers$Marker, c('Marker', 'AGI')]
  
  ## Link marker <-> Ath gene
  result_df = merge( ath_hits,
                     selectedMarkers,
                     by = "Marker")
  
  ## Select only important columns
  result_df = result_df[, c(
    "Marker",
    'Chr',
    "AGI",
    "Allele1",
    "Allele2",
    "log10P"
  )]
  
  result_df = result_df[ order(result_df$log10P, decreasing = TRUE), ]
  
  return(result_df)
}

####################################################################################################
####################################################################################################



model_type = select.list( choices = c(MODEL.GLM, MODEL.MLM), title = 'Please select the model type...', graphics = FALSE)
outputName = readline(prompt = "Please enter output identifier...   ") 

filenames = getDataFiles( model_type = model_type )

print( "Processing..." )

graph_df = getGrapherDirections( direction_file = filenames$direction_file )
result_df = load_data( model_type = model_type, tassel_effects = filenames$tassel_effects, tassel_stats = filenames$tassel_stats)

trait_select = readline(prompt = "Please enter trait you wish to graph... ")

trait_df = getTRAIT( result_df = result_df, trait_select = trait_select )
tassel_df = getTasselData( trait_df = trait_df )
chromo_df = getChromoInfo( graph_df = graph_df )
annos_df = getAnnotationData( annotation_file = filenames$annotation_file )

## merge and sort...
tassel_df = merge(tassel_df, graph_df, by.x = "Marker", na.rm = TRUE)
tassel_df = merge(tassel_df, chromo_df, by.x = "Chr")
tassel_df = tassel_df[with(tassel_df, order(sort)),]

ylim = c(0, max(tassel_df$log10P))

createGenomeManhattan(
  filename = paste0 (outputName, "Hemis_ManhattanStyle_", ".pdf"),
  df = tassel_df,
  ylim = ylim
)

createChromosomeManhattan(
  df = tassel_df,
  ylim = ylim,
  choices = chromo_df$Chr,
  annos_df = annos_df
)

print( paste("Thank you for using TASSEL Grapher version", VERSION) )
