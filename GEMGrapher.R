rm(list = ls())

require('tidyr')

TEST = FALSE

VERSION = "1.0.0"

require('rstudioapi')

script_path = rstudioapi::getActiveDocumentContext()$path
setwd(dirname(script_path))

getDataFiles <- function( mode ) {
  metal_file = tcltk::tk_choose.files(caption = paste("Select", "output file from Regress"))
  annotation_file = tcltk::tk_choose.files(caption = "Select Arabidopsis annotation file")
  direction_file = tcltk::tk_choose.files(caption = "Select C graph direction file")
 
  l = list( "metal_file" = metal_file,
            "annotation_file" = annotation_file,
            "direction_file" = direction_file )
  
  return(l)
}

getGrapherDirections <- function( direction_file, mode ) {
  ## Read GrapherDirections file with graphing orders
  graph_df = read.delim( direction_file, stringsAsFactors = FALSE )
  
    graph_df = graph_df[, c('Unigene', 'Chr', 'sort.GEM')]
    colnames(graph_df)[ colnames(graph_df) == 'Unigene' ] <- 'Marker'
    colnames(graph_df)[ colnames(graph_df) == 'sort.GEM' ] <- 'sort'
    graph_df = graph_df[ !duplicated(graph_df), ]

  graph_df = graph_df %>% drop_na(sort)
  
  return (graph_df)
}

getMetalData <- function( metal_file ) {
  ## Get the output file from METAL and augment with "-log10(p)"
  
  metal_df = read.delim( metal_file )
  #metal_df = metal_df[ !is.nan(metal_df$P.value), ]
  
  ##rename file columns
  names(metal_df)[names(metal_df) == "unigene"] <- "Unigene"
  names(metal_df)[names(metal_df) == "C.Chr"] <- "Chr"
  names(metal_df)[names(metal_df) == "Pvalue"] <- "P.value"
  
  
  
  #metal_df$log10P = -log10(metal_df$P.value)
  
  metal_df = separate(
    metal_df,
    col = "Unigene",
    into = c("Marker"),
    sep = ",",
    extra = 'drop'
  )
  
  metal_df = metal_df[, c("Marker", "log10P")] #Allele1 Allele2
  
  return(metal_df)
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

    annos_df = annos_df[, c('Unigene', 'Chr', 'AGI')]
    colnames(annos_df)[ colnames(annos_df) == 'Unigene' ] <- 'Marker'
    annos_df = annos_df[ !duplicated(annos_df), ]
  
  return(annos_df)
}

createGenomeManhattan <- function(filename, df, ylim) {
  ## Saves output as jpg
  ##
  jpeg(
    filename = filename,
    quality = 1500,
    width = 2000,
    height = 800,
    units = "px",
    bg = "white",
    res = NA
  )
  
  ## Sets up a split screen
  #par(mfrow = c(2, 1))
  
  ## Plot graphs
  
  for (genome in c("C")) {
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
      pch = 19,
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
      filename = paste0 (outputName, "_Edited_", choice, ".jpeg")
      dev.copy(jpeg, filename)
      dev.off()
      
      print( paste("Plot saved to", filename ))
    }
    
    selectedMarkers = chrom_df[idx$ind,]
    results = getHitsAndEffects( selectedMarkers = selectedMarkers,
                                 annos_df = annos_df )
    print(results)
    
    csv_file = paste0 (outputName, "_Edited_", choice, ".tsv")
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
    #"Allele1",
    #"Allele2",
    "log10P"
  )]
  
  result_df = result_df[ order(result_df$log10P, decreasing = TRUE), ]
  
  return(result_df)
}

####################################################################################################
####################################################################################################

outputName = readline(prompt = "Please enter output identifier...   ")

filenames = getDataFiles()

print( "Processing..." )

metal_df = getMetalData( metal_file = filenames$metal_file )
graph_df = getGrapherDirections( direction_file = filenames$direction_file, mode = mode )
chromo_df = getChromoInfo( graph_df = graph_df )
annos_df = getAnnotationData( annotation_file = filenames$annotation_file, mode = mode )

## merge and sort...
metal_df = merge(metal_df, graph_df, by.x = "Marker", na.rm = TRUE, all.y=T)
metal_df = merge(metal_df, chromo_df, by.x = "Chr")
metal_df = metal_df[with(metal_df, order(sort)),]

ylim = c(0, 0.5+max(metal_df$log10P[!(is.na(metal_df$log10P))]))

TMP=sort(metal_df$log10P)
createGenomeManhattan(
  filename = paste0 (outputName, "Hemis_ManhattanStyle_", ".jpeg"),
  df = metal_df,
  ylim = ylim
)

createChromosomeManhattan(
  df = metal_df,
  ylim = ylim,
  choices = chromo_df$Chr,
  annos_df = annos_df
)

print( paste("Thank you for using Grapher version", VERSION) )
