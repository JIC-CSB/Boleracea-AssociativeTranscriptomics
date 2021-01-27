require('data.table')

## Read in delimited file with sequence identifiers and trait values
traitDataFile <- tcltk::tk_choose.files(default = "", caption = "Please select one trait data file")
traitData <- read.delim(traitDataFile, header=TRUE, row.names=1, check.names = FALSE, na.strings = "-999") ## -999 is TASSEL output

## Read in delimited files with sequence names and rpkm values
rpkmC.temp <- read.table ("Bol_RPKMTASSEL.txt", header=TRUE, sep= "\t", row.names=1) # can use the check.names=FALSE argument if there's "unusual" symbols in the column headers  2018-01-15

outputName=""
outputName=readline(prompt="Please enter output identifier...   ")
print("reading files...")

traitLabel = colnames(traitData)

meansCutOff = 0.4

## remove low rpkm means
print(paste0("Removing markers with a mean RPKM < ", meansCutOff))
delrpkmC = rpkmC.temp[rowMeans(rpkmC.temp) >= meansCutOff, ]

rpkmC <- t(delrpkmC)

## Delete missing varieties from rpkm files
rpkmmergeC <- merge(traitData,rpkmC, by="row.names")
rownames(rpkmmergeC) <- rpkmmergeC[,1]

## Function for apply()ing. Does a linear model between trait data and genes
performLinearRegression <- function(exp_vector) {
  model <- lm(rpkmmergeC[,2] ~ exp_vector)
  return(as.numeric(anova(model)[1,]))
}

## Do the linear models. Each gene is in a column hence apply() over
## columns. rpkmmergeA contains lines in column 1 and the trait data
## in column 2 hence these are dropped in the apply() call. The beauty
## of this way is that it makes the row names the genes in the final
## lmResults output.
lmResults <- t(apply(rpkmmergeC[, -c(1,2)], 2, performLinearRegression))
colnames(lmResults) <- c("Df", "SumSq", "MeanSq", "Fvalue", "Pvalue")

## calculate log10P
log10P <- -log10(lmResults[,"Pvalue"])
print("Calculated p-values")

## Read B. oleracea to Arabidopsis matches 
codesFile = read.csv("Marker_to_At_Bol.csv", header=TRUE, na.strings = ".", fileEncoding="UTF-8-BOM")

## Get everything together for writing to a file. Made an explicit
## rownames column so don't have to use row.names in next regress
## plotter step. Not a problem because I use row.names=FALSE in
## write.table. Note needed to use as.data.frame otherwise it all came
## out as a list of character vectors i.e. in quotes which obviously
## can't be order()ed.
getGeneModels = function(results) {
  unigenes = unlist(strsplit(sub("_", "~#~", results$unigene), "~#~"))
  ArabidopsisHits = codesFile[ codesFile$unigene %in% unigenes, ]
  results$unigenes = unigenes
  r = data.table( results )
  ah = data.table( ArabidopsisHits )
  resultsWithAGI = merge(r, ah, by="unigene")  ## likely in the future can just merge by="unigene" as identical(as.character(results$unigene), results$tmp) == TRUE
  return(resultsWithAGI)
}

BoleraceaUnigene = rownames(lmResults)
traitAndUnigene = as.data.frame(cbind(trait = traitLabel, unigene = BoleraceaUnigene))

print("Producing final results")

finalResults <- cbind(traitAndUnigene, log10P, lmResults)
finalResults = getGeneModels(finalResults)
finalResults <- finalResults[order(finalResults$log10P, decreasing = TRUE),]

## Write these columns for now. In future if the apply() is quick
## enough probably combine this script and grapher one together.
finalResults = finalResults[,c("trait", "unigene", "C.Chr", "AGI", "log10P", "Df", "SumSq", "MeanSq", "Fvalue", "Pvalue")]
resultsFile = paste("Reg-", paste(outputName, collapse="-"), ".txt", sep="")
write.table(finalResults, resultsFile, quote=FALSE, sep="\t" ,row.names=FALSE, col.names=TRUE)



print ("Phew! Complete...")







