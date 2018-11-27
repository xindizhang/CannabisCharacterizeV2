library(dplyr)
library(ggplot2)
library(readr) #for fast table reading
library(here)

#==========================Helper functions =====================================================
#read in expression data (exon array expression)
loadBrainSpanExpression <- function(pathToBrainSpanExpression, doRankTransform=TRUE) {
  rowLabels <- read.table(pathToBrainSpanExpression, header=T,row.names=1,sep=",",stringsAsFactors = F)$gene_symbol
  print("Loading expression, please wait")
  expression <- read_csv(ExpressionMatrix, col_names=F, progress=T)
  expression <- as.data.frame(expression)
  rownames(expression) <- expression[,1] #probably not needed
  expression[,1] <- NULL
  
  inboth <- c()
  backgroundGene <- read.table(here("CannabisCharacterizeV2", "geneLists", "Custom.CannabisGWAS.human.Table2.background.doi.org_10.1038_s41593-018-0206-1.txt"), stringsAsFactors = F)[,"V1"]
  
  intersection <- intersect(rowLabels, backgroundGene)
  lostbackGene <- setdiff(backgroundGene, intersection)
  for (rowLabel in intersection) {
    locations <- which(rowLabels == rowLabel)
    inboth <- c(inboth, locations)
  }
  expression <- expression[inboth,]
  rowLabels <- rowLabels[inboth]
  
  # Select duplicated genes
  dupRowLabels <- rowLabels[duplicated(rowLabels)]
  correlations <- c()
  dupedIndices <- c()
  if (doRankTransform) {
    # rank by columns
    expression <- apply(expression,2,rank) #use ranked expression
  }
  
  # Remove duplicated rows
  for(rowLabel in dupRowLabels) {
    # find the row number where the gene is duplicated
    locations <- which(rowLabels == rowLabel)
    #most are pefectly correlated except one gene (r = 0.69), with a correlation of 1
    # if it has a perfect correlation(r = 1), duplicated genes, can remove, if not .....loose
    # about data here
    correlations <- c(correlations, cor(as.numeric(expression[locations[1],]),as.numeric(expression[locations[2],])))
    dupedIndices <- c(dupedIndices, locations[2:length(locations)]) #some have more than one
  }
  rowsBefore <- nrow(expression)
  expression <- expression[-dupedIndices,]
  rowLabels <- rowLabels[-dupedIndices]
  rownames(expression) <- rowLabels #unqiue gene symbols as row names
  print(paste("Duplicate gene symbol rows removed:" ,rowsBefore-nrow(expression)))
  
  columnData <- getExpressionColumnData(ExpressionColumn)
  colnames(expression) <- rownames(columnData)
  expression
}

#get the meta data associated with the brainspan expression samples
getExpressionColumnData <- function(ExpressionColumn) {
  #now set the column names using donor and region name
  columnData <- read.table(ExpressionColumn, header=T,row.names=1,sep=",",stringsAsFactors = F)
  #update the name of a single region if we are using an old expression dataset (508 samples)
  sum(columnData$structure_name=="posterior (caudal) superior temporal cortex (area TAc)") #count of regions
  columnData[columnData$structure_name=="posterior (caudal) superior temporal cortex (area TAc)", "structure_name"] <- "posterior (caudal) superior temporal cortex (area 22c)"
  
  rownames(columnData) <- paste(columnData$donor_name, "region", columnData$structure_name,sep=".")
  #create age in Months column
  # negative acge means pwc
  convertAgeToMonths <- function(ageString) {
    ageDigit <- as.numeric(gsub(" .*", "", ageString))
    ageUnit <- gsub(".* ", "", ageString)
    if (ageUnit == "pcw") age <- (ageDigit-70)/12
    if (ageUnit == "yrs") age <- ageDigit*12
    if (ageUnit == "mos") age <- ageDigit
    age
  }
  columnData$AgeInMonths <- sapply(columnData$age, convertAgeToMonths)
  columnData
}

#convert age in months to ranked age values - gives ties the same rank
convertToRanksForAge <- function(ageVariable) {
  # rank returns the order of a object
  monthsToRankedAge <- c(rank(unique(ageVariable)))
  names(monthsToRankedAge) <- unique(ageVariable)
  ageVariableRanked <- monthsToRankedAge[as.character(ageVariable)]
  # two columns, age, rank
  ageVariableRanked
}

#=========================================== code =================================================
here()
#load BrainSpan expression
setwd(here("CannabisCharacterizeV2", "R"))

useRNAseq <- T

if (useRNAseq) {
  pathToBrainSpanExpression <- here("CannabisCharacterizeV2", "genes_matrix_csv", "rows_metadata.csv")
  ExpressionColumn <- here("CannabisCharacterizeV2", "genes_matrix_csv", "columns_metadata.csv")
  ExpressionMatrix <- here("CannabisCharacterizeV2", "genes_matrix_csv", "expression_matrix.csv")
} else {
  pathToBrainSpanExpression <- here("CannabisCharacterizeV2", "gene_array_matrix_csv") 
}
 
#load expression (slow), there are 16986 genes in total
expression <- loadBrainSpanExpression(pathToBrainSpanExpression, doRankTransform = FALSE)

genesOfInterest <- read.table(here("CannabisCharacterizeV2", "geneLists", "Custom.CannabisGWAS.human.Table2.doi.org_10.1038_s41593-018-0206-1.txt"),stringsAsFactors = F)[,"V1"]
length(genesOfInterest)


genesOfInterestFiltered <- intersect(genesOfInterest, rownames(expression))
lostGene <- setdiff(genesOfInterest, genesOfInterestFiltered)
print(paste("Lost Genes:", lostGene))

donorMetaData <- getExpressionColumnData(ExpressionColumn)


#filter regions to remove early non sampled regions
counts <- dplyr::summarise(group_by(donorMetaData, structure_name), n = n())
regionsWithFewSamples <- subset(counts, n < 7)$structure_name

samplesToPlot <- subset(donorMetaData, !(structure_name %in% regionsWithFewSamples))


samplesToPlot$AgeRank <- convertToRanksForAge(samplesToPlot$AgeInMonths)

expressionToPlot <- expression[,rownames(samplesToPlot)]
dim(expressionToPlot)

if (useRNAseq) {
  expressionToPlot <- log(expressionToPlot + 1)
  yaxisLabel <- "log(Expression) (z-score)"
} else {
  yaxisLabel <- "Expression (z-score)"
}
#scale across rows
#z-score across all regions
expressionToPlot <- t(scale(t(expressionToPlot)))

#compute wilcoxon for each column

genesOfNoInterest <- setdiff(rownames(expressionToPlot), genesOfInterest)
pValueColumns <- data.frame(pValue = numeric(), stringsAsFactors = F)
for(column in colnames(expressionToPlot)) { #speed it up somehow?
  pValue <- wilcox.test(expressionToPlot[genesOfInterest,column], expressionToPlot[genesOfNoInterest,column], alternative="greater")$p.value
  pValueColumns[column, "pValue"] <- pValue
}
pValueColumns$pValueAdjusted <- p.adjust(pValueColumns$pValue)

regionOrder <- rev(c("ventrolateral prefrontal cortex","dorsolateral prefrontal cortex","orbital frontal cortex","anterior (rostral) cingulate (medial prefrontal) cortex","posterior (caudal) superior temporal cortex (area 22c)","inferolateral temporal cortex (area TEv, area 20)","posteroventral (inferior) parietal cortex","primary motor cortex (area M1, area 4)","primary somatosensory cortex (area S1, areas 3,1,2)","primary auditory cortex (core)","primary visual cortex (striate cortex, area V1/17)","hippocampus (hippocampal formation)","amygdaloid complex","mediodorsal nucleus of thalamus","striatum","cerebellar cortex"))
regionToAcro <- unique(donorMetaData[,c("structure_name","structure_acronym")])
rownames(regionToAcro) <- regionToAcro$structure_name
regionOrderTable <- data.frame(structure_name=regionOrder, structure_acronym = regionToAcro[regionOrder, "structure_acronym"], stringsAsFactors = F)
textAcro <- "Brain regions include:"
for(row in rev(rownames(regionOrderTable))) {
  textAcro <- paste(textAcro, regionOrderTable[row, "structure_name"])
  textAcro <- paste0(textAcro, " (" , regionOrderTable[row, "structure_acronym"], "),")
  
}
textAcro

#add region and donor and age IDs
pValueColumns$age <- donorMetaData[rownames(pValueColumns),"age"]
pValueColumns$structure_name <- donorMetaData[rownames(pValueColumns),"structure_name"]
pValueColumns$AgeRank <- samplesToPlot[rownames(pValueColumns),"AgeRank"]
pValueColumns$donor_name <- donorMetaData[rownames(pValueColumns),"donor_name"]
pValueColumns$structure_acronym <- donorMetaData[rownames(pValueColumns),"structure_acronym"]


agePairs <- unique(pValueColumns[,c("AgeRank","age", "donor_name" )])
agePairs <- agePairs[order(agePairs$AgeRank,agePairs$donor_name),]

pValueColumns$pValueAdjustedBins <- cut(pValueColumns$pValueAdjusted,breaks=c(1,0.05,0.005,0.0005,0))

dplyr::summarise(pValueColumns,
                 regions = n_distinct(structure_name),
                 donors = n_distinct(donor_name),
                 observations = n_distinct(age)
)
nrow(pValueColumns)

regionOrderTable$structure_nameWithoutParam <- trimws(gsub("[(][^)]*[)]","", regionOrderTable$structure_name))

base_size <- 14
(heatmap <- ggplot(pValueColumns, aes(donor_name, structure_name)) +
    geom_tile(aes(fill = pValueAdjustedBins)) +
    scale_fill_manual(drop=FALSE, values=colorRampPalette(c("yellow","red", "black"))(4), na.value="#EEEEEE", name= "Corrected p-value", labels=c("< 0.0005", "< 0.005", "< 0.05", "> 0.05")) +
    theme_grey(base_size = base_size) + labs(x = "", y = "") +
    scale_x_discrete(expand = c(0, 0), limits=agePairs$donor_name, labels=agePairs$age) +
    scale_y_discrete(expand = c(0, 0), limits=regionOrderTable$structure_name, labels=regionOrderTable$structure_nameWithoutParam) + #switch to structure_acronym
    theme(axis.ticks.y = element_blank(), axis.text.x = element_text(size = 11, angle = 330, hjust = 0, colour = "black")) +
    theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
    theme(legend.position="bottom"))

# [END]