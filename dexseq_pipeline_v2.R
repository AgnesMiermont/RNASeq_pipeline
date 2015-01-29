library(DEXSeq)
library(BiocParallel)


getArgs <- function() {
  myargs.list <- strsplit(grep("=",gsub("--","",commandArgs()),value=TRUE),"=")
  myargs <- lapply(myargs.list,function(x) x[2] )
  names(myargs) <- lapply(myargs.list,function(x) x[1])
  return (myargs)
}



########################## read arguments
dexseq.compute <- TRUE

BPPARAM = MulticoreParam(workers=8)

annotation.file <- '/cluster/project8/vyp/vincent/Software/RNASeq_pipeline/bundle/Tc1_mouse/tc1_annotations.tab'
iFolder <- '/scratch2/vyp-scratch2/IoN_RNASeq/Frances/processed'
support.frame <- 'data/RNASeq_AD_Tc1J20.tab'
code <- 'Zanda_AD_Tc1J20'
gff <- '/cluster/project8/vyp/vincent/Software/pipeline/RNASeq/bundle/Tc1_mouse/GTF/Tc1.gff'
keep.dups <- FALSE

myArgs <- getArgs()
if ('support.frame' %in% names(myArgs)) support.frame <- myArgs[['support.frame']]
if ('code' %in% names(myArgs)) code <- myArgs[['code']]
if ('iFolder' %in% names(myArgs)) iFolder <- myArgs[['iFolder']]
if ('annotation.file' %in% names(myArgs)) annotation.file <- myArgs[['annotation.file']]
if ('gff.file' %in% names(myArgs)) gff <- myArgs[['gff.file']]
if ('keep.dups' %in% names(myArgs)) keep.dups <- as.logical(myArgs[['keep.dups']])





###check input files and data frame
message('Now reading ', support.frame)
support <- read.table(support.frame, header = TRUE, stringsAsFactors = FALSE)
my.ids <- support$sample
list.conditions <- grep(names(support), pattern = '^condition.*', value  = TRUE)
annotation <- read.table(annotation.file, header = TRUE, sep = '\t', na.string = c('NA', ''))

files <- paste(iFolder, '/', my.ids, '/dexseq/', my.ids, '_dexseq_counts.txt', sep = '')
if (sum(!file.exists(files)) > 0) stop('Some input files are missing')




### dexseq output folders
dexseq.folder <- paste(iFolder, '/dexseq', sep = '')
dexseq.counts <- paste(dexseq.folder, '/dexseq_counts_', code, '.RData', sep = '')  ##this contains the key data
if (!file.exists(dexseq.folder)) dir.create(dexseq.folder)




my.ids <- support$sample
if (!keep.dups) countFiles <- paste(iFolder, '/', my.ids, '/dexseq/', my.ids, '_dexseq_counts.txt', sep = '')
if (keep.dups) countFiles <- paste(iFolder, '/', my.ids, '/dexseq/', my.ids, '_dexseq_counts_keep_dups.txt', sep = '')


for (condition in list.conditions) {

  message('Condition ', condition)
  support.loc <- support

  ##handle the type variable
  support.loc$condition <- factor(support[, condition])
  loc.countFiles <- countFiles[ !is.na(support.loc$condition) ]
  support.loc <-  support.loc[ !is.na(support.loc$condition), ]

  
  ##handle the type variable
  type.loc <- gsub(x = condition, pattern = 'condition', replacement = 'type')
  if ( (! type.loc %in% names(support.loc)) & ('type' %in% names(support.loc))) {type.loc <- 'type'}  ##if only "type" is present, use it
  if (type.loc %in% names(support)) {
    support.loc$type <- factor(support.loc[, type.loc])
    support.loc <- subset(support.loc, !is.na(type))
  }
  
  loc.code <-  paste(unique(support.loc$condition), collapse = '_')
  message('Support data frame', loc.code)
  print(support.loc)

  
  ################### create the appropriate folders
  loc.dexseq.folder <- paste(iFolder, '/dexseq/', loc.code, sep = '')
  dexseq.figs <- paste(loc.dexseq.folder, '/figs', sep = '')
  dexseq.data <- paste(loc.dexseq.folder, '/dexseq_', code, '_', loc.code, '.RData', sep = '')  ##this will contain the final output of dexseq
  
  for (folder in c(loc.dexseq.folder, dexseq.figs)) {
    if (! file.exists(folder)) dir.create(folder)
  }

  
  if (dexseq.compute) {  
    #load(dexseq.counts)  ##object mycounts is key
    #DexSeqExons <- subset(DexSeqExons, c(rep(TRUE, 300), rep(FALSE, nrow(counts(DexSeqExons)) - 300))) 
    

    use.covariate <- FALSE
    if ('type' %in% names(support.loc)) {
      if (length(unique(as.character(support.loc$type))) > 1) {
        use.covariate <- TRUE
      }
    }

    if (use.covariate) {
      formuladispersion <- ~ sample + (condition + type) * exon
      formula0 <-  ~ sample + type * exon + condition
      formula1 <-  ~ sample + type * exon + condition * exon
      my.design <- support.loc[, c('type', 'condition')]
      my.design$type <- factor(my.design$type) ## probably not needed
      my.design$condition <- factor(my.design$condition)  ## probably not needed
      my.design.loc <- my.design  ##just to print basically
    } else {
      formuladispersion <-  ~ sample + condition * exon
      formula0 <-  ~ sample + exon + condition
      formula1 <-  ~ sample + exon + condition * exon
      my.design <- factor(support.loc[, c('condition')])
      my.design.loc <- support.loc[, c('condition'), drop = FALSE]  ##just to print basically
    }

    row.names(my.design.loc) <- factor(support.loc$sample)
    
    DexSeqExons.loc <- DEXSeqDataSetFromHTSeq(loc.countFiles,
                                              sampleData= my.design.loc,
                                              design=  formula1,
                                              flattenedfile= gff)
    
    write.table(x = my.design.loc, file = paste(loc.dexseq.folder, '/design.tab', sep = ''), row.names = TRUE, quote = FALSE, sep = '\t')

    
    
    #message('Updating the dexseq object')
    #DexSeqExons.loc <- DEXSeqDataSet(countData= featureCounts(mycounts),
    #                                 sampleData = my.design.loc,
    #                                 design= formula1,
    #                                 groupID=geneIDs(mycounts),
    #                                 featureID=exonIDs(mycounts),
    #                                 transcripts = DexSeqExons@featureData@data$transcripts,
    #                                 featureRanges = GRanges(DexSeqExons@featureData@data$chr, IRanges (start = DexSeqExons@featureData@data$start, end = DexSeqExons@featureData@data$end)) )


    ##DexSeqExons.loc <- DexSeqExons.loc[1:300,]  ##VP
    
    message('Starting the computations')
    DexSeqExons.loc <- estimateSizeFactors(DexSeqExons.loc)

    
    message('Here is the part that takes a lot of time')
    DexSeqExons.loc <- DEXSeq::estimateDispersions(DexSeqExons.loc, BPPARAM=BPPARAM)
    message('Done with estimateDispersions')    
    DexSeqExons.loc <- DEXSeq::testForDEU(DexSeqExons.loc, BPPARAM=BPPARAM)
    message('Done with testDEU')    
    DexSeqExons.loc <- DEXSeq::estimateExonFoldChanges(DexSeqExons.loc, BPPARAM=BPPARAM)
    message('Done with estimateFoldChange')
    
######################### output basic table
    res <- DEXSeq::DEXSeqResults (DexSeqExons.loc)
    logname <- grep(names(res), pattern = 'log2fold', value = TRUE)
    res.clean <- as(res[, c('groupID', 'featureID', 'exonBaseMean', logname, 'dispersion', 'stat', 'pvalue')], 'data.frame')
    names(res.clean)<- c("EnsemblID", "exonID", "meanBase", "log2FoldChange", "dispersion", "stat", "pvalue")

    res.clean$FDR <- p.adjust(res.clean$pvalue, method = 'fdr')    
    res.clean$chromosome <- as.character(seqnames( res$genomicData))
    res.clean$exon.start <- start(res$genomicData)
    res.clean$exon.end <- end(res$genomicData)


    res.clean$external_gene_id <- annotation$external_gene_id[ match(res.clean$EnsemblID, table = annotation$EnsemblID) ]
    if ('strand' %in% names(annotation)) res.clean$strand <- annotation$strand[ match(res.clean$EnsemblID, table = annotation$EnsemblID) ] ## add strand if available
    #res.clean <- merge(res.clean, annotation[, c('EnsemblID', 'external_gene_id', 'chromosome_name', 'start_position', 'end_position', 'strand')], by = 'EnsemblID')
    res.clean <- res.clean[ order(res.clean$pvalue),]  
    
    write.csv(x = res.clean,
              file=paste(loc.dexseq.folder, "/", code, "_", loc.code, "_SignificantExons.csv", sep = ''),
              row.names = FALSE)
    
    message('Saving results in ', dexseq.data)
    save(list = c('res.clean', 'DexSeqExons.loc'), file = dexseq.data)
  } else {
    load(dexseq.data)
  }


########################## Now plot a subset
  file.remove(list.files(dexseq.figs, pattern = 'DEXSeq*', full.names = TRUE)) ##remove the old plots

  n.sig <- sum(res.clean$FDR < 0.01, na.rm = TRUE)
  if (n.sig > 50) {
    res.cleanSigs <- subset(res.clean, FDR<0.01)
  } else res.cleanSigs <- res.clean[1:50,]


  genes.to.plot <- unique(res.cleanSigs$EnsemblID)
  pretty.gene.names <- as.character(annotation$external_gene_id[ match(genes.to.plot, table = annotation$EnsemblID) ])
  
  for (i in 1:length(genes.to.plot)) {
    gene <- as.character(genes.to.plot[i])
    gene.pretty <- as.character(pretty.gene.names[ i ])
    
    message(i, ' ', gene, ' ', gene.pretty)
    output.pdf <- paste(dexseq.figs, '/DEXSeq-', gene.pretty, '.pdf', sep = '')
    pdf(output.pdf, width = 8, height = 4.9)
    try(plotDEXSeq(res,
                   geneID = gene,  ##I suspect it has to be gene, otherwise it crashes
                   cex.axis = 1.2,
                   cex=1.3,
                   lwd=2,
                   legend=TRUE,
                   displayTranscripts = TRUE,
                   names = TRUE,
                   main = gene.pretty)
        )
    dev.off()
    print(output.pdf)
  }
  
  
#################################### plot some graphs
  
  pdf(file=paste(dexseq.figs, '/DEXSeq-MeanVsDispPoints.pdf', sep = ''))
  plotDispEsts (DexSeqExons.loc)
  dev.off()
  
    
  pdf(file=paste(dexseq.figs, '/DEXSeq-MeanVsDispCircles.png', sep = ''))
  plotMA(data.frame(baseMean = res.clean[,6],log2FoldChange = res.clean[,7], padj = res.clean[,5] < 0.1),
         ylim=c(-4,4), cex=0.8)
  dev.off()
  
  message('Done with ', condition)
  rm('DexSeqExons.loc')
  gc()

}

warnings()

sessionInfo()
