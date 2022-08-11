CTL_genes = c('CD8A','CD8B','GZMA','GZMB','PRF1')
Merck_genes = c('TIGIT', 'CD27','CD8A','PDCD1LG2','LAG3','CD274','CXCR6','CMKLR1','NKG7','CCL5','PSMB10','IDO1','CXCL9','HLA-DQA1','CD276','STAT1','HLA-DRB1','HLA-E')
Antigen_pre_genes = c('B2M','HLA-A','HLA-B','HLA-C','HLA-E',
            'HLA-F','HLA-G','HLA-F','HLA-G','HLA-H','TAP2')
Szabo_Inflammation_genes = c('CD274', 'CD8A', 'LAG3','STAT1')
LKB1_loss_genes = c('AVPI1', 'BAG1', 'CPS1', 'DUSP4', 'FGA','GLCE', 
              'HAL', 'IRS2', 'MUC5AC', 'PDE4D', 'PTP4A1', 
              'RFK', 'SIK1','TACC2', 'TESC', 'TFF1')
NRF2_score_genes = c('ABCB6', 'AKR1C1', 'AKR1C2', 'AKR1C3', 
               'TXNRD1', 'GCLM', 'SRXN1', 'PGD', 'TRIM16L', 'TRIM16', 
               'G6PD', 'OSGIN1', 'TALDO1', 'NQO1', 'NR0B1', 
               'UGDH', 'GSR', 'CYP4F3', 'AKR1B10')
Neuroendocrine_genes = c('AGXTs2L1', 'COL25A1', 'SLC14A2', 'ODC1', 
                   'ENO3', 'ASCL1', 'CALCA', 'CALCB', 'RET', 
                   'ZMAT4', 'MTMR7', 'SLC38A8', 'CELF3', 
                   'KLK11', 'KLK13', 'KLK12', 'KLK14', 
                   'SEC11C', 'BAALC', 'ERO1LB', 'CTNND2', 'UGT3A1')
IFNG_genes =  c('IFNG','STAT1','IDO1','CXCL10','CXCL9','HLA-DRA')

Surface_recepter_genes <- read.table('gene_markers/Surface_gene_receptor.list', sep=',', header = FALSE, stringsAsFactors = FALSE)
Surface_ligand_genes <- read.table('gene_markers/Surface_gene_ligand.list', sep=',', header = FALSE, stringsAsFactors = FALSE)
Antigen_genes <- read.table('gene_markers/Antigen_processing_presentation.list', sep=',', header = FALSE, stringsAsFactors = FALSE)
P53_pathway_genes <- read.table('gene_markers/P53_pathway.list', sep=',', header = FALSE, stringsAsFactors = FALSE)
DNA_repair_genes <- read.table('gene_markers/DNA_repair.list', sep=',', header = FALSE, stringsAsFactors = FALSE)





get_markers <- function(exprsn) {
  #output = list()
  CTL <- colMeans(exprsn[intersect(rownames(exprsn),CTL_genes),])
  Merck <- colMeans(exprsn[intersect(rownames(exprsn),Merck_genes),])
  Inflammation <- colMeans(exprsn[intersect(rownames(exprsn),Szabo_Inflammation_genes),])
  LKB1_loss <- colMeans(exprsn[intersect(rownames(exprsn),LKB1_loss_genes),])
  NRF2 <- colMeans(exprsn[intersect(rownames(exprsn),NRF2_score_genes),])
  Neuroendocrine <- colMeans(exprsn[intersect(rownames(exprsn),Neuroendocrine_genes),])
  IFNG <- colMeans(exprsn[intersect(rownames(exprsn),IFNG_genes ),])
  Surface_recepter <- colMeans(exprsn[intersect(rownames(exprsn),Surface_recepter_genes$V1),])
  Surface_ligand <- colMeans(exprsn[intersect(rownames(exprsn),Surface_ligand_genes$V1),])
  Antigen_processing <- colMeans(exprsn[intersect(rownames(exprsn),Antigen_genes$V1),])
  Antigen_presentation <- colMeans(exprsn[intersect(rownames(exprsn),Antigen_pre_genes),])
  P53 <- colMeans(exprsn[intersect(rownames(exprsn),P53_pathway_genes$V1),])
  DNA_repair <- colMeans(exprsn[intersect(rownames(exprsn),DNA_repair_genes$V1),])
  
  
  output <- data.frame(CTL,Merck,Inflammation,LKB1_loss,NRF2,Neuroendocrine,IFNG,Surface_recepter,
                       Surface_ligand,Antigen_processing,Antigen_presentation,P53,DNA_repair)
  return (output)
}


suppressMessages(library(optparse))
suppressMessages(library (plyr))
suppressMessages(library(pROC))
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="list input files", metavar="character"),
  make_option(c("-d", "--design"), type="character", default=NULL,
              help="design column", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL,
              help="output prefix name", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

input <- opt$input
design <- opt$design
outdir <- opt$outdir

#input <- "/Users/yangliu/Dropbox (Partners HealthCare)/CIDC Trial Analyses/10026/RNAseq_updates/10026_cohort/matrix/Noreplicates_False_tpm.genesymbol.csv"

exprsn <- read.table(file = input, sep='\t', header = TRUE, stringsAsFactors = FALSE, row.names = 1)
exprsn <- log2(exprsn  + 1)
exprsn_norm <- exprsn - rowMeans(exprsn)

output <- get_markers(exprsn_norm)
outfile <- paste(outdir,design,"_signature_output.txt",sep='')
write.table(output,outfile,sep='\t',quote=FALSE,row.names = T,col.names = T)



