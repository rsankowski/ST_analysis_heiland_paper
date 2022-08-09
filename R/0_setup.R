dir.create("data")
dir.create("R")
dir.create("plots")

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils', 'EBImage'))

devtools::install_github(repo = "kueckelj/confuns")
devtools::install_github(repo = "theMILOlab/SPATA2")
devtools::install_github(repo = "theMILOlab/SPATA")
devtools::install_github("theMILOlab/SPATAData")

