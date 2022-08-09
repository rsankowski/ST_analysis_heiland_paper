library(SPATAData)
library(SPATA2)


source(file.path("R","miscellaneous.R"))

sourceDataFrame()

downloadSpataObjects(sample_names = c("275_T", "334_T"), folder = "objects")

launchSpataData()

