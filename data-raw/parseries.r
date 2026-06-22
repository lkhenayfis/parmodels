library(arrow)
library(data.table)

DIR <- "/home/lkhenayfis/Documents/ONS/MODELOS/NEWAVE/modelos-comuns/dados/data/vazoes"

cods    <- c("1", "6", "8")
objects <- paste0("vazoes-codigo=", cods, ".parquet.gzip")
objects <- file.path(DIR, objects)

dummyseries <- lapply(objects, read_parquet)
dummyseries <- lapply(dummyseries, function(d) ts(d$incremental, start = c(1931, 1), freq = 12))
dummyseries <- do.call(cbind, dummyseries)
colnames(dummyseries) <- c("serie1", "serie2", "serie3")

usethis::use_data(dummyseries, overwrite = TRUE)
