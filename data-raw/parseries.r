library(aws.s3)
library(arrow)

cods    <- c("1", "6", "8")
objects <- paste0("vazoes-codigo=", cods, ".parquet.gzip")
objects <- paste0("modelos/cenarios/modelos-comuns/dados/vazoes/", objects)

BUCKET <- "s3://ons-pem-historico"

dummyseries <- lapply(objects, function(obj) s3read_using(read_parquet, object = obj, bucket = BUCKET))
dummyseries <- lapply(dummyseries, function(d) ts(d$incremental, start = c(1931, 1), freq = 12))
dummyseries <- do.call(cbind, dummyseries)
colnames(dummyseries) <- c("serie1", "serie2", "serie3")

usethis::use_data(dummyseries, overwrite = TRUE)
