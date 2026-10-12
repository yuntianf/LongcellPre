args <- commandArgs(trailingOnly=TRUE)
mode <- "updated"
root <- normalizePath(args[2])
stopifnot(!length(find.package("Longcellsrc", quiet=TRUE)))
library(LongcellPre)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(429)
iso <- c('1,10','30,70','1,10|30,70','2,10|30,69','1,10|50,70','100,200')
out <- lapply(1:50,function(i) LongcellPre:::isoset_mid_diff(sample(iso,12,TRUE),sample(iso,16,TRUE),
                       sample(0:3,1),sample(c(0,.5,1),1),0,'|',','))
  stopifnot(identical(out,readRDS(args[2])))
  cat('PASS: 50 exact isoform-match table comparisons.\n')
