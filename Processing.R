library(dplyr)
library(tidyr)
library(purrr)

resdir <- "./res/"

resfiles <- list.files(resdir)

loadres <- function(fpath, cols) {
  res <- read.delim(fpath, sep=" ", header=F)
  colnames(res) <- cols
  return(res)
}

tsdf <- map(resfiles[grepl("^twosamp", resfiles)], 
                ~loadres(paste(resdir, .x, sep=""),
                         c("job", "i", "dim", "target", "d", "ns", "wvar", "gfun", "method", "Tobs", "pval"))) %>%
  list_rbind() %>% 
  as_tibble()

print("Done two sample testing")

cpdf <- map(resfiles[grepl("^CP", resfiles)], 
             ~loadres(paste(resdir, .x, sep=""),
                      c("job", "i", "dist", "n", "nU0", "setting", "Delta", "wvar", "gfun", "dprop", "ntest", "Tobs", "pval", "B"))) %>%
  list_rbind() %>% 
  as_tibble()

print("Done changepoint")

clustdf <- map(resfiles[grepl("^zipclust", resfiles)], 
            ~loadres(paste(resdir, .x, sep=""),
                     c("job", "i", "n", "p", "wvar", "dist", "Delta", "gfun", "ntest", "method", "Tobs", "pval"))) %>%
  list_rbind() %>% 
  as_tibble()

print("Done clustering")

save(cpdf, tsdf, clustdf, file="Orthresults.rda")
