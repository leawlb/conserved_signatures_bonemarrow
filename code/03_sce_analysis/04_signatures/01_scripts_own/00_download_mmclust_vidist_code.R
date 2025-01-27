

# download R script containing code for a function needed for reclustering
# score calculation

# the function is called vi.dist and is from the mccluste package
# https://doi.org/10.32614/CRAN.package.mcclust
# https://github.com/cran/mcclust

set.seed(37)
options(timeout=1000)


utils::download.file(
  url = "https://github.com/cran/mcclust/blob/master/R/vi.dist.R",
  destfile = snakemake@output[["path_output"]])


utils::download.file(
  url = "https://github.com/cran/mcclust/archive/refs/heads/master.zip",
  destfile = "/home/l012t/test/")
