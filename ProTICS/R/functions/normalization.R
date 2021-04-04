# normalization
#
# This is normalization function mapping the data to (0,1).
#
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'


normalization<-function(x) {
  min_v  <- min(x)
  max_v <- max(x)
  A<-x-replicate(dim(x)[2],min_v)
  B<-replicate(dim(x)[2],(max_v-min_v))
  return(A/B)
}
