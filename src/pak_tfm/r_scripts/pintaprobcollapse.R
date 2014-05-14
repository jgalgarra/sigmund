pintaprobcollapse <- function(directory, filea) {
  ficheroa <- paste0(getwd(),"/",directory,"/",filea)
  a <- read.csv(ficheroa)
  nspecies <- ncol(a)-1
  colores<-rainbow(nspecies)
  for (i in 1:nspecies){
    colnames(a)[i+1]<-paste0("species_",as.character(i))
  }
  plot(a, type="l", lwd=2, xlab='',ylab='', main=filea, font.main=20, cex.main=0.7)
  mtext(side = 2, "Prob. network collapse", line = 2.5, cex = 0.8)
  mtext(side = 1, "Blossom probability",line = 2.0, cex = 0.8)
  return(a)
}