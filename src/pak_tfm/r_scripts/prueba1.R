prueba1 <- function(directory, filea, fileb="", ylegend, titulo, position,grays="no") {
  ficheroa <- paste0(getwd(),"/",directory,"/",filea)
  #a<-list.files(directorio)
  #print(fichero)
  maximob <- 0
  minimob <- 0
  a <- read.csv(ficheroa)
  maximoa <- max(a)
  minimoa <- min(a)
  if (nchar(fileb)>0)  {
    ficherob <- paste0(getwd(),"/",directory,"/",fileb)
    b <- read.csv(ficherob)
    maximob <- max(b)
    minimob <- min(b)
  }
  maximo <-max(c(maximoa,maximob))
  minimo <-min(c(minimoa,minimob))
  if (minimo>0){
    minimo <- 0
  }
  nspecies <- ncol(a)
  if (grays=="no"){
    colores<-rainbow(nspecies)
  }
  else {
    colfunc <- colorRampPalette(c("grey15", "grey90"))
    colores <- colfunc(nspecies)
  }
  #colores=gray(seq(0.3,0.7,length=nspecies))
  
  for (i in 1:nspecies){
    colnames(a)[i]<-paste0("species_",as.character(i))
  }
  plot(a$species_1, col=colores[1], type="l", lwd=2, ylim=c(minimo,maximo),xlab='Days',ylab=ylegend, main=titulo, font.main=20, cex.main=1.5, cex.axis=0.8)
  intervv = 7
  abline(h=seq((minimo*10)%/%10,maximo,length=intervv), col="grey60", lty="dotted")
  intervh = 5
  periodoh = round(nrow(a)/intervh)
  abline(v=seq(0,round(nrow(a)),length=intervh), col="grey60", lty="dotted")
  #axis(1, lwd = 2, cex.axis=1)
  if (nspecies > 1){
    for (i in 1:nspecies){
      lines(a[i],col=colores[i], lwd=2)
    }
  }
  if (nchar(fileb)>0)  {
    for (i in 1:nspecies){
      lines(b[i],col=colores[i], lwd=2, lty="dotted")
    }
  }
  
  axis(1, lwd = 1, cex.axis=0.8, xaxs="r", font=2)
  axis(2, lwd = 1, cex.axis=0.8, xaxs="r", font=2)
  #legend("topleft", c(1,nspecies), cex=0.8, col=colores,lty=1:3, lwd=2, bty="n")
  legend(position,horiz=FALSE,as.character(seq(1,nspecies)),cex=0.8,fill=c(colores))
  return(a)
}