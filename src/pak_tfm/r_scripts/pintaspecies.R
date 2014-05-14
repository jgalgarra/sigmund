pintaspecies <- function(directory, filea, fileb="", ylegend, titulo, position = "",grays="no",mmaximo="",mminimo="",wlinea="",hr=5,vr=7,horizlegend="no",xtitle="",cexlegend=0.6,caxis=0.75) {
  ficheroa <- paste0(getwd(),"/",directory,"/",filea)
  if (wlinea=="")
    wl=3
  else
    wl=wlinea
  a <- read.csv(ficheroa)
  if (nchar(fileb)>0) {
    ficherob <- paste0(getwd(),"/",directory,"/",fileb)
    b <- read.csv(ficherob)
  }
  if (mmaximo=="") {
    maximob <- 0
    minimob <- 0
    maximoa <- max(a)
    minimoa <- min(a)
    if (nchar(fileb)>0)  {
      
      maximob <- max(b)
      minimob <- min(b)
    }
    maximo <-max(c(maximoa,maximob))
    minimo <-min(c(minimoa,minimob))
    if (minimo>0){
      minimo <- 0
    }
  }
  else {
    maximo <- mmaximo
    minimo <- mminimo
  }
  nspecies <- ncol(a)
  if (grays=="no"){
    colores<-rainbow(nspecies)
  }
  else {
    colfunc <- colorRampPalette(c("grey15", "grey70"))
    colores <- colfunc(nspecies)
  }
  #colores=gray(seq(0.3,0.7,length=nspecies))
  
  for (i in 1:nspecies){
    colnames(a)[i]<-paste0("species_",as.character(i))
  }
  par(mar=c(2,4,2,2))
  plot(a$species_1, col=colores[1], type="l", lwd=wl, ylim=c(minimo,maximo),xlab='',ylab='', main=titulo, font.main=20, cex.main=1.2, cex.axis=caxis, xaxt="n")
  mtext(side = 2, ylegend, line = 2.5, cex = 0.9)
  #mtext(side = 1, adj=1, xtitle,line = 0, cex = 0.9)
  intervv = vr
  dr <-min(c(ceiling(abs(log(maximo,10))),ceiling(abs(log(abs(minimo),10)))))
  if (maximo>=1){
    d<-as.integer(log(maximo,10))
    if (maximo != (maximo%/%(10^d))*10^d)
      vmax<-(1+maximo%/%(10^d))*10^d
    else
      vmax<-maximo
  }
  else
  {
    #d <- ceiling(abs(log(maximo,10)))
    vmax <- floor(abs(maximo)*100^dr)/100^dr
  }
  if (minimo >=0)
    vmin=0
  else
  {
    #d <- ceiling(abs(log(abs(minimo),10)))
    vmin <- -1*round(abs(minimo)*10^dr)/10^dr
  }
  valh <- seq(vmin,vmax,length=intervv)
  #print(valh)
  abline(h = valh, col="grey60", lty="dotted")
  abline(col="grey60", lty="dotted")
  intervh = hr
  tickhr = round(nrow(a))%/%hr
  posticks_x = seq(0,tickhr*hr,length=hr)
  
  labels_x = round(posticks_x/365)
  plast_label <- length(labels_x)
  labels_x[plast_label] <- paste(labels_x[plast_label],xtitle)
  abline(v=seq(0,round(nrow(a)),length=intervh), col="grey60", lty="dotted")
  if (nspecies > 1){
    for (i in 1:nspecies){
      lines(a[i],col=colores[i], lwd=wl)
    }
  }
  if (nchar(fileb)>0)  {
    for (i in 1:nspecies){
      lines(b[i],col=colores[i], lwd=wl, lty="dotted")
    }
  }
  
  axis(1, at = posticks_x, labels = labels_x, lwd = 1, cex.axis=caxis, xaxs="r", font=1)
  axis(2, lwd = 1, cex.axis=caxis, xaxs="r", font=1)
  #legend("topleft", c(1,nspecies), cex=0.8, col=colores,lty=1:3, lwd=2, bty="n")
  if (nchar(position)>0){
    legend(position,inset = 0.025,x.intersp=1,y.intersp=1.0,bty="o",box.lwd=0,box.col="white",horiz=(horizlegend=="yes"),as.character(seq(1,nspecies)),cex=cexlegend,fill=c(colores),border=c(colores))  }
  return(a)
}