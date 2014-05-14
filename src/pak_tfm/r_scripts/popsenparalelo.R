popsenparalelo <- function (directory,file1,position,rs="no",grays="no",wl="",intervh=6, intervvpops=5, intervvrs =5, forcermax="", forcermin="",vcexlegend=0.6,vcaxis=0.75)
{
  buscamax <- function (directory,filea,fileb)
  {
    ficheroa <- paste0(getwd(),"/",directory,"/",filea)
    ficherob <- paste0(getwd(),"/",directory,"/",fileb)
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
    return(c(maximo,minimo))
  }
  #print(file1)
  file2 = gsub("_a_","_b_",file1,fixed=TRUE)
  mx <- buscamax (directory,file1,file2)
  #print(file2)
  par(mfrow=c(2,2))
  pospops = position
  if (rs=="yes")
  {
    file3 = gsub("_populations_","_rs_",file1,fixed=TRUE)
    file4 = gsub("_a_","_b_",file3,fixed=TRUE)
    mxr <- buscamax (directory,file3,file4)
    #print(file3)
    #print(file4)
    par(mfrow=c(2,2))
    pospops <- ""
  }
  s<-pintaspecies(directory,file1,"","Population","Plants",pospops,grays,mmaximo=mx[1],mminimo=mx[2],wlinea=wl,hr=intervh,vr=intervvpops,xtitle="Years",caxis=vcaxis,cexlegend=vcexlegend)
  s<-pintaspecies(directory,file2,"","","Pollinators",pospops,grays,mmaximo=mx[1],mminimo=mx[2],wlinea=wl,hr=intervh,vr=intervvpops,xtitle="Years",caxis=vcaxis,cexlegend=vcexlegend)
  if (rs=="yes")
  {
    if (forcermax!="")
      mxr[1]<-as.numeric(forcermax)
    if (forcermin!="")
      mxr[2]<-as.numeric(forcermin)
    s<-pintaspecies(directory,file3,"","Effective rates","",position,grays,mmaximo=mxr[1],mminimo=mxr[2],wlinea=wl,hr=intervh,vr=intervvrs,xtitle="Years",caxis=vcaxis,cexlegend=vcexlegend)
    s<-pintaspecies(directory,file4,"","","",position,grays,mmaximo=mxr[1],mminimo=mxr[2],wlinea=wl,hr=intervh,vr=intervvrs,xtitle="Years",caxis=vcaxis,cexlegend=vcexlegend)
  }
}