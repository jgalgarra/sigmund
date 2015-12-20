popsenvertical <- function (directory,file1,position,rs="no",grays="no",wl="",intervh=6, 
                            intervvpops=5,vcexlegend=0.6,vcaxis=0.75, legy = 1.0, legx = 1.0, lgpch = "no", legins= 0.025)
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
  par(mfrow=c(2,1))
  pospops = position
  s<-pintaspecies(directory,file1,"","Population","Plants",pospops,grays,mmaximo=mx[1],mminimo=mx[2],wlinea=wl,hr=intervh,vr=intervvpops,
                  horizlegend="no",xtitle="Years",caxis=vcaxis,cexlegend=vcexlegend, leginset = legins,
                  legy.intersp= legy, legx.intersp = legx, legpch=lgpch)
  s<-pintaspecies(directory,file2,"","Population","Pollinators",pospops,grays,mmaximo=mx[1],mminimo=mx[2],wlinea=wl,hr=intervh,vr=intervvpops,
                  horizlegend="no",xtitle="Years",caxis=vcaxis,cexlegend=vcexlegend, leginset = legins,
                  legy.intersp= legy, legx.intersp = legx,legpch=lgpch)
}