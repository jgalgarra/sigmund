comparaoutcomes <- function (directory,file1,leyenda="",position="",rs="no",grays="no")
{
  file2 = gsub("_abs_","_u_",file1,fixed=TRUE)
  file3 = gsub("_populations_","_rs_",file1,fixed=TRUE)
  file4 = gsub("_a_","_b_",file3,fixed=TRUE)
  print(file1)
  print(file2)
  print(file3)
  print(file4)
  par(mfrow=c(2,1))
  mx<-pintaspecies(directory,file1,file2,"Individuals",leyenda,position,grays)
  print(mx)
  s<-pintaspecies(directory,file3,file4,"rs",leyenda,position,grays)
}