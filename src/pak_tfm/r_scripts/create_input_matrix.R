library(readr)



create_interaction_mtx <- function(file="",Factorbij = 0.0004, alfa = 0.0006, ci = 0.00002, r_b = 0, 
                                   r_d = 0.08, min_mort = 0.12, ini_pop = 10000, transposeinp = FALSE)
{

options(scipen=999)
  
  create_aux <- function(x_data)
  {
    nrow_aux <- 5
    sum_cols <- colSums(x_data)#/sum(x_data)
    promedio_cols <- mean(sum_cols)
    am <- matrix( rep(0,ncol(x_data)*nrow_aux),nrow=nrow_aux,ncol=ncol(x_data))
    # alpha
    am[2,] <- as.numeric(alfa)
    # ci
    for (j in 1:ncol(am))
      am[3,j] <- ci*sum_cols[j]/promedio_cols
    # r_b
    am[4,] <- r_b
    # initial populations
    am[1,] <- ini_pop
    # r_d
    for (j in 1:ncol(am))
      #am[5,] <- max(min_mort,as.numeric(sum_cols*r_d/promedio_cols))
      am[5,j] <- min(1,max(min_mort,as.numeric(sum_cols[j]*r_d/promedio_cols)))
    return(am)
  }
    
  input_data <- read_csv(paste0("../input/csvs/",file,".csv"))
  matrix_data <- unname(as.matrix(input_data[1:nrow(input_data),2:ncol(input_data)]))
  if (transposeinp)
    matrix_data <- t(matrix_data)
  num_a	<- nrow(matrix_data)
  num_b <- ncol(matrix_data)
  
  
  a_data <- matrix_data
  a_data <- log10(1+a_data)*Factorbij
  b_data <- t(matrix_data)
  b_data <- log10(1+b_data)*Factorbij
  
  aux_a_data <- create_aux(a_data)
  matrix_a <- rbind(a_data,aux_a_data)
  write.table(matrix_a,paste0("../input/",file,"_AUTO_a.txt"),col.names = FALSE, row.names=FALSE, sep = "\t")
  
  aux_b_data <- create_aux(b_data)
  matrix_b <- rbind(b_data,aux_b_data)
  write.table(matrix_b,paste0("../input/",file,"_AUTO_b.txt"),col.names = FALSE, row.names=FALSE, sep = "\t")
}

#create_interaction_mtx(file = "burkle", transposeinp = TRUE)
create_interaction_mtx(file = "M_PL_059")