  
#' Determine the transition matrix at each h given a Guassian process
#' 
#' @param Ef GP mean
#' @param V GP vector of variances
#' @param x_grid x grid over which to define transition matrices
#' @param h_grid the harvest grid over which to compute.  Assumed to match the stock grid x_grid if not given.
#' @return a list of matrices corresponding to the state transition matrix at each harvest level h
#' @export
gp_transition_matrix <- function(Ef, V, x_grid, h_grid=NULL){  
  if(is.null(h_grid))
    h_grid <- x_grid
  lapply(h_grid, gp_calc_F, Ef, V, x_grid)
}



gp_calc_F <- function(h, Ef, V, x_grid){
  j <- which.min(abs(x_grid - h))
  indices <- 1:length(Ef)
  shift <- pmax(indices-j,1)
  mu <- Ef[shift]
  mu[shift == 1] <- 0
  F_ <- t(sapply(1:length(mu), function(i){    
    if(snap_to_grid(mu[i],x_grid) < x_grid[2]){ # 
      out <- numeric(length(x_grid))
      out[1] <- 1
      out
    } else {
      out <- dlnorm(x_grid/mu[i], 0, V[i])
    }
  }))
  F_ <- rownorm(F_)
}

#' calculate transition matrix at each harvest level based on parameterized function
#' 
#' @param f the transtion function, a function of state x, harvest h, and parameters p
#' @param p the parameters p of the transition function
#' @param x_grid the grid over which to solve 
#' @param h_grid the harvest grid over which to compute.  Assumed to match the stock grid gp$X if not given.
#' @param sigma_g the scale of the noise
#' @return a list of matrices corresponding to the state transition matrix at each harvest level h
#' @export
f_transition_matrix <- function(f, p, x_grid, h_grid = NULL, sigma_g){
  if(is.null(h_grid))
    h_grid <- x_grid
  lapply(h_grid, par_F, f, p, x_grid, sigma_g)
}

#' calculate transition matrix based on parameterized function
#' 
#' @param f the transtion function, a function of state x, harvest h, and parameters p
#' @param p the parameters p of the transition function
#' @param x_grid the grid over which to solve 
#' @param h_grid the harvest grid over which to compute.  Assumed to match the stock grid gp$X if not given.
#' @param sigma_g the scale of the noise
#' @return a list of matrices corresponding to the state transition matrix at each harvest level h
#' @export
par_F <- function(h, f, p, x_grid, sigma_g){
  mu <- sapply(x_grid, f, h, p)
  F_true <- t(sapply(mu, function(m){    
    if(snap_to_grid(m,x_grid) < x_grid[2]){ # 
      out <- numeric(length(x_grid))
      out[1] <- 1
      out
    } else {
      out <- dlnorm(x_grid/m, 0, sigma_g)
    }
  }))
  F_true <- rownorm(F_true)
}
  
                       
snap_to_grid <- function(x, grid) sapply(x, function(x) grid[which.min(abs(grid - x))])   

                       
# internal helper function
rownorm <- function(M) 
  t(apply(M, 1, function(x){ 
    if(sum(x)>0){
      x/sum(x)
    } else {
      out <- numeric(length(x))
      out[1] <- 1
      out
    }
  }))
