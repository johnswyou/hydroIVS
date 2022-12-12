#-------------------------------------------------------------------------------
knn_max_dist <- function(X,k) {

  # Finds the maximum k nearest neighbour distance to each data point
  # (excluding itself)
  max_dists = RANN::nn2(X,k=k,treetype="kd",searchtype="priority")$nn.dists[,k]

  return(max_dists)
}
#-------------------------------------------------------------------------------
knn_radius <- function(X,r) {

  # Finds which data points are within a given radius of a particular data
  # point
  X = as.matrix(X)
  np = nrow(X)

  # preallocate space
  npoints_radius = vector(mode="numeric",length=np)

  for(i in 1:np){

    # use fact that indices of points NOT within raidus are zero
    npoints_radius[i] = sum( RANN::nn2(X,t(as.matrix(X[i,])), k=np,
                                       treetype="kd",searchtype = "radius",
                                       radius=r[i])$nn.idx > 0)
  }

  return(npoints_radius)
}
