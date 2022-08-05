
mru.start <- function( X, G, m = 2, start = c( "random", "da" ) )
# return starting values for B and V for the multi-nomial restricted unfolding model
# random = random starting values from uniform distribution
# da     = discriminant analysis start
{
  n <- nrow( X )
  p <- ncol( X )
  c <- ncol( G )
  if( start == "random" ) {
    B <- matrix( runif( p * m ), p, m )
    V <- matrix( runif( c * m ), c, m )
  }
  if ( start == "da" ) {
    tGGinv <- solve( t( G ) %*% G )
    U <- t( X ) %*% G %*% tGGinv
    e <- eigen( ( 1 / n ) * t( X ) %*% X )
    Tmp <- e$vectors %*% diag( sqrt( 1 / e$values ) )
    A <- t( U ) %*% Tmp
    s = svd( A, nu = m, nv = m )
    sV <- matrix( s$v, p, m )
    B <- Tmp %*% sV 
    V <- tGGinv %*% t( G ) %*% X %*% B
  }
  return( list( B = B, V = V ) )
} # mru.start

mru.finish <- function( X, B, V )
# return multi-nomial restricted unfolding solution rotated to principal axes  
{
  U <- X %*% B
  A <- t( U ) %*% U
  e <- eigen( A )
  R <- e$vectors
  B <- B %*% R
  sb = sign(B[1, ])
  B = outer(rep(1, nrow(B)), sb) * B
  V <- V %*% R
  V = outer(rep(1, nrow(V)), sb) * V
  return( list( B = B, V = V ) )
} # mru.finish

mru.random <- function( X, G, m = 2 )
{
  s <- mru.start( X, G, m = m, start = "random" )
  r <- fastmru( G, X, B = s$B, V = s$V)
  # r <- fastmru( G, X, B = s$B, V = s$V, MAXINNER = 5)
  f <- mru.finish( X, r$B, r$V )
  return( list( B = f$B, V = f$V, deviance = r$deviance ) )
} # mru.random

mru.user <- function( X, G, m = 2, B.start, V.start )
{
  r <- fastmru( G, X, B = B.start, V = V.start)
  f <- mru.finish( X, r$B, r$V )
  return( list( B = f$B, V = f$V, deviance = r$deviance ) )
} # mru.user

mru.da <- function( X, G, m = 2 )
{
  s <- mru.start( X, G, m = m, start = "da" )
  r <- fastmru( G, X, B = s$B, V = s$V)
  f <- mru.finish( X, r$B, r$V )
  return( list( B = f$B, V = f$V, deviance = r$deviance ) )
} # mru.da

mru.plot2 <- function( X, G, out )
{
  B <- out$B
  V <- out$V
  
  NN <- as.data.frame( X %*% B )
  colnames( NN ) <- c( "Dim1", "Dim2" )
  
  VV <- as.data.frame( V )
  colnames( VV ) <- c("Dim1", "Dim2")
  rownames( VV ) <- colnames( G )

  library( ggplot2 )  
  
  p <- ggplot() +
    geom_point( data = VV, aes( x = Dim1, y = Dim2 ), colour = "darkgreen", size = 3 ) +
    geom_point( data = NN, aes( x = Dim1, y = Dim2, color = y ), size = 1 ) +
    xlab( "Dimension 1" ) +
    ylab( "Dimension 2" ) +
    xlim( -12, 12 ) + ylim( -12, 12 ) +
    coord_fixed()
  
  #p <- p + geom_text( data = VV, aes( x = Dim1, y = Dim2 ),
  #                    label = rownames( V ), vjust = 0, nudge_y = -0.5 )

  p <- p + geom_abline( intercept = 0, slope = B[,2] / B[,1], colour = "lightskyblue" )
  
  idx1 <- apply( abs( B ), 1, which.max )
  s <- rep( NA, 3 )
  t <- rep( NA, 3 )
  
  for( pp in 1:3 ) {
    t[(pp)] = 12.4 / ( abs( B[pp,idx1[pp]] ) ) * B[pp,-idx1[pp]]
    s[(pp)] = sign( B[pp,idx1[pp]] )
  }
  
  CC <- cbind( idx1, t, s )
  bottom <- which( CC[, "idx1"] == 2 & CC[, "s"] == -1 )
  top <-  which( CC[, "idx1"] == 2 & CC[, "s"] == 1 )
  right <- which( CC[, "idx1"] == 1 & CC[, "s"] == 1 )
  left <- which( CC[, "idx1"] == 1 & CC[, "s"] == -1 )
  
  p <- p + scale_x_continuous( limits = c(-12,12), 
                               breaks = CC[bottom, "t"], 
                               labels = xnames[bottom],
                               sec.axis = sec_axis( trans ~ ., breaks = CC[top, "t"], labels = xnames[top] ) )
  
  p <- p + scale_y_continuous( limits = c(-12,12), 
    breaks = CC[left, "t"], 
    labels = xnames[left],
    sec.axis = sec_axis( trans ~ ., breaks = CC[right, "t"], labels = xnames[right] ) )
  
  plot( p )
} # mru.plot2
