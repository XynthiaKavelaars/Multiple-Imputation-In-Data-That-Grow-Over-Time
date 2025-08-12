# sigma.function checks positive definiteness of the correlation matrix, and 
# and makes non-positive definite matrices positive definite
sigma.function <- function(w1, w2, b1, b2, c1, c2){
  sigma <- matrix(c(1, w1, b1, c1,
                    w1, 1, c2, b2,
                    b1, c2, 1, w2,
                    c1, b2, w2, 1),
                  ncol = 4, nrow = 4, byrow=TRUE)
  e.up <- e.down <- eigen(sigma)$values
  c1.up <- c1.down <- c2.up <- c2.down <- c1
  while(any(e.down <= 0) && c1.down > -1){
    c1.down <- c2.down <- c1.down - 0.01
    sigma.down <- matrix(c(1, w1, b1, c1.down,
                           w1, 1, c2.down, b2,
                           b1, c2.down, 1, w2,
                           c1.down, b2, w2, 1),
                         ncol = 4, nrow = 4)
    e.down <- eigen(sigma.down)$values}
  while(any(e.up <= 0) && c1.up < 1){
    c1.up <- c2.up <- c1.up + 0.01
    sigma.up <- matrix(c(1, w1, b1, c1.up,
                         w1, 1, c2.up, b2,
                         b1, c2.up, 1, w2,
                         c1.up, b2, w2, 1),
                       ncol = 4, nrow = 4)
    
    e.up <- eigen(sigma)$values}
  if (!any(e.down <= 0)){
    c1 <- c1.down
    c2 <- c2.down}
  if (!any(e.up <= 0)){
    c1 <- c1.up
    c2 <- c2.up}
  sigma <- matrix(c(1, w1, b1, c1,
                    w1, 1, c2, b2,
                    b1, c2, 1, w2,
                    c1, b2, w2, 1),
                  ncol = 4, nrow = 4)
  if (any(e.up <= 0) & any(e.down <= 0)){
    sigma <- as.matrix(nearPD(sigma, corr = TRUE)$mat)}
  return(sigma)}

