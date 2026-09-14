clip_q <- function(x, q = 0.95, method = "up"){
  
  stopifnot(method %in% c("up","down"))
  
  x_q = quantile(x, q)
  
  y = x
  
  if (method == "up"){
    y[y > x_q ] <- x_q
  } else {
    y[y < x_q ] <- x_q
  }
  
  return(y)
}
