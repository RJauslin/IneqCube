
x <- seq(1,6,1)
n <- length(x)
m <- 3
a <- seq_len(m)
e <- 0
h <- m
len.r <- length(r <- x[a])
count <- as.integer(round(choose(n, m)))
out <- matrix(r, nrow = len.r, ncol = count)
if (m > 0) {
  i <- 2
  nmmp1 <- n - m + 1
  # print(nmmp1)
  # while (a[1] != nmmp1) {
  while(i < 4){
    print(e)
    print(n-h)
    if (e < n - h) {
      # print(h)# only one element to change and it is the last element
      h <- 1
      e <- a[m]
      j <- 1
      # print(e + j)
      # print(m - h + j)
      a[m - h + j] <- e + j
    }
    else {
      print("kjsdf")
      e <- a[m - h]
      h <- h + 1
      j <- 1:h
      # print(j);cat('\n')
      a[m - h + j] <- e + j
    }

    # a[m - h + j] <- e + j
    # print(m - h + j)
    # print(a[m - h + j])
    # print(m - h + j)
    r <- x[a]
    out[, i] <- r
    i <- i + 1
  }
}
