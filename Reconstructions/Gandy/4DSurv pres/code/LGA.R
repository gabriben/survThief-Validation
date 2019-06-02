library(imager)
b <- load.image("Benrashid original.jpg")
b <- b[,,1,]

# reshape image into a matrix
m <- matrix(c(c(b[,,1]), c(b[,,2]), c(b[,,3])),
            nrow = prod(dim(b)[-3]))

#### Linear Grouping Analysis
library(lga)
s <- m[sample(1:nrow(m), 10000),]
l <- lga(s, 3, niter = 20, biter = 300)
plot(l) #does not seem to find three lines

rl <- rlga(s, 3, alpha=0.85) # runs into an error
