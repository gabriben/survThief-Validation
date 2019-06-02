library(imager)

path <- "/Users/gabrielbenedict/Google_Drive/docs/UNIS/KU Leuven/Courses/Master Thesis/Reconstructions/reproducing survival meta analysis/Reconstructed/Benrashid"
b <- load.image(paste0(path, "/Benrashid original.jpg"))
b <- b[,,1,]




# reshape image into a data frame
df = data.frame(
  red = matrix(b[,,1], ncol=1),
  green = matrix(b[,,2], ncol=1),
  blue = matrix(b[,,3], ncol=1)
)

m <- matrix(c(c(b[,,1]), c(b[,,2]), c(b[,,3])),
            nrow = prod(dim(b)[-3]))

# K-MEANS 

# compute the k-means clustering
K = kmeans(df,4)
df$label = K$cluster

# Replace the color of each pixel in the image with the mean 
# R,G, and B values of the cluster in which the pixel resides:

# get the coloring
colors = data.frame(
  label = 1:nrow(K$centers), 
  R = K$centers[,"red"],
  G = K$centers[,"green"],
  B = K$centers[,"blue"]
)

# merge color codes on to df
# IMPORTANT: we must maintain the original order of the df after the merge!
df$order = 1:nrow(df)
df = merge(df, colors)
df = df[order(df$order),]
df$order = NULL


# get mean color channel values for each row of the df.
R = matrix(df$R, nrow=dim(b)[1])
G = matrix(df$G, nrow=dim(b)[1])
B = matrix(df$B, nrow=dim(b)[1])

# reconstitute the segmented image in the same shape as the input image
b.segmented = array(dim=dim(b))
b.segmented[,,1] = R
b.segmented[,,2] = G
b.segmented[,,3] = B

# save(b.segmented, file = "segmentedImage.rds")
# save.image(as.cimg(b.segmented[,,1,]), "segmentedImage.png")
# save(df, file = "dfSegmented.rds")



### DBSCAN
library(dbscan)
s <- df[sample(1:nrow(df), 20000),]

# choosing eps
dbscan::kNNdistplot(s, k = 5)
abline(h = 0.005, lty = 2)
abline(h = 0.015, lty = 2)
abline(h = 0.02, lty = 2)

K = dbscan(s, eps = 0.35)
print(K)

s$centerRed <- 0
s$centerGreen <- 0
s$centerBlue <- 0

for(i in 1:max(K$cluster)){
  group <- s[K$cluster == i,1:3]
  a <- array(rep(colMeans(group), nrow(group)), dim = c(3, nrow(group)))
  a <- aperm(a)  
  s[K$cluster == i, colnames(s) %in% 
      c("centerRed", "centerGreen", "centerBlue")] <- a
}

save(s, K, file = "dbscan.rds")

library(plotly)
plot_ly(data.frame(s), size = 0.1, x = ~red, y = ~green, z = ~blue, 
        marker = list(color = ~rgb(s$centerRed, s$centerGreen, s$centerBlue),
                      line = list(width = 0))) %>%
        # color = K$cluster) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'R'),
                      yaxis = list(title = 'G'),
                      zaxis = list(title = 'B'))) %>%
  layout(paper_bgcolor='rgb(254, 247, 234)')

### SPECTRAL CLUSTERING
library(kernlab)
s <- m[sample(1:nrow(m), 2000),]
system.time(
  sc <- specc(s, centers=4)
)

setwd("/Users/gabrielbenedict/Google_Drive/docs/UNIS/KU Leuven/Courses/Master Thesis/Reconstructions/reproducing survival meta analysis/")
save(s, sc, file = "sc.rds")

sc@centers

library(plotly)
# original
plot_ly(data.frame(s), size = 0.1, x = ~X1, y = ~X2, z = ~X3, 
        # marker = list(color = ~rgb(dfSmall$red, dfSmall$green, dfSmall$blue),
        #               line = list(width = 0))) %>%
  color = sc@.Data) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'R'),
                      yaxis = list(title = 'G'),
                      zaxis = list(title = 'B'))) %>%
  layout(paper_bgcolor='rgb(254, 247, 234)')


plot(m, col=sc, pch=4)            # estimated classes (x)
points(my.data, col=obj$classes, pch=5) # true classes (<>)

#### trimmed Robust Clustering


library(tclust)
s <- m[sample(1:nrow(m), 100000),]
tc <- tclust(s, k = 3, alpha = 0.001, restr = "eigen", restr.fact = 1)

centers <- cbind(tc$centers, c(0.2,0.8,0.2))
tc$cluster[tc$cluster == 0] <- ncol(tc$centers) + 1

library(plotly)
plot_ly(data.frame(s), size = 0.1, x = ~X1, y = ~X2, z = ~X3, 
            marker = list(color = ~rgb(centers[1, tc$cluster], 
                                       centers[2, tc$cluster], 
                                       centers[3, tc$cluster]),
                          line = list(width = 0))) %>%
      add_markers() %>%
      layout(scene = list(xaxis = list(title = 'R'),
                          yaxis = list(title = 'G'),
                          zaxis = list(title = 'B'))) %>%
      layout(paper_bgcolor='rgb(254, 247, 234)')


#### Linear Grouping Analysis

library(lga)
s <- m[sample(1:nrow(m), 10000),]
l <- lga(s, 3, niter = 20, biter = 300)
plot(l)

rl <- rlga(s, 3, alpha=0.85)
