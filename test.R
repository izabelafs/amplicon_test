#Izabela Silva
#30/06/2022

library(ggpubr)
library(factoextra)
library(patchwork)

##Question 1 
# 1. Write a function that generates N million random DNA fragments of length L base pairs.
#The function takes as entry the size @N number of desired sequences, @L the lenght of those sequences.

bases = c ("A", "T", "C", "G") #, "N", "a", "t", "c", "g", "n")
prob = c(0.2, 0.2, 0.2, 0.4) #, 0.05, 0.01, 0.01, 0.01, 0.01, 0.01)
N <- 100
set.seed(1000) 
L <- 100

generate_seq <- function(bases, prob, N, L) {
  for (i in 1:N) {
    sequence[i] = paste(sample(bases,
                            L,
                            replace = TRUE),
                     collapse = "")
  }
  return(sequence)
}

sequences <- generate_seq(bases,prob,N,L)


#Question 2
#Define function to count nucleotides
# 2. Write a function that takes as input the vector created in question 1 and plots the
# distribution of A, C, G and T frequencies over the N million fragments.
# 
return_n_count <- function(seq = sequences) {
  df <-
    data.frame(
      fq_a = 1:N,
      fq_c = 1:N,
      fq_g = 1:N,
      fq_t = 1:N
    )
  for (i in 1:N) {
    sub <- strsplit(seq[i], "")
      df$fq_a[i] <- str_count(sub,"A") / nchar(seq[i])
      df$fq_c[i] <- str_count(sub,"C") / nchar(seq[i])
      df$fq_g[i] <- str_count(sub,"G") / nchar(seq[i])
      df$fq_t[i] <- str_count(sub,"T") / nchar(seq[i])
    }
    return(df)
}

return_n_count(seq = sequences)

#Question 3 
# 3. Generate a matrix with 10,000 rows (gene) and 6 columns (time-points), each row value
# (score) follows a uniform distribution between 0 and 1. Plot the average gene score time-
#   series.
# 
generate_time_point <- function(N, n) {
  n_points = matrix(runif(N*n, 0,1),ncol = n)
  return(n_points)
}

n_time <- generate_time_point(N = 10000,6)

get_mean <- function(n){
for(i in 1:n){
m[i] = mean(n_time[,i])
}
return(m)
}

mean <- get_mean(n = 6)
dmean <- as.data.frame(mean) %>% rownames_to_column()
colnames(dmean) <- c("n_time","mean")

plot1 <- dmean %>%  ggplot(aes(n_time, mean)) + geom_point(na.rm=TRUE, color="blue", size=1)
  
 # 4. Bonus question – only start when you’ve finished questions 1.-3.
# Using the matrix generated previously, use your favourite method to cluster genes with
# similar time-series evolution into 6 groups. Generate a graph containing 6 sub-graphs for
# which the time-series scores of all genes of each group is plotted.

n.kmeans <- kmeans(n_time, 6,nstart = 25)

# Dimension reduction using PCA
res.pca <- prcomp(n_time,  scale = TRUE)
# Coordinates of individuals
ind.coord <- as.data.frame(get_pca_ind(res.pca)$coord)
# Add clusters obtained using the K-means algorithm
ind.coord$cluster <- factor(n.kmeans$cluster)

# Data inspection
eigenvalue <- round(get_eigenvalue(res.pca), 1)
variance.percent <- eigenvalue$variance.percent

#cluster plotting
plot2 <- ggscatter(
  ind.coord, x = "Dim.1", y = "Dim.2", 
  color = "cluster", palette = "npg", ellipse = TRUE, ellipse.type = "convex",
  size = 1.5,  legend = "right", ggtheme = theme_bw(),
  xlab = paste0("Dim 1 (", variance.percent[1], "% )" ),
  ylab = paste0("Dim 2 (", variance.percent[2], "% )" )
) +
  stat_mean(aes(color = cluster), size = 4) + plot1

print(plot2)

  