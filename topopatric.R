##This script generates plots for the outfiles from the
####program topopatric.f90
library(animation) # needed in sections 2, 3 and 5

#Color housekeeping
library(RColorBrewer)
##Old papers - colrd - colorblind friendly (http://colrd.com/palette/22444/)
colors <- c("#043227", "#097168", "#ffcc88", "#fa482e",  "#f4a32e")
color.function <- function(palette,n) colorRampPalette(palette)(n)

##List of parameters
#T		N		NF		deltat
#mut		diff		m		MNPM
#S		G		B
#IREAD		INIS		INISBIT         DMI
#INDEPENDENT_LOCI window		nparameter       replica
###############################################################
## SECTION 1
## Read parameter values
pars <- read.table('pop-SIM001.001.dat', nrows = 5)
N <- 1000#pars$V2[1]
S <- pars$V1[3]
B <- 1000#pars$V3[3]
relativeS <- round((as.numeric(as.character(pars$V1[3])))**2/(as.numeric(as.character(pars$V3[1])))**2,4)*100
deltat <- 50#as.numeric(as.character(pars$V4[1]))

###############################################################
## SECTION 2
##Distribution of number of compatible mates (=degree)
degree <- read.table("degree_distplotSIM007.001.dat")
filesize <- dim(degree)[1]
saveGIF({
  for (i in seq(N,filesize,N))
    {
      hist(degree[(i-999):i,1], xlab="number of compatible partners",
        main = paste(paste(paste('mating area(%) = ',relativeS,sep=''),', T = ',sep=''),(i/N)*deltat, sep=''),
        xlim=c(0,N), breaks = seq(0,N,10), col="gray",
        ylim = c(0,N))
    }
})

###############################################################
## SECTION 3
##Number of offspring vs. number of compatible mates (=degree)
saveGIF({
  for (i in seq(N,filesize,N))
  {
    plot(degree[(i-999):i,1],degree[(i-999):i,2], xlab="number of compatible partners",
         ylab = "number of offspring",
         main = paste(paste(paste('mating area(%) = ',relativeS,sep=''),', T = ',sep=''),(i/N)*deltat, sep=''),
         xlim=c(0,1000), ylim = c(0,10), pch=19)
    fit = lm(degree[(i-999):i,2] ~ degree[(i-999):i,1])
    color = "blue"
    if(i > 1) if(fit$coefficients[2] < 0) color = "red"
    if(i > 1) abline(fit,lwd=1.5,lty=2,col=color)
    text(900,10,"fit = ax+b")
    text(900,9,paste("a=",round(fit$coefficients[2],5), sep=''),col=color)
  }
})
dev.off()
prob_mating <- c(0,0,0)
for (sim in 1:7)
{
  filename <- paste("degree_distplotSIM00",sim,".001.dat", sep="")
  degree <- read.table(filename)
  filesize <- dim(degree)[1]
  i=filesize
  final_degree <- degree[(i-999):i,]
  for (i in seq(50,(N-100),100))
  {
    subset <- final_degree[which(final_degree$V1>i & final_degree$V1<=(i+50)),]
    if(dim(subset)[1] > 10)
    {
      prob <- 1-(length(which(subset$V2 == 0))/dim(subset)[1])
      prob_mating <- rbind(prob_mating, c(sim,(i+i+50)/2,prob))
    }
  }
}
prob_mating <- prob_mating[-1,]
plot(prob_mating[,1],prob_mating[,2],pch=19, xlim=c(0,N), ylim=c(0,1))
points(prob_mating[,1],prob_mating[,2],pch=19, xlim=c(0,N), ylim=c(0,1), col='red')
prob_mating <- as.data.frame(prob_mating)
names(prob_mating) <- c('c','k','prob')

library(ggplot2)
c.discrete <- factor(prob_mating$c, 
                     levels=unique(as.character(prob_mating$c)), 
                     unique(as.character(prob_mating$c)))
ggplot(prob_mating, aes(x=k, y=prob, color=c.discrete)) + 
  geom_line(aes(group=c)) +
  labs(color="number of attempts") + 
  xlab('Number of compatible mates') + 
  ylab('Probability of finding a compatible mate') + 
  xlim(0,N) + ylim(0,1.02) +
  geom_hline(yintercept = 1, color="black",
             linetype="dashed", size=0.8) +
  geom_vline(xintercept = 0, color="black",
             linetype="dashed", size=0.8) +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
        legend.position = "bottom", 
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA)
  )
###############################################################
## SECTION 4
##Fst
#Calculate the threshold for Fst outliers
fstNull <- read.table("FstNullSIM001.001.dat")
hist(as.numeric(fstNull[1,]))
outlier <- tail(sort(as.numeric(fstNull[1,])),(length(fstNull[1,])*0.01))[1]
#Plot Fst values along the genotype
fst <- read.table("FstSIM001.001.dat")
hist(as.numeric(fst[1,]))
plot(seq(1,B),fst[1,1:B],pch=19, xlab="position in genome", ylab="Fst")
abline(h=outlier,lty=2)

###############################################################
## SECTION 5
##Distribution of genetic distances
dist <- read.table('distSIM007.001.dat')
nvalues <- N*(N-1)/2
filesize <- dim(dist)[1]
saveGIF({
  for (i in seq(nvalues,filesize,nvalues))
  {
    hist(sample(dist[(i-(nvalues-1)):i,1],10000), xlab="genetic distance", breaks=B,
         main = paste(paste(paste('mating area(%) = ',relativeS,sep=''),', T = ',sep=''),(i/nvalues)*deltat, sep=''),
         col="gray", xlim=c(0,B))
  }
})

###############################################################
## SECTION 6
##Plot individuals in space colored by species
speciesplot <- read.table('speciesplotSIM005.001.dat')
nspecies <- max(speciesplot$V1)
sspcolors <- color.function(colors, nspecies)
plot(speciesplot$V1, speciesplot$V2, col=sspcolors[speciesplot$V3], pch=16, cex=1.1, axes=FALSE,asp=1, xlab="", ylab="")

ssp1g <- (speciesplot[which(speciesplot[,1]==1),4:1003])
ssp2g <- (speciesplot[which(speciesplot[,1]==2),4:1003])
ssp3g <- (speciesplot[which(speciesplot[,1]==3),4:1003])
ssp1m <- (speciesplot[which(speciesplot[,1]==1),1004:2003])
ssp2m <- (speciesplot[which(speciesplot[,1]==2),1004:2003])
ssp3m <- (speciesplot[which(speciesplot[,1]==3),1004:2003])


hist(rowSums(ssp2m),xlim=c(0,150),col="gray",ylim=c(0,50))
hist(rowSums(ssp2g), col="blue",add=TRUE)

ssp1gP <- colSums(ssp1g)/dim(ssp1g)[1]
length(which(ssp1gP == 0 | ssp1gP == 1))
ssp1mP <- colSums(ssp1m)/dim(ssp1m)[1]
length(which(ssp1mP == 0 | ssp1mP == 1))

ssp2gP <- colSums(ssp2g)/dim(ssp2g)[1]
Pn <- length(which(ssp2gP == 0 | ssp2gP == 1))
ssp2mP <- colSums(ssp2m)/dim(ssp2m)[1]
Ps <- length(which(ssp2mP == 0 | ssp2mP == 1))

ssp3gP <- colSums(ssp3g)/dim(ssp3g)[1]
length(which(ssp3gP == 0 | ssp3gP == 1))
ssp3mP <- colSums(ssp3m)/dim(ssp3m)[1]
length(which(ssp3mP == 0 | ssp3mP == 1))

Dn <- length(which(abs(ssp2gP - ssp3gP) == 1))
Ds <- length(which(abs(ssp2mP - ssp3mP) == 1))

(Pn/Ps)/(Dn/Ds)

(Dn/(Dn+Ds)) - (Pn/(Pn+Ps))
hist(test, add=TRUE, col="blue")
test2 <- rowSums(speciesplot[which(speciesplot[,3]==1),1004:2003])
hist(test2, col="gray" )

##Bootstrap example
# Bootstrap 95% CI for R-Squared
library(boot)
# function to obtain R-Squared from the data 
rsq <- function(formula, data, indices) {
  d <- data[indices,] # allows boot to select sample 
  fit <- lm(formula, data=d)
  return(summary(fit)$r.square)
} 
# bootstrapping with 1000 replications 
results <- boot(data=mtcars, statistic=rsq, 
                R=1000, formula=mpg~wt+disp)

# view results
results 
plot(results)

# get 95% confidence interval 
boot.ci(results, type="bca")
