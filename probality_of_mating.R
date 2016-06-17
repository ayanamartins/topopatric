#Expression to calculate mating probability as functions of the number
## of compatible matings and the number of trials
rm(list = ls())
single.term <- function(i) {
  t1 <- k/(n-(i-1))
  if(i==1)
    prod.terms <- 1
  if(i>1)
    prod.terms <- sapply(1:(i-1),
                         function(j) 1-k/(n-(j-1)))
  t1*prod(prod.terms)
}

n <- 5
k <- 2
c <- 2

sum <- 0
for(i in 1:c)
{
  print(single.term(i))
  sum <- sum+single.term(i)
}
sum

n <- 10000
cvalues <- c(1,2,5,10,20,100,200,500)
out <- c(0,0,0)
for (c in cvalues)
{
  for(k in 1:n)
  {
    sum <- 0
    for(i in 1:c)
    {
      sum <- sum+single.term(i)
    }
    out <- rbind(out,c(c,k,sum))
  }
}
out <- as.data.frame(out[-1,])
names(out) <- c('c','k','prob')


library(ggplot2)
c.discrete <- factor(out$c, 
                     levels=unique(as.character(out$c)), 
                     unique(as.character(out$c)))
ggplot(out, aes(x=k, y=prob, color=c.discrete)) + 
  geom_line(aes(group=c)) +
  labs(color="number of attempts") + 
  xlab('Number of compatible mates') + 
  ylab('Probability of finding a compatible mate') + 
  xlim(0,n) + ylim(0,1.02) +
  geom_hline(yintercept = 1, color="black",
             linetype="dashed", size=0.8) +
  geom_vline(xintercept = 0, color="black",
             linetype="dashed", size=0.8) +
  annotate("text", x = 6000, y = 0.5, size=3.5, label = "blind case", color='red') +
  annotate("text", x = 1200, y = 1.02, size=3.5, label = "fully-aware case")
+
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
        legend.position = "bottom", 
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA)
        )
save.image(file="matprob.RData")
