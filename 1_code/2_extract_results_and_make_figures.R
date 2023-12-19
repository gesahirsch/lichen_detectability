library(jagsUI)
library(tidyverse)
library(scales)
library(ggplot2)
library(ggridges)
library(viridis)
library(gridExtra)
library(LaplacesDemon)

load('0_data/data.RData')
load('2_output/output.jags_model.RData')

# explore output
print(output,2)
paste("running time:",round(output$mcmc.info$elapsed.mins,0),"min =",round(output$mcmc.info$elapsed.mins/60,1),"hours")
ifelse(sum(unlist(output$Rhat) > 1.1)==0, paste("Successful convergence based on Rhat values (all < 1.1)."), paste("Rhat values indicate convergence failure."))
# traceplot(output)

# extract chains
out <- data.frame(output$sims.list)
params.titles <- colnames(out)

# parameters that did not converge
params.titles[unlist(output$Rhat) > 1.1]

# MCMC settings & traced params
nc=output$mcmc.info$n.chains; ni=output$mcmc.info$n.iter; nb=output$mcmc.info$n.burnin; nt=output$mcmc.info$n.thin; na=output$mcmc.info$n.adapt
npost <- nc*(ni-nb)/nt # number of posterior samples
paste('burn-in:',nb,' total number of iterations:',ni)
params <- output$parameters



#  Table 1 - Point estimates and CRI of main parameters  ####

# Occupancy process
output$summary[1:15,c(1,3,7)] %>% round(2)

# Detection process
output$summary[16:32,c(1,3,7)] %>% round(2)

table1 <- output$summary[1:32,c(1,3,7)] %>% round(2)

# export
write.table(table1, file="2_output/Table.1_Posterior_estimates_main_parameters.csv",sep=";", row.names=TRUE)




#  Figure 1 - Map (ArcGIS)  ####



#  Figure 2 - Ridgeplot with observer-specific mean and sd detection probability  ####

start <- match("mu.beta1.1", colnames(out))
end <- match("mu.beta1.7", colnames(out))

order.obs <- out[,start:end] %>% 
  apply(2, median) %>% 
  order()
p <- out[,start:end][,order.obs] %>% 
  unlist() %>% 
  plogis()
df <- data.frame(p=p, observer=rep(paste('observer',order.obs), each=npost))
df$observer <- factor(x=df$observer, levels=paste('observer',order.obs))
head(df)

start <- match("sd.beta1.1", colnames(out))
end <- match("sd.beta1.7", colnames(out))
sd <- out[,start:end][,order.obs] %>% 
  unlist()
df$spread <- p %>% 
  qlogis() %>% 
  + rnorm(length(sd), mean=0, sd=sd) %>% 
  plogis()
df$variance <- out[,start:end]^2 %>%
  unlist


mapPalette <- colorRampPalette(c('khaki1','khaki3','khaki4','mediumorchid1','mediumorchid3','mediumorchid4'))

# Mean detection probabilities per observer
plot1 <- ggplot(df, aes(x=p, y=observer, fill=observer)) +
  xlab('Estimated detection probability') + ylab('') +
  geom_density_ridges(alpha=0.8) +
  scale_fill_manual(values=mapPalette(7)) +
  theme_ridges() +
  theme(legend.position='none'
        ,axis.title=element_text(size=65) # 20
        ,axis.title.x=element_text(margin = unit(c(10, 0, 0, 0), "mm"))
        ,axis.text=element_text(size=60) # 15
        ,plot.tag=element_text(size=80)
  ) +
  labs(tag='A')

# Detection variance across species per observer
plot2 <- ggplot(df, aes(x=variance, y=observer, fill=observer)) +
  xlab('Variance across species') + ylab('') +
  geom_density_ridges(alpha=0.8) +
  scale_fill_manual(values=mapPalette(7)) +
  theme_ridges() +
  theme(legend.position='none'
        ,axis.title=element_text(size=65) # 20
        ,axis.title.x=element_text(margin = unit(c(10, 0, 0, 0), "mm"))
        ,axis.text=element_text(size=60) # 15
        ,plot.tag=element_text(size=80)
  ) +
  scale_x_continuous(#breaks=seq(0,1,0.2),
                     limits=c(NA,3)) +
  labs(tag='B')

# dev.off()
tiff(paste('2_output/Fig.2_Observer-specific_detectability_and_variance.tif'), width=2286, height=1440, units='px', compression='lzw')
grid.arrange(plot1, plot2, ncol=2)
dev.off()




#  Figure 3 - Line plot with effect of identifiability, experience, conspicuousness  ####

start <- match("mu.beta1.1", colnames(out))
end <- match("mu.beta1.7", colnames(out))

x <- c(-1, -0.5, 0, 0.5, 1)
y <- array(data=NA, dim=c(2, length(x), 2, nrow(out)),
           dimnames=list(c('consp0','consp1'),
                         paste('identif',x,sep=''),
                         c('exper0','exper1'),
                         NULL)); dimnames(y)

# dev.off()
tiff(paste('2_output/Fig.3_Effect_of_identifiability_experience_and_conspicuousness.tif'), width=1200, height=800, units='px', compression='lzw')
par(mfrow=c(1,1), mar=c(8,10,2,1))
plot(x=x, y=x, col='white', xlim=c(-1.05, 1.05), ylim=c(0.2, 0.83), axes=F, xlab='', ylab='')
for(c in 1:2){
  for(e in 1:2){
    for(i in 1:length(x)){
      y[c,i,e,] <- out[,start:end] %>%
        apply(1, mean) %>% 
        + (out[,'beta2']*(c-1)) %>% 
        + (out[,'beta3']*x[i]) %>% 
        + (out[,'beta4']*(e-1)) %>% 
        plogis()
    }
    polygon(x=c(x, rev(x)),
            y=c(apply(y[c,,e,],1,quantile,probs=0.025), rev(apply(y[c,,e,],1,quantile,probs=0.975))),
            border=NA, col=alpha(c('khaki2','mediumorchid1')[e], c(0.2,0.1)[e]))
    lines(x=x, y=apply(y[c,,e,],1,quantile,probs=0.5), lwd=3,
          col=c('khaki3','mediumorchid')[e],
          lty=c(2,1)[c])
  }
}
axis(side=1, at=seq(-1,1,0.5), pos=0.2, cex.axis=2.5, padj=0.5)
axis(side=2, at=seq(0.2, 0.8, 0.1), pos=-1.1, cex.axis=2.5, padj=0.5, las=1)
mtext('Identifiability', side=1, line=4, cex=3)
mtext('Detectability estimate', side=2, line=6, cex=3)
legend(x=-1.06, y=0.84, bty='n', lwd=3, cex=2.2,
       lty=c(1,2,1,2), col=rep(c('mediumorchid','khaki3'), each=2),
       legend=c('with experience, conspicuous species',
                'with experience, insconspicuous species',
                'no experience, conspicuous species',
                'no experience, inconspicuous species'))
dev.off()





#  Figure 4 - Scatter plot of detectability vs. occupancy  ####

# detectability
beta1 <- output$sims.list$beta1; dim(beta1)
beta2 <- output$sims.list$beta2; str(beta2)
beta3 <- output$sims.list$beta3; str(beta3)

# detectability per species*observer based on that species' conspicuousness & identifiability
#    and ASSUMING NO EXPERIENCE
detectability <- array(NA, dim=c(nspecies,nobservers,dim(beta1)[1]),
                       dimnames=list(dimnames(data$y)[[3]],
                                     paste('obs',1:7,sep=''),
                                     NULL)); dim(detectability)
for(k in 1:nspecies){
  for(o in 1:nobservers){
    detectability[k,o,] <- plogis(beta1[,o,k] + beta2*data$conspicuousness[k] + beta3*data$identifiability[k]) #  + beta4*1
  }
}
help <- apply(detectability, c(1,3), mean)
detectability.per.species <- data.frame(mean=apply(help, 1, mean),
                                        LCL=apply(help, 1, quantile, probs=0.025),
                                        UCL=apply(help, 1, quantile, probs=0.975))

# occupancy
out <- data.frame(output$sims.list)
params.titles <- colnames(out)
start.occ <- match("occ.pred.1", colnames(out))
end.occ <- match("occ.pred.373", colnames(out))


# dev.off()
tiff(paste('2_output/Fig.4_Detection_vs_occupancy.tif'), width=1700, height=1200, units='px', compression='lzw')
par(mfrow=c(1,1), mar=c(8,10,2,1))
plot(x=output$summary[start.occ:end.occ,'mean'], y=detectability.per.species$mean,
     xlab='', ylab='', axes=F,
     xlim=c(0,1), ylim=c(0.1,0.8), las=1, frame=F, pch=20, cex=2.6,
     col=alpha(c('khaki3','mediumorchid3'),0.3)[data$conspicuousness+1])
axis(side=1, at=seq(0, 1, 0.2), pos=0.1, cex.axis=2.5, padj=0.5)
axis(side=2, at=seq(0.1, 0.8, 0.1), pos=-0.04, cex.axis=2.5, padj=0.5, las=1)
mtext('Occupancy estimate', side=1, line=4, cex=3)
mtext('Detectability estimate', side=2, line=6, cex=3)
segments(x0=output$summary[start.occ:end.occ,3], x1=output$summary[start.occ:end.occ,7],
         y0=detectability.per.species$mean, # y1=y0 by default: y1=detectability.per.species$mean,
         col = alpha(c('khaki2','mediumorchid2'),c(0.5,0.2))[data$conspicuousness+1], lwd = 2)
segments(x0=output$summary[start.occ:end.occ,1], # x1=x0 by default output$summary[start.occ:end.occ,1],
         y0=detectability.per.species$LCL, y1=detectability.per.species$UCL,
         col = alpha(c('khaki2','mediumorchid2'),c(0.5,0.2))[data$conspicuousness+1], lwd = 2)
points(x=output$summary[start.occ:end.occ,'mean'], y=detectability.per.species$mean,
       pch=20, cex=2.6,
       col=alpha(c('khaki3','mediumorchid3'),c(0.6,0.4))[data$conspicuousness+1])
legend(x=0.7, y=0.2, xjust=0, yjust=0.5,pch=20, bty='n', cex=2.6, pt.cex=2.6,
       col=alpha(c('mediumorchid3','khaki3'),0.6),
       legend=c('conspicuous species','inconspicuous species'))

# arrows to individual species
selection <- c('Myriolecis persimilis','Graphis scripta','Pseudevernia furfuracea','Violella fucata','Lecanora barkmaniana')
colours <- c('goldenrod1','hotpink','black','forestgreen','blue')
for(i in 1:length(selection)){
  occ <- output$summary[start.occ:end.occ,'mean'][which(dimnames(data$y)[[3]]==selection[i])]
  det <- detectability.per.species$mean[which(dimnames(data$y)[[3]]==selection[i])]
  if(i==1){
    arrows(x0=occ-0.07,
           y0=det-0.07,
           x1=occ-0.01,
           y1=det-0.01,
           length=0.3, col=colours[i], lwd=3)
    text(x=occ-0.09,
         y=det-0.09,
         labels=expression(paste(italic('Mper'))),
         col=colours[i], cex=2.3)
  }
  if(i==2){
    arrows(x0=occ-0.07,
           y0=det+0.07,
           x1=occ-0.01,
           y1=det+0.01,
           length=0.3, col=colours[i], lwd=3)
    text(x=occ-0.09,
         y=det+0.09,
         labels=expression(paste(italic('Gscr'))),
         col=colours[i], cex=2.3)
  }
  if(i==3){
    arrows(x0=occ+0.07,
           y0=det+0.07,
           x1=occ+0.01,
           y1=det+0.01,
           length=0.3, col=colours[i], lwd=3)
    text(x=occ+0.09,
         y=det+0.09,
         labels=expression(paste(italic('Pfur'))),
         col=colours[i], cex=2.3)
  }
  if(i==4){
    arrows(x0=occ+0.07,
           y0=det-0.07,
           x1=occ+0.01,
           y1=det-0.01,
           length=0.3, col=colours[i], lwd=3)
    text(x=occ+0.09,
         y=det-0.09,
         labels=expression(paste(italic('Vfuc'))),
         col=colours[i], cex=2.3)
  }
  if(i==5){
    arrows(x0=occ,
           y0=det-0.1,
           x1=occ,
           y1=det-0.015,
           length=0.3, col=colours[i], lwd=3)
    text(x=occ,
         y=det-0.12,
         labels=expression(paste(italic('Lbar'))),
         col=colours[i], cex=2.3)
  }
}
dev.off()




#  Figure 5 - Barplot observed/estimated number of occurrences per species  ####

start <- match("N.est.1", colnames(out))
end <- match("N.est.416", colnames(out))
N.est <- apply(out[,start:end], 2, mean)

start <- match("n.occ.est.1", colnames(out))
end <- match("n.occ.est.373", colnames(out))
n.occ.est <- apply(out[,start:end], 2, mean)

N.obs <- apply(data$y[,1,], 1, sum)
n.occ.obs <- apply(data$y[,1,], 2, sum)

# number of species with n.occ.obs > 200
sum(n.occ.obs > 200)

# number of species with n.occ.est > 200
sum(n.occ.est > 200)

# dev.off()
tiff(paste('2_output/Fig.5_Occurrences_per_species.tif'), width=2000, height=1000, units='px', compression='lzw')
par(mfrow=c(1,1), mar=c(8,10,1,1))
plot(x=1:nspecies, y=sort(n.occ.obs, decreasing=T),
     type='h', lwd=3,
     xlab='', ylab='', axes=F,
     xlim=c(0, 400),
     ylim=c(0, max(n.occ.est)+20),
     frame=F, las=1)
axis(side=1, at=seq(0,400, 100), pos=-15, cex.axis=2.5, padj=0.5)
axis(side=2, at=seq(0, 300, 50), pos=-15, cex.axis=2.5, padj=0.5, las=1)
mtext('Species (ranked by observed frequency)', side=1, line=5, cex=3)
mtext('Number of occurrences', side=2, line=6, cex=3)
points(x=1:nspecies, y=n.occ.est[order(n.occ.obs, decreasing=T)],
       type='h', lwd=3, col='grey')
points(x=1:nspecies, y=sort(n.occ.obs, decreasing=T),
       type='h', lwd=3)
legend(x=320, y=300, xjust=1, yjust=1, lty=1, lwd=3, bty='n', cex=2.6,
       col=c('grey','black'),
       legend=c('estimated','observed'))

# add arrows to individual species
species.ordered <- dimnames(data$y)[[3]][order(n.occ.obs, decreasing=T)]
selection <- c('Myriolecis persimilis','Graphis scripta','Pseudevernia furfuracea','Violella fucata','Lecanora barkmaniana')
colours <- c('goldenrod1','hotpink','black','forestgreen','blue')
for(i in 1:length(selection)){
  arrows(x1=which(species.ordered==selection[i]),
         y1=n.occ.est[which(dimnames(data$y)[[3]]==selection[i])]+8,
         x0=which(species.ordered==selection[i]),
         y0=n.occ.est[which(dimnames(data$y)[[3]]==selection[i])]+80,
         length=0.3, col=colours[i], lwd=3)
  text(x=which(species.ordered==selection[i]), cex=2.6,
       y=n.occ.est[which(dimnames(data$y)[[3]]==selection[i])]+98,
       labels=c(expression(paste(italic(Mper))),
                expression(paste(italic(Gscr))),
                expression(paste(italic(Pfur))),
                expression(paste(italic(Vfuc))),
                expression(paste(italic(Lbar))))[i],
       col=colours[i])
}
dev.off()



#  Appendix 1 - Prior sensitivity analysis  ####

##  Illustrate priors themselves  ##

# priors 1
tiff('2_output/Appendix1_Fig.S1_Illustration_of_prior1.tif', width=1200, height=400, units='px')
par(mfrow=c(1,3))
x0 <- rcauchy(100000, scale=2.5)
hist(x0, main='prior used\nrcauchy(scale=2.5)',
     xlim=c(-30, 30),
     xlab='intercept',
     col=alpha('black',0.2),
     breaks=seq(min(x0)-2,(max(x0)+2), 0.8))

x1 <- rnorm(100000, mean=0, sd=10)
hist(x1, main='less informative prior\nrnorm(sd=10)',
     xlim=c(-30,30), xlab='intercept',
     col=alpha('goldenrod1',0.3),
     breaks=seq(min(x1)-2,(max(x1)+2), 0.8))

hist(x1, main='comparison', col=alpha('goldenrod1',0.3),
     xlim=c(-30, 30), ylim=c(0, 10000),
     xlab='intercept',
     breaks=seq(min(x0,x1)-2,(max(x0,x1)+2), 0.8))
hist(x0, add=T, col=alpha('black',0.2),
     breaks=seq(min(x0,x1)-2,(max(x0,x1)+2), 0.8))
dev.off()


# priors 2
tiff('2_output/Appendix1_Fig.S2_Illustration_of_prior2.tif', width=1200, height=400, units='px')
par(mfrow=c(1,3))
x0 <- rhalfcauchy(100000, scale=2.25)
hist(x0, main='prior used\nrcauchy(scale=2.25)',
     xlim=c(0, 30),
     xlab='intercept',
     col=alpha('black',0.2),
     breaks=seq(0,(max(x0)+2), 0.8))

x1 <- rhalfcauchy(100000, scale=25)
hist(x1, main='less informative prior\nrhalfcauchy(scale=25)',
     xlim=c(0,30), xlab='intercept',
     col=alpha('mediumorchid2',0.3),
     breaks=seq(0,(max(x1)+2), 0.8))

hist(x1, main='comparison', col=alpha('mediumorchid2',0.3),
     xlim=c(0, 30), ylim=c(0, 22000),
     xlab='intercept',
     breaks=seq(0,(max(x0,x1)+2), 0.8))
hist(x0, add=T, col=alpha('black',0.2),
     breaks=seq(0,(max(x0,x1)+2), 0.8))
dev.off()


##  Compare posteriors of main parameters  ##

load('2_output/output.jags_model.RData')
out0 <- as.data.frame(output$sims.list)[,1:32]
out0$prior <- 'used'

out0[,1:32] # are the main parameters
colnames(out0[,1:32]) # are their names

load('2_output/output.jags_model_priors1.RData')
out1 <- as.data.frame(output$sims.list)[,1:32]
out1$prior <- 'prior1'

load('2_output/output.jags_model_priors2.RData')
out2 <- as.data.frame(output$sims.list)[,1:32]
out2$prior <- 'prior2'

df <- rbind(out0, out1, out2)
params <- colnames(out0[,1:32])

# dev.off()
tiff('2_output/Appendix1_Fig.S3_Prior_sensitivity_assessment.tif', width=1900, height=1100, units='px', compression='lzw')
par(mfrow=c(4,8), mar=c(4,4,3,1))
for(i in 1:length(params)){
  boxplot(df[,params[i]] ~ df$prior, las=1, pch=20, cex.lab=1.2, cex.axis=2,
          col=alpha(c('darkgrey','goldenrod1','mediumorchid2'),0.6),
          main='', xlab='', ylab='', frame=F, xaxt='n')
  axis(1, at=1:3, col='white', col.ticks='white', labels=unique(df$prior), padj=-0.7, cex.axis=1.8)
  mtext(params[i], side=3, line=0, cex=1.8)
}
dev.off()

rm(out0, out1, out2)
load('2_output/output.jags_model.RData')




#  Appendix 2 - Observer*species-specific detectabilities  ####

beta1 <- output$sims.list$beta1; dim(beta1)
beta2 <- output$sims.list$beta2; str(beta2)

# detectability per species*observer based on that species conspicuousness
#    and ASSUMING NO EXPERIENCE and identifiability=0 (i.e. average identifiability)
detectability <- array(NA, dim=c(nspecies,nobservers,dim(beta1)[1]),
                       dimnames=list(dimnames(data$y)[[3]],
                                     paste('obs',1:7,sep=''),
                                     NULL)); dim(detectability)
for(k in 1:nspecies){
  for(o in 1:nobservers){
    detectability[k,o,] <- beta1[,o,k] + beta2*data$conspicuousness[k] #  + beta4*1
  }
}
detectability[1:10,1:3,1:2]

# DETECTABILITIES ARE SUMMARIZED ACROSS ALL OBSERVERS
mean.detectability <- apply(detectability, c(1,2), mean) %>% 
  plogis() %>% 
  round(3) %>% 
  as.data.frame()

# export
write.table(mean.detectability, file="2_output/Appendix2_Table.S1_Observer-species-specific-detectabilities.csv",sep=";", row.names=TRUE)

mean.detectability$mean <- round(apply(mean.detectability[,1:7], 1, mean),3)
mean.detectability$diff <- apply(mean.detectability[,1:7], 1, max) - apply(mean.detectability[,1:7], 1, min)




#  Appendix 3 - Species richness per site  ####

start <- match("N.est.1", colnames(out))
end <- match("N.est.416", colnames(out))
N.est <- apply(out[,start:end], 2, mean)

start <- match("n.occ.est.1", colnames(out))
end <- match("n.occ.est.373", colnames(out))
n.occ.est <- apply(out[,start:end], 2, mean)

N.obs <- apply(data$y[,1,], 1, sum)
n.occ.obs <- apply(data$y[,1,], 2, sum)

# Species richness per site  ##
# dev.off()
tiff(paste('2_output/Appendix3_Fig.S1_Species_richness_per_site.tif'), width=2000, height=1000, units='px', compression='lzw')
par(mfrow=c(1,1), mar=c(8,8,2,1))
plot(x=1:nsites, y=sort(N.obs, decreasing=T),
     type='h', lwd=3,
     xlab='', ylab='', axes=F,
     xlim=c(0, 420),
     ylim=c(0, max(N.est)),
     frame=F, las=1)
axis(side=1, at=seq(0, 400, 100), pos=-3, cex.axis=2.5, padj=0.5)
axis(side=2, at=seq(0, 80, 20), pos=-10, cex.axis=2.5, padj=0.5, las=1)
mtext('Sites (ranked by observed species richness)', side=1, line=5, cex=3)
mtext('Species richness per site', side=2, line=3.5, cex=3)
points(x=1:nsites, y=N.est[order(N.obs, decreasing=T)],
       type='h', lwd=3, col='grey')
points(x=1:nsites, y=sort(N.obs, decreasing=T),
       type='h', lwd=3)
legend(x=320, y=80, xjust=0, yjust=0.5, bty='n',lty=1, lwd=3, cex=2.6,
       col=c('grey','black'),
       legend=c('estimated','observed'))
dev.off()




#  Appendix 4 - Posterior predictive checks (Bayesian p-values)  ####

load('2_output/output.jags_model_gof.RData')
ifelse(sum(unlist(output$Rhat) > 1.1)==0, paste("Successful convergence based on Rhat values (all < 1.1)."), paste("Rhat values indicate convergence failure."))
# print(output, 2)

out <- data.frame(output$sims.list)
params.titles <- colnames(out)


# set threshold
# threshold <- 0.05 # often deemed to restrictive
threshold <- 0.1  # more conservative and reliable


# p.mean
start <- match("p.mean.1", colnames(out))
end <- match("p.mean.373", colnames(out))
p.means <- output$summary[start:end, "mean"]

par(mfrow=c(1,2), mar=c(3,5,5,1))
hist(p.means, breaks=51, xlim=c(0,1), las=1, cex.lab=1.2, cex.axis=1.1, xlab='',
     main="Bayesian p-value for \nmean number of detections per species")
abline(v=c(threshold,(1-threshold)), col="red", lwd=2)
text(x=0, y=2, labels=sum(p.means<=threshold), pos=3, cex=1.1, col="red")
text(x=1, y=2, labels=sum(p.means>=(1-threshold)), pos=3, cex=1.1, col="red")

# which species lie outside the acceptable range?
which(p.means<threshold|p.means>=(1-threshold))
dimnames(data$y)[[3]][which(p.means<=threshold|p.means>=(1-threshold))]

# how many times were these species observed?
n.occ.obs <- apply(data$y, 3, sum, na.rm=T)
# n.occ.obs[which(p.means<=threshold|p.means>=(1-threshold))]
# n.occ.obs[which(p.means<=threshold)]
enframe(n.occ.obs, name="var")[which(p.means<=threshold),]
enframe(n.occ.obs, name="var")[which(p.means>=(1-threshold)),]


# p.cv
start <- match("p.cv.1", colnames(out))
end <- match("p.cv.373", colnames(out))
p.cvs <- output$summary[start:end, "mean"]

hist(p.cvs, breaks=50, xlim=c(0,1), las=1, cex.lab=1.2, cex.axis=1.1, xlab='',
     main="Bayesian p-value for \nsd in number of detections per species")
abline(v=c(threshold,(1-threshold)), col="red", lwd=2)
text(x=0, y=2, labels=sum(p.cvs<=threshold), pos=3, cex=1.1, col="red")
text(x=1, y=2, labels=sum(p.cvs>=(1-threshold)), pos=3, cex=1.1, col="red")

# which species lie outside the acceptable range?
which(p.cvs<threshold|p.cvs>=(1-threshold))
dimnames(data$y)[[3]][which(p.cvs<=threshold|p.cvs>=(1-threshold))]
dimnames(data$y)[[3]][which(p.cvs<=threshold)]
dimnames(data$y)[[3]][which(p.cvs>=(1-threshold))]

n.occ.obs <- apply(data$y, 3, sum, na.rm=T)
# n.occ.obs[which(p.cvs<=threshold|p.cvs>=(1-threshold))]
# n.occ.obs[which(p.cvs<=threshold)]
enframe(n.occ.obs, name="var")[which(p.cvs<=threshold),]
enframe(n.occ.obs, name="var")[which(p.cvs>=(1-threshold)),]




#  Summary numbers for manuscript  ####

# RESULTS

# occupancy intercept mean = without suitable substrate
mean(plogis(out$mu.alpha0)) %>% round(3)
quantile(plogis(out$mu.alpha0), probs=c(0.025,0.975)) %>% round(3)

# average occupancy with suitable substrate
mean(plogis(out$mu.alpha0 + out$alpha1)) %>% round(3)
quantile(plogis(out$mu.alpha0 + out$alpha1), probs=c(0.025,0.975)) %>% round(3)

# observer-specific detectabilities as transformed from table 1
#   this is the intercept only: conspicuousness=0 and experience=0
output$summary[16:22,c(1,3,7)] %>% 
  plogis() %>% 
  round(3)

# species with the smallest difference between observers
#   taking species conspicuousness into account, but assuming experience=0
mean.detectability %>% filter(diff == min(mean.detectability$diff))
mean.detectability %>% filter(diff < 0.20) %>% arrange(diff)

# species with the greatest difference between observers
#   taking species conspicuousness into account, but assuming experience=0
mean.detectability %>% filter(diff == max(mean.detectability$diff))
mean.detectability %>% filter(diff > 0.66) %>% arrange(-diff)

# average detectabilities accounting for conspicuousness, identifiability & experience
load('0_data/experience_gain.RData')
dim(experience.gained)
experience.gained[1:10,,]

beta1 <- output$sims.list$beta1; dim(beta1)
beta2 <- output$sims.list$beta2; str(beta2)
beta3 <- output$sims.list$beta3; str(beta3)
beta4 <- output$sims.list$beta4; str(beta4)

p.over.time <- array(NA, dim=c(nspecies, nobservers, 2, length(beta2)),
                     dimnames=list(dimnames(data$y)[[3]],
                                   paste('obs',1:7,sep=''),
                                   c('start','end'),
                                   NULL)); dim(p.over.time)

for(o in 1:7){
  for(k in 1:373){
    p.over.time[k,o,'start',] <- plogis(beta1[,o,k]
                                        + beta2*data$conspicuousness[k]
                                        + beta3*data$identifiability[k]
                                        + beta4*experience.gained[k,o,'start'])
    p.over.time[k,o,'end',] <- plogis(beta1[,o,k]
                                      + beta2*data$conspicuousness[k]
                                      + beta3*data$identifiability[k]
                                      + beta4*experience.gained[k,o,'end'])
  }
}


# across the 'regular' first 6 observers
help1 <- apply(p.over.time[,1:6,,], c(1,3,4), mean); dim(help1)
help2 <- apply(help1, c(2,3), mean); dim(help2)

# at the start of the study (less experience)
mean(help2['start',]) %>% round(3)
quantile(help2['start',], probs=c(0.025, 0.5, 0.975)) %>% round(3)

# at the end of the study (more experience)
mean(help2['end',]) %>% round(3)
quantile(help2['end',], probs=c(0.025, 0.5, 0.975)) %>% round(3)

# across both periods
help3 <- apply(help2, 2, mean); str(help3)
mean(help3) %>% round(3)
quantile(help3, probs=c(0.025, 0.5, 0.975)) %>% round(3)



# ABSTRACT

# species with the lowest detectability across observers
mean.detectability %>% filter(mean==min(mean.detectability$mean)) %>% select(mean)

# species with the highest detectability across observers
mean.detectability %>% filter(mean==max(mean.detectability$mean)) %>% select(mean)

# average detection probability of conspicuous species
mean(help1[which(data$conspicuousness==1),,]) %>% round(3)
quantile(help1[which(data$conspicuousness==1),,], probs=c(0.025, 0.975)) %>% round(3)

# average detection probability of inconspicuous species
mean(help1[which(data$conspicuousness==0),,]) %>% round(3)
quantile(help1[which(data$conspicuousness==0),,], probs=c(0.025, 0.975)) %>% round(3)

# observer detectability means across species
apply(mean.detectability[,1:7], 2, mean) %>% round(3)

