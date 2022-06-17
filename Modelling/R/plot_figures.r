# June 2022 - plot all figures 

# Set directory to R data folder in repo. Then load data and libraries with:
library(ggplot2)
library(RColorBrewer)
load("Figures.RData")

# Note: The prediction examples needs to be run separately, not via sourcing entire file. 

##########################  
# Figure: Log Z
##########################  
pdf("logZ.pdf", width=11, height=8)
par(mfrow=c(2,2))
boxplot(as.numeric(Figures_L[,5]),as.numeric(Figures_L_4k[,5]),as.numeric(Figures_L_4theta[,5]),as.numeric(Figures_L_4k4theta[,5]),
        main= "Lysis data, n specimens",#main = "Lysate data, n specimens",
        at = c(1,2,3,4),
        names = c(expression(paste("1k 1",theta)),expression(paste("4k 1",theta)),expression(paste("1k 4",theta)),expression(paste("4k 4",theta))),     
        col=brewer.pal(n = 4, name = "Blues"),ylab="log Z",ylim=c(-1800,-1740)) 
title("(a)", adj = 0, line = 1.5)
boxplot(as.numeric(Figures_L_w[,5]),as.numeric(Figures_L_4k_w[,5]),as.numeric(Figures_L_4theta_w[,5]),as.numeric(Figures_L_4k4theta_w[,5]),
        main= "Lysis data, weights",#main = "Lysate data, weights",
        at = c(1,2,3,4),
        names = c(expression(paste("1k 1",theta)),expression(paste("4k 1",theta)),expression(paste("1k 4",theta)),expression(paste("4k 4",theta))),
        col=brewer.pal(n = 4, name = "Blues"),ylab="log Z",ylim=c(-1800,-1740)) 
title("(b)", adj = 0, line = 1.5)
boxplot(as.numeric(Figures_H[,5]),as.numeric(Figures_H_4k[,5]),as.numeric(Figures_H_4theta[,5]),as.numeric(Figures_H_4k4theta[,5]),
        main= "Homogenate data, n specimens",# main = "Homogenate data, n specimens",
        at = c(1,2,3,4),
        names = c(expression(paste("1k 1",theta)),expression(paste("4k 1",theta)),expression(paste("1k 4",theta)),expression(paste("4k 4",theta))),
        col=brewer.pal(n = 4, name = "Blues"),ylab="log Z",ylim=c(-669,-656)) 
title("(c)", adj = 0, line = 1.5)
boxplot(as.numeric(Figures_H_w[,5]),as.numeric(Figures_H_4k_w[,5]),as.numeric(Figures_H_4theta_w[,5]),as.numeric(Figures_H_4k4theta_w[,5]),
        main= "Homogenate data, weights",#main = "Homogenate data, weights",
        at = c(1,2,3,4),
        names = c(expression(paste("1k 1",theta)),expression(paste("4k 1",theta)),expression(paste("1k 4",theta)),expression(paste("4k 4",theta))),
        col=brewer.pal(n = 4, name = "Blues"),ylab="log Z",ylim=c(-669,-656))
title("(d)", adj = 0, line = 1.5)
dev.off()

############################## 
# Figure: d and CV - n
############################## 
pdf("Densitites_d_CV_n.pdf", width=11, height=8)
par(mfrow=c(2,2))
#Lysate
plot(density(as.numeric(Figures_L[,6])*as.numeric(Figures_L[,7]),weights=as.numeric(as.numeric(Figures_L[,5]))/sum(as.numeric(Figures_L[,5])),na.rm=T,bw=0.75),type='l',bty='n',main="Lysis data",col="Blue",xlim=c(0,20),ylim=c(0,0.65),xlab="Mean yield (d)",lwd = 1.8)
points(density(as.numeric(Figures_L_4k_k[,4]*as.numeric(Figures_L_4k[,7])),weights=as.numeric(as.numeric(Figures_L_4k[,5]))/sum(as.numeric(Figures_L_4k[,5])),na.rm=T,bw=0.5),type='l',bty='n',main="4 k * theta",col="Black",xlim=c(0,20),ylim=c(0,0.65),xlab="",lwd = 1.8)
points(density(as.numeric(Figures_L_4k_k[,1])*as.numeric(Figures_L_4k[,7]),weights=as.numeric(as.numeric(Figures_L_4k[,5]))/sum(as.numeric(Figures_L_4k[,5])),na.rm=T,bw=0.75),type='l',bty='n',col="Dark green",lwd = 1.8)
points(density(as.numeric(Figures_L_4k_k[,2]*as.numeric(Figures_L_4k[,7])),weights=as.numeric(as.numeric(Figures_L_4k[,5]))/sum(as.numeric(Figures_L_4k[,5])),na.rm=T,bw=0.75),type='l',bty='n',col="Dark red",lwd = 1.8)
points(density(as.numeric(Figures_L_4k_k[,3]*as.numeric(Figures_L_4k[,7])),weights=as.numeric(as.numeric(Figures_L_4k[,5]))/sum(as.numeric(Figures_L_4k[,5])),na.rm=T,bw=0.75),type='l',bty='n',col="Dark orange",lwd = 1.8)
legend("topright", legend=c("Simple model","S. lateralis","G. bimaculatus","G. sigillatus","D. yakuba"),col = c("Blue","Dark green","Dark red","Dark orange","Black"),pch=15,cex=1,y.intersp=1.0,bty = "n")
title("(a)", adj = 0, line = 1)
#Homogenate
plot(density(as.numeric(Figures_H[,6])*as.numeric(Figures_H[,7]),weights=as.numeric(as.numeric(Figures_H[,5]))/sum(as.numeric(Figures_H[,5])),na.rm=T,bw=0.8),type='l',bty='n',main="Homogenate data",col="Blue",xlim=c(0,20),ylim=c(0,0.25),xlab="Mean yield (d)",lwd = 1.8)
points(density(as.numeric(Figures_H_4k_k[,4]*as.numeric(Figures_H_4k[,7])),weights=as.numeric(as.numeric(Figures_H_4k[,5]))/sum(as.numeric(Figures_H_4k[,5])),na.rm=T,bw=0.8),type='l',bty='n',main="4 k * theta",col="Black",xlim=c(0,20),ylim=c(0,0.25),xlab="",lwd = 1.8)
points(density(as.numeric(Figures_H_4k_k[,1])*as.numeric(Figures_H_4k[,7]),weights=as.numeric(as.numeric(Figures_H_4k[,5]))/sum(as.numeric(Figures_H_4k[,5])),na.rm=T,bw=0.8),type='l',bty='n',col="Dark green",lwd = 1.8)
points(density(as.numeric(Figures_H_4k_k[,2]*as.numeric(Figures_H_4k[,7])),weights=as.numeric(as.numeric(Figures_H_4k[,5]))/sum(as.numeric(Figures_H_4k[,5])),na.rm=T,bw=0.8),type='l',bty='n',col="Dark red",lwd = 1.8)
points(density(as.numeric(Figures_H_4k_k[,3]*as.numeric(Figures_H_4k[,7])),weights=as.numeric(as.numeric(Figures_H_4k[,5]))/sum(as.numeric(Figures_H_4k[,5])),na.rm=T,bw=0.8),type='l',bty='n',col="Dark orange",lwd = 1.8)
legend("topright", legend=c("Simple model","S. lateralis","G. bimaculatus","G. sigillatus","D. yakuba"),col = c("Blue","Dark green","Dark red","Dark orange","Black"),pch=15,cex=1,y.intersp=1.0,bty = "n")
title("(b)", adj = 0, line = 1)
#Lysate
plot(density(1/sqrt(as.numeric(Figures_L[,6])),weights=as.numeric(as.numeric(Figures_L[,5]))/sum(as.numeric(Figures_L[,5])),na.rm=T),type='l',bty='n',xlab="CV",col="Blue",xlim=c(0,1.5),ylim=c(0,20),lwd = 1.8,main="Lysis data")
lines(density(1/sqrt(as.numeric(Figures_L_4k_k[,4])),weights=as.numeric(as.numeric(Figures_L_4k[,5]))/sum(as.numeric(Figures_L_4k[,5])),na.rm=T),type='l',bty='n',col="Black",xlim=c(0,1.5),ylim=c(0,20),xlab="k",lwd = 1.8)
lines(density(1/sqrt(as.numeric(Figures_L_4k_k[,1])),weights=as.numeric(as.numeric(Figures_L_4k[,5]))/sum(as.numeric(Figures_L_4k[,5])),na.rm=T),type='l',bty='n',col="Dark green",lwd = 1.8)
lines(density(1/sqrt(as.numeric(Figures_L_4k_k[,2])),weights=as.numeric(as.numeric(Figures_L_4k[,5]))/sum(as.numeric(Figures_L_4k[,5])),na.rm=T),type='l',bty='n',col="Dark red",lwd = 1.8)
lines(density(1/sqrt(as.numeric(Figures_L_4k_k[,3])),weights=as.numeric(as.numeric(Figures_L_4k[,5]))/sum(as.numeric(Figures_L_4k[,5])),na.rm=T),type='l',bty='n',col="Dark orange",lwd = 1.8)
legend("topright", legend=c("Simple model","S. lateralis","G. bimaculatus","G. sigillatus","D. yakuba"),col = c("Blue","Dark green","Dark red","Dark orange","Black"),pch=15,cex=1,y.intersp=1.0,bty = "n")
title("(c)", adj = 0, line = 1.5)
#Homogenate
plot(density(1/sqrt(as.numeric(Figures_H[,6])),weights=as.numeric(as.numeric(Figures_H[,5]))/sum(as.numeric(Figures_H[,5])),na.rm=T),type='l',bty='n',xlab="CV",col="Blue",xlim=c(0,1.5),ylim=c(0,20),lwd = 1.8,main="Homogenate data")
lines(density(1/sqrt(as.numeric(Figures_H_4k_k[,4])),weights=as.numeric(as.numeric(Figures_H_4k[,5]))/sum(as.numeric(Figures_H_4k[,5])),na.rm=T),type='l',bty='n',col="Black",xlim=c(0,1.5),ylim=c(0,20),xlab="k",lwd = 1.8)
lines(density(1/sqrt(as.numeric(Figures_H_4k_k[,1])),weights=as.numeric(as.numeric(Figures_H_4k[,5]))/sum(as.numeric(Figures_H_4k[,5])),na.rm=T),type='l',bty='n',col="Dark green",lwd = 1.8)
lines(density(1/sqrt(as.numeric(Figures_H_4k_k[,2])),weights=as.numeric(as.numeric(Figures_H_4k[,5]))/sum(as.numeric(Figures_H_4k[,5])),na.rm=T),type='l',bty='n',col="Dark red",lwd = 1.8)
lines(density(1/sqrt(as.numeric(Figures_H_4k_k[,3])),weights=as.numeric(as.numeric(Figures_H_4k[,5]))/sum(as.numeric(Figures_H_4k[,5])),na.rm=T),type='l',bty='n',col="Dark orange",lwd = 1.8)
legend("topright", legend=c("Simple model","S. lateralis","G. bimaculatus","G. sigillatus","D. yakuba"),col = c("Blue","Dark green","Dark red","Dark orange","Black"),pch=15,cex=1,y.intersp=1.0,bty = "n")
title("(d)", adj = 0, line = 1.5)
dev.off()

##################################
# Supplemental figure: d and CV - weights
##################################
#Lysate-w
pdf("Densitites_d_CV_w.pdf", width=11, height=8)
par(mfrow=c(2,2))
plot(density(as.numeric(Figures_L_w[,6])*as.numeric(Figures_L_w[,7]),weights=as.numeric(as.numeric(Figures_L_w[,5]))/sum(as.numeric(Figures_L_w[,5])),na.rm=T),type='l',bty='n',main="Lysis data (weights)",col="Blue",xlim=c(0,60),ylim=c(0,0.65),xlab="Mean yield (d)",lwd = 1.8)
points(density(as.numeric(Figures_L_4k_w_k[,4]*as.numeric(Figures_L_4k_w[,7])),weights=as.numeric(as.numeric(Figures_L_4k_w[,5]))/sum(as.numeric(Figures_L_4k_w[,5])),na.rm=T),type='l',bty='n',col="Black",xlim=c(0,60),ylim=c(0,0.65),xlab="",lwd = 1.8)
points(density(as.numeric(Figures_L_4k_w_k[,1])*as.numeric(Figures_L_4k_w[,7]),weights=as.numeric(as.numeric(Figures_L_4k_w[,5]))/sum(as.numeric(Figures_L_4k_w[,5])),na.rm=T),type='l',bty='n',col="Dark green",lwd = 1.8)
points(density(as.numeric(Figures_L_4k_w_k[,2]*as.numeric(Figures_L_4k_w[,7])),weights=as.numeric(as.numeric(Figures_L_4k_w[,5]))/sum(as.numeric(Figures_L_4k_w[,5])),na.rm=T),type='l',bty='n',col="Dark red",lwd = 1.8)
points(density(as.numeric(Figures_L_4k_w_k[,3]*as.numeric(Figures_L_4k_w[,7])),weights=as.numeric(as.numeric(Figures_L_4k_w[,5]))/sum(as.numeric(Figures_L_4k_w[,5])),na.rm=T),type='l',bty='n',col="Dark orange",lwd = 1.8)
legend("topright", legend=c("Simple model","S. lateralis","G. bimaculatus","G. sigillatus","D. yakuba"),col = c("Blue","Dark green","Dark red","Dark orange","Black"),pch=15,cex=1,y.intersp=1.0,bty = "n")
title("(a)", adj = 0, line = 1)
#Homogenate-w
plot(density(as.numeric(Figures_H_w[,6])*as.numeric(Figures_H_w[,7]),weights=as.numeric(as.numeric(Figures_H_w[,5]))/sum(as.numeric(Figures_H_w[,5])),na.rm=T),type='l',bty='n',main="Homogenate data (weights)",col="Blue",xlim=c(0,60),ylim=c(0,0.25),xlab="Mean yield (d)",lwd = 1.8)
points(density(as.numeric(Figures_H_4k_w_k[,4]*as.numeric(Figures_H_4k_w[,7])),weights=as.numeric(as.numeric(Figures_H_4k_w[,5]))/sum(as.numeric(Figures_H_4k_w[,5])),na.rm=T),type='l',bty='n',col="Black",xlim=c(0,60),ylim=c(0,0.25),xlab="",lwd = 1.8)
points(density(as.numeric(Figures_H_4k_w_k[,1])*as.numeric(Figures_H_4k_w[,7]),weights=as.numeric(as.numeric(Figures_H_4k_w[,5]))/sum(as.numeric(Figures_H_4k_w[,5])),na.rm=T),type='l',bty='n',col="Dark green",lwd = 1.8)
points(density(as.numeric(Figures_H_4k_w_k[,2]*as.numeric(Figures_H_4k_w[,7])),weights=as.numeric(as.numeric(Figures_H_4k_w[,5]))/sum(as.numeric(Figures_H_4k_w[,5])),na.rm=T),type='l',bty='n',col="Dark red",lwd = 1.8)
points(density(as.numeric(Figures_H_4k_w_k[,3]*as.numeric(Figures_H_4k_w[,7])),weights=as.numeric(as.numeric(Figures_H_4k_w[,5]))/sum(as.numeric(Figures_H_4k_w[,5])),na.rm=T),type='l',bty='n',col="Dark orange",lwd = 1.8)
legend("topright", legend=c("Simple model","S. lateralis","G. bimaculatus","G. sigillatus","D. yakuba"),col = c("Blue","Dark green","Dark red","Dark orange","Black"),pch=15,cex=1,y.intersp=1.0,bty = "n")
title("(b)", adj = 0, line = 1)
#Lysate-w
plot(density(1/sqrt(as.numeric(Figures_L_w[,6])),weights=as.numeric(as.numeric(Figures_L_w[,5]))/sum(as.numeric(Figures_L_w[,5])),na.rm=T),type='l',bty='n',xlab="CV",col="Blue",xlim=c(0,1.5),ylim=c(0,20),lwd = 1.8,main="Lysis data (weights)")
lines(density(1/sqrt(as.numeric(Figures_L_4k_w_k[,4])),weights=as.numeric(as.numeric(Figures_L_4k_w[,5]))/sum(as.numeric(Figures_L_4k_w[,5])),na.rm=T),type='l',bty='n',col="Black",xlim=c(0,1.5),ylim=c(0,20),lwd = 1.8)
lines(density(1/sqrt(as.numeric(Figures_L_4k_w_k[,1])),weights=as.numeric(as.numeric(Figures_L_4k_w[,5]))/sum(as.numeric(Figures_L_4k_w[,5])),na.rm=T),type='l',bty='n',col="Dark green",lwd = 1.8)
lines(density(1/sqrt(as.numeric(Figures_L_4k_w_k[,2])),weights=as.numeric(as.numeric(Figures_L_4k_w[,5]))/sum(as.numeric(Figures_L_4k_w[,5])),na.rm=T),type='l',bty='n',col="Dark red",lwd = 1.8)
lines(density(1/sqrt(as.numeric(Figures_L_4k_w_k[,3])),weights=as.numeric(as.numeric(Figures_L_4k_w[,5]))/sum(as.numeric(Figures_L_4k_w[,5])),na.rm=T),type='l',bty='n',col="Dark orange",lwd = 1.8)
legend("topright", legend=c("Simple model","S. lateralis","G. bimaculatus","G. sigillatus","D. yakuba"),col = c("Blue","Dark green","Dark red","Dark orange","Black"),pch=15,cex=1,y.intersp=1.0,bty = "n")
title("(c)", adj = 0, line = 1.5)
#Homogenate-w
plot(density(1/sqrt(as.numeric(Figures_H_w[,6])),weights=as.numeric(as.numeric(Figures_H_w[,5]))/sum(as.numeric(Figures_H_w[,5])),na.rm=T),type='l',bty='n',xlab="CV",col="Blue",xlim=c(0,1.5),ylim=c(0,20),lwd = 1.8,main="Homogenate data (weights)")
lines(density(1/sqrt(as.numeric(Figures_H_4k_w_k[,4])),weights=as.numeric(as.numeric(Figures_H_4k_w[,5]))/sum(as.numeric(Figures_H_4k_w[,5])),na.rm=T),type='l',bty='n',col="Black",xlim=c(0,1.5),ylim=c(0,20),lwd = 1.8)
lines(density(1/sqrt(as.numeric(Figures_H_4k_w_k[,1])),weights=as.numeric(as.numeric(Figures_H_4k_w[,5]))/sum(as.numeric(Figures_H_4k_w[,5])),na.rm=T),type='l',bty='n',col="Dark green",lwd = 1.8)
lines(density(1/sqrt(as.numeric(Figures_H_4k_w_k[,2])),weights=as.numeric(as.numeric(Figures_H_4k_w[,5]))/sum(as.numeric(Figures_H_4k_w[,5])),na.rm=T),type='l',bty='n',col="Dark red",lwd = 1.8)
lines(density(1/sqrt(as.numeric(Figures_H_4k_w_k[,3])),weights=as.numeric(as.numeric(Figures_H_4k_w[,5]))/sum(as.numeric(Figures_H_4k_w[,5])),na.rm=T),type='l',bty='n',col="Dark orange",lwd = 1.8)
legend("topright", legend=c("Simple model","S. lateralis","G. bimaculatus","G. sigillatus","D. yakuba"),col = c("Blue","Dark green","Dark red","Dark orange","Black"),pch=15,cex=1,y.intersp=1.0,bty = "n")
title("(d)", adj = 0, line = 1.5)
dev.off()

####################################
# Figure: Predictions - n
####################################
#Plot n Predictions lysate and homogenate
pdf("Predictions_hist_n.pdf", width=6, height=8)
par(mfrow=c(2,1))
foo<-hist(Figures2_LPred_n[,4]-Figures2_LPred_n[,3],breaks = -11:12,xaxt="n",xlab="MAP - true n",main="",col="lightskyblue")
axis(side=1,at=foo$mids,labels=seq(-10,12))
text(x = -6, y = 60,                
     "Underestimate",
     col = "ivory4",
     cex = 1)
text(x = 4, y = 60,                
     "Overestimate",
     col = "ivory4",
     cex = 1)
title("(a) Predictions for lysis data", adj = 0, line = 1.5)
foo<-hist(Figures2_HPred_n[,4]-Figures2_HPred_n[,3],breaks = -11:12,xaxt="n",xlab="MAP - true n",main="",col="lightsteelblue1")
axis(side=1,at=foo$mids,labels=seq(-10,12))
text(x = -6, y = 12,                
     "Underestimate",
     col = "ivory4",
     cex = 1)
text(x = 4, y = 12,                
     "Overestimate",
     col = "ivory4",
     cex = 1)
title("(c) Predictions for homogenate data", adj = 0, line = 1.5)
dev.off()

####################################
# Supplemental figure: Predictions - weights
####################################
pdf("Predictions_hist_w.pdf", width=6, height=8)
par(mfrow=c(2,1))
hist(log(Figures2_LPred_w[,6]),xlab="MAP / true biomass - log scale",main="",xlim=c(-10,10),breaks = -11:10,col="lightskyblue")
title("(a) Predictions for lysis data (weights)", adj = 0, line = 1.5)
text(x = -4, y = 40,                
     "Underestimate",
     col = "ivory4",
     cex = 1)
text(x = 7, y = 40,                
     "Overestimate",
     col = "ivory4",
     cex = 1)
hist(log(Figures2_HPred_w[,6]),xlab="MAP / true biomass - log scale",main="",xlim=c(-10,10),breaks = -11:10,col="lightsteelblue1")
title("(c) Predictions for homogenate data (weights)", adj = 0, line = 1.5)
text(x = -6, y = 12,                
     "Underestimate",
     col = "ivory4",
     cex = 1)
text(x = 5, y = 12,                
     "Overestimate",
     col = "ivory4",
     cex = 1)
dev.off()

####################################
# Figure: Predictions examples - n
####################################
#Homogenate - n
pdf("Predictions_examples_H_n.pdf", width=5, height=4)
suppressWarnings(ggplot(df1, aes(x=Species, y=Mean)) + 
  scale_y_continuous(limits = c(0, 10), breaks = c(0, 2, 4, 6, 8, 10))+
  ylab("Number of specimens (n)") +
  geom_point(size = 2,aes(colour="Predicted n"))+
  geom_errorbar(aes(ymin=Mean-sd, ymax=Mean+sd), width=0,colour="black",
                position=position_dodge(0.05))+
  geom_point(aes(x=1,y=1,colour="True n"),size = 3,alpha=0.05)+
  geom_point(aes(x=2,y=1,colour="True n"),size = 3,alpha=0.05)+
  geom_point(aes(x=3,y=1,colour="True n"),size = 3,alpha=0.05)+
  geom_point(aes(x=5,y=5,colour="True n"),size = 3,alpha=0.05)+
  geom_point(aes(x=6,y=8,colour="True n"),size = 3,alpha=0.05)+
  geom_point(aes(x=7,y=8,colour="True n"),size = 3,alpha=0.05)+
  geom_point(aes(x=9,y=1,colour="True n"),size = 3,alpha=0.05)+
  geom_point(aes(x=10,y=1,colour="True n"),size = 3,alpha=0.05)+
  geom_point(aes(x=11,y=1,colour="True n"),size = 3,alpha=0.05)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = a))+
  scale_color_manual(name = "",
                     breaks = c("Predicted n", "True n"),
                     values = c("Predicted n" = "black", "True n" = adjustcolor( "red", alpha.f = 0.2)) )+
  labs(title="(d) Species examples, homogenate")+
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  theme(plot.title = element_text(face="bold"))+
  theme(plot.title = element_text(hjust = 0.1))+
  theme(plot.title = element_text(vjust = 5)) )
dev.off()
################################
#Lysate - n
pdf("Predictions_examples_L_n.pdf", width=5, height=4)
suppressWarnings(ggplot(df2, aes(x=Species, y=Mean)) + 
  scale_y_continuous(limits = c(0, 10), breaks = c(0, 2, 4, 6, 8, 10))+
  ylab("Number of specimens (n)") +
  geom_point(size = 2,aes(colour="Predicted n"))+
  geom_errorbar(aes(ymin=Mean-sd, ymax=Mean+sd), width=0,colour="black",
                position=position_dodge(0.05))+
  geom_point(aes(x=1,y=1,colour="True n"),size = 3,alpha=0.05)+
  geom_point(aes(x=2,y=1,colour="True n"),size = 3,alpha=0.05)+
  geom_point(aes(x=3,y=2,colour="True n"),size = 3,alpha=0.05)+
  geom_point(aes(x=5,y=2,colour="True n"),size = 3,alpha=0.05)+
  geom_point(aes(x=6,y=2,colour="True n"),size = 3,alpha=0.05)+
  geom_point(aes(x=7,y=9,colour="True n"),size = 3,alpha=0.05)+
  geom_point(aes(x=9,y=1,colour="True n"),size = 3,alpha=0.05)+
  geom_point(aes(x=10,y=1,colour="True n"),size = 3,alpha=0.05)+
  geom_point(aes(x=11,y=1,colour="True n"),size = 3,alpha=0.05)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = a))+
  scale_color_manual(name = "",
                     breaks = c("Predicted n", "True n"),
                     values = c("Predicted n" = "black", "True n" = adjustcolor( "red", alpha.f = 0.2)) )+
  labs(title="(b) Species examples, lysis")+
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  theme(plot.title = element_text(face="bold"))+
  theme(plot.title = element_text(hjust = 0.1))+
  theme(plot.title = element_text(vjust = 5)) )
dev.off()
####################################
# Supplemental figure: Predictions examples - weights
####################################
#Homogenate - w
pdf("Predictions_examples_H_w.pdf", width=5, height=4)
suppressWarnings(ggplot(df3, aes(x=Species, y=Mean)) + 
  scale_y_continuous(limits = c(0, 40))+#, breaks = c(0, 2, 4, 6, 8, 10))+
  ylab("Biomass (mg)") +
  geom_point(size = 2,aes(colour="Predicted mg"))+
  geom_errorbar(aes(ymin=Mean-sd, ymax=Mean+sd), width=0,colour="black",
                position=position_dodge(0.05))+
  geom_point(aes(x=1,y=3.65,colour="True mg"),size = 3,alpha=0.05)+
  geom_point(aes(x=2,y=3.65,colour="True mg"),size = 3,alpha=0.05)+
  geom_point(aes(x=3,y=3.65,colour="True mg"),size = 3,alpha=0.05)+
  geom_point(aes(x=5,y=38.48,colour="True mg"),size = 3,alpha=0.05)+
  geom_point(aes(x=6,y=19.24,colour="True mg"),size = 3,alpha=0.05)+
  geom_point(aes(x=7,y=19.24,colour="True mg"),size = 3,alpha=0.05)+
  geom_point(aes(x=9,y=0.98,colour="True mg"),size = 3,alpha=0.05)+
  geom_point(aes(x=10,y=0.98,colour="True mg"),size = 3,alpha=0.05)+
  geom_point(aes(x=11,y=0.98,colour="True mg"),size = 3,alpha=0.05)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = a))+
  scale_color_manual(name = "",
                     breaks = c("Predicted mg", "True mg"),
                     values = c("Predicted mg" = "black", "True mg" = adjustcolor( "red", alpha.f = 0.2)) )+
  labs(title="(d) Species examples, homogenate")+
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  theme(plot.title = element_text(face="bold"))+
  theme(plot.title = element_text(hjust = 0.1))+
  theme(plot.title = element_text(vjust = 5)) )
dev.off()
################################
#Lysate - w 
pdf("Predictions_examples_L_w.pdf", width=5, height=4)
suppressWarnings(ggplot(df4, aes(x=factor(Species, level = level_order), y=Mean)) + 
  scale_y_continuous(limits = c(-0.5, 8), breaks = c(0, 2, 4, 6, 8))+
  ylab("Biomass (mg)") +
  xlab("Species") +
  geom_point(size = 2,aes(colour="Predicted mg"))+
  geom_errorbar(aes(ymin=Mean-sd, ymax=Mean+sd), width=0,colour="black",
                position=position_dodge(0.05))+
  geom_point(aes(x=1,y=0.17,colour="True mg"),size = 3,alpha=0.05)+
  geom_point(aes(x=2,y=0.34,colour="True mg"),size = 3,alpha=0.05)+
  geom_point(aes(x=3,y=0.34,colour="True mg"),size = 3,alpha=0.05)+
  geom_point(aes(x=5,y=3.65,colour="True mg"),size = 3,alpha=0.05)+
  geom_point(aes(x=6,y=3.65,colour="True mg"),size = 3,alpha=0.05)+
  geom_point(aes(x=7,y=7.3,colour="True mg"),size = 3,alpha=0.05)+
  geom_point(aes(x=9,y=0.003,colour="True mg"),size = 3,alpha=0.05)+
  geom_point(aes(x=10,y=0.0015,colour="True mg"),size = 3,alpha=0.05)+
  geom_point(aes(x=11,y=0.18,colour="True mg"),size = 3,alpha=0.05)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = a))+
  scale_color_manual(name = "",
                     breaks = c("Predicted mg", "True mg"),
                     values = c("Predicted mg" = "black", "True mg" = adjustcolor( "red", alpha.f = 0.2)) )+
  labs(title="(b) Species examples, lysis")+
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  theme(plot.title = element_text(face="bold"))+
  theme(plot.title = element_text(hjust = 0.1))+
  theme(plot.title = element_text(vjust = 5)) )
dev.off()
################################################################


