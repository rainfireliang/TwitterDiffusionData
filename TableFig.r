library(lme4)
library(sjPlot)
library(sjstats)
library(eoffice)
library(tidyr)
library(dplyr)

#######
load("data_model.Rdata")

fm1 = update(fm_t[[1]],.~.+p_inter_s+p_nexpose_s+p_sim_s+t_inter_s+t_nexpose_s+t_sim_s,data=data)
fm2 = update(fm_t[[2]],.~.+p_inter_s+p_nexpose_s+p_sim_s+t_inter_s+t_nexpose_s+t_sim_s,data=data)
fm3 = update(fm_t[[3]],.~.+p_inter_s+p_nexpose_s+p_sim_s+t_inter_s+t_nexpose_s+t_sim_s,data=data)

fm_igm = list(fm1,fm2,fm3)
save(fm_igm,file="fm_igm.Rdata")
tab_model(fm_igm[[1]],fm_igm[[2]],fm_igm[[3]],show.est = T,show.ci = F,show.icc = F,show.se = T,show.aic = T,transform = NULL,digits=3)

####### Interaction Plot
load("fm_igm.Rdata")
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

library(interplot) #https://cran.r-project.org/web/packages/interplot/vignettes/interplot-vignette.html
A = interplot(fm_igm[[3]],var1 = 'inter_s',var2 = 'depth_s',xmin=-1.253899,xmax=0.8868936) + 
  labs(x = "Cascade depth", y = "Contagion (exposure) effect")
B = interplot(fm_igm[[3]],var1 = 'nexpose_s',var2 = 'depth_s',xmin=-1.253899,xmax=0.8868936) + 
  labs(x = "Cascade depth", y = "Contagion (exposure) effect") 
C = interplot(fm_igm[[3]],var1 = 'sim_s',var2 = 'depth_s',xmin=-1.253899,xmax=0.8868936) + 
  labs(x = "Cascade depth", y = "Interest similarity effect") 


a = interplot(fm_igm[[3]],var1 = 'inter_s',var2 = 'depth_s',xmin=mean(data$depth_s[data$depth==1]),xmax=mean(data$depth_s[data$depth==7])) + 
  theme_bw()+geom_hline(yintercept = 0, linetype = "dashed")+
  scale_x_continuous(breaks=sapply(1:7,function(i) mean(data$depth_s[data$depth==i])),
                   labels=paste(c(1:7)))+
  scale_y_continuous(breaks=c(0,1.2558597,0.7087153),# 1.2558597 [1.1956456,1.3172100];0.7087153 [0.6481587,0.7700426]
                     labels=c("0.00","1.25","0.71"))+
  labs(x = "Cascade depth", y = "Contagion (interaction) effect") 

b = interplot(fm_igm[[3]],var1 = 'nexpose_s',var2 = 'depth_s',xmin=mean(data$depth_s[data$depth==1]),xmax=mean(data$depth_s[data$depth==7])) + 
  theme_bw()+geom_hline(yintercept = 0, linetype = "dashed")+
  scale_x_continuous(breaks=sapply(1:7,function(i) mean(data$depth_s[data$depth==i])),
                     labels=paste(c(1:7)))+
  scale_y_continuous(breaks=c(0,-2.4970234,-0.6212904),# -2.4970234 [-2.6347713,-2.3614676];-0.6212904 [-0.7189364,-0.5262261]
                     labels=c("0.00","-2.50","-0.62"))+
  labs(x = "Cascade depth", y = "Contagion (exposure) effect")

c = interplot(fm_igm[[3]],var1 = 'sim_s',var2 = 'depth_s',xmin=mean(data$depth_s[data$depth==1]),xmax=mean(data$depth_s[data$depth==7])) + 
  theme_bw()+geom_hline(yintercept = 0, linetype = "dashed")+
  scale_x_continuous(breaks=sapply(1:7,function(i) mean(data$depth_s[data$depth==i])),
                     labels=paste(c(1:7)))+
  scale_y_continuous(breaks=c(0,0.5188291,0.8658957),# 0.5188291 [0.4365918,0.5995283];0.8658957 [0.7891052,0.9429618]
                     labels=c("0.00","0.52","0.87"))+
  labs(x = "Cascade depth", y = "Interest similarity effect")

##
x = get_model_data(fm_igm[[3]],type = "pred",terms=c('inter_s[-0.5867,0.3199]',"depth_s[-1.253899,0.8868936]"))
y = get_model_data(fm_igm[[3]],type = "pred",terms=c('nexpose_s[-0.6317,0.145]',"depth_s[-1.253899,0.8868936]"))
z = get_model_data(fm_igm[[3]],type = "pred",terms=c('sim_s[-0.4371,2.787166]',"depth_s[-1.253899,0.8868936]"),ci.lvl = NA)

e = plot_model(fm_igm[[3]],type = "pred",
               terms=c("inter_s[-0.5867,0.3199]","depth_s[-1.253899,0.8868936]"),
               title="",axis.title = "",axis.labels = "",
               legend.title = "Depth",colors = "bw")+
  scale_x_continuous(breaks=c(-0.5867,0.3199),#0,2
                     labels=c("0","2"))+
  scale_y_continuous(breaks=x$predicted,
                     labels = c("1.9%","1.5%","5.7%","2.8%"))+
  labs(x = "Interaction frequency", y = "Predicted probability of\nretweeting")+
  theme_bw() + theme(legend.position = c(0.3, 0.8), legend.direction = "horizontal",
                     plot.title = element_blank(),plot.margin = unit(c(5.5,5.5,12,5.5),"pt"))

f = plot_model(fm_igm[[3]],type = "pred",terms=c("nexpose_s[-0.6317,0.145]","depth_s[-1.253899,0.8868936]"),title="",show.legend = F,colors = "bw")+
  scale_x_continuous(breaks=c(-0.6317,0.145),#1,2
                     labels=c("1","2"))+
  scale_y_continuous(breaks=y$predicted,
                     labels = c("1.6%","3.3%","2.8%","2.1%"))+
  labs(x = "Number of exposures", y = "Predicted probability of\nretweeting")+
  theme_bw()+theme(plot.title = element_blank(),plot.margin = unit(c(5.5,5.5,12,5.5),"pt"))

g = plot_model(fm_igm[[3]],type = "pred",terms=c("sim_s[-0.4371,2.787166]","depth_s[-1.253899,0.8868936]"),title="",show.legend = F,colors = "bw")+
  scale_x_continuous(breaks=c(-0.4371,2.787166),#0,0.5
                     labels=c("0","0.5"))+
  scale_y_continuous(breaks=z$predicted,
                     labels = c("3.1%","1.5%","14.8%","20.4%"))+
  labs(x = "Interaction similarity", y = "Predicted probability of\nretweeting")+
  theme_bw()+theme(plot.title = element_blank(),plot.margin = unit(c(5.5,5.5,12,5.5),"pt"))

multiplot(a,b,c,e,f,g,cols = 2)
topptx(filename = "inter_marginal.pptx")


####### SENSITIVITY ANALYSIS RTs
load("fm_igm.Rdata")
d = data%>%group_by(seed_user)%>%summarise(n=length(sharers))
data = data%>%group_by(tid)%>%mutate(ntids = sum(selected))

#
senr = list()
i=1
for (n in c(100,95,90,85,80,75,70)){
  fmx = update(fm_igm[[3]],data=data[data$ntids<=n,])
  senr[[i]] = fmx
  i=i+1
  print(i)
}
save(senr,file='senr_v2.Rdata')

#
load("senr_v2.Rdata")
civ = list()
j=1
for(i in senr){
  se <- sqrt(diag(vcov(i)))
  tab <- cbind(Est = fixef(i), LL = fixef(i) - 2.58 * se, UL = fixef(i) + 2.58 * se)
  civ[[j]] = tab
  j=j+1
}

library(ggplot2)
fgs = list()
j=1
for (i in c(2,3,4,5,6,8,17,24,25,26)){
  tmpd = data.frame(t(sapply(1:7,function(x)civ[[x]][i,])))
  tmpd$RTs = c(100,95,90,85,80,75,70)
  fg = ggplot(tmpd, aes(x=RTs, y=Est)) + geom_pointrange(aes(ymin=LL, ymax=UL))
  fgs[[j]] = fg
  j=j+1
}

multiplot(fgs[[1]],fgs[[2]],fgs[[3]],fgs[[4]],fgs[[5]],fgs[[6]],fgs[[7]],fgs[[8]],fgs[[9]],fgs[[10]],cols=2)
topptx(file="sensitivity_rts.pptx")

######## SENSITIVITY 2 (likes)
tids = unique(data$tid)
tweets = tweets_later[tweets_later$tweetid%in%tids,]
tweets$favorite_count = as.numeric(tweets$favorite_count)
tweets$retweet_count = as.numeric(tweets$retweet_count)

data = merge(data,tweets[,c("tweetid","favorite_count")],by.x='tid',by.y='tweetid',all.x=T)

#
load("senr2_v2.Rdata")
civ = list()
j=1
for(i in senr2){
  se <- sqrt(diag(vcov(i)))
  tab <- cbind(Est = fixef(i), LL = fixef(i) - 2.58 * se, UL = fixef(i) + 2.58 * se)
  civ[[j]] = tab
  j=j+1
}

library(ggplot2)
fgs = list()
j=1
for (i in c(2,3,4,5,6,8,17,24,25,26)){
  tmpd = data.frame(t(sapply(1:6,function(x)civ[[x]][i,])))
  tmpd$Likes = c(10,30,50,70,90,110)
  fg = ggplot(tmpd, aes(x=Likes, y=Est)) + geom_pointrange(aes(ymin=LL, ymax=UL))
  fgs[[j]] = fg
  j=j+1
}

multiplot(fgs[[1]],fgs[[2]],fgs[[3]],fgs[[4]],fgs[[5]],fgs[[6]],fgs[[7]],fgs[[8]],fgs[[9]],fgs[[10]],cols=2)
topptx(file="sensitivity2_lks.pptx")

######## CORRELATION ########
library("Hmisc")
rcorr(as.matrix(data[data$selected==1,c("inter","sim","nexpose","sr",'reciprocity','rtb','depth')]),type='spearman')
rcorr(as.matrix(data[data$selected==0,c("inter","sim","nexpose","sr",'reciprocity','rtb','depth')]),type='spearman')

##
library(ggplot2)
library(scales)

distr = data%>%group_by(tid)%>%summarise(depth=max(depth),size=length(tid))
dpt = as.data.frame(table(distr$depth),stringsAsFactors = F)
dpt$Depth = as.numeric(dpt$Var1)
dpt$Probability = rev(cumsum(rev(dpt$Freq/sum(dpt$Freq))))

sze = as.data.frame(table(distr$size),stringsAsFactors = F)
sze$Size = as.numeric(sze$Var1)
plot(dpt,log="xy")
plot(sze,log="xy")
plot(distr$size,distr$depth,log="xy")

p0 <- ggplot(distr, aes(x = size, y = depth))+ geom_point()+ 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) 
  theme_bw() + xlab("Cascade size") + ylab("depth") 

p1 <- ggplot(distr, aes(x=depth)) +
  geom_point(aes(y = 1 - ..y..), stat='ecdf') +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))
  theme_bw()  + xlab("Cascade depth") + ylab("CCDF") 

p2 <- ggplot(distr, aes(x=size)) +
  geom_point(aes(y = 1 - ..y..), stat='ecdf') +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + 
  theme_bw() + xlab("Cascade size") + ylab("CCDF") 

#multiplot(p1,p2,p0,cols = 2)
library(ggpubr)
ggarrange(
  ggarrange(p1, p2, ncol = 2, labels = c("A", "B")), 
  p0,
  nrow = 2, 
  labels =c("","C")       # Label of the line plot
) 

topptx(file="distr.pptx")

## case control ratio
table(data$selected)
