setwd("~")
options(scipen = 999)

library(dplyr)
library(plyr)

#
load("sampledRTs.Rdata")
sampledRTs$t = as.POSIXct(sampledRTs$created_at,format='%a %b %d %H:%M:%S %z %Y', tz='GMT')
sampledRTs = sampledRTs[order(sampledRTs$retweeted_tweetid,sampledRTs$t),]
retws = sampledRTs[,c("id","retweeted_userid","t","retweeted_tweetid")]
colnames(retws) = c("sharers","seed_user","time","tid")
uids = unique(c(retws$seed_user,retws$sharers))

#
followings = readRDS("friends.rds") # dowload here: https://www.dropbox.com/s/fna98ccv5ooaalj/friends.rds?dl=0
colnames(followings) = c("egos","followees")
followings <- followings %>%
  filter(egos%in%uids) %>%
  filter(followees%in%uids)%>%
  collect(n=Inf)

tids = unique(retws$tid) #988 [21,329]

# reconstruct diffusion path function
recontruct = function(retws,followings){
  retws = retws[retws$seed_user!=retws$sharers,]
  if (nrow(followings)>0){
    dn = merge(followings,retws[,c('sharers','time')],by.x='egos',by.y='sharers')
    dn = merge(dn,retws[,c('sharers','time')],by.x='followees',by.y='sharers')
    
    dn = dn[which(dn$time.x>dn$time.y),]
    
    if (nrow(dn)>0){
      dn = plyr::ddply(dn,.(egos),mutate,n=length(followees),mint=max(time.y))
      dn = dn[which(dn$time.y==dn$mint),]
      dn = dn[,c('egos','followees','time.x','n')]
      colnames(dn) = c('sharers','seed_user','time','n')
      
      rest = retws[which(!retws$sharers %in% dn$sharers),c('sharers','seed_user','time')]
      rest$n=1
      
      tres = rbind(dn,rest)
      tres = tres[order(tres$time),]
    } else {
      tres = retws[,c('sharers','seed_user','time')]
      tres$n=1
    }
    
  } else {
    tres = retws[,c('sharers','seed_user','time')]
    tres$n=1
  }
  return (tres) # n indicate the number of exposures
}

# test: return a dataframe with 4 colums: sharers, shared from, time, number of exposures
# diffnet = recontruct(retws[retws$tid==tids[1],1:3],followings)

#
diffNet <- data.frame()
jj <- 0
for (tid in tids){
  y = retws[retws$tid==tid,1:3]
  y = y%>% distinct(sharers,.keep_all = T)
  x = recontruct(y,followings)
  x$tid = tid
  diffNet = bind_rows(diffNet,x)
  gc()
  
  jj = jj+1
  print (jj)
}
diffNet = diffNet %>% distinct(sharers,tid,.keep_all = T)
save(diffNet,file="diffNet.Rdata") #21188

#
library(igraph)
g = graph.data.frame(diffNet[diffNet$tid == tids[886],c('seed_user','sharers')],directed = T)
plot(g,edge.arrow.size=0.3,layout=layout_as_tree)

# depth
depth = function(i){
  n = diffNet[diffNet$tid==i,c("seed_user","sharers")] # from,to
  g = simplify(graph.data.frame(n,directed = T))
  ids=distances(g,mode='out')[1,]
  df=cbind(i,names(ids),ids)
  colnames(df) = NULL
  return(df)
}
diffD=matrix(ncol = 3)
tids=unique(diffNet$tid)

for (i in 1:length(tids)){
  diffD=rbind(diffD,depth(tids[i]))
  print(i)
}
colnames(diffD)=c('tid','sharers','depth')

diffD=as.data.frame(diffD,stringsAsFactors = F)
diffD=diffD[!is.na(diffD$sharers),]
diffD$depth=as.numeric(paste(diffD$depth))
save(diffD,file='diffD.Rdata')

diffD$idds=paste(diffD$tid,diffD$sharers)
diffNet$idds=paste(diffNet$tid,diffNet$sharers)

diffNet=merge(diffNet,diffD[,c('idds','depth')],by='idds')
diffNet=diffNet[order(diffNet$tid,diffNet$time),]
save(diffNet,file='diffNet.Rdata')


