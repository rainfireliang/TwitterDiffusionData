setwd("D:/MISC/Random Tweets/Data and Analyses")
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

###
load("diffNet.Rdata")
load("UserInf.Rdata")

diffNet = diffNet[diffNet$sharers!=diffNet$seed_user&diffNet$depth>0,]#21188
diffNet$idds =NULL

userInf$id_str = paste(userInf$id_str)
diff_data = merge(diffNet, userInf[,c("id_str","followers_count","friends_count","location")],by.x="sharers",by.y = 'id_str',all.x=T)
diff_data = merge(diff_data, userInf[,c("id_str","followers_count","friends_count","location")],by.x="seed_user",by.y = 'id_str',all.x=T)

diff_data = diff_data[,c("tid","sharers","seed_user","time","n",'depth',"followers_count.x","friends_count.x","location.x","followers_count.y","friends_count.y","location.y")]
colnames(diff_data) = gsub(".x",".sharers",colnames(diff_data))
colnames(diff_data) = gsub(".y",".seed",colnames(diff_data))
diff_data = diff_data[order(diff_data$tid,diff_data$time),]
save(diff_data,file = "diff_data.Rdata") #data in use


####

diff_data$s_d = ifelse(diff_data$followers_count.seed>diff_data$followers_count.sharers,1,0)

diff_data$s_d = (diff_data$followers_count.seed-diff_data$followers_count.sharers)/diff_data$followers_count.sharers

x = diff_data %>% dplyr::group_by(depth) %>% dplyr::summarise(n=median(s_d,na.rm = T))
x$p = x$n/x$m
