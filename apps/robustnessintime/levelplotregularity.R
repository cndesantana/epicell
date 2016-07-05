library(latticeExtra)

folders<-system("ls",intern=TRUE)
nfolders<-(length(folders))
i<-1
ncells<-array(0,nfolders);
lambdasd<-array(0,nfolders);
gammasd<-array(0,nfolders);
reg<-array(0,nfolders);

for(i in 1:nfolders){
#  celldivfile<-paste("./",folders[i],"/tissue_celldiv.sta",sep="")
#  celldiv<-read.csv(celldivfile,sep=" ",header=FALSE)
#  ncells[i]<-max(celldiv$V2);
#  lambdasd[i]<-celldiv$V6[1]
#  gammasd[i]<-celldiv$V7[1]
  
  
  regularityfile<-paste("./",folders[i],"/regularness_celldiv.reg",sep="")
  regularity<-read.csv(regularityfile,sep=" ",header=FALSE)
  reg[i]<-regularity$V6[length(regularity$V6)]
}

png("levelplot_regularity.png",width=1980,height=1280)
grid<-expand.grid(x=unique(sort(lambdasd)),y=unique(sort(gammasd)))
grid$z<-reg
levelplot(reg ~ lambdasd * gammasd,data=grid,col.regions=terrain.colors(100),
          xlab=expression(Lambda),ylab=expression(Gamma),cuts=5,
          log="xy",panel = function(...){panel.levelplot(...)})
dev.off()



png("levelplot_area.png",width=1980,height=1280)
grid<-expand.grid(x=unique(sort(lambdasd)),y=unique(sort(gammasd)))
grid$z<-reg
levelplot(log(reg) ~ lambdasd * gammasd,data=grid,col.regions=terrain.colors(100),
          xlab=expression(Lambda),ylab=expression(Gamma),cuts=5,
          log="xy",panel = function(...){panel.levelplot(...)})
dev.off()
