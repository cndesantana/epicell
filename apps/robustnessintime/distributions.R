scenarios<-system("ls",intern=TRUE)
nscenarios<-length(scenarios)
scenarios<-scenarios[-c(nscenarios,nscenarios-1,nscenarios-2)];#removing the script files from the list of files
nscenarios<-length(scenarios);#after removing the extra files
for(i in 1:nscenarios){
   folder<-scenarios[i];
   #folder<-"./1_10_1_-0.8_0.12_5_5_0.0_0.0"
   lambda<-unlist(strsplit(folder,split='_'))[4]
   gamma<-unlist(strsplit(folder,split='_'))[5]
   sdlambda<-unlist(strsplit(folder,split='_'))[8]
   sdgamma<-unlist(strsplit(folder,split='_'))[9]
   cat(paste("(sdGamma,sdLambda) -------> (",sdgamma,",",sdlambda,")",sep=""),sep="\n");
   picture<-paste("pictures/dist_G_",gamma,"_sdG_",sdgamma,"_L_",lambda,"_sdL_",sdlambda,".png",sep="")
   
   dist<-read.csv(paste(folder,"/distributionofshapes.dist",sep=""),sep=" ",header=FALSE)
   steps<-unique(sort(dist$V1))
   nsteps<-length(steps)
   nshapes<-length(unique(sort(dist[1,])))-1
   meanshape<-array(0,c(nsteps,nshapes))
   sdshape<-array(0,c(nsteps,nshapes))
   
   for(j in 1:nshapes){
       for(i in 1:nsteps){
           pos<-which(dist$V1==steps[i])
           meanshape[i,j]<-mean(dist[pos,j+1],na.rm=TRUE)
           sdshape[i,j]<-sd(dist[pos,j+1],na.rm=TRUE)
       }
   }
   coeffvarlambda<-abs(as.numeric(sdlambda)/as.numeric(lambda))
   coeffvargamma<-abs(as.numeric(sdgamma)/as.numeric(gamma))
   png(picture,width=1980,height=1280,res=300);
   plot(meanshape[,1],type="b",col="black",xlim=c(0,150),ylim=c(-0.1,1),xlab="Time",ylab="Frequency",main=paste("Distribution of shapes\n(sdLambda,sdGamma)  =  (",coeffvarlambda,",",coeffvargamma,")",sep=""))
   points(meanshape[,2],type="b",col=1)
   points(meanshape[,3],type="b",col=2)
   points(meanshape[,4],type="b",col=3)
   points(meanshape[,5],type="b",col=4)
   points(meanshape[,6],type="b",col=5)
   points(meanshape[,7],type="b",col=6)
   legend("topright",legend=c("<4sides","5sides","6sides","7sides","8sides",">9sides"),col=c(1,2,3,4,5,6),pch=1,lwd=1);
   dev.off();
}
