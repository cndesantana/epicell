library(colorRamps)
library(lattice);
library(latticeExtra);

getColor<-function(matdata,sdlambda){
    allcolors<-blue2red(7);
    allsdlambda<-unique(sort(matdata[,2]));
    mycolorix<-which(allsdlambda%in%sdlambda);
    if(length(mycolorix) == 0){
        mycolorix<-max(which(allsdlambda < sdlambda))+1;
    }
    mycolor<-(allcolors[mycolorix])
    return(mycolorix);
}

getMAXTS<-function(dat,ref){
	MAXTS<-dat$V1[which(dat$V2 > ref)[1]];
	return(MAXTS);
}

getMatData<-function(workdir){
    scenarios<-system("ls",intern=TRUE)
    nscenarios<-length(scenarios)
    scenarios<-scenarios[-c(nscenarios,nscenarios-1,nscenarios-2)];
    nscenarios<-length(scenarios)
    variables<-c("sdgamma", "sdlambda", "meanreg","sdreg","color","ncells","meanarea","meanperimeter","meanshape")
    matdata<-matrix(0,nscenarios,length(variables))
    colnames(matdata)<-variables
    for(i in 1:nscenarios){
        file1<-paste(scenarios[i],"/regularness_celldiv.reg",sep="")
        dat1<-read.csv(file1,sep=" ",header=FALSE)
        file2<-paste(scenarios[i],"/tissue_celldiv.sta",sep="")
        dat2<-read.csv(file2,sep=" ",header=FALSE)
        MAXTS<-getMAXTS(dat2,230);
	sdgamma<-dat1$V4[1]
        sdlambda<-dat1$V5[1]

	cat(paste("Folder ",i," --> MAXTS = ",MAXTS,sep=""),sep="\n");	
	if(!is.na(MAXTS)){
		#From REG file
        	posmax<-which(dat1$V1 == MAXTS)
        	meanreg<-mean(dat1$V6[posmax])
        	sdreg<-sd(dat1$V6[posmax])
        	mycolor<-getColor(matdata,sdlambda);
		#From STA file
        	posmax2<-which(dat2$V1 == MAXTS)
        	ncells<-mean(dat2$V2[posmax2])
        	meanarea<-mean(dat2$V8[posmax2])
        	meanperimeter<-mean(dat2$V9[posmax2])
        	meanshape<-mean(dat2$V10[posmax2])
        	matdata[i,]<-c(sdgamma,sdlambda,meanreg,sdreg,mycolor,ncells,meanarea,meanperimeter,meanshape)
    	}else{
        	matdata[i,]<-c(sdgamma,sdlambda,NA,NA,NA,NA,NA,NA,NA)
	}
    }

    return(matdata)
}

plotPerturbations<-function(matdata,currdir){    
    seqgamma<-order(matdata[,1])
    seqlambda<-order(matdata[,2])
    
#Regularity
    png("./gammaperturb_regularity.png",width=1980,height=1280,res=300)
    plot(matdata[seqgamma,1],matdata[seqgamma,3],log="x",type="b",xlab="Perturbation in gamma",ylab="Regularity",main="Perturbations vs. Regularity",pch=matdata[seqgamma,5],col=matdata[seqgamma,5],ylim=c(0.2,0.9))
    dev.off()
}

args<-commandArgs(TRUE);
#workdir<-args[1];
workdir<-"home/cdesantana/Data/EpiphysX/epicell/apps/robustnessintime/tmp/model_2/gammaN_0.12_lambdaN_-0.8";
currdir<-"home/cdesantana/Data/EpiphysX/epicell/apps/robustnessintime";
#setwd(workdir);
matdata<-getMatData(workdir);
#setwd(currdir);
plotPerturbations(matdata,currdir)

meanlambda<--0.8;
meangamma<-0.12;
lambdas<-sort(unique(matdata[,2]))
gammas<-sort(unique(matdata[,1]))
mygrid<-expand.grid(x=lambdas,y=gammas)
ngridpoints<-length(mygrid$x)

regularity<-array(NA,ngridpoints)
for(i in 1:ngridpoints){
	regdata<-matdata[which((matdata[,1]==mygrid$y[i]) & (matdata[,2]==mygrid$x[i])),3];
	if(length(regdata) > 0){
		regularity[i]<-regdata;
	}
}

ncells<-array(NA,ngridpoints)
for(i in 1:ngridpoints){
	ncellsdata<-matdata[which((matdata[,1]==mygrid$y[i]) & (matdata[,2]==mygrid$x[i])),6];
	if(length(ncellsdata) > 0){
		ncells[i]<-ncellsdata;
	}
}


area<-array(NA,ngridpoints)
for(i in 1:ngridpoints){
	areadata<-matdata[which((matdata[,1]==mygrid$y[i]) & (matdata[,2]==mygrid$x[i])),7];
	if(length(areadata) > 0){
		area[i]<-areadata;
	}
}

perimeter<-array(NA,ngridpoints)
for(i in 1:ngridpoints){
	perimeterdata<-matdata[which((matdata[,1]==mygrid$y[i]) & (matdata[,2]==mygrid$x[i])),8];
	if(length(perimeterdata) > 0){
		perimeter[i]<-perimeterdata;
	}
}

shape<-array(NA,ngridpoints)
for(i in 1:ngridpoints){
	shapedata<-matdata[which((matdata[,1]==mygrid$y[i]) & (matdata[,2]==mygrid$x[i])),9];
	if(length(shapedata) > 0){
		shape[i]<-shapedata;
	}
}

cat("Saving...\n");
mygrid$regularity<-regularity;
mygrid$shape<-shape;
mygrid$ncells<-ncells;
mygrid$area<-area;
mygrid$perimeter<-perimeter;
save(mygrid,file="./mygrid.Rdat");
mygrid$x<-abs(mygrid$x/meanlambda);
mygrid$y<-abs(mygrid$y/meangamma);

cat("Plotting...\n");
myxlab<-substitute(paste(sigma,"(",Lambda,")",sep=""));
myylab<-expression(Gamma);
png("./levelplot_shape.png",width=1980,height=1280,res=300);
levelplot(shape ~ x * y,data=mygrid,col.regions=terrain.colors(10),
          xlab=myxlab,ylab=myylab,main="Mean Shape",cuts=3,
          log="xy",panel = panel.levelplot.raster, interpolate=TRUE,
          xlim=c(0,max(mygrid$x)),ylim=c(0,max(mygrid$y)))
dev.off();

png("./levelplot_ncells.png",width=1980,height=1280,res=300);
levelplot(ncells ~ x * y,data=mygrid,col.regions=terrain.colors(10),
          xlab=myxlab,ylab=myylab,main="Number of Cells",cuts=3,
          log="xy",panel = panel.levelplot.raster, interpolate=TRUE,
          xlim=c(0,max(mygrid$x)),ylim=c(0,max(mygrid$y)))
dev.off();

png("./levelplot_area.png",width=1980,height=1280,res=300);
levelplot(area ~ x * y,data=mygrid,col.regions=terrain.colors(10),
          xlab=myxlab,ylab=myylab,main="Mean Area",cuts=3,
          log="xy",panel = panel.levelplot.raster, interpolate=TRUE,
          xlim=c(0,max(mygrid$x)),ylim=c(0,max(mygrid$y)))
dev.off();

png("./levelplot_perimeter.png",width=1980,height=1280,res=300);
levelplot(perimeter ~ x * y,data=mygrid,col.regions=terrain.colors(10),
          xlab=myxlab,ylab=myylab,main="Mean Perimeter",cuts=3,
          log="xy",panel = panel.levelplot.raster, interpolate=TRUE,
          xlim=c(0,max(mygrid$x)),ylim=c(0,max(mygrid$y)))
dev.off();

png("./levelplot_regularity.png",width=1980,height=1280,res=300);
levelplot(regularity ~ x * y,data=mygrid,col.regions=terrain.colors(10),
          xlab=myxlab,ylab=myylab,main="Regularity",cuts=3,
          log="xy",panel = panel.levelplot.raster, interpolate=TRUE,
          xlim=c(0,max(mygrid$x)),ylim=c(0,max(mygrid$y)))
dev.off();
