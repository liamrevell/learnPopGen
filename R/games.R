hawk.dove<-function(p=c(0.01,0.99),M=NULL,time=100){
	if(is.null(M)) M<-matrix(c(0.6,1.5,0.5,1.0),2,2,byrow=T)
	if(is.null(dimnames(M))) dimnames(M)<-list(c("hawk","dove"),c("hawk","dove"))
	interactors<-rownames(M)
	message("Pay-off matrix:")
	print(M)
	if(sum(p)!=1) p<-p/sum(p)
	hawk<-p[1]; dove<-p[2]
	wbar<-hawk*hawk*M[1,1]+hawk*dove*M[1,2]+dove*hawk*M[2,1]+dove*dove*M[2,2]
	whawk<-hawk*M[1,1]+dove*M[1,2]
	wdove<-hawk*M[2,1]+dove*M[2,2]
	for(i in 2:time){
		hawk[i]<-hawk[i-1]*whawk[i-1]/wbar[i-1]
		dove[i]<-dove[i-1]*wdove[i-1]/wbar[i-1]
		whawk[i]<-hawk[i]*M[1,1]+dove[i]*M[1,2]
		wdove[i]<-hawk[i]*M[2,1]+dove[i]*M[2,2]
		wbar[i]<-hawk[i]*whawk[i]+dove[i]*wdove[i]
	}
	time<-1:time
	layout(c(1,2))
	par(mar=c(5.1,4.1,1.1,2.1))
	plot(time,dove,type="l",ylim=c(0,1),ylab="frequency",col="blue",lwd=2)
	lines(time,hawk,type="l",col="red",lwd=2)
	legend(max(time),1,interactors,col=c("red","blue"),bg="gray90",lwd=2,xjust=1,yjust=1)
	plot(time,wbar,type="l",lwd=2,ylim=c(0,max(M)),ylab="mean fitness")
	lines(time,wdove,type="l",col="blue",lwd=2)
	lines(time,whawk,type="l",col="red",lwd=2)
	legend(max(time),max(M),c("overall",interactors),col=c("black","red","blue"),bg="gray90",lwd=2,xjust=1,yjust=1)
}
