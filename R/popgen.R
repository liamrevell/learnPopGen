## functions by Liam Revell 2012 (some small updates 2017, 2018)

selection<-function(p0=0.01,w=c(1.0,0.9,0.8),time=100,show="p",pause=0,...){
	if(hasArg(add)) add<-list(...)$add
	else add<-FALSE
	if(hasArg(color)) color<-list(...)$color
	else color<-"black"
	if(hasArg(equil)) equil<-list(...)$equil
	else equil<-FALSE
	if(hasArg(lwd)) lwd<-list(...)$lwd
	else lwd<-2
	if(hasArg(lty)) lty<-list(...)$lty
	else lty<-"solid"
	## compute equilibrium
	if(w[2]>w[1]&&w[2]>w[3]){ ## overdominance
		W<-w/w[2]
		eq<-(1-W[3])/((1-W[1])+(1-W[3]))
	} else if(w[2]<w[1]&&w[2]<w[3]){ ## underdominance
		W<-w/w[2]
		eq<-(1-W[3])/((1-W[1])+(1-W[3]))
	} else if(1%in%which(w==max(w))){ 
		eq<-1
	} else eq<-0
	if(show=="surface"){
		p<-0:100/100
		wbar<-p^2*w[1]+2*p*(1-p)*w[2]+(1-p)^2*w[3]
		plot(p,wbar,type="l",ylim=c(0,1.1*max(w)),
			main=expression(paste("mean fitness (",
			bar(w),")",sep="")),ylab=expression(bar(w)),col=color)
		if(equil) abline(v=eq,lty="dotted")
	}
	else if(show=="deltap"){
		p<-0:100/100
		wbar<-p^2*w[1]+2*p*(1-p)*w[2]+(1-p)^2*w[3]
		deltap<-(p/wbar)*(p*w[1]+(1-p)*w[2]-wbar)
		plot(p,deltap,type="l",main=expression(paste(Delta,"p as a function of p",
			sep="")),ylab=expression(paste(Delta,"p",sep="")),col=color)
		lines(c(0,1),c(0,0),lty=2)
		if(equil) abline(v=eq,lty="dotted")
	} else {
		if(show=="cobweb"){
			p<-0:100/100
			wbar<-p^2*w[1]+2*p*(1-p)*w[2]+(1-p)^2*w[3]
			p2<-(p/wbar)*(p*w[1]+(1-p)*w[2]-wbar)+p
			plot(p,p2,type="l",xlab=expression(p[t]),ylab=expression(p[t+1]),
				main=expression(paste(p[t+1]," as a function of ",p[t],sep="")),
				col=color)
			lines(c(0,1),c(0,1),lty=2)
			if(equil){ 
				abline(v=eq,lty="dotted")
				abline(h=eq,lty="dotted")
			}
			dev.flush()
		}
		p<-wbar<-vector()
		p[1]<-p0
		wbar[1]<-p[1]^2*w[1]+2*p[1]*(1-p[1])*w[2]+(1-p[1])^2*w[3]
		for(i in 2:time){
			p[i]<-p[i-1]
			p[i]<-(p[i]^2*w[1]+p[i]*(1-p[i])*w[2])/wbar[i-1]
			wbar[i]<-p[i]^2*w[1]+2*p[i]*(1-p[i])*w[2]+(1-p[i])^2*w[3]
			ii<-(i-1):i
			if(show=="p"){
				if(i==2 && !add) plot(1:i,p,type="l",xlim=c(0,time),ylim=c(0,1),xlab="time",
					main="frequency of A",col=color,lwd=lwd,lty=lty)
				else lines(ii,p[ii],type="l",col=color,lwd=lwd,lty=lty)
			} else if(show=="q" && !add){
				if(i==2) plot(1:i,1-p,type="l",xlim=c(0,time),ylim=c(0,1),xlab="time",
					ylab="q",main="frequency of a",col=color,lwd=lwd,lty=lty)
				else lines(ii,1-p[ii],type="l",col=color,lwd=lwd,lty=lty)
			} else if(show=="fitness" && !add){
				if(i==2) plot(1:i,wbar,type="l",xlim=c(0,time),
					ylim=c(0,1.1*max(w)),
					xlab="time",main=expression(paste("mean fitness (",
					bar(w),")",sep="")),ylab=expression(bar(w)),col=color)
				else lines(ii,wbar[ii],type="l",col=color)
			} else if(show=="cobweb"){
				lines(c(p[i-1],p[i-1]),c(p[i-1],p[i]),col=color)
				lines(c(p[i-1],p[i]),c(p[i],p[i]),col=color)
			} else {
				message("not a recognized option")
				break
			}
			if(equil) {
				if(show=="p") abline(h=eq,lty="dotted")
				if(show=="q") abline(h=1-eq,lty="dotted")
				if(show=="fitness") abline(h=eq^2*w[1]+2*eq*(1-eq)*w[2]+(1-eq)^2*w[3],
					lty="dotted")
			}
			dev.flush()
			Sys.sleep(pause)
		}
	}
}

freqdep<-function(p0=0.01,s=0,time=100,show="p",pause=0,...){
	if(hasArg(color)) color<-list(...)$color
	else color<-"black"
	p<-0:100/100
	f11<-p^2
	f12<-2*p*(1-p)
	f22<-(1-p)^2
	wbar<-f11*(1-3*f12+3*f22)+f12*(1-s*f12)+f22*(1-3*f12+3*f11)
	wmax<-max(wbar)
	if(show=="surface"){
		plot(p,wbar,type="l",ylim=c(0,max(1,wmax)),
			main=expression(paste("mean fitness (",
			bar(w),")",sep="")),ylab=expression(bar(w)),col=color)
	}
	else if(show=="deltap"){
		p<-0:100/100
		f11<-p^2
		f12<-2*p*(1-p)
		f22<-(1-p)^2
		wbar<-f11*(1-3*f12+3*f22)+f12*(1-s*f12)+f22*(1-3*f12+3*f11)
		w11<-1-3*f12+3*f22
		w12<-1-s*f12
		w22<-1-3*f12+3*f11
		deltap<-(p/wbar)*(p*w11+(1-p)*w12-wbar)
		plot(p,deltap,type="l",main=expression(paste(Delta,"p as a function of p",
			sep="")),ylab=expression(paste(Delta,"p",sep="")),col=color)
		lines(c(0,1),c(0,0),lty=2)
	} else {
		if(show=="cobweb"){
			p<-0:100/100
			f11<-p^2; f12<-2*p*(1-p); f22<-(1-p)^2
			wbar<-f11*(1-3*f12+3*f22)+f12*(1-s*f12)+f22*(1-3*f12+3*f11)
			w11<-1-3*f12+3*f22
			w12<-1-s*f12
			w22<-1-3*f12+3*f11
			p2<-(p/wbar)*(p*w11+(1-p)*w12-wbar)+p
			plot(p,p2,type="l",xlab=expression(p[t]),ylab=expression(p[t+1]),
				main=expression(paste(p[t+1]," as a function of ",p[t],sep="")))
			lines(c(0,1),c(0,1),lty=2)
			dev.flush()
		}
		p<-wbar<-vector()
		p[1]<-p0
		f11<-p[1]^2
		f12<-2*p[1]*(1-p[1])
		f22<-(1-p[1])^2
		wbar[1]<-f11*(1-3*f12+3*f22)+f12*(1-s*f12)+f22*(1-3*f12+3*f11)
		for(i in 2:time){
			p[i]<-p[i-1]
			w11<-1-3*f12+3*f22
			w12<-1-s*f12
			w22<-1-3*f12+3*f11
			p[i]<-(p[i]^2*w11+p[i]*(1-p[i])*w12)/wbar[i-1]
			f11<-p[i]^2; f12<-2*p[i]*(1-p[i]); f22<-(1-p[i])^2
			wbar[i]<-f11*(1-3*f12+3*f22)+f12*(1-s*f12)+f22*(1-3*f12+3*f11)
			ii<-(i-1):i
			if(show=="p"){
				if(i==2) plot(1:i,p,type="l",xlim=c(0,time),ylim=c(0,1),xlab="time",
					main="frequency of A",col=color)
				else lines(ii,p[ii],type="l",col=color)
			} else if(show=="q"){
				if(i==2) plot(1:i,1-p,type="l",xlim=c(0,time),ylim=c(0,1),xlab="time",
					ylab="q",main="frequency of a",col=color)
				else lines(ii,1-p[ii],type="l",col=color)
			} else if(show=="fitness"){
				if(i==2) plot(1:i,wbar,type="l",xlim=c(0,time),ylim=c(0,max(1,wmax)),
					xlab="time",
					main=expression(paste("mean fitness (",bar(w),")",sep="")),
					ylab=expression(bar(w)),col=color)
				else lines(ii,wbar[ii],type="l",col=color)
			} else if(show=="cobweb"){
				lines(c(p[i-1],p[i-1]),c(p[i-1],p[i]),col=color)
				lines(c(p[i-1],p[i]),c(p[i],p[i]),col=color)
			} else {
				message("not a recognized option")
				break
			}
			dev.flush()
			Sys.sleep(pause)
		}
	}
}

sexratio<-function(p0=0.01,time=40,show="p",pause=0,sex.Aa=c(0.5,0.5)){
	p<-fm<-ff<-wm<-wf<-wbar<-vector()
	p[1]<-p0
	fm[1]<-p[1]^2+2*p[1]*(1-p[1])*sex.Aa[1]
	ff[1]<-(1-p[1])^2+2*p[1]*(1-p[1])*sex.Aa[2]
	wm[1]<-0.5/fm[1]
	wf[1]<-0.5/ff[1]
	wbar[1]<-fm[1]*wm[1]+wf[1]*ff[1]
	for(i in 2:time){
		p[i]<-p[i-1]
		w<-c(wm[i-1],sex.Aa[1]*wm[i-1]+sex.Aa[2]*wf[i-1],wf[i-1])
		p[i]<-(p[i]^2*w[1]+p[i]*(1-p[i])*w[2])/wbar[i-1]
		fm[i]<-p[i]^2+2*p[i]*(1-p[i])*sex.Aa[1]
		ff[i]<-(1-p[i])^2+2*p[i]*(1-p[i])*sex.Aa[2]
		wm[i]<-0.5/fm[i]
		wf[i]<-0.5/ff[i]
		wbar[i]<-fm[i]*wm[i]+wf[i]*ff[i]
		ii<-(i-1):i
		if(show=="p"){
			if(i==2) plot(1:i,p,type="l",xlim=c(1,time),ylim=c(0,1),xlab="time",
				main="frequency of A")
			else lines(ii,p[ii],type="l")
		} else if(show=="sex-ratio"){
			if(i==2){ 
				plot(1:i,fm,type="l",xlim=c(1,time),ylim=c(0,1),xlab="time",
					main="frequency of each sex",col=make.transparent("blue",0.5),
					lwd=2,ylab="relative frequencies of each sex")
				lines(1:i,ff,lwd=2,col=make.transparent("red",0.5))
				legend(x="topright",c("males","females"),lwd=2,
					col=c(make.transparent("blue",0.5),
					make.transparent("red",0.5)))
			} else { 
				lines(ii,fm[ii],type="l",lwd=2,col=make.transparent("blue",0.5))
				lines(ii,ff[ii],type="l",lwd=2,col=make.transparent("red",0.5))
			}
		} else if(show=="fitness"){
			if(i==2){
				plot(1:i,wbar,type="l",xlim=c(1,time),ylim=c(min(wbar,wm,wf),max(wbar,wm,wf)),
					xlab="time",main="fitness",lwd=2,col=make.transparent("grey",1/3),
					log="y")
				lines(1:i,wm,lwd=2,col=make.transparent("blue",1/3))
				lines(1:i,wf,lwd=2,col=make.transparent("red",1/3))
				legend(x="topright",c("average","males","females"),lwd=2,
					col=c(make.transparent("grey",1/3),
					make.transparent("blue",1/3),
					make.transparent("red",1/3)))
			} else { 
				lines(ii,wbar[ii],lwd=2,col=make.transparent("grey",1/3))
				lines(ii,wm[ii],lwd=2,col=make.transparent("blue",1/3))
				lines(ii,wf[ii],lwd=2,col=make.transparent("red",1/3))
			}
		}
		Sys.sleep(pause)
		dev.flush()
	}
}

mutation.selection<-function(p0=1.0,w=c(1,1),u=0.001,time=100,show="q",
	pause=0,ylim=c(0,1)){
	p<-wbar<-vector(); p[1]<-p0
	wbar[1]<-p[1]^2*1.0+2*p[1]*(1-p[1])*w[1]+(1-p[1])^2*w[2]
	for(i in 2:time){
		p[i]<-p[i-1]
		p[i]<-(1-u)*(p[i]^2*1.0+p[i]*(1-p[i])*w[1])/wbar[i-1]
		wbar[i]<-p[i]^2*1.0+2*p[i]*(1-p[i])*w[1]+(1-p[i])^2*w[2]
		ii<-(i-1):i
		if(show=="p"){
			if(i==2) plot(1:i,p,type="l",xlim=c(0,time),ylim=c(0,1),xlab="time",
				main="frequency of A")
			else lines(ii,p[ii],type="l")
		} else if(show=="q"){
			if(i==2) plot(1:i,1-p,type="l",xlim=c(0,time),ylim=ylim,xlab="time",
				ylab="q",main="frequency of a")
			else lines(ii,1-p[ii],type="l")
		} else if(show=="fitness"){
			if(i==2) plot(1:i,wbar/max(w),type="l",xlim=c(0,time),ylim=c(0,1),
				xlab="time",main="mean fitness")
			else lines(ii,wbar[ii]/max(w),type="l")
		} else {
			message("not a recognized option")
			break
		}
		dev.flush()
		Sys.sleep(pause)
	}
}

genetic.drift<-function(p0=0.5,Ne=20,nrep=10,time=100,show="p",pause=0.1,
	...){
	if(hasArg(colors)) colors<-list(...)$colors
	else colors<-rep("black",nrep)
	if(hasArg(lwd)) lwd<-list(...)$lwd
	else lwd<-1
	freq<-repMatrix(matrix(0,time,3),nrep)
	p<-matrix(0,time,nrep)
	hbar<-vector()
	genotypes<-list()
	for(i in 1:nrep) genotypes[[i]]<-matrix(sample(c(rep(1,round(2*p0*Ne)),
		rep(0,2*Ne-round(2*p0*Ne)))),Ne,2)
	for(i in 1:nrep) freq[[i]][1,]<-hist(rowSums(genotypes[[i]]),
		c(-0.5,0.5,1.5,2.5),plot=FALSE)$density
	for(i in 1:nrep) p[1,i]<-mean(genotypes[[i]])
	X<-matrix(NA,nrep,3,dimnames=list(NULL,c("aa","Aa","AA")))
	for(i in 1:nrep) X[i,]<-freq[[i]][1,]
	hbar[1]<-mean(X[,2])
	if(show=="genotypes") barplot(X,ylim=c(0,1),main="genotype frequencies",beside=TRUE)
	if(show=="fixed"){
		fixedA<-sum(sapply(genotypes,sum)==(2*Ne))
		fixeda<-sum(sapply(genotypes,sum)==0)
		barplot(c(fixeda,fixedA)/nrep,ylim=c(0,1),names.arg=c("a","A"),
			main="populations fixed",ylab="frequency")
	}	
	for(i in 2:time){
		new.gen<-repMatrix(matrix(NA,Ne,2),nrep)
		for(j in 1:nrep) for(k in 1:Ne) new.gen[[j]][k,]<-sample(genotypes[[j]],size=2)
		genotypes<-new.gen
		for(j in 1:nrep){ 
			freq[[j]][i,]<-hist(rowSums(genotypes[[j]]),c(-0.5,0.5,1.5,2.5),
				plot=FALSE)$density
			p[i,j]<-mean(genotypes[[j]])
		}
		for(j in 1:nrep) X[j,]<-freq[[j]][i,]
		if(show=="genotypes") barplot(X,ylim=c(0,1),main="genotype frequencies",
			beside=TRUE,col=if(all(colors!="black")) colors else NULL)
		else if(show=="p"){ 
			if(i==2){ 
				plot(1:i,p[1:i,1],type="l",ylim=c(0,1),xlim=c(1,time),
					main="frequency of A",xlab="time",ylab="p",col=colors[1],
					lwd=lwd)
				if(p0<=0.5) text(paste("N =",Ne,sep=" "),x=0,y=1,pos=4)
				else text(paste("N =",Ne,sep=" "),x=0,y=0,pos=4)
				if(nrep>1) for(j in 2:nrep) lines(1:i,p[1:i,j],col=colors[j],
					lwd=lwd)
			} else for(j in 1:nrep) lines((i-1):i,p[(i-1):i,j],col=colors[j],
				lwd=lwd)
		}
		else if(show=="heterozygosity"){
			hbar[i]<-mean(X[,2])
			if(i==2) plot(1:time,2*p0*(1-p0)*(1-1/(2*Ne))^(0:(time-1)),type="l",
				lwd=2,col="red",main="heterozygosity",xlab="time",ylab="f(Aa)",ylim=c(0,0.6))
			lines((i-1):i,hbar[(i-1):i])
		}
		else if(show=="fixed"){
			fixedA<-sum(sapply(genotypes,sum)==(2*Ne))
			fixeda<-sum(sapply(genotypes,sum)==0)
			barplot(c(fixeda,fixedA)/nrep,ylim=c(0,1),names.arg=c("a","A"),
				main="populations fixed",ylab="frequency")
		}
		dev.flush()
		Sys.sleep(pause)
	}
}

repMatrix<-function(X,times){
	Z<-list()
	for(i in 1:times) Z[[i]]<-X
	return(Z)
}

founder.event<-function(p0=0.5,Ne=1000,Nf=10,ttime=100,etime=50,show="p",...){
	if(hasArg(ltype)) ltype<-list(...)$ltype
	else ltype<-"s"
	genotypes<-matrix(sample(c(rep(1,round(2*p0*Ne)),
		rep(0,2*Ne-round(2*p0*Ne)))),Ne,2)
	p<-v<-vector()
	p[1]<-mean(genotypes)
	v[1]<-p[1]*(1-p[1])
	for(i in 2:ttime){
		if(i%in%etime){
			if(i==etime[1]){ 
				founder<-genotypes[sample(1:Ne,size=Nf),]
				genotypes<-founder
			}
			if((i+1)%in%etime){
				new.gen<-matrix(NA,Nf,2)
				for(j in 1:Nf) new.gen[j,]<-sample(genotypes,size=2,replace=TRUE)
			} else {
				new.gen<-matrix(NA,Ne,2)
				for(j in 1:Ne) new.gen[j,]<-sample(genotypes,size=2,replace=TRUE)
			}
		} else {
			new.gen<-matrix(NA,Ne,2) 
			for(j in 1:Ne) new.gen[j,]<-sample(genotypes,size=2,replace=TRUE)
		}
		genotypes<-new.gen
		p[i]<-mean(genotypes)
		v[i]<-p[i]*(1-p[i])
	}
	if(show=="p")
		plot(1:ttime,p,main="frequency of A",xlab="time",ylim=c(0,1),type=ltype,
			ylab=expression(p[A]))
	else if(show=="var")
		plot(1:ttime,v,main="genetic variation",xlab="time",ylim=c(0,0.3),type=ltype,
			ylab=expression(p[A]*(1-p[A])))
	if(length(etime)==1) abline(v=etime,col=make.transparent("grey",0.5),lwd=6)
	else rect(xleft=etime[1],ybottom=par()$usr[3],
			xright=etime[length(etime)],ytop=par()$usr[4],
			col=make.transparent("grey",0.5),border=NA)
}
