# functions by Liam Revell 2012 (some small updates 2017, 2018)

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
	if(equil){
		if(w[2]>w[1]&&w[2]>w[3]){ ## overdominance
			W<-w/w[2]
			eq<-(1-W[3])/((1-W[1])+(1-W[3]))
		} else if(w[2]<w[1]&&w[2]<w[3]){ ## underdominance
			W<-w/w[2]
			eq<-(1-W[3])/((1-W[1])+(1-W[3]))
		} else if(1%in%which(w==max(w))){ 
			eq<-1
		} else eq<-0
	}
	if(show=="surface"){
		p<-0:100/100
		wbar<-p^2*w[1]+2*p*(1-p)*w[2]+(1-p)^2*w[3]
		plot(p,wbar,type="l",ylim=c(0,1),main=expression(paste("mean fitness (",
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
				if(i==2) plot(1:i,wbar/max(w),type="l",xlim=c(0,time),ylim=c(0,1),
					xlab="time",main="mean fitness",col=color)
				else lines(ii,wbar[ii]/max(w),type="l",col=color)
			} else if(show=="cobweb"){
				lines(c(p[i-1],p[i-1]),c(p[i-1],p[i]),col=color)
				lines(c(p[i-1],p[i]),c(p[i],p[i]),col=color)
			} else {
				message("not a recognized option")
				break
			}
			if(equil) abline(h=eq,lty="dotted")
			dev.flush()
			Sys.sleep(pause)
		}
	}
}

freqdep<-function(p0=0.01,s=0,time=100,show="p",pause=0){
	if(show=="surface"){
		p<-0:100/100
		f11<-p^2; f12<-2*p*(1-p); f22<-(1-p)^2
		wbar<-f11*(1-3*f12+3*f22)+f12*(1-s*f12)+f22*(1-3*f12+3*f11)
		plot(p,wbar,type="l",ylim=c(0,1),main="mean fitness")
	}
	else if(show=="deltap"){
		p<-0:100/100
		f11<-p^2; f12<-2*p*(1-p); f22<-(1-p)^2
		wbar<-f11*(1-3*f12+3*f22)+f12*(1-s*f12)+f22*(1-3*f12+3*f11)
		w11<-1-3*f12+3*f22
		w12<-1-s*f12
		w22<-1-3*f12+3*f11
		deltap<-(p/wbar)*(p*w11+(1-p)*w12-wbar)
		plot(p,deltap,type="l",main="delta p")
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
			plot(p,p2,type="l",xlab="p(t)",ylab="p(t+1)")
			lines(c(0,1),c(0,1),lty=2)
			dev.flush()
		}
		p<-wbar<-vector(); p[1]<-p0
		f11<-p[1]^2; f12<-2*p[1]*(1-p[1]); f22<-(1-p[1])^2
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
					main="frequency of A")
				else lines(ii,p[ii],type="l")
			} else if(show=="q"){
				if(i==2) plot(1:i,1-p,type="l",xlim=c(0,time),ylim=c(0,1),xlab="time",
					ylab="q",main="frequency of a")
				else lines(ii,p[ii],type="l")
			} else if(show=="fitness"){
				if(i==2) plot(1:i,wbar,type="l",xlim=c(0,time),ylim=c(0,1),xlab="time",
					main="mean fitness")
				else lines(ii,wbar[ii],type="l")
			} else if(show=="cobweb"){
				lines(c(p[i-1],p[i-1]),c(p[i-1],p[i]))
				lines(c(p[i-1],p[i]),c(p[i],p[i]))
			} else {
				message("not a recognized option")
				break
			}
			dev.flush()
			Sys.sleep(pause)
		}
	}
}

sexratio<-function(p0=0.01,time=100,show="p",pause=0){
	p<-wbar<-vector(); p[1]<-p0
	fm<-p[1]^2+p[1]*(1-p[1]); ff<-(1-p[1])^2+p[1]*(1-p[1])
	wm<-0.5/fm; wf<-0.5/ff
	wbar[1]<-fm*wm+wf*ff
	for(i in 2:time){
		p[i]<-p[i-1]
		w<-c(wm,(wm+wf)/2,wf)
		p[i]<-(p[i]^2*w[1]+p[i]*(1-p[i])*w[2])/wbar[i-1]
		fm<-p[i]^2+p[i]*(1-p[i]); ff<-(1-p[i])^2+p[i]*(1-p[i])
		wm<-0.5/fm; wf<-0.5/ff
		wbar[i]<-fm*wm+wf*ff
		ii<-(i-1):i
		if(show=="p"){
			if(i==2) plot(1:i,p,type="l",xlim=c(0,time),ylim=c(0,1),xlab="time",
				main="frequency of A")
			else lines(ii,p[ii],type="l")
		} else if(show=="fitness"){
			if(i==2) plot(1:i,wbar,type="l",xlim=c(0,time),ylim=c(0,1),xlab="time",
				main="mean fitness")
			else lines(ii,wbar[ii],type="l")
		}
		dev.flush()
	}
}

mutation.selection<-function(p0=1.0,w=c(1,1),u=0.001,time=100,show="q",pause=0,ylim=c(0,1)){
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

founder.event<-function(p0=0.5,Ne=1000,Nf=10,ttime=100,etime=50,show="p"){
	genotypes<-matrix(sample(c(rep(1,round(2*p0*Ne)),rep(0,2*Ne-round(2*p0*Ne)))),Ne,2)
	p<-v<-vector()
	p[1]<-mean(genotypes)
	v[1]<-p[1]*(1-p[1])
	for(i in 2:ttime){
		new.gen<-matrix(NA,Ne,2)
		if(i==etime){
			founder<-genotypes[sample(1:Ne,size=Nf),]
			for(j in 1:Ne) new.gen[j,]<-sample(founder,size=2)
		} else for(j in 1:Ne) new.gen[j,]<-sample(genotypes,size=2)
		genotypes<-new.gen
		p[i]<-mean(genotypes)
		v[i]<-p[i]*(1-p[i])
	}
	if(show=="p"){ 
		plot(1:ttime,p,type="l",main="frequency of A",xlab="time",ylim=c(0,1))
		lines(c(etime,etime),c(0,1),lty=2)
	} else if(show=="var"){ 
		plot(1:ttime,v,type="l",main="genetic variation",xlab="time",ylim=c(0,0.3))
		lines(c(etime,etime),c(0,0.3),lty=2)
	}
}

