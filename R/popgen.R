## functions by Liam Revell 2012 (some small updates 2017, 2018, 2019)

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
	if(hasArg(main)) main<-list(...)$main
	else main<-NULL
	if(hasArg(xlim)) xlim<-list(...)$xlim
	else xlim<-NULL
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
		plot(p,wbar,type="l",xlim=xlim,ylim=c(0,1.1*max(w)),
			main=if(is.null(main)) expression(paste("mean fitness (",
			bar(w),")",sep="")) else main,ylab=expression(bar(w)),col=color)
		if(equil) abline(v=eq,lty="dotted")
	}
	else if(show=="deltap"){
		p<-0:100/100
		wbar<-p^2*w[1]+2*p*(1-p)*w[2]+(1-p)^2*w[3]
		deltap<-(p/wbar)*(p*w[1]+(1-p)*w[2]-wbar)
		plot(p,deltap,type="l",xlim=xlim,main=if(is.null(main)) 
			expression(paste(Delta,"p as a function of p",sep="")) else main,
			ylab=expression(paste(Delta,"p",sep="")),col=color)
		lines(c(0,1),c(0,0),lty=2)
		if(equil) abline(v=eq,lty="dotted")
	} else {
		if(show=="cobweb"){
			p<-0:100/100
			wbar<-p^2*w[1]+2*p*(1-p)*w[2]+(1-p)^2*w[3]
			p2<-(p/wbar)*(p*w[1]+(1-p)*w[2]-wbar)+p
			plot(p,p2,type="l",xlim=xlim,xlab=expression(p[t]),ylab=expression(p[t+1]),
				main=if(is.null(main)) expression(paste(p[t+1]," as a function of ",
				p[t],sep="")) else main,col=color)
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
				if(i==2 && !add) plot(1:i,p,type="l",xlim=if(is.null(xlim)) c(0,time) else xlim,
					ylim=c(0,1),xlab="time",main=if(is.null(main)) "frequency of A" else main,
					col=color,lwd=lwd,lty=lty)
				else lines(ii,p[ii],type="l",col=color,lwd=lwd,lty=lty)
			} else if(show=="q" && !add){
				if(i==2) plot(1:i,1-p,type="l",xlim=if(is.null(xlim)) c(0,time) else xlim,
					ylim=c(0,1),xlab="time",ylab="q",
					main=if(is.null(main)) "frequency of a" else main,
					col=color,lwd=lwd,lty=lty)
				else lines(ii,1-p[ii],type="l",col=color,lwd=lwd,lty=lty)
			} else if(show=="fitness" && !add){
				if(i==2) plot(1:i,wbar,type="l",xlim=if(is.null(xlim)) c(0,time) else xlim,
					ylim=c(0,1.1*max(w)),
					xlab="time",main=if(is.null(main)) 
					expression(paste("mean fitness (",bar(w),")",sep="")) else main,
					ylab=expression(bar(w)),col=color)
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
	object<-list(p0=p0,w=w,time=time,
		p=if(show%in%c("surface","deltap")) NULL else p,
		equilibrium=eq)
	class(object)<-"selection"
	invisible(object)
}

print.selection<-function(x,...){
	cat("\nObject of class \"selection\" normally consisting of the expected allele")
	cat("\nfrequencies through time under frequency independent selection with genotype")
	cat("\nfitnesses as follows:\n")
	cat(paste("    W(AA) =",round(x$w[1],2),"\n"))
	cat(paste("    W(Aa) =",round(x$w[2],2),"\n"))
	cat(paste("    W(aa) =",round(x$w[3],2),"\n"))
	cat("\nTo plot enter plot(\'object_name\') at the command line\ninterface.\n\n")
}

plot.selection<-function(x,...){
	if(hasArg(color)) color<-list(...)$color
	else color<-"black"
	if(hasArg(lwd)) lwd<-list(...)$lwd
	else lwd<-2
	if(hasArg(lty)) lty<-list(...)$lty
	else lty<-"solid"
	if(hasArg(main)) main<-list(...)$main
	else main<-NULL
	if(hasArg(xlim)) xlim<-list(...)$xlim
	else xlim<-NULL
	if(hasArg(show)) show<-list(...)$show
	else show<-"p"
	if(hasArg(type)) type<-list(...)$type
	else type<-"l"
	w<-x$w
	eq<-x$equilibrium
	if(show%in%c("surface","deltap")){
		p<-0:100/100
		wbar<-p^2*w[1]+2*p*(1-p)*w[2]+(1-p)^2*w[3]
		if(show=="surface"){
			plot(p,wbar,type="l",xlim=xlim,ylim=c(0,1.1*max(w)),
				main=if(is.null(main)) expression(paste("mean fitness (",
				bar(w),")",sep="")) else main,ylab=expression(bar(w)),
				col=color)
			abline(v=eq,lty="dotted")
		} else if(show=="deltap"){
			deltap<-(p/wbar)*(p*w[1]+(1-p)*w[2]-wbar)
			plot(p,deltap,type="l",xlim=xlim,main=if(is.null(main)) 
				expression(paste(Delta,"p as a function of p",sep="")) else main,
				ylab=expression(paste(Delta,"p",sep="")),col=color)
			lines(c(0,1),c(0,0),lty=2)
			abline(v=eq,lty="dotted")
		}
	} else if(show%in%c("p","q","fitness","cobweb")&&!is.null(x$p)){
		if(show=="cobweb"){
			p<-0:100/100
			wbar<-p^2*w[1]+2*p*(1-p)*w[2]+(1-p)^2*w[3]
			p2<-(p/wbar)*(p*w[1]+(1-p)*w[2]-wbar)+p
			plot(p,p2,type="l",xlim=xlim,xlab=expression(p[t]),ylab=expression(p[t+1]),
				main=if(is.null(main)) expression(paste(p[t+1]," as a function of ",
				p[t],sep="")) else main,col=color)
			lines(c(0,1),c(0,1),lty=2)
			abline(v=eq,lty="dotted")
			abline(h=eq,lty="dotted")
		}
		p<-x$p
		wbar<-x$wbar
		if(show=="p"){
			plot(1:x$time,p,type=type,xlim=if(is.null(xlim)) c(0,x$time) else xlim,
				ylim=c(0,1),xlab="time",main=if(is.null(main)) "frequency of A" else main,
				col=color,lwd=lwd,lty=lty)
		} else if(show=="q"){
			plot(1:x$time,1-p,type=type,xlim=if(is.null(xlim)) c(0,x$time) else xlim,
				ylim=c(0,1),xlab="time",ylab="q",
				main=if(is.null(main)) "frequency of a" else main,
				col=color,lwd=lwd,lty=lty)
		} else if(show=="fitness"){
			plot(1:x$time,wbar,type="l",xlim=if(is.null(xlim)) c(0,x$time) else xlim,
				ylim=c(0,1.1*max(w)),
				xlab="time",main=if(is.null(main)) 
				expression(paste("mean fitness (",bar(w),")",sep="")) else main,
				ylab=expression(bar(w)),col=color)
		} else if(show=="cobweb"){
			for(i in 2:x$time){
				lines(c(p[i-1],p[i-1]),c(p[i-1],p[i]),col=color)
				lines(c(p[i-1],p[i]),c(p[i],p[i]),col=color)
			}
		}
		if(show=="p") abline(h=eq,lty="dotted")
		if(show=="q") abline(h=1-eq,lty="dotted")
		if(show=="fitness") abline(h=eq^2*w[1]+2*eq*(1-eq)*w[2]+(1-eq)^2*w[3],
			lty="dotted")
	} else cat("\nThis option not available for your input object as\ncurrently configured.\n\n")
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
	object<-list(p0=p0,
		s=s,
		time=time,
		p=if(show%in%c("surface","deltap")) NULL else p)
	class(object)<-"freqdep"
	invisible(object)
}

print.freqdep<-function(x,...){
	cat("\nAn object of class \"freqdep\" normally consisting of the expected")
	cat("\ngene frequency of the A allele through time under the following")
	cat("\nfrequency dependent selection model:\n")
	cat("    w(AA)=1-3*f(Aa)+3*f(aa)\n")
	if(x$s>0)
		cat(paste("    w(Aa)=1-",round(x$s,2),"*f(Aa)\n",sep=""))
	else
		cat(paste("    w(Aa)=1+",round(-x$s,2),"*f(Aa)\n",sep=""))
	cat("    w(aa)=1-3*f(Aa)+3*f(AA)\n\n")
}

plot.freqdep<-function(x,...){
	if(hasArg(type)) type<-list(...)$type
	else type<-"l"
	if(hasArg(lwd)) lwd<-list(...)$lwd
	else lwd<-1
	if(is.null(x$p)) cat("\nNo allele frequencies in \"freqdep\" object.\n\n")
	else {
		if(hasArg(color)) color<-list(...)$color
		else color<-par()$fg
		plot(1:length(x$p),x$p,type=type,xlim=c(0,x$time),ylim=c(0,1),
			xlab="time",ylab="p",main="frequency of A",col=color,lwd=lwd)
	}
}
	
sexratio<-function(p0=0.01,time=40,show="p",pause=0,sex.Aa=c(0.5,0.5)){
	p<-f.AA<-f.Aa<-f.aa<-fm<-ff<-wf<-wm<-wbar<-vector()
	p[1]<-p0
	f.AA[1]<-p[1]^2
	f.Aa[1]<-2*p[1]*(1-p[1])
	f.aa[1]<-(1-p[1])^2
	fm[1]<-f.AA[1]+sex.Aa[1]*f.Aa[1]
	ff[1]<-sex.Aa[2]*f.Aa[1]+f.aa[1]
	wm[1]<-0.5/fm[1]
	wf[1]<-0.5/ff[1]
	wbar[1]<-fm[1]*wm[1]+wf[1]*ff[1]
	t<-0:time
	for(i in 2:length(t)){
		## M(AA) x F(Aa)
		p_AA_Aa=f.AA[i-1]*f.Aa[i-1]*sex.Aa[2]
		## M(AA) x F(aa)
		p_AA_aa=f.AA[i-1]*f.aa[i-1]
		## M(Aa) x F(Aa)
		p_Aa_Aa=f.Aa[i-1]*sex.Aa[1]*f.Aa[i-1]*sex.Aa[2]
		## M(Aa) + F(aa)
		p_Aa_aa=f.Aa[i-1]*sex.Aa[1]*f.aa[i-1]
		## normalize
		sump<-p_AA_Aa+p_AA_aa+p_Aa_Aa+p_Aa_aa
		p_AA_Aa<-p_AA_Aa/sump
		p_AA_aa<-p_AA_aa/sump
		p_Aa_Aa<-p_Aa_Aa/sump
		p_Aa_aa<-p_Aa_aa/sump
		f.AA[i]<-f.Aa[i]<-f.aa[i]<-0
		## from M(AA) x F(Aa)
		f.AA[i]<-f.AA[i]+0.5*p_AA_Aa
		f.Aa[i]<-f.Aa[i]+0.5*p_AA_Aa
		## from M(Aa) x F(aa)
		f.Aa[i]<-f.Aa[i]+p_AA_aa
		## from M(Aa) x F(Aa)
		f.AA[i]<-f.AA[i]+0.25*p_Aa_Aa
		f.Aa[i]<-f.Aa[i]+0.5*p_Aa_Aa
		f.aa[i]<-f.aa[i]+0.25*p_Aa_Aa
		## from M(Aa) x F(aa)
		f.Aa[i]<-f.Aa[i]+0.5*p_Aa_aa
		f.aa[i]<-f.aa[i]+0.5*p_Aa_aa
		## compute frequencies
		p[i]<-f.AA[i]+0.5*f.Aa[i]
		fm[i]<-f.AA[i]+sex.Aa[1]*f.Aa[i]
		ff[i]<-sex.Aa[2]*f.Aa[i]+f.aa[i]
		## compute fitnesses
		wm[i]<-0.5/fm[i]
		wf[i]<-0.5/ff[i]
		wbar[i]<-fm[i]*wm[i]+wf[i]*ff[i]
	}
	for(i in 2:length(t)){
		ii<-c(i-1,i)
		if(show=="p"){
			if(i==2) plot(t[1:i],p[1:i],type="l",xlim=c(0,time),ylim=c(0,1),xlab="time",
				ylab="p",main="frequency of A",bty="l",col="darkgrey",lwd=2)
			else lines(t[ii],p[ii],type="l",lwd=2,col="darkgrey")
		} else if(show=="sex-ratio"){
			if(i==2){ 
				plot(t[1:i],fm[1:i],type="l",xlim=c(0,time),ylim=c(0,1),xlab="time",
					main="frequency of each sex",col=make.transparent("blue",0.5),
					lwd=2,ylab="relative frequencies of each sex",bty="l")
				lines(t[1:i],ff[1:i],lwd=2,col=make.transparent("red",0.5))
				legend(x="topright",c("males","females"),lwd=2,
					col=c(make.transparent("blue",0.5),
					make.transparent("red",0.5)))
			} else { 
				lines(t[ii],fm[ii],type="l",lwd=2,col=make.transparent("blue",0.5))
				lines(t[ii],ff[ii],type="l",lwd=2,col=make.transparent("red",0.5))
			}
		} else if(show=="fitness"){
			if(i==2){
				plot(t[1:i],wbar[1:i],type="l",xlim=c(0,time),ylim=c(min(wbar,wm,wf),max(wbar,wm,wf)),
					xlab="time",main=expression(paste("mean fitness (",bar(w),")",sep="")),
					ylab=expression(bar(w)),lwd=2,col=make.transparent("grey",1/3),
					log="y",bty="l")
				lines(t[1:i],wm[1:i],lwd=2,col=make.transparent("blue",1/3))
				lines(t[1:i],wf[1:i],lwd=2,col=make.transparent("red",1/3))
				legend(x="topright",c("average","males","females"),lwd=2,
					col=c(make.transparent("grey",1/3),
					make.transparent("blue",1/3),
					make.transparent("red",1/3)))
			} else { 
				lines(t[ii],wbar[ii],lwd=2,col=make.transparent("grey",1/3))
				lines(t[ii],wm[ii],lwd=2,col=make.transparent("blue",1/3))
				lines(t[ii],wf[ii],lwd=2,col=make.transparent("red",1/3))
			} 
		} else if(show=="genotypes"){
			if(i==2){
				plot(t[1:i],f.AA[1:i],type="l",xlim=c(0,time),ylim=c(0,1),xlab="time",
					main="frequency of each sex",col=make.transparent("red",0.5),
					lwd=2,ylab="relative frequency of each genotype",bty="l")
				lines(t[1:i],f.Aa[1:i],lwd=2,col=make.transparent("purple",0.5))
				lines(t[1:i],f.aa[1:i],lwd=2,col=make.transparent("blue",0.5))
				legend(x="topright",c("AA","Aa","aa"),lwd=2,
					col=c(make.transparent("red",0.5),
					make.transparent("purple",0.5),make.transparent("blue",0.5)))
			} else {
				lines(t[ii],f.AA[ii],lwd=2,col=make.transparent("red",0.5))
				lines(t[ii],f.Aa[ii],lwd=2,col=make.transparent("purple",0.5))
				lines(t[ii],f.aa[ii],lwd=2,col=make.transparent("blue",0.5))
			}
		}
		Sys.sleep(pause)
		dev.flush()
	}
	object<-list(p=p,fm=fm,ff=ff,
		f.AA=f.AA,f.Aa=f.Aa,f.aa=f.aa,
		wm=wm,wf=wf,wbar=wbar,
		sex.Aa=sex.Aa)
	class(object)<-"sexratio"
	invisible(object)
}

print.sexratio<-function(x,...){
	cat("\nObject of class \"sexratio\" consisting of the expected frequency of the two")
	cat("\nalternative alleles (A & a) at a sex-determination locus in which heterozygotes")
	cat(paste("\nare male with a probability of ",round(x$sex.Aa[1],2),
		"; and female with a probability of ",round(x$sex.Aa[2],2),".\n",sep=""))
	cat("\nTo plot enter plot(\'object_name\') at the command line\ninterface.\n\n")
}

plot.sexratio<-function(x,...){
	if(hasArg(show)) show<-list(...)$show
	else show<-"p"
	if(hasArg(lwd)) lwd<-list(...)$lwd
	else lwd<-2
	if(hasArg(type)) type<-list(...)$type
	else type<-"l"
	t<-1:length(x$p)-1
	if(show=="p"){
		plot(t,x$p,type=type,xlim=c(0,max(t)),ylim=c(0,1),xlab="time",
			ylab="p",main="frequency of sex-determination allele A",lwd=lwd,bty="l",col="darkgrey")
	} else if(show=="sex-ratio"){
		plot(t,x$fm,type=type,xlim=c(0,max(t)),ylim=c(0,1),xlab="time",
			main="frequency of each sex",col=make.transparent("blue",0.5),
			ylab="relative frequencies of each sex",lwd=lwd,bty="l")
		lines(t,x$ff,lwd=lwd,col=make.transparent("red",0.5),
			type=type)
		legend(x="topright",c("males","females"),lwd=lwd,
			col=c(make.transparent("blue",0.5),
			make.transparent("red",0.5)))
	} else if(show=="fitness"){
		plot(t,x$wbar,type=type,xlim=c(0,max(t)),ylim=c(min(x$wbar,x$wm,x$wf),
			max(x$wbar,x$wm,x$wf)),xlab="time",main=expression(paste("mean fitness (",bar(w),
			")",sep="")),ylab=expression(bar(w)),lwd=lwd,col=make.transparent("grey",1/3),
			log="y",bty="l")
		lines(t,x$wm,lwd=lwd,col=make.transparent("blue",1/3),type=type)
		lines(t,x$wf,lwd=lwd,col=make.transparent("red",1/3),type=type)
		legend(x="topright",c("average","males","females"),lwd=lwd,
			col=c(make.transparent("grey",1/3),
			make.transparent("blue",1/3),
			make.transparent("red",1/3)))
	} else if(show=="genotypes"){
		plot(t,x$f.AA,type=type,xlim=c(0,max(t)),ylim=c(0,1),xlab="time",
			main="frequency of each genotype",col=make.transparent("red",0.5),
			lwd=lwd,ylab="relative frequency of each genotype",bty="l")
		lines(t,x$f.Aa,lwd=lwd,col=make.transparent("purple",0.5),type=type)
		lines(t,x$f.aa,lwd=lwd,col=make.transparent("blue",0.5),type=type)
		legend(x="topright",c("AA","Aa","aa"),lwd=lwd,
			col=c(make.transparent("red",0.5),
			make.transparent("purple",0.5),make.transparent("blue",0.5)))
	}
}

mutation.selection<-function(p0=1.0,w=c(1,0),u=0.001,time=100,show="q",
	pause=0,ylim=c(0,1)){
	p<-wbar<-vector()
	p[1]<-p0
	wbar[1]<-p[1]^2*1.0+2*p[1]*(1-p[1])*w[1]+(1-p[1])^2*w[2]
	for(i in 2:time){
		p[i]<-p[i-1]
		p[i]<-(1-u)*(p[i]^2*1.0+p[i]*(1-p[i])*w[1])/wbar[i-1]
		wbar[i]<-p[i]^2*1.0+2*p[i]*(1-p[i])*w[1]+(1-p[i])^2*w[2]
		ii<-(i-1):i**
		if(show=="p"){
			if(i==2) plot(1:i,p,type="l",xlim=c(0,time),ylim=c(0,1),xlab="time",
				main="frequency of A")
			else lines(ii,p[ii],type="l")
		} else if(show=="q"){
			if(i==2) plot(1:i,1-p,type="l",xlim=c(0,time),ylim=ylim,xlab="time",
				ylab="q",main="frequency of a")
			else lines(ii,1-p[ii],type="l")
		} else if(show=="fitness"){
			if(i==2) plot(1:i,wbar,type="l",xlim=c(0,time),ylim=c(0,1),
				xlab="time",main=expression(paste("mean fitness (",
				bar(w),")",sep="")),ylab=expression(bar(w)))
			else lines(ii,wbar[ii],type="l")
		} else {
			message("not a recognized option")
			break
		}
		dev.flush()
		Sys.sleep(pause)
	}
	object<-list(p0=p0,w=w,u=u,time=time,
		p=p,wbar=wbar)
	class(object)<-"mutation.selection"
	invisible(object)
}

print.mutation.selection<-function(x,...){
	cat("\nObject of class \"mutation.selection\" containing the exected frequency\n")
	cat("of the A allele through time under mutation & selection, as specified by\n")
	cat("the user.\n\n")
	cat("\nTo plot enter plot(\'object_name\') at the command line interface.\n\n")
}

plot.mutation.selection<-function(x,...){
	if(hasArg(show)) show<-list(...)$show
	else show<-"q"
	if(hasArg(color)) color<-list(...)$color
	else color<-"blue"
	if(hasArg(lwd)) lwd<-list(...)$lwd
	else lwd<-2
	if(hasArg(ylim)) ylim<-list(...)$ylim
	else ylim<-c(0,1)
	if(hasArg(type)) type<-list(...)$type
	else type<-"l"
	if(show=="p"){
		plot(1:x$time,x$p,type=type,xlim=c(0,x$time),ylim=ylim,xlab="time",
			main="frequency of A",col=color,lwd=lwd)
	} else if(show=="q"){
		plot(1:x$time,1-x$p,type=type,xlim=c(0,x$time),ylim=ylim,xlab="time",
			ylab="q",main="frequency of a",col=color,lwd=lwd)
	} else if(show=="fitness"){
		plot(1:x$time,x$wbar,type=type,xlim=c(0,x$time),ylim=ylim,
			xlab="time",main="mean fitness",col=color,lwd=lwd,ylab="fitness")
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
	attr(p,"p0")<-p0
	attr(p,"Ne")<-Ne
	class(p)<-"genetic.drift"
	invisible(p)
}

repMatrix<-function(X,times){
	Z<-list()
	for(i in 1:times) Z[[i]]<-X
	return(Z)
}

print.genetic.drift<-function(x,...){
	cat("\nObject of class \"genetic.drift\" consisting of allele frequencies")
	cat(paste("\nfrom ",ncol(x)," independent genetic drift simulation(s) each initiated",sep=""))
	cat(paste("\nwith a starting allele frequency of ",round(attr(x,"p0"),2)," and an effective population",sep=""))
	cat(paste("\nsize of ",attr(x,"Ne"),".\n\n",sep=""))
	cat("\nTo plot enter plot(\'object_name\') at the command line interface.\n\n")
}

plot.genetic.drift<-function(x,...){
	if(hasArg(show)) show<-list(...)$show
	else show<-"p"
	if(hasArg(pause)) pause<-list(...)$pause
	else pause<-0.05
	if(show=="p"){
		if(hasArg(colors)) colors<-list(...)$colors
		else colors<-rep("black",ncol(x))
		if(hasArg(lwd)) lwd<-list(...)$lwd
		else lwd<-1
		if(hasArg(type)) type<-list(...)$type
		else type<-"l"
		plot(1:nrow(x),x[,1],type=type,col=colors[1],
			lwd=lwd,main="frequency of A",xlab="time",
			ylab="p",ylim=c(0,1))
		if(attr(x,"p0")<=0.5) text(paste("N =",attr(x,"Ne"),
			sep=" "),x=0,y=1,pos=4)
		else text(paste("N =",attr(x,"Ne"),sep=" "),x=0,y=0,
			pos=4)
		for(i in 2:ncol(x)) lines(1:nrow(x),x[,i],col=colors[i],
			lwd=lwd,type=type)
	} else if(show=="genotypes"){
		AA<-unclass(x)^2
		Aa<-2*unclass(x)*(1-unclass(x))
		aa<-(1-unclass(x))^2
		for(i in 1:nrow(x)){
			X<-cbind(aa[i,],Aa[i,],AA[i,])
			colnames(X)<-c("aa","Aa","AA")
			barplot(X,ylim=c(0,1),main="genotype frequencies",beside=TRUE)
			Sys.sleep(pause)
		}
	} else if(show=="fixed"){
		for(i in 1:nrow(x)){
			fixedA<-mean(x[i,]==1)
			fixeda<-mean(x[i,]==0)
			barplot(c(fixeda,fixedA),ylim=c(0,1),names.arg=c("a","A"),
				main="populations fixed",ylab="frequency")
			Sys.sleep(pause)
		}
	} else if(show=="heterozygosity"){
		if(hasArg(lwd)) lwd<-list(...)$lwd
		else lwd<-1
		if(hasArg(type)) type<-list(...)$type
		else type<-"l"
		plot(1:nrow(x),2*attr(x,"p0")*(1-attr(x,"p0"))*(1-1/(2*attr(x,"Ne")))^(0:(nrow(x)-1)),
			lwd=2,col="red",main="heterozygosity",xlab="time",ylab="f(Aa)",
			ylim=c(0,0.6))
		lines(1:nrow(x),rowMeans(2*x*(1-x)),lwd=lwd,type=type)
	}
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
	class(p)<-"founder.event"
	attr(p,"Ne")<-Ne
	attr(p,"Nf")<-Nf
	attr(p,"ttime")<-ttime
	attr(p,"etime")<-etime
	invisible(p)
}

print.founder.event<-function(x,...){
	if(length(attr(x,"etime"))==1) duration<-1 else duration<-abs(diff(attr(x,"etime")))
	cat(paste("\nObject of class \"founder.event\" consisting of ",attr(x,"ttime"),sep=""))
	cat(paste("\ntotal generations, including a bottleneck of ",duration,sep=""))
	cat("\ngenerations.\n")
	cat(paste("The simulation had an effective population size of ",attr(x,"Ne"),"\n",sep=""))
	cat(paste("and a founder population size of ",attr(x,"Nf"),".\n",sep=""))
	cat("\nTo plot enter plot(\'object_name\') at the command line\ninterface.\n\n")
}

plot.founder.event<-function(x,...){
	if(hasArg(show)) show<-list(...)$show
	else show<-"p"
	if(hasArg(ltype)) ltype<-list(...)$ltype
	else ltype<-"s"
	Ne<-attr(x,"Ne")
	Nf<-attr(x,"Nf")
	ttime<-attr(x,"ttime")
	etime<-attr(x,"etime")
	if(show=="p")
		plot(1:ttime,x,main="frequency of A",xlab="time",ylim=c(0,1),type=ltype,
			ylab=expression(p[A]))
	else if(show=="var"){
		v<-x*(1-x)
		plot(1:ttime,v,main="genetic variation",xlab="time",ylim=c(0,0.3),type=ltype,
			ylab=expression(p[A]*(1-p[A])))
	}
	if(length(etime)==1) abline(v=etime,col=make.transparent("grey",0.5),lwd=6)
	else rect(xleft=etime[1],ybottom=par()$usr[3],
			xright=etime[length(etime)],ytop=par()$usr[4],
			col=make.transparent("grey",0.5),border=NA)
}
