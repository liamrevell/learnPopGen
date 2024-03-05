## function to compute multiallelic Hardy-Weinberg frequencies
## written by Liam J. Revell 2018

hardy.weinberg<-function(p=c(0.5,0.5),alleles=c("A","a"),as.matrix=FALSE){
	if(sum(p)!=1) p<-p/sum(p)
	nalleles<-max(length(p),length(alleles))
	if(nalleles!=length(p)) p<-rep(1/nalleles,nalleles)
	if(nalleles!=length(alleles)) {
		if(length(p)<=26) alleles<-LETTERS[1:length(p)]
		else alleles<-paste("A",1:length(p),sep="")
	}
	p<-setNames(p,alleles)
	HW<-as.matrix(p)%*%t(as.matrix(p))
	if(as.matrix) return(HW)
	else {
		hw<-vector()
		k<-1
		for(i in 1:length(alleles)){
			for(j in i:length(alleles)){
				hw[k]<-if(i==j) HW[i,j] else 2*HW[i,j]
				names(hw)[k]<-paste(alleles[c(i,j)],collapse="")
				k<-k+1
			}
		}
		return(hw)
	}
}

## compute multilocus Hardy-Weinberg in which each individual locus is biallelic
## written by Liam J. Revell 2018

multilocus.hw<-function(nloci=2,p=NULL){
	if(is.null(p)) p<-rep(0.5,nloci)
	if(length(p)!=nloci) nloci<-length(p)
	genotypes<-t(apply(cbind(p,1-p),1,hardy.weinberg))
	COMBN<-permutations(n=3,r=nloci,set=T,repeats.allowed=T)
	ALLELESp<-LETTERS[1:nloci]
	ALLELESq<-letters[1:nloci]
	GENOTYPE<-vector()
	FREQ<-rep(1,nrow(COMBN))
	for(i in 1:nrow(COMBN)){
		for(j in 1:nloci){
			FREQ[i]<-FREQ[i]*genotypes[j,COMBN[i,j]]
			gtype<-if(COMBN[i,j]==1) paste(rep(ALLELESp[j],2),collapse="")
				else if(COMBN[i,j]==2) paste(c(ALLELESp[j],ALLELESq[j]),collapse="")
				else if(COMBN[i,j]==3) paste(rep(ALLELESq[j],2),collapse="")
			GENOTYPE[i]<-if(j==1) gtype else paste(GENOTYPE[i],gtype,sep="")
		}
	}
	hw<-setNames(FREQ,GENOTYPE)
	hw
}

## function to compute relative frequencies of a phenotypic trait
## for a polygenic trait under a simple additive genetic model
## written by Liam J. Revell 2018, 2019

phenotype.freq<-function(nloci=6,p=NULL,effect=1/nloci){
	if(is.null(p)) p<-rep(0.5,nloci)
	if(length(effect)==1) effect<-rep(effect,nloci)
	else if(length(effect)>1&&length(effect)!=nloci){
		cat("The length of \'effect\' should match \'nloci\'.\n")
		cat("Setting to recycle first value of \'effect\'.\n")
		effect<-rep(effect[1],nloci)
	}
	genotypes<-t(apply(cbind(p,1-p),1,hardy.weinberg))
	COMBN<-permutations(n=3,r=nloci,set=T,repeats.allowed=T)
	PHEN<-round(-rowSums((COMBN-2)%*%effect),12)
	FREQ<-rep(1,nrow(COMBN))
	for(i in 1:nrow(COMBN)){
		for(j in 1:nloci)	FREQ[i]<-FREQ[i]*genotypes[j,COMBN[i,j]]
	}
	phen<-sort(unique(PHEN))
	if(length(phen)<0.2*length(PHEN)){
		type<-"discrete"
		freq<-rep(0,length(phen))
		for(i in 1:length(phen)) freq[i]<-sum(FREQ[which(PHEN==phen[i])])
		plot(phen,freq,type="b",pch=21,bg="grey",
			xlab="phenotypic trait value",ylab="relative frequency",
			cex=1.5,xlim=range(phen)+c(-0.5,0.5)*mean(effect),
			ylim=c(0,max(freq)))
		by<-min(rowSums(cbind(-phen[2:length(phen)-1],phen[2:length(phen)])))
		for(i in 1:length(phen))
			rect(phen[i]-0.4*by,0,phen[i]+0.4*by,
			freq[i],border="grey",
			col=make.transparent("blue",0.2))
	} else {
		type<-"continuous"
		h<-hist(PHEN,breaks=seq(min(PHEN),max(PHEN),
			by=diff(range(PHEN))/20),plot=FALSE)
		phen<-h$mids
		freq<-h$counts/sum(h$counts)
		plot(phen,freq,type="b",pch=21,bg="grey",
			xlab="phenotypic trait value",ylab="relative frequency",
			cex=1.5,xlim=range(phen)+c(-0.5,0.5)*mean(effect),
			ylim=c(0,max(freq)))
		by<-phen[2]-phen[1]
		for(i in 1:length(phen))
			rect(phen[i]-0.5*by,0,phen[i]+0.5*by,
				freq[i],border="grey",
				col=make.transparent("blue",0.2))
	}
	object<-list(phenotype=phen,frequency=freq,
		nloci=nloci,p=p,effect=effect,
		type=type)
	class(object)<-"phenotype.freq"
	invisible(object)
}

print.phenotype.freq<-function(x,...){
	cat("\nObject of class \"phenotype.freq\" containing the frequencies of each value")
	cat(paste("\nof a hypothetical polygenic trait determined by the additive effect of",x$nloci))
	cat("\ngenetic loci with the following frequencies,\n")
	cat(paste("  p(A): ",paste(round(x$p,2),collapse=", "),"\n",sep=""))
	cat("and the following additive effect of allelic subsitution,\n")
	cat(paste("  a: ",paste(round(x$effect,2),collapse=", "),"\n\n",sep=""))
	cat("\nTo plot enter plot(\'object_name\') at the command line interface.\n\n")
}

plot.phenotype.freq<-function(x,...){
	phen<-x$phenotype
	freq<-x$frequency
	if(x$type=="discrete"){
		plot(phen,freq,type="b",pch=21,bg="grey",
			xlab="phenotypic trait value",ylab="relative frequency",
			cex=1.5,xlim=range(phen)+c(-0.5,0.5)*mean(x$effect),
			ylim=c(0,max(freq)))
		by<-min(rowSums(cbind(-phen[2:length(phen)-1],
			phen[2:length(phen)])))
		for(i in 1:length(phen))
			rect(phen[i]-0.4*by,0,phen[i]+0.4*by,
			freq[i],border="grey",
			col=make.transparent("blue",0.2))
	} else {
		plot(phen,freq,type="b",pch=21,bg="grey",
			xlab="phenotypic trait value",ylab="relative frequency",
			cex=1.5,xlim=range(phen)+c(-0.5,0.5)*mean(x$effect),
			ylim=c(0,max(freq)))
		by<-phen[2]-phen[1]
		for(i in 1:length(phen))
			rect(phen[i]-0.5*by,0,phen[i]+0.5*by,
				freq[i],border="grey",
				col=make.transparent("blue",0.2))
	}
}

## function to conduct numerical analysis of phenotypic selection on 
## a polygenic quantitative trait
## written by Liam J. Revell 2018

phenotype.selection<-function(nloci=6,p=NULL,effect=1/nloci,beta=0.1,ngen=20,...){
	if(hasArg(sleep)) sleep<-list(...)$sleep
	else sleep<-0.1
	if(is.null(p)) p<-rep(0.5,nloci)
	COMBN<-permutations(n=3,r=nloci,set=T,repeats.allowed=T)
	PHEN<--rowSums(COMBN-2)*effect
	w<-beta*(PHEN-min(PHEN))+1
	het<-as.matrix(apply(COMBN,2,function(x) which(x==2)))
	homo<-as.matrix(apply(COMBN,2,function(x) which(x==1)))
	for(i in 1:ngen){
		genotypes<-t(apply(cbind(p,1-p),1,hardy.weinberg))
		FREQ<-rep(1,nrow(COMBN))
		for(j in 1:nrow(COMBN)){
			for(k in 1:nloci)	FREQ[j]<-FREQ[j]*genotypes[k,COMBN[j,k]]
		}
		phen<-unique(PHEN)
		freq<-rep(0,length(phen))
		for(j in 1:length(phen)) freq[j]<-sum(FREQ[which(PHEN==phen[j])])
		dev.hold()
		plot(phen,freq,type="b",pch=21,bg="grey",
			xlab="phenotypic trait value",ylab="relative frequency",
			cex=1.5,xlim=range(phen)+c(-0.5,0.5)*effect,ylim=c(0,max(freq)))
		for(j in 1:length(phen))
			rect(phen[j]-0.4*effect,0,phen[j]+0.4*effect,
			freq[j],border="grey",
			col=make.transparent("blue",0.2))
		dev.flush()
		FREQp<-FREQ*w/sum(FREQ*w)
		pp<-vector()
		for(j in 1:nloci) pp[j]<-sum(FREQp[het[,j]]/2)+sum(FREQp[homo[,j]])
		p<-pp
		Sys.sleep(sleep)
	}
}

## plot a coalescent genealogy
## written by Liam J. Revell 2018, 2019

coalescent.plot<-function(n=10,ngen=20,colors=NULL,...){
	if(hasArg(sleep)) sleep<-list(...)$sleep
	else sleep<-0.2
	if(hasArg(lwd)) lwd<-list(...)$lwd
	else lwd<-2
	if(hasArg(mar)) mar<-list(...)$mar
	else mar<-c(2.1,4.1,2.1,1.1)
	if(is.null(colors)){ 
		colors<-rainbow(n=n)
		if(hasArg(col.order)) col.order<-list(...)$col.order
		else col.order<-"sequential"
		if(col.order=="alternating"){
			if(n%%2==1) 
				ii<-as.vector(rbind(1:((n+1)/2),1:((n+1)/2)+(n+1)/2))
			else
				ii<-as.vector(rbind(1:(n/2),n/2+1:(n/2)))
			colors<-colors[ii]
		}
	}
	popn<-matrix(NA,ngen+1,n)
	parent<-matrix(NA,ngen,n)
	popn[1,]<-1:n
	for(i in 1:ngen){
		parent[i,]<-sort(sample(1:n,replace=TRUE))
		popn[i+1,]<-popn[i,parent[i,]]
	}
	plot.new()
	par(mar=mar)
	plot.window(xlim=c(0.5,n+0.5),ylim=c(ngen,0))
	axis(2)
	title(ylab="time (generations)")
	cx.pt<-2*25/max(n,ngen)
	points(1:n,rep(0,n),bg=colors,pch=21,cex=cx.pt)
	for(i in 1:ngen){
		dev.hold()
		for(j in 1:n){
			lines(c(parent[i,j],j),c(i-1,i),lwd=lwd,
				col=colors[popn[i+1,j]])
		}
		points(1:n,rep(i-1,n),col="grey",bg=colors[popn[i,]],pch=21,
			cex=cx.pt)
		points(1:n,rep(i,n),col="grey",bg=colors[popn[i+1,]],pch=21,
			cex=cx.pt)
		dev.flush()
		Sys.sleep(sleep)
	}
	object<-list(allele=popn,parent=parent)
	class(object)<-"coalescent.plot"
	invisible(object)
}

print.coalescent.plot<-function(x,...){
	cat("\nObject of class \"coalescent.plot\" consisting of a simulated process of allelic\n")
	cat(paste("coalescence over",nrow(x$allele)-1,"generations, in a population containing",
		ncol(x$allele),"individuals.\n\n"))
	cat("The object consists of:\n")
	cat(paste("  (1) a ",nrow(x$allele)," x ",ncol(x$allele),
		" numeric matrix containing the unique \'alleles\'\n",sep=""))
	cat(paste("      present in the population from time=0 to time=",nrow(x$allele)-1,".\n",sep=""))
	cat(paste("  (2) a ",nrow(x$parent)," x ",ncol(x$parent),
		" numeric matrix giving the parent/offspring\n",sep=""))
	cat(paste("      relationships across all ",nrow(x$parent),
		" generations of the simulation.\n\n",sep=""))
	cat("To re-plot type plot(object_name) at the command line.\n\n")
}

plot.coalescent.plot<-function(x,...){
	popn<-x$allele
	parent<-x$parent
	n<-ncol(popn)
	ngen<-nrow(parent)
	if(hasArg(sleep)) sleep<-list(...)$sleep
	else sleep<-0.2
	if(hasArg(lwd)) lwd<-list(...)$lwd
	else lwd<-2
	if(hasArg(mar)) mar<-list(...)$mar
	else mar<-c(2.1,4.1,2.1,1.1)
	if(hasArg(colors)) colors<-list(...)$colors
	else colors<-NULL
	if(is.null(colors)){ 
		colors<-rainbow(n=n)
		if(hasArg(col.order)) col.order<-list(...)$col.order
		else col.order<-"sequential"
		if(col.order=="alternating"){
			if(n%%2==1) 
				ii<-as.vector(rbind(1:((n+1)/2),1:((n+1)/2)+(n+1)/2))
			else
				ii<-as.vector(rbind(1:(n/2),n/2+1:(n/2)))
			colors<-colors[ii]
		}
	}
	plot.new()
	par(mar=mar)
	plot.window(xlim=c(0.5,n+0.5),ylim=c(ngen,0))
	axis(2)
	title(ylab="time (generations)")
	cx.pt<-2*25/max(n,ngen)
	points(1:n,rep(0,n),bg=colors,pch=21,cex=cx.pt)
	for(i in 1:ngen){
		dev.hold()
		for(j in 1:n){
			lines(c(parent[i,j],j),c(i-1,i),lwd=lwd,
				col=colors[popn[i+1,j]])
		}
		points(1:n,rep(i-1,n),col="grey",bg=colors[popn[i,]],pch=21,
			cex=cx.pt)
		points(1:n,rep(i,n),col="grey",bg=colors[popn[i+1,]],pch=21,
			cex=cx.pt)
		dev.flush()
		Sys.sleep(sleep)
	}
}

## function to simulate drift & selection

drift.selection<-function(p0=0.5,Ne=100,w=c(1,1,1),ngen=400,nrep=10,
	colors=NULL,...){
	if(is.null(colors)) colors<-rainbow(nrep)
	if(hasArg(xlabel)) xlabel=list(...)$xlabel
	else xlabel<-"Time (generations)"
	if(hasArg(ylabel)) ylabel=list(...)$ylabel
	else ylabel<-"f(A)"
	w<-(w/max(w))[3:1]
	gametes<-rep(0,2*Ne)
	if(p0>0) gametes[1:round(p0*2*Ne)]<-1
	gametes<-replicate(nrep,gametes,simplify=FALSE)
	p<-lapply(gametes,mean)
	for(i in 1:ngen){
		genotypes<-lapply(gametes,function(x) matrix(sample(x),
			length(x)/2,2))
		fitness<-lapply(genotypes,function(x,w) w[rowSums(x)+1],w=w)
		selected<-lapply(fitness,function(prob,x) 
			cbind(sample(x,prob=prob,replace=TRUE),
			sample(x,prob=prob,replace=TRUE)),x=Ne)
		copy<-replicate(nrep,matrix(sample(1:2,2*Ne,replace=TRUE),
			Ne,2),simplify=FALSE)
		gametes<-mapply(function(g,s,c) c(diag(g[s[,1],][,c[,1]]),
			diag(g[s[,2],][,c[,2]])),
			g=genotypes,s=selected,c=copy,SIMPLIFY=FALSE)
		for(j in 1:nrep) p[[j]][i+1]<-mean(gametes[[j]])
	}
	plot(0:ngen,p[[1]],type="l",col=colors[1],lwd=2,ylim=c(0,1),
		xlab=xlabel,ylab=ylabel)
	if(nrep>1)
		nulo<-mapply(lines,x=replicate(nrep-1,0:ngen,simplify=FALSE),
			y=p[2:nrep],col=colors[2:nrep],lwd=2)
	class(p)<-"drift.selection"
	attr(Ne,"Ne")<-Ne
	attr(p,"w")<-w[3:1]
	invisible(p)
}

print.drift.selection<-function(x,...){
	ngen<-length(x[[1]])-1
	nsim<-length(x)
	cat("\nObject of class \"drift.selection\" consisting of")
	cat(paste("\nallele frequencies through time from ",nsim," simulation(s) of\n",sep=""))
	cat(paste(ngen," generations of genetic drift & natural selection",sep=""))
	cat("\nwith normalized fitnesses of:\n")
	cat(paste("    W(AA) =",round(attr(x,"w")[1],2),"\n"))
	cat(paste("    W(Aa) =",round(attr(x,"w")[2],2),"\n"))
	cat(paste("    W(aa) =",round(attr(x,"w")[3],2),"\n"))
	cat("\nTo plot enter plot(\'object_name\') at the command line\ninterface.\n\n")
}

plot.drift.selection<-function(x,...){
	ngen<-length(x[[1]])-1
	nrep<-length(x)
	if(hasArg(colors)) colors<-list(...)$colors
	else colors<-rainbow(nrep)
	if(hasArg(type)) type<-list(...)$type
	else type<-"l"
	if(hasArg(lwd)) lwd<-list(...)$lwd
	else lwd<-2
	if(hasArg(xlabel)) xlabel=list(...)$xlabel
	else xlabel<-"Time (generations)"
	if(hasArg(ylabel)) ylabel=list(...)$ylabel
	else ylabel<-"f(A)"
	plot(0:ngen,x[[1]],type=type,col=colors[1],lwd=lwd,ylim=c(0,1),
		xlab=xlabel,ylab=ylabel)
	if(nrep>1)
		nulo<-mapply(lines,x=replicate(nrep-1,0:ngen,simplify=FALSE),
			y=x[2:nrep],col=colors[2:nrep],type=type,lwd=lwd)
}

as.data.frame.drift.selection<-function(x,...){
	ngen<-length(x[[1]])-1
	nsim<-length(x)
	x<-unclass(x)
	as.data.frame(x,row.names=0:ngen,
		col.names=paste("sim_",1:nsim,sep=""))
}

## function to simulate migration, selection, & drift

msd<-function(p0=c(0.5,0.5),Ne=c(100,100),
	w=list(c(1,1,1),c(1,1,1)),m=c(0.01,0.01),ngen=400,
	colors=c("red","blue"),...){
	if(hasArg(show.legend)) show.legend=list(...)$show.legend
	else show.legend<-TRUE
	if(hasArg(xlabel)) xlabel=list(...)$xlabel
	else xlabel<-"Time (generations)"
	if(hasArg(ylabel)) ylabel=list(...)$ylabel
	else ylabel<-"f(A)"
	w<-lapply(w,function(w) (w/max(w))[3:1])
	gametes<-lapply(Ne,function(Ne) rep(0,2*Ne))
	gametes<-mapply(function(p0,g,N){
		g[1:round(p0*2*N)]<-1
		g},p0=p0,g=gametes,N=Ne,SIMPLIFY=FALSE)
	p<-lapply(gametes,mean)
	for(i in 1:ngen){
		genotypes<-lapply(gametes,function(x) matrix(sample(x),
			length(x)/2,2))
		migrants<-mapply(function(N,m) which(runif(N)<=m),N=Ne,
			m=m,SIMPLIFY=FALSE)
		for(j in 1:length(genotypes)){
			to<-if(j==1) 2 else 1
			genotypes[[to]]<-rbind(genotypes[[to]],
				genotypes[[j]][migrants[[j]],])
		}
		for(j in 1:length(genotypes)){
			if(length(migrants[[j]])>0)
				genotypes[[j]]<-genotypes[[j]][-migrants[[j]],]
		}
		fitness<-mapply(function(x,w) w[rowSums(x)+1],x=genotypes,
			w=w,SIMPLIFY=FALSE)
		selected<-mapply(function(prob,N,Ne) 
			cbind(sample(N,Ne,prob=prob,replace=TRUE),
			sample(N,Ne,prob=prob,replace=TRUE)),prob=fitness,
			N=sapply(genotypes,nrow),Ne=Ne,SIMPLIFY=FALSE)
		copy<-lapply(Ne,function(Ne)
			matrix(sample(1:2,2*Ne,replace=TRUE),Ne,2))
		gametes<-mapply(function(g,s,c) c(diag(g[s[,1],][,c[,1]]),
			diag(g[s[,2],][,c[,2]])),
			g=genotypes,s=selected,c=copy,SIMPLIFY=FALSE)
		for(j in 1:2) p[[j]][i+1]<-mean(gametes[[j]])
	}
	plot(0:ngen,p[[1]],type="l",col=colors[1],lwd=2,ylim=c(0,1),
		xlab=xlabel,ylab=ylabel)
	lines(x=0:ngen,y=p[[2]],col=colors[2],lwd=2)
	if(show.legend) legend(x="topright",legend=1:2,lty=1,col=colors,
		lwd=2,bg=make.transparent("white",0.8))
	attr(p,"p0")<-p0
	attr(p,"Ne")<-Ne
	attr(p,"w")<-lapply(w,function(w) w[3:1])
	attr(p,"m")<-m
	class(p)<-"msd"
	invisible(p)
}

print.msd<-function(x,...){
	cat("\nObject of class \"msd\" containing the results from a numerical simulation")
	cat("\nof drift, selection, and migration within and between two populations.")
	cat("\n\nThe following parameters were used in the simulation:\n")
	cat(paste("\t",paste("Ne[",1:2,"]",sep="",collapse="\t"),"\n",sep=""))
	cat(paste("\t",paste(attr(x,"Ne"),collapse="\t")))
	cat("\n")
	cat("\tm[1->2]\tm[2->1]\n")
	cat(paste("\t",paste(attr(x,"m"),collapse="\t")))
	cat("\n")
	print(matrix(c(attr(x,"w")[[1]],attr(x,"w")[[2]]),2,3,byrow=TRUE,
		dimnames=list(c("pop'n 1","pop'n 2"),c("w(AA)","w(Aa)","w(aa)"))))
	cat("\n")
	cat("\nTo plot enter plot(\'object_name\') at the command line interface.\n\n")
}

plot.msd<-function(x,...){
	if(hasArg(show.legend)) show.legend=list(...)$show.legend
	else show.legend<-TRUE
	if(hasArg(colors)) colors<-list(...)$colors
	else colors<-c("red","blue")
	if(hasArg(lwd)) lwd<-list(...)$lwd
	else lwd<-2
	if(hasArg(type)) type<-list(...)$type
	else type<-"l"
	ngen<-length(x[[1]])-1
	plot(0:ngen,x[[1]],type=type,col=colors[1],lwd=lwd,ylim=c(0,1),
		xlab="time (generations)",ylab="f(A)")
	lines(x=0:ngen,y=x[[2]],col=colors[2],lwd=lwd,type=type)
	if(show.legend) legend(x="topright",legend=1:2,lty=1,col=colors,
		lwd=lwd,bg=make.transparent("white",0.8))
}

as.data.frame.msd<-function(x,...){
	ngen<-length(x[[1]])-1
	x<-unclass(x)
	as.data.frame(x,row.names=0:ngen,
		col.names=paste("sim_",1:nsim,sep=""))
}

## function to illustrate the central limit theorem (CLT)

clt<-function(nvar=1,nobs=1000,df=c("normal","uniform","exponential","binomial"),
	theta=NULL,breaks="Sturges",show=c("sum","mean")){
	df<-df[1]
	show<-show[1]
	if(is.null(theta))
		theta<-if(df%in%c("normal","uniform","exponential")) 1 else 0.5
	foo<-if(df=="normal") function(n,theta){
		rnorm(n,sd=sqrt(theta))
	} else if(df=="uniform") function(n,theta){
		runif(n,0,theta)
	} else if(df=="exponential") function(n,theta){
		rexp(n,rate=theta)
	} else if(df=="binomial") function(n,theta){
		rbinom(n,1,theta)
	}
	X<-matrix(foo(nvar*nobs,theta),nobs,nvar)
	x<-if(show=="sum") rowSums(X) else rowMeans(X)
	if(is.integer(breaks)) breaks<-seq(min(x),max(x),by=diff(range(x))/(breaks-1))
	obj<-hist(x,border="darkgrey",col=phytools::make.transparent("blue",0.1),
		main=paste("CLT: the",show,"of",nvar,df,"distribution(s)"),
		breaks=breaks)
	lines(obj$mids,obj$counts,type="b",pch=21,bg="grey")
	object<-list(dist=obj,data=X,df=df,show=show)
	class(object)<-"clt"
	invisible(object)
}

## S3 methods for clt

print.clt<-function(x,...){
	cat("\nObject of class \"clt\" consisting of:\n")
	cat(paste("  (1) ",ncol(x$data)," ",x$df,
		"ly distributed independent random variables, each with\n",sep=""))
	cat(paste("  (2)",nrow(x$data),"observations (stored in x$data), and\n"))
	cat(paste("  (3) a histogram giving the distribution of their observation-wise ",
		x$show,".\n\n",sep=""))
}

plot.clt<-function(x,...){
	plot(x$dist,border="darkgrey",col=phytools::make.transparent("blue",0.1),
		main=paste("CLT: the",x$show,"of",ncol(x$data),x$df,"distribution(s)"))
	lines(x$dist$mids,x$dist$counts,type="b",pch=21,bg="grey")
}

as.data.frame.clt<-function(x,...){
	x<-x$data
	as.data.frame(x,row.names=1:nrow(x),
		col.names=paste("V",1:ncol(x),sep=""))
}