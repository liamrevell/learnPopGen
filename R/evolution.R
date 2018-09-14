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
## written by Liam J. Revell 2018

phenotype.freq<-function(nloci=6,p=NULL,effect=1/nloci){
	if(is.null(p)) p<-rep(0.5,nloci)
	genotypes<-t(apply(cbind(p,1-p),1,hardy.weinberg))
	COMBN<-permutations(n=3,r=nloci,set=T,repeats.allowed=T)
	PHEN<--rowSums(COMBN-2)*effect
	FREQ<-rep(1,nrow(COMBN))
	for(i in 1:nrow(COMBN)){
		for(j in 1:nloci)	FREQ[i]<-FREQ[i]*genotypes[j,COMBN[i,j]]
	}
	phen<-unique(PHEN)
	freq<-rep(0,length(phen))
	for(i in 1:length(phen)) freq[i]<-sum(FREQ[which(PHEN==phen[i])])
	plot(phen,freq,type="b",pch=21,bg="grey",
		xlab="phenotypic trait value",ylab="relative frequency",
		cex=1.5,xlim=range(phen)+c(-0.5,0.5)*effect,ylim=c(0,max(freq)))
	for(i in 1:length(phen))
		rect(phen[i]-0.4*effect,0,phen[i]+0.4*effect,
		freq[i],border="grey",
		col=make.transparent("blue",0.2))
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
## written by Liam J. Revell 2018

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
}

drift.selection<-function(p0=0.5,Ne=100,w=c(1,1,1),ngen=400,nrep=10,
	colors=NULL,...){
	if(is.null(colors)) colors<-rainbow(nrep)
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
		xlab="time (generations)",ylab="f(A)")
	if(nrep>1)
		nulo<-mapply(lines,x=replicate(nrep-1,0:ngen,simplify=FALSE),
			y=p[2:nrep],col=colors[2:nrep],lwd=2)
	invisible(p)
}

msd<-function(p0=c(0.5,0.5),Ne=c(100,100),
	w=list(c(1,1,1),c(1,1,1)),m=c(0.01,0.01),ngen=400,
	colors=c("red","blue"),...){
	if(hasArg(show.legend)) show.legend=list(...)$show.legend
	else show.legend<-TRUE
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
		xlab="time (generations)",ylab="f(A)")
	lines(x=0:ngen,y=p[[2]],col=colors[2],lwd=2)
	if(show.legend) legend(x="topright",legend=1:2,lty=1,col=colors,
		lwd=2,bg=make.transparent("white",0.8))
	invisible(p)
}

clt<-function(nvar=1,nobs=1000,df=c("normal","uniform","exponential"),theta=1,
	breaks="Sturges"){
	df<-df[1]
	foo<-if(df=="normal") function(n,theta){
		rnorm(n,sd=sqrt(theta))
	} else if(df=="uniform") function(n,theta){
		runif(n,0,theta)
	} else if(df=="exponential") function(n,theta){
		rexp(n,rate=theta)
	}
	X<-matrix(foo(nvar*nobs,theta),nobs,nvar)
	x<-rowSums(X)
	if(is.integer(breaks)) breaks<-seq(min(x),max(x),by=diff(range(x))/(breaks-1))
	obj<-hist(x,border="darkgrey",col=phytools::make.transparent("blue",0.1),
		main=paste("CLT: the sum of",nvar,df,"distribution(s)"),
		breaks=breaks)
	lines(obj$mids,obj$counts,type="b",pch=21,bg="grey")
}
