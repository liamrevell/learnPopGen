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
	if(is.null(colors)) colors<-rainbow(n=n)
	popn<-matrix(NA,ngen+1,n)
	parent<-matrix(NA,ngen,n)
	popn[1,]<-1:n
	for(i in 1:ngen){
		parent[i,]<-sort(sample(1:n,replace=TRUE))
		popn[i+1,]<-popn[i,parent[i,]]
	}
	plot.new()
	par(mar=c(2.1,4.1,2.1,1.1))
	plot.window(xlim=c(0.5,n+0.5),ylim=c(ngen,0))
	axis(2)
	title(ylab="time (generations)")
	cx.pt<-2*25/max(n,ngen)
	points(1:n,rep(0,n),bg=colors,pch=21,cex=cx.pt)
	for(i in 1:ngen){
		dev.hold()
		for(j in 1:n){
			lines(c(parent[i,j],j),c(i-1,i),lwd=2,
				col=colors[popn[i+1,j]])
		}
		points(1:n,rep(i-1,n),bg=colors[popn[i,]],pch=21,
			cex=cx.pt)
		points(1:n,rep(i,n),bg=colors[popn[i+1,]],pch=21,
			cex=cx.pt)
		dev.flush()
		Sys.sleep(sleep)
	}
}
