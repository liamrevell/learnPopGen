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

## function to compute relative frequencies of a phenotypic trait
## for a polygenic trait under a simple additive genetic model
## written by Liam J. Revell 2018

phenotype.freq<-function(nloci=10,p=NULL,effect=1){
	if(is.null(p)) p<-rep(0.5,nloci)
	genotypes<-t(apply(cbind(p,1-p),1,hardy.weinberg))
	COMBN<-permutations(n=3,r=nloci,set=T,repeats.allowed=T)
	PHEN<-rowSums(COMBN-2)*effect
	FREQ<-rep(1,nrow(COMBN))
	for(i in 1:nrow(COMBN)){
		for(j in 1:nloci)	FREQ[i]<-FREQ[i]*genotypes[j,COMBN[i,j]]
	}
	phen<-unique(PHEN)
	freq<-rep(0,length(phen))
	for(i in 1:length(phen)) freq[i]<-sum(FREQ[which(PHEN==phen[i])])
	plot(phen,freq,type="b",pch=21,bg="grey",
		xlab="phenotypic trait value",ylab="relative frequency",
		cex=1.5,xlim=range(phen)+c(-0.5,0.5)*effect)
	for(i in 1:length(phen))
		rect(phen[i]-0.4*effect,0,phen[i]+0.4*effect,
		freq[i],border="grey",
		col=make.transparent("blue",0.2))
}

