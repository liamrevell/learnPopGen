# function performs individual based simulations of reproductive character displacement
# written by Liam J. Revell 2012

rcd<-function(nsp=3,nindivs=c(700,400,100),w_t=10,gen=c(500,500),figs="on",pf=100){
	# set simulation control conditions
	burnin<-gen[1]
	ngen<-gen[2]
	nloci.pref<-10
	nloci.trait<-10
	w.trait<-w_t
	w.pref<-1.0
	u<-0.02
	mutvar<-0.01
	# initiate population
	species<-list()
	for(i in 1:nsp){
		species[[i]]<-list()
		species[[i]]$individual<-list()
		for(j in 1:nindivs[i]){
			species[[i]]$individual[[j]]<-list()
			species[[i]]$individual[[j]]$genome.pref<-matrix(rep(0,nloci.pref*2),nloci.pref,2)
			species[[i]]$individual[[j]]$genome.trait<-matrix(rep(0,nloci.trait*2),nloci.trait,2)
			species[[i]]$individual[[j]]$pref<-sum(species[[i]]$individual[[j]]$genome.pref)
			species[[i]]$individual[[j]]$trait<-sum(species[[i]]$individual[[j]]$genome.trait)
			species[[i]]$individual[[j]]$sex<-r01() # assign random 0,1 sex
		}
	}
	nextgen.species<-species
	# evolve separately for burnin generations
	meanpref<-meantrait<-matrix(0,burnin+ngen,nsp)
	colnames(meanpref)<-colnames(meantrait)<-paste("sp",1:nsp)
	tt<-ts<-vector()
	for(i in 1:burnin){
		if(i%%pf==0) message(paste("burn-in generation",i))
		for(j in 1:nsp){
			for(k in 1:nindivs[j]){
				mated<-0
				# pick mother & father
				while(mated==0){
					# female at random
					female<-rInt(nindivs[j])
					while(species[[j]]$individual[[female]]$sex!=0) female<-rInt(nindivs[j])
					# a male is nearby
					male<-rInt(nindivs[j])
					while(species[[j]]$individual[[male]]$sex!=1) male<-rInt(nindivs[j])
					# does she seem him?
					pseen<-exp(-species[[j]]$individual[[male]]$trait^2/w.trait)
					# does she like him?
					plike<-exp(-(species[[j]]$individual[[male]]$trait-species[[j]]$individual[[female]]$pref)^2/w.pref)
					# does she accept him?
					if((pseen*plike)>runif(n=1)) mated=1
				}
				# now mate them
				chrom<-matrix(r01(n=nloci.pref*2)+1,nloci.pref,2)
				for(l in 1:nloci.pref){
					mut=0.0
					if(u>runif(n=1)) mut=rnorm(n=1,sd=sqrt(mutvar))
					nextgen.species[[j]]$individual[[k]]$genome.pref[l,1]<-species[[j]]$individual[[female]]$genome.pref[l,chrom[l,1]]+mut
				}
				for(l in 1:nloci.pref){
					mut=0.0
					if(u>runif(n=1)) mut=rnorm(n=1,sd=sqrt(mutvar))
					nextgen.species[[j]]$individual[[k]]$genome.pref[l,2]<-species[[j]]$individual[[male]]$genome.pref[l,chrom[l,2]]+mut
				}
				nextgen.species[[j]]$individual[[k]]$pref<-sum(nextgen.species[[j]]$individual[[k]]$genome.pref)
				meanpref[i,j]<-meanpref[i,j]+nextgen.species[[j]]$individual[[k]]$pref/nindivs[j]
				chrom<-matrix(r01(n=nloci.trait*2)+1,nloci.trait,2)
				for(l in 1:nloci.trait){
					mut=0.0
					if(u>runif(n=1)) mut=rnorm(n=1,sd=sqrt(mutvar))
					nextgen.species[[j]]$individual[[k]]$genome.trait[l,1]<-species[[j]]$individual[[female]]$genome.trait[l,chrom[l,1]]+mut
				}
				for(l in 1:nloci.trait){
					mut=0.0
					if(u>runif(n=1)) mut=rnorm(n=1,sd=sqrt(mutvar))
					nextgen.species[[j]]$individual[[k]]$genome.trait[l,2]<-species[[j]]$individual[[male]]$genome.trait[l,chrom[l,2]]+mut
				}
				nextgen.species[[j]]$individual[[k]]$trait<-sum(nextgen.species[[j]]$individual[[k]]$genome.trait)
				meantrait[i,j]<-meantrait[i,j]+nextgen.species[[j]]$individual[[k]]$trait/nindivs[j]
				# assign sex randomly
				nextgen.species[[j]]$individual[[k]]$sex<-r01()
			}
		}
		species<-nextgen.species
		tt[i]<-i
	}
	# burn-in over
	# now start evolving the species together
	mismatings<-totalmatings<-matrix(rep(0,ngen*nsp),ngen,nsp)
	colnames(mismatings)<-colnames(totalmatings)<-paste("sp",1:nsp)
	for(i in 1:ngen){
		if(i%%pf==0) message(paste("simulation generation",i))
		for(j in 1:nsp){
			for(k in 1:nindivs[j]){
				mated<-0
				keepmom<-0
				# pick mother & father
				while(mated==0){
					if(keepmom==0){
						# take a female from species j at random
						female<-rInt(nindivs[j])
						while(species[[j]]$individual[[female]]$sex!=0)
							female<-rInt(nindivs[j])
					}
					# a male is nearby
					male<-rInt(sum(nindivs))
					# get male species
					malesp<-getSp(male,nindivs)
					if(malesp>1) male<-male-sum(nindivs[1:(malesp-1)])
					while(species[[malesp]]$individual[[male]]$sex!=1){
						male<-rInt(sum(nindivs))
						malesp<-getSp(male,nindivs)
						if(malesp>1) male<-male-sum(nindivs[1:(malesp-1)])
					}
					# does she seem him?
					pseen<-exp(-species[[malesp]]$individual[[male]]$trait^2/w.trait)
					# does she like him?
					plike<-exp(-(species[[malesp]]$individual[[male]]$trait-species[[j]]$individual[[female]]$pref)^2/w.pref)
					# does she accept him?
					if(pseen>runif(n=1)){
						if(plike>runif(n=1)){
							if(j==malesp){
								mated<-1
								keepmom<-0
							} else {
								mismatings[i,j]<-mismatings[i,j]+1
								keepmom<-0
							}
							totalmatings[i,j]<-totalmatings[i,j]+1
						} else {
							if(j!=malesp) keepmom<-1
							else keepmom<-0
						}
					}
				}
				# now mate them
				chrom<-matrix(r01(n=nloci.pref*2)+1,nloci.pref,2)
				for(l in 1:nloci.pref){
					mut=0.0
					if(u>runif(n=1)) mut=rnorm(n=1,sd=sqrt(mutvar))
					nextgen.species[[j]]$individual[[k]]$genome.pref[l,1]<-species[[j]]$individual[[female]]$genome.pref[l,chrom[l,1]]+mut
				}
				for(l in 1:nloci.pref){
					mut=0.0
					if(u>runif(n=1)) mut=rnorm(n=1,sd=sqrt(mutvar))
					nextgen.species[[j]]$individual[[k]]$genome.pref[l,2]<-species[[j]]$individual[[male]]$genome.pref[l,chrom[l,2]]+mut
				}
				nextgen.species[[j]]$individual[[k]]$pref<-sum(nextgen.species[[j]]$individual[[k]]$genome.pref)
				meanpref[i+burnin,j]<-meanpref[i+burnin,j]+nextgen.species[[j]]$individual[[k]]$pref/nindivs[j]
				chrom<-matrix(r01(n=nloci.trait*2)+1,nloci.trait,2)
				for(l in 1:nloci.trait){
					mut=0.0
					if(u>runif(n=1)) mut=rnorm(n=1,sd=sqrt(mutvar))
					nextgen.species[[j]]$individual[[k]]$genome.trait[l,1]<-species[[j]]$individual[[female]]$genome.trait[l,chrom[l,1]]+mut
				}
				for(l in 1:nloci.trait){
					mut=0.0
					if(u>runif(n=1)) mut=rnorm(n=1,sd=sqrt(mutvar))
					nextgen.species[[j]]$individual[[k]]$genome.trait[l,2]<-species[[j]]$individual[[male]]$genome.trait[l,chrom[l,2]]+mut
				}
				nextgen.species[[j]]$individual[[k]]$trait<-sum(nextgen.species[[j]]$individual[[k]]$genome.trait)
				meantrait[i+burnin,j]<-meantrait[i+burnin,j]+nextgen.species[[j]]$individual[[k]]$trait/nindivs[j]
				# assign sex randomly
				nextgen.species[[j]]$individual[[k]]$sex<-r01()
			}
		}
		species<-nextgen.species
		tt[i+burnin]<-i+burnin
		ts[i]<-i
	}
	# ok simulation over
	ninety<-sqrt(-w.trait*log(0.90))
	tt<-c(0,tt)
	meantrait<-rbind(rep(0,nsp),meantrait)
	meanpref<-rbind(rep(0,nsp),meanpref)
	rownames(meantrait)<-rownames(meanpref)<-c(paste("b",0:burnin,sep=""),paste("s",1:ngen,sep=""))
	# create plots
	if(figs=="on"||figs=="minimal"){
		if(nsp==1){
			yMax<-max(1.1*c(meantrait,meanpref,ninety))
			yMin<-min(1.1*c(meantrait,meanpref,-ninety))
			if(figs=="on") x11()
			plot(tt,meantrait[,1],"l",ylim=c(yMin,yMax),col="blue",xlab="generation",ylab="signal trait",lwd=2)
			lines(tt,meanpref[,1],col="red",lwd=2)
			lines(c(tt[1],tt[length(tt)]),c(0,0))
			lines(c(tt[1],tt[length(tt)]),c(ninety,ninety),lty=2)
			lines(c(tt[1],tt[length(tt)]),-c(ninety,ninety),lty=2)
		} else {
			if(nsp==2){
				yMax<-max(1.1*c(meantrait,meanpref,ninety))
				yMin<-min(1.1*c(meantrait,meanpref,-ninety))
				if(figs=="on") x11()
				plot(tt,meantrait[,1],"l",ylim=c(yMin,yMax),col="blue",xlab="generation",ylab="signal trait",lwd=2)
				lines(tt,meanpref[,1],col="red",lwd=2)
				lines(tt,meantrait[,2],col="green",lwd=2)
				lines(tt,meanpref[,2],col="yellow",lwd=2)			
				lines(c(tt[1],tt[length(tt)]),c(0,0))
				lines(c(tt[1],tt[length(tt)]),c(ninety,ninety),lty=2)
				lines(c(tt[1],tt[length(tt)]),-c(ninety,ninety),lty=2)
				if(figs=="on"){ 
					x11(); plot(ts,(mismatings[,1]/totalmatings[,1])/(1-nindivs[1]/sum(nindivs)),type="l",xlab="generation",ylab="relative mismating",col="blue",lwd=2,xlim=c(0,max(ts)),ylim=c(0,1.2))
					lines(ts,(mismatings[,2]/totalmatings[,2])/(1-nindivs[2]/sum(nindivs)),type="l",col="green",lwd=2)
				}
			} else {
				yMax<-max(1.1*c(meantrait,meanpref,ninety))
				yMin<-min(1.1*c(meantrait,meanpref,-ninety))
				if(figs=="on") x11()
				plot(tt,meantrait[,1],"l",ylim=c(yMin,yMax),col="blue",xlab="generation",ylab="signal trait",lwd=2)
				lines(tt,meanpref[,1],col="red",lwd=2)
				lines(tt,meantrait[,2],col="green",lwd=2)
				lines(tt,meanpref[,2],col="yellow",lwd=2)
				lines(tt,meantrait[,3],col="black",lwd=2)
				lines(tt,meanpref[,3],col="purple",lwd=2)		
				lines(c(tt[1],tt[length(tt)]),c(0,0))
				lines(c(tt[1],tt[length(tt)]),c(ninety,ninety),lty=2)
				lines(c(tt[1],tt[length(tt)]),-c(ninety,ninety),lty=2)
				if(figs=="on"){ 
					x11(); plot(ts,(mismatings[,1]/totalmatings[,1])/(1-nindivs[1]/sum(nindivs)),type="l",ylab="relative mismating",xlab="generation",col="blue",lwd=2,xlim=c(0,max(ts)),ylim=c(0,1.2))
					lines(ts,(mismatings[,2]/totalmatings[,2])/(1-nindivs[2]/sum(nindivs)),type="l",col="green",lwd=2)
					lines(ts,(mismatings[,3]/totalmatings[,3])/(1-nindivs[3]/sum(nindivs)),type="l",col="black",lwd=2)
				}
			}
		}
	}
	Sys.sleep(1)
	return(list(trait=meantrait,pref=meanpref))
}

# function random 0 or 1
r01<-function(n=1) rbinom(n,size=1,prob=0.5)

# function to pull random integer
rInt<-function(max) round(runif(n=1)*max+0.5)

# return the species of an individual
getSp<-function(id,nindivs){
	sp<-1
	nsp<-length(nindivs)
	for(i in 2:nsp) if((id>sum(nindivs[1:i-1]))&&(id<=sum(nindivs[1:i]))) sp<-i
	return(sp)
}
