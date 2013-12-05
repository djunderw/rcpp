p_lik<-function(dat, x, zz){
	cov_num<-dim(zz)[2]
	rr=x[1:cov_num]	
	aa=x[(cov_num+1)]; bb=x[(cov_num+2)]
	
	zz<-as.matrix(zz)
	dat<-as.matrix(dat)
	nr<-dim(dat)[1]#subj
	nc<-dim(dat)[2]#events & cs
		
	length_ni<-apply(dat, 1, function(x) length(na.omit(x)))
	
	num<-lapply(1:nr, Mnum_derv, length_ni=length_ni, dat=dat, x=x, zz=zz)
	den<-lapply(1:nr, Mden_derv, length_ni=length_ni, dat=dat, x=x, zz=zz)
	
	ll_ij<-c()
	for (i in 1:nr){
		if(is.null(num[[i]])){ll_ij1<-NA}
		else{ll_ij1<-log(num[[i]][1,])-log(den[[i]][1,])}		
		ll_ij<-c(ll_ij, ll_ij1)
	}	
	l<-sum(ll_ij, na.rm=TRUE)
	return(-l)
}#fn

derv_fns<-function(dat, x, zz){	
	cov_num<-dim(zz)[2]
	rr=x[1:cov_num]	
	aa=x[(cov_num+1)]; bb=x[(cov_num+2)]
		
	dat<-as.matrix(dat)
	nr<-dim(dat)[1]#subj
	nc<-dim(dat)[2]#events & cs
		
	length_ni<-apply(dat, 1, function(x) length(na.omit(x)))
	
	num<-lapply(1:nr, Mnum_derv, length_ni=length_ni, dat=dat, x=x, zz=zz)
	den<-lapply(1:nr, Mden_derv, length_ni=length_ni, dat=dat, x=x, zz=zz)
	
	temp_sc_r<-c()
	temp_sc_a<-c()
	temp_sc_b<-c()
	
	temp_sd_rr<-c()
	temp_sd_aa<-c()
	temp_sd_bb<-c()
	temp_sd_ra<-c()
	temp_sd_rb<-c()
	temp_sd_ab<-c()
	
	temp_sd_rr2<-c()
	
	comb<-combn(seq(1:cov_num),2)
	comb_num<-dim(comb)[2]

	for (i in 1:nr){
		if(is.null(num[[i]])){sc_r1<-rep(NA,cov_num); sc_a1<-NA; sc_b1<-NA 
			sd_rr1<-rep(NA,cov_num); sd_aa1<-NA; sd_bb1<-NA; sd_ra1<-rep(NA,cov_num); sd_rb1<-rep(NA,cov_num); sd_ab1<-NA
			sd_rr2<-rep(NA,comb_num)}
		else{
			sc_r1<-zz[i,]-(den[[i]][2:(1+cov_num),]/den[[i]][1,])
			sc_a1<-num[[i]][2,]-(den[[i]][(2+cov_num),]/den[[i]][1,])
			sc_b1<-num[[i]][3,]-(den[[i]][(3+cov_num),]/den[[i]][1,]) #score funs
			
			sd_rr1<-((den[[i]][(4+cov_num):(3+2*cov_num),]*den[[i]][1,])-(den[[i]][2:(1+cov_num),])^2)/(den[[i]][1,]^2)
			sd_aa1<-((den[[i]][(4+2*cov_num),]*den[[i]][1,])-(den[[i]][(2+cov_num),])^2)/(den[[i]][1,]^2)
			sd_bb1<-(-num[[i]][4,])+(((den[[i]][(5+2*cov_num),]+den[[i]][(6+2*cov_num),])*den[[i]][1,])-((den[[i]][(3+cov_num),])^2))/(den[[i]][1,]^2)
			
			sd_ra1<-((den[[i]][(7+2*cov_num):(6+3*cov_num),]*den[[i]][1,])-(den[[i]][2:(1+cov_num),]*den[[i]][(2+cov_num),]))/(den[[i]][1,]^2)
			sd_rb1<-((den[[i]][(7+3*cov_num):(6+4*cov_num),]*den[[i]][1,])-(den[[i]][2:(1+cov_num),]*den[[i]][(3+cov_num),]))/(den[[i]][1,]^2)
			sd_ab1<-(-num[[i]][5,])+(((den[[i]][(7+4*cov_num),]+den[[i]][(8+4*cov_num),])*den[[i]][1,])-(den[[i]][(2+cov_num),]*den[[i]][(3+cov_num),]))/(den[[i]][1,]^2) 
			
			den_rr2<-c()
			for (jj in 1:comb_num){
			pre_den_rr2<-den[[i]][2:(1+cov_num),][comb[1,jj]]*den[[i]][2:(1+cov_num),][comb[2,jj]]
			den_rr2<-c(den_rr2, pre_den_rr2)
			}

			den[[i]][2:(1+cov_num),]
			sd_rr2<-((den[[i]][(9+4*cov_num):(8+4*cov_num+comb_num),]*den[[i]][1,])-(den_rr2))/(den[[i]][1,]^2)
			#-second derv funs
			}#else				
		temp_sc_r<-cbind(temp_sc_r, sc_r1)
		temp_sc_a<-c(temp_sc_a, sc_a1)
		temp_sc_b<-c(temp_sc_b, sc_b1)
		
		temp_sd_rr<-cbind(temp_sd_rr, sd_rr1)
		temp_sd_aa<-c(temp_sd_aa, sd_aa1)
		temp_sd_bb<-c(temp_sd_bb, sd_bb1)
		temp_sd_ra<-cbind(temp_sd_ra, sd_ra1)
		temp_sd_rb<-cbind(temp_sd_rb, sd_rb1)
		temp_sd_ab<-c(temp_sd_ab, sd_ab1)
		
		temp_sd_rr2<-cbind(temp_sd_rr2, sd_rr2)
		}#for i	
	sc_r<-rowSums(temp_sc_r, na.rm=TRUE)
	sc_a<-sum(temp_sc_a, na.rm=TRUE)
	sc_b<-sum(temp_sc_b, na.rm=TRUE)
	sc<-c(sc_r, sc_a, sc_b)
	
	sd_rr<-rowSums(temp_sd_rr, na.rm=TRUE)
	sd_aa<-sum(temp_sd_aa, na.rm=TRUE)
	sd_bb<-sum(temp_sd_bb, na.rm=TRUE)
	sd_ra<-rowSums(temp_sd_ra, na.rm=TRUE)
	sd_rb<-rowSums(temp_sd_rb, na.rm=TRUE)
	sd_ab<-sum(temp_sd_ab, na.rm=TRUE)
	
	sd_rr2<-rowSums(temp_sd_rr2, na.rm=TRUE)
	
	sd_mat1<-diag(sd_rr); sd_mat1[lower.tri(sd_mat1)]<-sd_rr2
	sd_mat1<-sd_mat1+t(sd_mat1); diag(sd_mat1)<-diag(sd_mat1)/2
	sd_mat2<-matrix(c(sd_aa, sd_ab, sd_ab, sd_bb),nrow=2, ncol=2, byrow=TRUE)
	sd_mat3<-matrix(c(sd_ra, sd_rb),nrow=cov_num, ncol=2)
	sd_mat4<-cbind(sd_mat1, sd_mat3)
	sd_mat5<-cbind(t(sd_mat3), sd_mat2)
	sd_mat<-rbind(sd_mat4, sd_mat5)
	
	res_derv<-list(sc,sd_mat)
	return(res_derv)
}#fn

