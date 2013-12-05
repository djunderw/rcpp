#cov:mat, rr:vector
haz<-function(ti, covi, rr, aa, bb, t1) {
	n_ti<-max(which(t1<ti)) 
	ifelse(t1[n_ti]==0, exp(covi%*%rr), exp(covi%*%rr + sum(aa*exp(-bb*(ti-t1[2:n_ti])))))
	}
	
haz2<-function(ti, bb, t1) {
	n_ti<-max(which(t1<ti)) 
	ifelse(t1[n_ti]==0, 0, sum(exp(-bb*(ti-t1[2:n_ti]))))
	}

haz3<-function(ti, aa, bb, t1) {
	n_ti<-max(which(t1<ti)) 
	ifelse(t1[n_ti]==0, 0, sum(aa*exp(-bb*(ti-t1[2:n_ti]))*(t1[2:n_ti]-ti)))
	}
	
haz7<-function(ti, aa, bb, t1) {
	n_ti<-max(which(t1<ti)) 
	ifelse(t1[n_ti]==0, 0, sum(aa*exp(-bb*(ti-t1[2:n_ti]))*((t1[2:n_ti]-ti)^2)))
	}
	
haz8<-function(ti, bb, t1) {
	n_ti<-max(which(t1<ti)) 
	ifelse(t1[n_ti]==0, 0, sum(exp(-bb*(ti-t1[2:n_ti]))*(t1[2:n_ti]-ti)))
	}

Mnum_i_derv<-function(i, j, length_ni, dat, x, zz){
	cov_num<-dim(zz)[2]
	rr=x[1:cov_num]; aa=x[(cov_num+1)]; bb=x[(cov_num+2)]

	ni<-length_ni[i]
	if(ni>2){ 
	t_i<-dat[i,j]
	num_ij<-haz(ti=t_i, covi= zz[i,], rr=rr, aa=aa, bb=bb, t1=dat[i,])
	num_ij1<-haz2(ti=t_i, bb=bb, t1=dat[i,])
	num_ij2<-haz3(ti=t_i, aa=aa, bb=bb, t1=dat[i,])
	num_ij3<-haz7(ti=t_i, aa=aa, bb=bb, t1=dat[i,])
	num_ij4<-haz8(ti=t_i, bb=bb, t1=dat[i,])
	}#if
	res_num_ij<-c(num_ij,num_ij1,num_ij2,num_ij3,num_ij4)
	return(res_num_ij)
}

Mnum_derv<-function(i, length_ni, dat, x, zz){	
	nc<-dim(dat)[2]
	temp_num_i1<-rep(NA, (nc-1))
	temp_num_i<-c()
	ni<-length_ni[i]
	if(ni>2){
	#temp_num_i[2:(ni-1)]<-sapply(2:(ni-1), Mnum_i_derv, i=i, length_ni=length_ni, dat=dat, x=x, zz=zz)
		for (k in 2:(ni-1)){
			temp_num_i1<-sapply(k, Mnum_i_derv, i=i, length_ni=length_ni, dat=dat, x=x, zz=zz)
			temp_num_i<-cbind(temp_num_i, temp_num_i1)
		}
	}#if
	return(temp_num_i)
}
		
Mden_i_derv<-function(i, j, length_ni, dat, x, zz){
	cov_num<-dim(zz)[2]
	rr=x[1:cov_num]; aa=x[(cov_num+1)]; bb=x[(cov_num+2)]
	comb<-combn(seq(1:cov_num),2)
	comb_num<-dim(comb)[2]
	
	comb_zz<-c()
	for (jj in 1:comb_num){
		pre_comb_zz<-zz[,comb[1,jj]]*zz[,comb[2,jj]]
		comb_zz<-cbind(comb_zz, pre_comb_zz)
	}
	
	nr<-dim(dat)[1]
	nc<-dim(dat)[2]
	ni<-length_ni[i]
	
	if(ni>2){
	t_i<-dat[i,j]
		
	temp_Aij<-c()
	temp_Aij1<-c()	
	temp_Aij2<-c()
	temp_Aij3<-c()
	temp_Aij4<-c()
	temp_Aij5<-c()
	temp_Aij6<-c()
	temp_Aij7<-c()	
	temp_Aij8<-c()
	temp_Aij9<-c()
	temp_Aij10<-c()
	temp_Aij11<-c()	
	temp_Aij12<-c()
	
	res_Aij<-rep(NA,(8+4*cov_num+comb_num))
	
	for(k in 1:nr){
		#print(k)
		if (t_i<=dat[k,nc]) { #note: tie data <=
			Aijk<-haz(covi= zz[k,], ti=t_i, rr=rr, aa=aa, bb=bb, t1=dat[k,])
			Aijk2<-haz2(ti=t_i, bb=bb, t1=dat[k,])
			Aijk3<-haz3(ti=t_i, aa=aa, bb=bb, t1=dat[k,])
			Aijk7<-haz7(ti=t_i, aa=aa, bb=bb, t1=dat[k,])
			Aijk8<-haz8(ti=t_i, bb=bb, t1=dat[k,])
			temp_Aij<-c(temp_Aij, Aijk)
			temp_Aij1<-cbind(temp_Aij1, Aijk*zz[k,])
			temp_Aij2<-c(temp_Aij2, Aijk*Aijk2)
			temp_Aij3<-c(temp_Aij3, Aijk*Aijk3)
			temp_Aij4<-cbind(temp_Aij4, Aijk*(zz[k,]^2))
			temp_Aij5<-c(temp_Aij5, Aijk*(Aijk2^2))
			temp_Aij6<-c(temp_Aij6, Aijk*(Aijk3^2))
			temp_Aij7<-c(temp_Aij7, Aijk*Aijk7)
			temp_Aij8<-cbind(temp_Aij8, Aijk*zz[k,]*Aijk2)
			temp_Aij9<-cbind(temp_Aij9, Aijk*zz[k,]*Aijk3)
			temp_Aij10<-c(temp_Aij10, Aijk*Aijk2*Aijk3)
			temp_Aij11<-c(temp_Aij11, Aijk*Aijk8)
			temp_Aij12<-cbind(temp_Aij12, Aijk*comb_zz[k,])
			}#if
		}#k
	res_Aij[1]<-sum(temp_Aij,na.rm=TRUE)
	res_Aij[2:(1+cov_num)]<-rowSums(temp_Aij1,na.rm=TRUE)
	res_Aij[(2+cov_num)]<-sum(temp_Aij2,na.rm=TRUE)
	res_Aij[(3+cov_num)]<-sum(temp_Aij3,na.rm=TRUE)
	res_Aij[(4+cov_num):(3+2*cov_num)]<-rowSums(temp_Aij4,na.rm=TRUE)
	res_Aij[(4+2*cov_num)]<-sum(temp_Aij5,na.rm=TRUE)
	res_Aij[(5+2*cov_num)]<-sum(temp_Aij6,na.rm=TRUE)
	res_Aij[(6+2*cov_num)]<-sum(temp_Aij7,na.rm=TRUE)
	res_Aij[(7+2*cov_num):(6+3*cov_num)]<-rowSums(temp_Aij8,na.rm=TRUE)
	res_Aij[(7+3*cov_num):(6+4*cov_num)]<-rowSums(temp_Aij9,na.rm=TRUE)
	res_Aij[(7+4*cov_num)]<-sum(temp_Aij10,na.rm=TRUE)
	res_Aij[(8+4*cov_num)]<-sum(temp_Aij11,na.rm=TRUE)
	res_Aij[(9+4*cov_num):(8+4*cov_num+comb_num)]<-rowSums(temp_Aij12,na.rm=TRUE)
	}
	return(res_Aij)
}

Mden_derv<-function(i, length_ni, dat, x, zz){
	nc<-dim(dat)[2]
	#temp_Ai<-rep(NA, (nc-1))
	temp_Ai1<-rep(NA, (nc-1))
	temp_Ai<-c()
	ni<-length_ni[i]
	if(ni>2){
		#temp_Ai[2:(ni-1)]<-sapply(2:(ni-1), Mden_i, i=i, length_ni=length_ni, dat=dat, x=x, zz=zz)
		for (k in 2:(ni-1)){
		temp_Ai1<-sapply(k, Mden_i_derv, i=i, length_ni=length_ni, dat=dat, x=x, zz=zz)
		temp_Ai<-cbind(temp_Ai, temp_Ai1)
		}
	}#if
	return(temp_Ai)
}
