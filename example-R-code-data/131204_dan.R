ex1_dat<-read.table("ex1_dat.txt",header=TRUE)
ex1_zz<-read.table("ex1_zz.txt",header=TRUE)

ex1_zz<-as.matrix(ex1_zz)

source("sim2_num_den1.R")
source("sim2_pl_derv1.R")

derv_fns(dat=ex1_dat, zz=ex1_zz, x=c(-0.3,1,0.5,10))
#results
[[1]]
            X           zz2                             
  -0.53964222 -153.57897821   -2.66116346    0.07557093 

[[2]]
              [,1]      [,2]          [,3]         [,4]
[1,] 29.1745104316 -1.578662  0.0003082929  0.006879055
[2,] -1.5786617970 52.217884 -1.0134314844  0.132611004
[3,]  0.0003082929 -1.013431  9.9628832465 -0.373769324
[4,]  0.0068790545  0.132611 -0.3737693238  0.027106078


library(numDeriv)
hessian(p_lik, dat=ex1_dat, zz=ex1_zz, x=c(-0.3,1,0.5,10))
#results
             [,1]         [,2]       [,3]        [,4]
[1,]  31.22300744 -10.10535429 -0.8921878  0.02983206
[2,] -10.10535429 101.50350657 -0.6590392  0.04627048
[3,]  -0.89218782  -0.65903919  9.9628832 -0.37376932
[4,]   0.02983206   0.04627048 -0.3737693  0.02710608

#If you see these two matrices then except 2*2 matrix in the right below part the other values are not matched.
#I need to fix these errors. They are caused by derivative calculation, so if the 2*2 matrix values are correct then it will be fine! 