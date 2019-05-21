!***************************************************************************!
!			  calculates E(I( sum(\tilde Yij = k) )| Zj=0) 					!
!***************************************************************************!
      subroutine eszj(res,pj,Sp,cj,Sej)
      implicit none
      integer cj
      double precision res(cj),pj(cj),Sp
      double precision Sej(cj)
	  
      double precision prb(cj+1)  
      integer k,s
      double precision prob, pz0
  
      do 60 k = 1, cj	  
         do 65 s = 0, cj      
            call dbinf100(prob,s,pj,cj,cj) 
            prb(s+1) = prob   	  
65       continue

         pz0 = 0.0
         do 70 s = 1, (cj+1) 	  	 
            if(s .eq. 1) then
               pz0 = pz0 + prb(s)*Sp
            else
               pz0 = pz0 + prb(s)*(1.0 - Sej(s-1))
            end if
70       continue

         res(k) = (1.0 - Sej(k))*prb(k+1)/pz0
60    continue	 
 
      return
      end   
	   
!*****************************************************!
! Binomial probability with non-identical
! success probability
!   
      subroutine dbinf100(res,n,p,NN,cj)
      implicit none  
      integer n, NN, cj
      double precision res, p(NN)  
	  
      double precision ppos(cj), pfrac(cj), pcomp(cj), temp(cj-n+1)
      double precision summ, prod, sm
      integer i, j, k, ind
	  
      do 80 i=1, cj	  
         ppos(i) = p(i)
         pfrac(i) = p(i)/(1.0-p(i))
         pcomp(i) = 1.0 - p(i)
80    continue

      if(n .eq. 0) then
         res = prod(pcomp, cj)
      end if
 
      if(n .gt. 0) then
          if(n .eq. 1) then
             res = summ(pfrac, cj)*prod(pcomp, cj)
          end if
       		 
          if(n .gt. 1) then
           	  sm = 0.0
              do 85 j=cj,n,-1
                 sm = sm + pfrac(j)
                 temp(j-(n-1)) = sm
85            continue	 

              if(n .gt. 2) then
                 do 90 i=1,(n-2) 
                    ind = 1
                    do 95 j = (n-i),(cj-i)
                       temp(ind) = pfrac(j)*temp(ind)
				       ind = ind+1
95                  continue
                    sm = 0.0
                    do 100 k=(cj-n+1),1,-1
				       sm = sm + temp(k)
				       temp(k) = sm
100                 continue  
90               continue	
              end if
              sm = 0.0
              do 105 i=1,(cj-n+1)
			     sm = sm + pfrac(i)*temp(i)
105           continue
              res = sm*prod(pcomp, cj)
           end if
       end if
       return
       end	  
	   
!**************** sum function ***************************!	  	  
      function summ(vec, N)
      implicit none
      integer N
      double precision vec(N)
      double precision summ
      integer i
	  
      summ = 0.0
      do 920 i = 1, N
         summ = summ + vec(i)
920    continue
      return
      end	     
	  
!************* product function **************************!  	  
      function prod(vec, N)
      implicit none
      integer i, N
      double precision vec(N)
      double precision prod
      prod = 1.0
      do 930 i = 1, N
         prod = prod*vec(i)
930    continue
      return
      end	       	
 