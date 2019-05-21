!*******************************************************************************!
! 	This subroutine calculates log likelihood function for			            !
! 	dilution model under master pool testing 						            !
! 	Last updated: 11/25/2016								            		!
!*******************************************************************************!
      subroutine Dil_MPT_LogLik(ll,pt,gamma,Z,Sp,cvec,J,N,Se_mat,cmax)
      implicit none
      integer J,N,Se_row,cmax
      double precision ll,pt(N),gamma,Sp
      double precision Se_mat(cmax,cmax)
      integer Z(J),cvec(J)
	  
      double precision ptemp(cmax),PHjk(cmax+1)
      double precision prob,pz1,temp1,temp2
      integer ind(J+1),i,iter,s,cj,temp0

      ind(1) = 0
      temp0 = 0
      do 72 i=1,J
         temp0 = temp0 + cvec(i)
         ind(i+1) = temp0		 
72    continue  
	  
      ll = 0
      do 3 i=1,J
        cj = cvec(i)
        iter = 1
        do 5 s=(ind(i)+1),ind(i+1)	
           ptemp(iter) = pt(s)
           iter = iter + 1
5       continue	

        do 10 s = 0, cj      
           call dbin(prob,s,ptemp,cj,cmax) 
	       PHjk(s+1) = prob   	  
10    	continue

      pz1 = 0.0	
      do 20 s = 1, (cj+1) 
	     if(s .eq. 1) then
		     temp1 = PHjk(s)*(1.0-Sp)
         else
		     temp1 = PHjk(s)*Se_mat(cj,s-1)
         end if
         pz1 = pz1 + temp1
20    continue
	  
      if(Z(i) .gt. 0) then
          temp2 = log(pz1)    
      else
          temp2 = log(1.0 - pz1)	  
      end if  
      ll = ll + temp2
3     continue
      return
      end

!*******************************************************************!
!	This subroutine calculates probabilities for binomial 
!	distribution with non-identical success probability
!
      subroutine dbin(res,n,p,cj,cmax)
      implicit none
      integer n,cj,cmax
      double precision res, p(cmax)  
	  
      double precision ppos(cj), pfrac(cj), pcomp(cj), temp(cj-n+1)
      double precision summ, prod, sm
      integer i, j, k, ind
	  
      do 610 i=1, cj	  
	     ppos(i) = p(i)
	     pfrac(i) = p(i)/(1.0-p(i))
	     pcomp(i) = 1.0 - p(i)
610   continue

      if(n .eq. 0) then
        res = prod(pcomp, cj)
      end if
 
      if(n .gt. 0) then
          if(n .eq. 1) then
             res = summ(pfrac, cj)*prod(pcomp, cj)
          end if
       		 
          if(n .gt. 1) then
           	  sm = 0.0
              do 620 j=cj,n,-1
                sm = sm + pfrac(j)
                temp(j-(n-1)) = sm
620           continue	 

              if(n .gt. 2) then
                do 630 i=1,(n-2) 
                  ind = 1
                  do 640 j = (n-i),(cj-i)
                    temp(ind) = pfrac(j)*temp(ind)
				    ind = ind+1
640               continue
                  sm = 0.0
                  do 650 k=(cj-n+1),1,-1
				    sm = sm + temp(k)
				    temp(k) = sm
650               continue  
630             continue	
              end if
              sm = 0.0
              do 660 i=1,(cj-n+1)
			    sm = sm + pfrac(i)*temp(i)
660           continue
            res = sm*prod(pcomp, cj)
       end if
       end if
       return
       end

!*******************************************************************!	   
!	sum function 
!	  	  
      double precision function summ(vec, vlen)
      implicit none
      integer i, vlen
      double precision vec(vlen)
      summ = 0.0
      do 700 i = 1, vlen
         summ = summ + vec(i)
700   continue
      return
      end	       

!*******************************************************************!	  
!	product function   
!	  
      double precision function prod(vec, vlen)
      implicit none
      integer i, vlen
      double precision vec(vlen)
      prod = 1.0
      do 800 i = 1, vlen
         prod = prod*vec(i)
800   continue
      return
      end	      
