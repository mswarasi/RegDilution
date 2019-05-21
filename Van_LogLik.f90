!*******************************************************************!
! 	This subroutine calculates log likelihood function for			!
! 	Vandsteelandts et al. (2000) model, where Se/Sp is constant.	!
! 	Last modification date: 11/27/2016								!
!*******************************************************************!
      subroutine Van_LogLik(ll,pt,Z,Sp,Se,cvec,J,N) 
      implicit none
      integer J,N
      double precision ll,pt(N),Sp,Se
      integer cvec(J),Z(J)

      double precision pz1,prd,temp
      integer s,i,temp0,ind(J+1)

      ind(1) = 0
      temp0 = 0
      do 7 i=1,J
         temp0 = temp0 + cvec(i)
         ind(i+1) = temp0		 
7     continue
  
      ll = 0
      do 10 i=1,J
         prd = 1
         do 15 s=(ind(i)+1),ind(i+1)
            prd = prd*(1-pt(s))
15       continue
         pz1 = Se+(1-Se-Sp)*prd

         if(Z(i) .eq. 1) then
            temp = log(pz1)
         end if
         
         if(Z(i) .eq. 0) then         
            temp = log(1-pz1)	
         end if	
         
         ll = ll + temp	         
10    continue
	  
      return
      end
	  
	  