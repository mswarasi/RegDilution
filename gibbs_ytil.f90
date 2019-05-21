!*********************************************************************
!     Gibbs sampler: generates Ytil_{ij}|{Y,Z} for all i, j  
!     Ytil is returened via res, an N by GI matrix
!
      subroutine gibbs_ytil(res,dataf,p,SeI,Sp,N,cmax,GI,burn,seed,Se_mat)
      IMPLICIT NONE
      INTEGER N, cmax, GI, burn, seed
      INTEGER res(N,GI),dataf(N,(4+cmax))
      double precision p(N),Sp,SeI,Se_mat(cmax,cmax)
	  
      INTEGER cj, SY, ind, g, i, s
      double precision zeta0, zeta1
      INTEGER temp
      double precision U 	

      call srand(seed)	  
      DO 100 g = 1,(burn+GI)
         DO 200 i=1,N
            dataf(i,1) = 0
            cj = dataf(i,4)
		  
            SY = 0
	        DO 300 s =1,cj
		       ind = dataf(i,4+s)
		       SY = SY + dataf(ind,1)
300         CONTINUE	 

!   success probability
            IF(dataf(i,3) .gt. 0) THEN
               zeta1 = p(i)*Se_mat(cj,SY+1)
	           zeta1 = zeta1*(SeI**dataf(i,2))*((1-SeI)**(1-dataf(i,2)))
            ELSE
		       zeta1 = p(i)*(1-Se_mat(cj,SY+1))
            END IF
	  
!   failure probability
           IF(dataf(i,3) .gt. 0) THEN
              IF(SY .eq. 0) THEN
		         zeta0 = (1-p(i))*(1-Sp)*((1-Sp)**dataf(i,2))*(Sp**(1-dataf(i,2)))
              ELSE
		         zeta0 = (1-p(i))*Se_mat(cj,SY)
                 zeta0 = zeta0*((1-Sp)**dataf(i,2))*(Sp**(1-dataf(i,2)))
              END IF
           ELSE
	          IF(SY .eq. 0) THEN
		         zeta0 = (1-p(i))*Sp
              ELSE
		         zeta0 = (1-p(i))*(1-Se_mat(cj,SY))
              END IF
           END IF
		   
!     now generate the Bernoulli random variable		
          zeta0 = zeta0/(zeta0+zeta1)
		  
          U = rand(0)
          IF(U .gt. zeta0) THEN
             temp = 1
             dataf(i,1) = 1
          ELSE
             temp = 0
             dataf(i,1) = 0
          END IF
	  
          IF(g .gt. burn) THEN
             res(i,(g-burn)) = temp
          END IF
	  
200     CONTINUE	 
100   CONTINUE	
      RETURN
      END
