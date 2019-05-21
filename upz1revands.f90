!***********************************************************************!
! 		Dorfman decoding: Vandsteelandt et al. (2000)					!
! 		that is, with constant Se/Sp 									!		
!   	This function calculates: 									    !
!   	Pr(Zj = 1) = P(Zj,Y1j = y1,...,Ycjj = ycj)						!
! 		Last modified: 08/05/2015										!
!***********************************************************************!
      subroutine upz1revands(sm,Yj,pj,Sp,cj,ssp,row,col,SeI)
      implicit none
      integer cj, row, col
      integer Yj(cj), ssp(row,col)
      double precision sm, pj(cj)
      double precision Sp
	  
      double precision Sej(cj), SeI
      integer s, k
      double precision prob, psyz(cj+1)
      double precision Sek
	    
      sm = 0.0
      do 10 s = 0, cj
          call psyzj(prob,s,Yj,pj,Sej,Sp,SeI,cj,cj,ssp,row,col)
          psyz(s+1) = prob
          sm = sm + prob
10    continue
      
      return
      end
!********************************************************************	  
! 	Calculates: pr(Y1j=y1,...,Ycj,  Zj=1,  \tilY1j+...+\tilYcjj = k)
!
      subroutine psyzj(res,k,Yj,pj,Sej,Sp,SeI,cj,cmax,ssp,row,col)
      implicit none
      integer k, cj, cmax, row, col
      integer Yj(cmax), ssp(row,col)	  
      double precision pj(cmax), Sej(cmax)
      double precision res, Sp, SeI	  
      
      integer i, s, t
      integer sm
      integer rsum(row) 
      double precision  prd1, prod1, prd2, prod2	  
      integer SY, ytil, y
      double precision tst, lk
      
      do 10 i = 1, row	  
	      sm = 0	 
          do 15 s = 1, col
              sm = sm + ssp(i,s)	  
15        continue
          rsum(i) = sm
10    continue
      
      res = 0.0
      do 20 i = 1, row	 
          if((rsum(i) .gt. (k-2)) .and. (rsum(i) .lt. (k+1))) then
              prod1 = 1.0	
	          SY = 0
	          prod2 = 1.0
              do 25 t = 1, (col+1)
	              if(t .eq. 1) then
                      ytil = k - rsum(i)
                  else
                      ytil = ssp(i,(t-1))
                  end if
                  y = Yj(t)				
				  
                  if((ytil .eq. 1) .and. (y .eq. 1)) then  ! for term 1
		              prd1 = SeI
                  else if((ytil .eq. 1) .and. (y .eq. 0)) then 	
					  prd1 = 1.0 - SeI
                  else if((ytil .eq. 0) .and. (y .eq. 1)) then 	
					  prd1 = 1.0 - Sp	
                  else
					  prd1 = Sp			
                  end if
                  prod1 = prod1*prd1
				  
                  SY = SY + ytil  			! for term 2
			  
		          if(ytil .gt. 0) then 		! for term 3
		              prd2 = pj(t)
                  else
                      prd2 = 1.0 - pj(t)
                  end if
                  prod2 = prod2*prd2
25            continue			  
		  
              if(SY .eq. 0) then
                  tst = 1.0 - Sp
              else
                  tst = SeI
              end if		  
           
              lk = prod1*tst*prod2		   
          else
		      lk = 0.0
          end if
		  
          res = res + lk
20    continue
	  
      return
      end	  
 	  