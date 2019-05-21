!*******************************************************************!
!   This subroutine calculates E(Yij_til|Zj=1,Yj=yj), for i=1,..cj	!	
!   under Dorfman retesting											!
!*******************************************************************!
      subroutine eyijtilz1(eyj,pj,cj,yj,ssp,row,col,SeI,Sp,Sej)
      implicit none
      integer cj,row,col,sj
      integer yj(cj),ssp(row,col)
      double precision eyj(cj),pj(cj)
      double precision SeI,Sp,Sej(cj)

      integer i,r,j,yjc(cj-1)
      double precision sm,pjc(cj-1)
      double precision rsl
	  
      call pz1yj(sm,yj,pj,Sp,cj,ssp,row,col,SeI,Sej)
	  
      do 105 i = 1,cj
         call icomp_real(pjc,pj,cj,i)	
         call icomp_int(yjc,yj,cj,i)
         call pyzjyij(rsl,pj,pjc,cj,i,yj,yjc,ssp,row,col,SeI,Sp,Sej)
         eyj(i) = rsl/sm
105   continue
      return
      end

!********************************************************************!
!   This subroutine calculates P(Zj = 1, Yj = yj, Ytil_ij = 1)
!
      subroutine pyzjyij(res,pj,pjc,cj,sj,yj,yjc,ssp,row,col,SeI,Sp,Sej)
      implicit none	  
      integer cj,row,col,sj
      integer yj(cj),yjc(cj-1),ssp(row,col)
      double precision res,pj(cj),pjc(cj-1)
      double precision SeI,Sp	  
      double precision Sej(cj)
	  
      integer sm,i,s,rsum(row)
      integer s11,s10,s01,s00
      double precision prd,h,temp1,temp2
	  
      do 10 i = 1, row	  
	      sm = 0	 
          do 15 s = 1, col
              sm = sm + ssp(i,s)	  
15        continue
          rsum(i) = sm
10    continue		  

      res = 0
      do 20 i = 1, row
	      s11 = 0
          s10 = 0
          s01 = 0
          s00 = 0
          prd = 1.0		  
          h = Sej(rsum(i)+1)
          do 25 s = 1, col
             s11 = s11 + yjc(s)*ssp(i,s)
			 s10 = s10 + yjc(s)*(1-ssp(i,s))
			 s01 = s01 + (1-yjc(s))*ssp(i,s)
			 s00 = s00 + (1-yjc(s))*(1-ssp(i,s))
			 prd = prd*(pjc(s)**ssp(i,s))*((1-pjc(s))**(1-ssp(i,s)))
25        continue
          temp1 = (SeI**s11)*((1-SeI)**s01)*((1-Sp)**s10)*(Sp**s00)
          temp2 = (SeI**yj(sj))*((1-SeI)**(1-yj(sj)))
          res = res + h*temp1*temp2*prd*pj(sj)		  
20    continue
      return
      end

!*************************************************************!
!   This subroutine calculates: Pr(Zj = 1, Yj = yj)
!   This is the probability that a pool tests positively
!   under Dorfman decoding algorithm
!*************************************************************!
      subroutine pz1yj(sm,Yj,pj,Sp,cj,ssp,row,col,SeI,Sej)
      implicit none	  
      integer cj, row, col
      integer Yj(cj), ssp(row,col)
      double precision sm, pj(cj)
      double precision Sp, SeI	  
      double precision Sej(cj)
	  
      integer s, k
      double precision prob, psyz(cj+1)
      double precision Sek
	    
      sm = 0.0
      do 10 s = 0, cj
          call psyzj(prob,s,Yj,pj,Sej,Sp,SeI,cj,ssp,row,col)
          psyz(s+1) = prob
          sm = sm + prob
10    continue
      
      return
      end    
	  
!**************************************************************!	  
! 	Calculates: pr(Y1j=y1,...,Ycj,Zj=1,\tilY1j+...+\tilYcjj = k)
!
      subroutine psyzj(res,k,Yj,pj,Sej,Sp,SeI,cj,ssp,row,col)
      implicit none
      integer k,cj,row,col
      integer Yj(cj),ssp(row,col)	  
      double precision pj(cj),Sej(cj)
      double precision res,Sp,SeI	  
      
      integer i,s,t
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
                  tst = Sej(SY)
              end if		  
           
              lk = prod1*tst*prod2		   
          else
		      lk = 0.0
          end if
		  
          res = res + lk
20    continue
	  
      return
      end      

!***************************************************	  
      subroutine icomp_int(newvec,vec,vlen,i)
      implicit none	  
      integer vlen,i 	  
      integer newvec(vlen),vec(vlen)
 
      integer j,r 
      j = 1	  
      do 110 r = 1,vlen
         if(r .ne. i) then		 
	         newvec(j) = vec(r)
             j = j + 1
         end if			 
110   continue

      return
      end
	  
!***************************************************	  
      subroutine icomp_real(newvec,vec,vlen,i)
      implicit none	  
      integer vlen,i 	  
      double precision newvec(vlen),vec(vlen)
 
      integer j,r 
      j = 1	  
      do 110 r = 1,vlen
         if(r .ne. i) then		 
	         newvec(j) = vec(r)
             j = j + 1
         end if			 
110   continue

      return
      end
	  