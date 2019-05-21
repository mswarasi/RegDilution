
!*******************************************************************!
!   This subroutine calculates E(Yij_til|Zj=0), for i=1,..cj		!
!   under Dorfman retesting											!
!*******************************************************************!
      subroutine eyijtilz0(eyj,pj,cj,Sp,Sej)
      implicit none
      integer cj
      double precision eyj(cj),pj(cj),Sp
      double precision Sej(cj)
	  
      integer i,k
      double precision psub(cj-1),prob,pz0 
      double precision PHjki(cj,cj),PHjk(cj+1)
      double precision sm
	  
      do 25 i = 1,cj
         call icomp(psub,pj,i,cj)
         do 30 k = 0,(cj-1)
            call dbinf(prob,k,psub,cj-1)
            PHjki(i,k+1) = prob	
30       continue
25    continue
                              
      do 35 k = 0,cj
         call dbinf100(prob,k,pj,cj)
            PHjk(k+1) = prob
35    continue	  

      pz0 = 0.0
      do 40 i = 1, (cj+1)
         if(i .eq. 1) then
            pz0 = pz0 + Sp*PHjk(i)	  
         else
            pz0 = pz0 + (1.0 - Sej(i-1))*PHjk(i)	  
         end if
40    continue

      do 45 i = 1, cj
         sm = 0.0
         do 50 k = 1, cj
            sm = sm + PHjki(i,k)*(1.0 - Sej(k))	  
50       continue
         eyj(i) = sm*pj(i)/pz0
45    continue
		 
      return
      end

!--------------- not iid binomial probability ---------------------!  
      subroutine dbinf100(res,n,p,cj)
      implicit none 
      integer n,cj
      double precision res,p(cj)  
	  
      double precision ppos(cj),pfrac(cj),pcomp(cj),temp(cj-n+1)
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
	   
!--------------- not iid binomial probability ---------------------!
      subroutine dbinf(res,n,p,cj)
      implicit none
      integer n, cj
      double precision res, p(cj)  
	  
      double precision ppos(cj), pfrac(cj), pcomp(cj), temp(cj-n+1)
      double precision summ, prod, sm
      integer i, j, k, ind
	  
      do 110 i=1, cj	  
         ppos(i) = p(i)
	     pfrac(i) = p(i)/(1.0-p(i))
	     pcomp(i) = 1.0 - p(i)
110   continue

      if(n .eq. 0) then
         res = prod(pcomp, cj)
      end if
 
      if(n .gt. 0) then
          if(n .eq. 1) then
             res = summ(pfrac, cj)*prod(pcomp, cj)
          end if
       		 
          if(n .gt. 1) then
           	  sm = 0.0
              do 120 j=cj,n,-1
                 sm = sm + pfrac(j)
                 temp(j-(n-1)) = sm
120           continue	 

              if(n .gt. 2) then
                 do 130 i=1,(n-2) 
                    ind = 1
                    do 140 j = (n-i),(cj-i)
                       temp(ind) = pfrac(j)*temp(ind)
				       ind = ind+1
140                 continue
                    sm = 0.0
                    do 150 k=(cj-n+1),1,-1
				       sm = sm + temp(k)
				       temp(k) = sm
150                 continue  
130              continue	
              end if
              sm = 0.0
              do 160 i=1,(cj-n+1)
			     sm = sm + pfrac(i)*temp(i)
160           continue
              res = sm*prod(pcomp, cj)
           end if
       end if
       return
       end		      	   
	   
!-------- calculates p(-i) ----------------
      subroutine icomp(res,p,i,cj)
      implicit none  
      integer i,cj
      double precision res(cj-1),p(cj)
      integer ind, j
	  
      ind = 1
      do 165 j = 1, cj
         if(j .ne. i) then        
            res(ind) = p(j)
            ind = ind + 1
         end if
165   continue
      return 
      end           
	  
!------ sum function ---------------------	  	  
      function summ(vec, vlen)
      implicit none
      integer i, vlen
      double precision summ, vec(vlen)
      summ = 0.0
      do 170 i = 1, vlen
         summ = summ + vec(i)
170   continue
      return
      end	      	
	  
!------ product function ------------------  	  
      function prod(vec, vlen)
      implicit none
      integer i, vlen
      double precision prod, vec(vlen)
      prod = 1.0
      do 180 i = 1, vlen
         prod = prod*vec(i)
180   continue
      return
      end	     	

