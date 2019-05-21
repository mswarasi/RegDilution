!*****************************************************************!
!   This routine calculates the "correlation" part of the
!   variance-covariance matrix using Louis's method. 
!	Note: This one works with only logit link 'g' for the 
!   individuals' regression model, but with any submodel 'h' 
! 	which should must be specified through R. 
! 	Date: 12/17/2016
!
      subroutine cov_dl_glgt(cov,p,Z,X,N,J,GI,c,row,col,Ymat,plen,Se_mat,cmax,dSe_mat)
      implicit none
      integer N, J, GI, row, col, plen, cmax
      integer c(J),Z(J)
      double precision cov(row,col), p(N), X(N,plen-1)
      double precision Se_mat(cmax,cmax),dSe_mat(cmax,cmax)
      integer Ymat(N,GI)

      double precision dldb(plen,GI), dlbdlb(plen*plen,GI)
      double precision dldb_sm(plen), dlbdlb_sm(plen*plen)  
	  
      integer g,i,SY,jj,k
      integer Ijk,ind(J+1)
      double precision sm
      integer h,idk,temp0
      double precision temp	  
      double precision term1
	  
      ind(1) = 0
      temp0 = 0
      do 72 i=1,J
         temp0 = temp0 + c(i)
         ind(i+1) = temp0		 
72    continue
	  
      do 10 g = 1, GI	

	  
      do 5 h = 1, (plen-1)	
         temp = 0
		 do 11 i = 1, N
            temp = temp + (Ymat(i,g)*(1-p(i))-(1-Ymat(i,g))*p(i))*X(i,h)	  
11       continue
         dldb(h,g) = temp
5     continue
		 	  
      temp = 0
      do 12 jj = 1, J

         SY = 0		 
         do 13 i = (ind(jj)+1), ind(jj+1)
            SY = SY + Ymat(i,g)	  
13       continue
	  
         sm = 0
         do 14 k = 1, c(jj)
            Ijk = 0
            if(SY .eq. k) then
               Ijk = 1           
            end if
			if(Z(jj) .gt. 0) then
			    term1 = dSe_mat(c(jj),k)/Se_mat(c(jj),k)
			else
			    term1 = -dSe_mat(c(jj),k)/(dble(1.0)-Se_mat(c(jj),k))
			end if
			sm = sm + term1*Ijk
14       continue
		 temp = temp + sm 
12    continue

 	  dldb(plen,g) = temp

      idk = 1
      do 15 i = 1, plen
         do 20 jj = 1, plen
            dlbdlb(idk,g) = dldb(i,g)*dldb(jj,g)
			idk = idk + 1
20       continue    		 
15    continue	  

10    continue

      do 25 i = 1, plen
         temp = 0
		 do 30 g = 1, GI
		    temp = temp + dldb(i,g)
30       continue	
		 dldb_sm(i) = temp
25    continue	  
	  
	  do 35 i = 1, (plen*plen)
	     temp = 0
		 do 40 g = 1, GI
		    temp = temp + dlbdlb(i,g)
40       continue		 
         dlbdlb_sm(i) = temp
35    continue	  
	  
      idk = 1
	  do 45 i = 1, plen
         do 50 jj=1, plen
		    cov(i,jj) = dlbdlb_sm(idk)/dble(GI) - (dldb_sm(i)/dble(GI))*(dldb_sm(jj)/dble(GI))
			idk = idk + 1
50       continue		 
45    continue	  
	  
      return
      end
  