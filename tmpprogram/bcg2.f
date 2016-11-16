      subroutine bicgstab2 (okprint,l, n, x, rhs, matvec, nonzero, tol,
     &      typestop,mxmv, work, ldw, info)
C 
C Simple BiCGstab(\ell) iterative method, \ell <= 2
C By M.A.Botchev, Jan.'98 
C report bugs to botchev@cwi.nl or botchev@excite.com
C
C Copyright (c) 1998 by M.A.Botchev
C Permission to copy all or part of this work is granted,
C provided that the copies are not made or distributed
C for resale, and that the copyright notice and this
C notice are retained.
C
C This is the "vanilla" version of BiCGstab(\ell) as described
C in PhD thesis of D.R.Fokkema, Chapter 3.  It includes two enhancements 
C to BiCGstab(\ell) proposed by G.Sleijpen and H.van der Vorst in
C 1) G.Sleijpen and H.van der Vorst "Maintaining convergence 
C    properties of BiCGstab methods in finite precision arithmetic",
C    Numerical Algorithms, 10, 1995, pp.203-223
C 2) G.Sleijpen and H.van der Vorst "Reliable updated residuals in
C    hybrid BiCG methods", Computing, 56, 1996, pp.141-163
C
C {{ This code based on:
C subroutine bistbl v1.0 1995
C
C Copyright (c) 1995 by D.R. Fokkema.
C Permission to copy all or part of this work is granted,
C provided that the copies are not made or distributed
C for resale, and that the copyright notice and this
C notice are retained.  }}
C
C okprint == (input) LOGICAL. If okprint=.true. residual norm
C            will be printed to *
C l       == (input) INTEGER BiCGstab's dimension <= 2
C            Set l=2 for highly nonsymmetric problems
C n       == (input) INTEGER size of the system to solve 
C x       == (input/output) DOUBLE PRECISION array dimension n
C            initial guess on input, solution on output
C rhs     == (input) DOUBLE PRECISION array dimension n
C            right-hand side (rhs) vector
C matvec  == (input) EXTERNAL name of matrix vector subroutine
C            to deliver y:=A*x by CALL matvec(n,x,y)
C nonzero == (input) LOGICAL tells
C            BiCGstab(\ell) if initial guess in x is zero or not. 
C            If nonzero is .FALSE., one MATVEC call is saved.
C tol     == (input/output) DOUBLE PRECISION tolerance for all possible
C            stopping criteria (see the 'typestop' parameter)
C            On output, if info=0 or 1, tol is actually achieved
C            residual reduction or whatever (see the 'typestop' parameter)
C typestop== (input) CHARACTER*3 stopping criterion (||.|| denotes 
C            the 2-norm):
C            typestop='rel' -- relative stopping crit.: ||res|| < tol*||res0||
C            typestop='abs' -- absolute stopping crit.: ||res||<tol
C            typestop='max' -- maximum  stopping crit.: max(abs(res))<tol
C NOTE(for typestop='rel' and 'abs'): To save computational work, the value of 
C            residual norm used to check the convergence inside the main iterative 
C            loop is computed from 
C            projections, i.e. it can be smaller than the true residual norm
C            (it may happen when e.g. the 'matrix-free' approach is used).
C            Thus, it is possible that the true residual does NOT satisfy
C            the stopping criterion ('rel' or 'abs').
C            The true residual norm (or residual reduction) is reported on 
C            output in parameter TOL -- this can be changed to save 1 MATVEC
C            (see comments at the end of the subroutine)
C mxmv   ==  (input/output) INTEGER.  On input: maximum number of matrix 
C            vector multiplications allowed to be done.  On output: 
C            actual number of matrix vector multiplications done
C work   ==  (workspace) DOUBLE PRECISION array dimension (n,2*l+5))
C ldw    ==  (input) INTEGER size of work, i.e. ldw >= n*(2*l+5)
C info   ==  (output) INTEGER.  info = 0 in case of normal computations
C            and 
C            info = -m (<0) - means paramater number m has an illegal value
C            info = 1 - means no convergence achieved (stopping criterion
C            is not fulfilled)
C            info = 2 - means breakdown of the algorithm (try to enlarge
C            parameter l=\ell to get rid of this)
C ----------------------------------------------------------
      implicit none
      external matvec
      integer l, n, mxmv, ldw, info
      double precision  x(n), rhs(n), tol
      logical   okprint,nonzero
      character*3   typestop
      double precision  work(n,5+2*l)
C     -----------------------------------------
      integer   lmax
      parameter(lmax=2)
      double precision rwork(lmax+1,3+2*(lmax+1))

      logical GoOn, rcmp, xpdt
      integer ii, i1, jj, kk, nmv
      double precision alpha,beta,hatgamma,kappa0, kappal,maxval1,
     &     mxnrmr,mxnrmx,omega,rho0,rho1,rnrm0,rnrm,rnrmMax,
     &     sigma,sum1,varrho
      integer z, zz, y0, yl, y
      integer rr, r, u, xp, bp

      double precision    zero, one, delta
      parameter(zero=0d0,one=1d0,delta=1d-2)

      info = 0
      
      if (l.gt.lmax .or. l.lt.1) info = -2
      if (tol.le.zero) info = -8
      if (mxmv.lt.0) info = -10
      
      rr = 1
      r = rr+1
      u = r+(l+1)
      xp = u+(l+1)
      bp = xp+1
      if (bp*n.gt.ldw) info = -12
      
      z = 1
      zz = z+(l+1)
      y0 = zz+(l+1)
      yl = y0+1
      y = yl+1
      
      if (info.ne.0) return

C      Initialize first residual
      if (nonzero) then
         call matvec (n, x, work(1,r) )
         do ii=1,n
            work(ii,r) = rhs(ii) - work(ii,r)
         enddo
         nmv = 1
      else
         do ii=1,n
            work(ii,r) = rhs(ii)
         enddo
         nmv = 0
      endif

C     Initialize iteration loop
      sum1 = zero
      do ii=1,n
         work(ii,rr) = work(ii,r)
         work(ii,bp) = work(ii,r)
         work(ii,xp) = x(ii)
         x(ii) = zero
         sum1=sum1+ work(ii,r)**2 
      enddo

      rnrm0 = sqrt( sum1 )
      rnrm = rnrm0
      if (typestop.eq.'max') then
         maxval1=zero
         do ii=1,n
            maxval1=max( maxval1, abs( work(ii,r) ) )
         enddo
         rnrmMax = maxval1
      endif

      mxnrmx = rnrm0
      mxnrmr = rnrm0
      rcmp = .false.
      xpdt = .false.

      alpha = zero
      omega = one
      sigma = one
      rho0 =  one

C     Iterate
      if(typestop.eq.'rel')then
         GoOn = rnrm.ge.tol*rnrm0 .and. nmv.lt.mxmv
      else if(typestop.eq.'abs')then
         GoOn = rnrm.ge.tol       .and. nmv.lt.mxmv
      else if(typestop.eq.'max')then
         GoOn = rnrmMax.ge.tol    .and. nmv.lt.mxmv
      else
         info = -9
         return
      end if

      do while (GoOn)
C     =====================
C     --- The BiCG part ---
C     =====================
      	 rho0 = -omega*rho0
      	 do kk=1,l
            sum1 = zero
            do ii=1,n
               sum1=sum1+ work(ii,rr)*work(ii,r+kk-1) 
            enddo
            rho1 = sum1

      	    if (rho0.eq.zero) then
      	       info = 2
      	       return
      	    endif
      	    beta = alpha*(rho1/rho0)
      	    rho0 = rho1
      	    do jj=0,kk-1
               do ii=1,n
                  work(ii,u+jj) = work(ii,r+jj) - beta*work(ii,u+jj)
               enddo
      	    enddo

            call matvec(n, work(1,u+kk-1), work(1,u+kk))
      	    nmv = nmv+1

            sum1 = zero
            do ii=1,n
               sum1=sum1+ work(ii,rr)*work(ii,u+kk) 
            enddo
            sigma = sum1
            
      	    if (sigma.eq.zero) then
      	       info = 2
      	       return
      	    endif

      	    alpha = rho1/sigma
            do ii=1,n
               x(ii) = alpha*work(ii,u) + x(ii)
            enddo
            
      	    do jj=0,kk-1
               do ii=1,n
                  work(ii,r+jj) = -alpha*work(ii,u+jj+1) + work(ii,r+jj)
               enddo
      	    enddo

	    call matvec (n, work(1,r+kk-1), work(1,r+kk))
      	    nmv = nmv+1

            sum1 = zero
            do ii=1,n
               sum1=sum1+ work(ii,r)**2 
            enddo
            rnrm = sqrt( sum1 )
            if (typestop.eq.'max') then
               maxval1=zero
               do ii=1,n
                  maxval1=max( maxval1, abs( work(ii,r) ) )
               enddo
               rnrmMax = maxval1
            endif
      	    mxnrmx = max (mxnrmx, rnrm)
      	    mxnrmr = max (mxnrmr, rnrm)
      	 enddo

C     ==================================
C     --- The convex polynomial part ---
C     ================================== 

C        Z = R'R
      	 do i1=1,l+1 
            do jj=i1-1,l
               sum1 = zero
               do ii=1,n
                  sum1=sum1+ work(ii,r+jj)*work(ii,r+i1-1) 
               enddo
               rwork(jj+1,z+i1-1) = sum1
               rwork(z+i1-1,jj+1) = rwork(jj+1,z+i1-1) 
            enddo
      	 enddo 

         do i1=zz,zz+l
            do ii=1,l+1
               rwork(ii,i1)   = rwork(ii,i1+(z-zz)) 
            enddo
         enddo
C        tilde r0 and tilde rl (small vectors)

       	 rwork(1,y0) = -one
         rwork(2,y0) = rwork(2,z) / rwork(2,zz+1)
      	 rwork(l+1,y0) = zero

      	 rwork(1,yl) = zero
         rwork(2,yl) = rwork(2,z+l) / rwork(2,zz+1)
      	 rwork(l+1,yl) = -one

C        Convex combination
         do ii=1,l+1
            rwork(ii,y) = zero
         enddo
         do jj=1,l+1
            do ii=1,l+1
               rwork(ii,y) = rwork(ii,y) + rwork(jj,yl)*
     &              rwork(ii,z+jj-1)
            enddo
         enddo
         sum1 = zero
         do ii=1,l+1
            sum1=sum1+ rwork(ii,yl)*rwork(ii,y) 
         enddo
         kappal = sqrt( sum1 )

         do ii=1,l+1
            rwork(ii,y) = zero
         enddo
         do jj=1,l+1
            do ii=1,l+1
               rwork(ii,y) = rwork(ii,y) + rwork(jj,y0)*
     &              rwork(ii,z+jj-1)
            enddo
         enddo
         sum1 = zero
         do ii=1,l+1
            sum1=sum1+ rwork(ii,y0)*rwork(ii,y)  
         enddo
         kappa0 = sqrt( sum1 )

         sum1 = zero
         do ii=1,l+1
            sum1=sum1+ rwork(ii,yl)*rwork(ii,y) 
         enddo
         varrho = sum1  
         varrho = varrho / (kappa0*kappal)
         
      	 hatgamma = sign(1d0,varrho)*max(abs(varrho),7d-1) *
     &        (kappa0/kappal)

         do ii=1,l+1
            rwork(ii,y0) = -hatgamma*rwork(ii,yl) +      rwork(ii,y0)
         enddo

C        Update
      	 omega = rwork(l+1,y0)
         do jj=1,l
            do ii=1,n
               work(ii,u) = work(ii,u) - rwork(1+jj,y0)*work(ii,u+jj)
               x(ii)      = x(ii)      + rwork(1+jj,y0)*work(ii,r+jj-1)
               work(ii,r) = work(ii,r) - rwork(1+jj,y0)*work(ii,r+jj)
            enddo
         enddo

         do ii=1,l+1
            rwork(ii,y) = zero
         enddo
         do jj=1,l+1
            do ii=1,l+1
               rwork(ii,y) = rwork(ii,y) + rwork(jj,y0)*
     &              rwork(ii,z+jj-1)
            enddo
         enddo

         sum1 = zero
         do ii=1,l+1
            sum1=sum1+ rwork(ii,y0)*rwork(ii,y) 
         enddo
         rnrm = sqrt( sum1 )

C     ================================
C     --- The reliable update part ---
C     ================================
      	 mxnrmx = max (mxnrmx, rnrm)
      	 mxnrmr = max (mxnrmr, rnrm)
      	 xpdt = (rnrm.lt.delta*rnrm0.and.rnrm0.lt.mxnrmx)
      	 rcmp = ((rnrm.lt.delta*mxnrmr.and.rnrm0.lt.mxnrmr) .or.xpdt)
      	 if (rcmp) then
      	    call matvec (n, x, work(1,r) )
            nmv = nmv + 1
            do ii=1,n
               work(ii,r) = work(ii,bp) - work(ii,r)
            enddo
      	    mxnrmr = rnrm
      	    if (xpdt) then
               do ii=1,n
                  work(ii,xp) = x(ii) + work(ii,xp)
                  x(ii) = zero
                  work(ii,bp) = work(ii,r)
               enddo
      	       mxnrmx = rnrm
      	    endif
      	 endif

         if(typestop.eq.'rel')then
            GoOn = rnrm.ge.tol*rnrm0 .and. nmv.lt.mxmv
            if(okprint) print *, nmv,' ',rnrm/rnrm0
         else if(typestop.eq.'abs')then
            GoOn = rnrm.ge.tol       .and. nmv.lt.mxmv
            if(okprint) print *, nmv,' ',rnrm
         else if(typestop.eq.'max')then
            maxval1=zero
            do ii=1,n
               maxval1=max(maxval1, abs( work(ii,r) ) )
            enddo
            rnrmMax = maxval1
            GoOn = rnrmMax.ge.tol    .and. nmv.lt.mxmv
            if(okprint) print *, nmv,' ',rnrmMax
         end if

      enddo

C     =========================
C     --- End of iterations ---
C     =========================

      do ii=1,n
         x(ii) = work(ii,xp) + x(ii)
      enddo

C --------------------- One matvec can be saved by commenting out this:
C                       (this is to compute the true residual)
      call matvec (n, x, work(1,r) )
      do ii=1,n
         work(ii,r) = rhs(ii) - work(ii,r)   
      enddo

      sum1 = zero
      do ii=1,n
         sum1=sum1+ work(ii,r)**2 
      enddo
      rnrm = sqrt( sum1 )
C --------------------- One matvec can be saved by commenting out this^

      if(typestop.eq.'rel')then
         if (rnrm.gt.tol*rnrm0) info = 1
         tol = rnrm/rnrm0
      else if(typestop.eq.'abs')then
         if (rnrm.gt.tol) info = 1
         tol = rnrm
      else if(typestop.eq.'max')then
         if (rnrmMax.gt.tol) info = 1
         tol = rnrmMax
      end if

      mxmv = nmv

      return
      end
