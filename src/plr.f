      subroutine plr ( n, xn, zsmall, lenz, inform )

      implicit           double precision (a-h,o-z)
      integer            n, lenz, inform
      double precision   xn(n+1), zsmall(lenz)
*     ------------------------------------------------------------------
      integer            nwcore
      parameter          (nwcore = 10000000)
*     ------------------------------------------------------------------
      integer            m, nb, ne, nname,
     &                   nncon, nnobj, nnjac, iobj,  
     &                   mincor, ns, ninf, 
     &                   ka(n+1), name1, name2,
     &                   iprint, isumm, ispecs, i
      integer*4          ha(n), hs(n+1)
      double precision   objadd, sinf, obj,
     &                   a(n), bl(n+1), bu(n+1),
     &                   pi(1), rc(n+1), z(nwcore)
      character*8        names(5)
*     zero, one, infinity-----------------------------------------------
      double precision   zero,             one
      parameter         (zero   = 0.0d+0,  one    = 1.0d+0)
      double precision   bplus,            bminus
      parameter         (bplus  = 1.0d+20, bminus = -bplus)
*     ------------------------------------------------------------------
      iprint = 0   ! The MINOS PRINT   file.
      isumm  = 0   ! The MINOS SUMMARY file.
      ispecs = 0   ! The MINOS SPECS   file.
      call mistart( iprint, isumm, ispecs )  ! Initialize MINOS and open
*     ------------------------------------------------------------------
*     User workspace: 1  +  1  +  1  + 7 + (p+1)*nobs + nobs
*                   (nobs) (p)  (lam) (b)  (x, y data)  (w)  
*     ------------------------------------------------------------------      
      call miopti( 'Workspace (user) ', lenz, 0, 0, inform )
      call miopti( 'LOG FREQUENCY ', 0, 0, 0, inform )
      call miopti( 'PRINT LEVEL ', 0, 0, 0, inform )
      call miopti( 'SUMMARY FILE ', 0, 0, 0, inform )
      call miopti( 'SUMMARY FREQUENCY ', 0, 0, 0, inform )
      call miopti( 'SUPERBASICS LIMIT ', (n+10), 0, 0, inform )
*     ------------------------------------------------------------------
*     Now set parameters for moniss
*     ------------------------------------------------------------------
      do i = 1, lenz
         z(i) = zsmall(i)
      end do
      m = 1
      nb = n+1
      ne = n
      nname = 1
      nncon = 0
      nnobj = n
      nnjac = 0
      iobj = 0
      objadd = zero
      do i = 1, n
         a(i) = one
         ha(i) = 1
         ka(i) = i
      end do
      ka(n+1) = n+1
      do i = 1, nb
         bl(i) = bminus
         bu(i) = bplus
         hs(i) = 0
      end do
      pi(1) = zero
      call minoss( 'Cold', m, n, nb, ne, nname,
     $             nncon, nnobj, nnjac,
     $             iobj, objadd, names,
     $             a, ha, ka, bl, bu, name1, name2,
     $             hs, xn, pi, rc, 
     $             inform, mincor, ns, ninf, sinf, obj,
     $             z, nwcore )
      end ! plr solution

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine funobj( mode, n, x, f, g, nstate, nprob, z, nwcore )

      implicit           double precision (a-h,o-z)
      integer            mode, n, nstate, nprob, nwcore
      double precision   x(n), g(n), f, z(nwcore) 
      integer            nobs
*     ------------------------------------------------------------------
*     User workspace: 1  +  1  +  1  + 7 + (p+1)*nobs + nobs
*                   (nobs) (p)  (lam) (b)  (x, y data)  (w)  
*     ------------------------------------------------------------------      
      mode = mode
      nstate = nstate
      nprob = nprob
      nobs = int(z(1))
      call subfunobj( n, x, f, g, z, nwcore, nobs )
      return
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine subfunobj( n, x, f, g, z, nwcore, nobs )

      integer            n, nwcore, nobs
      double precision   x(n), f, g(n), z(nwcore)
      integer            i, j, ii, p
      double precision   eta(nobs), mu(nobs), wt(nobs), 
     &                   y(nobs), resid(nobs), xi(nobs), 
     &                   lam, bnorm, loglik, ddot
      double precision   zero,          one,          two
      parameter         (zero = 0.0d+0, one = 1.0d+0, two = 2.0d+0)
*     ------------------------------------------------------------------
      p = int(z(2))
      lam = z(3)
      ii = 10 
      do 200 i = 1, nobs
         eta(i) = zero
         do 100 j = 1, p
            eta(i) = eta(i) + x(j)*z(ii+(j-1)*nobs+i)
 100     continue
         y(i) = z(ii+p*nobs+i)
 200  continue
      ii = 10+(p+1)*nobs
      loglik = zero
      do 300 i = 1, nobs
         wt(i) = z(ii+i)
         mu(i) = one/(one+exp(-eta(i)))
         loglik = loglik + wt(i)*(y(i)*eta(i)-log(one+exp(eta(i))))
         resid(i) = wt(i)*(y(i)-mu(i))
 300  continue 
      bnorm = zero
      do 350 j = 2, p
         bnorm = bnorm + x(j)**2
 350  continue
      f = -loglik + lam*bnorm
      do i = 1, nobs
        xi(i) = -z(10+i)
      end do
      g(1) = ddot(nobs,xi,1,resid,1)
      do 500 j = 2, p
         ii = 10+(j-1)*nobs
         do 400 i = 1, nobs
            xi(i) = -z(ii+i)
 400     continue
         g(j) = ddot(nobs,xi,1,resid,1) + two*lam*x(j)
 500  continue
      return
      end

