      subroutine solveplr(x, n, z, lenz, nobs, lam, status)

      implicit none

      external myvalue, mygrad

      integer n,lenz,nobs,status
      double precision x(n),z(lenz),lam

      integer iter,nfunc,ngrad
      double precision eps,gnorm,f,d(n),g(n),xtemp(n),gtemp(n)

      eps = 1.d-8

      call cg_descent(eps, x, n, myvalue, mygrad, status,
     &                gnorm, f, iter, nfunc, ngrad, d, g, xtemp, gtemp,
     &                z, lenz, nobs, lam)

      end

      subroutine myvalue(f, x, n, z, lenz, nobs, lam)

      implicit none

C     INPUT VARIABLES
      integer n,lenz,nobs
      double precision f,x(n),z(lenz),lam

C     LOCAL VARIABLES
      integer i,j,iy,iw,io,hasoffset
      double precision y(nobs),w(nobs),offset(nobs),eta(nobs),
     &                 loglik,norm

      hasoffset = int(z(lenz))
      if (hasoffset .gt. 0) then
         io = (n+2)*nobs
         do i = 1,nobs
            offset(i) = z(io+i)
         end do
      else
         do i = 1,nobs
            offset(i) = 0.0
         end do
      endif

      iy = n*nobs
      iw = (n+1)*nobs
      do i = 1,nobs
         y(i) = z(iy+i)
         w(i) = z(iw+i)
         eta(i) = offset(i)
         do j = 1,n
            eta(i) = eta(i) + x(j)*z((j-1)*nobs+i)
         end do
      end do

      loglik = 0.0
      do i = 1,nobs
         loglik = loglik + w(i)*(y(i)*eta(i)-log(1.0+exp(eta(i))))
      end do

      norm = 0.0
      if (n .gt. 1) then
         do j = 2,n
            norm = norm + x(j)**2
         end do
      end if

      f = -loglik + lam*norm

      end

      subroutine mygrad(g, x, n, z, lenz, nobs, lam)

      implicit none

C     INPUT VARIABLES
      integer n,lenz,nobs
      double precision g(n),x(n),z(lenz),lam

C     LOCAL VARIABLES
      integer i,j,iy,iw,ip,io,hasoffset
      double precision y(nobs),w(nobs),offset(nobs),
     &                 eta(nobs),mu(nobs),resid(nobs),xi(nobs),
     &                 ddot

      hasoffset = int(z(lenz))
      if (hasoffset .gt. 0) then
         io = (n+2)*nobs
         do i = 1,nobs
            offset(i) = z(io+i)
         end do
      else
         do i = 1,nobs
            offset(i) = 0.0
         end do
      endif

      iy = n*nobs
      iw = (n+1)*nobs
      do i = 1,nobs
         y(i) = z(iy+i)
         w(i) = z(iw+i)
         eta(i) = offset(i)
         do j = 1,n
            eta(i) = eta(i) + x(j)*z((j-1)*nobs+i)
         end do
      end do

      do i = 1,nobs
         mu(i) = 1.0/(1.0+exp(-eta(i)))
         resid(i) = w(i)*(y(i)-mu(i))
      end do

      do j = 1,n
         ip = (j-1)*nobs
         do i = 1,nobs
            xi(i) = -z(ip+i)
         end do
         g(j) = ddot (nobs, xi, 1, resid, 1)
      end do

      if (n. gt. 1) then
         do j = 2,n
            g(j) = g(j) + 2.0*lam*x(j)
         end do
      end if

      end

      double precision function ddot(n,dx,incx,dy,incy)

      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n

      ddot = 0.0d0
      dtemp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20

      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      ddot = dtemp
      return

   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +
     *   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue
   60 ddot = dtemp
      return
      end
