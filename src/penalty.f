      subroutine penalty(n,m,q,x,y,bnd,tlist,tlptr,tlend,rax,jax,ned,
&     eps,ierr)
      integer n,m,q,lp,lpl,ned,ierr
      integer bnd(n),tlist(q),tlptr(q),tlend(n),n4(4),p4(4),jax(m)
      double precision x(n),y(n),rax(m),eps
      double precision x4(4),y4(4),g4(4)
      ned = 0
      do 23000 i=1,n
      lpl = tlend(i)
      lp = lpl
23002 continue
      lp = tlptr(lp)
      j = iabs(tlist(lp))
      if(.not.(j .gt. i))goto 23005
      n4(1) = i
      n4(2) = j
      call fadjs(n4,n,q,tlist,tlptr,tlend)
      if(.not.(bnd(i)*bnd(j) .eq. 0))goto 23007
      ned = ned + 1
      do 23009 k = 1,4
      x4(k) = x(n4(k))
      y4(k) = y(n4(k))
23009 continue
      call ggap(x4,y4,g4,eps,ierr)
      if(.not.(ierr .eq. 1))goto 23011
      return
23011 continue
      call srtpai(n4,1,p4,1,4)
      do 23013 k = 1,4
      rax((ned - 1)*4 + k) = g4(p4(k))
      jax((ned - 1)*4 + k) = n4(p4(k))
23013 continue
      if(.not.(ned*4 .gt. m))goto 23015
      return
23015 continue
23007 continue
23005 continue
      if(.not.(lp .eq. lpl))goto 23017
      goto 23004
23017 continue
23003 goto 23002
23004 continue
23000 continue
      return
      end
      subroutine fadjs(n4,n,q,tlist,tlptr,tlend)
      integer n,q,vp,vpl,v,v0,match
      integer n4(4),tlist(q),tlptr(q),tlend(n)
      match = 0
      vpl = tlend(n4(1))
      vp = vpl
      k = 0
23019 continue
      k = k+1
      vp = tlptr(vp)
      v = tlist(vp)
      if(.not.(k.gt.1 .and. iabs(v) .eq. n4(2)))goto 23022
      n4(3) = iabs(v0)
      match = 1
      goto 23020
23022 continue
      if(.not.(match .gt. 0))goto 23024
      n4(4) = iabs(v)
      goto 23021
23024 continue
      v0 = v
23020 goto 23019
23021 continue
      return
      end
      subroutine ggap(x,y,g,eps,ierr)
      double precision x(4),y(4),g(4),w(2,4),h(2),d1,d2,eps
      d1 = -x(2) * y(1) + x(3) * y(1) + x(1) * y(2) -x(3) * y(2) - x(1) 
&     * y(3) + x(2) * y(3)
      d2 = -x(2) * y(1) + x(4) * y(1) + x(1) * y(2) -x(4) * y(2) - x(1) 
&     * y(4) + x(2) * y(4)
      if(.not.(dabs(d1) .lt. eps .or. dabs(d2) .lt. eps))goto 23026
      ierr = 1
      return
23026 continue
      h(1) = -(y(1) - y(2))
      h(2) = (x(1) - x(2))
      w(1, 1) = (y(2) - y(3))/d1 - (y(2) - y(4))/d2
      w(2, 1) = (x(3) - x(2))/d1 - (x(4) - x(2))/d2
      w(1, 2) = (y(3) - y(1))/d1 - (y(4) - y(1))/d2
      w(2, 2) = (x(1) - x(3))/d1 - (x(1) - x(4))/d2
      w(1, 3) = (y(1) - y(2))/d1
      w(2, 3) = (x(2) - x(1))/d1
      w(1, 4) = (y(2) - y(1))/d2
      w(2, 4) = (x(1) - x(2))/d2
      do 23028 i = 1,4
      g(i) = h(1)*w(1,i)+h(2)*w(2,i)
23028 continue
      ierr = 0
      return
      end
