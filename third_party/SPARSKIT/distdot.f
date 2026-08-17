       function distdot(n,x,ix,y,iy)
       integer*8 n, ix, iy
       real*8 distdot, x(*), y(*), ddot
       external ddot
       distdot = ddot(n,x,ix,y,iy)
       return
       end