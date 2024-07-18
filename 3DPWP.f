      program Ideal
      
          
      parameter (nt=20000, nx=851, ny=1805, nz=92, nts=121, 
     x nn=19000, dx=12.e3, dt=60, nxr=83)

c     nt为总时间间隔，nx、ny、nz为x、y和z方向格点数
c     nts和nxr为输入台风数据时所用，nn为这么多步长输出一次结果，
c     dx为x和y方向间隔，dt为时间间隔，需符合cfl条件，dt<dx/sqrt(g*H)

      real rr ideal_t0 ideal_t1 ideal_t2 ideal_h0 ideal_h1 ideal_h2
      dimension zob(nts), sob(nts), tob(nts), uob(nts), vob(nts)
      dimension z(nz), dz(nz), z1(nz)
      dimension tinit(nz), sinit(nz)
      dimension alpha(nz), beta(nz)
      dimension x(nx),y(ny),yp(ny) 
      dimension u(nx,ny,nz),v(nx,ny,nz),t(nx,ny,nz),s(nx,ny,nz)
      dimension uo(nx,ny,nz),vo(nx,ny,nz),to(nx,ny,nz),so(nx,ny,nz)
      dimension tham(nx,ny,nz),tvam(nx,ny,nz),tentm(nx,ny,nz)      
      dimension dtham(nx,ny,nz),dtvam(nx,ny,nz),dtentm(nx,ny,nz)
      dimension timd(nt)
      dimension xi(nxr*nxr),yi(nxr*nxr),uspe(nxr*nxr)
     x,vspe(nxr*nxr),tauxi(nxr*nxr),tauyi(nxr*nxr)
      dimension taux(nx,ny),tauy(nx,ny)
C     dimension tautheta(nx,ny),taur(nx,ny)
C     dimension taumax(nx,ny)

      dimension tt(nx,ny,nz),ut(nx,ny,nz),vt(nx,ny,nz),st(nx,ny,nz)
      dimension t9(nx,ny,nz),u9(nx,ny,nz),v9(nx,ny,nz),s9(nx,ny,nz)
      dimension w(nx,ny,nz),pp(nx,ny,nz)
      dimension tb5(nz),sb5(nz),ub5(nz),vb5(nz),tb4(nz),sb4(nz),ub4(nz)
     x,vb4(nz),db4(nz)
      dimension ashe(nx,ny),dmlm(nx,ny),dtlm(nx,ny)
      dimension sx(nz), sy(nz)

          taumax= 2.293 ! 设置理想台风的最大风应力
          factor = 1 ! 此为用台风移速修正风速时的系数，2为hspd，1为0.5hspd
          windmax=(70*0.5144444-0.5*hspd)*0.92*0.8
          rmax= 130.  ! 设置理想台风的最大风速半径，单位km
          rhlat = 18.7           ! 设定纬度
          hspd = 8.5             ! 设定移动速度
          ntss = 2  ! leapfrog-trapezoidal two step iterative method，两次迭代
          b=0.5   ! 兰金涡旋的系数
          
c  these flags control the dynamics of the model:
      
c  realwind = 1 use wind from file stressfield.dat
c           = 2 use ideal wind of Rankine
C           = 3 use ideal wind of SLOSH
c           = 4 use ideal wind of Shay
c  zlevel= 0 use origional z level
c          = 1 use an unequal z level
C          = 2 use parabolic z level
C          = 3 use e exponential z level      
c  iupw  = 1 use upwind scheme in u,v calculation
c        = 0 use central scheme in u,v calculation
c  iupwz = 1 use upwind scheme in w calculation
c        = 0 use central scheme in w calculation
c  noasym = 1 shuts off the asymmetry of the storm winds
c         = 2 also shuts off radial wind component of storm winds
c  nonlcl = 1 shuts off horizontal advection, 
c         = 2 shuts off horizontal and vertical advection,
c         = 3 above plus applies spatially constant wind stress
c                (useful for testing 1-d dynamics)
c  nodeep = 1 shuts off ml deepening by entrainment
c  nopg = 1 shuts off pressure gradients
c  idia = 1 turns on diagnostics
c  ts_switch = 1 use Temperature and Salinity from txt document
c            = 2 use ideal Temperature and Salinity
c
c  set flags to default values
      
      pi = 3.1415926
      
      ts_switch = 1
      realwind = 3 
      zlevel = 1
      inorb = 0 ! 如果inorb为1，则垂向等间隔，间隔值为dzfix1      
      iupw = 0
      iupwz = 0     
      
      noasym = 0
      nonlcl = 0
      nodeep = 0
      nopg = 0
c assume that ep,q0,qsw=0
      ashe=0
      rsi=0
      emp=0
      
      idia = 0

      iglo = 0
      
      pscale = 1.0
      if(nopg.eq.1) pscale = 0.
c
      hascale = 1.0
      if(nonlcl.ge.1) hascale = 0.
c
      vascale = 1.0
      if(nonlcl.ge.2) then
      vascale = 0.
      hascale = 0.
      end if
      
      if (ts_switch.eq.1) then
          open(unit=9, file='TS_Haiou_new.txt')
C         open(unit=9,file='TSinitprofiles.dat') 
      elseif (ts_switch.eq.2) then
          idael_t0=30.0     !   海表温度
          ideal_t1=15.0     !   温跃层底  
          ideal_t2=2.0      !   海底温度
          ideal_h0=30.0     !   混合层厚度
          ideal_h1=300.0    !   温跃层厚度
          ideal_h2=4000.0   !   海洋深度
      endif
      open(unit=18, file='G:\fortran\ourput\dml.dat',
     xstatus='unknown')
      open(unit=19, file='G:\fortran\ourput\dtl.dat',
     xstatus='unknown')
      open(unit=20, file='G:\fortran\ourput\tham.dat',
     xstatus='unknown')
      open(unit=21, file='G:\fortran\ourput\tvam.dat',
     xstatus='unknown')
      open(unit=22, file='G:\fortran\ourput\tentm.dat',
     xstatus='unknown')
      open(unit=23, file='G:\fortran\ourput\dtham.dat',
     xstatus='unknown')
      open(unit=24, file='G:\fortran\ourput\dtvam.dat',
     xstatus='unknown')
      open(unit=25, file='G:\fortran\ourput\dtentm.dat',
     xstatus='unknown')
      
      open(unit=61, file='G:\fortran\ourput\u.dat',
     x status='unknown')
      open(unit=62, file='G:\fortran\ourput\v.dat',
     x status='unknown')
      open(unit=63, file='G:\fortran\ourput\w.dat',
     x status='unknown')
      open(unit=64, file='G:\fortran\ourput\t.dat',
     x status='unknown')
      open(unit=65, file='G:\fortran\ourput\s.dat',
     x status='unknown')  

      open(unit=71, file='G:\fortran\ourput\u.surface'
     x, status='unknown')
      open(unit=72, file='G:\fortran\ourput\v.surface'
     x, status='unknown')
      open(unit=73, file='G:\fortran\ourput\w.surface'
     x, status='unknown')
      open(unit=74, file='G:\fortran\ourput\t.surface',
     x status='unknown')
      open(unit=75, file='G:\fortran\ourput\s.surface',
     x status='unknown')     
      open(unit=76, file='G:\fortran\ourput\taux.dat',
     x status='unknown')
      open(unit=77, file='G:\fortran\ourput\tauy.dat',
     x status='unknown')   

c  set grid interval, meters
      ncx = (nx-1)/2.0    
      ncy = (ny-1)/2.0
      
      dy = dx  
      
      dt7 = dt
      dtm = dt
      
      y0 = -ncy*dy-3*rmax*1000.    ! 设定起始点
C     y0 = 0                ! 设定起始点
      
      nstep = nt            !  run this many time steps 
      nh = 0   
                                           
      iarb = 0               ! iarb=1则用给定（观测）的台风路径
      
      nxm=nx-1
      nym=ny-1
      
c  you can halve (or more) the time step by setting ihalf = 2,3,4 etc 

      ihalf = 1   ! 将时间间隔分割，2则为1/2
      naw = nstep                 !  start this no of steps before center 
      tzero = -dt*float(naw)/8.64e4   ! 按天算的时间（before center）
      dt29 = dt/float(ihalf)  ! 分割后的时间      
            
c  set several constants.  cp is speed used in radiation
c  boundary condition. 

      cp = 1.5
      crbc = -cp/dx
      cadv = -1./(2.*dx)
      ddx = 1./(2.*dx)
      ddy = 1./(2.*dy)
      g = 9.80665
      
      beta1 = 0.6
      beta2 = 20.
      udrag = 9999.
      rg = 0.25
      rb = 0.65
      
      d2r = 3.1415/180. ! 角度变为弧度的参数
      
      f = 2.*7.292e-5*sin(rhlat*3.141/180.)                                      
      ti = 2.*3.141/f

      do i=1,nx
          do j=1,ny
              taux(i,j)=0;
              tauy(i,j)=0;
          enddo
      enddo
      do k=1,nts
        read(9,*,end=653) zob(k),tob(k),sob(k)
        uob(k) = 0
        vob(k) = 0
C      write(*,*) zob(k),tob(k),sob(k),uob(k),vob(k)        
      enddo
      close(9)      
  653 continue
           
      if (zlevel.eq.1) then
          dzfix1 = 10.  ! 较浅处垂向间隔
          dzfix2 = 20.  ! 较深处垂向间隔
          dzfix3 = 100.  ! 最深处垂向间隔
C         dzfix4 = 50.  ! 最深处垂向间隔
          
          nz1 = 10
          nz2 = 18
C         nz3 = 23

          do 542 j=1,nz
          dz(j) = dzfix1
          if(inorb.eq.1) go to 6606
          if(j.gt.nz1) dz(j) = dzfix2
          if(j.gt.nz2) dz(j) = dzfix3
C         if(j.gt.nz3) dz(j) = dzfix4
          go to 6607
 6606     continue      
 6607     continue          
          z(1)=dzfix1
          if(j.gt.1) z(j) = z(j-1) + dz(j)  
          z1(j) = z(j)  + dz(j)/2.   ! z1  is the mid depth of a grid cell
  542     continue  
      elseif (zlevel.eq.2) then
          do j=1,nz
              interval=ideal_h2/nz+1
              z(j) = (interval*j)**2/ideal_h2
          enddo
              dz(1)=z(2)/2
          do j=2,nz-1
              dz(j)=(z(j+1)-z(j-1))/2
          enddo
              dz(nz)=z(nz)-z(nz-1)
      elseif (zlevel.eq.3) then
          interval=ideal_h2/nz
          z(j) = (interval*j)**3.71828/ideal_h2**2.71828        
          dz(1)=z(2)/2
          do j=2,nz-1
              dz(j)=(z(j+1)-z(j-1))/2
          enddo
          dz(nz)=z(nz)-z(nz-1)
      elseif (zlevel.eq.0) then
          z(1)=0.81
          z(2)=2.43
          z(3)=4.05
          z(4)=5.70
          z(5)=7.37
          z(6)=9.07
          z(7)=10.81
          z(8)=12.63
          z(9)=14.54
          z(10)=16.56
          z(11)=18.76
          z(12)=21.19
          z(13)=23.93
          z(14)=27.11
          z(15)=30.88
          z(16)=35.46
          z(17)=41.14
          z(18)=48.34
          z(19)=57.60
          z(20)=69.69
          z(21)=85.65
          z(22)=106.88
          z(23)=135.32
          z(24)=173.63
          z(25)=225.42
          z(26)=295.63
          z(27)=391.02
          z(28)=520.84
          z(29)=697.70
          z(30)=938.87
          z(31)=1267.9
          z(32)=1717.1
          dz(1)=z(2)/2  ! 计算间隔  
          do j=2,nz-1
              dz(j)=(z(j+1)-z(j-1))/2
          enddo
          dz(nz)=z(nz)-z(nz-1)        
      endif
          

          do j=1,nz 
          if (ts_switch.eq.1) then      
              call lintrp(zob,tob,nts,z(j),tinit(j))
              call lintrp(zob,sob,nts,z(j),sinit(j))
C             write(*,*) z(j),tinit(j)
          elseif (ts_switch.eq.2) then 
              sinit(j)=35 ! 认为盐度为常数35
              if (z(j).le.ideal_h0) then
                  tinit(j)=ideal_t0;
              elseif (z(j).le.ideal_h1) then
                  tinit(j)=ideal_t0-(ideal_t1-ideal_t0)/
     x(ideal_h1-ideal_h0)*(z(j)-ideal_h0)
              else    
                  tinit(j)=ideal_t1-(ideal_t2-ideal_t1)/
     x(ideal_h2-ideal_h1)*(z(j)-ideal_h1)
              endif
          endif
          enddo

C         do j=1,nz
CC         z(1) = dzfix1/2. ! 因为是C网格，u、v值和T、S值间开
CC         if(j.gt.1) z(j) = z(j-1) + dz(j) 
C  
C         z(j) = (z1(j-1)+z1(j))/2 
C         call lintrp(zob,tob,nts,z(j),tinit(j))
C         call lintrp(zob,sob,nts,z(j),sinit(j))
C     enddo
  

          do j=1,nz    
              call lintrp(zob,tob,nts,z(j),tinit(j))
              call lintrp(zob,sob,nts,z(j),sinit(j))
          enddo
          
c  initialization of z, t, s, d is finished.

c  calculate the thermal and haline expansion coefficients      

      do 1180 l=1,nz
      t0 = tinit(l)
      s0 = sinit(l)
      d0 = sgt(t0,s0,gg)
      alpha(l) = sgt(t0+0.5,s0,gg) - sgt(t0-0.5,s0,gg)
      beta(l) = sgt(t0,s0+0.5,gg) - sgt(t0,s0-0.5,gg)   
 1180 continue
      
      
      if(nodeep.eq.1) then
      rg = 0.
      rb = 0.
      end if      
      
      call pwpint(nz,rhlat,tob(1),sob(1),rg,rb)
      
c  initialize variables

      do 1 l=1,nz
      do 1 k=1,ny                                                             
      do 1 j=1,nx
      u(j,k,l) = 0.
      v(j,k,l) = 0.
      t(j,k,l) = tinit(l)  ! 所有位置取相同的温度和盐度          
      s(j,k,l) = sinit(l)
                                                                           
      to(j,k,l) = t(j,k,l)                                                          
      so(j,k,l) = s(j,k,l)                                                          
      uo(j,k,l) = 0.
      vo(j,k,l) = 0.
      dtvam(j,k,l)=0.
      dtham(j,k,l)=0.
      dtentm(j,k,l)=0.
    1 continue                     
      
C     do 8844 i=1,nstep
C     timesec(i) = tzero*8.64e4/ti + float(i)*dt7/ti
C     timd(i) = timesec(i)*ti/8.64e4      ! 时间（单位/天）
C8844 continue      
      
      
c     iforcing = 1    compute stress and heat flux fields  
c     iforcing = 2    read stress and heat flux data from a file  

      iforcing = 2
      if (realwind.eq.1) then
        open(unit=51, file='stressfield.dat', status='old')
          do i=1,nxr*nxr
              read(51,*) xi(i),yi(i),uspe(i),vspe(i),tauxi(i),tauyi(i)
          enddo
        close(51)
          do k=1,nx 
              x(k) = xi(k)
              y(k) = yi(1+(k-1)*nx)
C             xoc(k) = x(k)
C             yoc(k) = y(k) 
C             toc(k) = (yoc(j)/hspd)/8.64e4   !注意toc是按天的
          enddo
      elseif (realwind.gt.1) then
          do j=1,nx
              x(j)=-1.*dx*ncx/1000+(j-1)*dx/1000
          enddo
          do k=1,ny
              y(k)=-1.*dy*ncy/1000+(k-1)*dy/1000
          enddo
          
      endif
      
      
      do 10 i=1,nstep
      
      write(*,*) 'time step: ', i

c      timesec(i) = tzero*8.64e4/ti + float(i)*dt7/ti
          timd(i) = dt*i/8.64e4
          yhe = y0 + dt*i*hspd     !  move it along at hspd
          yhkm = yhe/1000.
          xhe = 0.
          yp = y+yhkm

      if (realwind.eq.1) then
        do j=1,nx
          do k=1,ny
              call lintrp(yp,tauxi(j:j+(nx-1)*nx:nx),nx,y(k),taux(j,k)) 
              call lintrp(yp,tauyi(j:j+(nx-1)*nx:nx),nx,y(k),tauy(j,k))
          enddo
        enddo
      
      elseif (realwind.eq.2) then  ! 兰金模型
          do j=1,nx
              do k=1,ny
                  rr = sqrt(x(j)**2+(y(k)-yhkm)**2)
                  taux(j,k)=0
                  tauy(j,k)=0
                  if (rr.le.rmax) then
                      wind=windmax*(rr/rmax)**b 
                      theta=(10*(1+rr/rmax))/180.0*pi
                  else
                      wind=windmax*(rmax/rr)**b
                      if (rr.lt.(1.2*rmax)) then
                          theta=(20+25*(rr/rmax-1))/180.0*pi
                      else
                          theta=25/180.0*pi
                      endif
                  endif
                  
                  windt=wind*cos(theta)
                  windr=-wind*sin(theta)
                  angle=atan2(y(k)-yhkm,x(j))
                  winds=rmax*rr/(rmax**2+rr**2)*hspd*factor
                  windu=-windt*sin(angle)+windr*cos(angle)
                  windv=windt*cos(angle)+windr*sin(angle)+winds ! 加上台风移速影响
                  wind=sqrt(windu**2+windv**2) ! 对风速进行修正，0.8为修正至10m，8-10分钟平均风
                  
                  if (wind.le.11) then   ! 用Oey方法计算拖曳系数
                      cd=1.2e-3
                  elseif (wind.le.19) then
                      cd=0.49e-3+0.065e-3*wind
                  else                  
                      cd=1.364e-3+0.0234e-3*wind-2e-7*wind**2
                  endif   
                  
                  taux(j,k)=1.22*cd*windu*wind
                  tauy(j,k)=1.22*cd*windv*wind
                  
              enddo
          enddo
      
      elseif (realwind.eq.3) then ! SlOSH模型
          do j=1,nx
              do k=1,ny
                  rr = sqrt(x(j)**2+(y(k)-yhkm)**2)
                  taux(j,k)=0
                  tauy(j,k)=0
                  wind=windmax*(2*rmax*rr)/(rmax**2+rr**2)
                  if (rr.le.rmax) then 
                      theta=(10*(1+rr/rmax))/180.0*pi
                  elseif (rr.lt.(1.2*rmax)) then
                      theta=(20+25*(rr/rmax-1))/180.0*pi
                  else
                      theta=25/180.0*pi
                  endif


                  windr=wind*cos(theta)
                  windt=-wind*sin(theta)
                  angle=atan2(y(k)-yhkm,x(j))
                  winds=rmax*rr/(rmax**2+rr**2)*hspd*factor
                  windu=-windr*sin(angle)+windt*cos(angle)
                  windv=windr*cos(angle)+windt*sin(angle)+winds  ! 加上台风移速影响
                  wind=sqrt(windu**2+windv**2) ! 对风速进行修正，0.8为修正至10m，8-10分钟平均风
                  
                  if (wind.le.11) then   ! 用Oey方法计算拖曳系数
                      cd=1.2e-3
                  elseif (wind.le.19) then
                      cd=0.49e-3+0.065e-3*wind
                  else                  
                      cd=1.364e-3+0.0234e-3*wind-2e-7*wind**2
                  endif                 
                  
                  taux(j,k)=1.22*cd*windu*wind
                  tauy(j,k)=1.22*cd*windv*wind
              enddo
          enddo
      endif

 
      if(i.eq.nstep) go to 10   ! 不计算最后一次循环

      do 1010 i79 = 1,ihalf  
      do 1010 ii = 1,ntss
          if(ii.eq.1) then
              dti = 2.*dt
C             dti = 1.*dt
      
              do 80 l=1,nz
              do 80 k=1,ny  
              do 80 j=1,nx 

                  tt(j,k,l) = 0.
                  ut(j,k,l) = 0.
                  vt(j,k,l) = 0.
                  st(j,k,l) = 0.                                               
                  t9(j,k,l) = t(j,k,l) 
                  s9(j,k,l) = s(j,k,l)
                  u9(j,k,l) = u(j,k,l)
                  v9(j,k,l) = v(j,k,l) 
 
                  w(j,k,l) = 0.
                  pp(j,k,l) = 0.
                  tham(j,k,l) = 0.
                  tvam(j,k,l) = 0.
                  tentm(j,k,l) = 0.
                     
   80         continue
C     write(*,*) t(1,1,1),s(1,1,1),u(1,1,1),v(1,1,1) 
          elseif (ii.eq.2) then
   
c  note that dti is half dt. this accounts for summing the tendencies.   
   
              dti = dt/2.                                               
              do 81 l=1,nz
              do 81 k=1,ny
              do 81 j=1,nx
                  to(j,k,l) = t9(j,k,l)
                  so(j,k,l) = s9(j,k,l)
                  uo(j,k,l) = u9(j,k,l)            
                  vo(j,k,l) = v9(j,k,l)
                  w(j,k,l) = 0.
                  pp(j,k,l) = 0.      
   81         continue 
          endif

c  estimate the vertical velocity by integrating divergence from 
c  the surface downward   
      
      do 635 k=2,nym        
      do 635 j=2,nxm
          
          jp = j + 1
          jm = j - 1
          kp = k + 1
          km = k - 1
          
      do 6381 l=1,nz
          divu = dz(l)*ddx*(u(jp,k,l) - u(jm,k,l))
          divv = dz(l)*ddy*(v(j,kp,l) - v(j,km,l)) 
          w(j,k,l) = divu + divv 
 6381 continue

      do 6321 l=2,nz
          w(j,k,l) = w(j,k,l) + w(j,k,l-1)
 6321 continue

  635 continue
      
      iq9 = 0  ! 这个语句用来设置边界条件，非1则为无梯度边界条件
      if(iq9.eq.1) go to 8833      
          do 6328 l=1,nz
              do 633 j=1,nx
                  w(j,1,l) = w(j,2,l)
  633             w(j,ny,l) = w(j,ny-1,l)
              do 634 k=1,ny
                  w(1,k,l) = w(2,k,l)
                  w(nx,k,l) = w(nx-1,k,l)
  634         continue

c  do the corners, just to make sure.

              w(1,1,l) = w(2,2,l)
              w(nx,1,l) = w(nx-1,2,l)
              w(1,ny,l) = w(2,ny-1,l)
              w(nx,ny,l) = w(nx-1,ny-1,l)
 6328     continue
 8833 continue
 
c  compute the pressure perturbation by integrating from depth upwards 

      nzm = nz-1
      do 6327 j=1,nx
      do 6327 k=1,ny
          pp(j,k,nz) = -g*dz(nz)*
     c (alpha(nz)*(t(j,k,nz) - tinit(nz)) + 
     c   beta(nz)*(s(j,k,nz) - sinit(nz))) 

      do 632 l=1,nzm 
          ll = nzm - (l-1)

      pp(j,k,ll) = -g*dz(ll)*
     c (alpha(ll)*(t(j,k,ll) - tinit(ll)) + 
     c   beta(ll)*(s(j,k,ll) - sinit(ll))) 
     c  + pp(j,k,ll+1)    ! 此处认为最深处压强为0

  632 continue
 6327 continue
 
c     write(*,*)'*'
c  compute tendencies on interior points

      do 410 l=1,nz
      do 410 k=2,nym
          kp = k + 1
          km = k - 1                                                                
      do 410 j=2,nxm 
          jp = j + 1 
          jm = j - 1                                                        
                                                                             
          tha = 0. 
          sha = 0.   
          uha = 0.
          vha = 0.  

          uva = 0.
          vva = 0.
          tva = 0.
          sva = 0.
                                                                 
          pgx = 0. 
          pgy = 0.

          if(iupw.eq.1) then
c  use first upwind differencing
              vs = sign(1.,t(j,k,l)) ! sign(x,y)为取x的值，y的符号
              kvs = -1*nint(vs) ! nint(x),将n转为整数，四舍五入
              us = sign(1.,u(j,k,l))
              kus = -1*nint(us)

              dtdy = vs*(t(j,k,l) - t(j,k+kvs,l))/dy ! 若v为正，则取k+1，负则取k-1
              dtdx = us*(t(j,k,l) - t(j+kus,k,l))/dx

              dsdy = vs*(s(j,k,l) - s(j,k+kvs,l))/dy
              dsdx = us*(s(j,k,l) - s(j+kus,k,l))/dx

              dudy = vs*(u(j,k,l) - u(j,k+kvs,l))/dy
              dudx = us*(u(j,k,l) - u(j+kus,k,l))/dx

              dvdy = vs*(v(j,k,l) - v(j,k+kvs,l))/dy
              dvdx = us*(v(j,k,l) - v(j+kus,k,l))/dx
 
              tha = -(v(j,k,l)*dtdy + u(j,k,l)*dtdx)
              sha = -(v(j,k,l)*dsdy + u(j,k,l)*dsdx)
              uha = -(v(j,k,l)*dudy + u(j,k,l)*dudx)
              vha = -(v(j,k,l)*dvdy + u(j,k,l)*dvdx)
 
          elseif (iupw.eq.0) then
c     use central differencing          
          tha = cadv*(u(j,k,l)*(t(jp,k,l) - t(jm,k,l)) + 
     c            v(j,k,l)*(t(j,kp,l) - t(j,km,l)))  ! 即u*dtemp/dx+v*dtemp/dy
          sha = cadv*(u(j,k,l)*(s(jp,k,l) - s(jm,k,l)) +
     c            v(j,k,l)*(s(j,kp,l) - s(j,km,l)))
          uha = cadv*(u(j,k,l)*(u(jp,k,l) - u(jm,k,l)) +
     c            v(j,k,l)*(u(j,kp,l) - u(j,km,l)))
          vha = cadv*(u(j,k,l)*(v(jp,k,l) - v(jm,k,l)) +
     c            v(j,k,l)*(v(j,kp,l) - v(j,km,l)))
      end if 
	 
c  estimate the tendencies due to vertical advection
c  skip if l = 1  (assume there is no vertical advect at surface)

      if(l.eq.1) go to 3322
C         ! 以下求t,s,u,v的垂向梯度
          wavg = 0.5*(w(j,k,l-1) + w(j,k,l)) ! 第l-1层垂向速度

          !  seems better with iupwz = 0;  check the scheme?????  

          if(iupwz.eq.0) then
          lm = l-1
          lp = l+1
          if(l.eq.nz) lp = l  
          dtdz = (t(j,k,lm) - t(j,k,lp))/(z(lm) - z(lp))
          dsdz = (s(j,k,lm) - s(j,k,lp))/(z(lm) - z(lp))
          dudz = (u(j,k,lm) - u(j,k,lp))/(z(lm) - z(lp))
          dvdz = (v(j,k,lm) - v(j,k,lp))/(z(lm) - z(lp))

C  TRY UPWIND DIFFERENCING FOR THE VERTICAL
C  这里给出了垂向的迎风格式，这个格式似乎不稳定

          elseif(iupwz.eq.1) then
              ws = 1.
              lp = l + 1                 !  W POSITIVE
              if(wavg.le.0.) then      
                  lp = l - 1  !  IN CASE W IS DOWN
                  ws = 1.
              end if
              dtdz = ws*(t(j,k,lp) - t(j,k,l))/(z(lp) - z(l))
              stdz = ws*(s(j,k,lp) - s(j,k,l))/(z(lp) - z(l))
              utdz = ws*(u(j,k,lp) - u(j,k,l))/(z(lp) - z(l))
              dvdz = ws*(v(j,k,lp) - v(j,k,l))/(z(lp) - z(l))
          end if

          tva = wavg*dtdz
          sva = wavg*dsdz
          uva = wavg*dudz
          vva = wavg*dvdz

 3322 continue

c  evaluate the pressure gradient 

          pgx = ddx*(pp(jp,k,l) - pp(jm,k,l))
          pgy = ddy*(pp(j,kp,l) - pp(j,km,l))
                                                
          tt(j,k,l) = tt(j,k,l) + tha*hascale + tva*vascale  
          st(j,k,l) = st(j,k,l) + sha*hascale + sva*vascale
          ut(j,k,l) = ut(j,k,l) + uha*hascale + uva*vascale 
     c        - pgx*pscale/1023.               
          vt(j,k,l) = vt(j,k,l) + vha*hascale + vva*vascale 
     c        - pgy*pscale/1023.

          tham(j,k,l) = tham(j,k,l) + tha*hascale/2.
          tvam(j,k,l) = tvam(j,k,l) + tva*vascale/2.
      
  410 continue      
  
C     write(*,*) '%'
      
c  compute the tendency on boundaries by applying
c    a radiation boundary condition

      do 790 l=1,nz
      do 790 k=1,ny,nym  !只有南北边界处
          ka = 1
          if(k.eq.ny) ka = -1
      do 790 j=1,nx
     
          tt(j,k,l) = tt(j,k,l) + crbc*(t(j,k,l) - t(j,k+ka,l))
          st(j,k,l) = st(j,k,l) + crbc*(s(j,k,l) - s(j,k+ka,l))
          ut(j,k,l) = ut(j,k,l) + crbc*(u(j,k,l) - u(j,k+ka,l))
          vt(j,k,l) = vt(j,k,l) + crbc*(v(j,k,l) - v(j,k+ka,l))
 
  790 continue
  
      do 793 l=1,nz
      do 793 k=1,ny
      do 793 j=1,nx,nxm  ! 只有东西边界处                               
          ja = 1 
          if(j.eq.nx) ja = -1 
          
          tt(j,k,l) = tt(j,k,l) + crbc*(t(j,k,l) - t(j+ja,k,l))
          st(j,k,l) = st(j,k,l) + crbc*(s(j,k,l) - s(j+ja,k,l))
          ut(j,k,l) = ut(j,k,l) + crbc*(u(j,k,l) - u(j+ja,k,l))
          vt(j,k,l) = vt(j,k,l) + crbc*(v(j,k,l) - v(j+ja,k,l))   
          
  793 continue
 
      do 510 j=1,nx
      do 510 k=1,ny	
          if(ii.eq.1) then
          do 511 l=1,nz
              tb5(l) = t(j,k,l)
              sb5(l) = s(j,k,l)
              ub5(l) = u(j,k,l)
              vb5(l) = v(j,k,l)    
              tb4(l) = tb5(l)
              sb4(l) = sb5(l)
              ub4(l) = ub5(l)
              vb4(l) = vb5(l)
c	write(*,*)'ll',l
  511         continue

          elseif(ii.eq.2) then
          do 561 l=1,nz
              tb5(l) = to(j,k,l)
              sb5(l) = so(j,k,l)
              ub5(l) = uo(j,k,l)
              vb5(l) = vo(j,k,l)    
              tb4(l) = tb5(l)
              sb4(l) = sb5(l)
              ub4(l) = ub5(l)
              vb4(l) = vb5(l)
  561     continue
          end if
      
C     write(*,*)'#'
      
      ashe(j,k) = 0.          ! as a stop gap，海表热通量
      
      call pwpgo (taux(j,k), tauy(j,k), ashe(j,k), dtm, tb4, sb4,
     x db4, ub4, vb4, z, dz, dml, dtl, sx, sy)
     
C     write(*,*)'@'

c  add  mixing tendencies j
    
      do 512 l=1,nz
          ut(j,k,l) = ut(j,k,l) + (ub4(l) - ub5(l))/dtm 
          tt(j,k,l) = tt(j,k,l) + (tb4(l) - tb5(l))/dtm 
          st(j,k,l) = st(j,k,l) + (sb4(l) - sb5(l))/dtm
          vt(j,k,l) = vt(j,k,l) + (vb4(l) - vb5(l))/dtm
          tentm(j,k,l) = tentm(j,k,l) + 0.5*(tb4(l) - tb5(l))/dtm 
  512 continue
 
          dmlm(j,k) = dml
          dtlm(j,k) = dtl

  510 continue
 
c  apply all tendencies     
  
      do 24 l=1,nz
      do 24 k=1,ny 
      do 24 j=1,nx
          t(j,k,l) = to(j,k,l) + dti*tt(j,k,l)
          s(j,k,l) = so(j,k,l) + dti*st(j,k,l) 
          u(j,k,l) = uo(j,k,l) + dti*ut(j,k,l)
          v(j,k,l) = vo(j,k,l) + dti*vt(j,k,l) 
   24 continue

   
C     write(*,*)'!'
 1010 continue
 
      do j=1,nx
      do k=1,ny
      do l=1,nz
          dtham(j,k,l)=dtham(j,k,l)+tham(j,k,l)*dt
          dtvam(j,k,l)=dtvam(j,k,l)+tvam(j,k,l)*dt
          dtentm(j,k,l)=dtentm(j,k,l)+tentm(j,k,l)*dt
      enddo
      enddo
      enddo
      
      if (mod(i,nn) .eq. 0) then   ! 每nn步输出一个结果
          do l=1,nz
C             do k=1,ny
              write(20,*) ((tham(j,k,l),j=1,nx),k=1,ny)
              write(21,*) ((tvam(j,k,l),j=1,nx),k=1,ny)
              write(22,*) ((tentm(j,k,l),j=1,nx),k=1,ny)
              write(23,*) ((dtham(j,k,l),j=1,nx),k=1,ny)
              write(24,*) ((dtvam(j,k,l),j=1,nx),k=1,ny)
              write(25,*) ((dtentm(j,k,l),j=1,nx),k=1,ny)
              
              write(61,*) ((u(j,k,l),j=1,nx),k=1,ny)
              write(62,*) ((v(j,k,l),j=1,nx),k=1,ny)
              write(63,*) ((w(j,k,l),j=1,nx),k=1,ny)
              write(64,*) ((t(j,k,l),j=1,nx),k=1,ny)
              write(65,*) ((s(j,k,l),j=1,nx),k=1,ny)
C             enddo
          enddo
     
 
          do k=1,ny               ! 只输出表层的值
              write(18,*) (dmlm(j,k),j=1,nx)
              write(19,*) (dtlm(j,k),j=1,nx)
              write(71,*) (u(j,k,1),j=1,nx)
              write(72,*) (v(j,k,1),j=1,nx)
              write(73,*) (w(j,k,1),j=1,nx)
              write(74,*) (t(j,k,1),j=1,nx)
              write(75,*) (s(j,k,1),j=1,nx)
              write(76,*) (taux(j,k),j=1,nx)
              write(77,*) (tauy(j,k),j=1,nx)
          enddo
 
      endif
	
   10 continue   
      
      call exit(1)
      close(19)
      close(20)
      close(21)
      close(22)
      close(23)
      close(24)
      close(25)
      close(26)
      close(27)
      close(28)
      close(29)
      close(30)
      close(31)
      close(32)
      close(33)
      close(34)

      close(51)
      close(52)
      close(53)
      close(54)
      close(55)

      close(61)
      close(62)
      close(63)
      close(64)
      close(65)

      close(71)
      close(72)
      close(73)
      close(74)
      close(75)
      close(76)
      close(77)  
      
      end
      
      
      subroutine maxspd(u,v,nt,ny,x,t,y,hs,us,vs,
     ch2s,u2s,v2s,h3s,u3s,v3s)
c
c  find and list the maximum current speeds and layer depths
c 
      dimension hs(nt,ny),us(nt,ny),vs(nt,ny),u(nt,ny),v(nt,ny),
     ch2s(nt,ny),u2s(nt,ny),v2s(nt,ny),h3s(nt,ny),u3s(nt,ny),
     cv3s(nt,ny),x(nt),t(nt),y(ny)
c
      spmx = 0.
      jmx = 1
      kmx = 1
c
      do 1 j=1,nt
      do 1 k=1,ny
      sp = sqrt(u(j,k)**2 + v(j,k)**2)
      if(sp.gt.spmx) go to 2
      go to 1
    2 continue
      spmx = sp
      jmx = j
      kmx = k
    1 continue
c
      ymx = y(kmx)/1.e3
      xmx = x(jmx)/1.e3
      write (6,3) spmx,xmx,ymx,t(jmx)
      write (1,3) spmx,xmx,ymx,t(jmx)
    3 format (1x,/,1x,'speed max(m/s) =',f6.2,',   at x,y(km) =',2f6.1,
     c',   at time(days) = ',f5.2)
c
      write (6,4) hs(jmx,kmx),us(jmx,kmx),vs(jmx,kmx)
      write (1,4) hs(jmx,kmx),us(jmx,kmx),vs(jmx,kmx)
    4 format (1x,'layer 1; h, u, v, are',3f8.2)
      write (6,5) h2s(jmx,kmx),u2s(jmx,kmx),v2s(jmx,kmx)
      write (1,5) h2s(jmx,kmx),u2s(jmx,kmx),v2s(jmx,kmx)
    5 format (1x,'layer 2; h, u, v, are',3f8.2)
      write (6,6) h3s(jmx,kmx),u3s(jmx,kmx),v3s(jmx,kmx)
      write (1,6) h3s(jmx,kmx),u3s(jmx,kmx),v3s(jmx,kmx)
    6 format (1x,'layer 3; h, u, v, are',3f8.2)
c
      return
      end
c
      subroutine smooth(a,nx,ny)
c
c  apply some smoothing along the grid edges
c
      dimension a(nx,ny),b(60,60)
c  the dimensions of b must match or exceed those of a
      do 3 j=1,nx
      do 3 k=1,ny
    3 b(j,k) = a(j,k)
      do 1 k=1,ny
      b(1,k) = 0.5*(a(1,k) + a(2,k))
      b(nx,k) = 0.5*(a(nx,k) + a(nx-1,k))
    1 continue
      do 2 k=1,ny
      do 2 j=1,5
      j1 = j+1
      j2 = nx - j
      b(j1,k) = 0.25*a(j1-1,k) + 0.5*a(j1,k) + 0.25*a(j1+1,k)
      b(j2,k) = 0.25*a(j2-1,k) + 0.5*a(j2,k) + 0.25*a(j2+1,k)
    2 continue
      do 4 j=1,nx
      do 4 k=1,ny
    4 a(j,k) = b(j,k)
      return
      end
c
c
      subroutine smooth3(a,nx,ny,nz)
c
c  apply some smoothing to the interior only of array a
c
      dimension a(nx,ny,nz),b(60,60)
c
c  the dimensions of b must match or exceed those of a
c
      do 99 l=1,nz
c
c
      nxm = nx-1
      nym = ny-1

      do 1 j=1,nx
      do 1 k=1,ny
c
      jp = j+1
      if(jp.gt.nx) jp=nx-1
      jm = j-1
      if(jm.lt.1) jm=2
      kp = k+1
      if(kp.gt.ny) kp=ny-1
      km = k-1
      if(km.lt.1) km=2
c
      b(j,k) = (0.5*a(j,k,l) + 0.25*(a(jm,k,l) + a(jp,k,l) + 
     x  a(j,km,l) + a(j,kp,l)))/1.5
c
    1 continue
c
      do 4 j=1,nx
      do 4 k=1,ny
    4 a(j,k,l) = b(j,k)
c
   99 continue
      return
      end
c
c
      subroutine hurri (idtau,x,y,xhspd,yhspd,noasym,wspd,hsize,
     c taux,tauy,
     c tdew,tdry,aspd,iglo)

c  this subroutine specifies the hurricane wind stress,
c  and air temperaures, and wind speed as functions of
c  radial distance form (0,0) the storm center.
c  x,y are assumed in meters on input, r is computed in km.
c  hspd is the hurricane translation speed, wspd is the maximum
c  wind speed, both in meters per second
c  hsize is the hurricane size given as radius to maximum winds, km
c  taux,tauy are in mks kinematic units (m**2/sec**2)
c  tdew,tdry are in celsius, aspd is in m/sec
c  noasym is a flag, = 1 shuts off asymmetry of winds due to translation
c                    = 2 shuts off radial wind component also                                        
c  this array (vs) gives ths wind speed in a stationary storm
c  as a function of radial distance
c  data below are from figs. 13.6 and 13.11 of noaa tech
c  report nws 23, meteorological criteria for standard
c  project hurricane .......
c
      dimension vs(15),vsr(15)
      dimension phi(10),phir(10)
c      data vs/0.0,-0.1,-0.5,-0.8,-0.95,-1.0,
c     c-0.97,-0.72,-0.54,-0.44,-0.4,-0.36,
c     c-0.27,-0.23,0.0/
      data vs/0.0,0.1,0.5,0.8,0.95,1.0,
     c 0.97,0.72,0.54,0.44,0.4,0.36,
     c 0.27,0.23,0.0/
c
      data vsr/0.,0.4,0.7,0.8,0.95,1.0,
     c1.35,2.7,4.05,5.4,6.75,8.10,
     c10.8,13.5,40.0/
      nr = 15
c
c  these arrays set the inflow angle, from section 14
c  of the noaa hurricane report
c
      data phi/0.,4.0,7.0,12.0,15.0,
     c22.0,25.5,23.5,21.0,20.0/
      data phir/0.,0.7,1.0,1.25,1.5,
     c2.0,3.0,4.5,6.5,27.0/
      np = 10
c
      rraw = 1.23e-3
c  rraw is the ratio density air/water
c
      taux = 0.                                                                 
      tauy = 0.                                                                 
      aspd = 0.                                                                 
      r = (sqrt(x**2 + y**2))/1.e3                                              
      if(r.gt.1000.) go to 5
c
      rn = r/hsize
      call intx (vsr,vs,nr,rn,vsi)
      call intx (phir,phi,np,rn,phii)
c
      if(noasym.eq.2) phii = 0.
c  the above shuts off all inflow in noasym = 2
c
c  as a test for gloria only, try doubling the inflow angle on the 
c  right side of the storm (to better match the ow stress field)
c  this is better.仅作为对 Gloria 的测试，尝试将风暴右侧的流入角加倍（以更好地匹配 OW 应力场）
c
      if(iglo.eq.1) then
      if(y.gt.-20000.) phii = phii*1.7
      end if
c
      ut = wspd*vsi*cos(phii/57.3)                                       
      ur = -wspd*vsi*sin(phii/57.3)
c 
c      if(idtau.eq.1) then
c      write (3,2210) noasym, x, y, xhspd, yhspd, rn, vsi, phii
c 2210 format (1x,'noasym,  x, y, xhspd, yhspd, rn, vsi, phii',/,
c     x 1x, i5, 7f10.3)
c      write (6,2210) noasym, x, y, xhspd, yhspd, rn, vsi, phii
c      end if
c
c                                                                               
      a = atan2 (y,x)                                                           
      u = ur*cos(a) - ut*sin(a)                                                 
      v = ur*sin(a) + ut*cos(a)   
c
c      if(idtau.eq.1) write (3,2211) ut, ur, y, x, a, u, v
c 2211 format (1x,'ut, ur, y, x, a, u, v',/,1x, 7f11.3)
c      if(idtau.eq.1) write (6,2211)  ut, ur y, x, a, u, v
                                       
c
c
c  add asymmetry to wind field if noasym not equal 1
c
      if(noasym.lt.1) then
c
      uadd = xhspd/2.
      vadd = yhspd/2.
c
c
      u = u + uadd
      v = v + vadd

c
      end if
c
c
      aspd = sqrt(u**2 + v**2)                               
c  drag coefficient is computed according to large and pond
c
      cd = (0.5 + .06*aspd)*1.e-3  
c
c  try constant cd for icd = 1, just to check sensitivity
c
c      icd = 1
c      if(icd.eq.1) cd = 1.3e-3
c
c                                           
      taux = rraw*cd*aspd*u
      tauy = rraw*cd*aspd*v
c
    5 continue 
c
      tdew = 24.
      tdry = 26. 
      if(iglo.eq.1) then
c  below is to get some air/sea exchange for gloria.
      tdew = 23.
      tdry = 24.
      end if
c
      return
      end
c
c
      subroutine thrmo (x,y,timd,sst,aspd,tdew,tdry,ashe,rsi,emp)
c
c  this subroutine computes air/sea heat flux, ashe, and mass
c  flux, emp, and solar insolation rsi.  inputs are sst, air 
c  temperatures, tdew and tdry and wind speed, aspd.
c  units are mks.  fluxes are in mks kinematic units.
c
      rraw = 1.23e-3
      cq = 1.3e-3
      cql = rraw*cq*0.70
      cqs = rraw*cq*0.24
c
      ql = cql*aspd*(tdew - sst)
      qs = cqs*aspd*(tdry - sst)
      ashe = ql + qs
c
c  set the evaporation minus precipitation
c
      p = 0.
      e = -ql/560.
      emp = e - p
      emp = 0.
c
c  set the solar insolation
c
      rsi = 0.
c
      return
      end
c
c
c  interpolate the array s(r) to get value at r = ri
c
      subroutine intx (r,s,n,ri,ui)
      dimension r(10),s(10)                                                       
      nm = n -  1                                                               
      do 2 k=1,nm                                                               
      if(ri.le.r(k+1).and.ri.gt.r(k)) go to 1                                   
    2 continue                                                                  
      return                                                                    
    1 continue                                                                  
      ui = (s(k)*(r(k+1) - ri) + s(k+1)*(ri - r(k)))/(r(k+1) - r(k))            
      return                                                                    
      end                                                                       
c
                                                                    
       subroutine owinit(rlat,rlon)
c
c  this subr reads in the hurricane stress data provided
c  by oceanweather.  call as owinit before the first use
c
c
      dimension tx(43,43),ty(43,43)
      dimension rlat(43),rlon(43)

      dimension a(2,47,45)

      character*30 ifile
      ro = 1.025e3
c
      write (6,3)
    3 format (1x,'enter file name of stress data')
      read (5,4) ifile
    4 format (a20)
c
      open (unit=30,file=ifile,form= 'unformatted',
     cstatus='old')
c
c
      read (30) a
 
c
c
      do 80 j=1,43
      jp = 43 + 2 - j
      rlat(j) = a(1,10,jp)
      kp = j + 2
      rlon(j) = a(2,kp,10)
c
      write (6,81) j,jp,kp,rlat(j),rlon(j)


   81 format (1x,'lat, lon',3i5,2f10.3)
c
   80 continue
c
c
c
      return
c
c  use this entry point for one snap shot of data
c
      entry owget(tx,ty)
c
c
      read (30) a
c
      do 5 j=1,43
      do 5 k=1,43
c
      jp = 43 + 2 - j
      kp = k+2
c
      tx(j,k) = -1.*a(2,kp,jp)/ro
      ty(j,k) = a(1,kp,jp)/ro
c
    5 continue
      txq = 1.023e3*tx(16,22)
      tyq = 1.023e3*ty(16,22)
      write (6,89) txq,tyq,rlat(16),rlon(22)
   89 format (1x,'tx, ty near center',4f10.4)
c
c
      return
      end
      subroutine rhot(ang,x,y)
c
c  rotate x,y into a unknown corrdinate system which
c  is at an angle ang to old system (radians, ccw > 0)
c
      ct = cos(ang)
      st = sin(ang)
      xs = x
      ys = y
      x = xs*ct + st*ys
      y = ys*ct - st*xs
      return
      end
      
      subroutine lintrp(r,s,n,ri,si)
      dimension r(120),s(120)
      nm = n - 1
      do 2 k=1,nm
      if(ri.lt.r(k+1).and.ri.ge.r(k)) go to 1
    2 continue
C     write (6,3) r(1),r(n),ri
C   3 format (1x,'times are bad in interp; rf,rl,ri',3e12.5)
      return
    1 continue
      si = (s(k)*(r(k+1) - ri) + s(k+1)*(ri - r(k)))/(r(k+1) - r(k))
      return
      end
c
      subroutine spdr(u,v,nx,ny,spd,dir)
c
c  calculate speed and direction from current components
c
      dimension u(nx,ny),v(nx,ny),spd(nx,ny),dir(nx,ny)
c
      do 1 j=1,nx
      do 1 k=1,ny
c
      spd(j,k) = sqrt(u(j,k)**2 + v(j,k)**2)
      ujk = -u(j,k)
      if(abs(ujk).lt.1.e-6) ujk = ujk + sign(1.e6,ujk)
c
      dir(j,k) = atan2(v(j,k),ujk)*57.3
      if(dir(j,k).lt.0.) dir(j,k) = dir(j,k) + 360.
c
    1 continue
      return
      end
c
      subroutine pwph
c
c  this subroutine is an implementation of the price, weller,
c  pinkel upper ocean model described in jgr 91, c7 8411-8427
c  july 15, 1986. this version was coded and documented by
c  jim price in february 92 to use in the 3-d hurrricane model.
c  please direct questions and comments to jim price, whoi,
c  woods hole, ma, 02543, usa, tel. 508-548-1400, x2526.
c
c  all units within the model are mks.
c
c  the array absrb (gone) is internal to the subroutine, and
c  must have dimension larger than nz, the number of points
c  in the t,s,u,v,d arrays. note that t,s,u,v,d are dimensioned 
c  in the main or calling program.
c
      dimension t(40), s(40), u(40), v(40), d(40), dz(40), z(40)
      dimension ub4(40), vb4(40), sx(40), sy(40)
c
      save
c
c  come here once to initialize the model.
c
      entry pwpint(nz4,rlat4,tr4,sr4,rg4,rb4)
c
      nz = nz4
      rlat = rlat4
      sr = sr4
      tr = tr4
      rg = rg4
      rb = rb4
c
c
c  nz is the number of grid points in the vertical.
c  dz is the grid interval (meters) (an array).
c  beta1 and beta2 are the longwave and shortwave extinction
c  coefficients for solar insolation (meters).
c  rlat is the latitude (deg).
c  tr, sr are the reference values of temperature and 
c  salinity used in the density equation (c, and ppt).
c  rg is the critical gradient richardson number (= 0.25, 
c  in the pwp model).
c  rb is the critical bulk richardson number (= 0.65, 
c  in the pwp model).
c
c  set a few constants that do not change during the run.
c
      d2r = 2.*3.14159/360.
      g = 9.8
      ro = 1.024e3
      cpw = 3990. 
      f = 2.*7.29e-5*sin(rlat*d2r)
c
c  evaluate alpha and beta, the thermal and haline expansion
c  coefficients used in a linear equation of state.
c
      dr = sgt(tr,sr,gg)
      alpha = sgt(tr+.5,sr,gg) - sgt(tr-.5,sr,gg)
      beta = sgt(tr,sr+.5,gg) - sgt(tr,sr-.5,gg)
c
      write (6,3344) alpha,beta,rb,rg,em0
 3344 format (1x,'alpha, beta, rb, rg, em0',2e10.2,3f9.2)
c
c
c  finished with initialization.
c
      return
c
c
c  come here when ready to time step ahead by an amount dt.
c
      entry pwpgo (tx,ty,qn,dt,t,s,d,u,v,z,dz,dml,dtl,sx,sy)
c
c
c  input to the subroutine is the following:
c
c  tx and ty are the east and north stress components (pa).
c  qn is the net heaqt flux (watts/m*m)
c  dt is the time step (sec).
c  t, s, d, u, v are arrays that store temperature, salinity,
c  density and east and north current -  these arrays also return
c  the unknown values.
c  z is the mid depth of a grid cell, dz is the thickness of the
c  cell, thus the total depth or thickness to cell j is z(j) + dz(j)/2.
c
c  output only are dml and dtl, the depths of the mixed layer
c  and the transition layer (meters).
c
c  the heat and salt budgets used to be here. gone now.
c
c  compute the density (sigma units), and relieve static
c  instability, if it occurs.
c
      t(1) = t(1) + dt*qn/(dz(1)*4.e6)
      do 71 j=1,nz
      d(j) = dr + (t(j) - tr)*alpha + (s(j) - sr)*beta
c     if(j.lt.5) write (6,776) j,t(j),s(j),d(j)
c  776 format (1x,'density ; t,s,d',i6,3f10.4)
   71 continue
      
      call mldep(d,nz,ml)
      ds1 = z(ml)  + dz(ml)/2.
      
      call si(d,dz,nz,ml)
      
      if(ml.gt.1) call mix5(t,s,u,v,d,dz,ml)
      ds2 = z(ml) + dz(ml)/2.
            
c  at this point the density profile is statically stable.
c
c  time step the momentum equation.
c
c  rotate the current throughout the water column through an
c  angle equal to inertial rotation for half a time step.
c
      ang = f*dt/2.
      sa = sin(ang)
      ca = cos(ang)
 
      do 32 j=1,nz
      call rot(u(j),v(j),sa,ca)
      ub4(j) = u(j)
      vb4(j) = v(j)
      sx(j) = 0.
      sy(j) = 0.
   32 continue

c  apply the wind stress to the mixed layer only.
c
c  find the ml depth.

      call mldep(d,nz,ml)      ! 这句话是否多余？
      dml = z(ml) + dz(ml)/2.
      
      du = (tx/(dml*ro))*dt
      dv = (ty/(dml*ro))*dt
      do 33 j=1,ml
      u(j) = u(j) + du
      v(j) = v(j) + dv
   33 continue

c  cause the mixed-layer to deepen.
c
c  will skip this section if rb = rg = 0.
c
      nzm = nz - 1
      call mldep(d,nz,ml)     ! 是否又多余了？
      dml1 = z(ml) + dz(ml)/2.
      
      if(rb.le.1.e-4) go to 81
c
c  come here to do the bulk richardson number instability
c  form of mixing (as in pwp).
c
      rvc = rb
      do 80 j=ml,nzm
      delr = (d(j+1) - d(j))/ro
      ds2 = (u(j+1) - u(j))**2 + (v(j+1) - v(j))**2
      ds2 = ds2 + 1.e-8       ! 为防止0值出现
      h1 = z(j) + dz(j)/2.
      rv = g*delr*h1/ds2
c     if(j.eq.ml) write (6,778) delr,h1,ds2,rv
C 778 format (1x,'delr,h1,ds2,rv,',4f10.2)
      if(rv.gt.rvc) go to 81
c
c  rv is critical, come here to entrain.
c
      jp = j + 1
      
      call mix5(t,s,u,v,d,dz,jp)
      
   80 continue
   81 continue
      
c  relieve gradient richardson number instability if it occurs.
c
      if(rg.gt.1.e-4) call rimix(rg,u,v,d,t,s,dz,nz,nmix,kcd)
      
c  estimate the stress profile by integrating the acceleration 
c  from depth upward. 
c 
      do 90 j = 1,nzm
      jz = nzm - (j-1)
      sx(jz) = sx(jz+1) + ro*dz(jz)*(u(jz) - ub4(jz))/dt  ! 这里需要最后一层力为0
      sy(jz) = sy(jz+1) + ro*dz(jz)*(v(jz) - vb4(jz))/dt
   90 continue
c
c  rotate the currents through the rest of the time step.
c
      do 34 j = 1,nz
      call rot(u(j),v(j),sa,ca)
   34 continue
c
      call mldep(d,nz,ml)
      if(kcd.eq.0) kcd = 1
      if(ml.eq.0) ml = 1
      dml = z(ml) + dz(ml)/2.
      dtl = z(kcd) + dz(kcd)/2.
      if(dtl.lt.dml) dtl = dml
c
      return
      end
c
c
c
      subroutine rimix(rg,u,v,d,t,s,dz,nz,nmix,kcd)
c
c  this subroutine performs the gradient richardson number
c  relaxation by mixing adjacent cells just enough to bring
c  them to a unknown richardson number = runknown (see subroutine
c  stir for details of this).
c
c  all of u, v, d, t, s, dz, nz are exactly as in the pwp subroutine.
c  nmix is the number of iterations required to achieve the relaxation,
c  kcd is the deepest grid level reached (the base of the transition
c  layer if the mixing is ordinary wind-driven mixing).
c
      dimension u(nz), v(nz), d(nz), t(nz), s(nz), dz(nz)
      dimension r(nz)
      rc = rg
      g = 9.8
      d0 = 1.023e3
      nzm = nz - 1
      kcd = 0
      nmix = 0
c
      j1 = 1
      j2 = nzm
   10 continue
c
c  compute the gradient richardson number, taking care
c  to avoid dividing by zero in the mixed-layer. these
c  limits are arbitrary, and have caused problems in the
c  past in some applications. 
c 
      do 1 j = j1,j2
      dd = d(j+1) - d(j)
      if(dd.lt.1.e-3) dd = 1.e-3
      dv = (u(j+1) - u(j))**2 + (v(j+1) - v(j))**2
      if(dv.lt.1.e-6) dv = 1.e-6
      dza = (dz(j) + dz(j+1))/2.
      r(j) = g*dza*dd/(dv*d0)
    1 continue
c
c  find the smallest value of r in profile.
c
      rs = r(1)
      js = 1
      do 2 j=2,nzm
      if(r(j).lt.rs) go to 3
      go to 2
    3 continue
      rs = r(j)
      js = j
    2 continue
c
c  check to see whether the smallest r is critical or not.
c
      if(rs.gt.rc) go to 99
c
c  mix the cells js and js+1 that had the smallest richardson number.
c
      if(js.ge.kcd) kcd = js + 1
      call stir(rs,t,dz,js)
      call stir(rs,s,dz,js)
      call stir(rs,d,dz,js)
      call stir(rs,u,dz,js)
      call stir(rs,v,dz,js)
      nmix = nmix + 1
c
c  recompute the richardson number over the part of the profile
c  that has changed.
c
      j1 = js - 2
      if(j1.lt.1) j1 = 1
      j2 = js + 2
      if(j2.gt.nzm) j2 = nzm
      go to 10
c
   99 continue
      return
      end
c
c
      subroutine stir (r,a,dz,j)
c
c  this subroutine mixes cells j and j+1 just enough so that
c  the richardson number after the mixing is brought up to 
c  the value runknown. in order to have this mixing process
c  converge, runknown must exceed the critical value of the 
c  richardson number where mixing is presumed to start. if
c  r critical is 1/4, then runknown = 0.3 would be reasonable.
c 
c  r is the richardson number before mixing. 
c  a is the array to be mixed.
c  j is the cell where r was found to be critical.
c  b is the factor that accounts for variabel dz.
c
      dimension a(40),dz(40)
      runknown = 0.3
c
      f = 1. - r/runknown
      da = (a(j+1) - a(j))*f/2.
      b = 2./(1. + dz(j)/dz(j+1))
      a(j+1) = a(j+1) - b*da*dz(j)/dz(j+1)
      a(j) = a(j) + b*da
      return
      end
c
c
      subroutine mldep(d,nz,ml)
c
c  this subroutine scans through the density array d to find
c  the index of the cell containing the surface mixed-layer.
c
      dimension d(nz)
c
      deps = 1.e-4
      nzm = nz - 1
      do 1 j=1,nzm
      dd = abs(d(j+1) - d(j))
      if(dd.gt.deps) go to 2
    1 continue
    2 continue
      ml = j
      return
      end
c
c
      subroutine mix1(d,dz,k)
c
c  this subroutine completely mixes the array d down to level k.
c
      dimension d(40),dz(40)
      ds = 0.
      zs = 0.
      do 1 j=1,k
      ds = ds + d(j)*dz(j)
      zs = zs + dz(j)
    1 continue
      ds = ds/zs
      do 2 j=1,k
      d(j) = ds
    2 continue
      return
      end
c
c
      subroutine mix5(t,s,u,v,d,dz,k)
c
c  this subroutine mixes arrays t, s, u, v, d down to level k.
c
      dimension t(40), s(40), u(40), v(40), d(40), dz(40)
      ts = 0.
      ss = 0.
      us = 0.
      vs = 0.
      ds = 0.
      zs = 0.
      do 1 j=1,k
      ts = ts + t(j)*dz(j)
      ss = ss + s(j)*dz(j)
      us = us + u(j)*dz(j)
      vs = vs + v(j)*dz(j)
      ds = ds + d(j)*dz(j)
      zs = zs + dz(j)
    1 continue
      do 2 j=1,k
      u(j) = us/zs
      v(j) = vs/zs
      t(j) = ts/zs
      s(j) = ss/zs
      d(j) = ds/zs
    2 continue
      return
      end
c
c
      subroutine rot(u,v,sa,ca)
c
c  this subroutine rotates the vector (u,v) through an
c  angle whose sine and cosine are sa, ca.
c
      u0 = u
      v0 = v
      u = u0*ca + v0*sa
      v = v0*ca - u0*sa
      return
      end
c
c
      subroutine si(d,dz,nz,ml)        
      dimension d(nz),dz(nz)
      dimension di(nz)
      
c  find and relieve static instability that may occur in the
c  density array d.
c  ml is the depth of the surface mixed layer after adjustment.
c  sipe is the change in potential energy due to free convection.

      do 10 j=1,nz
      di(j) = d(j)
   10 continue
      
      mml = 1
      do 1 j=2,nz
      if(d(j).gt.d(j-1)) go to 2
      if(d(j).lt.d(j-1)) then
      mml = j
      call mix1(d,dz,mml)
      end if
    1 continue
    2 continue
      ml = mml 
      return
      end

      function sg0(s) 
c 
c  a sigma-0 subroutine neede by the sigma-t subroutine
c  taken from seaprop.
c
c  sigma-0 knudsen 
      sg0 = ((6.76786136e-6*s-4.8249614e-4)*s+0.814876577)*s  
     x -0.0934458632
      return
      end
c
c
      function sgt(t,s,sg)
c
c  a sigma-t subroutine taken from seaprop.
c
c sigma-t knudsen 
      sg = sg0(s) 
   20 sgt = ((((-1.43803061e-7*t-1.98248399e-3)*t-0.545939111)*t
     x +4.53168426)*t)/(t+67.26)+((((1.667e-8*t-8.164e-7)*t 
     x +1.803e-5)*t)*sg+((-1.0843e-6*t+9.8185e-5)*t-4.7867e-3)*t
     x +1.0)*sg 
      return
      end


      