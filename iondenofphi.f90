! Read in an ion distribution function in the chiofv format
! Or form the distribution using the gpar parameters.
! Calculate n(phi) for a steady potential fixed in a frame (hole)
! moving at speed vf, by integrating over the distribution. 
! Form the corresponding force ion on the hole by integrating
! over a hole shape. 

module idenofphi
integer, parameter :: nv=600,nvp=100
integer, parameter ::  npsimax=20,nvsmax=60
real, dimension(0:nv+2) :: f,v
real, dimension(1:nv+2) :: fprime,vprime
real, dimension(1:nv+1,npsimax) :: forcea,delfrcdelx,fva,fvthresh
real :: psi,psimax=.5,vh=0.,den,delden,v0,dv,force,fperx,f0
integer :: nvin,icross
real, dimension(nvp) :: dena,vpa,deldena,x,phia,dena2
integer, parameter :: npar=4, ngmax=10
integer :: ng=-1
real, dimension(npar,ngmax) :: gpar
real :: theta=0.,vt2=1.
real, dimension(nvsmax,npsimax) :: fperxa
real, dimension(nvsmax) :: vsa
real, dimension(npsimax) :: psia
integer, dimension(npsimax) :: icrsa,jtha
integer :: nvs=nvsmax,npsi=npsimax,imultiscan=0,nframes=2
real :: vt2set2=1.,denset2=.5,vt2set=1.,Tthresh=0.,vsmax=2.5
logical :: lplotfv=.false.,lstab=.false.,lscanplot=.true.
logical :: lmain=.true.,lequil=.false.,ldebug=.false.,lsech4=.true.
character*30 :: label=''
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine denprimeint(vp,ldebug)
! Integrate to get density of vp using fprime and parts expression.
! The velocity vp has magnitude sqrt(2.*phi), and sign equal to the
! sign of the incoming velocity that reflects at the point
! considered. That is, vp is the velocity relative to the hole that
! reflects at the point. vp+vh is therefore the incoming velocity that
! reflects there. The distant potential is zero on both sides.
    logical ldebug
    logical endcor
    endcor=.false. ! endcor breaks symmetric cases
!    endcor=.true. ! and does not seem to do much for asymmetric
    vpsi=sqrt(2.*psi)                    ! relative speed that reaches psi
    vksm=max((-vpsi+vh-vprime(1))/dv+1,1.)  ! real caret index of -vpsi
    vksp=min((vpsi+vh-vprime(1))/dv+1,nv+0.)! +vpsi. +1 for rounding down.
    vkp =min((abs(vp)+vh-vprime(1))/dv,float(nv))+1     ! index of +|vp|
    vkpm=max(-abs(vp)+vh-vprime(1),0.)/dv+1             ! -|vp|
    vkfm=max(-vpsi+vh-v(0),0.)/dv+1            ! real caret f-index of -vpsi
    vkfp=min((vpsi+vh-v(0))/dv-1.,float(nv))   ! real caret f-index of +vpsi
    if(vkfm.lt.0.or.vkfp.gt.nv+1)write(*,'(10f8.3)')vkfp,vkfm,vpsi,vh,v(0)
    ! Interpolate integrated f(+-vpsi) terms and add
    fvpsim=f(int(vkfm))+(vkfm-int(vkfm))*(f(int(vkfm)+1)-f(int(vkfm)))
    fvpsip=f(int(vkfp))+(vkfp-int(vkfp))*(f(int(vkfp)+1)-f(int(vkfp)))
    if(ldebug)write(*,*)'vkfp   fvpsip     vkfm     fvpsim,    vh    &
         &  v(0)      vpsi    vp'
    if(ldebug)write(*,'(8f9.4)')vkfp,fvpsip,vkfm,fvpsim,vh,v(0),vpsi,vp
    delden=0.
    den=0.
    den=sign(1.,vp)*(fvpsip-fvpsim)*sqrt(abs(vpsi**2-vp**2))
    do i=1,int(vksm)                    !Transmitted/Unreflected
       vrel=vprime(i)-vh
       delden=delden+fprime(i)*(vrel-sign(sqrt(max(vrel**2-vp**2,0.)),vrel))*dv
       den=den+fprime(i)*sqrt(max(vrel**2-vp**2,0.))*dv
    enddo
    if(endcor) den=den+(vksm-int(vksm)-0.5)*(fprime(int(vksm))&
         *sqrt(max((vprime(int(vksm))-vh)**2-vp**2,0.)))*dv ! End correction
    if(.not.abs(delden).ge.0)stop  'delden NAN 1'
    do i=int(vksp)+1,nvin+1             !Unreflected/Transmitted
       vrel=vprime(i)-vh
       delden=delden+fprime(i)*(vrel-sign(sqrt(max(vrel**2-vp**2,0.)),vrel))*dv
       den=den-fprime(i)*sqrt(max(vrel**2-vp**2,0.))*dv
    enddo
    if(endcor) den=den+(vksp-int(vksp)-0.5)*(fprime(int(vksp+1))&
         *sqrt(max((vprime(int(vksp+1))-vh)**2-vp**2,0.)))*dv ! End correction
    if(.not.abs(delden).ge.0)stop  'delden NAN 2'
    if(vp.gt.0)then            
       do i=int(vkp)+1,min(int(vksp),nvin+1)   !Reflected
          vrel=vprime(i)-vh
          delden=delden+2.*fprime(i)*(vrel-sign(sqrt(max(vrel**2-vp&
               &**2,0.)),vrel))*dv
          den=den-2.*fprime(i)*sqrt(max(vrel**2-vp**2,0.))*dv
       enddo
       if(endcor) den=den+(vkp-int(vkp)-0.5)*(fprime(int(vkp+1))&
         *sqrt(max((vprime(int(vkp+1))-vh)**2-vp**2,0.)))*dv ! End correction
       if(endcor) den=den+(vksp-int(vksp)-0.5)*(fprime(min(int(vksp),nvin+1))&
            *sqrt(max((vprime(min(int(vksp),nvin+1))-vh)**2-vp**2,0.)))*dv
    else
       do i=int(vksm)+1,int(vkpm)       !Reflected
          vrel=vprime(i)-vh
          delden=delden+2.*fprime(i)*(vrel-sign(sqrt(max(vrel**2-vp&
               &**2,0.)),vrel))*dv
          den=den+2.*fprime(i)*sqrt(max(vrel**2-vp**2,0.))*dv
       enddo
       if(endcor) den=den+(vksm-int(vksm)-0.5)*(fprime(int(vksm+1))&
         *sqrt(max((vprime(int(vksm+1))-vh)**2-vp**2,0.)))*dv ! End correction
       if(endcor) den=den+(vkpm-int(vkpm)-0.5)*(fprime(int(vkpm))&
            *sqrt(max((vprime(int(vkpm))-vh)**2-vp**2,0.)))*dv ! End correction
    endif
    do i=int(vkpm)+1,int(vkp)           !Absent background.
       vrel=vprime(i)-vh
       delden=delden+fprime(i)*vrel*dv
    enddo
    if(.not.abs(delden).ge.0)stop  'delden NAN'
    if(ldebug)then
       write(*,*)'denprimeint vp=',vp,' den',den
       write(*,*)'    1,vksm,vksm+1,vkpm,  vkp, vksp,vksp+1,nvin+1'
       write(*,'(10i6)')1,int(vksm),int(vksm)+1,int(vkpm)&
            &,int(vkp),min(int(vksp),nvin+1),int(vksp)+1,nvin+1
       write(*,'(10f6.3)')f(1),f(int(vksm)),f(int(vksm)+1)&
            &,f(int(vkpm)) ,f(int(vkp)) ,f(min(int(vksp),nvin +1))&
            &,f(int(vksp)+1) ,f(nvin+1)
       call autoplot(v,f,nv)
       call axlabels('v','f')
       call axis2
       call polymark([v(1),v(int(vksm)),v(int(vksp)+1),v(nvin+1)],&
            & [f(1),f(int(vksm)),f(int(vksp)+1) ,f(nvin+1)],4,1)
       call polymark([v(max(int(vksm),1)+1),v(int(vkpm))&
            &,v(int(vkp)),v(min(int(vksp),nvin+1))],&
            & [f(max(int(vksm),1)+1),f(int(vkpm)),f(int(vkp)) &
            &,f(min(int(vksp),nvin +1))],4,2)
       call pltend
    endif
!write(*,'(a,4f8.4,8f6.1)')'denprimeint',den,vp,vh,vpsi,vksm,vksp,vkpm,vkfm,vkfp
  end subroutine denprimeint
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine denint(vp,ldebug)
    logical ldebug
! Find the density at potential corresponding to vp by integrating over
! the distribution function f, given peak psi supposed at a v(i)location. 
! Leave result in den.
    vpsi=sqrt(2.*psi)
    ivpsim=int(max((-vpsi+vh-v(0))/dv+.5,0.     ))
    if(v(ivpsim)**2.le.vp**2)ivpsim=ivpsim-1         !Hack
    ivpsip=int(min(( vpsi+vh-v(0))/dv+.5,nvin+0.))
    if(v(ivpsip)**2.le.vp**2)ivpsip=ivpsip+1
    vpk=min(max((vp-v(0))/dv,0.),float(nvin))
! ivp is the first index in the |v|>vp region. Its contribution is obtained
! by integrating from vp to [ivp+0.5]*dv using the constant value f.
! That is f*sqrt((vp+dv*(0.5+dvp))**2-vp**2). The singularity addition.
! The central f-value to be used is between ivp-0.5 and ivp+0.5
! i.e. f=[f(ivp)+dvp*f(ivp-1)+(1-dvp)*f(ivp+1)]/2. But arguably it is 
! weighted toward the low end, and might be [dvp*f(ivp-1)+(1-dvp)*f(ivp)]
    if(vp.ge.0)then
       ivp=int(vpk+1.)
       dvp=ivp-vpk
    else
       ivp=int(vpk)
       dvp=vpk-ivp
    endif
    den=0.
    do i=0,ivpsim                    !Transmitted/Unreflected
       if(i.lt.ivpsim)then
          den=den+f(i)/sqrt(max((v(i)-vh)**2-vp**2,0.))*abs(v(i))*dv
       else
          den=den+0.5*f(i)/sqrt(max((v(i)-vh)**2-vp**2,0.))*abs(v(i))*dv
       endif
    enddo
    do i=ivpsip,nvin                  !Unreflected/Transmitted
       if(i.gt.ivpsip)then
          den=den+f(i)/sqrt(max((v(i)-vh)**2-vp**2,0.))*abs(v(i))*dv
       else
          den=den+0.5*f(i)/sqrt(max((v(i)-vh)**2-vp**2,0.))*abs(v(i))*dv
       endif
    enddo
    if(ldebug)write(*,*)den,nvin,dv
    if(ldebug)write(*,*)den
    if(vp.ge.0)then            
       do i=ivp,ivpsip   !Reflected
          if(i.eq.ivp)then  ! Not much to choose between these.
!             den=den+2.*(dvp*f(ivp-1)+f(ivp)+(1.-dvp)*f(ivp+1))/2.& 
             den=den+2.*(dvp*f(ivp-1)+(1.-dvp)*f(ivp)) &
             *sqrt((abs(vp)+dv*(0.5+dvp))**2-vp**2) ! Singular end
          elseif(i.lt.ivpsip) then
             den=den+2.*f(i)/sqrt(max((v(i)-vh)**2-vp**2,0.))*abs(v(i))*dv
          else
             den=den+f(i)/sqrt(max((v(i)-vh)**2-vp**2,0.))*abs(v(i))*dv
          endif
       enddo
    else
       do i=ivp,ivpsim,-1       !Reflected
          if(i.eq.ivp)then
             den=den+2.*(dvp*f(ivp+1)+(1.-dvp)*f(ivp)) &
                  *sqrt((abs(vp)+dv*(0.5+dvp))**2-vp**2)  ! Singular end
          elseif(i.gt.ivpsim)then
             den=den+2.*f(i)/sqrt(max((v(i)-vh)**2-vp**2,0.))*abs(v(i))*dv
          else
             den=den+f(i)/sqrt(max((v(i)-vh)**2-vp**2,0.))*abs(v(i))*dv
          endif
       enddo
    endif
    if(.not.abs(den).lt.1.e30)then
       write(*,*)'denint ERROR',den
       write(*,*)ivpsim,ivpsip,ivp,i
       write(*,*)vp,v(i),f(i)
       stop
    endif
    if(ldebug)then
       if(ldebug)write(*,*)den
       write(*,*)'denint      vp=',vp,' den',den
       write(*,'(3i6,5f9.3)')ivpsim,ivp,ivpsip,vpk,dvp
       write(*,'(10f6.3)')f(ivpsim),f(ivp),f(ivpsip)
       call autoplot(v,f,nv)
       call axlabels('v','f')
       call axis2
       call polymark([v(ivpsim),v(ivp),v(ivpsip)],&
            & [f(ivpsim),f(ivp),f(ivpsip)],3,1)
!       call polymark(v(ivpalt),f(ivpalt),1,2)
       call pltend
    endif

  end subroutine denint
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine readfvin
    open(13,file='COPTICverif.dat',status='old',err=101)
    read(13,*) nvin
    read(13,*) dv
    read(13,*) v1
    if(nvin.gt.nv)then
       write(*,*) 'Input f longer than expected.','nvin=',nvin
       call exit(1)
    endif
    read(13,*) (f(i),i=1,nvin)
    f(0)=0.
    f(1)=0.
    f(nvin+1)=0.
    f(nvin+2)=0.
    do i=0,nvin+2
       v(i)=(v1+(i-1.)*dv)
    enddo
    close(13)
! Calculate derivative fprime and the corresponding vprime
    do i=1,nvin+2
       fprime(i)=(f(i)-f(i-1))/dv
       vprime(i)=(v(i)+v(i-1))/2.
    enddo
    return
101 write(*,*)'Error opening file'
    close(14,status='delete')
    call exit(1)
  end subroutine readfvin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!********************************************************************
! Return the 1-D distribution function in array f(0:nv+2) for nv+3
! velocity points from v0 to v0+(nv+2)*dv being the sum of ng gaussians
! whose parameters are in gpar(*,ig)=(den, vel, vtperp, vtparallel) 
! At an angle of velocity component to drift velocity, theta.
! vt2h= the temperature in scaled units of the Gaussian component.
  subroutine fvinit
! Parameters; den, vel, vtperp, vtparallel.          
    nvin=nv
    do i=0,nv+2
       f(i)=0.
    enddo
    ct=cos(theta)
    st=sin(theta)
    do k=1,ng
       vt2h=gpar(3,k)*st**2+gpar(4,k)*ct**2
       fnorm=1/(sqrt(2.*3.141593*vt2h))
       do i=0,nv+2
          vi=v0+i*dv
          if(k.eq.1)v(i)=vi
          f(i)=f(i)+gpar(1,k)*fnorm*exp(-(vi-gpar(2,k)*ct)**2/(2.*vt2h))
       enddo
    enddo
! Calculate derivative fprime and the corresponding vprime
    do i=1,nv+2
       fprime(i)=(f(i)-f(i-1))/dv
       vprime(i)=(v(i)+v(i-1))/2.
    enddo
  end subroutine fvinit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ForceOfX(psiin,vhin,ioff)  
! Calculate density dena(i) and integrate times dphidx to get force.
! -ioff is density-index-shift (cf potential). This is equivalent to
! moving the entire density profile forward by ioff*dx, or the entire
! potential profile by delta x=-dx*ioff. 
! The xwidth of the profile is chosen to make \int phi dx =\psi
  real :: psiin,vhin  
  psi=psiin
  if(lsech4)then
     xwidth=3/4.
     xmax=2.5*xwidth
  else
     xwidth=1/sqrt(3.1415926)
     xmax=3.*xwidth
  endif
  vh=vhin
  force=0.
  f0=0.
  dx=2.*xmax/(nvp-1.)
  do i=1,nvp
     x(i)=xmax*(-1.+2.*(i-1.)/(nvp-1.))
     if(lsech4)then
        phia(i)=psi/cosh(x(i)/xwidth)**4
     else
        phia(i)=psi*exp(-(x(i)/xwidth)**2)
     endif
     vpa(i)=-sign(sqrt(2.*phia(i)),x(i))
     call denprimeint(vpa(i),.false.)
     dena(i)=den
  enddo
  do i=max(3,2+ioff),min(nvp,nvp+1+ioff)
     dphidx=(phia(i)-phia(i-2))/(2.*dx)  ! Equal x space
     force=force+(1.-dena(i-1-ioff))*dphidx*dx ! Accumulate offset force on ions
     f0=f0+(1.-dena(i-1))*dphidx*dx            ! Accumulate force on ions f0
!     if(ldebug)write(*,*)i,x(i-1),dphidx,1.-dena(i-1-ioff),force
  enddo
  fperx=(force-f0)/(-dx*ioff-1.e-6) ! Force difference per phi-displacement.
end subroutine ForceOfX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine parseargs
    character(30) :: argument
! ng indicates: -1 internal iteration, 0 read in file, >0 cmdline set
    vw=0.
    vdmax=0.
    vdmin=0.
    vwfac=5.5
! Defaults:
    gpar=1.
    gpar(2,:)=0.
    psi=psimax
    call pfset(3)    ! Default save plot files.
! Read commandline arguments and adjust parameters.
    do ia=1,iargc()
       call getarg(ia,argument)
       if(argument(1:2).eq.'-g')then
          if(ng.lt.0)ng=0                         ! No iteration  
          read(argument(3:),*,end=106)gpar(1,ng+1)! Null -g ignore.
          ng=ng+1                                 ! Read g
          read(argument(3:),*,err=103,end=105)(gpar(k,ng),k=1,4)
105       continue
          if(gpar(2,ng).lt.vdmin)vdmin=gpar(2,ng)
          if(gpar(2,ng).gt.vdmax)vdmax=gpar(2,ng)
          if(vwfac*gpar(4,ng).gt.vw)vw=max(vw,vwfac*gpar(4,ng))
          dv=(vdmax-vdmin+2.*vw)/(nv+1)
          v0=vdmin-vw
106       continue
       endif
       if(argument(1:2).eq.'-p')read(argument(3:),*)psi
       if(argument(1:2).eq.'-h')goto 103
       if(argument(1:2).eq.'-s')read(argument(3:),*)nvs
       if(argument(1:2).eq.'-r')read(argument(3:),*)npsi
       if(argument(1:2).eq.'-T')read(argument(3:),*)vt2set
       if(argument(1:2).eq.'-t')read(argument(3:),*)vt2set2
       if(argument(1:2).eq.'-d')read(argument(3:),*)denset2
       if(argument(1:2).eq.'-v')read(argument(3:),*)vh
       if(argument(1:2).eq.'-i')read(argument(3:),*)imultiscan
       if(argument(1:2).eq.'-m')read(argument(3:),*)nframes
       if(argument(1:2).eq.'-w')read(argument(3:),*)vsmax
       if(argument(1:2).eq.'-l')read(argument(3:),'(a)')label
       if(argument(1:2).eq.'-c')call pfset(-3)
       if(argument(1:2).eq.'-E')lplotfv=.not.lplotfv
       if(argument(1:2).eq.'-M')lmain=.not.lmain
       if(argument(1:2).eq.'-U')lequil=.not.lequil
       if(argument(1:2).eq.'-C')lstab=.not.lstab
       if(argument(1:2).eq.'-S')lsech4=.not.lsech4
       if(argument(1:3).eq.'-DE')ldebug=.true.
       nvs=min(nvsmax,nvs)
       npsi=min(npsimax,npsi)
    enddo
    if(ng.gt.0)then
       call fvinit
    elseif(ng.eq.0)then
       call readfvin
       write(*,*)'Read in',nvin,dv,v(1),v(nvin+1)
    endif

    return
103 continue
    write(*,*)' Switches:'
    write(*,'(a,i2,a)')' -g<n>,<vs>,<Tpp>,<Tpl> no iteration <add&
         & gaussian> [',ng,']'
    do i=1,ng
       write(*,'(a,4f7.3,a)') '    [',gpar(:,i),']'
    enddo
    write(*,*)'-M  toggle force plot [',lmain
    write(*,*)'-U  equil scan        [',lequil
    write(*,*)'-C  ion-ion contours  [',lstab
    write(*,*)'-E  extra fv plots    [',lplotfv
    write(*,*)'-S  toggle sech4      [',lsech4
    write(*,*)'-DEBUG debug          [',ldebug
    write(*,*)'-c continuous running without stopping for plot viewing'
    write(*,'(a,i4)')' -i<imu> multiscan  [',imultiscan
    write(*,'(a,i4)')' -m<nfr> frames     [',nframes
    write(*,'(a,i4)')' -s<nvs>  read nvs  [',nvs
    write(*,'(a,i4)')' -r<npsi> read npsi [',npsi
    write(*,'(a,f6.2)')' -p<psi>  read psi  [',psi
    write(*,'(a,f6.2)')' -v<vh>   read vh   [',vh
    write(*,'(a,f6.2,a)')' -T<T1>   read T1   [',vt2set,' +v beam width'
    write(*,'(a,f6.2,a)')' -t<T2>   read T2   [',vt2set2,' -v beam width'
    write(*,'(a,f6.2)')' -d<n2>   read n2   [',denset2
    write(*,'(a,f6.2)')' -w<vsm>  read vsmax[',vsmax
    write(*,'(a,a)')' -l<char> read label[',label
    call exit
  end subroutine parseargs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  subroutine nofx(psithresh,ioff,ncase,vscase,iptype)
    real vscase(ncase)
    character*10 str2
    write(*,'(''For psi='',f6.3,'' Vshift cases:'',6f7.3)')psithresh,vscase
!       ldebug=.true.
    call ForceOfX(psithresh,vhin,ioff)
    call dcharsize(.02)
    call multiframe(2,1,3)
    if(iptype.eq.2)then
       call pltinit(x(1),x(nvp),0.,psithresh)
       call polyline(x(1+ioff),phia,nvp-ioff)
!    call autoplot(x(1+ioff),phia,nvp-ioff)
       call axis
       call axis2
       call axlabels('','!Af!@')
       call legendline(.9,.9,258,'(a)')
    endif
    do k=1,ncase
       vsthresh=vscase(k)
       gpar(2,1)=vsthresh
       gpar(2,2)=-vsthresh
       write(*,*)'gpar',gpar(:,1)
       call fvinit
       vhin=0.
!       ldebug=.true.
       call ForceOfX(psithresh,vhin,ioff)
       ldebug=.false.
       if(abs(dena(1)-1.).gt.0.001)then
          write(*,'(a,/,8f8.4)')'nofx call of ForceOfX gave dena not 1&
               & at edge',psithresh,vsthresh,vhin,dena(1),dena(nvp)
       endif
!       write(*,'(10f8.4)')'dena=',dena
       write(*,*)'f0,force,vhin=',f0,force,vhin
       if(k.eq.1)then
          call minmax(dena,nvp,dmin,dmax)
          dd=dmax-dmin
          call pltinit(x(1),x(nvp),-1.5*dd,1.5*dd)
          call axis; call axis2
!        call autoplot(x,dena-1.,nvp)
          call polyline(x,dena-1.,nvp)
          if(iptype.eq.2)then
             call legendline(.9,.9,258,'(b)') 
             call axlabels('','n-1')
          endif
          if(iptype.eq.3)then
             call legendline(.9,.9,258,'(c)')
             call axlabels('x','n-1')
          endif
          call fwrite(psithresh,iwidth,3,str2)
          call legendline(.1,.9,258,'!Ay!@='//str2(1:iwidth))
          call winset(.true.)
       else 
          call color(k-1)
          call polyline(x,dena-1.,nvp)
       endif
    enddo
    call pltend
    call multiframe(0,0,0)
  end subroutine nofx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
end module idenofphi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Ion-ion stability calculation
subroutine ionionstab
  use idenofphi
  real, dimension(ngmax) :: fstat
  real, dimension(nv) :: cvimag,cvreal
  real :: cvr0=0.
! Calculate the chi along just the real axis for thresholds
  call stationary(nv,f,ngmax,fstat,nstat)
  if(nstat.eq.3)then
! Find fractional depth assuming that the middle one is the local minimum.
     fmax=max(fstat(1),fstat(3))
     fmid=min(fstat(1),fstat(3))
     depth=(fmid-fstat(2))/fmid
!     write(*,*)nstat,(fstat(i),i=1,nstat),depth
  else
     write(*,*)'Incorrect number of stationaries',nstat,fstat
     call autoplot(v,f,nv)
     call pltend
  endif
! Here we are replacing this call. Not need for npts != nv
!  call chiion(v0,dv,f,nv,   fv,vp,npts,fmax                       &
!       &     ,nimax,0,0,vpi,vpimax,smfac,cvreal,cvimag)
  npts=nv
  do j=1,npts
     vp=v0+dv*j
     vi=1.+(j-1.)
     i=int(vi)
     xi=vi-i
     if(i.gt.nv)write(*,*)'Too large i',i,j
! Interpolated values of f(v) and f'(v)
!     fv(j)=(f(i)*(1.-xi)+f(i+1)*xi)
!     if(fv(j).gt.fmax)fmax=fv(j)
! Derivative indices
     ixd=nint(xi)
     id=i+ixd
     xid=xi+0.5-ixd
!--------------------
! Real values of vp.
! For a range of vp's, calculate the pvintegral and imaginary part.
     cvimag(j)=-3.1415926*(                                        &
          &        (f(id)-f(id-1))/(dv)*(1.-xid)                              &
          &        +(f(id+1)-f(id))/(dv)*xid)
! The ion principal value integral
     cvreal(j)=-pvint(v0,dv,f,nv,vp)
  enddo
!  call autoplot(v(1),cvimag,nv-1)
!  call polyline(v(1),cvreal,nv-1)
!  call polyline(v(1),f(1),nv-1)
!  call pltend
  do i=nv/10,nv-1  ! Find stability threshold.
     if(cvimag(i)*cvimag(i+1).le.0.and.cvreal(i).le.0)then
        ri=cvimag(i)/(cvimag(i)-cvimag(i+1))
        cvr0=ri*cvreal(i+1)+(1-ri)*cvreal(i)
        cvi0=ri*cvimag(i+1)+(1-ri)*cvimag(i)
        vp0=ri*(v0+dv*(i+1))+(1-ri)*(v0+dv*i)
!               write(*,*)i,ri,vp0,cvr0
        write(*,'(a,f8.4,a,f8.4,a,f8.4,a,f8.4,a,f8.4)')&
             &'shift',abs(gpar(2,1)-gpar(2,2))/2.,' n2',gpar(1,2),'&
             & T2/T1' ,gpar(3,2)/gpar(3,1), ' -T1/cvr0', -gpar(3,1)&
             &/cvr0,' depth',depth
        goto 1
     endif
  enddo
1 continue
  Tthresh=-gpar(3,1)/cvr0
end subroutine ionionstab
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine stationary(nv,f,ngmax,fstat,nstat)
  real f(nv),fstat(ngmax)
! Finds stationary values of f and return nstat values in fstat
  nstat=0
  do i=2,nv-1
     if((f(i+1)-f(i))*(f(i)-f(i-1)).lt.0..or.f(i+1)-f(i).eq.0)then
        nstat=nstat+1
        fstat(nstat)=f(i)
     endif
     if(nstat.eq.ngmax)return
  enddo
end subroutine stationary
!*********************************************************************
!      real function pvint(v0,dv,f,nv,vp)
! Integrate a distribution function f to obtain the susceptibility
! principal value integral defined as I(v_p)= \int {df(v)\over
! dv}/(v-v_p) dv. Here v_p is the phase velocity omega/k. f is the 1-d
! distribution function along v, and the resultant must be multiplied by
! (\omega_{pi}/k)^2/n_i to yield susceptibility.
! 
! The principal value integral is implemented by integrating by parts
! twice to obtain the form I = \int (v-v_p)(ln|v-v_p|-1)
! f'''(v)dv. (Which requires that f', f'' are zero at the extrema of the
! integration.  Then the value of f''' is evaluated by natural finite
! differences on a uniform mesh:
!
!  f'''_{i+1/2} = ([f_{i+2}+f_i-2f_{i+1}]- [f_{i+1}+f_{i-1}-2f_i])/dv^3
! 
! This is then integrated with I = \sum dv x value at i+1/2.  The
! limitations arise from the high order difference and rounding.  But
! the velocity function is non-singular once integrated, so no big
! difficulty arises in evaluating the integral through the former pole.
! f(v) ought to be zero for at least two values at the ends of its
! range.  The integral is done only over the non-zero values. The
! velocity grid must be uniform, and only v_0,dv is passed. 2:nv is the
! non-zero range. Convergence is second order in dv, i.e. in 1/nv.
!
! If velocities are in units of sqrt(2T_i/m_i), and the total i-th species
! density is given by n_i = \int f dv, then the i-susceptibility is
! \chi_i(\omega,k)= [Z_i^2 T_e/T_i] [-I(v_p)/2]/(n_e k^2\lambda_{De}^2)
!
! The imaginary part of the susceptibility in these units comes from the
! above equation with I(v_p) replaced by \pi df/dv giving half the residue.
! 
! Alternatively, if velocities are in units of sqrt(T_e/m_i), then
! \chi_i(\omega,k)= [Z_i^2] [-I(v_p)]/(n_e k^2\lambda_{De}^2)
!
! In either case, if f_i is normalized, one must multiply by n_i.
!
! Return pvint(=I) for the distribution function f, which is an array of
! values on a velocity grid starting at v0 (at index 0), spaced by dv,
! and extending to index nv+2. The first two and last two values would
! normally be zero, so as to satisfy boundary conditions.
! 
      real function pvint(v0,dv,f,nv,vp)
      integer nv
      real v0,vp,dv
      real f(0:nv+2)
      parameter (small=1.e-20)
      pvint=0.
      do i=1,nv
         v=v0+(i+0.5)*dv
         vd=v-vp
         if(abs(vd).lt.small)then 
            vfn=0.
         else
            vfn=vd*(alog(abs(vd))-1)
         endif
         pvint=pvint+vfn*((f(i+2)-f(i-1))-3.*(f(i+1)-f(i)))
      enddo
      pvint=pvint/dv**2
      end
!*********************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine testreadfvin
  use idenofphi
  if(ng.le.0)then
     call readfvin
     write(*,*)'Read in',nvin,dv,v(1),v(nvin+1)
  endif
  if(.true.)then
     call autoplot(v(0),f(0),nvin+3)
     call pltend
     call autoplot(vprime,fprime,nvin+2)
     call pltend
  endif
  phi=psi
  call denprimeint(0.,.false.)
  dp=den
  call denint(0.,.false.)
  write(*,*)'For vp=0 denprimeint gives',dp,' denint gives',den
  vpmax=sqrt(2.*phi)
  vpfac=.9
  call denprimeint(vpfac*vpmax,.true.)
  call denint(vpfac*vpmax,.true.)
  do i=1,nvp
     vpa(i)=-vpmax+2.*vpmax*(i-1)/(nvp-1.)
     call denprimeint(vpa(i),.false.)
     deldena(i)=delden
     dena(i)=den
!     write(*,*)'calling denint',i,nvp,vpa(i)
     call denint(vpa(i),.false.)
     dena2(i)=den
  enddo
  write(*,*)'dena2(1),dena2(nvp)',dena2(1),dena2(nvp)
  if(vh.eq.0)then
  write(*,*)'exp(-vpa(nvp)**2/2.)-1.=',exp(-vpa(nvp)**2/2.)-1.,' Error='&
       ,(deldena(nvp)-exp(-vpa(nvp)**2/2.)+1.)/(exp(-vpa(nvp)**2/2.)-1.)
     call charsize(.02,.02)
     call autoplot(vpa**2,dena-1.,nvp)
     call axlabels('vp**2','!Ad!@n')
     call color(1)
     call dashset(1)
!     call polyline(vpa**2,exp(-vpa**2/2.)-1.,nvp)
     call dashset(2)
     call color(4)
!     call polyline(vpa**2,deldena,nvp)
!     call dashset(4)
!     call polyline(vpa**2,dena2-1.,nvp)
     call dashset(0)
     call pltend
  else
     call autoplot(vpa,deldena,nvp)
     call axlabels('vp','!Ad!@n')
     call polyline([vpa(1),vpa(nvp)],[0.,0.],2)
     call pltend
  endif
  call autoplot(vpa,dena,nvp)
  call axis2
  call axlabels('vp','n')
  call color(5)
  call polyline(vpa,dena2,nvp)
  call pltend
  call autoplot(psi-vpa**2/2.,dena,nvp)
  call axlabels('!Ay!@-v!dp!d!u2!u/2','n')
  call pltend
  stop
end subroutine testreadfvin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine mainidenofphi
use idenofphi
character(20) :: string
integer :: jthresh=0
if(ng.ge.0)then
   if(ng.eq.0)call readfvin
   if(ng.gt.0)call fvinit
   call ForceOfX(.1,0.,1)
   write(*,'(a,1f6.3,a,g11.4e1,a,g11.4e1,a)')'psi=',psi  &
        ,' Force=',force,' dForce/dx=',fperx,' (-ve stable)'
   call charsize(.02,.02)
   call autoplot(v(0),f(0),nvin+3)
   call pltend
   call autoplot(vprime,fprime,nvin+2)
   call pltend
elseif(ng.eq.-1)then
   if(nvs.eq.0)nvs=nvsmax
   if(npsi.eq.0)npsi=npsimax
   ioff=1
   vsmax=2.
   psimax=psi
   v0=-vsmax-4.*sqrt(vt2)                    
   dv=2.*abs(v0)/(nv+3)
   write(*,'(a,f7.3,a,e9.3,a,f6.3,a,f6.3)')'Iteration over v0=',v0,'&
        & dv=',dv ,' vsmax=',vsmax,' vt2=',vt2
   do j=1,nvs            ! Cycle over Gaussian velocity shifts
      vs=vsmax*j/nvs
      vsa(j)=vs
      ng=2
      gpar(:,1)=[.5, vs,vt2,vt2]
      gpar(:,2)=[.5,-vs,vt2,vt2]
      if(lplotfv)then
         if(j.gt.3)then
            do k=1,npsi,npsi-1
               if(fperxa(j-2,k).gt.0..and.fperxa(j-1 ,k).le.0)then
!Mark the threshold
!            write(*,*)'jthresh+1,j',jthresh+1,k
                  call color(15-(15-4)*min(k-1,1))
                  call polymark(v,f,nvin+3,3)
                  call fwrite(psia(k),iwidth,2,string)
                  if(k.eq.1)call legendline(.03,.95,258,'Threshold case')
                  call legendline(-0.05,0.87-0.07*k/npsi,3, &
                       ' !Ay!@='//string(1:iwidth))
               endif
            enddo
         endif
      endif
      call fvinit
      if(lplotfv)then
         if(j.eq.1)then
            call charsize(.02,.02)
            call autoplot(v(0),f(0),nvin+3)
            call axis2
            call axlabels('!Bv!@','!Bf!@(!Bv!@)')
         endif
         call color(mod(j,15)+1)
         call polyline(v(0),f(0),nvin+3)   ! Plot the f
         if(j.eq.nvs)call pltend
      endif
      do i=1,npsi
         psiin=psimax*i/npsi
         call ForceOfX(psiin,0.,ioff)
         psia(i)=psiin
         fperxa(j,i)=fperx/psia(i)**2
!        if(j.eq.1)write(*,'(a,2f6.3,a,g11.4e1,a,g11.4e1,a)')'psi,vs=',psiin,vs&
!              ,' Force=',force,' dForce/dx=',fperxa(j,i),' (-ve stable)'
      enddo
      if(j.gt.1)then
      do i=1,npsi,npsi-1
         if(fperxa(j,i).le.0..and.fperxa(j-1,i).gt.0)then
            fmax=0. ; fmin=0.
            do k=1,nv
               fmax=max(fmax,f(k))
               if(v(k).le.dv/2)fmin=f(k)
            enddo
            jthresh=j-1
            jtha(i)=jthresh
            write(*,'(a,f7.4,a,f7.3,i4,a,f6.3)')'Stability threshold&
                 & for psi=' ,psia(i),' : !Bv!@!db!d=',vs,j,' f-dip='&
                 &,fmin/fmax
         endif
      enddo
      endif
   enddo
   if(nvs.gt.1)then
      call charsize(.02,.02)
      call autoinit(vsa,fperxa(1,1),nvs)
      call axis; call axis2
      call axlabels('!Bv!@!db!d','!Ad!@F/!Ad!@x/!Ay!@!u2!u')
      call polyline([vsa(1),vsa(nvs)],[0.,0.],2)
      
      jan=nint(0.25*nvs)
      call jdrwstr(wx2nx(vsa(jan)),wy2ny(fperxa(1,1))-.04,'!Ay!@=',0.)
      call jdrwstr(wx2nx(vsa(jan)),wy2ny(0)-.02,'stable',0.)
      call jdrwstr(wx2nx(vsa(jan)),wy2ny(0)+.02,'unstable',0.)
      do i=1,npsi
         call fwrite(psia(i),iwidth,2,string)
         call color(mod(i-1,15)+1)
         call polyline(vsa,fperxa(1,i),nvs)
         if(i.eq.1)call jdrwstr(wx2nx(vsa(jan)),wy2ny(fperxa(jan,i)),&
              string(1:lentrim(string)),1.)
         if(i.eq.npsi)call jdrwstr(wx2nx(vsa(jan)),wy2ny(fperxa(jan,i)),&
              string(1:lentrim(string)),-1.)
!   call labeline(vsa,fperxa(1,i),nvs,string,iwidth)
      enddo
      
      call pltend
      if(lplotfv)then
         do i=1,npsi,npsi-1
            jthresh=jtha(i)
            write(*,*)jthresh,psia(i)
            call nofx(psia(i),ioff,4,&
                 [(5+jthresh)*vsmax/nvs,(jthresh)*vsmax/nvs,&
                 &(-1.+jthresh)*vsmax/nvs,(-5.+jthresh)*vsmax/nvs],2)
            call nofx(psia(i),ioff,4,&
                 [(5+jthresh)*vsmax/nvs,(jthresh)*vsmax/nvs,&
                 &(-1.+jthresh)*vsmax/nvs,(-5.+jthresh)*vsmax/nvs],3)
         enddo
      endif
   endif
endif
end subroutine mainidenofphi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Given the gpar parameters, scan the hole velocity calculating steady
! force.
subroutine scanvhvs
use idenofphi
!use ISO_FORTRAN_ENV 
integer STDERR
!STDERR=ERROR_UNIT ! Use this if the following is broken.
STDERR=0
if(psi.eq.0)psi=.1
nvs=min(10,int(npsimax/2.))
vsmin=0.
vsmaxin=vsmax
denset1=1.-denset2
jthresh=0
vhcross=0.
fperx=0.
fperxm=0.
joff=0
v0=-(vsmax+3.0*sqrt(vt2set+vt2set2))           ! for fvinit
dv=2.*abs(v0)/(nv+3)
write(STDERR,*)
niterk=4    ! Number of iterative refinements default 4.
do k=1,niterk
   write(STDERR,'(a,3f8.5)')'Finding threshold',vsmin,vsmax,(vsmax-vsmin)/nvs
   jthresh=0
   do j=1,nvs
      vs=vsmin+(vsmax-vsmin)*j/nvs
      write(STDERR,'(i3,f7.4,$)')j,vs
      ng=2
      gpar(:,1)=[denset1, vs,vt2set,vt2set]
      gpar(:,2)=[denset2,-vs,vt2set2,vt2set2]
      call fvinit
      fva(:,j+joff)=f(1:nv+1)
      do i=1,nv+1
         vhin=v(i)
         fperxm=fperx
         call ForceOfX(psi,vhin,1)
!   write(*,'(a,1f6.3,a,g11.4e1,a,g11.4e1,a,g11.4e1)')'psi=',psi  &
!        ,' Force=',f0,' vhin=',vhin,' dForce/dx=',fperx
         forcea(i,j+joff)=f0/psi**2
         delfrcdelx(i,j+joff)=fperx/psi**2
         if(i.gt.2)then
            if(forcea(i,j+joff).le.0.and.forcea(i-1,j+joff).gt.1.e-5)then
               if(fperx.lt.0.and.jthresh.eq.0)then
                  write(STDERR,'(/,a,f7.3,a,f7.3,i4,a,2g10.3e1,a,2g10.3e1)')&
                       & 'vshift=',vs ,' vhole=',vhin,i,', F''s='&
                       &,forcea(i-1,j+joff) ,forcea(i,j+joff),', dFdx='&
                       &,fperx
                  jthresh=j
                  icross=i
                  vhcross=vhin
!                  write(STDERR,*)'joff,jthresh',joff,jthresh,fperxm,fperx
                  if(joff.gt.0.and.k.lt.niterk)exit ! Shortcuts
               endif
            endif
         endif
      enddo
     ! Always complete first vs scan, but then don't bother when joff>0.
      if(joff.gt.0.and.jthresh.gt.0)exit
   enddo
! Refine by setting vsmax=vthresh and vsmin=vs(jthresh-1), and doing over.
   vmx=vsmax
   vsmax=vsmin+(vsmax-vsmin)*(jthresh+.001)/float(nvs)
   vsmin=vsmin+(vmx-vsmin)*(jthresh-1.001)/float(nvs)
   if(joff.eq.0)write(STDERR,*)
   joff=nvs
enddo
vs=vsmin+(vsmax-vsmin)*(jthresh)/float(nvs)
jthresh=jthresh+joff
vsmax=vsmaxin

if(lplotfv)then
   call ForceOfX(psi,(v(icross)+v(icross-1))/2.,1)
   call autoplot(x,dena-1.,nvp)
   call axis2
   call boxtitle('Density change for threshold case')
   call axlabels('x','n-1')
   call pltend
endif

if(jthresh.gt.1+nvs)then    !Put interpolated threshold case in nvs+1
   ffp=abs(fperxm)/(abs(fperxm)+abs(fperx))
   fva(:,nvs+1)=(ffp*fva(:,jthresh)+(1.-ffp)*fva(:,jthresh-1))
   forcea(:,nvs+1)=(ffp*forcea(:,jthresh)+(1.-ffp)*forcea(:,jthresh-1))
   delfrcdelx(:,nvs+1)=(ffp*delfrcdelx(:,jthresh)+(1.-ffp)&
        &*delfrcdelx(:,jthresh-1))
elseif(jthresh.gt.0)then   ! Put threshold case in nvs+1
   fva(:,nvs+1)=fva(:,jthresh)
   forcea(:,nvs+1)=forcea(:,jthresh)
   delfrcdelx(:,nvs+1)=delfrcdelx(:,jthresh)
endif
if(denset1.eq.denset2.and.vt2set.eq.vt2set2)vhcross=0.

write(*,'(a,f6.3,a,f6.3,a,f7.3,a,f7.3)')'Threshold for psi=',psi,'&
     & stable vh=',vhcross,' vb/sqrt(T1)=',gpar(2,1)/sqrt(gpar(4,1))
write(STDERR,'(a,8f8.4)')'gpars:',gpar(:,1),gpar(:,2)

! Do plot:
if(lscanplot)then
call ionionstab

call dcharsize(.022)
call multiframe(3,1,1)
!call manautoinit(v(1),fva(1,1),nv+1,3,v(1),v(nv+1),symin,symax)
call minmax(fva(1,1),nv,fmin,fmax)
call pltinit(v(1),v(nv+1),fmin,1.05*fmax)
call axis; call axis2
call axlabels('','!Bf!@(!Bv!@)')
call legendline(.48,1.07,258,'!Bv!@')   ! Removed for PRL
do j=1,nvs
   call color(mod(j,14)+1)
   call polyline(v(1),fva(1,j),nv+1)
enddo
call color(15)
call winset(.true.)
call polyline([vhcross,vhcross],[0.,20.],2)
if(jthresh.gt.0)call polymark(v(1),fva(1,nvs+1),nv+1,3)
call legendline(.94,.9,258,'(i)')

call minmax(forcea(1,1),nv,fmin,fmax)
call pltinit(v(1),v(nv+1),1.1*fmin,1.1*fmax)
call axis; call axis2
call axlabels('','F/!Ay!@!u2!u')
!call axlabels('','F/!AM!@!d0!d!u2!u')         !PRL
do j=1,nvs
   call color(mod(j,14)+1)
   call polyline(v(1),forcea(1,j),nv+1)
enddo
call color(15)
call winset(.true.)
call polyline([vhcross,vhcross],[-30.,30.],2)
if(jthresh.gt.0)call polymark(v(1),forcea(1,nvs+1),nv+1,3)
call polyline([v(1),v(nv+1)],[0.,0.],2)
call legendline(.93,.9,258,'(ii)')

call minmax(delfrcdelx(1,1),nv,fmin,fmax)
call pltinit(v(1),v(nv+1),1.3*fmin,1.05*fmax)
call axis; call axis2
call axlabels('!Bv!@!dh!d'//label(1:lentrim(label)),'!Ad!@F/!Ad!@x/!Ay!@!u2!u')
!call axlabels('!Bv!@!dh!d'//label(1:lentrim(label)),'!Ad!@F/!Ad!@x/!AM!@!d0!d!u2!u')!PRL
do j=1,nvs
   call color(mod(j,14)+1)
   call polyline(v(1),delfrcdelx(1,j),nv+1)
enddo
call color(15)
call winset(.true.)
call polyline([vhcross,vhcross],[-30.,30.],2)
if(jthresh.gt.0)call polymark(v(1),delfrcdelx(1,nvs+1),nv+1,3)
call polyline([v(1),v(nv+1)],[0.,0.],2)
call legendline(.92,.9,258,'(iii)')
call multiframe(0,0,0)
call accisflush
endif

end subroutine scanvhvs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine multiscanvhvs(nscan)
use idenofphi
integer, parameter :: nT2=20,nn2=20,ncmax=20
real, dimension(nT2) :: T2a
real, dimension(nn2) :: n2a
real :: T2max=1.,n2max=0.5
character*20 string
!call pfset(-3) ! If you don't want to plot
lscanplot=.false.
call multiframe(nframes,nframes,0)
call dcharsize(.024,.024)
do k=1,nframes
   T2a(k)=T2max*(0.2+0.8/max(nframes-1,1)*(k-1))
   do j=1,nframes
      write(*,*)k,j
      n2a(j)=n2max*(0.2+0.8/max(nframes-1,1)*(nframes-j))
      denset2=n2a(j)
      denset1=1.-denset2
      vt2set=1.
      vt2set2=T2a(k)
do i=1,nscan
   psi=i*psimax/nscan
   if(nscan.eq.2)psi=psimax/10.**(i-1.)
   call scanvhvs
   if(psi.ne.i*psimax/nscan)write(*,*)'psi CHANGED. Do not use -p switch'
   fvthresh(:,i)=fva(:,nvs+1)
   icrsa(i)=icross
enddo
!call pltinit(v(1),v(nv),0.,.49)
call pltinit(-5.,5.,0.,.49)
call axis; call axis2
if(k.eq.1)call axlabels('','!Bf!@(!Bv!@)')
if(j.eq.nframes)call axlabels('!Bv!@','')
if(j.eq.nframes)call legendline(.05,.9,258,'+!Bv!@!dh!d !Ay!@=')
!call charsize(.017,.017)
hs=.08
do i=1,nscan
   psi=i*psimax/nscan
   if(nscan.eq.2)psi=psimax/10.**(i-1.)
   call fwrite(psi,iwidth,2,string)
   call color(i)
   call winset(.true.)
   if(j.eq.nframes)call legendline(hs,.9-i*hs,0,' '//string(1:iwidth))
   call polyline(v(1),fvthresh(1,i),nv+1)
   call polymark(v(icrsa(i)),fvthresh(icrsa(i),i),1,10)
   call color(15)
enddo
call fwrite(T2a(k),iwidth,2,string)
call legendline(hs,.9-i*hs,258,'T!d2!d='//string(1:iwidth))
call fwrite(n2a(j),iwidth,2,string)
i=i+1
call legendline(hs,.9-i*hs,258,'n!d2!d='//string(1:iwidth))
enddo
enddo
call pltend
lscanplot=.true.
end subroutine multiscanvhvs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ThreshCont
use idenofphi
integer, parameter :: nT2=20,nn2=20,ncmax=20
real, dimension(nT2,nn2) :: Tthresha,cworka
real, dimension(nT2) :: T2a
real, dimension(nn2) :: n2a
real, dimension(ncmax) :: zclv
real :: T2max=1.,n2max=0.5
!character*10 string
if(psi.lt.0.3)then
   thmax=20.
else
   thmax=10.
endif
do i=1,nT2
   T2a(i)=T2max*i/nT2
   do j=1,nn2
      n2a(j)=n2max*j/nn2
      denset2=n2a(j)
      denset1=1.-denset2
      vt2set=1.
      vt2set2=T2a(i)
      if(.false.)then
         Tthresh=i+.2*j
      else
         call scanvhvs
         call prtend(' ')  ! We don't pause for viewing here.
      endif
!      write(*,*)Tthresh
      Tthresha(i,j)=min(Tthresh,thmax)
   enddo
enddo
!call autocontour(Tthresha,nT2,nT2,nn2)
icl=0
zclv(1)=10.
icsw=1
call charsize(.02,.02)
call pltinit(T2a(1),T2a(nT2),n2a(1),n2a(nn2))
call contourl(Tthresha,cworka,nT2,nT2,nn2,zclv,icl,T2a,n2a,icsw)
call axis
call axis2
call axlabels('T!d2!d','n!d2!d')
call pltend

end subroutine ThreshCont
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use idenofphi
call parseargs
if(ldebug)call testreadfvin   
if(lmain)call mainidenofphi
if(lequil)then; call scanvhvs; call pltend
endif
if(imultiscan.gt.0)call multiscanvhvs(imultiscan)
if(lstab)call ThreshCont
end
