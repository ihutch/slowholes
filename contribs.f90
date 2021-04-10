integer, parameter :: npts=100
real, dimension(npts) :: v, f

vrange=8.
dv=vrange/(npts-1.)
voff=+1.
v0=-vrange/2.
vpsip=2.2
vpsim=-2.2
vphi=.7

do i=1,npts
   v(i)=v0+(i-1.)*dv
   f(i)=exp(-(v(i)-voff)**2/2.)
enddo
call pfset(3)
call dcharsize(.02)
call multiframe(2,1,2)
call autoplot(v,f,npts)
call axis2
call axlabels('!Bv!@!d!A;!@!d','!Bf!@!d!A;!@!d')
call polyline([0.,0.],[0.,1.],2)
call dashset(2)
call polyline([vpsim,vpsim],[0.,1.],2)
call jdrwstr(wx2nx(vpsim),wy2ny(1.1),'!As!@!d!A;!@!d!Bv!@!d!Ay!@!d',0.)
call polyline([vpsip,vpsip],[0.,1.],2)
call jdrwstr(wx2nx(vpsip),wy2ny(1.1),'!As!@!d!A;!@!d!Bv!@!d!Ay!@!d',0.)
call dashset(3)
call color(3)
call polyline([vphi,vphi],[0.,1.],2)
call jdrwstr(wx2nx(vphi),wy2ny(1.1),'!Bv!@!d!Af!@!d',0.)
call color(4)
call polyline([-vphi,-vphi],[0.,1.],2)
call jdrwstr(wx2nx(-vphi),wy2ny(1.1),'!Bv!@!d!Af!@!d',0.)
call dashset(5)
call color(3)
call polyline([vphi,vrange/2],[.25,.25],2)
call polyline([vpsim,-vrange/2],[.25,.25],2)
call polyline([vphi,vpsip],[.2,.2],2)
call jdrwstr(wx2nx((vphi+vpsip)/2.),wy2ny(.15),'Reflected',0.)
call jdrwstr(wx2nx(vpsip),wy2ny(.3),'Incoming',0.)
call jdrwstr(wx2nx((vpsim-vrange/2.)/2.),wy2ny(.3),'Transmitted',0.)
call jdrwstr(wx2nx((vpsim+vphi)/2.),wy2ny(.3),'No Contribution',0.)
call jdrwstr(wx2nx((vpsim-vrange/2.)/2.),wy2ny(.1),'!As!@!d!A;!@!d positive',0.)
call color(4)
call jdrwstr(wx2nx((vpsim-vrange/2.)/2.),wy2ny(.9),'!As!@!d!A;!@!d negative',0.)
call polyline([-vphi,-vrange/2],[.75,.75],2)
call polyline([vpsip,vrange/2],[.75,.75],2)
call polyline([-vphi,vpsim],[.8,.8],2)
call jdrwstr(wx2nx((-vphi+vpsim)/2.),wy2ny(.85),'Reflected',0.)
call jdrwstr(wx2nx(vpsim),wy2ny(.7),'Incoming',0.)
call jdrwstr(wx2nx((vpsip+vrange/2.)/2.),wy2ny(.7),'Transmitted',0.)
call jdrwstr(wx2nx((vpsip-vphi)/2.),wy2ny(.7),'No Contribution',0.)

call pltend
end
