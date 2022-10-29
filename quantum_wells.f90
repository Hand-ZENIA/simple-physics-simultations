program wells
implicit none
integer, parameter:: N=2000,Nw=11,LWORK=3*N-1,NE=1000
double precision, parameter:: V0=300,W=1.5,d=0.5,dext=3*d,Vmax=1000,echelle=0.07357,eta=0.5
double precision::x(N),V(N),H(N,N),WORK(LWORK),E(N),eps(NE),dos(NE)
integer:: i,j
double precision:: alpha,dx,de
character:: INFO
double precision, external:: lorentz

call makewells(N,Nw,W,d,dext,V0,Vmax,V,x)

open(unit=20,file='pot.dat')
do i=1,N
   write(20,*) x(i),"  ",v(i)
enddo
close(20)

dx=x(2)-x(1)
alpha=1/dx**2/echelle

H=0
do i=1,N
   H(i,i)=v(i)+2*alpha
   if(i<N) H(i,i+1)=-alpha
   if(i>1) H(i,i-1)=-alpha
enddo

call dsyev('V','U',N,H,N,E,WORK,LWORK,INFO)
!E=-E*echelle
write(*,*) E(1),E(2),E(3),E(4)
open(unit=21,file='wfs.dat')
do i=1,N
   write(21,*) x(i),"  ",H(i,1),"  ",H(i,2),"  ",H(i,3),"  ",H(i,4)
enddo
close(21)


! compute DOS
de=V0/(NE-1)
do i=1,NE
   eps(i)=-V0+(i-1)*de
enddo

open(unit=22,file='dos.dat')
do i=1,NE
   dos(i)=0
   do j=1,N
     dos(i)=dos(i)+lorentz(eps(i)-E(j),eta)
   enddo
   write(22,*) eps(i),"  ",dos(i)
enddo
close(22)

end program wells

subroutine makewells(N,Nw,w,d,dext,v0,vmax,v,x)
implicit none
integer, intent(in)::N,Nw
double precision, intent(in):: v0,w,d,dext,vmax
double precision, intent(out)::x(N),v(N)
!!! internal variables
integer :: i,j
double precision:: xmax,dx

xmax=Nw*w+(Nw-1)*d+2*dext
dx=xmax/(N-1)

do i=1,N
 x(i)=(i-1)*dx
 v(i)=0
 do j=0,Nw-1
    if((x(i)>(3*d+j*(w+d))) .and. (x(i)<(3*d+j*(w+d)+w))) v(i)=-v0
 enddo
enddo
v(1)=vmax
v(N)=v(1)

end subroutine makewells

double precision function lorentz(x,eta)
implicit none
double precision, intent(in)::x,eta
double precision, parameter:: pi=3.141592653589793
lorentz=eta**2/(x**2+eta**2)/pi
end function lorentz
