module constantes
implicit none
save
real*8, parameter :: pi=3.14159274 ! pi
real*8, parameter :: mu=3.98600E5 ! GMearth (km^3/s^2)
real*8, parameter :: Re=6378.14 ! Earth radius (km)
real*8, parameter :: deg2rad = pi/180 
end module constantes

program main
use constantes
implicit none
real*8 :: sma,ecc,inc,Omega,om,Ma,epoch
real*8 :: E,X,Y,Z,x1,y1,z1,x2,y2,z2,x3,y3,z3
write(*,*)pi
call readTLE(sma,ecc,inc,Omega,om,Ma,epoch)
call eccano(Ma,ecc,E)
call XYZ (sma,E,ecc,X,Y,Z)
call rotation1 (om,X,Y,Z,x1,y1,z1)
call rotation2 (inc,x1,y1,z1,x2,y2,z2)
call rotation3 (Omega,x2,y2,z2,x3,y3,z3)
end program main

subroutine readTLE(sma,ecc,inc,Omega,om,Ma,epoch)
use constantes
implicit none
integer ::y,e1

real*8 ::d,sma,ecc,inc,Omega,om,epoch,AoP,Ma,Mm,Mm1
read(*,"(18x,i2,f12.8)")y,d
if (y==20) then
epoch=d
endif
if (y==21) then
epoch=d+366.0
endif
print*, "The epoch is :",epoch
read(*,"(9x,f7.4,x,f8.4,4x,i4,x,f8.4,x,f8.4,x,f12.9)")inc,Omega,e1,AoP,Ma,Mm1
write(*,"(a25,f7.4)")"Inclination (degrees) : ",inc
inc = inc*deg2rad
print*, "Inclination (radian) : ",inc
write(*,"(a47,f8.4)")"Longitude of the ascending node Ω (degrees) : ",Omega
Omega = Omega*deg2rad
print*,"Longitude of the ascending node Ω (radian) : ",Omega 
write(*,"(a6,i4)")" e1 = ",e1
ecc=e1*1E-7
write(*,"(a40,f9.7)") "Eccentricity (decimal point assumed) : ",ecc
write(*,"(a22,f8.4)")"Argument of Perigee :",AoP
AoP = AoP*deg2rad
print*,"Argument of Perigee (in radian) :",AoP
write(*,"(a17,f8.4)")"Mean anomaly M : ",Ma
Ma = Ma*deg2rad
print*,"Mean anomaly M (radian) : ",Ma
print*,"The revolution per day is :",Mm1
Mm=((Mm1*2*pi)/86400)
write(*,"(a20,f18.14)")"Mean motion Mm is :",Mm
sma= (mu**(1.0/3.0))/(Mm**(2.0/3.0))
print*, "The semi major axis is:", sma

end subroutine readTLE

subroutine eccano(Ma,ecc,E)
use constantes
implicit none
real*8, intent(in) :: Ma,ecc
real*8, intent(out) :: E
integer :: n
n = 0
E = Ma
print*,"E",n,"=",E
n = 1 
do while(n<=5)
E = E-(E-ecc*sin(E)-Ma)/(1-ecc*cos(E))
print*,"E",n,"=",E
n = n+1
enddo
end subroutine eccano

subroutine XYZ (sma,E,ecc,X,Y,Z)
implicit none
real*8, intent(out) :: X,Y,Z
real*8, intent(in):: sma,E,ecc
X = sma*(cos(E)-ecc)
Y = sma*sqrt(1-ecc**2)*sin(E)
Z = 0
print*,"the position X,Y,Z are : ", X,Y,Z
end subroutine XYZ

subroutine rotation1 (AoP,X,Y,Z,x1,y1,z1)
implicit none
real*8, intent(in) :: X,Y,Z,AoP
real*8, intent(out):: x1,y1,z1
!call XYZ (sma,E,ecc,X,Y,Z)
x1 = cos(AoP)*X - sin(AoP)*Y
y1 = sin(AoP)*X + cos(AoP)*Y
z1 = Z
print*,"x1,y1,z1 are equal to : ", x1,y1,z1
end subroutine rotation1

subroutine rotation2 (inc,x1,y1,z1,x2,y2,z2)
real*8, intent(out) :: x2,y2,z2
real*8, intent(in):: inc,x1,y1,z1
x2 = x1
y2 = cos(inc)*y1 - sin(inc)*z1
z2 = sin(inc)*y1 + cos(inc)*z1
print*,"x2,y2,z2 are equal to : ", x2,y2,z2
end subroutine rotation2

subroutine rotation3 (omega,x2,y2,z2,x3,y3,z3)
real*8, intent(out) :: x3,y3,z3
real*8, intent(in):: omega,x2,y2,z2
x3 = cos(omega)*x2 - sin(omega)*y2
y3 = sin(omega)*x2 + cos(omega)*y2
z3 = z2
print*,"x3,y3,z3 are equal to : ", x3,y3,z3
end subroutine rotation3

print,"x3,y3,z3 are equal to : ", x3,y3,z3
end subroutine rotation3

subroutine geographic (x,w,y,z,r,α,φ′)
real8, intent(in) :: x,w,y,z
real8, intent(out) :: α,φ′
r = sqrt(x2 + y2 + z**2)
α = artan(y/x)
φ′ = arcsin(z/r)
print,"The geocentric distance is : ",r 
print, "the ",α
print, "the geocentric latitude :",φ′
end subroutine geographic