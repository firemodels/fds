program pipe_exact

implicit none
real h,k(3),r(4),T(0:5),q,pi,TT,rr
integer i,iz

pi = atan(1.)*4.
h = 10.
k(1) = 0.2
k(2) = 50.
k(3) = 0.2
r(1) = 0.01
r(2) = r(1)+0.02
r(3) = r(2)+0.01
r(4) = r(3)+0.02
T(0) = 20.
T(5) = 480.
q = (T(5)-T(0))/( (1./(2*pi*r(1)*h)) + log(r(2)/r(1))/(2*pi*k(1)) + log(r(3)/r(2))/(2*pi*k(2)) + log(r(4)/r(3))/(2*pi*k(3)) +&
                  (1./(2*pi*r(4)*h)) )
T(1) = T(0) + q/(2*pi*r(1)*h)
T(2) = T(1) + log(r(2)/r(1))*q/(2*pi*k(1))
T(3) = T(2) + log(r(3)/r(2))*q/(2*pi*k(2))
T(4) = T(3) + log(r(4)/r(3))*q/(2*pi*k(3))

write(10,'(a)') 'Depth,Temp'  
do i=100,0,-1
rr = r(1) + i*(r(4)-r(1))/100.
if (rr<r(2)) iz=2
if (r(2)<=rr .and. rr<r(3)) iz=3
if (r(3)<=rr .and. rr<=r(4)) iz=4
TT = T(iz) + log(rr/r(iz))*(T(iz-1)-T(iz))/log(r(iz-1)/r(iz))
write(10,'(f7.4,a,f7.2)') r(4)-rr,',',TT
enddo

end program
