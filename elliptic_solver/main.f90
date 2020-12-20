program main
    use precision
    use input_functions
    use elliptic_equation_solver
implicit none
integer :: N, M
real(mp) :: lx, ly !endpoints of x and y intervals, respectively
real(mp) :: hx, hy
real(mp) :: c1, c2, d1, d2
real(mp) :: eps, rho
real(mp), allocatable :: x(:), y(:), residuals(:,:)
real(mp), allocatable :: u_init_guess(:,:), u(:,:)
integer :: i,j ! counter
integer :: file_id
N = 20
M = 20
lx = 1.0_mp
ly = 1.0_mp
c1 = 1.0_mp
c2 = 1.0_mp
d1 = 1.0_mp
d2 = 1.0_mp
eps = 1.0d-3
allocate(x(0:N))
allocate(y(0:M))
allocate(u_init_guess(0:N,0:M), u(0:N,0:M))
call coord_grid(N,M,lx,ly,x,y,hx,hy)

u_init_guess = 0.0_mp

!open(unit=1,file='C:\Users\Marta\source\repos\elliptic_solver\elliptic_solver\iterative_method_5.txt',status='replace')
!open(unit=2,file='C:\Users\Marta\source\repos\elliptic_solver\elliptic_solver\iterative_method_10.txt',status='replace')
open(unit=3,file='C:\Users\Marta\source\repos\elliptic_solver\elliptic_solver\iterative_method_20.txt',status='replace')
file_id = 3
call iterative_method(p, q, f, mu, u_init_guess, x, y, hx, hy, lx, ly, c1, c2, d1, d2, eps, u, file_id)

write(file_id,*)'U(x,y):'
write(file_id,'("y \ x     ",6("    &    ",f10.3))') x(0), x(1), x(2), x(3), x(4), x(5)
do j=0,M
    write(file_id,'(f10.3,6("    &    ",f10.5))') y(j), u(0,j), u(1,j), u(2,j), u(3,j), u(4,j), u(5,j)
enddo
write(file_id,*)

write(file_id,*)'U*(x,y):'
write(file_id,'("y \ x     ",6("    &    ",f10.3))') x(0), x(1), x(2), x(3), x(4), x(5)
do j=0,M
    write(file_id,'(f10.3,6("    &    ",f10.5))') y(j), reference_solution(x(0),y(j)), reference_solution(x(1),y(j)), reference_solution(x(2),y(j)), reference_solution(x(3),y(j)), reference_solution(x(4),y(j)), reference_solution(x(5),y(j)) !reference_solution(x(0),y(j))
enddo
close(file_id)
deallocate(x,y,u_init_guess,u)   
end program