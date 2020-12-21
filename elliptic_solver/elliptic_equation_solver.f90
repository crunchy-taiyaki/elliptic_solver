module elliptic_equation_solver
    use precision
    use constants
    use input_functions, only : reference_solution
    implicit none
    contains
    
    subroutine coord_grid(N,M,lx,ly,x,y,hx,hy)
    integer, intent(in) :: N, M ! number of grid points for x and y axis
    real(mp), intent(in) :: lx, ly !endpoints of x and y intervals, respectively
    real(mp), intent(out) :: x(0:N), y(0:M)
    real(mp), intent(out) :: hx, hy ! size of grid step
    integer :: i,j
    hx = lx/N
    hy = ly/M
    forall(i=0:N) x(i) = i*hx
    forall(j=0:M) y(j) = j*hy    
    end subroutine coord_grid
    
    function eigenvalues_interval(lx,ly,hx,hy,c1,c2,d1,d2) result(delta)
    real(mp), intent(in) :: lx, ly !endpoints of x and y intervals, respectively
    real(mp), intent(in) :: hx, hy ! the grid steps
    real(mp), intent(in) :: c1, c2, d1, d2! strictly positive constants for which the conditions (c1.le.p(x,y)).and.(c2.ge.p(x,y)), (d1.le.q(x,y)).and.(d2.ge.q(x,y)) is satisfied
    real(mp) :: delta(1:2) !this is endpoints of interval which consist eigenvalues
    delta(1) = c1*4.0_mp*sin(pi*hx/2.0_mp/lx)**2/hx**2 + d1*4.0_mp*sin(pi*hy/2.0_mp/ly)**2/hy**2
    delta(2) = c2*4.0_mp*cos(pi*hx/2.0_mp/lx)**2/hx**2 + d2*4.0_mp*cos(pi*hy/2.0_mp/ly)**2/hy**2
    end function eigenvalues_interval
    
    function matrix_norm(A)
        real(mp), intent(in) :: A(0:,0:)
        real(mp) :: matrix_norm
        integer :: N, M, i, j
        N = size(A,1) - 1
        M = size(A,2) - 1
        matrix_norm = maxval(abs(A(1:N-1,1:M-1)))
    end function matrix_norm
    
    subroutine iterative_method(p, q, f, mu, u0, x, y, hx, hy, lx, ly, c1, c2, d1, d2, eps, u, file_id)
    interface    
    pure function p(x,y)
    real(8), intent(in) :: x,y
    real(8) :: p
    end function p
    
    pure function q(x,y)
    real(8), intent(in) :: x,y
    real(8) :: q
    end function q
    
    pure function mu(x,y)
    real(8), intent(in) :: x,y
    real(8) :: mu
    end function mu    
    
    pure function f(x,y)
    real(8), intent(in) :: x,y
    real(8) :: f
    end function f
    end interface
    
    real(mp), intent(in) :: x(0:), y(0:), u0(0:,0:)
    real(mp), intent(in) :: hx, hy ! the grid steps
    real(mp), intent(in) :: lx, ly !endpoints of x and y intervals, respectively
    real(mp), intent(in) :: c1, c2, d1, d2
    real(mp), intent(in) :: eps
    integer, intent(in) :: file_id
    real(mp), intent(out) :: u(0:,0:)
    real(mp) :: rho ! radius of method convergence
    integer  :: N, M
    integer :: max_iter ! maximum of iteration number
    integer :: k ! iteration counter
    integer :: i,j ! grid nodes counters
    real(mp) :: delta(1:2), xi
    real(mp), allocatable :: u_prev(:,:)
    real(mp), allocatable :: Lh_u(:,:), Lh_u0(:,:),f_ij(:,:), residual_k(:), abs_residual_k(:), u_ref(:,:), u_dif_k(:), dash_rho_k(:)!dummy
    real(mp), allocatable :: rho_dep_residual(:)!dummy
    real(mp) :: residual_u0, abs_residual_u0 !dummy

    delta = eigenvalues_interval(lx,ly,hx,hy,c1,c2,d1,d2)
    xi = delta(1)/delta(2)
    max_iter = int(log(1.0_mp/eps)/(2.0_mp*xi)) + 1
    rho = (delta(2)-delta(1))/(delta(2)+delta(1))
    N = size(x)-1
    M = size(y)-1
    allocate(u_prev(0:N,0:M))
    allocate(Lh_u(0:N,0:M),Lh_u0(0:N,0:M),f_ij(0:N,0:M),residual_k(0:max_iter),abs_residual_k(0:max_iter))!dummy
    allocate(u_ref(0:N,0:M),u_dif_k(0:max_iter),rho_dep_residual(0:max_iter),dash_rho_k(0:max_iter))!dummy
    u = u0
    
    forall (j=0:M) u(0,j) = mu(0.0_mp,y(j))
    forall (j=0:M) u(N,j) = mu(lx,y(j))
    forall (i=1:N-1) u(i,0) = mu(x(i),0.0_mp)
    forall (i=1:N-1) u(i,M) = mu(x(i),ly)
    k = 0
    u_dif_k(0) = 1.0_mp !dummy
    do while (k.lt.max_iter)
        k = k+1
        u_prev = u
        do i=1,N-1
            do j=1,M-1
                u(i,j) = (p(x(i)-hx/2.0_mp,y(j))*u_prev(i-1,j)/(hx**2) + p(x(i)+hx/2.0_mp,y(j))*u_prev(i+1,j)/(hx**2) +&
                        & q(x(i),y(j)-hy/2.0_mp)*u_prev(i,j-1)/(hy**2) + q(x(i),y(j)+hy/2.0_mp)*u_prev(i,j+1)/(hy**2) +&
                        & f(x(i),y(j)))  /  (p(x(i)-hx/2.0_mp,y(j))/(hx**2) + p(x(i)+hx/2.0_mp,y(j))/(hx**2) +&
                        & q(x(i),y(j)-hy/2.0_mp)/(hy**2) + q(x(i),y(j)+hy/2.0_mp)/(hy**2))
                
                Lh_u(i,j) = p(x(i)+hx/2.0_mp,y(j))*(u(i+1,j)-u(i,j))/hx**2 - &  !dummy
                    & p(x(i)-hx/2.0_mp,y(j))*(u(i,j)-u(i-1,j))/hx**2 + &
                    & q(x(i),y(j)+hy/2.0_mp)*(u(i,j+1)-u(i,j))/hy**2 - &
                    & q(x(i),y(j)-hy/2.0_mp)*(u(i,j)-u(i,j-1))/hy**2
                
                f_ij(i,j) = f(x(i),y(j)) !dummy
                u_ref(i,j) = reference_solution(x(i),y(j)) !dummy
            enddo
        enddo
        residual_k(k) = matrix_norm(Lh_u + f_ij) !dummy
        abs_residual_k(k) = matrix_norm(u - u_ref)!dummy
        u_dif_k(k) = matrix_norm(u - u_prev)!dummy
        rho_dep_residual(k) = rho*u_dif_k(k)/(1.0_mp-rho)!dummy
        if (k.lt.2) then !dummy
            dash_rho_k(k) = sqrt(-1.d0) !dummy
        else !dummy
            dash_rho_k(k) = sqrt(u_dif_k(k)/u_dif_k(k-2))!dummy
        endif !dummy
    enddo
    !dummy part!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i=1,N-1
        do j=1,M-1
            Lh_u0(i,j) = p(x(i)+hx/2.0_mp,y(j))*(u0(i+1,j)-u0(i,j))/hx**2 - &
                    & p(x(i)-hx/2.0_mp,y(j))*(u0(i,j)-u0(i-1,j))/hx**2 + &
                    & q(x(i),y(j)+hy/2.0_mp)*(u0(i,j+1)-u0(i,j))/hy**2 - &
                    & q(x(i),y(j)-hy/2.0_mp)*(u0(i,j)-u0(i,j-1))/hy**2
        enddo
    enddo

    residual_u0 = matrix_norm(Lh_u0 + f_ij)
    abs_residual_u0 = matrix_norm(u0 - u_ref)
    write(file_id,*)'1) residual u*:',residual_k(max_iter)
    write(file_id,*)'2) residual u0:',residual_u0
    write(file_id,*)'3) eps:',eps,'m:',max_iter
    write(file_id,*)'4) rho:', rho
    write(file_id,'("    k    &    ||F-A*Uk||    &     rel.d        &     ||Uk-u*||    &     rel.error    &  ||Uk - Uk-1||   &    apost.er.     &    dash_rho_k")')   
    do k=1, max_iter
        write(file_id,'(i5,7("    &    ",e10.3))') k, residual_k(k), residual_k(k)/residual_u0, abs_residual_k(k), abs_residual_k(k)/abs_residual_u0,u_dif_k(k),rho_dep_residual(k),dash_rho_k(k)
    enddo
    write(file_id,*)
    deallocate(Lh_u,Lh_u0,f_ij,residual_k,abs_residual_k)!dummy
    deallocate(u_ref,u_dif_k,rho_dep_residual,dash_rho_k)!dummy
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    deallocate(u_prev)
    
    end subroutine iterative_method
            
    subroutine gauss_seidel_method(p, q, f, mu, u0, x, y, hx, hy, lx, ly, c1, c2, d1, d2, eps, u, file_id)
    interface    
    pure function p(x,y)
    real(8), intent(in) :: x,y
    real(8) :: p
    end function p
    
    pure function q(x,y)
    real(8), intent(in) :: x,y
    real(8) :: q
    end function q
    
    pure function mu(x,y)
    real(8), intent(in) :: x,y
    real(8) :: mu
    end function mu    
    
    pure function f(x,y)
    real(8), intent(in) :: x,y
    real(8) :: f
    end function f
    end interface
    real(mp), intent(in) :: x(0:), y(0:), u0(0:,0:)
    real(mp), intent(in) :: hx, hy ! the grid steps
    real(mp), intent(in) :: lx, ly !endpoints of x and y intervals, respectively
    real(mp), intent(in) :: c1, c2, d1, d2
    real(mp), intent(in) :: eps
    integer, intent(in) :: file_id
    real(mp), intent(out) :: u(0:,0:)
    real(mp) :: rho ! radius of method convergence
    integer  :: N, M
    integer :: max_iter ! maximum of iteration number
    integer :: k ! iteration counter
    integer :: i,j ! grid nodes counters
    real(mp) :: delta(1:2), xi
    real(mp), allocatable :: u_prev(:,:)
    real(mp), allocatable :: Lh_u(:,:), Lh_u0(:,:),f_ij(:,:), residual_k(:), abs_residual_k(:), u_ref(:,:), u_dif_k(:), dash_rho_k(:)!dummy
    real(mp), allocatable :: rho_dep_residual(:)!dummy
    real(mp) :: residual_u0, abs_residual_u0 !dummy
    
    delta = eigenvalues_interval(lx,ly,hx,hy,c1,c2,d1,d2)
    xi = delta(1)/delta(2)
    max_iter = int(log(1.0_mp/eps)/(4.0_mp*xi)) + 1
    rho = (delta(2)-delta(1))/(delta(2)+delta(1))
    N = size(x)-1
    M = size(y)-1
    allocate(u_prev(0:N,0:M))
    allocate(Lh_u(0:N,0:M),Lh_u0(0:N,0:M),f_ij(0:N,0:M),residual_k(0:max_iter),abs_residual_k(0:max_iter))!dummy
    allocate(u_ref(0:N,0:M),u_dif_k(0:max_iter),rho_dep_residual(0:max_iter),dash_rho_k(0:max_iter))!dummy
    u = u0

    forall (j=0:M) u(0,j) = mu(0.0_mp,y(j))
    forall (j=0:M) u(N,j) = mu(lx,y(j))
    forall (i=1:N-1) u(i,0) = mu(x(i),0.0_mp)
    forall (i=1:N-1) u(i,M) = mu(x(i),ly)
    k = 0
    u_dif_k(0) = 1.0_mp
    do while (k.lt.max_iter)
        k = k+1
        u_prev = u
        do i=1,N-1
            do j=1,M-1
                u(i,j) = (p(x(i)-hx/2.0_mp,y(j))*u(i-1,j)/hx**2 + p(x(i)+hx/2.0_mp,y(j))*u_prev(i+1,j)/hx**2 +&
                        & q(x(i),y(j)-hy/2.0_mp)*u(i,j-1)/hy**2 + q(x(i),y(j)+hy/2.0_mp)*u_prev(i,j+1)/hy**2 +&
                        & f(x(i),y(j)))  /  (p(x(i)-hx/2.0_mp,y(j))/hx**2 + p(x(i)+hx/2.0_mp,y(j))/hx**2 +&
                        & q(x(i),y(j)-hy/2.0_mp)/hy**2 + q(x(i),y(j)+hy/2.0_mp)/hy**2)
                
                Lh_u(i,j) = p(x(i)+hx/2.0_mp,y(j))*(u(i+1,j)-u(i,j))/hx**2 - &  !dummy
                    & p(x(i)-hx/2.0_mp,y(j))*(u(i,j)-u(i-1,j))/hx**2 + &
                    & q(x(i),y(j)+hy/2.0_mp)*(u(i,j+1)-u(i,j))/hy**2 - &
                    & q(x(i),y(j)-hy/2.0_mp)*(u(i,j)-u(i,j-1))/hy**2
                
                f_ij(i,j) = f(x(i),y(j)) !dummy
                u_ref(i,j) = reference_solution(x(i),y(j)) !dummy
            enddo
        enddo
        residual_k(k) = matrix_norm(Lh_u + f_ij) !dummy
        abs_residual_k(k) = matrix_norm(u - u_ref)!dummy
        u_dif_k(k) = matrix_norm(u - u_prev)!dummy
        rho_dep_residual(k) = rho*u_dif_k(k)/(1.0_mp-rho)!dummy
        if (k.lt.2) then !dummy
            dash_rho_k(k) = sqrt(-1.d0) !dummy
        else !dummy
            dash_rho_k(k) = sqrt(u_dif_k(k)/u_dif_k(k-2))!dummy
        endif !dummy
    enddo
    !dummy part!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i=1,N-1
        do j=1,M-1
            Lh_u0(i,j) = p(x(i)+hx/2.0_mp,y(j))*(u0(i+1,j)-u0(i,j))/hx**2 - &
                    & p(x(i)-hx/2.0_mp,y(j))*(u0(i,j)-u0(i-1,j))/hx**2 + &
                    & q(x(i),y(j)+hy/2.0_mp)*(u0(i,j+1)-u0(i,j))/hy**2 - &
                    & q(x(i),y(j)-hy/2.0_mp)*(u0(i,j)-u0(i,j-1))/hy**2
        enddo
    enddo
    residual_u0 = matrix_norm(Lh_u0 + f_ij)
    abs_residual_u0 = matrix_norm(u0 - u_ref)
    write(file_id,*)'1) residual u*:',residual_k(max_iter)
    write(file_id,*)'2) residual u0:',residual_u0
    write(file_id,*)'3) eps:',eps,'m:',max_iter
    write(file_id,*)'4) rho:', rho
    write(file_id,'("    k    &    ||F-A*Uk||    &     rel.d        &     ||Uk-u*||    &     rel.error    &  ||Uk - Uk-1||   &    apost.er.     &    dash_rho_k")')   
    do k=1, max_iter
        write(file_id,'(i5,7("    &    ",e10.3))') k, residual_k(k), residual_k(k)/residual_u0, abs_residual_k(k), abs_residual_k(k)/abs_residual_u0,u_dif_k(k),rho_dep_residual(k),dash_rho_k(k)
    enddo
    write(file_id,*)
    deallocate(Lh_u,Lh_u0,f_ij,residual_k,abs_residual_k)!dummy
    deallocate(u_ref,u_dif_k,rho_dep_residual,dash_rho_k)!dummy
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    deallocate(u_prev)
    end subroutine gauss_seidel_method
            
            
    subroutine successive_over_relaxation_method(p, q, f, mu, u0, x, y, hx, hy, lx, ly, c1, c2, d1, d2, eps, u, file_id)
    interface    
    pure function p(x,y)
    real(8), intent(in) :: x,y
    real(8) :: p
    end function p
    
    pure function q(x,y)
    real(8), intent(in) :: x,y
    real(8) :: q
    end function q
    
    pure function mu(x,y)
    real(8), intent(in) :: x,y
    real(8) :: mu
    end function mu    
    
    pure function f(x,y)
    real(8), intent(in) :: x,y
    real(8) :: f
    end function f
    end interface
    real(mp), intent(in) :: x(0:), y(0:), u0(0:,0:)
    real(mp), intent(in) :: hx, hy ! the grid steps
    real(mp), intent(in) :: lx, ly !endpoints of x and y intervals, respectively
    real(mp), intent(in) :: c1, c2, d1, d2
    real(mp), intent(in) :: eps
    integer, intent(in) :: file_id
    real(mp), intent(out) :: u(0:,0:)
    real(mp) :: rho, S ! radius of method convergence and spectral radius, respectively
    integer  :: N, M
    integer :: max_iter ! maximum of iteration number
    integer :: k ! iteration counter
    integer :: i,j ! grid nodes counters
    real(mp) :: delta(1:2), xi, omega
    real(mp), allocatable :: u_prev(:,:)
    real(mp), allocatable :: Lh_u(:,:), Lh_u0(:,:),f_ij(:,:), residual_k(:), abs_residual_k(:), u_ref(:,:), u_dif_k(:), dash_rho_k(:)!dummy
    real(mp), allocatable :: rho_dep_residual(:)!dummy
    real(mp) :: residual_u0, abs_residual_u0 !dummy

    
    delta = eigenvalues_interval(lx,ly,hx,hy,c1,c2,d1,d2)
    xi = delta(1)/delta(2)
    max_iter = int(log(1.0_mp/eps)/sqrt(xi)) + 1
    rho = (delta(2)-delta(1))/(delta(2)+delta(1))
    omega = 2.0_mp/(1.0_mp + sqrt(1.0_mp-rho**2))
    S = omega - 1.0_mp
    N = size(x)-1
    M = size(y)-1
    allocate(u_prev(0:N,0:M))
    allocate(Lh_u(0:N,0:M),Lh_u0(0:N,0:M),f_ij(0:N,0:M),residual_k(0:max_iter),abs_residual_k(0:max_iter))!dummy
    allocate(u_ref(0:N,0:M),u_dif_k(0:max_iter),rho_dep_residual(0:max_iter),dash_rho_k(0:max_iter))!dummy

    u = u0
    forall (j=0:M) u(0,j) = mu(0.0_mp,y(j))
    forall (j=0:M) u(N,j) = mu(lx,y(j))
    forall (i=1:N-1) u(i,0) = mu(x(i),0.0_mp)
    forall (i=1:N-1) u(i,M) = mu(x(i),ly)
    k = 0
    u_dif_k(0) = 1.0_mp
    do while (k.lt.max_iter)
        k = k+1
        u_prev = u
        do i=1,N-1
            do j=1,M-1
                u(i,j) = u_prev(i,j) + omega*( f(x(i),y(j)) + p(x(i)+hx/2.0_mp,y(j))*(u_prev(i+1,j)-u_prev(i,j))/hx**2 - &
                        & p(x(i)-hx/2.0_mp,y(j))*(u_prev(i,j)-u(i-1,j))/hx**2 + &
                        & q(x(i),y(j)+hy/2.0_mp)*(u_prev(i,j+1)-u_prev(i,j))/hy**2 - &
                        & q(x(i),y(j)-hy/2.0_mp)*(u_prev(i,j)-u(i,j-1))/hy**2 ) /&
                        & ( p(x(i)-hx/2.0_mp,y(j))/hx**2 + p(x(i)+hx/2.0_mp,y(j))/hx**2 +&
                        & q(x(i),y(j)-hy/2.0_mp)/hy**2 + q(x(i),y(j)+hy/2.0_mp)/hy**2 )
  
                Lh_u(i,j) = p(x(i)+hx/2.0_mp,y(j))*(u(i+1,j)-u(i,j))/hx**2 - &  !dummy
                    & p(x(i)-hx/2.0_mp,y(j))*(u(i,j)-u(i-1,j))/hx**2 + &
                    & q(x(i),y(j)+hy/2.0_mp)*(u(i,j+1)-u(i,j))/hy**2 - &
                    & q(x(i),y(j)-hy/2.0_mp)*(u(i,j)-u(i,j-1))/hy**2
                
                f_ij(i,j) = f(x(i),y(j)) !dummy
                u_ref(i,j) = reference_solution(x(i),y(j)) !dummy
            enddo
        enddo
        residual_k(k) = matrix_norm(Lh_u + f_ij) !dummy
        abs_residual_k(k) = matrix_norm(u - u_ref)!dummy
        u_dif_k(k) = matrix_norm(u - u_prev)!dummy
        rho_dep_residual(k) = rho*u_dif_k(k)/(1.0_mp-rho)!dummy
        if (k.lt.2) then !dummy
            dash_rho_k(k) = sqrt(-1.d0) !dummy
        else !dummy
            dash_rho_k(k) = sqrt(u_dif_k(k)/u_dif_k(k-2))!dummy
        endif !dummy
    enddo
    !dummy part!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i=1,N-1
        do j=1,M-1
            Lh_u0(i,j) = p(x(i)+hx/2.0_mp,y(j))*(u0(i+1,j)-u0(i,j))/hx**2 - &
                    & p(x(i)-hx/2.0_mp,y(j))*(u0(i,j)-u0(i-1,j))/hx**2 + &
                    & q(x(i),y(j)+hy/2.0_mp)*(u0(i,j+1)-u0(i,j))/hy**2 - &
                    & q(x(i),y(j)-hy/2.0_mp)*(u0(i,j)-u0(i,j-1))/hy**2
        enddo
    enddo

    residual_u0 = matrix_norm(Lh_u0 + f_ij)
    abs_residual_u0 = matrix_norm(u0 - u_ref)
    write(file_id,*)'1) residual u*:',residual_k(max_iter)
    write(file_id,*)'2) residual u0:',residual_u0
    write(file_id,*)'3) eps:',eps,'m:',max_iter
    write(file_id,*)'4) rho:', rho
    write(file_id,'("    k    &    ||F-A*Uk||    &     rel.d        &     ||Uk-u*||    &     rel.error    &  ||Uk - Uk-1||   &    apost.er.     &    dash_rho_k")')   
    do k=1, max_iter
        write(file_id,'(i5,7("    &    ",e10.3))') k, residual_k(k), residual_k(k)/residual_u0, abs_residual_k(k), abs_residual_k(k)/abs_residual_u0,u_dif_k(k),rho_dep_residual(k),dash_rho_k(k)
    enddo
    write(file_id,*)
    deallocate(Lh_u,Lh_u0,f_ij,residual_k,abs_residual_k)!dummy
    deallocate(u_ref,u_dif_k,rho_dep_residual,dash_rho_k)!dummy
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    deallocate(u_prev)
    end subroutine successive_over_relaxation_method
            
            
    subroutine triangle_matrix_method(p, q, f, mu, u0, x, y, hx, hy, lx, ly, c1, c2, d1, d2, eps, u, file_id)
    interface    
    pure function p(x,y)
    real(8), intent(in) :: x,y
    real(8) :: p
    end function p
    
    pure function q(x,y)
    real(8), intent(in) :: x,y
    real(8) :: q
    end function q
    
    pure function mu(x,y)
    real(8), intent(in) :: x,y
    real(8) :: mu
    end function mu    
    
    pure function f(x,y)
    real(8), intent(in) :: x,y
    real(8) :: f
    end function f
    end interface
    real(mp), intent(in) :: x(0:), y(0:), u0(0:,0:)
    real(mp), intent(in) :: hx, hy ! the grid steps
    real(mp), intent(in) :: lx, ly !endpoints of x and y intervals, respectively
    real(mp), intent(in) :: c1, c2, d1, d2
    real(mp), intent(in) :: eps
    integer, intent(in) :: file_id
    real(mp), intent(out) :: u(0:,0:)
    real(mp) :: rho, S ! radius of method convergence and spectral radius, respectively
    integer  :: N, M
    integer :: max_iter ! maximum of iteration number
    integer :: k ! iteration counter
    integer :: i,j ! grid nodes counters
    real(mp) :: k1, k2
    real(mp) :: eta, delta(1:2), omega, gamma1, gamma2, tau
    real(mp) :: Lh_u_prev
    real(mp), allocatable :: u_prev(:,:)
    real(mp), allocatable :: dash_omega(:,:), omega_k(:,:)
    real(mp), allocatable :: Lh_u(:,:), Lh_u0(:,:),f_ij(:,:), residual_k(:), abs_residual_k(:), u_ref(:,:), u_dif_k(:), dash_rho_k(:)!dummy
    real(mp), allocatable :: rho_dep_residual(:)!dummy
    real(mp) :: residual_u0, abs_residual_u0 !dummy
    delta(1) = c1*4.0_mp*sin(pi*hx/2.0_mp/lx)**2/hx**2 + d1*4.0_mp*sin(pi*hy/2.0_mp/ly)**2/hy**2
    delta(2) = c2*4.0_mp/hx**2 + d2*4.0_mp/hy**2
    eta = delta(1)/delta(2)
    rho = (1.0_mp-eta)/(1.0_mp+eta)
    gamma1 = delta(1)/(2.0_mp + 2.0_mp*sqrt(eta))
    gamma2 = delta(1)/(4.0_mp*sqrt(eta))
    tau = 2.0_mp/(gamma1+gamma2)
    omega = 2.0_mp/sqrt(delta(1)*delta(2))
    k1 = omega/hx**2
    k2 = omega/hy**2
    max_iter = int(log(1.0_mp/eta)/log(1.0_mp/rho)) + 1
    N = size(x)-1
    M = size(y)-1
    allocate(u_prev(0:N,0:M))
    allocate(dash_omega(0:N,0:M), omega_k(0:N,0:M))
    allocate(Lh_u(0:N,0:M),Lh_u0(0:N,0:M),f_ij(0:N,0:M),residual_k(0:max_iter),abs_residual_k(0:max_iter))!dummy
    allocate(u_ref(0:N,0:M),u_dif_k(0:max_iter),rho_dep_residual(0:max_iter),dash_rho_k(0:max_iter))!dummy   

    u = u0
    omega_k = 0.0_mp
    dash_omega = 0.0_mp
    
    forall (j=0:M) u(0,j) = mu(0.0_mp,y(j))
    forall (j=0:M) u(N,j) = mu(lx,y(j))
    forall (i=1:N-1) u(i,0) = mu(x(i),0.0_mp)
    forall (i=1:N-1) u(i,N) = mu(x(i),ly)
    k = 0
    u_dif_k(0) = 1.0_mp !dummy
    do while (k.lt.max_iter)
        k = k+1
        u_prev = u
        do i=1,N-1
            do j=1,M-1
                Lh_u_prev = p(x(i)+hx/2.0_mp,y(j))*(u_prev(i+1,j)-u_prev(i,j))/hx**2 - &
                            & p(x(i)-hx/2.0_mp,y(j))*(u_prev(i,j)-u_prev(i-1,j))/hx**2 + &
                            & q(x(i),y(j)+hy/2.0_mp)*(u_prev(i,j+1)-u_prev(i,j))/hy**2 - &
                            & q(x(i),y(j)-hy/2.0_mp)*(u_prev(i,j)-u_prev(i,j-1))/hy**2
                
                dash_omega(i,j) = ( k1*p(x(i)-hx/2.0_mp,y(j))*dash_omega(i-1,j) + &
                        & k2*q(x(i),y(j)-hy/2.0_mp)*dash_omega(i,j-1) + &
                        & Lh_u_prev + f(x(i),y(j)) ) / ( 1.0_mp + k1*p(x(i)-hx/2.0_mp,y(j)) + &
                        & k2*q(x(i),y(j)-hy/2.0_mp) )
                
                Lh_u(i,j) = p(x(i)+hx/2.0_mp,y(j))*(u(i+1,j)-u(i,j))/hx**2 - &  !dummy
                    & p(x(i)-hx/2.0_mp,y(j))*(u(i,j)-u(i-1,j))/hx**2 + &
                    & q(x(i),y(j)+hy/2.0_mp)*(u(i,j+1)-u(i,j))/hy**2 - &
                    & q(x(i),y(j)-hy/2.0_mp)*(u(i,j)-u(i,j-1))/hy**2
                
                f_ij(i,j) = f(x(i),y(j)) !dummy
                u_ref(i,j) = reference_solution(x(i),y(j)) !dummy
            enddo
        enddo
        do i=N-1,1,-1
            do j=M-1,1,-1
                omega_k(i,j) = ( k1*p(x(i)+hx/2.0_mp,y(j))*omega_k(i+1,j) + &
                               & k2*q(x(i),y(j)+hy/2.0_mp)*omega_k(i,j+1) + &
                               & dash_omega(i,j) ) /&
                                & ( 1.0_mp + k1*p(x(i)+hx/2.0_mp,y(j)) + &
                                & k2*q(x(i),y(j)+hy/2.0_mp) )
            enddo
        enddo
        u = u_prev + tau*omega_k
        residual_k(k) = matrix_norm(Lh_u + f_ij) !dummy
        abs_residual_k(k) = matrix_norm(u - u_ref)!dummy
        u_dif_k(k) = matrix_norm(u - u_prev)!dummy
        rho_dep_residual(k) = rho*u_dif_k(k)/(1.0_mp-rho)!dummy
        if (k.lt.2) then !dummy
            dash_rho_k(k) = sqrt(-1.d0) !dummy
        else !dummy
            dash_rho_k(k) = sqrt(u_dif_k(k)/u_dif_k(k-2))!dummy
        endif !dummy
    enddo
   !dummy part!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i=1,N-1
        do j=1,M-1
            Lh_u0(i,j) = p(x(i)+hx/2.0_mp,y(j))*(u0(i+1,j)-u0(i,j))/hx**2 - &
                    & p(x(i)-hx/2.0_mp,y(j))*(u0(i,j)-u0(i-1,j))/hx**2 + &
                    & q(x(i),y(j)+hy/2.0_mp)*(u0(i,j+1)-u0(i,j))/hy**2 - &
                    & q(x(i),y(j)-hy/2.0_mp)*(u0(i,j)-u0(i,j-1))/hy**2
        enddo
    enddo

    residual_u0 = matrix_norm(Lh_u0 + f_ij)
    abs_residual_u0 = matrix_norm(u0 - u_ref)
    write(file_id,*)'1) residual u*:',residual_k(max_iter)
    write(file_id,*)'2) residual u0:',residual_u0
    write(file_id,*)'3) eps:',eps,'m:',max_iter
    write(file_id,*)'4) rho:', rho
    write(file_id,'("    k    &    ||F-A*Uk||    &     rel.d        &     ||Uk-u*||    &     rel.error    &  ||Uk - Uk-1||   &    apost.er.     &    dash_rho_k")')   
    do k=1, max_iter
        write(file_id,'(i5,7("    &    ",e10.3))') k, residual_k(k), residual_k(k)/residual_u0, abs_residual_k(k), abs_residual_k(k)/abs_residual_u0,u_dif_k(k),rho_dep_residual(k),dash_rho_k(k)
    enddo
    write(file_id,*)
    deallocate(Lh_u,Lh_u0,f_ij,residual_k,abs_residual_k)!dummy
    deallocate(u_ref,u_dif_k,rho_dep_residual,dash_rho_k)!dummy
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    deallocate(u_prev)
    deallocate(dash_omega, omega_k)
    end subroutine triangle_matrix_method
    
end module elliptic_equation_solver