module input_functions
    use precision
    implicit none
    contains
     
    ! -Lu = f(x,y), where (x,y) :: {0 .le. x .le. lx, 0 .le. y .le. ly}
    !               and Lu = d(p(x,y)*du/dx)/dx + d(q(x,y)*du/dy)/dy
    !               also 0 .lt. c1 .le. p(x,y) .le. c2, and 0 .lt. d1 .le q(x,y) .le. d2,
    !               where c1, c2, d1, d2 - constants.
    !mu(x,y)  - boundary condition
    
    pure function reference_solution(x,y) result(u)
    real(mp), intent(in) :: x,y
    real(mp) :: u
    u = x**3*y + x*y**2
    end function reference_solution    
    
    pure function p(x,y)
    real(mp), intent(in) :: x,y
    real(mp) :: p
    p = 1.0_mp
    end function p 
    
    pure function q(x,y)
    real(mp), intent(in) :: x,y
    real(mp) :: q
    q = 1.0_mp
    end function q
    
    pure function mu(x,y)
    real(mp), intent(in) :: x,y
    real(mp) :: mu
    mu = reference_solution(x,y)
    end function mu
    
    pure function f(x,y)
    real(mp), intent(in) :: x,y
    real(mp) :: f
    f = -(6.0_mp*x*y + 2.0_mp*x)
    end function f
    
end module input_functions
    