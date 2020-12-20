module function_interfaces
    implicit none
    
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
    
end module function_interfaces