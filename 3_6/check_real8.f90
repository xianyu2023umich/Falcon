program check_real8
    implicit none
    real(8) :: x, eps, tiny_val, small1, small2, small3
    
    eps = epsilon(x)
    tiny_val = tiny(x)
    
    write(*,*) '=== real(8) properties ==='
    write(*,*) 'Maximum value (huge):', huge(x)
    write(*,*) 'Minimum positive normalized (tiny):', tiny_val
    write(*,*) 'Machine epsilon (epsilon):', eps
    write(*,*) 'Decimal precision (digits):', digits(x)
    write(*,*) ''
    write(*,*) '=== Testing values smaller than epsilon ==='
    write(*,*) ''
    
    ! Values much smaller than epsilon but still representable
    small1 = 1.0d-100  ! Much smaller than epsilon (2.22e-16)
    small2 = 1.0d-200  ! Even smaller
    small3 = 1.0d-300  ! Very small but still > tiny
    
    write(*,*) 'small1 = 1.0e-100 =', small1
    write(*,*) '  Is it representable?', (small1 > 0.0d0)
    write(*,*) '  Is it < epsilon?', (small1 < eps)
    write(*,*) ''
    
    write(*,*) 'small2 = 1.0e-200 =', small2
    write(*,*) '  Is it representable?', (small2 > 0.0d0)
    write(*,*) '  Is it < epsilon?', (small2 < eps)
    write(*,*) ''
    
    write(*,*) 'small3 = 1.0e-300 =', small3
    write(*,*) '  Is it representable?', (small3 > 0.0d0)
    write(*,*) '  Is it < epsilon?', (small3 < eps)
    write(*,*) '  Is it > tiny?', (small3 > tiny_val)
    write(*,*) ''
    
    write(*,*) '=== What epsilon actually means ==='
    write(*,*) '1.0 + epsilon =', 1.0d0 + eps
    write(*,*) '1.0 + (epsilon/2) =', 1.0d0 + eps/2.0d0
    write(*,*) '  (epsilon/2) < epsilon, but 1.0 + (epsilon/2) == 1.0?', &
               (1.0d0 + eps/2.0d0 == 1.0d0)
    write(*,*) ''
    write(*,*) 'But you CAN represent values smaller than epsilon:'
    write(*,*) '  Just not as additions to 1.0!'
    
end program check_real8