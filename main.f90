module global
    implicit none
    real(kind=8)                    ::     data    
    real(kind=8)                    ::     dx, lu_value
contains   
end module global

! Boundary value processing
subroutine Boundary_value(u)
    real(kind=8), dimension(1:N+7)  ::     u
    u(4) = u(N+4)
    u(3) = u(N+3)
    u(2) = u(N+2)
    u(1) = u(N+1)
    u(N+5) = u(5)
    u(N+6) = u(6)
    u(N+7) = u(7)
end subroutine Boundary_value

subroutine Lu(uold, j, umax)
    use global
    real(kind=8), dimension(1:N+7)  ::     uold
    real(kind=8)                    ::     flax_data, diffision_data, umax
    integer                         ::     j
    flax_data = 1.0d0
    write(*,*) flax_data
    call flax(uold, j, umax, flax_data)
    write(*,*) flax_data

    diffision_data = 1.0d0
    write(*,*) diffision_data
    call diffision(uold, j, diffision_data)
    write(*,*) diffision_data
end subroutine Lu

   
subroutine flax(uold, j, umax, flax_data)
    use global
    real(kind=8)                    ::     upp, upm, umm, ump, flax1, flax2, flax_data, umax 
    real(kind=8), dimension(1:N+7)  ::     uold
    integer                         ::     j
    call weno(uold(j+3), uold(j+2), uold(j+1), uold(j), uold(j-1))
    upp = data
    call weno(uold(j-2), uold(j-1), uold(j), uold(j+1), uold(j+2))
    upm = data
    call weno(uold(j-3), uold(j-2), uold(j-1), uold(j), uold(j+1))
    umm = data
    call weno(uold(j+2), uold(j+1), uold(j), uold(j-1), uold(j-2))
    ump = data
    flax1 = 1.0d0/2*((upm/2.0d0)**2+(upp/2.0d0)**2-umax*(upp-pm))
    flax2 = 1.0d0/2*((umm/2.0d0)**2+(ump/2.0d0)**2-umax*(ump-mm))
    flax_data = -(flax1-flax2)/dx
end subroutine flax

subroutine diffision(uold, j, diffision_data)
    real(kind=8), dimension(1:N+7)  ::     uold
    real(kind=8)                    ::     diffision_data
    integer                         ::     j
    diffision_data = 1.0d0/12*uold(j-1) - 5.0d0/4*uold(j) + 5.0d0/4*uold(j+1) - 1.0d0/12*uold(j+2)

end subroutine diffision

! for x=j+1/2                                               for x=j-1/2
! minal a:u(j-2) b:u(j-1) c:u(j) d:u(j+1) e:u(j+2)          minal a:u(j-3) b:u(j-2) c:u(j-1) d:u(j) e:u(j+1)
! plus  a:u(j+3) b:u(j+2) c:u(j+1) d:u(j) e:u(j-1)          plus  a:u(j+2) b:u(j+1) c:u(j) d:u(j-1) e:u(j-2)
subroutine weno(a,b,c,d,e)
    use global 
    real(kind=8)                    ::     u1, u2, u3, beta1, beta2, beta3, w1, w2, w3, wf1, wf2, wf3, a, b, c, d, e
    beta1 = 13.0d0/12*(a-2*b+c)**2 + 1.0d0/4*(a-4*b+3*c)**2
    beta2 = 13.0d0/12*(b-2*c+d)**2 + 1.0d0/4*(b-d)**2
    beta3 = 13.0d0/12*(c-2*d+e)**2 + 1.0d0/4*(3*c-4*d+e)**2
    wf1   = 1.0d0/10/(beta1+1E-06)**2
    wf2   = 3.0d0/5/(beta2+1E-06)**2
    wf3   = 3.0d0/10/(beta3+1E-06)**2
    w1    = wf1/(wf1+wf2+wf3)
    w2    = wf2/(wf1+wf2+wf3)
    w3    = wf3/(wf1+wf2+wf3)
    u1    = 1.0d0/3*a   - 7.0d0/6*b  + 11.0d0/6*c
    u2    = -1.0d0/6*b  + 5.0d0/6*c  + 1.0d0/3*d
    u3    = 1.0d0/3*c   + 5.0d0/6*d  - 1.0d0/6*e
    data  = w1*u1 + w2*u2 + w3*u3

end subroutine weno

program main
    use global
    implicit none
    integer, parameter              ::     N = 200, L = 1
    integer                         ::     j, it
    real(kind=8), dimension(1:N+7)  ::     uold, un, u1, u2
    real(kind=8)                    ::     dt, t, pi, umax

    dx = 1.0d0*L/N
    dt = 0.005d0
    pi = 4.0d0*atan(1.0d0)
    t  = 0.0d0
    
! initial numerical
    do j = 5, N+4
        un(j) = cos(2.0d0*pi*(j-2.5d0))/(2*pi) + j - 2.5 - (cos(2.0d0*pi*(j-3.5d0))/(2*pi) + j - 3.5)
        un(j) = un(j) / dx
      !  write(*,*) (j-2.5)*dx, (j-3.5)*dx
    end do
    call Boundary_value(un)
    uold = un
    write(*,*) data
    do it = 1, 50
        umax = maxval(uold)
        do j = 5, N+4        
            call Lu(uold, j, umax)    
            u1(j) = uold(j) + dt*lu_value  
        end do
        call Boundary_value(u1)

        umax = maxval(u1)
        do j = 5, N+4        
            call Lu(u1, j, umax)    
            u2(j) = 3.0d0/4*uold(j) + 1.0d0/4*u1(j) + 1.0d0/4*dt*lu_value  
        end do
        call Boundary_value(u2)

        umax = maxval(u2)
        do j = 5, N+4        
            call Lu(u2, j, umax)    
            un(j) = 1.0d0/3*uold(j) + 2.0d0/3*u2(j) + 2.0d0/3*dt*lu_value  
        end do
        call Boundary_value(u2)

        uold = un
    
    end do


end program main