subroutine reconduction(uold)
    real(kind=8), dimension(1:N+5)  ::     uold
    write(*,*) uold(1)

end subroutine reconduction

subroutine weno(u, j)
    real(kind=8), dimension(1:N+5)  ::     u
    integer                         ::     j
    real(kind=8)                    ::     beta1, beta2, beta3, w1, w2, w3, wf1, wf2, wf3
    beta1 = 13/12*(u(j-2)-2*u(j-1)+u(j))**2 + 1/4*(u(j-2)-4*u(j-1)+3*u(j))**2
    beta2 = 13/12*(u(j-1)-2*u(j)+u(j+1))**2 + 1/4*(u(j-1)-u(j+1))**2
    beta3 = 13/12*(u(j)-2*u(j+1)+u(j+2))**2 + 1/4*(3*u(j)-4*u(j+1)+u(j+2))**2
    wf1   = 1/10/(beta1+1E-06)**2
    wf2   = 3/5/(beta2+1E-06)**2
    wf3   = 3/10/(beta3+1E-06)**2

end subroutine weno

program main
    implicit none
    integer, parameter        ::     N = 200, L = 1
    integer                   ::     j, it
    real(kind=8), dimension(1:N+5)  ::     uold, un
    real(kind=8)                    ::     dx, dt, x, t, pi

    dx = 1.0d0*L/N
    dt = 0.005d0
    pi = 4.0d0*atan(1.0d0)
    t  = 0.0d0

! initial numerical
    do j = 4, N+3
        un(j) = cos(2.0d0*pi*(j-2.5d0))/(2*pi) + j - 2.5 - (cos(2.0d0*pi*(j-3.5d0))/(2*pi) + j - 3.5)
        write(*,*) un(j)
        un(j) = un(j) / dx
        write(*,*) (j-2.5)*dx, (j-3.5)*dx
    end do
    un(3) = un(N+3)
    un(2) = un(N+2)
    un(1) = un(N+1)
    un(N+4) = un(4)
    un(N+5) = un(5)
    uold = un
    call reconduction(uold)
    do it = 1, 50

    end do


end program main