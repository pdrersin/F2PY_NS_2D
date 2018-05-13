!Fortran subroutines to interface with python


!Function 1 - This specifies initial conditions to the domain
subroutine initialize_ic_bc_fort(omega,psi,dx,dy,nx,ny,kappa)
implicit none

double precision, intent(inout), dimension(0:nx-1,0:ny-1) :: omega
double precision, intent(inout), dimension(0:nx-1,0:ny-1) :: psi
double precision, intent(in) :: dx, dy, kappa
integer, intent(in) :: nx, ny

integer :: i,j
double precision :: x, y


do j = 0,ny-1
do i = 0,nx-1

x = dfloat(i)*dx
y = dfloat(j)*dy

omega(i,j) = 2.0*kappa*dcos(kappa*x)*dcos(kappa*y)
psi(i,j) = 1.0/kappa*dcos(kappa*x)*dcos(kappa*y)

end do
end do

end subroutine initialize_ic_bc_fort



subroutine euler_tstep_fort(omega,psi,dt,dx,dy,nx,ny,Re_N)
implicit none

double precision, intent(inout), dimension(0:nx-1,0:ny-1) :: omega
double precision, intent(in), dimension(0:nx-1,0:ny-1) :: psi
double precision, intent(in) :: dt, dx, dy, Re_N
integer, intent(in) :: nx, ny


double precision :: jj1, jj2, jj3, d2wdy2, d2wdx2
double precision, dimension(:,:), allocatable :: omega_new, psi_new
integer :: i, j

allocate(omega_new(-1:nx,-1:ny)) !For ghost points
allocate(psi_new(-1:nx,-1:ny)) !For ghost points

do j = 0, ny-1
do i = 0, nx-1

omega_new(i,j) = omega(i,j)
psi_new(i,j) = psi(i,j)

end do
end do


!BC update
do j = 0,ny-1
omega_new(-1,j) = omega_new(nx-1,j)
omega_new(nx,j) = omega_new(0,j)

psi_new(-1,j) = psi_new(nx-1,j)
psi_new(nx,j) = psi_new(0,j)
end do

do i = -1,nx
omega_new(i,-1) = omega_new(i,ny-1)
omega_new(i,ny) = omega_new(i,0)

psi_new(i,-1) = psi_new(i,ny-1)
psi_new(i,ny) = psi_new(i,0)
end do


!Find Arakawa coefficients
do j = 0,ny-1
do i = 0,nx-1


jj1 = 1.0/(4.0*dx*dy) * ((omega_new(i+1,j)-omega_new(i-1,j)) * (psi_new(i,j+1) - psi_new(i,j-1)) &
			- (omega_new(i,j+1)-omega_new(i,j-1)) * (psi_new(i+1,j) - psi_new(i-1,j)))

jj2 = 1.0 / (4.0 * dx * dy) * (omega_new(i+1, j) * (psi_new(i+1, j+1) - psi_new(i+1, j-1)) &
                                         - omega_new(i-1, j) * (psi_new(i-1, j+1) - psi_new(i-1, j-1)) &
                                         - omega_new(i, j+1) * (psi_new(i+1, j+1) - psi_new(i-1, j+1)) &
                                         + omega_new(i, j-1) * (psi_new(i+1, j-1) - psi_new(i-1, j-1)) &
                                          )

jj3 = 1.0 / (4.0 * dx * dy) * (omega_new(i+1, j+1) * (psi_new(i, j+1) - psi_new(i+1, j)) &
                                        -  omega_new(i-1, j-1) * (psi_new(i-1, j) - psi_new(i, j-1)) &
                                        -  omega_new(i-1, j+1) * (psi_new(i, j+1) - psi_new(i-1, j)) &
                                        +  omega_new(i+1, j-1) * (psi_new(i+1, j) - psi_new(i, j-1)) &
                                          )

d2wdy2 = (omega_new(i, j+1) + omega_new(i, j-1) - 2.0 * omega_new(i, j)) / (dy * dy)
d2wdx2 = (omega_new(i+1, j) + omega_new(i-1, j) - 2.0 * omega_new(i, j)) / (dx * dx)

omega_new(i, j) = omega(i, j) + dt * (-(jj1 + jj2 + jj3) + 1.0 / Re_n * (d2wdy2 + d2wdx2))

end do
end do


do j = 0,ny-1
do i = 0,nx-1
omega(i,j) = omega_new(i,j)
end do
end do

deallocate(omega_new,psi_new)



end subroutine euler_tstep_fort


subroutine gauss_siedel_fort(psi,omega,nx,ny,dx,dy,gs_tol)
implicit none

integer, intent(in) :: nx, ny
double precision, intent(in) :: dx, dy, gs_tol
double precision, intent(in), dimension(0:nx-1,0:ny-1) :: omega
double precision, intent(inout), dimension(0:nx-1,0:ny-1) :: psi

double precision, dimension(:,:), allocatable :: psi_update_1, psi_update_2
integer :: not_converged, i, j, gsiter

allocate(psi_update_1(-1:nx,-1:ny))
allocate(psi_update_2(-1:nx,-1:ny))

do j = 0, ny-1
do i = 0, nx-1

psi_update_1(i,j) = psi(i,j)
psi_update_2(i,j) = psi(i,j)

end do
end do

!BC update
do j = 0,ny-1
psi_update_1(-1,j) = psi_update_1(nx-1,j)
psi_update_1(nx,j) = psi_update_1(0,j)

psi_update_2(-1,j) = psi_update_2(nx-1,j)
psi_update_2(nx,j) = psi_update_2(0,j)
end do

do i = -1,nx
psi_update_1(i,-1) = psi_update_1(i,ny-1)
psi_update_1(i,ny) = psi_update_1(i,0)

psi_update_2(i,-1) = psi_update_2(i,ny-1)
psi_update_2(i,ny) = psi_update_2(i,0)
end do

!GS Iterations
not_converged = 1
gsiter = 0
do while (not_converged == 1)

do j = 0, ny - 1
do i = 0, nx - 1

psi_update_2(i, j) = -0.25*(-omega(i, j) * dx * dx - (psi_update_2(i+1, j) + psi_update_2(i-1, j) &
			+ psi_update_2(i, j+1) + psi_update_2(i, j-1)))

end do
end do

gsiter = gsiter + 1


!Update BC
do j = 0,ny-1
psi_update_2(-1,j) = psi_update_2(nx-1,j)
psi_update_2(nx,j) = psi_update_2(0,j)
end do

do i = -1,nx
psi_update_2(i,-1) = psi_update_2(i,ny-1)
psi_update_2(i,ny) = psi_update_2(i,0)
end do

if (maxval(psi_update_2-psi_update_1) < gs_tol) then

do j = 0, ny-1
do i = 0, nx-1
psi(i,j) = psi_update_2(i,j)
end do
end do

!print*,'Converged: ',gsiter,' Iterations of GS , Residual: ',maxval(psi_update_2-psi_update_1)

deallocate(psi_update_1,psi_update_2)
exit

else


do j = -1,ny
do i = -1,nx
psi_update_1(i,j) = psi_update_2(i,j)
end do
end do

end if

end do

















end subroutine gauss_siedel_fort












