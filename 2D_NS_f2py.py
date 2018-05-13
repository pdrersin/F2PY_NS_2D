import numpy as np
import matplotlib.pyplot as plt
import time
from numpy import f2py

import Fortran_functions

#Python functions
def euler_tstep(omega,psi):

    omega_new = np.copy(omega)

    for i in range(0,nx):
        for j in range(0,ny):

            im1 = i - 1
            ip1 = i + 1
            jm1 = j - 1
            jp1 = j + 1

            if im1 < 0 :
                im1 = im1 + nx

            if jm1 < 0 :
                jm1 = jm1 + ny

            if ip1 > nx - 1:
                ip1 = ip1 - nx

            if jp1 > ny - 1:
                jp1 = jp1 - ny

            # FTCS
            # dsdy = (psi[i, jp1] - psi[i, jm1]) / (2.0 * dy)
            # dsdx = (psi[ip1, j] - psi[im1, j]) / (2.0 * dx)
            #
            # dwdy = (omega[i, jp1] - omega[i, jm1]) / (2.0 * dy)
            # dwdx = (omega[ip1, j] - omega[im1, j]) / (2.0 * dx)

            # Arakawa
            jj1 = 1.0 / (4.0 * dx * dy) * ((omega[ip1, j] - omega[im1, j]) * (psi[i, jp1] - psi[i, jm1])
                                         - (omega[i, jp1] - omega[i, jm1]) * (psi[ip1, j] - psi[im1, j]))

            jj2 = 1.0 / (4.0 * dx * dy) * (omega[ip1, j] * (psi[ip1, jp1] - psi[ip1, jm1])
                                         - omega[im1, j] * (psi[im1, jp1] - psi[im1, jm1])
                                         - omega[i, jp1] * (psi[ip1, jp1] - psi[im1, jp1])
                                         + omega[i, jm1] * (psi[ip1, jm1] - psi[im1, jm1])
                                          )

            jj3 = 1.0 / (4.0 * dx * dy) * (omega[ip1, jp1] * (psi[i, jp1] - psi[ip1, j])
                                        -  omega[im1, jm1] * (psi[im1, j] - psi[i, jm1])
                                        -  omega[im1, jp1] * (psi[i, jp1] - psi[im1, j])
                                        +  omega[ip1, jm1] * (psi[ip1, j] - psi[i, jm1])
                                          )

            d2wdy2 = (omega[i, jp1] + omega[i, jm1] - 2.0 * omega[i, j]) / (dy * dy)
            d2wdx2 = (omega[ip1, j] + omega[im1, j] - 2.0 * omega[i, j]) / (dx * dx)

            omega_new[i, j] = omega[i, j] + dt * (-(jj1 + jj2 + jj3) + 1.0 / Re_n * (d2wdy2 + d2wdx2))


    for i in range(nx):
        for j in range(ny):
            omega[i,j] = omega_new[i,j]

    del omega_new

def gauss_seidel(psi,omega):

    psi_update_1 = np.copy(psi)
    psi_update_2 = np.copy(psi)

    not_converged = True

    while not_converged:
        for i in range(nx):
            for j in range(ny):

                im1 = i - 1
                ip1 = i + 1
                jm1 = j - 1
                jp1 = j + 1

                if im1 < 0:
                    im1 = im1 + nx

                if jm1 < 0:
                    jm1 = jm1 + ny

                if ip1 > nx - 1:
                    ip1 = ip1 - nx

                if jp1 > ny - 1:
                    jp1 = jp1 - ny

                psi_update_2[i, j] = -0.25*(-omega[i, j] * dx * dx - (
                            psi_update_2[ip1, j] + psi_update_2[im1, j] + psi_update_2[i, jp1] + psi_update_2[i, jm1]))

        if np.max(psi_update_2-psi_update_1)<gs_tol:

            for i in range(nx):
                for j in range(ny):
                    psi[i,j] = psi_update_2[i,j]

            del psi_update_1, psi_update_2

            break

        else:
            for i in range(nx):
                for j in range(ny):
                    psi_update_1[i,j] = psi_update_2[i,j]

def initialize_ic_bc(omega,psi):

    for i in range(nx):
        for j in range(ny):

            x = float(i) * dx
            y = float(j) * dy

            omega[i, j] = 2.0 * kappa * np.cos(kappa * x) * np.cos(kappa * y)
            psi[i, j] = 1.0 / kappa * np.cos(kappa * x) * np.cos(kappa * y)

def init_domain():

    global nx, ny, kappa, Re_n, lx, ly, lt, nt, dt, dx, dy, gs_tol
    global use_fortran
    nx = 64
    ny = 64
    kappa = 2.0
    Re_n = 1.0
    lx = 2.0 * np.pi
    ly = 2.0 * np.pi
    dx = lx / float(nx)
    dy = ly / float(ny)

    gs_tol = 1.0e-4


    use_fortran = 0#if 0 - python functions, 1 - fortran

    lt = 0.1
    dt = 5.0e-4
    nt = int(lt/dt)

    omega = np.zeros(shape=(nx,ny),dtype='double', order='F')
    psi = np.zeros(shape=(nx, ny), dtype='double', order='F')

    return omega, psi

def post_process(omega):

    fig, ax = plt.subplots(nrows=2,ncols=1)

    levels = np.linspace(-4,4,10)

    ax[0].set_title("Numerical Solution")
    plot1 = ax[0].contourf(omega[:, :],levels=levels)
    plt.colorbar(plot1, format="%.2f",ax=ax[0])

    omega_true = np.zeros(shape=(nx,ny),dtype='double')

    for i in range(nx):
        for j in range(ny):

            x = float(i) * dx
            y = float(j) * dy

            omega_true[i, j] = 2.0 * kappa * np.cos(kappa * x) * np.cos(kappa * y)*np.exp(-2.0*kappa*kappa*lt/Re_n)

    ax[1].set_title("Exact Solution")
    plot2 = ax[1].contourf(omega_true[:, :],levels=levels)
    plt.colorbar(plot2, format="%.2f",ax=ax[1])
    plt.show()

def main_func():

    #Initialize my domain and constants
    omega, psi = init_domain()

    if use_fortran == 0:
    	#Initial and boundary conditions - python
        initialize_ic_bc(omega,psi)
    else:
    	# Initial and boundary conditions - Fortran
        Fortran_functions.initialize_ic_bc_fort(omega, psi, dx, dy, kappa)


    t = 0.0

    clock_time_init = time.clock()

    for tstep in range(nt):

        t = t + dt
        print(t)

        if use_fortran == 0:

            #FTCS update - python
            euler_tstep(omega,psi)
            #Gauss-Siedel for Vort-Stream Poisson
            gauss_seidel(psi,omega)
        else:
        	#FTCS update - Fortran
            Fortran_functions.euler_tstep_fort(omega,psi,dt,dx,dy,Re_n)
            # Gauss-Siedel Prediction - Fortran
            Fortran_functions.gauss_siedel_fort(psi, omega, dx, dy, gs_tol)
       

    total_clock_time = time.clock()-clock_time_init

    print('Total Clock Time = ',total_clock_time)

    post_process(omega)


##### RUN HERE #####
main_func()

