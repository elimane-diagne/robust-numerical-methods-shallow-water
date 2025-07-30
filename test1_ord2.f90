program hydrostatic_reconstruction_second_order
  implicit none

  ! Paramètres globaux
  integer, parameter :: n_cells = 100, n_steps = 1
  real, parameter :: L = 2.4, g = 9.81, cfl = 1
  integer :: t_step
  real :: dt, s_max, dx
  real :: h(n_cells), u(n_cells), z(n_cells), x(n_cells)
  real :: h_new(n_cells), u_new(n_cells)
  real :: h_pred(n_cells), u_pred(n_cells)
  dx = L / (n_cells - 1)

  ! Initialisation des variables
  call initialize_test1(h, u, z, x, dx, n_cells)

  ! Boucle temporelle
  do t_step = 1, n_steps

     ! Calcul du pas de temps
     call calculate_max_wave_speed(h, u, z, n_cells, g, s_max)
     dt = cfl * dx / s_max

     ! Étape de prédiction
     call reconstruction_and_flux(h, u, z, n_cells, dx, dt, g, h_pred, u_pred)

     ! Étape de correction (Heun)
     call reconstruction_and_flux(h_pred, u_pred, z, n_cells, dx, dt, g, h_new, u_new)

     ! Moyenne des prédictions pour monter à l'ordre 2
     h_new = 0.5 * (h + h_new)
     u_new = 0.5 * (u + u_new)

     ! Mise à jour des valeurs
     h = h_new
     u = u_new

     call write_results(h, u, z, n_cells, t_step)
  end do

contains

! Fonction pour convertir un entier en chaîne de caractères
  function itoa(num) result(str)
    integer, intent(in) :: num
    character(len=32) :: str
    write(str, '(I32)') num
  end function itoa

  subroutine initialize_test1(h, u, z, x, dx, n_cells)
      integer, intent(in) :: n_cells
      real, intent(in) :: dx
      real, intent(out) :: h(n_cells), u(n_cells), z(n_cells), x(n_cells)
      integer :: i
      real :: pos

      do i = 1, n_cells
        pos = (i - 1) * dx
        x(i) = pos

        ! Topographie 
        if (pos >= 1.4 .and. pos <= 1.6) then
          z(i) = 0.25 * (cos(3.14159 * (pos - 1.5) / 0.1) + 1.0)
        else
          z(i) = 0.0
        end if

        ! Perturbation dans l'intervalle donné
        if (pos >= 1.1 .and. pos <= 1.2) then
          h(i) = 1.001  
        else
          h(i) = 1.0 - z(i)  
        end if

        ! Vitesse initiale nulle
        u(i) = 0.0
      end do
    end subroutine initialize_test1

  subroutine calculate_max_wave_speed(h, u, z, n_cells, g, s_max)
    integer, intent(in) :: n_cells
    real, intent(in) :: h(n_cells), u(n_cells), z(n_cells), g
    real, intent(out) :: s_max
    integer :: i
    real :: wave_speed

    s_max = 0.0
    do i = 1, n_cells
      wave_speed = abs(u(i)) + sqrt(g * max(h(i), 1e-6))
      s_max = max(s_max, wave_speed)
    end do
  end subroutine calculate_max_wave_speed

  subroutine reconstruction_and_flux(h, u, z, n_cells, dx, dt, g, h_new, u_new)
    integer, intent(in) :: n_cells
    real, intent(in) :: h(n_cells), u(n_cells), z(n_cells), dx, dt, g
    real, intent(out) :: h_new(n_cells), u_new(n_cells)

    integer :: i
    real :: h_l, h_r, z_l, z_r, z_interface
    real :: grad_h, grad_z
    real :: F_h_plus, F_h_minus, F_hu_plus, F_hu_minus, lambda_max
    real :: h_interface_l, h_interface_r
    real :: delta_z, sci_hu, si_hu
    real :: S_ip1_half_minus, S_im1_half_plus

    ! Conditions aux limites : 
    h_new(1) = h(1)
    u_new(1) = 0.0
    h_new(n_cells) = h(n_cells)
    u_new(n_cells) = 0.0

    do i = 2, n_cells - 1
        ! Calcul des gradients avec minmod
        grad_h = minmod((h(i+1) - h(i)) / dx, (h(i) - h(i-1)) / dx)
        grad_z = minmod((z(i+1) - z(i)) / dx, (z(i) - z(i-1)) / dx)

        ! Reconstruction au second ordre
        h_l = h(i) - 0.5 * dx * grad_h
        h_r = h(i) + 0.5 * dx * grad_h
        z_l = z(i) - 0.5 * dx * grad_z
        z_r = z(i) + 0.5 * dx * grad_z

        ! Reconstruction hydrostatique
        z_interface = max(z_l, z_r)
        h_interface_l = max(0.0, h_l + z_l - z_interface)
        h_interface_r = max(0.0, h_r + z_r - z_interface)

        ! Calcul de ∆z pour le terme source centré Sci
        delta_z = z_r - z_l
        sci_hu = g * (h_l + h_r) / 2.0 * delta_z

        ! Calcul des flux (Rusanov)
        lambda_max = max(abs(u(i)) + sqrt(g * h_interface_l), abs(u(i+1)) + sqrt(g * h_interface_r))

        F_h_plus = 0.5 * (h_interface_l * u(i) + h_interface_r * u(i+1)) - &
                   0.5 * lambda_max * (h_interface_r - h_interface_l)

        F_hu_plus = 0.5 * (h_interface_l * u(i)**2 + 0.5 * g * h_interface_l**2 + &
                           h_interface_r * u(i+1)**2 + 0.5 * g * h_interface_r**2) - &
                    0.5 * lambda_max * (h_interface_r * u(i+1) - h_interface_l * u(i))

        lambda_max = max(abs(u(i-1)) + sqrt(g * h_interface_l), abs(u(i)) + sqrt(g * h_interface_r))

        F_h_minus = 0.5 * (h_interface_l * u(i-1) + h_interface_r * u(i)) - &
                    0.5 * lambda_max * (h_interface_r - h_interface_l)

        F_hu_minus = 0.5 * (h_interface_l * u(i-1)**2 + 0.5 * g * h_interface_l**2 + &
                            h_interface_r * u(i)**2 + 0.5 * g * h_interface_r**2) - &
                     0.5 * lambda_max * (h_interface_r * u(i) - h_interface_l * u(i-1))

        ! Calcul des termes sources distribués 
        S_ip1_half_minus = g * 0.5 * (h_interface_l**2 - h_r**2)
        S_im1_half_plus = g * 0.5 * (h_l**2 - h_interface_r**2)

        si_hu = S_ip1_half_minus + S_im1_half_plus

        ! Appliqué dans les flux
        F_hu_plus = F_hu_plus + S_im1_half_plus
        F_hu_minus = F_hu_minus + S_ip1_half_minus

        ! Mise à jour des variables
        h_new(i) = h(i) - dt / dx * (F_h_plus - F_h_minus)
        if (h_new(i) > 1e-6) then
            u_new(i) = (h(i) * u(i) - dt / dx * ((F_hu_plus - F_hu_minus) + sci_hu )) / h_new(i)
        else
            u_new(i) = 0.0
        end if
    end do
  end subroutine reconstruction_and_flux

  function minmod(a, b)
    real, intent(in) :: a, b
    real :: minmod

    if (a * b > 0.0) then
      minmod = sign(1.0, a) * min(abs(a), abs(b))
    else
      minmod = 0.0
    end if
  end function minmod

  subroutine write_results(h, u, z, n_cells, t_step)
    integer, intent(in) :: n_cells, t_step
    real, intent(in) :: h(n_cells), u(n_cells), z(n_cells)
    
    integer :: i
    character(len=50) :: filename
     write(filename, "(A,I0)") "results_t1_order2_step_", t_step
    open(unit=10, file=filename, status="unknown")
    write(10, *) "Cell | h |  z| h +z"
    do i = 1, n_cells
      write(10, '(I5,5F15.6)') i, (i - 1) * (2.4 / (n_cells - 1)), h(i), z(i), h(i) + z(i)
    end do
    close(10)
  end subroutine write_results

end program hydrostatic_reconstruction_second_order
