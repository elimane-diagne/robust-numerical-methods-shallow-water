program hydrostatic_reconstruction_second_order
  implicit none

  ! Paramètres globaux
  integer, parameter :: n_cells = 350, n_steps = 20
  real, parameter :: dx = 1.0 / n_cells, g = 9.81, cfl = 1
  integer :: t_step
  real :: dt, s_max
  real :: h(n_cells), u(n_cells), z(n_cells), x(n_cells)
  real :: h_new(n_cells), u_new(n_cells)
  real :: h_pred(n_cells), u_pred(n_cells)

  ! Initialisation des variables
  call initialize(h, u, z, x, dx, n_cells)
  call save_initial_conditions(x, h, z, n_cells)

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

     ! Sauvegarde des résultats
     if (mod(t_step, 10) == 0) then
        call save_results(t_step, x, h, z, n_cells)
     end if
  end do

contains

  subroutine initialize(h, u, z, x, dx, n_cells)
    integer, intent(in) :: n_cells
    real, intent(in) :: dx
    real, intent(out) :: h(n_cells), u(n_cells), z(n_cells), x(n_cells)
    integer :: i

    ! Initialisation des conditions spécifiques au test 3
    do i = 1, n_cells
        x(i) = (i - 1) * dx
        z(i) = 0.5 * (1.0 - 0.5 * (cos(3.14159 * (x(i) - 0.5) / 0.5) + 1.0))
        h(i) = max(0.0, 0.4 - z(i) + 0.04 * sin((x(i) - 0.5) / 0.25) - max(0.0, -0.4 + z(i)))
        u(i) = 0.0
    end do
  end subroutine initialize

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
    real :: delta_z, source_term_hu, sci_hu, si_hu
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

        ! Hydrostatic reconstruction
        z_interface = max(z_l, z_r)
        h_interface_l = max(0.0, h_l + z_l - z_interface)
        h_interface_r = max(0.0, h_r + z_r - z_interface)

        ! Calcul de ∆z pour le terme source centré Sci
        delta_z = z_r - z_l
        sci_hu = g * (h_l + h_r) / 2.0 * delta_z

        ! Calcul des flux de Rusanov
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

        ! Calcul des termes sources distribués Si
        S_ip1_half_minus = g * 0.5 * (h_interface_l**2 - h_r**2)
        S_im1_half_plus = g * 0.5 * (h_l**2 - h_interface_r**2)
        si_hu = S_ip1_half_minus + S_im1_half_plus
        F_hu_plus = F_hu_plus + S_im1_half_plus
        F_hu_minus = F_hu_minus + S_ip1_half_minus

        
        ! Mise à jour des variables avec flux et termes sources
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

  subroutine save_initial_conditions(x, h, z, n_cells)
    integer, intent(in) :: n_cells
    real, intent(in) :: x(n_cells), h(n_cells), z(n_cells)
    integer :: i
    real :: h_plus_z

    open(unit=10, file="initial_conditions_order2.txt", status="replace")
    write(10, "(A)") "x   z(x)   h(0,x)   h(0,x) + z(x)"
    do i = 1, n_cells
      h_plus_z = h(i) + z(i)
      write(10, "(F8.5, 3F12.5)") x(i), z(i), h(i), h_plus_z
    end do
    close(10)
  end subroutine save_initial_conditions

  subroutine save_results(t_step, x, h, z, n_cells)
    integer, intent(in) :: t_step, n_cells
    real, intent(in) :: x(n_cells), h(n_cells), z(n_cells)
    integer :: i
    real :: h_plus_z
    character(len=50) :: filename

    write(filename, "(A,I0)") "results_order2_step_", t_step
    open(unit=20, file=trim(filename) // ".txt", status="replace")
    write(20, "(A)") "x   z(x)   h   h + z"
    do i = 1, n_cells
      h_plus_z = h(i) + z(i)
      write(20, "(F8.5, 3F12.5)") x(i), z(i), h(i), h_plus_z
    end do
    close(20)
  end subroutine save_results

end program hydrostatic_reconstruction_second_order
