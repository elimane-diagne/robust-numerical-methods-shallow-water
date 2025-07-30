program hydrostatic_reconstruction
  implicit none

  ! Paramètres globaux
  integer, parameter :: n_cells = 200, n_steps = 1
  real, parameter :: dx = 0.1 , g = 9.81, cfl = 1
  integer :: t_step
  real :: dt, s_max
  real :: h(n_cells), u(n_cells), z(n_cells)
  real :: h_new(n_cells), u_new(n_cells)
  real :: x 

  ! Initialisation des variables
  call initialize(h, u, z, dx , n_cells)

  ! Boucle temporelle
  do t_step = 1, n_steps

     ! Calcul de la vitesse maximale d'onde pour le pas de temps
     call calculate_max_wave_speed(h, u, z, n_cells, g, s_max)
     dt = cfl * dx / s_max

     ! Reconstruction hydrostatique et calcul des flux
     call reconstruction_and_flux(h, u, z, n_cells, dx, dt, g, h_new, u_new)

     ! Copier les nouvelles valeurs dans les anciennes
     h = h_new
     u = u_new

     ! Affichage des résultats
     call print_results(h, u, n_cells, z)
  end do

contains

  subroutine initialize(h, u, z, dx , n_cells)
    integer, intent(in) :: n_cells
    real, intent(in) :: dx
    real, intent(out) :: h(n_cells), u(n_cells), z(n_cells)
    integer :: i
    real :: x 

    do i = 1, n_cells
        z(i) = 0.5
        h(i) = 1.0 - z(i)  
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
    implicit none
    integer, intent(in) :: n_cells
    real, intent(in) :: h(n_cells), u(n_cells), z(n_cells), dx, dt, g
    real, intent(out) :: h_new(n_cells), u_new(n_cells)

    integer :: i
    real :: h_i_plus_half_G, h_i_plus_half_D, u_i_plus_half_G, u_i_plus_half_D
    real :: h_i_minus_half_G, h_i_minus_half_D, u_i_minus_half_G, u_i_minus_half_D
    real :: F_h_i_plus_half, F_hu_i_plus_half
    real :: F_h_i_minus_half, F_hu_i_minus_half
    real :: z_max_i_plus_half, z_max_i_minus_half
    real :: source_term_i_plus_half, source_term_i_minus_half

    ! Conditions aux limites : murs solides
    h_new(1) = h(1)
    u_new(1) = 0.0
    h_new(n_cells) = h(n_cells)
    u_new(n_cells) = 0.0

    ! Calcul des flux et des termes sources
    do i = 2, n_cells - 1
        ! Indices \(i+1/2\)
        z_max_i_plus_half = max(z(i), z(i+1))
        h_i_plus_half_G = max(0.0, h(i) + z(i) - z_max_i_plus_half)  ! Gauche
        h_i_plus_half_D = max(0.0, h(i+1) + z(i+1) - z_max_i_plus_half)  ! Droite
        u_i_plus_half_G = u(i)
        u_i_plus_half_D = u(i+1)

        ! Indices \(i-1/2\)
        z_max_i_minus_half = max(z(i-1), z(i))
        h_i_minus_half_G = max(0.0, h(i-1) + z(i-1) - z_max_i_minus_half)  ! Gauche
        h_i_minus_half_D = max(0.0, h(i) + z(i) - z_max_i_minus_half)  ! Droite
        u_i_minus_half_G = u(i-1)
        u_i_minus_half_D = u(i)

        ! Terme source pour \(i+1/2\) et \(i-1/2\)
        source_term_i_plus_half = g * (h_i_plus_half_G + h_i_plus_half_D) * (z(i+1) - z(i)) / 2.0
        source_term_i_minus_half = g * (h_i_minus_half_G + h_i_minus_half_D) * (z(i) - z(i-1)) / 2.0

        ! Flux à \(i+1/2\)
        F_h_i_plus_half = 0.5 * (h_i_plus_half_G * u_i_plus_half_G + h_i_plus_half_D * u_i_plus_half_D) - &
                        0.5 * max(abs(u_i_plus_half_G) + sqrt(g * h_i_plus_half_G), &
                                    abs(u_i_plus_half_D) + sqrt(g * h_i_plus_half_D)) * &
                                    (h_i_plus_half_D - h_i_plus_half_G)
                                    
        F_hu_i_plus_half = 0.5 * (h_i_plus_half_G * u_i_plus_half_G**2 + 0.5 * g * h_i_plus_half_G**2 + &
                                h_i_plus_half_D * u_i_plus_half_D**2 + 0.5 * g * h_i_plus_half_D**2) - &
                        0.5 * max(abs(u_i_plus_half_G) + sqrt(g * h_i_plus_half_G), &
                                    abs(u_i_plus_half_D) + sqrt(g * h_i_plus_half_D)) * &
                                    (h_i_plus_half_D * u_i_plus_half_D - h_i_plus_half_G * u_i_plus_half_G)

        ! Flux à \(i-1/2\)
        F_h_i_minus_half = 0.5 * (h_i_minus_half_G * u_i_minus_half_G + h_i_minus_half_D * u_i_minus_half_D) - &
                        0.5 * max(abs(u_i_minus_half_G) + sqrt(g * h_i_minus_half_G), &
                                    abs(u_i_minus_half_D) + sqrt(g * h_i_minus_half_D)) * &
                                    (h_i_minus_half_D - h_i_minus_half_G)
        F_hu_i_minus_half = 0.5 * (h_i_minus_half_G * u_i_minus_half_G**2 + 0.5 * g * h_i_minus_half_G**2 + &
                                h_i_minus_half_D * u_i_minus_half_D**2 + 0.5 * g * h_i_minus_half_D**2) - &
                            0.5 * max(abs(u_i_minus_half_G) + sqrt(g * h_i_minus_half_G), &
                                    abs(u_i_minus_half_D) + sqrt(g * h_i_minus_half_D)) * &
                                    (h_i_minus_half_D * u_i_minus_half_D - h_i_minus_half_G * u_i_minus_half_G)

        F_hu_i_plus_half = F_hu_i_plus_half + source_term_i_plus_half
        F_hu_i_minus_half = F_hu_i_minus_half + source_term_i_minus_half

        ! Mise à jour des états avec le terme source
        h_new(i) = h(i) - dt / dx * (F_h_i_plus_half - F_h_i_minus_half)
        if (h_new(i) > 1e-6) then
        u_new(i) = (h(i) * u(i) - dt / dx * (F_hu_i_plus_half - F_hu_i_minus_half) ) / h_new(i)
                    
        else
        u_new(i) = 0.0
        end if
    end do
 end subroutine reconstruction_and_flux

 subroutine print_results(h, u, n_cells, z)
    integer, intent(in) :: n_cells
    real, intent(in) :: h(n_cells), u(n_cells), z(n_cells)
    integer :: i

    print *, "Cell | h + z | h | z"
    do i = 1, n_cells
    print *, i, h(i) + z(i), h(i), z(i)
    end do
 end subroutine print_results

end program hydrostatic_reconstruction
