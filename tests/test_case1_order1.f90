program comparison_schemes
  implicit none

  ! Paramètres globaux
  integer, parameter :: n_cells = 100, n_steps = 1
  real, parameter :: L = 2.4, g = 9.81, cfl = 0.9
  real :: dx, dt, s_max, t
  integer :: t_step
  real :: h_standard(n_cells), u_standard(n_cells), z(n_cells), x(n_cells)
  real :: h_balanced(n_cells), u_balanced(n_cells)
  real :: h_new(n_cells), u_new(n_cells)

  ! Initialisation
  dx = L / (n_cells - 1)
  t = 0.0
  call initialize_test2(h_standard, u_standard, z, x, dx, n_cells)
  h_balanced = h_standard
  u_balanced = u_standard

  ! Boucle temporelle
  do t_step = 1, n_steps

    ! Calcul du pas de temps
    call calculate_max_wave_speed(h_standard, u_standard, z, n_cells, g, s_max)
    dt = cfl * dx / s_max
    t = t + dt

    ! Schéma standard
    call standard_flux(h_standard, u_standard, z, n_cells, dx, dt, g, h_new, u_new)
    h_standard = h_new
    u_standard = u_new

    ! Schéma bien équilibré
    call balanced_flux(h_balanced, u_balanced, z, n_cells, dx, dt, g, h_new, u_new)
    h_balanced = h_new
    u_balanced = u_new

    
    ! Écriture des résultats
    call write_results(h_standard, u_standard, z, n_cells, t_step, "standard")
    call write_results(h_balanced, u_balanced, z, n_cells, t_step, "balanced")

  end do

contains


    ! Fonction pour convertir un entier en chaîne de caractères
    function itoa(num) result(str)
        integer, intent(in) :: num
        character(len=32) :: str
        write(str, '(I32)') num
    end function itoa


    subroutine initialize_test2(h, u, z, x, dx, n_cells)
      integer, intent(in) :: n_cells
      real, intent(in) :: dx
      real, intent(out) :: h(n_cells), u(n_cells), z(n_cells), x(n_cells)
      integer :: i
      real :: pos

      do i = 1, n_cells
        pos = (i - 1) * dx
        x(i) = pos

        ! Topographie conforme au document
        if (pos >= 1.4 .and. pos <= 1.6) then
          z(i) = 0.25 * (cos(3.14159 * (pos - 1.5) / 0.1) + 1.0)
        else
          z(i) = 0.0
        end if

        ! Perturbation dans l'intervalle donné
        if (pos >= 1.1 .and. pos <= 1.2) then
          h(i) = 1.001  ! Hauteur perturbée
        else
          h(i) = 1.0 - z(i)  ! Hauteur d'équilibre
        end if

        ! Vitesse initiale nulle
        u(i) = 0.0
      end do
    end subroutine initialize_test2

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

    subroutine standard_flux(h, u, z, n_cells, dx, dt, g, h_new, u_new)
      integer, intent(in) :: n_cells
      real, intent(in) :: h(n_cells), u(n_cells), z(n_cells), dx, dt, g
      real, intent(out) :: h_new(n_cells), u_new(n_cells)
      integer :: i
      real :: F_h_i_plus_half, F_hu_i_plus_half
      real :: F_h_i_minus_half, F_hu_i_minus_half

      ! Conditions aux limites
      h_new(1) = h(1)
      u_new(1) = 0.0
      h_new(n_cells) = h(n_cells)
      u_new(n_cells) = 0.0

      ! Flux numériques
      do i = 2, n_cells - 1
        F_h_i_plus_half = 0.5 * (h(i) * u(i) + h(i+1) * u(i+1)) - &
                          0.5 * max(abs(u(i) + sqrt(g * h(i))), abs(u(i+1) + sqrt(g * h(i+1)))) * &
                          (h(i+1) - h(i))
        F_hu_i_plus_half = 0.5 * (h(i) * u(i)**2 + 0.5 * g * h(i)**2 + &
                                  h(i+1) * u(i+1)**2 + 0.5 * g * h(i+1)**2) - &
                          0.5 * max(abs(u(i) + sqrt(g * h(i))), abs(u(i+1) + sqrt(g * h(i+1)))) * &
                          (h(i+1) * u(i+1) - h(i) * u(i))

        F_h_i_minus_half = 0.5 * (h(i-1) * u(i-1) + h(i) * u(i)) - &
                          0.5 * max(abs(u(i-1) + sqrt(g * h(i-1))), abs(u(i) + sqrt(g * h(i)))) * &
                          (h(i) - h(i-1))
        F_hu_i_minus_half = 0.5 * (h(i-1) * u(i-1)**2 + 0.5 * g * h(i-1)**2 + &
                                  h(i) * u(i)**2 + 0.5 * g * h(i)**2) - &
                            0.5 * max(abs(u(i-1) + sqrt(g * h(i-1))), abs(u(i) + sqrt(g * h(i)))) * &
                            (h(i) * u(i) - h(i-1) * u(i-1))

        h_new(i) = h(i) - dt / dx * (F_h_i_plus_half - F_h_i_minus_half)
        if (h_new(i) > 1e-6) then
          u_new(i) = (h(i) * u(i) - dt / dx * (F_hu_i_plus_half - F_hu_i_minus_half)) / h_new(i)
        else
          u_new(i) = 0.0
        end if
      end do
    end subroutine standard_flux


    ! Schéma bien équilibré (code similaire au standard mais avec termes sources)
    subroutine balanced_flux(h, u, z, n_cells, dx, dt, g, h_new, u_new)
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

        ! Conditions aux limites 
        h_new(1) = h(1)
        u_new(1) = 0.0
        h_new(n_cells) = h(n_cells)
        u_new(n_cells) = 0.0

        ! Calcul des flux et des termes sources (schéma bien équilibré)
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

          ! Flux à \(i+1/2\) (bien équilibré)
          F_h_i_plus_half = 0.5 * (h_i_plus_half_G * u_i_plus_half_G + h_i_plus_half_D * u_i_plus_half_D) - &
                            0.5 * max(abs(u_i_plus_half_G) + sqrt(g * h_i_plus_half_G), &
                                      abs(u_i_plus_half_D) + sqrt(g * h_i_plus_half_D)) * &
                                      (h_i_plus_half_D - h_i_plus_half_G)

          F_hu_i_plus_half = 0.5 * (h_i_plus_half_G * u_i_plus_half_G**2 + 0.5 * g * h_i_plus_half_G**2 + &
                                    h_i_plus_half_D * u_i_plus_half_D**2 + 0.5 * g * h_i_plus_half_D**2) - &
                            0.5 * max(abs(u_i_plus_half_G) + sqrt(g * h_i_plus_half_G), &
                                      abs(u_i_plus_half_D) + sqrt(g * h_i_plus_half_D)) * &
                                      (h_i_plus_half_D * u_i_plus_half_D - h_i_plus_half_G * u_i_plus_half_G)

          ! Flux à \(i-1/2\) (bien équilibré)
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
            u_new(i) = (h(i) * u(i) - dt / dx * (F_hu_i_plus_half - F_hu_i_minus_half)) / h_new(i)
          else
            u_new(i) = 0.0
          end if
        end do
        if (abs(h_new(i) + z(i) - h(i) - z(i)) > 1e-6) then
            print *, "Warning: Loss of equilibrium at cell ", i, " at time step ", t_step
        else
            print *, "bien equilibre"
        end if


    end subroutine balanced_flux


    ! les résultats dans un fichier
  subroutine write_results(h, u, z, n_cells, t_step, scheme)
    integer, intent(in) :: n_cells, t_step
    real, intent(in) :: h(n_cells), u(n_cells), z(n_cells)
    character(len=*), intent(in) :: scheme
    integer :: i
    character(len=50) :: filename
    write(filename, '(A)') "results_" // scheme // "_t" // trim(adjustl(itoa(t_step))) // ".tex"
    open(unit=10, file=filename, status="unknown")
    write(10, *) "Cell | h |  z| h + z"
    do i = 1, n_cells
      write(10, '(I5,5F15.6)') i, (i - 1) * (2.4 / (n_cells - 1)), h(i),  z(i), h(i) + z(i)
    end do
    close(10)
  end subroutine write_results



end program comparison_schemes
