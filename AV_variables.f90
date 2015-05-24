	module AV_variables

	implicit none

	! Paramètres du programme
	character(60) :: nom_run
	integer :: type_calcul, type_lpo, n_cristaux, n_dispersions, n_delta, n_angle
	! Paramètres de viscosité
	real(kind=8) :: K, mu, delta_cub, delta_1, delta_2, delta_3
	character(20) :: fichier_eta
	! Parametres de lpo
	real(kind=8) :: dispersion_perso
	character(20) :: fichier_lpo
	! Contraintes
	real(kind=8), dimension(6) :: sigma_perso, sigma_ps, sigma_ss
	
	! Tenseurs de viscosité
	real(kind=8), dimension(6,6) :: S0, S_perso, S_test, S_Reuss, S_Voigt
	
	! Angles d'Euler
	real(kind=8), dimension(:), allocatable :: phi1, theta, phi2, phi1_perso, theta_perso, phi2_perso
	
	! Variables générales
	integer :: i,j,i_delta, k1,k2,k3
	real(kind=8) :: Pi
	character(80) :: commande_shell

	! Calcul de eta=f(delta(s))
	real(kind=8), dimension(:), allocatable :: V_delta_cub, V_delta_1, V_delta_2, V_delta_3, mu0_Reuss, mu0_Voigt, &
			&Delta_Reuss_delta, Delta_Voigt_delta
	character(60) :: nom_output
	integer :: longueur_repertoire
	
	! Calcul de mu0, chi0 pour toutes les dispersions
	real(kind=8), dimension(:), allocatable :: dispersion
	real(kind=8), dimension(:,:), allocatable :: mu0_Reuss_dispersion, mu0_Voigt_dispersion, chi0_Reuss_dispersion, chi0_Voigt_dispersion
	real(kind=8), dimension(:), allocatable :: Delta_Reuss_dispersion, Delta_Voigt_dispersion
!	Tenseurs de stiffness
!	double precision, dimension(6,6) :: S0, S_temp, S_iso, eta_ij, eta_temp, S_iso_Reuss, S_iso_Voigt
!	double precision, dimension(6,6) :: S_ij_Voigt,S_ij_Reuss,S_ij_moyen
!	double precision, dimension(3,3,3,3) :: test3333
!	LPO
!	integer :: nangle
!	double precision :: k, dcostheta, dphi1, dphi2, denomi
!	double precision, dimension(:), allocatable :: k_out,Delta_total, delta_cub
!	double precision, dimension(:), allocatable :: theta, phi1, phi2
!	double precision, dimension(3,3) :: acs
!	sigma, mu0, chi0
!	double precision, dimension(6) :: sigma, sigma_ps, sigma_ss, sigma_cs
!	double precision, dimension(:), allocatable :: mu0_Voigt,mu0_Reuss,mu0_moyen,chi0_Voigt,chi0_Reuss,chi0_moyen, mu_iso_Voigt,&
!				&mu_iso_Reuss
! 	output
!	integer :: type_calcul
!	character(40) :: nom_output

	end module AV_variables
