	module AV_inout
	
	USE AV_variables
	IMPLICIT NONE
	
	CONTAINS
	
	SUBROUTINE message_accueil()
	
	! Affiche un message d'accueil précisant le numéro de version
	WRITE(*,*) "************************************"
	WRITE(*,*) " "
	WRITE(*,*) "       Anisotropic Viscosity        "
	WRITE(*,*) " 			of a"
	WRITE(*,*) " 	Polycrystalline Aggregate"
	WRITE(*,*) " "
	WRITE(*,*) "Version 0.9.8"
	WRITE(*,*) " "
	WRITE(*,*) "************************************"
	WRITE(*,*) " "
	WRITE(*,*) " "
	
	END SUBROUTINE message_accueil
	
	SUBROUTINE lecture_fichier_parametres()
	! lit les paramètres généraux du code
	WRITE(*,*) "Lecture des paramètes d'entrée"
	
	OPEN(18,file='parametres_generaux.in',status='old',action='READ')
	DO i=1,4
		READ(18,*)
	end do
	READ(18,*) nom_run
	READ(18,*) type_calcul
	READ(18,*) type_lpo
	READ(18,*) sigma_perso
	CLOSE(18)

	longueur_repertoire= len_trim(nom_run)
	commande_shell="mkdir " // nom_run(1:longueur_repertoire)
	CALL SYSTEM(commande_shell)
	commande_shell="cp parametres.in "//nom_run(1:longueur_repertoire)//"/"//nom_run(1:longueur_repertoire)//".in"
	CALL SYSTEM(commande_shell)
	if (type_calcul==0) then
		OPEN(15,file=fichier_eta,status='old',action='READ')
		READ(15,*) S_perso				
		CLOSE(15)
		commande_shell="cp "//fichier_eta//" "//nom_run(1:longueur_repertoire)//"/"//fichier_eta
		CALL SYSTEM(commande_shell)
	end if
	
	WRITE(*,*) "Terminé"	
	end subroutine lecture_fichier_parametres
	
	
	subroutine lecture_lpo(n_cristaux,fichier_lpo)
	
	integer :: n_cristaux
	character(20) :: fichier_lpo
	
	OPEN(12,file=fichier_lpo,status='old',action='READ')
	do i=1,n_cristaux
		READ(12,*) phi1_perso(i), theta_perso(i), phi2_perso(i)
	end do
	
	CLOSE(12)
	
	end subroutine lecture_lpo
	
	subroutine ecriture_eta_delta_cub(delta_out,mu_out_Reuss,mu_out_Voigt)
	integer :: i
	real(kind=8), dimension(:) :: delta_out,mu_out_Reuss,mu_out_Voigt
	
	nom_output= nom_run(1:longueur_repertoire) // "/eta_delta_cub"
	OPEN(20,FILE=nom_output,status='unknown')
	do i=0,n_delta
		WRITE(20,'(5F10.6)') delta_out(i),mu_out_Reuss(i),mu_out_Voigt(i),Delta_Reuss_delta(i),Delta_Voigt_delta(i)
	end do
	CLOSE(20)
	
	end subroutine ecriture_eta_delta_cub

	subroutine ecriture_mu0_chi0_dispersion()

	nom_output= nom_run(1:longueur_repertoire) // "/Delta_dispersion_RV"
	OPEN(12,file=nom_output,status='unknown')
	DO i=0,n_dispersions
		WRITE(12,'(5F12.8)') dispersion(i), Delta_Reuss_dispersion(i), Delta_Voigt_dispersion(i)
	END DO
	CLOSE(12)
	

	nom_output= nom_run(1:longueur_repertoire) // "/mu0_chi0_RV_ps"
	OPEN(13,file=nom_output,status='unknown')
	DO i=0,n_dispersions
		WRITE(13,'(5F12.8)') dispersion(i), mu0_Reuss_dispersion(i,1), chi0_Reuss_dispersion(i,1), mu0_Voigt_dispersion(i,1), chi0_Voigt_dispersion(i,1)
	END DO
	CLOSE(13)

	nom_output= nom_run(1:longueur_repertoire) // "/mu0_chi0_RV_ss"
	OPEN(14,file=nom_output,status='unknown')
	DO i=0,n_dispersions
		WRITE(14,'(5F12.8)') dispersion(i), mu0_Reuss_dispersion(i,2), chi0_Reuss_dispersion(i,2), mu0_Voigt_dispersion(i,2), chi0_Voigt_dispersion(i,2)
	END DO	
	CLOSE(14)
	nom_output= nom_run(1:longueur_repertoire) // "/mu0_chi0_RV_perso"
	OPEN(15,file=nom_output,status='unknown')
	DO i=0,n_dispersions
		WRITE(15,'(5F12.8)') dispersion(i), mu0_Reuss_dispersion(i,3), chi0_Reuss_dispersion(i,3), mu0_Voigt_dispersion(i,3), chi0_Voigt_dispersion(i,3)
	END DO	
	CLOSE(15)
	
	end subroutine 	ecriture_mu0_chi0_dispersion

	end module AV_inout
