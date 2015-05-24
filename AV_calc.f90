	MODULE AV_calc
	
	USE tensmod
	USE decmod
	USE AV_variables
	USE AV_inout
	
	CONTAINS
	
	SUBROUTINE initialisation()
	
	! affichage du message d'accueil
	CALL message_accueil()
	! Lecture du fichier d'entrée 
	CALL lecture_fichier_parametres
	! Affectation du tenseur du monocrytal
	
	END SUBROUTINE initialisation()
	
	SUBROUTINE parametres_stratifie()
	
	OPEN(10,file="parametres_stratifie.in")
	READ(10,*) mu_star, K_star
	READ(10,*) compressibility
	READ(10,*) Phi
	READ(10,*) n_x,n_y,n_z
	CLOSE(10)
	
	END SUBROUTINE parametres_stratifie
	
	SUBROUTINE parametres_polycristal
	
	OPEN(10,file="parametres_polycristal.in")
	READ(10,*) type_symetrie
	SELECT(type_symetrie)
		CASE(0) ! Cas cubique
		READ(10,*) K, mu, delta_cub
		
		CASE(1) ! Cas hexagonal
		READ(10,*) K, mu, delta_hex1, delta_hex2, delta_hex3
		
		CASE(2) ! Milieu donné
		DO i=1,6
			DO j=1,6
				READ(10,*) eta(i,j)
			END DO
		END DO
	END SELECT
	
	CLOSE(10)
	
	
	END SUBROUTINE parametres_polycrystal
	
	if (type_lpo==0) then
		CALL lecture_lpo(n_cristaux,fichier_lpo)
	end if
	n_angle=abs(real(n_cristaux)**(1./3.))+1
	
	ALLOCATE(phi1_perso(n_cristaux), theta_perso(n_cristaux), phi2_perso(n_cristaux))
	ALLOCATE(V_delta_cub(n_delta+1),mu0_Reuss(n_delta+1),mu0_Voigt(n_delta+1),phi1(n_angle), theta(n_angle), phi2(n_angle))
    ALLOCATE(V_delta_1(n_delta+1), V_delta_2(n_delta+1), V_delta_3(n_delta+1))
    ALLOCATE(Delta_Reuss_delta(n_delta+1),Delta_Voigt_delta(n_delta+1))
    ALLOCATE(dispersion(n_dispersions+1),mu0_Reuss_dispersion(n_dispersions+1,3),mu0_Voigt_dispersion(n_dispersions+1,3))
    ALLOCATE(chi0_Reuss_dispersion(n_dispersions+1,3),chi0_Voigt_dispersion(n_dispersions+1,3))
    ALLOCATE(Delta_Reuss_dispersion(n_dispersions+1),Delta_Voigt_dispersion(n_dispersions+1))

	sigma_ps=0.
	sigma_ss=0.
	sigma_ps(1)=-0.5 ; sigma_ps(2)=-0.5; sigma_ps(3)=1;
	sigma_ss(6)=sqrt(2.)
	
	Pi=acos(-1.0)
	end subroutine initialisation
	
	subroutine creation_S_perturbation(S,calcul,delta,delta1,delta2,delta3)
	
	integer :: calcul,i
	real(kind=8) :: delta, delta1, delta2, delta3
	real(kind=8),dimension(6,6) :: S
	SELECT CASE(type_calcul)
		CASE(0)
		! Tenseur personnalisé
		OPEN(10,file=fichier_eta)
		do i=1,6
			READ(10,*) S(i,1),S(i,2),S(i,3),S(i,4),S(i,5),S(i,6)
		end do
		CLOSE(10)
		CASE(1)
		! Tenseur cubique
		S=0
		do i=1,3
			S(i,i)=1/(2*mu+delta*3/5)
			S(i+3,i+3)=1/(2*mu-delta*2/5)
		end do
		CASE(2)
		S=0
		write(*,*) "cas hexagonal non implémenté"
	END SELECT
	end subroutine creation_S_perturbation
	
	
	
	subroutine calcul_eta_delta()
	
	CALL affectation_angles(1D0)
	SELECT CASE(type_calcul)
		CASE(1)
		V_delta_cub=0
		do i_delta=0,n_delta
			call SYSTEM("clear")
			write(*,*) "Calcul de mu isotrope en fonction de delta_cub pour un aggrégat uniformément distribué"
			write(*,'(F5.1,A2)') 100*real(i_delta)/real(n_delta),' %'
			V_delta_cub(i_delta)=-10./3.+(5.+10./3.)*i_delta/n_delta
			CALL creation_S_perturbation(S_test,type_calcul,V_delta_cub(i_delta),0D0,0D0,0D0)
			call moyennage_reuss(S_test,S_Reuss,phi1,theta,phi2)
			call DEC_CONT(CKEL(S_Reuss))
			call DEC_ISO(CKEL(S_Reuss))
			Delta_Reuss_delta(i_delta)=PERCA
			call moyennage_voigt(S_test,S_Voigt,phi1,theta,phi2)
			call DEC_CONT(CKEL(S_Voigt))
			call DEC_ISO(CKEL(S_Voigt))
			Delta_Voigt_delta(i_delta)=PERCA
			call calcul_viscosite(S_Reuss,sigma_ps,mu0_Reuss(i_delta))
			call calcul_viscosite(S_Reuss,sigma_ss,mu0_Reuss(i_delta))
			call calcul_viscosite(S_Voigt,sigma_ps,mu0_Voigt(i_delta))
			call calcul_viscosite(S_Voigt,sigma_ss,mu0_Voigt(i_delta))
		end do
		call ecriture_eta_delta_cub(V_delta_cub,mu0_Reuss,mu0_Voigt)		
	
		CASE(2)
		WRITE(*,*) "Calcul du tenseur de viscosité pour un aggrégat de glace"
	END SELECT


	end subroutine calcul_eta_delta

	SUBROUTINE affectation_angles(d)
	! Dispersion d
	integer :: j
	real(kind=8) :: d
	
	DO j=1,n_angle
     		phi1(j) = (REAL(j)/REAL(n_angle))*ACOS(-1d0)*d
			theta(j) = ACOS(1d0 - (REAL(j))/REAL(n_angle)*2d0*d)
	! Ligne PAS ISOTROPE.. theta(i,j) = (REAL(i))/REAL(nangle)*ACOS(-1d0)*k_out(j)
      		phi2(j) = (REAL(j))/REAL(n_angle)*ACOS(-1d0)*d
	END DO

	END SUBROUTINE affectation_angles
	
	subroutine moyennage_reuss(S0,S_ij_Reuss,phi1,theta,phi2)
	
	real(kind=8), dimension(6,6) :: S0, S_ij_Reuss, S_temp
	real(kind=8), dimension(:) :: phi1,theta,phi2
    	real(kind=8), dimension(3,3,3,3) :: test3333
	real(kind=8), dimension(3,3) :: acs

	S_ij_Reuss(:,:)=0
	do k1=1,n_angle
		do k2=1,n_angle
			do k3=1,n_angle
				acs=ROTMAT(180/Pi*phi1(k1),180/Pi*theta(k2),180/Pi*phi2(k3))
				test3333=T4C(CKEL(S0))
				test3333=RT4(test3333,acs)
				S_temp=CT4(test3333)
				S_temp=KELC(S_temp)
				do i=1,6 
					do j=1,6
						S_ij_Reuss(i,j)=S_ij_Reuss(i,j)+S_temp(i,j)
					end do
				end do
			end do
		end do
	end do
	S_ij_Reuss(:,:)=S_ij_Reuss(:,:)/n_angle**3
	
	end	subroutine moyennage_reuss

	subroutine moyennage_voigt(S0,S_ij_Voigt,phi1,theta,phi2)
	
	real(kind=8), dimension(6,6) :: S0, S_ij_Voigt, S_temp, eta_ij, eta_temp
	real(kind=8), dimension(:) :: phi1,theta,phi2
   	real(kind=8), dimension(3,3,3,3) :: test3333	
    	real(kind=8), dimension(3,3) :: acs
        
	eta_ij=0
	S_ij_Voigt(:,:)=0
	do k1=1,n_angle
		do k2=1,n_angle
			do k3=1,n_angle
				acs=ROTMAT(180/Pi*phi1(k1),180/Pi*theta(k2),180/Pi*phi2(k3))
				test3333=T4C(CKEL(S0))
				test3333=RT4(test3333,acs)
				S_temp=CT4(test3333)
				S_temp=KELC(S_temp)
				CALL LUINV(S_temp,6,eta_temp)
				do i=1,6 
					do j=1,6
						eta_ij(i,j)=eta_ij(i,j)+eta_temp(i,j)
					end do
				end do
			end do
		end do
	end do
	eta_ij=eta_ij/n_angle**3
	CALL LUINV(eta_ij,6,S_ij_Voigt(:,:))
	
	end	subroutine moyennage_voigt
	
	subroutine calcul_viscosite(S_ij,sigma_applique,mu0)
	real(kind=8), dimension(6,6) :: S_ij
	real(kind=8), dimension(6) :: sigma_applique, epsilon_vrai,epsilon_iso
	real(kind=8) :: s, denomi, mu0
	! Calcul de la déformation réelle
	epsilon_vrai=matmul(S_ij,sigma_applique)
	! Détermination de la meilleure viscosité
	s=0.
	denomi=0.
	do i=1,6
		s=s+epsilon_vrai(i)*sigma_applique(i)
		denomi=denomi+sigma_applique(i)*sigma_applique(i)
	end do
	s=s/denomi
	mu0=1/(2*s)


	end subroutine calcul_viscosite
	
	
	subroutine calcul_mu0_chi0(n_dispersion)
	integer :: i_dispersion
	
	dispersion=0
	do i_dispersion=0,n_dispersion
		dispersion(i_dispersion)=real(i_dispersion)/real(n_dispersion)
		call affectation_angles(dispersion(i_dispersion))
		call moyennage_reuss(S0,S_Reuss,phi1,theta,phi2)
		call DEC_CONT(CKEL(S_Reuss))
		call DEC_ISO(CKEL(S_Reuss))
		Delta_Reuss_dispersion(i_dispersion)=PERCA
		call moyennage_voigt(S0,S_Voigt,phi1,theta,phi2)
		call DEC_CONT(CKEL(S_Voigt))
		call DEC_ISO(CKEL(S_Voigt))
		Delta_Voigt_dispersion(i_dispersion)=PERCA
		call calcul_viscosite_erreur(S_Reuss,sigma_ps,mu0_Reuss_dispersion(i_dispersion,1),chi0_Reuss_dispersion(i_dispersion,1))
		call calcul_viscosite_erreur(S_Reuss,sigma_ss,mu0_Reuss_dispersion(i_dispersion,2),chi0_Reuss_dispersion(i_dispersion,2))
		call calcul_viscosite_erreur(S_Reuss,sigma_perso,mu0_Reuss_dispersion(i_dispersion,3),chi0_Reuss_dispersion(i_dispersion,3))
		call calcul_viscosite_erreur(S_Voigt,sigma_ps,mu0_Voigt_dispersion(i_dispersion,1),chi0_Voigt_dispersion(i_dispersion,1))
		call calcul_viscosite_erreur(S_Voigt,sigma_ss,mu0_Voigt_dispersion(i_dispersion,2),chi0_Voigt_dispersion(i_dispersion,2))
		call calcul_viscosite_erreur(S_Voigt,sigma_perso,mu0_Voigt_dispersion(i_dispersion,3),chi0_Voigt_dispersion(i_dispersion,3))
		call SYSTEM("clear")
		write(*,*) "Calcul de mu0 et chi0 en fonction de la dispersion"
		write(*,'(F5.1,A2,F10.3)') 100*real(i_dispersion)/real(n_dispersion), ' %'
	end do
		
	
	end subroutine calcul_mu0_chi0
	
	
	subroutine calcul_viscosite_erreur(S_ij,sigma_applique,mu0,chi0)
	real(kind=8), dimension(6,6) :: S_ij, S_iso
	real(kind=8), dimension(6) :: sigma_applique, epsilon_vrai,epsilon_iso
	real(kind=8) :: s, denomi, chi0, mu0
	
	call calcul_viscosite(S_ij,sigma_applique,mu0)
! Création du tenseur isotrope correspondant et calcul du tenseur de taux de déformation correspondant
	s=1/(2*mu0)
	S_iso=0
	do i=1,3
		S_iso(i,i)=s
		S_iso(i+3,i+3)=s
	end do
	epsilon_iso=matmul(S_iso,sigma_applique)
	epsilon_vrai=matmul(S_ij,sigma_applique)
	!calcul du chi0
	chi0 = 0.
	denomi=0.
	do i=1,6
		chi0 = chi0 + (epsilon_vrai(i)-epsilon_iso(i))**2
		denomi= denomi + epsilon_vrai(i)**2
	end do
	chi0 = sqrt(chi0/denomi)*100
	
	end subroutine calcul_viscosite_erreur
	end module AV_calc
