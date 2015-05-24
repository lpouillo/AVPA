	module AV_routines

	use AV_variables
	use tensmod
	use decmod

	implicit none

	contains
	
	SUBROUTINE initialisation()
!!! Subroutine permettant l'initialisation du programme
! * Affichage d'un message d'accueil
! * Lecture du fichier d'entrée
! * Vérification des paramètres

	WRITE(*,*) "************************************"
	WRITE(*,*) " "
	WRITE(*,*) "       Anisotropic Viscosity        "
	WRITE(*,*) " 			of a"
	WRITE(*,*) " 	Polycrystalline Aggregate"
	WRITE(*,*) " "
	WRITE(*,*) "Version 0.0.9"
	WRITE(*,*) " "
	WRITE(*,*) "************************************"
	WRITE(*,*) " "
	WRITE(*,*) " "

	! Lecture des paramètres
	WRITE(*,*) "Lecture du fichier d'entrée ..."
	OPEN(18,FILE='parametres.init')
	! Parametres du tenseur du monocristal cubique
	READ(18,*) mu, delta
	! nangle^3=nombre de cristaux, nk = résolution en dispersion pour les fichiers de données
	READ(18,*) nangle, nk, n_delta
	! Tenseur des contraintes appliques
	READ(18,*) sigma
    ! Lecture du nom de fichier de sortie
    READ(18,'(a15)') nom_output

    
	CLOSE(18);
	WRITE(*,*) "Terminé !"

	! Affichage des parmètres pour vérification
	WRITE(*,*) " "
	WRITE(*,'(A9,6F7.3,A1)') "sigma =(", sigma,")"
	WRITE(*,'(A9,F6.2)') "mu =",mu
	WRITE(*,'(A9,F6.2)') "delta =",delta
	WRITE(*,*) " "
	WRITE(*,'(A9,I5)') "nk =",nk
	WRITE(*,'(A9,I5)') "nangle =",nangle	
	WRITE(*,'(A9,I5)') "n_delta =",n_delta	
	WRITE(*,*) " "
        WRITE(*,*)
	WRITE(*,*) "Are this parameter correct ?" 
	WRITE(*,*) "Press Enter to continue, Ctrl+C to Abort"
	READ(*,*)
	
	! Création des tenseurs de contraintes standards
	sigma_ps=0.
	sigma_ss=0.
	sigma_cs=0.
	sigma_ps(1)=-0.5 ; sigma_ps(2)=-0.5; sigma_ps(3)=1;
	sigma_ss(6)=sqrt(2.)
	sigma_cs(1)=1 ; sigma_ps(2)=-1; sigma_ps(6)=sqrt(2.);
	
	! Calcul du tenseur du monocristal
	CALL calcul_S0(mu,delta)
	
	! Allocation des variables		
	ALLOCATE(mu0_Voigt(0:nk), chi0_Voigt(0:nk),mu0_Reuss(0:nk), chi0_Reuss(0:nk),mu0_moyen(0:nk), chi0_moyen(0:nk))
	ALLOCATE(theta(nangle),phi1(nangle),phi2(nangle),k_out(0:nk))
	ALLOCATE(delta_cub(n_delta),mu_iso_Voigt(n_delta),mu_iso_Reuss(n_delta))
	Pi=acos(-1.0)
	
	

	WRITE(*,*) "Initialisation terminée."
	WRITE(*,*) " "
	
	END SUBROUTINE initialisation

	SUBROUTINE calcul_S0(mu,delta)
	real(kind=8) :: mu,delta
	
	! Calcul du tenseur cubique du monocrystal
	S0=0
	do i=1,3
		S0(i,i)=1/(2*mu+delta*3/5)
		S0(i+3,i+3)=1/(2*mu-delta*2/5)
	end do
	!WRITE(*,'(6F12.4)') S0
	WRITE(*,*) "Terminé !"	
	WRITE(*,*) " "


	
	END SUBROUTINE calcul_S0


	
	

	SUBROUTINE affectation_angles(d)
! Dispersion d
	real(kind=8) :: d
	
	DO j=1,nangle
     		phi1(j) = (REAL(j)/REAL(nangle))*ACOS(-1d0)*d
			theta(j) = ACOS(1d0 - (REAL(j))/REAL(nangle)*2d0*d)
! Ligne PAS ISOTROPE.. theta(i,j) = (REAL(i))/REAL(nangle)*ACOS(-1d0)*k_out(j)
      		phi2(j) = (REAL(j))/REAL(nangle)*ACOS(-1d0)*d
	END DO
	
	END SUBROUTINE affectation_angles

	subroutine moyennage_reuss()
	
	S_ij_Reuss(:,:)=0
	do k1=1,nangle
		do k2=1,nangle
			do k3=1,nangle
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
	S_ij_Reuss(:,:)=S_ij_Reuss(:,:)/nangle**3
	
	end	subroutine moyennage_reuss

	subroutine moyennage_voigt()
	
	eta_ij=0
	S_ij_Voigt(:,:)=0
	do k1=1,nangle
		do k2=1,nangle
			do k3=1,nangle
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
	eta_ij=eta_ij/nangle**3
	CALL LUINV(eta_ij,6,S_ij_Voigt(:,:))
	
	end	subroutine moyennage_voigt

	subroutine calcul_viscosite(S_ij,sigma_applique,mu0)
	real(kind=8), dimension(6,6) :: S_ij
	real(kind=8), dimension(6) :: sigma_applique, epsilon_vrai,epsilon_iso
	real(kind=8) :: s, denomi, mu0
	! Calcul de la déformation réelle
	epsilon_vrai=matmul(S_ij,sigma_applique)
	write(*,'(6F10.6)')epsilon_vrai
	! Détermination de la meilleure viscosité
	s=0.
	denomi=0.
	do i=1,6
		s=s+epsilon_vrai(i)*sigma_applique(i)
		!write(*,*) s
		denomi=denomi+sigma(i)*sigma_applique(i)
		!write(*,*) denomi
	end do
	s=s/denomi
	mu0=1/(2*s)
	!write(*,*) mu0
	write(*,*) ' '

	end subroutine calcul_viscosite

	subroutine calcul_erreur(S_ij,sigma_applique,mu0,chi0)
	real(kind=8), dimension(6,6) :: S_ij
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
	
!calcul du chi0
	chi0 = 0.
	denomi=0.
	do i=1,6
		chi0 = chi0 + (epsilon_vrai(i)-epsilon_iso(i))**2
		denomi= denomi + epsilon_vrai(i)**2
	end do
	chi0 = sqrt(chi0/denomi)*100
	end subroutine calcul_erreur


	subroutine calcul_direct
	
	do a=0,nk
		WRITE(*,'(F7.3,A4)') 100*real(a)/real(nk),"/100"
		k_out(a)=real(a)/real(nk)
		call affectation_angles(k_out(a))
		! Moyennage de Reuss
		call moyennage_reuss()
		call calcul_erreur(S_ij_Reuss,sigma,mu0_Reuss(a),chi0_Reuss(a))
		! Moyennage de Voigt
		call moyennage_voigt()
		call calcul_erreur(S_ij_Voigt,sigma,mu0_Voigt(a),chi0_Voigt(a))
		! Tenseur moyen
		S_ij_moyen=0.5*(S_ij_Reuss+S_ij_Voigt)
		call calcul_erreur(S_ij_moyen,sigma,mu0_moyen(a),chi0_moyen(a))		
		
	end do
	write(*,*)'termine !'
	
	end subroutine calcul_direct

	subroutine calcul_isotrope()
	
	delta_cub=0
	do a=0,n_delta
		! Création du tenseur isotrope de référence
		delta_cub(a)=-10./3.+a*(5.+10./3.)/n_delta
		write(*,*) a,	delta_cub(a)
		call calcul_S0(mu,delta_cub(a))
		k=1
		call affectation_angles(k)
		call moyennage_reuss()
		call moyennage_voigt()
		S_iso_Reuss=S_ij_Reuss
		S_iso_Voigt=S_ij_Voigt
		CALL calcul_viscosite(S_iso_Reuss,sigma_ss,mu_iso_Reuss(a))
		CALL calcul_viscosite(S_iso_Voigt,sigma_ss,mu_iso_Voigt(a))
		WRITE(*,*) 'Simple shear : ', mu_iso_Reuss(a), mu_iso_Voigt(a)
		CALL calcul_viscosite(S_iso_Voigt,sigma_ps,mu_iso_Voigt(a))
		CALL calcul_viscosite(S_iso_Reuss,sigma_ps,mu_iso_Reuss(a))
		WRITE(*,*) 'Pure Shear : ',mu_iso_Reuss(a), mu_iso_Voigt(a)
		CALL calcul_viscosite(S_iso_Voigt,sigma_cs,mu_iso_Voigt(a))
		CALL calcul_viscosite(S_iso_Reuss,sigma_cs,mu_iso_Reuss(a))
		WRITE(*,*) 'Combined Stress : ',mu_iso_Reuss(a), mu_iso_Voigt(a)
	end do
	
	end subroutine calcul_isotrope
		
	SUBROUTINE ecriture_output()
!!! Ecriture des fichier de sortie nommé en fonction des paramètres du calcul
!sauvegardé dans le répertoire Resultat
	WRITE(*,*) "Ecriture du fichier de sortie"
        if (nom_output=="default") then
        	if(delta.lt.0) then
	        	WRITE(nom_output,'(A23,F3.1,A3,F3.1,A2,F4.1)') 'Results/mu0_chi0_VRM_k_',k,'_mu',mu,'_d',delta
        	else
	        	WRITE(nom_output,'(A23,F3.1,A3,F3.1,A2,F3.1)') 'Results/mu0_chi0_VRM_k_',k,'_mu',mu,'_d',abs(delta)
        	end if
        else
                WRITE(nom_output,'(A8,A15)')'Results/',nom_output
        end if
        WRITE(*,*) nom_output
       	OPEN(15,file=nom_output,action="write",status="unknown")
	do i=0,nk
		write(15,'(7F12.5)') k_out(i),mu0_Reuss(i),chi0_Reuss(i),mu0_Voigt(i),chi0_Voigt(i),mu0_moyen(i),chi0_moyen(i)
	end do
	CLOSE(15)
	WRITE(*,*)"Terminé !"

	END SUBROUTINE ecriture_output

	end module AV_routines
