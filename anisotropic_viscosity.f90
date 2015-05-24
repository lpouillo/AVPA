!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Programme permettant de déterminer le tenseur de viscosité anisotrope de la particule de fluide, selon deux manières :
!	* milieu stratifié avec un milieu 1 défini par mu1 et K1, et un milieu 2 défini par mu2 et K2
!	* aggrégats de 1 à 2 cristaux cubiques et/ou hexagonaux défini par leur paramètre de perturbation anisotrope
!
! Laurent Pouilloux (à partir de routines créées par Jules Browaeys
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	PROGRAM anisotropic_viscosity

	USE tensmod
	USE decmod
	USE AV_variables
	USE AV_inout
	USE AV_calc

	CALL initialisation()

	SELECT(type_materiau)
		CASE(0) ! milieu stratifié
			CALL parametres_stratifie()
		CASE(1) ! milieu polycrystallin
			CALL parametres_polycristal()
	END SELECT
	
	CALL creation_matrice()

	SELECT(type_materiau)
		CASE(0)
			CALL calcul_stratifie(eta,theta,phi1,phi2)
		CASE(1)
			CALL affectation_lpo()
			CALL moyennage_aggregat()
	END SELECT

	WRITE(*,*) " "
	WRITE(*,*) "Programme terminé normalement"
	END PROGRAM anisotropic_viscosity
