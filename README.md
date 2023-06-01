# UVLM-par-rotor
Unsteady vortex lattice method for helicopter rotor, parallel implementation Matlab-C++
/!\ Avant tout /!\ : 

Sur son propre ordinateur :	- Accéder au fichier "run_UVLM_own_computer.sh"
				- Changer la variable MATLABHOME par le chemin d'accès au dossier MATLAB	
				  regroupant les dossiers "appdata", "bin", "cefclient", "derived", "examples", "extern", "help", "interprocess", "java", "licences", etc..
				  
				  [Commandes d'execution du code] :
				  
				  	taper "bash run_UVLM_own_computer.sh <nombre_de_threads_alloués>" dans le terminal local de votre ordinateur 
				  	
Sur la grappe cedar       :	- Accéder au fichier "run_UVLM_cedar.sh" sur votre ordinateur
				- (Normalement, cette étape est déja effectuée) Changer la variable MATLABHOME par le chemin d'accès au dossier MATLAB	
				  regroupant les dossiers "appdata", "bin", "cefclient", "derived", "examples", "extern", "help", "interprocess", "java", "licences", etc..
				- Définir le nombre "n" de coeurs utilisés dans le fichier "run_UVLM_cedar.sh" :
					- #SBATCH --job-name=n
					- export OMP_NUM_THREADS=n
					- ./UVLM_2_Blade_Rotor_Axial n
 
				- Copier les fichiers du dossier local contenant le code (/home/guillaume/Bureau/UVLM_2BladeAxial_parallel) 
				  dans un dossier créé sur la grappe Cedar (gui13@cedar.computecanada.ca:projects/def-morency/gui13/UVLM_2BladeAxial_parallel) : 
				  	
				  	taper par exemple "scp /home/guillaume/Bureau/UVLM_2BladeAxial_parallel/* gui13@cedar.computecanada.ca:projects/def-morency/gui13/UVLM_2BladeAxial_parallel/" 
				  	dans le terminal de votre ordinateur
				 
				- Accéder à la grappe : 
				
					taper par exemple "ssh -Y <nom_utilisateur>@cedar.calculcanada.ca" dans le terminal de votre ordinateur
					
				- Accéder à votre dossier sur la grappe Cedar
				
				
				  [Commandes d'execution du code] :
				  
				  	taper "sbatch run_UVLM_cedar.sh" dans le terminal sur la grappe Cedar
				
				
				- Pour accéder à vos fichiers résultats, copier les depuis votre ordinateur local :
					
					taper par exemple "scp gui13@cedar.computecanada.ca:projects/def-morency/gui13/UVLM_2BladeAxial_parallel/* /home/guillaume/Bureau/UVLM_2BladeAxial_parallel/" 
					dans le terminal local de votre ordinateur
				


/i\ Il est très important de passer par ces commandes car elles comportent le chemin d'accès vers les librairies -lmx, -leng nécessaires au bon fonctionnement du couplage Matlab/C++ /i\

=====================================================================================================================


[Description des différents fichiers] :

[main.cpp]     	      : Fichier contenant le programme principal C++/Matlab et l'ensemble des fonctions nécessaires aux calculs des nouvelles positions des panneaux de vorticité à chaque 					itération.

[Documentation.html]         : Lien html vers la documentation du code et des formules mathématiques associées aux fonctions.

Sections résiduelles du code UVLM MATLAB :

[UVLM_2BladeAxial_Init.m]
[UVLM_2BladeAxial_Boucle1.m]
[UVLM_2BladeAxial_Boucle2.m]
[UVLM_2BladeAxial_avantParallel.m]
[UVLM_2BladeAxial_Parallel.m]
[UVLM_2BladeAxial_apresParallel.m]
[UVLM_2BladeAxial_Save.m]

=====================================================================================================================

[Sortie et commandes d'execution] :

[Thrust_Data.csv]            : [Azimuth]        : Tableau de double qui caractérise l'angle parcouru pour chaque itération - En degrés
				[NREV]           : Tableau de double qui caractérise le nombre de révolutions pour chaque itération
				[CT]             : Tableau caractérisant le coefficient de poussée adimensionnel à chaque itération.
				[CQ]             : Tableau caractérisant le coefficient de couple adimensionnel à chaque itération.
				[CL]             : Tableau caractérisant le coefficient de portance adimensionnel à chaque itération.
				[FM]             : Symbole de mérite à chaque itération.
				[IterationTime]  : Tableau de double qui caractérise le temps d'exécution de chaque itération.
				
[Initial_Data.csv]           : [C]              : Longueur de la corde.
				[B]              : Rayon total de la pale, y compris l'espacement entre le centre de rotation et la pale.
				[RootY]          : Plus grand espacement entre la pale et le centre de rotation dans la direction y.
				[IB]             : Nombre de panneaux au niveau d'une pale suivant la direction radial. 
				[JB]             : Nombre de panneaux au niveau d'une pale suivant le sens de la corde.
				[RPM]            : Vitesse angulaire (Nombre de révolutions par minute).
				[ALFA*180/pi]    : La pente (angle entre la vitesse de la pale et l’horizontale) 
				[AirfoilShape]   : Choix du profil des pales : 1 pour NACA 0012, 2 pour NACA 4412, 3 pour NACA63421 et 4 pour E374.
				[DTETA]          : Angle parcouru par les pales, pendant une itération de temps.
				[NREVSLOWSTART]  : Pendant ce nombre de révolutions du rotor, la vitesse angulaire du rotor augmente jusqu'à atteindre la vitesse angulaire définie par RPM.
				[SpeedofSound]   : Vitesse du son (m/s).
				[Solidity]       : Paramètre lié à la conception du système de rotor.
				[CH]             : Influence du sol (Ground clearence), définie comme la somme de B et RootY à un coefficient multicatif près.
				[ClimbRatio]     : Coefficient multiplicatif pour calculer la vitesse vertical de l'helicoptère par rapport à la vitesse des extrémités des pales.
				
[ReY.csv]                    : Tableau caractérisant le nombre de Reynolds de la section radiale d'une pale à chaque itération.

[CLY_Values.csv]             : Tableau caractérisant le coefficient de portance de la section radiale d'une pale, normalisé par toute la surface de la pale, à chaque itération.

[CLYJ_Values.csv]            : Tableau caractérisant le coefficient de portance de la section radiale d'une pale, normalisé par la surface de cette section radiale uniquement, à chaque itération.

[Wake_Shape_i.csv]           : Tableau caractérisant la localisation selon la direction i=1,2 ou 3 (x, y ou z) des panneaux de sillage à la dernière itération.

[Wake_ShapeA_i.csv]          : Tableau caractérisant la localisation selon la direction i=1,2 ou 3 (x, y ou z) des panneaux de sillage à la dernière itération pour la seconde pale.

[buffer.out]                 : Fichier contenant ce qui a été écrit dans le buffer (uniquement lors de l'utilisation de la grappe Cedar).
