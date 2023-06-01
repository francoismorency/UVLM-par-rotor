#include "stdio.h"
#include "string.h"
#include <stdlib.h>
#include <math.h>

#include "engine.h"
#include <omp.h>

/**
* @file main.cpp
*	@brief Parallélisation d'un code UVLM à l'aide d'un couplage MATLAB/C++
*
* Ce code fait partie intégrante d'une méthode de résolution UVLM, à priori, entièrement écrit en MATLAB. Cependant, ce code était de plus en plus lent d'une itération de temps à l'autre. La section la plus lente a donc été réécrite en C++, dans le but d'être parallélisée. Le fichier main.cpp contient la section parallélisable à coupler avec les autres sections MATLAB du code initiale. La section écrite en C++ a pour but de calculer le déplacement de chaque panneau du sillage en aval des pales à l'aide d'un schéma explicite.
*
*  @image html schemaaile.jpg "Figure 1. Répartition des panneaux au niveau d’une aile et au niveau de son sillage." width=870cm height=421cm
*
* \n
* Le calcul des nouvelles positions des panneaux du sillage est de plus en plus couteux car à chaque itération, de nouveaux panneaux sont créés (Figure 2). Ainsi, au temps \f$t_{0}\f$, il y a trois rangées de coins de panneau de taille JB (\f$t_{0}\f$, \f$t_{0}-DT\f$, \f$t_{0}-2 \cdot DT\f$) qui au temps \f$t_{0}+DT\f$ se déforment et à ces rangées s'ajoutent la rangée (\f$t_{0}+DT\f$).
*
* @image html helicotemps.jpg "Figure 2. Répartition des panneaux pour un rotor à deux pales et ajout de panneaux pendant une itération." width=850cm height=621cm
*
* \n
*
* Ainsi, pour réduire le temps d'exécution, le calcul du déplacement des panneaux est partionné sur les différents coeurs (threads) disponibles comme sur l'exemple suivant (Figure 3) :
*
* @image html parallelisation.jpg "Figure 3. Exemple de répartition des panneaux au niveau du sillage sur 3 coeurs." width=920cm height=480cm
*
* \n
* La parallélisation du code a été effectuée grâce à Openmp.
*
*	@author Guillaume DAMOUR
* @version 1.0
* @date 31/08/2022
*
*/
#define  BUFSIZE 256

/**
* @brief Blade_1_QF est le tableau contenant l'ensemble des coins des panneaux au niveau de l'une des ailes.
*
* \n
* Le partage de variables entre MATLAB/C++ contraint l'utilisateur à un agencement des données particulier (Blade_1_QF est un tableau à 3 dimensions dans MATLAB et à 1 dimension en C++)
* @param  int i              : 1ère dimension de Blade_1_QF : (IB1)
* @param  int j              : 2ème dimension de Blade_1_QF : (JB1)
* @param  int k               : 3ème dimension de Blade_1_QF : (3)
*/
#define  Blade_1_QF(i,j,k) Blade_1_QF[(k-1)*(IB+1)*(JB+1) + (i-1) + (j-1)*(IB+1)]
/**
* @brief Blade_2_QF est le tableau contenant l'ensemble des coins des panneaux au niveau de l'autre aile.
*
* \n
* Le partage de variables entre MATLAB/C++ contraint l'utilisateur à un agencement des données particulier (Blade_2_QF est un tableau à 3 dimensions dans MATLAB et à 1 dimension en C++)
* @param  int i              : 1ère dimension de Blade_2_QF : (IB1)
* @param  int j              : 2ème dimension de Blade_2_QF : (JB1)
* @param  int k               : 3ème dimension de Blade_2_QF : (3)
*/
#define  Blade_2_QF(i,j,k) Blade_2_QF[(k-1)*(IB+1)*(JB+1) + (i-1) + (j-1)*(IB+1)]
/**
* @brief QW est le tableau contenant l'ensemble des coins des panneaux liés au sillage derrière une des ailes.
*
* \n
* Le partage de variables entre MATLAB/C++ contraint l'utilisateur à un agencement des données particulier (QW est un tableau à 3 dimensions dans MATLAB et à 1 dimension en C++)
* @param  int i               : 1ère dimension de QW (NSTEPS)
* @param  int j               : 2ème dimension de QW (JB1)
* @param  int k               : 3ème dimension de QW (3)
*/
#define  QW(i,j,k) QW[(k-1)*NSTEPS*(JB+1) + (i-1) + (j-1)*NSTEPS]
/**
* @brief QW est le tableau contenant l'ensemble des coins des panneaux liés au sillage derrière l'autre aile.
*
* \n
* Le partage de variables entre MATLAB/C++ contraint l'utilisateur à un agencement des données particulier (QWA est un tableau à 3 dimensions dans MATLAB et à 1 dimension en C++)
* @param  int i               : 1ère dimension de QWA (NSTEPS)
* @param  int j               : 2ème dimension de QWA (JB1)
* @param  int k               : 3ème dimension de QWA (3)
*/
#define  QWA(i,j,k) QWA[(k-1)*NSTEPS*(JB+1) + (i-1) + (j-1)*NSTEPS]
/**
* @brief QWnplusun est le tableau contenant l'ensemble des coins des panneaux liés au sillage derrière une des ailes après un pas de temps, après déformations.
*
* \n
* Le partage de variables entre MATLAB/C++ contraint l'utilisateur à un agencement des données particulier (QWA est un tableau à 3 dimensions dans MATLAB et à 1 dimension en C++)
* @param  int i               : 1ère dimension de QWA (NSTEPS)
* @param  int j               : 2ème dimension de QWA (JB1)
* @param  int k               : 3ème dimension de QWA (3)
*/
#define  QWnplusun(i,j,k) QWnplusun[(k-1)*NSTEPS*(JB+1) + (i-1) + (j-1)*NSTEPS]
/**
* @brief VORTIC est le tableau contenant l'ensemble des forces \f$ \Gamma \f$ des segments des panneaux liés au sillage.
*
* \n
*  Le partage de variables entre MATLAB/C++ contraint l'utilisateur à un agencement des données particulier (VORTIC est un tableau à 2 dimensions dans MATLAB et à 1 dimension en C++)
* @param  int i               : 1ère dimension de VORTIC (NSTEPS+1)
* @param  int j              : 2ème dimension de VORTIC (JB)
*/
#define  VORTIC(i,j) VORTIC[(i-1) + (j-1)*(NSTEPS+1)]
/**
* @brief GAMA est le tableau contenant l'ensemble des forces \f$ \Gamma \f$ des segments des panneaux au niveau de l'aile.
*
* \n
* Le partage de variables entre MATLAB/C++ contraint l'utilisateur à un agencement des données particulier (GAMA est un tableau à 2 dimensions dans MATLAB et à 1 dimension en C++)
* @param  int i               : 1ère dimension de GAMA : (IB]
* @param  int j               : 2ème dimension de GAMA : (JB)
*/
#define  GAMA(i,j) GAMA[(i-1) + (j-1)*IB]
/**
*	@brief Vitesse est le tableau contenant l'ensemble des vitesses des coins des panneaux lié au sillage.
*
* \n
* Le partage de variables entre MATLAB/C++ contraint l'utilisateur à un agencement des données particulier (QW est un tableau à 3 dimensions dans MATLAB et à 1 dimension en C++)
* @param  int i               : 1ère dimension de QW (NSTEPS)
* @param  int j               : 2ème dimension de QW (JB1)
* @param  int k               : 3ème dimension de QW (3)
*/
#define  Vitesse(i,j,k) Vitesse[(k-1)*NSTEPS*(JB+1) + (i-1) + (j-1)*NSTEPS]

/**
* @brief Nombre d'itérations en temps
*/
int NSTEPS;
/**
* @brief Nombre d'itérations en temps
*/
double DT;
/**
* @brief Nombre de panneaux dans le sens de la corde
*/
int IB;
/**
* @brief Nombre de panneaux dans le sens de l'envergure (SPANWISE)
*/
int JB;
/**
* @brief NW est le nombre d'segments se déformant d'une itération à l'autre (ce paramètre est inutilisé dans la version original)
*/
int NW;
/**
* @brief Garde au sol (GROUND CLEARANCE)
*/
int CH;

void F_Wake_Rollup(double* U, double* V, double* W, double X, double Y, double Z, int IT, double* QW, double* QWA, double* VORTIC, double* GAMA, int CH, int IB, int JB, double T, double* Blade_1_QF, double* Blade_2_QF, double* Wake_Age);
void F_Wake_Influence(double* U, double* V, double* W, double X, double Y, double Z, int IT,  double* VORTIC, double* QW, int JB, double* Time);
void F_Wing_Influence_on_Wake(double* U, double* V, double* W, double X, double Y, double Z,  double* GAMA, double T, int IB, int JB, double* Blade_1_QF);
void F_Wake_Influence_opt(double* U, double* V, double* W, double X, double Y, double Z, int IT,  double* VORTIC, double* QW, int JB, double* Time);
void F_Wing_Influence_on_Wake_opt(double* U, double* V, double* W, double X, double Y, double Z,  double* GAMA, double T, int IB, int JB, double* Blade_1_QF);
void F_Vortex(double* U, double* V, double* W, double X, double Y, double Z, double X1, double Y1, double Z1, double X2, double Y2, double Z2, double GAMA, double T);
void F_Vortex_segment_commun_wing_influence(double* U, double* V, double* W, double X, double Y, double Z, double X1, double Y1, double Z1, double X2, double Y2, double Z2, double GAMA1, double GAMA2, double T);
void F_Vortex_segment_commun_wake_influence(double* U, double* V, double* W, double X, double Y, double Z, double X1, double Y1, double Z1, double X2, double Y2, double Z2, double GAMA1, double GAMA2, double T1, double T2);



int main(int argc, char *argv[])
{
	// Pour afficher ce qu'il y a dans le buffer
	char buffer[BUFSIZE];

	//------------------------------------ Ouverture de MATLAB --------------------------------------//
	Engine *ep;
	if (!(ep = engOpen("")))
	{
		fprintf(stderr, "\nCan't start MATLAB engine\n");
		return EXIT_FAILURE;
	}
	//----------------------------------------------------------------------------------------------//



	//-------------------------------------- Initialisation ----------------------------------------//
	engEvalString(ep, "UVLM_2BladeAxial_Init");
	NSTEPS = (int) mxGetPr(engGetVariable(ep, "NSTEPS"))[0];
	IB = (int) mxGetPr(engGetVariable(ep, "IB"))[0];
	JB = (int) mxGetPr(engGetVariable(ep, "JB"))[0];
	NW = (int) mxGetPr(engGetVariable(ep, "NW"))[0];
	CH = (int) mxGetPr(engGetVariable(ep, "CH"))[0];
	//----------------------------------------------------------------------------------------------//


	//-----------Variables n'ayant pas besoin d'être tout le temps recréé pour la mémoire ----------//
	mxArray *ITmat = mxCreateDoubleMatrix(1, 1, mxREAL);
	double *IT = mxGetPr(ITmat);
	/*       Sorties de la fonction F_Wake_Rollup        */
	mxArray *QWnplusunmat  = mxCreateDoubleMatrix(NSTEPS*(JB+1)*3, 1, mxREAL);
	double  *QWnplusun     = mxGetPr(QWnplusunmat)                           ;
	/*---------------------------------------------------*/
	//----------------------------------------------------------------------------------------------//


	//----------------------------------- Boucle sur le temps ------------------------------------- //

	for (int ITc=1; ITc<=NSTEPS; ITc++)
	{
		IT[0] = ITc;
		engPutVariable(ep, "IT", ITmat);

		/*  Partie de la booucle MATLAB avant parallélisation  */
		engEvalString(ep, "UVLM_2BladeAxial_Boucle1");
		DT = (double) mxGetPr(engGetVariable(ep, "DT"))[0];
		/*--------  -------------------------------------------*/

		if (ITc!=1)
		{
			//Pour l'économie de calcul, le paramètre IW peut être utilisé
			//NW est le nombre de segment de sillage déformant (dans le temps)
			int  IW;
			int  I1;

			IW=1;
			if (ITc >= NW)
			{
				IW = ITc-NW+1;
			}
			I1 = ITc-1;
			//int JS1=0; //présent dans le code MATLAB mais inutilisé
			//int JS2=0; //présent dans le code MATLAB mais inutilisé


			/*      Paramètres de la fonction F_Wake_Rollup      */
			mxArray *QWmat         = engGetVariable(ep, "QW")         ;
			mxArray *QWAmat        = engGetVariable(ep, "QWA")        ;
			mxArray *Blade_1_QFmat = engGetVariable(ep, "Blade_1_QF") ;
			mxArray *Blade_2_QFmat = engGetVariable(ep, "Blade_2_QF") ;
			mxArray *Wake_Agemat   = engGetVariable(ep, "Wake_Age")   ;
			mxArray *VORTICmat     = engGetVariable(ep, "VORTIC")     ;
			mxArray *GAMAmat       = engGetVariable(ep, "GAMA")       ;
			mxArray *Tmat          = engGetVariable(ep, "T")          ;
			double  *QW            = mxGetPr(QWmat)                   ;
			double  *QWA           = mxGetPr(QWAmat)                  ;
			double  *Blade_1_QF    = mxGetPr(Blade_1_QFmat)           ;
			double  *Blade_2_QF    = mxGetPr(Blade_2_QFmat)           ;
			double  *Wake_Age      = mxGetPr(Wake_Agemat)             ;
			double  *VORTIC        = mxGetPr(VORTICmat)               ;
			double  *GAMA          = mxGetPr(GAMAmat)                 ;
			double  *T             = mxGetPr(Tmat)                    ;
			/*---------------------------------------------------*/

			#pragma omp parallel default(shared)
			{
				double* Ubc=(double*)malloc(sizeof(double)), * Vbc=(double*)malloc(sizeof(double)), * Wbc=(double*)malloc(sizeof(double));
				int Ic, Jc;
				#pragma omp for schedule(static)
				for (int N=0; N<=(I1-IW+1)*(JB+1)-1;N++)
				{
					Ic = N%(I1-IW+1) + IW;
					Jc = N/(I1-IW+1) + 1;
					F_Wake_Rollup(Ubc,Vbc,Wbc,QW(Ic,Jc,1),QW(Ic,Jc,2),QW(Ic,Jc,3),ITc,QW,QWA,VORTIC,GAMA,CH,IB,JB,T[0],Blade_1_QF,Blade_2_QF,Wake_Age);
					QWnplusun(Ic,Jc,1)=QW(Ic,Jc,1)+DT*Ubc[0];
					QWnplusun(Ic,Jc,2)=QW(Ic,Jc,2)+DT*Vbc[0];
					QWnplusun(Ic,Jc,3)=QW(Ic,Jc,3)+DT*Wbc[0];
				}
				free(Ubc)            ;
				free(Vbc)            ;
				free(Wbc)            ;
			}

			engPutVariable(ep, "QWnplusun", QWnplusunmat);

			/* Partie de la booucle MATLAB après parallélisation */
			engEvalString(ep, "UVLM_2BladeAxial_apresParallel");
			/*---------------------------------------------------*/

			mxDestroyArray(QWmat)        ;
			mxDestroyArray(QWAmat)       ;
			mxDestroyArray(Blade_1_QFmat);
			mxDestroyArray(Blade_2_QFmat);
			mxDestroyArray(Wake_Agemat)  ;
			mxDestroyArray(VORTICmat)    ;
			mxDestroyArray(GAMAmat)      ;
			mxDestroyArray(Tmat)         ;
		}

		/*       Partie de la boucle valable pour ITc=1      */
		engEvalString(ep, "UVLM_2BladeAxial_Boucle2");
		/*---------------------------------------------------*/

		// affichage du buffer (à enlever pour gain de temps)
		engOutputBuffer(ep, buffer, BUFSIZE);
		printf("%s\n", buffer+3);
	}
	//--------------------------------------------------------------------------------------------- //

	mxDestroyArray(ITmat)        ;
	mxDestroyArray(QWnplusunmat) ;

	//--------------------------------------- Sauvegarde ----------------------------------------- //
	engEvalString(ep, "UVLM_2BladeAxial_Save");
	//-------------------------------------------------------------------------------------------- //
	engEvalString(ep, "close all; clear all;");

	//---------------------------- Nettoyage et fermeture de MATLAB ------------------------------- //
	engClose(ep);
	//-------------------------------------------------------------------------------------------- //

	return 0;

}




/** @brief Calcule le vecteur vitesse \f$(u,v,w)\f$ induit en un point \f$(x,y,z)\f$ par les panneaux situés au niveau des deux ailes (wing_influence) et derrière les deux ailes (wake_influence) en prenant en compte l'effet du sol. Les valeurs de U[0], V[0], W[0] sont modifiées par la subroutine. Cette subroutine sert principalement dans le code à calculer la vitesse de déplacement d'un côté d'un panneau derrière l'aile (sillage) pendant un pas de temps.
*
* Exemple : On considère \f$(x,y,z)=(QW(i,j,1),QW(i,j,2),QW(i,j,3))\f$. On applique cette fonction et on obtient le vecteur vitesse \f$(u,v,w)\f$. La nouvelle position du côté du panneau de sillage sera calculé ultérieurement au temps T :  \f$(x,y,z)=(x,y,z)+dt \cdot (u,v,w) \f$.
*
* \paragraph sec Détails mathématiques
* \f$(u,v,w)=\sum_{i=1}^{2}(u_{wake\_influence\_i},v_{wake\_influence\_i},w_{wake\_influence\_i})_{(x,y,z)}+(u_{wing\_influence\_i},v_{wing\_influence\_i},w_{wing\_influence\_i})_{(x,y,z)} + \delta_{CH<=100} \cdot ((u_{wake\_ground\_influence\_i},v_{wake\_ground\_influence\_i},w_{wake\_ground\_influence\_i})_{(x,y,z)}+(u_{wing\_ground\_influence\_i},v_{wing\_ground\_influence\_i},w_{wing\_ground\_influence\_i})_{(x,y,z)})\f$
*
* \n
* Où :
*  - \f$(u_{wake\_influence\_i},v_{wake\_influence\_i},w_{wake\_influence\_i})_{(x,y,z)}\f$ est le vecteur vitesse induit par tout les panneaux se situant derrière l'aile i. \n
*  - \f$(u_{wing\_influence\_i},v_{wing\_influence\_i},w_{wing\_influence\_i})_{(x,y,z)}\f$ est le vecteur vitesse induit par tout les panneaux se situant au niveau de l'aile i.
*  - \f$(u_{wake\_ground\_influence\_i},v_{wake\_ground\_influence\_i},w_{wake\_ground\_influence\_i})_{(x,y,z)}\f$ est le vecteur vitessse prenant en compte l'effet du sol induit par tout les panneaux se situant derrière l'aile i, grâce à la méthode d'effet mirroir : \f$(u_{wake\_ground\_influence\_i},v_{wake\_ground\_influence\_i},w_{wake\_ground\_influence\_i})_{(x,y,z)}=(u_{wake\_influence\_i},v_{wake\_influence\_i},w_{wake\_influence\_i})_{(x,y,-z)}\f$
*  - \f$(u_{wing\_ground\_influence\_i},v_{wing\_ground\_influence\_i},w_{wing\_ground\_influence\_i})_{(x,y,z)}\f$ est le vecteur vitessse prenant en compte l'effet du sol induit par tout les panneaux se situant au niveau de l'aile i, grâce à la méthode d'effet mirroir : \f$(u_{wing\_ground\_influence\_i},v_{wing\_ground\_influence\_i},w_{wing\_ground\_influence\_i})_{(x,y,z)})=(u_{wing\_ground\_influence\_i},v_{wing\_ground\_influence\_i},w_{wing\_ground\_influence\_i})_{(x,y,-z)}\f$
*
* \n
* @param[in,out]  double*  U : vitesse suivant x induite par les panneaux situés au niveau des deux ailes (wing_influence) et derrière les deux ailes (wake_influence) en prenant en compte l'effet du sol.
* @param[in,out] double*  V : vitesse suivant y induite par les panneaux situés au niveau des deux ailes (wing_influence) et derrière les deux ailes (wake_influence) en prenant en compte l'effet du sol
* @param[in,out] double*  W : vitesse suivant z induite par les panneaux situés au niveau des deux ailes (wing_influence) et derrière les deux ailes (wake_influence) en prenant en compte l'effet du sol.
* @param[in]  double  X : 1ère coordonnée du point où la vitesse est calculée.
* @param[in]  double Y : 2ème coordonnée du point où la vitesse est calculée.
* @param[in]  double Z : 3ème coordonnée du point où la vitesse est calculée.
* @param[in]  int IT : Numéro de l'itération actuelle.
* @param[in]  double* QW : panneaux liés aux sillages derrière la pale 1.
* @param[in]  double* QWA : panneaux liés aux sillages derrière la pale 2.
* @param[in]  double* VORTIC : tableau contenant l'ensemble des forces \f$ \Gamma \f$ des segments des panneaux liés au sillage.
* @param[in]  double* GAMA : tableau contenant l'ensemble des forces \f$ \Gamma \f$ des segments des panneaux au niveau de l'aile.
* @param[in]  int CH : Garde au sol (GROUND CLEARANCE)
* @param[in]  int IB : Nombre de panneaux dans le sens de la corde
* @param[in]  int JB : Nombre de panneaux dans le sens de l'envergure (SPANWISE)
* @param[in]  double T : temps actuel
* @param[in]  double* Blade\_1\_QF :  panneaux sur la pale 1
* @param[in]  double* Blade\_2\_QF  :  panneaux sur la pale 2
* @param[in]  double* Wake\_Age : tableau contenant les temps auquels les rangées de panneaux du sillage ont été créées.
*/
void F_Wake_Rollup(double* U, double* V, double* W, double X, double Y, double Z, int IT, double* QW, double* QWA, double* VORTIC, double* GAMA, int CH, int IB, int JB, double T, double* Blade_1_QF, double* Blade_2_QF, double* Wake_Age)
{
	double* U1,* V1,* W1;
	U1 = (double*)malloc(sizeof(double)*1),	V1 = (double*)malloc(sizeof(double)*1), W1 = (double*)malloc(sizeof(double)*1);
	double* U2,* V2,* W2;
	U2 = (double*)malloc(sizeof(double)*1),	V2 = (double*)malloc(sizeof(double)*1), W2 = (double*)malloc(sizeof(double)*1);

	F_Wake_Influence_opt(U1, V1, W1,X,Y,Z,IT,VORTIC,QW,JB,Wake_Age);
	F_Wake_Influence_opt(U2, V2, W2,X,Y,Z,IT,VORTIC,QWA,JB,Wake_Age);
	U[0]=U1[0]+U2[0];
	V[0]=V1[0]+V2[0];
	W[0]=W1[0]+W2[0];

	F_Wing_Influence_on_Wake_opt(U1, V1, W1,X,Y,Z,GAMA,T,IB,JB,Blade_1_QF);
	F_Wing_Influence_on_Wake_opt(U2, V2, W2,X,Y,Z,GAMA,T,IB,JB,Blade_2_QF);

	U[0]+=U1[0]+U2[0];
	V[0]+=V1[0]+V2[0];
	W[0]+=W1[0]+W2[0];

	if (CH <= 100.0)
	{
		F_Wake_Influence_opt(U1, V1, W1,X,Y,-Z,IT,VORTIC,QW,JB,Wake_Age);
		F_Wake_Influence_opt(U2, V2, W2,X,Y,-Z,IT,VORTIC,QWA,JB,Wake_Age);

		U[0]+=U1[0]+U2[0];
		V[0]+=V1[0]+V2[0];
		W[0]+=W1[0]+W2[0];

		F_Wing_Influence_on_Wake_opt(U1, V1, W1,X,Y,-Z,GAMA,T,IB,JB,Blade_1_QF);
		F_Wing_Influence_on_Wake_opt(U2, V2, W2,X,Y,-Z,GAMA,T,IB,JB,Blade_2_QF);

		U[0]+=U1[0]+U2[0];
		V[0]+=V1[0]+V2[0];
		W[0]+=W1[0]+W2[0];
	}
	free(U1);
	free(V1);
	free(W1);
	free(U2);
	free(V2);
	free(W2);
}




/**
* @brief Calcule le vecteur vitesse \f$(u_{wake\_influence\_i},v_{wake\_influence\_i},w_{wake\_influence\_i})\f$ induit en un point \f$(x,y,z)\f$ par les panneaux situés derrière l'aile au temps T. Les valeurs U[0], V[0] et W[0] sont modifiées par la subroutine.
*
* @image html panneauxsillage.jpg "Figure 3. Influence d'un panneau de sillage au point P"  width=350cm height=221cm
*
* \n
*
* \paragraph sec Détails mathématiques
* \f$(u_{wake\_influence\_i},v_{wake\_influence\_i},w_{wake\_influence\_i})=\sum_{i=1}^{IT-1} \sum_{j=1}^{JB} (u_{i,j},v_{i,j},w_{i,j})\f$
*
* \n
* Où :
*  - \f$(u_{i,j},v_{i,j},w_{i,j})\f$ est le vecteur vitesse induit par le panneau (i,j) : \f$(u_{i,j},v_{i,j},w_{i,j}) = (u_{1},v_{1},w_{1}) + (u_{2},v_{2},w_{2}) + (u_{3},v_{3},w_{3}) + (u_{4},v_{4},w_{4})\f$ \n
*  - \f$(u_{1},v_{1},w_{1})\f$ est le vecteur vitesse induit au point \f$(x,y,z)\f$ due à le segment 1 de force \f$\Gamma_{i,j}\f$ (par unité de longueur), pointant vers la direction \f$(x_{i+1,j},y_{i+1,j},z_{i+1,j})-(x_{i,j},y_{i,j},z_{i,j})\f$, au temps Time(i), calculé à partir de F_Vortex. \n
*  - \f$(u_{2},v_{2},w_{2})\f$ est le vecteur vitesse induit au point \f$(x,y,z)\f$ due à le segment 2 de force \f$\Gamma_{i,j}\f$ (par unité de longueur), pointant vers la direction \f$(x_{i+1,j+1},y_{i+1,j+1},z_{i+1,j+1})-(x_{i+1,j},y_{i+1,j},z_{i+1,j})\f$, au temps Time(i), calculé à partir de F_Vortex. \n
*  - \f$(u_{3},v_{3},w_{3})\f$ est le vecteur vitesse induit au point \f$(x,y,z)\f$ due à le segment 3 de force \f$\Gamma_{i,j}\f$ (par unité de longueur), pointant vers la direction \f$(x_{i,j+1}, y_{i,j+1}, z_{i,j+1})-(x_{i+1,j+1},y_{i+1,j+1},z_{i+1,j+1})\f$, au temps Time(i), calculé à partir de F_Vortex. \n
*  - \f$(u_{4},v_{4},w_{4})\f$ est le vecteur vitesse induit au point \f$(x,y,z)\f$ due à le segment 4 de force \f$\Gamma_{i,j}\f$ (par unité de longueur), pointant vers la direction \f$(x_{i,j},y_{i,j},z_{i,j})-(x_{i,j+1},y_{i,j+1},z_{i,j+1})\f$, au temps Time(i), calculé à partir de F_Vortex.
*
* \n
* @param[in,out] double* U : vitesse suivant x induite par les panneaux situés derrière l'aile au temps T
* @param[in,out] double* V : vitesse suivant y induite par les panneaux situés derrière l'aile au temps T
* @param[in,out] double* W : vitesse suivant z induite par les panneaux situés derrière l'aile au temps T
* @param[in] double X : 1ère coordonnée du point où la vitesse est calculée.
* @param[in] double Y : 2ème coordonnée du point où la vitesse est calculée.
* @param[in] double Z : 3ème coordonnée du point où la vitesse est calculée.
* @param[in] int IT : Numéro de l'itération actuelle.
* @param[in] double* VORTIC : tableau contenant l'ensemble des forces \f$ \Gamma \f$ des segments des panneaux liés au sillage.
* @param[in] double* QW : panneaux liés aux sillages derrière une pale.
* @param[in] int JB : Nombre de panneaux dans le sens de l'envergure (SPANWISE)
* @param[in] double* Time : tableau contenant les temps auquels les rangées de panneaux du sillage ont été créées.
*/
void F_Wake_Influence_opt(double* U, double* V, double* W, double X, double Y, double Z, int IT,  double* VORTIC, double* QW, int JB, double* Time)
{
	double* U1,* V1,* W1;
	U1 = (double*)malloc(sizeof(double)*1),	V1 = (double*)malloc(sizeof(double)*1), W1 = (double*)malloc(sizeof(double)*1);
	U1[0]=0.0, V1[0]=0.0, W1[0]=0.0;
	double VORTEK;
	int I1=IT-1;
	int I, J;
	for (int ip=1; ip<=I1; ip++)
	{
		I=ip;
		F_Vortex(U1,V1,W1,X,Y,Z,QW(I,1,1),QW(I,1,2),QW(I,1,3),QW(I+1,1,1),QW(I+1,1,2),QW(I+1,1,3),VORTIC(I,1),Time[I-1]);
		U[0]+=U1[0];
		V[0]+=V1[0];
		W[0]+=W1[0];
		for (int J=2; J<=JB; J++) //on se balade suivant les arrêtes suivant une ligne
		{
			F_Vortex_segment_commun_wake_influence(U1,V1,W1,X,Y,Z,QW(I,J,1),QW(I,J,2),QW(I,J,3),QW(I+1,J,1),QW(I+1,J,2),QW(I+1,J,3),VORTIC(I,J),VORTIC(I,J-1),Time[I-1],Time[I-1]);
			U[0]+=U1[0];
			V[0]+=V1[0];
			W[0]+=W1[0];
		}
		F_Vortex(U1,V1,W1,X,Y,Z,QW(I,JB+1,1),QW(I,JB+1,2),QW(I,JB+1,3),QW(I+1,JB+1,1),QW(I+1,JB+1,2),QW(I+1,JB+1,3),VORTIC(I,JB),Time[I-1]);
		U[0]-=U1[0];
		V[0]-=V1[0];
		W[0]-=W1[0];
	}
	for (int jp=1; jp<=JB; jp++)
	{
		J=jp;
		F_Vortex(U1,V1,W1,X,Y,Z,QW(1,J+1,1),QW(1,J+1,2),QW(1,J+1,3),QW(1,J,1),QW(1,J,2),QW(1,J,3),VORTIC(1,J),Time[0]);
		U[0]+=U1[0];
		V[0]+=V1[0];
		W[0]+=W1[0];
		for (int I=2; I<=I1; I++) //on se balade suivant les arrêtes suivant une colonne
		{
			F_Vortex_segment_commun_wake_influence(U1,V1,W1,X,Y,Z,QW(I,J+1,1),QW(I,J+1,2),QW(I,J+1,3),QW(I,J,1),QW(I,J,2),QW(I,J,3),VORTIC(I,J),VORTIC(I-1,J),Time[I-1],Time[I-2]);
			U[0]+=U1[0];
			V[0]+=V1[0];
			W[0]+=W1[0];
		}
		F_Vortex(U1,V1,W1,X,Y,Z,QW(I1+1,J+1,1),QW(I1+1,J+1,2),QW(I1+1,J+1,3),QW(I1+1,J,1),QW(I1+1,J,2),QW(I1+1,J,3),VORTIC(I1,J),Time[I1-1]);
		U[0]-=U1[0];
		V[0]-=V1[0];
		W[0]-=W1[0];
	}
	free(U1);
	free(V1);
	free(W1);
}




/**
* @brief Calcule le vecteur vitesse \f$(u_{wing\_influence\_i},v_{wing\_influence\_i},w_{wing\_influence\_i})\f$ induit en un point \f$(𝑥,𝑦,𝑧)\f$ par les panneaux situés au niveau de l'aile. Les valeurs de U[0], V[0], W[0] sont modifiées par la subroutine.
*
* @image html panneauxblade.jpg "Figure 4. Influence d'un panneau lié à l'aile au point P"  width=350cm height=221cm
*
* \n
*
* \paragraph sec Détails mathématiques
* \f$(u_{wing\_influ},v_{wing\_influ},w_{wing\_influ})=\sum_{i=1}^{IB} \sum_{j=1}^{JB} (u_{i,j},v_{i,j},w_{i,j})\f$
*
* \n
* Où :
*  - \f$(u_{i,j},v_{i,j},w_{i,j})\f$ est le vecteur vitesse induit par le panneau (i,j) : \f$(u_{i,j},v_{i,j},w_{i,j}) = (u_{1},v_{1},w_{1}) + (u_{2},v_{2},w_{2}) + (u_{3},v_{3},w_{3}) + (u_{4},v_{4},w_{4})\f$ \n
*  - \f$(u_{1},v_{1},w_{1})\f$ est le vecteur vitesse induit au point \f$(x,y,z)\f$ due à le segment 1 de force \f$\Gamma_{i,j}\f$ (par unité de longueur), pointant vers la direction \f$(x_{i+1,j},y_{i+1,j},z_{i+1,j})-(x_{i,j},y_{i,j},z_{i,j})\f$, calculé à partir de F_Vortex. \n
*  - \f$(u_{2},v_{2},w_{2})\f$ est le vecteur vitesse induit au point \f$(x,y,z)\f$ due à le segment 2 de force \f$\Gamma_{i,j}\f$ (par unité de longueur), pointant vers la direction \f$(x_{i+1,j+1},y_{i+1,j+1},z_{i+1,j+1})-(x_{i+1,j},y_{i+1,j},z_{i+1,j})\f$, calculé à partir de F_Vortex. \n
*  - \f$(u_{3},v_{3},w_{3})\f$ est le vecteur vitesse induit au point \f$(x,y,z)\f$ due à le segment 3 de force \f$\Gamma_{i,j}\f$ (par unité de longueur), pointant vers la direction \f$(x_{i,j+1}, y_{i,j+1}, z_{i,j+1})-(x_{i+1,j+1},y_{i+1,j+1},z_{i+1,j+1})\f$, calculé à partir de F_Vortex. \n
*  - \f$(u_{4},v_{4},w_{4})\f$ est le vecteur vitesse induit au point \f$(x,y,z)\f$ due à le segment 4 de force \f$\Gamma_{i,j}\f$ (par unité de longueur), pointant vers la direction \f$(x_{i,j},y_{i,j},z_{i,j})-(x_{i,j+1},y_{i,j+1},z_{i,j+1})\f$, calculé à partir de F_Vortex.
*
* \n
* @param[in,out] double* U : vitesse induite suivant x par les panneaux situés au niveau de l'aile.
* @param[in,out] double* V : vitesse induite suivant y par les panneaux situés au niveau de l'aile.
* @param[in,out] double* W : vitesse induite suivant z par les panneaux situés au niveau de l'aile.
* @param[in] double X : 1ère coordonnée du point où la vitesse est calculée.
* @param[in] double Y : 2ème coordonnée du point où la vitesse est calculée.
* @param[in] double Z : 3ème coordonnée du point où la vitesse est calculée.
* @param[in] double* GAMA : tableau contenant l'ensemble des forces \f$ \Gamma \f$ des segments des panneaux au niveau de l'aile.
* @param[in] double T : temps actuel.
* @param[in]  int IB : Nombre de panneaux dans le sens de la corde
* @param[in]  int JB : Nombre de panneaux dans le sens de l'envergure (SPANWISE)
* @param[in] double* Blade\_1\_QF : tableau contenant l'ensemble des coins des panneaux au niveau de l'une des ailes.
*/
void F_Wing_Influence_on_Wake_opt(double* U, double* V, double* W, double X, double Y, double Z,  double* GAMA, double T, int IB, int JB, double* Blade_1_QF)
{
	double* U1,* V1, * W1;
	U1 = (double*)malloc(sizeof(double)*1),	V1 = (double*)malloc(sizeof(double)*1), W1 = (double*)malloc(sizeof(double)*1);
	int I,J;
	U[0]=0.0;
	V[0]=0.0;
	W[0]=0.0;

	for (int ip=1; ip<=IB; ip++)
	{
		I=ip;
		F_Vortex(U1,V1,W1,X,Y,Z,Blade_1_QF(I+1,1,1),Blade_1_QF(I+1,1,2),Blade_1_QF(I+1,1,3),Blade_1_QF(I,1,1),Blade_1_QF(I,1,2),Blade_1_QF(I,1,3),GAMA(I,1),T);
		U[0]+=U1[0];
		V[0]+=V1[0];
		W[0]+=W1[0];
		for (int J=2; J<=JB; J++) //on se balade suivant les arrêtes suivant une ligne
		{
			F_Vortex_segment_commun_wing_influence(U1,V1,W1,X,Y,Z,Blade_1_QF(I+1,J,1),Blade_1_QF(I+1,J,2),Blade_1_QF(I+1,J,3),Blade_1_QF(I,J,1),Blade_1_QF(I,J,2),Blade_1_QF(I,J,3),GAMA(I,J),GAMA(I,J-1),T);
			U[0]+=U1[0];
			V[0]+=V1[0];
			W[0]+=W1[0];
		}
		F_Vortex(U1,V1,W1,X,Y,Z,Blade_1_QF(I+1,JB+1,1),Blade_1_QF(I+1,JB+1,2),Blade_1_QF(I+1,JB+1,3),Blade_1_QF(I,JB+1,1),Blade_1_QF(I,JB+1,2),Blade_1_QF(I,JB+1,3),GAMA(I,JB),T);
		U[0]-=U1[0];
		V[0]-=V1[0];
		W[0]-=W1[0];
	}
	for (int jp=1; jp<=JB; jp++)
	{
		J=jp;
		F_Vortex(U1,V1,W1,X,Y,Z,Blade_1_QF(1,J,1),Blade_1_QF(1,J,2),Blade_1_QF(1,J,3),Blade_1_QF(1,J+1,1),Blade_1_QF(1,J+1,2),Blade_1_QF(1,J+1,3),GAMA(1,J),T);
		U[0]+=U1[0];
		V[0]+=V1[0];
		W[0]+=W1[0];
		for (int I=2; I<=IB; I++) //on se balade suivant les arrêtes suivant une colonne
		{
			F_Vortex_segment_commun_wing_influence(U1,V1,W1,X,Y,Z,Blade_1_QF(I,J,1),Blade_1_QF(I,J,2),Blade_1_QF(I,J,3),Blade_1_QF(I,J+1,1),Blade_1_QF(I,J+1,2),Blade_1_QF(I,J+1,3),GAMA(I,J),GAMA(I-1,J),T);
			U[0]+=U1[0];
			V[0]+=V1[0];
			W[0]+=W1[0];
		}
		F_Vortex(U1,V1,W1,X,Y,Z,Blade_1_QF(IB+1,J,1),Blade_1_QF(IB+1,J,2),Blade_1_QF(IB+1,J,3),Blade_1_QF(IB+1,J+1,1),Blade_1_QF(IB+1,J+1,2),Blade_1_QF(IB+1,J+1,3),GAMA(IB,J),T);
		U[0]-=U1[0];
		V[0]-=V1[0];
		W[0]-=W1[0];
	}
	free(U1);
	free(V1);
	free(W1);
}




void F_Wake_Influence(double* U, double* V, double* W, double X, double Y, double Z, int IT,  double* VORTIC, double* QW, int JB, double* Time)
{
	double* U1,* V1,* W1;
	U1 = (double*)malloc(sizeof(double)*1),	V1 = (double*)malloc(sizeof(double)*1), W1 = (double*)malloc(sizeof(double)*1);
	U1[0]=0.0, V1[0]=0.0, W1[0]=0.0;
	double VORTEK;
	int I1=IT-1;
	for (int J=1; J<=JB; J++)
	{
		for (int I=1; I<=I1; I++)
		{
			VORTEK=VORTIC(I,J);

			F_Vortex(U1,V1,W1,X,Y,Z,QW(I,J,1),QW(I,J,2),QW(I,J,3),QW(I+1,J,1),QW(I+1,J,2),QW(I+1,J,3),VORTEK,Time[I-1]);
			U[0]+=U1[0];
			V[0]+=V1[0];
			W[0]+=W1[0];

			F_Vortex(U1,V1,W1,X,Y,Z,QW(I+1,J,1),QW(I+1,J,2),QW(I+1,J,3),QW(I+1,J+1,1),QW(I+1,J+1,2),QW(I+1,J+1,3),VORTEK,Time[I-1]);
			U[0]+=U1[0];
			V[0]+=V1[0];
			W[0]+=W1[0];

			F_Vortex(U1,V1,W1,X,Y,Z,QW(I+1,J+1,1),QW(I+1,J+1,2),QW(I+1,J+1,3),QW(I,J+1,1),QW(I,J+1,2),QW(I,J+1,3),VORTEK,Time[I-1]);
			U[0]+=U1[0];
			V[0]+=V1[0];
			W[0]+=W1[0];

			F_Vortex(U1,V1,W1,X,Y,Z,QW(I,J+1,1),QW(I,J+1,2),QW(I,J+1,3),QW(I,J,1),QW(I,J,2),QW(I,J,3),VORTEK,Time[I-1]);
			U[0]+=U1[0];
			V[0]+=V1[0];
			W[0]+=W1[0];

		}
	}
	free(U1);
	free(V1);
	free(W1);
}


void F_Wing_Influence_on_Wake(double* U, double* V, double* W, double X, double Y, double Z,  double* GAMA, double T, int IB, int JB, double* Blade_1_QF)
{
	double* U1,* V1, * W1;
	U1 = (double*)malloc(sizeof(double)*1),	V1 = (double*)malloc(sizeof(double)*1), W1 = (double*)malloc(sizeof(double)*1);
	int I, J;
	U[0]=0.0;
	V[0]=0.0;
	W[0]=0.0;

	for (I=1; I<=IB; I++)
	{
		for (J=1; J<=JB; J++)
		{
			F_Vortex(U1,V1,W1,X,Y,Z,Blade_1_QF(I,J,1),Blade_1_QF(I,J,2),Blade_1_QF(I,J,3),Blade_1_QF(I,J+1,1),Blade_1_QF(I,J+1,2),Blade_1_QF(I,J+1,3),GAMA(I,J),T);
			U[0]+=U1[0];
			V[0]+=V1[0];
			W[0]+=W1[0];
			F_Vortex(U1,V1,W1,X,Y,Z,Blade_1_QF(I,J+1,1),Blade_1_QF(I,J+1,2),Blade_1_QF(I,J+1,3),Blade_1_QF(I+1,J+1,1),Blade_1_QF(I+1,J+1,2),Blade_1_QF(I+1,J+1,3),GAMA(I,J),T);
			U[0]+=U1[0];
			V[0]+=V1[0];
			W[0]+=W1[0];
			F_Vortex(U1,V1,W1,X,Y,Z,Blade_1_QF(I+1,J+1,1),Blade_1_QF(I+1,J+1,2),Blade_1_QF(I+1,J+1,3),Blade_1_QF(I+1,J,1),Blade_1_QF(I+1,J,2),Blade_1_QF(I+1,J,3),GAMA(I,J),T);
			U[0]+=U1[0];
			V[0]+=V1[0];
			W[0]+=W1[0];
			F_Vortex(U1,V1,W1,X,Y,Z,Blade_1_QF(I+1,J,1),Blade_1_QF(I+1,J,2),Blade_1_QF(I+1,J,3),Blade_1_QF(I,J,1),Blade_1_QF(I,J,2),Blade_1_QF(I,J,3),GAMA(I,J),T);
			U[0]+=U1[0];
			V[0]+=V1[0];
			W[0]+=W1[0];
		}
	}
	free(U1);
	free(V1);
	free(W1);
}

/**
* @brief Calcule le vecteur vitesse \f$(u_k,v_k,w_k)\f$ induit en un point \f$(x_p,y_p,z_p)\f$ due à un segment \f$k\f$ de force \f$\Gamma\f$ (par unité de longueur), pointant vers la direction \f$(x_2,y_2,z_2)-(x_1,y_1,z_1)\f$. Les valeurs de U[0], V[0], W[0] sont modifiées par la subroutine.
*
* @image html influensegment.jpg "Figure 5. Influence d'un segment sur le vecteur vitesse au point P"  width=350cm height=221cm middle
*
* \n
*
* \paragraph sec Détails mathématiques
*
*  \f$\mathbf{r}_{0}=\left(x_{2}-x_{1}, y_{2}-y_{1}, z_{2}-z_{1}\right)\f$ ;
* \f$\quad \mathbf{r}_{1}=\left(x_{p}-x_{1}, y_{p}-y_{1}, z p-z_{1}\right)\f$ ;
* \f$\quad \mathbf{r}_{2}=\left(x_{p}-x_{2}, y_{p}-y_{2}, z_{p}-z_{2}\right)\f$ ;
* \f$\quad  a l f=1.25643 ;  \quad\mu=1,343 \times 10^{-5} ;\f$\n
* \f$ a_{1}=10^{-4} ; \quad R C U T=10^{-10} ;\f$
* \f$\quad RC=\sqrt{4 \cdot alf \cdot (1+a_{1} \cdot \frac{|\Gamma|}{\mu}) \cdot \mu \cdot T}\f$ ;
* \f$\quad K=\frac{\Gamma}{4 \pi\left|\mathbf{r}_{1} \times \mathbf{r}_{2}\right|^{2}}\left(\frac{\mathbf{r}_{0} \cdot \mathbf{r}_{1}}{|\mathbf{r}_{1}|}-\frac{\mathbf{r}_{0} \cdot \mathbf{r}_{2}}{|\mathbf{r}_{2}|}\right)\f$;
* \f$\quad LambOseen = \delta_{T=0}+\delta_{T\ne0} \cdot 1-\exp \left(-a l f \cdot (\frac{\left|\mathbf{r}_{1} \times \mathbf{r}_{2}\right|^{2}}{\left(|\mathbf{r}_{0}| \cdot RC\right)^{2}}\right))\f$
*
* \n
* Si \f$\quad (\mathbf{r}_0 < RCUT)\quad\f$ ou \f$\quad (\mathbf{r}_1<RCUT)\quad\f$ ou \f$\quad (\mathbf{r}_2<RCUT)\quad\f$ ou \f$\quad (\left|\mathbf{r}_{1} \times \mathbf{r}_{2}\right|^{2} <RCUT)\quad\f$ :
*
*  - \f$(u_{k},v_{k},w_{k}) = (0,0,0)\f$
*
* Sinon :
*
*  - \f$(u_{k},v_{k},w_{k}) = K \cdot LambOseen  \cdot (\left(\mathbf{r}_{1} \times \mathbf{r}_{2}\right)_{x}, \left(\mathbf{r}_{1} \times \mathbf{r}_{2}\right)_{y}, \left(\mathbf{r}_{1} \times \mathbf{r}_{2}\right)_{z})\f$
*
* \n
* @param[in,out] double* U : vitesse suivant x induite par le segment
* @param[in,out] double* V : vitesse suivant y induite par le segment
* @param[in,out] double* W : vitesse suivant z induite par le segment
* @param[in] double X : 1ère coordonnée du point où la vitesse est calculée.
* @param[in] double Y : 2ème coordonnée du point où la vitesse est calculée.
* @param[in] double Z : 3ème coordonnée du point où la vitesse est calculée.
* @param[in] double X1 : 1ère coordonnée du coin 1 du segment.
* @param[in] double Y1 : 2ème coordonnée du coin 1 du segment.
* @param[in] double Z1 : 3ème coordonnée du coin 1 du segment.
* @param[in] double X2 : 1ère coordonnée du coin 2 du segment.
* @param[in] double Y2 : 2ème coordonnée du coin 2 du segment.
* @param[in] double Z2 : 3ème coordonnée du coin 2 du segment.
* @param[in] double GAMA : force du segment.
* @param[in] double T : temps actuel.
*/
void F_Vortex(double* U, double* V, double* W, double X, double Y, double Z, double X1, double Y1, double Z1, double X2, double Y2, double Z2, double GAMA, double T)
{
	double alf=1.25643, u=1.343*pow(10,-5), a1=pow(10,-4), RCUT=pow(10,-10), pi=4*atan(1);
	double RC,R0R1, R0R2, COEF, LambOseen, R1R2X=(Y-Y1)*(Z-Z2)-(Z-Z1)*(Y-Y2), R1R2Y=-((X-X1)*(Z-Z2)-(Z-Z1)*(X-X2)), R1R2Z=(X-X1)*(Y-Y2)-(Y-Y1)*(X-X2), SQUARE=R1R2X*R1R2X+R1R2Y*R1R2Y+R1R2Z*R1R2Z,R0=sqrt((X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2)+(Z1-Z2)*(Z1-Z2)), R1=sqrt((X-X1)*(X-X1)+(Y-Y1)*(Y-Y1)+(Z-Z1)*(Z-Z1)), R2=sqrt((X-X2)*(X-X2)+(Y-Y2)*(Y-Y2)+(Z-Z2)*(Z-Z2));

	if ((R0<RCUT) || (R1<RCUT) || (R2<RCUT) || (SQUARE<RCUT))
	{
		U[0]=0.0;
		V[0]=0.0;
		W[0]=0.0;
	}
	else
	{
		R0R1=(X2-X1)*(X-X1)+(Y2-Y1)*(Y-Y1)+(Z2-Z1)*(Z-Z1),R0R2=(X2-X1)*(X-X2)+(Y2-Y1)*(Y-Y2)+(Z2-Z1)*(Z-Z2);
		COEF=GAMA/(4.0*pi*SQUARE)*(R0R1/R1-R0R2/R2),RC=sqrt(4*alf*(1+a1*abs(GAMA)/u)*u*T);
		if (T==0)
		{
			LambOseen=1;
		}
		else
		{
			LambOseen=(1-exp(-alf*(SQUARE/pow((R0*RC),2))));
		}
		COEF=COEF*LambOseen, U[0]=R1R2X*COEF, V[0]=R1R2Y*COEF, W[0]=R1R2Z*COEF;
	}
}



/**
* @brief Calcule le vecteur vitesse \f$(u_k,v_k,w_k)\f$ induit en un point \f$(x_p,y_p,z_p)\f$ due à un segment \f$k\f$ de force \f$\Gamma\f$ (par unité de longueur), pointant vers la direction \f$(x_2,y_2,z_2)-(x_1,y_1,z_1)\f$. Cette subroutine est semblable à F_Vortex mais elle prend en compte le fait qu'un segment peut apparaitre au niveau de deux panneaux en aval de l'aile (sillage) et avoir une double influence. Certains calculs sont alors redondants et peuvent être effectués une seule fois. Les valeurs de U[0], V[0], W[0] sont modifiées par la subroutine.
*
* @image html influensegmentblade.jpg "Figure 6. Influence d'un segment de l'aile sur le vecteur vitesse au point P"  width=350cm height=221cm middle
* @image html influen2.jpg "Figure 7. Placement des panneaux et mise en avant des segments communs au niveau de l'aile" width=850cm height=421cm
*
* \n
*
* \paragraph sec Détails mathématiques
*
*  \f$\mathbf{r}_{0}=\left(x_{2}-x_{1}, y_{2}-y_{1}, z_{2}-z_{1}\right)\f$ ;
* \f$\quad \mathbf{r}_{1}=\left(x_{p}-x_{1}, y_{p}-y_{1}, z p-z_{1}\right)\f$ ;
* \f$\quad \mathbf{r}_{2}=\left(x_{p}-x_{2}, y_{p}-y_{2}, z_{p}-z_{2}\right)\f$ ;
* \f$\quad  a l f=1.25643 ;  \quad\mu=1,343 \times 10^{-5} ;\f$\n
* \f$ a_{1}=10^{-4} ; \quad R C U T=10^{-10} ;\f$
* \f$\quad RC1=\sqrt{4 \cdot alf \cdot (1+a_{1} \cdot \frac{|\Gamma_{1}|}{\mu}) \cdot \mu \cdot T1}\f$ ;
* \f$\quad RC2=\sqrt{4 \cdot alf \cdot (1+a_{1} \cdot \frac{|\Gamma_{2}|}{\mu}) \cdot \mu \cdot T2}\f$ ;
* \f$\quad K=\frac{1}{4 \pi\left|\mathbf{r}_{1} \times \mathbf{r}_{2}\right|^{2}}\left(\frac{\mathbf{r}_{0} \cdot \mathbf{r}_{1}}{|\mathbf{r}_{1}|}-\frac{\mathbf{r}_{0} \cdot \mathbf{r}_{2}}{|\mathbf{r}_{2}|}\right)\f$; \n
* \f$LambOseen = \delta_{T=0}+\delta_{T\ne0} \cdot 1-\exp \left(-a l f \cdot (\frac{\left|\mathbf{r}_{1} \times \mathbf{r}_{2}\right|^{2}}{\left(|\mathbf{r}_{0}| \cdot RC1\right)^{2}}\right))\f$;
* \f$\quad LambOseen = \delta_{T=0}+\delta_{T\ne0} \cdot 1-\exp \left(-a l f \cdot (\frac{\left|\mathbf{r}_{1} \times \mathbf{r}_{2}\right|^{2}}{\left(|\mathbf{r}_{0}| \cdot RC2\right)^{2}}\right))\f$
*
* \n
* Si \f$\quad (\mathbf{r}_0 < RCUT)\quad\f$ ou \f$\quad (\mathbf{r}_1<RCUT)\quad\f$ ou \f$\quad (\mathbf{r}_2<RCUT)\quad\f$ ou \f$\quad (\left|\mathbf{r}_{1} \times \mathbf{r}_{2}\right|^{2} <RCUT)\quad\f$ :
*
*  - \f$(u_{k},v_{k},w_{k}) = (0,0,0)\f$
*
* Sinon :
*
*  - \f$(u_{k},v_{k},w_{k}) = (\Gamma_{1} \cdot LambOseen1 - \Gamma_{2} \cdot LambOseen2) \cdot K \cdot LambOseen  \cdot (\left(\mathbf{r}_{1} \times \mathbf{r}_{2}\right)_{x}, \left(\mathbf{r}_{1} \times \mathbf{r}_{2}\right)_{y}, \left(\mathbf{r}_{1} \times \mathbf{r}_{2}\right)_{z})\f$
*
* \n
* @param[in,out] double* U : vitesse suivant x induite par le segment
* @param[in,out] double* V : vitesse suivant y induite par le segment
* @param[in,out] double* W : vitesse suivant z induite par le segment
* @param[in] double X : 1ère coordonnée du point où la vitesse est calculée.
* @param[in] double Y : 2ème coordonnée du point où la vitesse est calculée.
* @param[in] double Z : 3ème coordonnée du point où la vitesse est calculée.
* @param[in] double X1 : 1ère coordonnée du coin 1 du segment.
* @param[in] double Y1 : 2ème coordonnée du coin 1 du segment.
* @param[in] double Z1 : 3ème coordonnée du coin 1 du segment.
* @param[in] double X2 : 1ère coordonnée du coin 2 du segment.
* @param[in] double Y2 : 2ème coordonnée du coin 2 du segment.
* @param[in] double Z2 : 3ème coordonnée du coin 2 du segment.
* @param[in] double GAMA1 : force du segment du premier panneau adjacent.
* @param[in] double GAMA2 : force du segment du deuxième panneau adjacent.
* @param[in] double T : temps actuel
*/
void F_Vortex_segment_commun_wing_influence(double* U, double* V, double* W, double X, double Y, double Z, double X1, double Y1, double Z1, double X2, double Y2, double Z2, double GAMA1, double GAMA2, double T)
{
	double alf=1.25643, u=1.343*pow(10,-5), a1=pow(10,-4), RCUT=pow(10,-10), pi=4*atan(1);
	double RC1, RC2, R0R1, R0R2, COEF, COEF1, COEF2, LambOseen1, LambOseen2, R1R2X=(Y-Y1)*(Z-Z2)-(Z-Z1)*(Y-Y2), R1R2Y=-((X-X1)*(Z-Z2)-(Z-Z1)*(X-X2)), R1R2Z=(X-X1)*(Y-Y2)-(Y-Y1)*(X-X2), SQUARE=R1R2X*R1R2X+R1R2Y*R1R2Y+R1R2Z*R1R2Z,R0=sqrt((X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2)+(Z1-Z2)*(Z1-Z2)), R1=sqrt((X-X1)*(X-X1)+(Y-Y1)*(Y-Y1)+(Z-Z1)*(Z-Z1)), R2=sqrt((X-X2)*(X-X2)+(Y-Y2)*(Y-Y2)+(Z-Z2)*(Z-Z2));

	if ((R0<RCUT) || (R1<RCUT) || (R2<RCUT) || (SQUARE<RCUT))
	{
		U[0]=0.0;
		V[0]=0.0;
		W[0]=0.0;
	}
	else
	{
		R0R1=(X2-X1)*(X-X1)+(Y2-Y1)*(Y-Y1)+(Z2-Z1)*(Z-Z1),R0R2=(X2-X1)*(X-X2)+(Y2-Y1)*(Y-Y2)+(Z2-Z1)*(Z-Z2);
		COEF=1.0/(4.0*pi*SQUARE)*(R0R1/R1-R0R2/R2),RC1=sqrt(4*alf*(1+a1*abs(GAMA1)/u)*u*T), RC2=sqrt(4*alf*(1+a1*abs(GAMA2)/u)*u*T);
		if (T==0)
		{
			LambOseen1=1;
			LambOseen2=1;
		}
		else
		{
			LambOseen1=(1-exp(-alf*(SQUARE/pow((R0*RC1),2))));
			LambOseen2=(1-exp(-alf*(SQUARE/pow((R0*RC2),2))));
		}
		COEF1=GAMA1*COEF*LambOseen1, COEF2=GAMA2*COEF*LambOseen2;
		U[0]=R1R2X*(COEF1-COEF2), V[0]=R1R2Y*(COEF1-COEF2), W[0]=R1R2Z*(COEF1-COEF2);
	}
}



/**
* @brief Calcule le vecteur vitesse \f$(u_k,v_k,w_k)\f$ induit en un point \f$(x_p,y_p,z_p)\f$ due à un segment \f$k\f$ de force \f$\Gamma\f$ (par unité de longueur), pointant vers la direction \f$(x_2,y_2,z_2)-(x_1,y_1,z_1)\f$. Cette subroutine est semblable à F_Vortex mais elle prend en compte le fait qu'un segment peut apparaitre au niveau de deux panneaux sur l'aile et avoir une double influence. Certains calculs sont alors redondants et peuvent être effectués une seule fois. Les valeurs de U[0], V[0], W[0] sont modifiées par la subroutine.
*
* @image html influensegmentsillage.jpg "Figure 8. Influence d'un segment du sillage sur le vecteur vitesse au point P"  width=350cm height=221cm middle
* @image html influen1.jpg "Figure 9. Placement des panneaux et mise en avant des segments communs au niveau du sillage" width=850cm height=421cm
*
* \n
*
* \paragraph sec Détails mathématiques
*
*  \f$\mathbf{r}_{0}=\left(x_{2}-x_{1}, y_{2}-y_{1}, z_{2}-z_{1}\right)\f$ ;
* \f$\quad \mathbf{r}_{1}=\left(x_{p}-x_{1}, y_{p}-y_{1}, z p-z_{1}\right)\f$ ;
* \f$\quad \mathbf{r}_{2}=\left(x_{p}-x_{2}, y_{p}-y_{2}, z_{p}-z_{2}\right)\f$ ;
* \f$\quad  a l f=1.25643 ;  \quad\mu=1,343 \times 10^{-5} ;\f$\n
* \f$ a_{1}=10^{-4} ; \quad R C U T=10^{-10} ;\f$
* \f$\quad RC1=\sqrt{4 \cdot alf \cdot (1+a_{1} \cdot \frac{|\Gamma_{1}|}{\mu}) \cdot \mu \cdot T1}\f$ ;
* \f$\quad RC2=\sqrt{4 \cdot alf \cdot (1+a_{1} \cdot \frac{|\Gamma_{2}|}{\mu}) \cdot \mu \cdot T2}\f$ ;
* \f$\quad K=\frac{1}{4 \pi\left|\mathbf{r}_{1} \times \mathbf{r}_{2}\right|^{2}}\left(\frac{\mathbf{r}_{0} \cdot \mathbf{r}_{1}}{|\mathbf{r}_{1}|}-\frac{\mathbf{r}_{0} \cdot \mathbf{r}_{2}}{|\mathbf{r}_{2}|}\right)\f$;\n
* \f$\ LambOseen1 = \delta_{T1=0}+\delta_{T1\ne0} \cdot 1-\exp \left(-a l f \cdot (\frac{\left|\mathbf{r}_{1} \times \mathbf{r}_{2}\right|^{2}}{\left(|\mathbf{r}_{0}| \cdot RC1\right)^{2}}\right))\f$;
* \f$\quad LambOseen2 = \delta_{T2=0}+\delta_{T2\ne0} \cdot 1-\exp \left(-a l f \cdot (\frac{\left|\mathbf{r}_{1} \times \mathbf{r}_{2}\right|^{2}}{\left(|\mathbf{r}_{0}| \cdot RC2\right)^{2}}\right))\f$
*
* \n
* Si \f$\quad (\mathbf{r}_0 < RCUT)\quad\f$ ou \f$\quad (\mathbf{r}_1<RCUT)\quad\f$ ou \f$\quad (\mathbf{r}_2<RCUT)\quad\f$ ou \f$\quad (\left|\mathbf{r}_{1} \times \mathbf{r}_{2}\right|^{2} <RCUT)\quad\f$ :
*
*  - \f$(u_{k},v_{k},w_{k}) = (0,0,0)\f$
*
* Sinon :
*
*  - \f$(u_{k},v_{k},w_{k}) = (\Gamma_{1} \cdot LambOseen1 - \Gamma_{2} \cdot LambOseen2) \cdot K \cdot (\left(\mathbf{r}_{1} \times \mathbf{r}_{2}\right)_{x}, \left(\mathbf{r}_{1} \times \mathbf{r}_{2}\right)_{y}, \left(\mathbf{r}_{1} \times \mathbf{r}_{2}\right)_{z})\f$
*
* \n
* @param[in,out] double* U : vitesse suivant x induite par le segment
* @param[in,out] double* V : vitesse suivant y induite par le segment
* @param[in,out] double* W : vitesse suivant z induite par le segment
* @param[in] double X : 1ère coordonnée du point où la vitesse est calculée.
* @param[in] double Y : 2ème coordonnée du point où la vitesse est calculée.
* @param[in] double Z : 3ème coordonnée du point où la vitesse est calculée.
* @param[in] double X1 : 1ère coordonnée du coin 1 du segment.
* @param[in] double Y1 : 2ème coordonnée du coin 1 du segment.
* @param[in] double Z1 : 3ème coordonnée du coin 1 du segment.
* @param[in] double X2 : 1ère coordonnée du coin 2 du segment.
* @param[in] double Y2 : 2ème coordonnée du coin 2 du segment.
* @param[in] double Z2 : 3ème coordonnée du coin 2 du segment.
* @param[in] double GAMA1 : force du segment du premier panneau adjacent.
* @param[in] double GAMA2 : force du segment du deuxième panneau adjacent.
* @param[in] double T1 : temps où le panneau 1 à été créé.
* @param[in] double T2 : temps où le panneau 2 à été créé. (le panneau 1 et 2 partagent un même segment)
*/
void F_Vortex_segment_commun_wake_influence(double* U, double* V, double* W, double X, double Y, double Z, double X1, double Y1, double Z1, double X2, double Y2, double Z2, double GAMA1, double GAMA2, double T1, double T2)
{
	double alf=1.25643, u=1.343*pow(10,-5), a1=pow(10,-4), RCUT=pow(10,-10), pi=4*atan(1);
	double RC1, RC2, R0R1, R0R2, COEF, COEF1, COEF2, LambOseen1, LambOseen2, R1R2X=(Y-Y1)*(Z-Z2)-(Z-Z1)*(Y-Y2), R1R2Y=-((X-X1)*(Z-Z2)-(Z-Z1)*(X-X2)), R1R2Z=(X-X1)*(Y-Y2)-(Y-Y1)*(X-X2), SQUARE=R1R2X*R1R2X+R1R2Y*R1R2Y+R1R2Z*R1R2Z,R0=sqrt((X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2)+(Z1-Z2)*(Z1-Z2)), R1=sqrt((X-X1)*(X-X1)+(Y-Y1)*(Y-Y1)+(Z-Z1)*(Z-Z1)), R2=sqrt((X-X2)*(X-X2)+(Y-Y2)*(Y-Y2)+(Z-Z2)*(Z-Z2));

	if ((R0<RCUT) || (R1<RCUT) || (R2<RCUT) || (SQUARE<RCUT))
	{
		U[0]=0.0;
		V[0]=0.0;
		W[0]=0.0;
	}
	else
	{
		R0R1=(X2-X1)*(X-X1)+(Y2-Y1)*(Y-Y1)+(Z2-Z1)*(Z-Z1),R0R2=(X2-X1)*(X-X2)+(Y2-Y1)*(Y-Y2)+(Z2-Z1)*(Z-Z2);
		COEF=1.0/(4.0*pi*SQUARE)*(R0R1/R1-R0R2/R2),RC1=sqrt(4*alf*(1+a1*abs(GAMA1)/u)*u*T1), RC2=sqrt(4*alf*(1+a1*abs(GAMA2)/u)*u*T2);
		if (T1==0)
		{
			LambOseen1=1;
		}
		else
		{
			LambOseen1=(1-exp(-alf*(SQUARE/pow((R0*RC1),2))));
		}
		if (T2==0)
		{
			LambOseen2=1;
		}
		else
		{
			LambOseen2=(1-exp(-alf*(SQUARE/pow((R0*RC2),2))));
		}
		COEF1=GAMA1*COEF*LambOseen1, COEF2=GAMA2*COEF*LambOseen2;
		U[0]=R1R2X*(COEF1-COEF2), V[0]=R1R2Y*(COEF1-COEF2), W[0]=R1R2Z*(COEF1-COEF2);
	}
}
