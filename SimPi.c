#include<stdio.h>
#include<math.h>


/* Period parameters */  
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* initializes mt[N] with a seed */
void init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] = 
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(unsigned long init_key[], int key_length)
{
    int i, j, k;
    init_genrand(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if init_genrand() has not been called, */
            init_genrand(5489UL); /* a default initial seed is used */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }
  
    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(void)
{
    return (long)(genrand_int32()>>1);
}

/* generates a random number on [0,1]-real-interval */
double genrand_real1(void)
{
    return genrand_int32()*(1.0/4294967295.0); 
    /* divided by 2^32-1 */ 
}

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void)
{
    return genrand_int32()*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3(void)
{
    return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(void) 
{ 
    unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6; 
    return(a*67108864.0+b)*(1.0/9007199254740992.0); 
} 

/*-------------------------------------------------------------------------------------------------------*/
/* simPy   calclule la valeur de pi ?? partir d'un certain nombre de points                               */                                                                     
/*																									     */
/* En entr??e : nbr_point est un long entier entr?? par l'utilisateur                                      */
/*																								   	     */			
/* En sortie : la fonction retourne la valeur de pi												         */
/*																									     */
/*-------------------------------------------------------------------------------------------------------*/

double simPi (long int nbr_point)
{
    long int j;
	double x,y;
	int onSurface = 0;                 //variable qui d??signe le nombre de points dans le cercle 

	for(j = 0; j<nbr_point;j++)
	{
		x = genrand_real1();          //abscisse dans un intervalle [0,1] g??n??r?? par MT
		y = genrand_real1();          //ordonn??e dans un intervalle [0,1] g??n??r?? par MT

		if( x*x + y*y <= 1)               //distance entre le point et le centre du cerle
		{
			onSurface++;                //on incr??mente si le point est dans le cercle(distance<=1)
		}
	}
	
	return 4 * (double) onSurface / nbr_point;	       //retourne pi
}

/*----------------------------------------------------------------------------------------------------------*/
/* estimPy  estimation de pi ?? partir des r??sultats de plusieurs exp??riences ind??pendantes                  */
/*                                                                                                          */
/* En entr??e : tab_exp un tableau pour stocker les r??sultats des exp??riences                                */
/*                     nbr_exp  nombre d'exp??riences                                                        */
/*                     nbr_point  nombre de points                                                          */
/*                             																	    		*/
/* En sortie  :  la valeur estim?? de pi				                      								    */
/*----------------------------------------------------------------------------------------------------------*/

double estimPy(double tab_exp[], int nbr_exp, long int nbr_point)
{
	int i;
	double somme = 0;
	
	for(i=0; i < nbr_exp; i++)
	{
		tab_exp[i] = simPi(nbr_point);          //on stocke les valeur de pi dans le tableau
		
		somme = somme + tab_exp[i];        //on calcule la somme des valeurs de pi
	}
	
	return (double) somme / nbr_exp;       //on retourne la moyenne des valeurs de pi
}

/*-----------------------------------------------------------------------------------------------------------------------------*/
/* conf_Int      calcule la marge d'erreur utilis??e dans le main pour le cacul de l'intervalle                                 */
/*                                                                                                                             */
/* En entr??e :  tab_exp le tableau ou on a stock?? les r??sultats des exp??riences                                                */
/*                     nbr_exp  nombre d'exp??riences                                                                           */
/*                     moyenne  estimation de pi                                                                               */
/*													                               											   */
/* En sortie :  la marge d'erreur																							   */
/*-----------------------------------------------------------------------------------------------------------------------------*/

double conf_Int(double tab_exp[], int nbr_exp, double moyenne)
{
	//tableau des valeurs de student law pour des experiences de 1 ?? 30(inclus)

	float t[30] = {12.706, 4.303, 3.182, 2.776, 2.571, 2.447, 2.365, 2.308, 2.262, 2.228,
								2.201, 2.179, 2.160, 2.145, 2.131, 2.120, 2.110, 2.101, 2.093, 2.086,
								2.080, 2.074, 2.069, 2.064, 2.060, 2.056, 2.052, 2.048, 2.045, 2.042};
	int i;
	double s = 0;  //variable repr??sentant l'estimation de la variance

	for(i = 0; i < nbr_exp; i++)
	{
	
		s = s + ((tab_exp[i] - moyenne) * (tab_exp[i] - moyenne));
	}
	
	s = s / (nbr_exp - 1);    //calcule de la variance

	return  t[nbr_exp - 1] * sqrt(s / nbr_exp);       //retourne la marge d'erreur
}

int main ()
{
    unsigned long init[4]={0x123, 0x234, 0x345, 0x456}, length=4;
    init_by_array(init, length);

    long int nbr_point;

	printf("entrer le nombre de points :");
    scanf("%ld", &nbr_point);

	int nbr_exp;

	printf("entrer le nombre d'exp??riences entre 1 et 30(inclus) :");
    scanf("%d", &nbr_exp);
	
	double tab_exp[nbr_exp];
    double moyenne = 0,
    				R = 0;
    
    moyenne = estimPy(tab_exp, nbr_exp,nbr_point);
    printf("estimation de pi : %10.8f\n",moyenne);

    R = conf_Int(tab_exp, nbr_exp, moyenne);
	printf("l'intervalle de confiance : [%10.8f, %10.8f]\n", moyenne - R, moyenne + R);

	init_by_array(init, length);
	printf("calcul de pi pour %ld points : %10.8f\n", nbr_point, simPi(nbr_point));
	
	return 0;
}

































