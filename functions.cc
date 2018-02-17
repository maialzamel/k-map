
/**
    K-MAP: In the k-mappability problem, we are given a string x of length n and integers m and k, and we are asked to count, for each length-m factor y of x, the number of other factors of length m of x that are at Hamming distance at most k from y.
    Copyright (C) 2017 Mai Alzamel and Solon P. Pissis. 

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
**/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
#include <sys/time.h>
#include <vector>
#include "mapdefs.h"


#include <divsufsort64.h>                                         // include header for suffix sort

#include <sdsl/rmq_support.hpp>
#include <sdsl/bit_vectors.hpp>                                      // include header for bit vectors

using namespace sdsl;
using namespace std;

double gettime( void )
{
    struct timeval ttime;
    gettimeofday( &ttime , 0 );
    return ttime.tv_sec + ttime.tv_usec * 0.000001;
};


unsigned int LCParray ( unsigned char *text, INT n, INT * SA, INT * ISA, INT * LCP )
{										
	INT i=0, j=0;

	LCP[0] = 0;
	for ( i = 0; i < n; i++ ) // compute LCP[ISA[i]]
		if ( ISA[i] != 0 ) 
		{
			if ( i == 0) j = 0;
			else j = (LCP[ISA[i-1]] >= 2) ? LCP[ISA[i-1]]-1 : 0;
			while ( text[i+j] == text[SA[ISA[i]-1]+j] )
				j++;
			LCP[ISA[i]] = j;
		}

	return ( 1 );
}

unsigned int compute_map ( unsigned char * seq,unsigned char * seq_id, struct TSwitch sw, INT m ,   double  * C )
{
    
    	unsigned char * seq_prime;
	INT * SA;
	INT * LCP;
	INT * invSA;
    	INT * SA_prime;
    	INT * LCP_prime;
    	INT * invSA_prime;

	INT n = strlen ( ( char * ) seq );
    	INT N = m + n;
    	INT L= floor (m/3.0);
    	INT B=0,w=0,B_i=0,B_j=0;

    	INT l,r,l_prime,r_prime ,r_1,r_2,l_1,l_2,u_1,u_2,L_1,L_2;
    	seq_prime = ( unsigned char * ) calloc ( ( n + 1 ) , sizeof ( unsigned char ) );
    
    	for( INT i=n;i>0;i--)
        	seq_prime[w++]=seq[i-1];
	
        /* Compute the suffix array and LCP array for the seq and inverse seq */
        SA = ( INT * ) malloc( ( n ) * sizeof( INT ) );
        SA_prime = ( INT * ) malloc( ( n ) * sizeof( INT ) );
        if( ( SA == NULL) )
        {
                fprintf(stderr, " Error: Cannot allocate memory for SA.\n" );
                return ( 0 );
        }
        if( ( SA_prime == NULL) )
        {
        fprintf(stderr, " Error: Cannot allocate memory for SA.\n" );
        return ( 0 );
        }

	#ifdef _USE_64
        if( divsufsort64( seq, SA,  n ) != 0 )
        {
                fprintf(stderr, " Error: SA computation failed.\n" );
                exit( EXIT_FAILURE );
        }
	#endif

	#ifdef _USE_32
        if( divsufsort( seq, SA,  n ) != 0 )
        {
                fprintf(stderr, " Error: SA computation failed.\n" );
                exit( EXIT_FAILURE );
        }
	#endif

	#ifdef _USE_64
    	if( divsufsort64( seq_prime, SA_prime,  n ) != 0 )
    	{
        	fprintf(stderr, " Error: SA computation failed.\n" );
        	exit( EXIT_FAILURE );
    	}
	#endif
    
	#ifdef _USE_32
    	if( divsufsort( seq_prime, SA_prime,  n ) != 0 )
    	{
       		fprintf(stderr, " Error: SA computation failed.\n" );
       	 	exit( EXIT_FAILURE );
    	}
	#endif
    
        invSA = ( INT * ) calloc( n , sizeof( INT ) );
        if( ( invSA == NULL) )
        {
                fprintf(stderr, " Error: Cannot allocate memory for invSA.\n" );
                return ( 0 );
        }

        for ( INT i = 0; i < n; i ++ )
        {
                invSA [SA[i]] = i;
        }
    
 
    	invSA_prime = ( INT * ) calloc( n , sizeof( INT ) );
    	if( ( invSA_prime == NULL) )
    	{
        	fprintf(stderr, " Error: Cannot allocate memory for invSA.\n" );
        	return ( 0 );
    	}
    
    	for ( INT i = 0; i < n; i ++ )
    	{
        	invSA_prime [SA_prime[i]] = i;
    	}
   

	LCP = ( INT * ) calloc  ( n, sizeof( INT ) );
        if( ( LCP == NULL) )
        {
                fprintf(stderr, " Error: Cannot allocate memory for LCP.\n" );
                return ( 0 );
        }

        /* Compute the LCP array */
        if( LCParray( seq, n, SA, invSA, LCP ) != 1 )
        {
                fprintf(stderr, " Error: LCP computation failed.\n" );
                exit( EXIT_FAILURE );
        }
    
    	//LCP REVERSED String
    	LCP_prime = ( INT * ) calloc  ( n, sizeof( INT ) );
    	if( ( LCP_prime == NULL) )
    	{
        	fprintf(stderr, " Error: Cannot allocate memory for LCP.\n" );
        	return ( 0 );
    	}

    	if( LCParray( seq_prime, n, SA_prime, invSA_prime, LCP_prime ) != 1 )
    	{
        	fprintf(stderr, " Error: LCP computation failed.\n" );
        	exit( EXIT_FAILURE );
    	}
    
    	int_vector<> v( n , 0 ); // create a vector of length n and initialize it with 0s
    
    	for ( INT i = 0; i < n; i ++ )
    	{
        	v[i] = LCP[i];
    	}
    
    	rmq_succinct_sct<> rmq(&v);
    
    	int_vector<> v_prime( n , 0 ); // create a vector of length n and initialize it with 0s
    
    	for ( INT i = 0; i < n; i ++ )
    	{
        	v_prime[i] = LCP_prime[i];
    	}
    
    	rmq_succinct_sct<> rmq_prime(&v_prime);
    
    
    
    
    
    
	/*1-MAPP ALGORTHIM */
    	INT index=0,j=0,i=0, alpha=0 , beta=0;
    	while ( i < n )
	{ // loop to on LCP table //all maximal sets of indices such that the (lcp) between any two of them is at least L (block)
    		if (LCP[i]>=L)
		{
			alpha=i;                          //store the start index in of the set in alpha
			j=i+1;                          // move j to next pos after start pos has been found
			while ( j < n && LCP[j]>= L)
			{            //move j counter until find LCP < L
				  // next index after i is >=
				j++;
			}
			beta=(j-1);
		
			//************************Processing each set between Alpha-1 and Beta ***********************//
			for(INT i=alpha-1;i<=beta;i++)
			{
				if (  SA[i] % L==0)        
				{                             //check if it is a starting position of a block.
					for (INT j=alpha-1; j<= beta;j++)
					{
						if (i!=j)
						{                                      // to avoid comparing the same pair
											//**finding poistions r_1 and r_2**/
							l = min ( i ,  j );
							r = max ( i ,  j);
							r_1 = LCP[rmq ( l+1 , r ) ]+SA[i];   //LCE to the right
							if(r_1<n-1 && (SA[j]+r_1-SA[i])<n-1)
							{                      //to check not reaching the end of the string
								l = min ( invSA[r_1+1], invSA[SA[j]+r_1-SA[i]+1]);
								r = max ( invSA[r_1+1], invSA[SA[j]+r_1-SA[i]+1]);
								r_2 = LCP[rmq ( l+1,r )] + r_1+1;
							}
							else                                                //check if it is not the last postion in the string
								r_2=r_1+1;
				
							//LCE to the Left
							u_1=n-SA[i]-1;     //find the poistion i in the reverse string
							u_2=n-SA[j]-1;    // find the poistion j in the reverse string
							l = min ( invSA_prime[u_1] ,  invSA_prime[u_2] );  // poistion in suffix array or LCP not the value of suffix array or LCP
							r = max ( invSA_prime[u_1] ,  invSA_prime[u_2]);  // poistion in suffix array or LCP not the value of suffix array or LCp
							L_1 = LCP_prime[rmq_prime( l+1 , r ) ]+u_1;
							l_1=n-L_1-1;                                     // find l_1 in the original string
				
							if(L_1<n-1 && (u_2+L_1-u_1)<n-1)
							{              //check if there is one poistion more for  L_2
								l = min ( invSA_prime[L_1+1], invSA_prime[u_2+L_1-u_1+1]);
								r = max ( invSA_prime[L_1+1], invSA_prime[u_2+L_1-u_1+1]);
								L_2 = LCP_prime[rmq_prime ( l+1,r )] + L_1+1;
							}
							else
								L_2=L_1+1;
				
							l_2=n-L_2-1; // find l_2 in the original string
							INT p=0;
							p=max(p,l_2+1);
							for (;p<=l_1;p++)
							{
								if (p+m-1>=SA[i]+L-1 && p+m-1< r_1 && p+m-1<n && p+SA[j]-SA[i]+m-1<n && p+SA[j]-SA[i]>=0 )
								{ // poistions between l_2 and l_1 and p+m-1 < r_1
									B_i=Counting_Blocks( p, L, l_1,  m);   // count the number of blocks
									B_j=Counting_Blocks( p+SA[j]-SA[i], L, l_1+SA[j]-SA[i],  m);   // count the number of blocks
									if ((B_i+B_j)!=0 )
									{
										C[p]+=1.0/(B_i+B_j);
										C[p+SA[j]-SA[i]]+=1.0/(B_i+B_j);
									}
								}
							}
			     
							//end poistions p
							B_j=0;
							B_i=0;
							INT q=0;
							q=max(q,l_1+1);
							r_2=min(r_2,n);
							for (;q<=SA[i];q++)
							{
								if (q+m-1>=r_1 && q+m-1 <r_2 && q+m-1<n && q+SA[j]-SA[i]+m-1<n && q+SA[j]-SA[i]>=0 )
								{
									B_i=Counting_Blocks( q, L, r_1,  m);
									B_j=Counting_Blocks( q+SA[j]-SA[i], L, r_1+SA[j]-SA[i],  m);
									if ((B_i+B_j)!=0)
									{
										C[q]+=1.0/(B_i+B_j);
										C[q+SA[j]-SA[i]]+=1.0/(B_i+B_j);
									}
								}
							}
			  
						}
						//end poistions q//
					}
				}
			}
			i=j;
    		}//end if
        	else
            		i++;
	}//end loop
    
        FILE * out_fd;
	if ( ! ( out_fd = fopen ( sw . output_filename, "a") ) )
	{
		fprintf ( stderr, " Error: Cannot open file %s!\n", sw . output_filename );
		return ( 1 );
        }
	fprintf ( out_fd, ">%s\n", ( char * ) seq_id );
	for ( INT i = 0; i < n; i++ )
		fprintf ( out_fd, "%4.00f ", C[i] );
	fprintf( out_fd, "\n" );

	if ( fclose ( out_fd ) )
	{
		fprintf( stderr, " Error: file close error!\n");
		return ( 1 );
	}
    
	free ( SA );
	free ( LCP );
    	free ( invSA );
    	util::clear(v);
    	util::clear(v_prime);
    	free(seq_prime);
    	free ( SA_prime );
    	free ( LCP_prime );
    	free ( invSA_prime );
	return ( 1 );
}

unsigned int compute_map_simple ( unsigned char * seq,unsigned char * seq_id, struct TSwitch sw, INT m ,   double  * C )
{
    
    	unsigned char * seq_prime;
	INT * SA;
	INT * LCP;
	INT * invSA;
    	INT * SA_prime;
    	INT * LCP_prime;
    	INT * invSA_prime;

	INT n = strlen ( ( char * ) seq );
    	INT N = m + n;
    	INT L= floor (m/3.0);
    	INT B=0,w=0,B_i=0,B_j=0;

    	INT l,r,l_prime,r_prime ,r_1,r_2,l_1,l_2,u_1,u_2,L_1,L_2;
    	seq_prime = ( unsigned char * ) calloc ( ( n + 1 ) , sizeof ( unsigned char ) );
    
    	for( INT i=n;i>0;i--)
        	seq_prime[w++]=seq[i-1];
	
        /* Compute the suffix array and LCP array for the seq and inverse seq */
        SA = ( INT * ) malloc( ( n ) * sizeof( INT ) );
        SA_prime = ( INT * ) malloc( ( n ) * sizeof( INT ) );
        if( ( SA == NULL) )
        {
                fprintf(stderr, " Error: Cannot allocate memory for SA.\n" );
                return ( 0 );
        }
        if( ( SA_prime == NULL) )
        {
        fprintf(stderr, " Error: Cannot allocate memory for SA.\n" );
        return ( 0 );
        }

	#ifdef _USE_64
        if( divsufsort64( seq, SA,  n ) != 0 )
        {
                fprintf(stderr, " Error: SA computation failed.\n" );
                exit( EXIT_FAILURE );
        }
	#endif

	#ifdef _USE_32
        if( divsufsort( seq, SA,  n ) != 0 )
        {
                fprintf(stderr, " Error: SA computation failed.\n" );
                exit( EXIT_FAILURE );
        }
	#endif

	#ifdef _USE_64
    	if( divsufsort64( seq_prime, SA_prime,  n ) != 0 )
    	{
        	fprintf(stderr, " Error: SA computation failed.\n" );
        	exit( EXIT_FAILURE );
    	}
	#endif
    
	#ifdef _USE_32
    	if( divsufsort( seq_prime, SA_prime,  n ) != 0 )
    	{
       		fprintf(stderr, " Error: SA computation failed.\n" );
       	 	exit( EXIT_FAILURE );
    	}
	#endif
    
        invSA = ( INT * ) calloc( n , sizeof( INT ) );
        if( ( invSA == NULL) )
        {
                fprintf(stderr, " Error: Cannot allocate memory for invSA.\n" );
                return ( 0 );
        }

        for ( INT i = 0; i < n; i ++ )
        {
                invSA [SA[i]] = i;
        }
    
 
    	invSA_prime = ( INT * ) calloc( n , sizeof( INT ) );
    	if( ( invSA_prime == NULL) )
    	{
        	fprintf(stderr, " Error: Cannot allocate memory for invSA.\n" );
        	return ( 0 );
    	}
    
    	for ( INT i = 0; i < n; i ++ )
    	{
        	invSA_prime [SA_prime[i]] = i;
    	}
   

	LCP = ( INT * ) calloc  ( n, sizeof( INT ) );
        if( ( LCP == NULL) )
        {
                fprintf(stderr, " Error: Cannot allocate memory for LCP.\n" );
                return ( 0 );
        }

        /* Compute the LCP array */
        if( LCParray( seq, n, SA, invSA, LCP ) != 1 )
        {
                fprintf(stderr, " Error: LCP computation failed.\n" );
                exit( EXIT_FAILURE );
        }
    
    	//LCP REVERSED String
    	LCP_prime = ( INT * ) calloc  ( n, sizeof( INT ) );
    	if( ( LCP_prime == NULL) )
    	{
        	fprintf(stderr, " Error: Cannot allocate memory for LCP.\n" );
        	return ( 0 );
    	}

    	if( LCParray( seq_prime, n, SA_prime, invSA_prime, LCP_prime ) != 1 )
    	{
        	fprintf(stderr, " Error: LCP computation failed.\n" );
        	exit( EXIT_FAILURE );
    	}
    
	/*1-MAPP ALGORTHIM */
    	INT index=0,j=0,i=0, alpha=0 , beta=0;
    	while ( i < n )
	{ // loop to on LCP table //all maximal sets of indices such that the (lcp) between any two of them is at least L (block)
    		if (LCP[i]>=L)
		{
			alpha=i;                          //store the start index in of the set in alpha
			j=i+1;                          // move j to next pos after start pos has been found
			while ( j < n && LCP[j]>= L)
			{            //move j counter until find LCP < L
				  // next index after i is >=
				j++;
			}
			beta=(j-1);
		
			//************************Processing each set between Alpha-1 and Beta ***********************//
			for(INT i=alpha-1;i<=beta;i++)
			{
				if (  SA[i] % L==0)        
				{                             //check if it is a starting position of a block.
					for (INT j=alpha-1; j<= beta;j++)
					{
						if (i!=j)
						{                                      // to avoid comparing the same pair
											//**finding poistions r_1 and r_2**/
							l = min ( i ,  j );
							r = max ( i ,  j );		//LCP(l+1,r)
							INT lce = 0; INT lSA = SA[l]; INT rSA = SA[r];
  							while ( lSA < n && rSA < n && seq[lSA++] == seq[rSA++]  )
								lce++;
							r_1 = lce + SA[i];   //LCE to the right
							if(r_1<n-1 && (SA[j]+r_1-SA[i])<n-1)
							{                      //to check not reaching the end of the string
								l = min ( invSA[r_1+1], invSA[SA[j]+r_1-SA[i]+1]);
								r = max ( invSA[r_1+1], invSA[SA[j]+r_1-SA[i]+1]);
								INT lce = 0; INT lSA = SA[l]; INT rSA = SA[r];	//LCP(l+1,r)
  								while ( lSA < n && rSA < n && seq[lSA++] == seq[rSA++]  )
									lce++;
								r_2 = lce + r_1 + 1;
							}
							else                                                //check if it is not the last postion in the string
								r_2=r_1+1;
				
							//LCE to the Left
							u_1=n-SA[i]-1;     //find the poistion i in the reverse string
							u_2=n-SA[j]-1;    // find the poistion j in the reverse string
							l = min ( invSA_prime[u_1] ,  invSA_prime[u_2] );  // poistion in suffix array or LCP not the value of suffix array or LCP
							r = max ( invSA_prime[u_1] ,  invSA_prime[u_2]);  // poistion in suffix array or LCP not the value of suffix array or LCp
							lce = 0; lSA = SA_prime[l]; rSA = SA_prime[r];
  							while ( lSA < n && rSA < n  && seq_prime[lSA++] == seq_prime[rSA++] )
								lce++;
							L_1 = lce + u_1;
							l_1=n-L_1-1;                                     // find l_1 in the original string
				
							if(L_1<n-1 && (u_2+L_1-u_1)<n-1)
							{              //check if there is one poistion more for  L_2
								l = min ( invSA_prime[L_1+1], invSA_prime[u_2+L_1-u_1+1]);
								r = max ( invSA_prime[L_1+1], invSA_prime[u_2+L_1-u_1+1]);
								INT lce = 0; INT lSA = SA_prime[l]; INT rSA = SA_prime[r];	//LCP(l+1,r)
  								while ( lSA < n && rSA < n  && seq_prime[lSA++] == seq_prime[rSA++] )
									lce++;
								L_2 = lce + L_1 + 1;
							}
							else
								L_2=L_1+1;
				
							l_2=n-L_2-1; // find l_2 in the original string
                            
                            
                          
                            
                            
							INT p=0;
							p=max(p,l_2+1);
							for (;p<=l_1;p++)
							{
								if (p+m-1>=SA[i]+L-1 && p+m-1< r_1 && p+m-1<n && p+SA[j]-SA[i]+m-1<n && p+SA[j]-SA[i]>=0 )
								{ // poistions between l_2 and l_1 and p+m-1 < r_1
									B_i=Counting_Blocks( p, L, l_1,  m);   // count the number of blocks
									B_j=Counting_Blocks( p+SA[j]-SA[i], L, l_1+SA[j]-SA[i],  m);   // count the number of blocks
									if ((B_i+B_j)!=0 )
									{
										C[p]+=1.0/(B_i+B_j);
										C[p+SA[j]-SA[i]]+=1.0/(B_i+B_j);
									}
								}
							}
			     
							//end poistions p
							B_j=0;
							B_i=0;
							INT q=0;
							q=max(q,l_1+1);
							r_2=min(r_2,n);
							for (;q<=SA[i];q++)
							{
								if (q+m-1>=r_1 && q+m-1 <r_2 && q+m-1<n && q+SA[j]-SA[i]+m-1<n && q+SA[j]-SA[i]>=0 )
								{
									B_i=Counting_Blocks( q, L, r_1,  m);
									B_j=Counting_Blocks( q+SA[j]-SA[i], L, r_1+SA[j]-SA[i],  m);
									if ((B_i+B_j)!=0)
									{
										C[q]+=1.0/(B_i+B_j);
										C[q+SA[j]-SA[i]]+=1.0/(B_i+B_j);
									}
								}
							}
			  
						}
						//end poistions q//
					}
				}
			}
			i=j;
    		}//end if
        	else
            		i++;
	}//end loop
    
 FILE * out_fd;
    if ( ! ( out_fd = fopen ( sw . output_filename, "a") ) )
    {
        fprintf ( stderr, " Error: Cannot open file %s!\n", sw . output_filename );
        return ( 1 );
    }
    fprintf ( out_fd, ">%s\n", ( char * ) seq_id );
    for ( INT i = 0; i < n; i++ )
        fprintf ( out_fd, "%4.00f ", C[i] );
    fprintf( out_fd, "\n" );
    
    if ( fclose ( out_fd ) )
    {
        fprintf( stderr, " Error: file close error!\n");
        return ( 1 );
    }
    
    

    

    
	free ( SA );
	free ( LCP );
    	free ( invSA );
    	free(seq_prime);
    	free ( SA_prime );
    	free ( LCP_prime );
    	free ( invSA_prime );
	return ( 1 );
}




unsigned int compute_At_Most_One_map_simple ( unsigned char * seq,unsigned char * seq_id, struct TSwitch sw, INT m ,   double  * C )
{
    
  
    INT * SA;
    INT * LCP;
    INT * invSA;

    
    INT n = strlen ( ( char * ) seq );
    INT N = m + n;
    INT L= floor (m/3.0);
    INT B=0,w=0,B_i=0,B_j=0;
    
    INT l,r,l_prime,r_prime ,r_1,r_2,l_1,l_2,u_1,u_2,L_1,L_2;

    int h=0;
 
    
    /* Compute the suffix array and LCP array for the seq and inverse seq */
    SA = ( INT * ) malloc( ( n ) * sizeof( INT ) );

    if( ( SA == NULL) )
    {
        fprintf(stderr, " Error: Cannot allocate memory for SA.\n" );
        return ( 0 );
    }
    
#ifdef _USE_64
    if( divsufsort64( seq, SA,  n ) != 0 )
    {
        fprintf(stderr, " Error: SA computation failed.\n" );
        exit( EXIT_FAILURE );
    }
#endif
    
#ifdef _USE_32
    if( divsufsort( seq, SA,  n ) != 0 )
    {
        fprintf(stderr, " Error: SA computation failed.\n" );
        exit( EXIT_FAILURE );
    }
#endif
    
    
    invSA = ( INT * ) calloc( n , sizeof( INT ) );
    if( ( invSA == NULL) )
    {
        fprintf(stderr, " Error: Cannot allocate memory for invSA.\n" );
        return ( 0 );
    }
    
    for ( INT i = 0; i < n; i ++ )
    {
        invSA [SA[i]] = i;
    }
    
    
    
    
    LCP = ( INT * ) calloc  ( n, sizeof( INT ) );
    if( ( LCP == NULL) )
    {
        fprintf(stderr, " Error: Cannot allocate memory for LCP.\n" );
        return ( 0 );
    }
    
    /* Compute the LCP array */
    if( LCParray( seq, n, SA, invSA, LCP ) != 1 )
    {
        fprintf(stderr, " Error: LCP computation failed.\n" );
        exit( EXIT_FAILURE );
    }
    
    
    
    
    
    /*0-Mapp Algorthim*/
    
    
 compute_ZeroMapp (n,  m , SA,LCP, invSA, C );

    
    /*1-MAPP ALGORTHIM */
    INT index=0,j=0,i=0, alpha=0 , beta=0;
    while ( i < n )
    {                                            // loop to on LCP table //all maximal sets of indices such that the (lcp) between any two of them is at least L (block)
        if (LCP[i]>=L)
        {
            alpha=i;                               //store the start index in of the set in alpha
            j=i+1;                                // move j to next pos after start pos has been found
            while ( j < n && LCP[j]>= L)
            {                                       //move j counter until find LCP < L
                                                   // next index after i is >=
                j++;
            }
            beta=(j-1);
            
            //************************Processing each set between Alpha-1 and Beta ***********************//
            for(INT i=alpha-1;i<=beta;i++)
            {
                if (  SA[i] % L==0)
                {                             //check if it is a starting position of a block.
                    for (INT j=alpha-1; j<= beta;j++)
                    {
                        if (i!=j)
                        {
                        // cout <<"i: "<<i <<" j: "<<j<<endl;// to avoid comparing the same pair
                            //**finding poistions r_1 and r_2**/
                            l = min ( i ,  j );
                            r = max ( i ,  j );		//LCP(l+1,r)
                            INT lce = 0; INT lSA = SA[l]; INT rSA = SA[r];
                            while ( lSA < n && rSA < n && seq[lSA++] == seq[rSA++]  )
                                lce++;
                            r_1 = lce + SA[i];   //LCE to the right
                            
                            
                            
                            
                            if(r_1<n-1 && (SA[j]+r_1-SA[i])<n-1)
                            {                                                       //to check not reaching the end of the string
                                l = min ( invSA[r_1+1], invSA[SA[j]+r_1-SA[i]+1]);
                                r = max ( invSA[r_1+1], invSA[SA[j]+r_1-SA[i]+1]);
                                
                                
                                INT lce = 0;
                                INT lSA = SA[l];
                                INT rSA = SA[r];
                                
                                while ( lSA < n && rSA < n && seq[lSA++] == seq[rSA++]  )
                                    lce++;
                                
                                r_2 = lce + r_1 + 1;
                            }
                            else                                                //check if it is not the last postion in the string
                                r_2=r_1+1;
                            
                            
                            //**finding poistions l_1 and l_2**/
                            l = min ( i ,  j );
                            r = max ( i ,  j );
                            
                            
                             lce = 0;
                             lSA = SA[l];
                             rSA = SA[r];
                            
                            
                            while ( lSA >= 0 && rSA >= 0 && seq[lSA--] == seq[rSA--]  )
                                lce++;
                            
                            
                            l_1 =  SA[i] - lce ;   //LCE to the left

                            if(l_1>0 && (abs(SA[i]-l_1-SA[j]))>0)  // to check there is enough pos berore reaching 0 to serach for L_2
                            {
                               INT lce = 0;
                                lSA= l_1-1; // to start from pos l_1-1 before i immeditly
                                rSA =SA[j]- abs(SA[i]-l_1)-1; // to start from the pos j-l_1-1
                           
                             
                                while ( lSA >= 0 && seq[lSA--] == seq[rSA--]  )
                                    lce++;
                                
                                l_2 = l_1 - lce -1 ;
                            }
                            
                            else                                                //check if it is not the last postion in the string
                                l_2=l_1-1;
                            

                            INT p=0;
                            p=max(p,l_2+1);
                
                            for (;p<=l_1;p++)
                            {
                           
                                if (p+m-1>=SA[i]+L-1 && p+m-1< r_1 && p+m-1<n && p+SA[j]-SA[i]+m-1<n && p+SA[j]-SA[i]>=0 )
                                { // poistions between l_2 and l_1 and p+m-1 < r_1
                                   
                                 /// cout << "p  S[I] "<<p<<endl;
                                   // cout << "mismatch "<< l_1<<endl;
                                    //cout << "blocks "<< B_i<<endl;
                                    B_i=Counting_Blocks( p, L, l_1,  m);   // count the number of blocks
                                    B_j=Counting_Blocks( p+SA[j]-SA[i], L, l_1+SA[j]-SA[i],  m);   // count the number of blocks
                                        //  cout <<"blocks "<<  B_i+      B_j<<endl;
                                    if ((B_i+B_j)!=0 )
                                    {
                                       
                                        C[p]+=1.0/(B_i+B_j);
                                        C[p+SA[j]-SA[i]]+=1.0/(B_i+B_j);
                                        
                              
                                    }
                       
                                }
                            }
 
                            //end poistions p
                            B_j=0;
                            B_i=0;
                            INT q=0;

                            q=max(q,l_1+1);
                            r_2=min(r_2,n);

                            for (;q<=SA[i];q++)
                            {

                                if (q+m-1>=r_1 && q+m-1 <r_2 && q+m-1<n && q+SA[j]-SA[i]+m-1<n && q+SA[j]-SA[i]>=0 )
                                {
                              
                            
                                // cout << "q  S[I] "<<q<<endl;
                                  //  cout << "mismatch "<< r_1<<endl;
                   
                                    B_i=Counting_Blocks( q, L, r_1,  m);
                                    B_j=Counting_Blocks( q+SA[j]-SA[i], L, r_1+SA[j]-SA[i],  m);
                                           // cout <<"blocks "<<  B_i+      B_j<<endl;
                                    if ((B_i+B_j)!=0)
                                    {
                                   
                                        C[q]+=1.0/(B_i+B_j);
                                        C[q+SA[j]-SA[i]]+=1.0/(B_i+B_j);
                                        
                                
                                    }
                                }
                            }
                            
                        }
                        //end poistions q//
                    }
                }
            }
            i=j;
        }//end if
        else
            i++;
    }//end loop

    FILE * out_fd;
    if ( ! ( out_fd = fopen ( sw . output_filename, "a") ) )
    {
        fprintf ( stderr, " Error: Cannot open file %s!\n", sw . output_filename );
        return ( 1 );
    }
    fprintf ( out_fd, ">%s\n", ( char * ) seq_id );
    for ( INT i = 0; i < n; i++ )
        fprintf ( out_fd, "%4.00f ", C[i] );
    fprintf( out_fd, "\n" );
    
    if ( fclose ( out_fd ) )
    {
        fprintf( stderr, " Error: file close error!\n");
        return ( 1 );
    }
    
    
    
   
    
    
    
    free ( SA );
    free ( LCP );
    free ( invSA );

   
    return ( 1 );
}


unsigned int compute_naive_map ( unsigned char * seq, unsigned char * seq_id, struct TSwitch sw, double * C, INT m , unsigned int k )
{
 
    unsigned char * substring_1=NULL;
    unsigned char * substring_2=NULL;

    substring_1 = ( unsigned char * ) realloc  ( substring_1,   ( m + ALLOC_SIZE ) * sizeof ( unsigned char ) );
    substring_2 = ( unsigned char * ) realloc  ( substring_2,   ( m + ALLOC_SIZE ) * sizeof ( unsigned char ) );
    INT n = strlen ( ( char * ) seq );
    
    for (INT i=0;i<n;i++){
        
     memmove(substring_1,seq+i,m); //subsrting with size m

        
        for (INT j=0;j<n;j++){
            if (i!=j){  // to avoid comparing the same subsrting in case k !=0 will not be noticed
            memmove(substring_2,seq+j,m);// substring with size mm
                if (strlen ( ( char * ) substring_2 )== m && strlen ( ( char * ) substring_1 )==m){
                if (Hamming_Dist(substring_1,substring_2, m)==k){
                    C[i]++;
                }//end if
                    
                }
            }
        }//end loop
    }//end loop


        FILE * out_fd;
	if ( ! ( out_fd = fopen ( sw . output_filename, "a") ) )
	{
		fprintf ( stderr, " Error: Cannot open file %s!\n", sw . output_filename );
		return ( 1 );
        }
	fprintf ( out_fd, ">%s\n", ( char * ) seq_id );
	for ( INT i = 0; i < n; i++ )
		fprintf ( out_fd, "%4.00f ", C[i] );
	fprintf( out_fd, "\n" );


	if ( fclose ( out_fd ) )
	{
		fprintf( stderr, " Error: file close error!\n");
		return ( 1 );
	}

     free(substring_1);
     free (substring_2);
        
}



unsigned int compute_naive_at_most_k_map ( unsigned char * seq, unsigned char * seq_id, struct TSwitch sw, double * C, INT m , unsigned int k )
{
  
    unsigned char * substring_1=NULL;
    unsigned char * substring_2=NULL;

    substring_1 = ( unsigned char * ) realloc  ( substring_1,   ( m + ALLOC_SIZE ) * sizeof ( unsigned char ) );
    substring_2 = ( unsigned char * ) realloc  ( substring_2,   ( m + ALLOC_SIZE ) * sizeof ( unsigned char ) );
    INT n = strlen ( ( char * ) seq );
    
    for (INT i=0;i<n;i++){
        
        memmove(substring_1,seq+i,m); //subsrting with size m
        cout << "substring"<<endl;
           cout << substring_1 <<endl;
        for (INT j=0;j<n;j++){
            
            if (i !=j){
    
             
            memmove(substring_2,seq+j,m);// substring with size m
            if (strlen ( ( char * ) substring_2 )== m && strlen ( ( char * ) substring_1 )==m){
                if (Hamming_Dist(substring_1,substring_2, m) <= k){
                    cout << substring_2 <<endl;
                    C[i]++;
                }//end if
                
            }// to avoid counting the same substrings
            }
        }//end loop
    }//end loop
    
    
    FILE * out_fd;
    if ( ! ( out_fd = fopen ( sw . output_filename, "a") ) )
    {
        fprintf ( stderr, " Error: Cannot open file %s!\n", sw . output_filename );
        return ( 1 );
    }
    fprintf ( out_fd, ">%s\n", ( char * ) seq_id );
    for ( INT i = 0; i < n; i++ )
        fprintf ( out_fd, "%4.00f ", C[i] );
    fprintf( out_fd, "\n" );
    
    
    if ( fclose ( out_fd ) )
    {
        fprintf( stderr, " Error: file close error!\n");
        return ( 1 );
    }
    
    free(substring_1);
    free (substring_2);
    
}







unsigned int Hamming_Dist ( unsigned char * sub_1, unsigned char * sub_2, INT m ){
    unsigned int dist=0;
    for (int i=0;i<m;i++){
        if (sub_1[i]!=sub_2[i])
            dist++;
        
    }
    return dist;
}







unsigned int compute_ZeroMapp (INT n, INT m ,  INT * SA, INT  * LCP, 	INT * invSA, double * C )
{
   

    /*0-MAPP ALGORTHIM */

    
    INT i=0, Counter=0,alpha=0,beta=0;
    bool newCluster=false;
    
    while ( i < n )
    {
        
        if (LCP[i]>=m && newCluster==false){
            newCluster=true;
        
            alpha=i-1;// save the sart poistion of  a new cluster
            
            
        }
        
        if (LCP[i]>=m)
        {
            
            Counter++; //calculating the size of the cluster
   
        
        }
        
        
        
        i++;
        
        
        if (LCP[i]<m && newCluster==true ){ // save the end poistion of the cluster
        beta=i-1;
        

        
        for (int j=alpha; j<=beta;j++) //update the array C with size of the clusters
        {
            C[SA[j]]= Counter; //+1 becuse the the poistion X-1 in the cluster
            
           
            
        }
              newCluster=false; // reset the value to start new cluster
            Counter=0;
            
        }//end if
        
        
        
      
    
    }//end loop



    
    

    return ( 1 );
}


unsigned int compute_At_Most_Two_map_simple ( unsigned char * seq,unsigned char * seq_id, struct TSwitch sw, INT m ,   double  * C )
{
    
    
    INT * SA;
    INT * LCP;
    INT * invSA;
    
    
    INT n = strlen ( ( char * ) seq );
    INT N = m + n;

    INT L= floor (m/4);
    INT B=0,w=0,B_i=0,B_j=0;
    
    INT l,r,l_prime,r_prime ,r_1,r_2,r_3,l_1,l_2,l_3,u_1,u_2,L_1,L_2;
    
    
    
    
    /* Compute the suffix array and LCP array for the seq and inverse seq */
    SA = ( INT * ) malloc( ( n ) * sizeof( INT ) );
    
    if( ( SA == NULL) )
    {
        fprintf(stderr, " Error: Cannot allocate memory for SA.\n" );
        return ( 0 );
    }
    
#ifdef _USE_64
    if( divsufsort64( seq, SA,  n ) != 0 )
    {
        fprintf(stderr, " Error: SA computation failed.\n" );
        exit( EXIT_FAILURE );
    }
#endif
    
#ifdef _USE_32
    if( divsufsort( seq, SA,  n ) != 0 )
    {
        fprintf(stderr, " Error: SA computation failed.\n" );
        exit( EXIT_FAILURE );
    }
#endif
    
    
    invSA = ( INT * ) calloc( n , sizeof( INT ) );
    if( ( invSA == NULL) )
    {
        fprintf(stderr, " Error: Cannot allocate memory for invSA.\n" );
        return ( 0 );
    }
    
    for ( INT i = 0; i < n; i ++ )
    {
        invSA [SA[i]] = i;
    }
    
    
    
    
    LCP = ( INT * ) calloc  ( n, sizeof( INT ) );
    if( ( LCP == NULL) )
    {
        fprintf(stderr, " Error: Cannot allocate memory for LCP.\n" );
        return ( 0 );
    }
    
    /* Compute the LCP array */
    if( LCParray( seq, n, SA, invSA, LCP ) != 1 )
    {
        fprintf(stderr, " Error: LCP computation failed.\n" );
        exit( EXIT_FAILURE );
    }
    INT index=0,j=0,i=0, alpha=0 , beta=0;
    while ( i < n )
    {                                            // loop to on LCP table //all maximal sets of indices such that the (lcp) between any two of them is at least L (block)
        if (LCP[i]>=L)
        {
            
            
            alpha=i;                               //store the start index in of the set in alpha
            j=i+1;                                // move j to next pos after start pos has been found
            while ( j < n && LCP[j]>= L)
            {                                       //move j counter until find LCP < L
                // next index after i is >=
                j++;
            }
            beta=(j-1);
            

            
            //************************Processing each set between Alpha-1 and Beta ***********************//
            for(INT i=alpha-1;i<=beta;i++)
            {
             
                if (  SA[i] % L==0)
                {                         //check if it is a starting position of a block.
                    for (INT j=alpha-1; j<= beta;j++)
                    {
                     
                        if (i!=j)
                        {
                        
                           // cout << "i: "<<i<<"j: "<<j<<endl;
                            // to avoid comparing the same pair
                            //**finding poistions r_1 and r_2**/
                            l = min ( i ,  j );
                            r = max ( i ,  j );		//LCP(l+1,r)
                            INT lce = 0; INT lSA = SA[l]; INT rSA = SA[r];
                            while ( lSA < n && rSA < n && seq[lSA++] == seq[rSA++]  )
                                lce++;
                            r_1 = lce + SA[i];   //LCE to the right
                            
                            
                            
                            //Finding r_2
                            if(r_1<n-1 && (SA[j]+r_1-SA[i])<n-1)
                            {                                                       //to check not reaching the end of the string
                                l = min ( invSA[r_1+1], invSA[SA[j]+r_1-SA[i]+1]);
                                r = max ( invSA[r_1+1], invSA[SA[j]+r_1-SA[i]+1]);
                                
                                
                                INT lce = 0;
                                INT lSA = SA[l];
                                INT rSA = SA[r];
                                
                                while ( lSA < n && rSA < n && seq[lSA++] == seq[rSA++]  )
                                    lce++;
                                
                                r_2 = lce + r_1 + 1;
                            }
                            else
                                r_2=r_1+1;
                            
                            //Finding r_3
                            
                            
                            if(r_2<n-1 && (SA[j]+r_2-SA[i])<n-1)
                            {                                                       //to check not reaching the end of the string
                                l = min ( invSA[r_2+1], invSA[SA[j]+r_2-SA[i]+1]);
                                r = max ( invSA[r_2+1], invSA[SA[j]+r_2-SA[i]+1]);
                                
                                
                                INT lce = 0;
                                INT lSA = SA[l];
                                INT rSA = SA[r];
                                
                                while ( lSA < n && rSA < n && seq[lSA++] == seq[rSA++]  )
                                    lce++;
                                
                                r_3 = lce + r_2 + 1;
                            }
                            else
                                r_3=r_2+1;
                            
                            
                            
                            
                 
                            
                            //********************************finding poistions l_1, l_2 and l_3***********************************************/
                            l = min ( i ,  j );
                            r = max ( i ,  j );
                            
                            
                            lce = 0;
                            lSA = SA[l];
                            rSA = SA[r];
                            
                            //l_1
                            while ( lSA >= 0 && rSA >= 0 && seq[lSA--] == seq[rSA--]  )
                                lce++;
                            
                            
                            l_1 =  SA[i] - lce ;   //LCE to the left
                            
                          
                            
                            
                            
                            //l_2
                            
                            if(l_1>0 && (abs(SA[i]-l_1-SA[j]))>0)  // to check there is enough pos berore reaching 0 to serach for L_2
                            {
                                
                                
                                
                                INT lce = 0;
                                lSA= l_1-1; // to start from pos l_1-1 before i immeditly
                                rSA =SA[j]- abs(SA[i]-l_1)-1; // to start from the pos j-l_1-1
                                
                                
                                while ( lSA >= 0 && seq[lSA--] == seq[rSA--]  )
                                    lce++;
                                
                                l_2 = l_1 - lce -1 ;
                            }
                            
                            else
                                l_2=l_1-1;
                            
                            
                            //l_3
                        
                            
                            if(l_2>0 && (abs(SA[i]-l_2-SA[j]))>0)  // to check there is enough pos berore reaching 0 to serach for L_2
                            {
                                
                                
                                
                                INT lce = 0;
                                
                                
                                
                                lSA= l_2-1; // to start from pos l_1-1 before i immeditly
                                rSA =SA[j]- abs(SA[i]-l_2)-1; // to start from the pos j-l_1-1
                                
                                
                                while ( lSA >= 0 && seq[lSA--] == seq[rSA--]  )
                                    lce++;
                                
                                l_3 = l_2 - lce -1 ;
                            }
                            
                            else
                                l_3=l_2-1;
                            
        
                          //  cout << r_1<< " "<< r_2<<" "<< r_3<<endl;
                           // cout << l_1<< " "<< l_2<<" "<< l_3<<endl;
                               /*****************processing each k's pos******************************************/
                              /**************Only pos with k=2******************************************************/
                   
                            //cout << "(i: "<< i<< " j: "<<j<<")"<<endl;
                            //cout << "r_1 "<<r_1<<" r_2 "<<r_2<<" r_3 "<<r_3<<endl;
                            //cout <<"l_1 "<< l_1<<" l_2 "<<l_2<<" l_3 "<<l_3<<endl;
                            INT k=0, u=0;
                          
                      
                                  k=max(k,l_3+1);
//l_3+=abs((SA[i]+L-1- m)-l_3); // to remove unwanted poistions in l_3
                         
                 //    cout <<"region 3"<<endl;
                     
                            for (;k<=l_2;k++)
                            {
                 
                         
                              
      
                                if ((k+m-1)>=(SA[i]+L-1) && (k+m-1)< r_1 && (k+m-1)<n && (k+SA[j]-SA[i]+m-1)<n && (k+SA[j]-SA[i])>=0 )
                                { // poistions between l_3 and l_2 and p+m-1 < r_1
                                   // cout << "pos  region 3 "<<k<<endl;
                                    //if ( k==7 ||(k+SA[j]-SA[i])==7 ){cout <<"loool"<<endl;}
                                    B_i=Counting_Blocks_for_Two_Mismatch( k, L, l_2,l_1,  m);
                                    B_j=Counting_Blocks_for_Two_Mismatch( k+SA[j]-SA[i], L,l_2+SA[j]-SA[i] ,l_1+SA[j]-SA[i],  m);
                     
                                    if ((B_i+B_j)!=0 )
                                    {
                                        
                                    
                 
                    
                                       
                                        C[k]+=1.0/(B_i+B_j);
                                        C[k+SA[j]-SA[i]]+=1.0/(B_i+B_j);
                               
                                    }
                                }
                           
                            }
                            
                            
                           /*****************processing each p's pos******************************************/
                          /**************Only pos with k=1 or k=2******************************************************/
                            B_j=0;
                            B_i=0;
                            INT p=0;
                          
                            p=max(p,l_2+1);
                       //  cout <<"region 2"<<endl;

                          l_2+=abs((SA[i]+L-1- m)-l_2); // to remove umwanted poistions in l_2
                           
                            for (;p<=l_1;p++)
                            {
                                //cout << "Region 2"<<endl;
                                //cout <<"l_pos[o-1] "<< l_2<<endl;
                                //cout << "(SA[i]+L-1)"<< (SA[i]+L-1)<<endl;
                                //cout <<"r_pos[k+1-o+1] "<< r_2<<endl;
                                if ((p+m-1)>=(SA[i]+L-1) && (p+m-1)< r_2 && (p+m-1)<n && (p+SA[j]-SA[i]+m-1)<n && (p+SA[j]-SA[i])>=0 )
                                
                                { // poistions between l_2 and l_1 and p+m-1 < r_1
                                      //  cout << "pos region 2 "<<p<<endl;
                                    //if ( p==7 ||(p+SA[j]-SA[i])==7 ){
                                      //  cout <<"loool"<<endl;}
                                    B_i=Counting_Blocks_for_Two_Mismatch( p, L, l_1,r_1,  m);   // count the number of blocks
                                    B_j=Counting_Blocks_for_Two_Mismatch( (p+SA[j]-SA[i]), L, (l_1+SA[j]-SA[i]),(r_1+SA[j]-SA[i]),  m);   // count the number of blocks
                                    
                          
                                    if ((B_i+B_j)!=0 )
                                    {
                                        C[p]+=1.0/(B_i+B_j);
                                        C[p+SA[j]-SA[i]]+=1.0/(B_i+B_j);
              
                                    }
                                }
                            }
                            
                            //end poistions p
                            B_j=0;
                            B_i=0;
                            INT q=0;
             
                            q=max(q,l_1+1);
                  
                           l_1+=abs((SA[i]+L-1- m)-l_1); // to remove unwanted poistions in l_1
                            
                            
                                /*****************processing each q's pos******************************************/
                                     /**************Only pos with K=0,k=1 and k=2******************************************************/
                            
                            
                           //    cout <<"region 1"<<endl;
                            for (;q<=SA[i];q++)
                            {
                               // cout << SA[i]<<endl;
                                //cout << q<<endl;
                                //cout <<r_3<<endl;
                                //cout << (SA[i]+L-1)<<endl;
                             if ((q+m-1)>=(SA[i]+L-1) && (q+m-1) <r_3 && (q+m-1)<n && (q+SA[j]-SA[i]+m-1)<n && (q+SA[j]-SA[i])>=0 )
                  
                                {
                                                   //   cout << "pos region 1 "<<q <<endl;
                                   // if ( q==7 ||(q+SA[j]-SA[i])==7){ cout <<"loool"<<endl;}
                                    B_i=Counting_Blocks_for_Two_Mismatch( q, L, r_1,r_2 , m);
                                   
                                    B_j=Counting_Blocks_for_Two_Mismatch( q+SA[j]-SA[i], L, r_1+SA[j]-SA[i],  r_2+SA[j]-SA[i],m);
                                   if ((B_i+B_j)!=0)
                                    {
                                        
                                        C[q]+=1.0/(B_i+B_j);
                                        C[q+SA[j]-SA[i]]+=1.0/(B_i+B_j);
                        
                                        
                                    }
                                }
                            }
                            
                            
                             //end poistions q//
                            
                            
                        }
                    
                        
                        
                        
                        
                        
                        
                    }
                }
            }
            
            
            
            i=j;
        }//end if
        else
            i++;
    }//end loop
    
    

    
    FILE * out_fd;
    if ( ! ( out_fd = fopen ( sw . output_filename, "a") ) )
    {
        fprintf ( stderr, " Error: Cannot open file %s!\n", sw . output_filename );
        return ( 1 );
    }
    fprintf ( out_fd, ">%s\n", ( char * ) seq_id );
    for ( INT i = 0; i < n; i++ )
        fprintf ( out_fd, "%4.00f ", C[i] );
    fprintf( out_fd, "\n" );
    
    if ( fclose ( out_fd ) )
    {
        fprintf( stderr, " Error: file close error!\n");
        return ( 1 );
    }
    
    
    free ( SA );
    free ( LCP );
    free ( invSA );
    

    return ( 1 );
}

unsigned int compute_At_Most_K_map_simple ( unsigned char * seq,unsigned char * seq_id, struct TSwitch sw, INT m ,   double  * C ,int k)
{
    INT L= floor ((m/(k+2.0))*1.0);
    int number_of_mismatch=0;
    INT n = strlen ( ( char * ) seq );
    //INT N = m + n;
    INT * SA;
    INT * LCP;
    INT * invSA;
    
    /* Compute the suffix array and LCP array for the seq and inverse seq */
    SA = ( INT * ) malloc( ( n ) * sizeof( INT ) );
    
    if( ( SA == NULL) )
    {
        fprintf(stderr, " Error: Cannot allocate memory for SA.\n" );
        return ( 0 );
    }
    
#ifdef _USE_64
    if( divsufsort64( seq, SA,  n ) != 0 )
    {
        fprintf(stderr, " Error: SA computation failed.\n" );
        exit( EXIT_FAILURE );
    }
#endif
    
#ifdef _USE_32
    if( divsufsort( seq, SA,  n ) != 0 )
    {
        fprintf(stderr, " Error: SA computation failed.\n" );
        exit( EXIT_FAILURE );
    }
#endif
    
    
    invSA = ( INT * ) calloc( n , sizeof( INT ) );
    if( ( invSA == NULL) )
    {
        fprintf(stderr, " Error: Cannot allocate memory for invSA.\n" );
        return ( 0 );
    }
    
    for ( INT i = 0; i < n; i ++ )
    {
        invSA [SA[i]] = i;
    }
    
    
    
    
    LCP = ( INT * ) calloc  ( n, sizeof( INT ) );
    if( ( LCP == NULL) )
    {
        fprintf(stderr, " Error: Cannot allocate memory for LCP.\n" );
        return ( 0 );
    }
    
    /* Compute the LCP array */
    if( LCParray( seq, n, SA, invSA, LCP ) != 1 )
    {
        fprintf(stderr, " Error: LCP computation failed.\n" );
        exit( EXIT_FAILURE );
    }
    
    INT index=0,j=0,i=0, alpha=0 , beta=0;
    //   omp_set_num_threads(sw . T);

    double t=0.0;
    while ( i < n )
    {                                            // loop to on LCP table //all maximal sets of indices such that the (lcp) between any two of them is at least L (block)
        if (LCP[i]>=L)
        {
            
            
            alpha=i;                               //store the start index in of the set in alpha
            j=i+1;                                // move j to next pos after start pos has been found
            while ( j < n && LCP[j]>= L)
            {                                       //move j counter until find LCP < L
                // next index after i is >=
                j++;
            }
            beta=(j-1);
            //************************Processing each set between Alpha-1 and Beta ***********************//

              // cout << alpha<< " "<<beta<<endl;
            //***************parallel goes here********************///
            
            //   #pragma omp parallel for shared( dic )
            for(INT i=alpha-1;i<=beta;i++)
            {
                
                int tid =0;
              
                
                
                INT * r_pos;
                INT * l_pos;
                INT * mismatches; INT * mismatches_j_pos ;
                
                
                
                INT B_i=0,B_j=0;
                
                INT l,r ,r_1,r_2,r_3,l_1,l_2,l_3,u_1,u_2,L_1,L_2;
                
                
                l_pos = ( INT * ) malloc( ( k+2 ) * sizeof( INT ) );
                r_pos = ( INT * ) malloc( ( k+2 ) * sizeof( INT ) );
                mismatches_j_pos  = ( INT * ) malloc( ( k ) * sizeof( INT ) );
                mismatches = ( INT * ) malloc( ( k ) * sizeof( INT ) );
                
                if (  SA[i] % L==0)
                {                                                                //check if it is a starting position of a block.
                    
                    for (INT j=alpha-1; j<= beta;j++)
                    {
                        
                        
                        if (i!=j) // to avoid comparing the same pair
                        {
                      
                            for(INT o=1;o<=k+1;o++){  //extending to right with number of mismatches+1
                                
                                if (o==1){ //if k <= 1
                                    l = min ( i ,  j );
                                    r = max ( i ,  j );        //LCP(l+1,r)
                                    INT lce = 0; INT lSA = SA[l]; INT rSA = SA[r];
                                    while ( lSA < n && rSA < n && seq[lSA++] == seq[rSA++]  )
                                        lce++;
                                    r_1 = lce + SA[i];   //LCE to the right
                                    //Finding r_2
                                    if(r_1<n-1 && (SA[j]+r_1-SA[i])<n-1)
                                    {                                                       //to check not reaching the end of the string
                                        l = min ( invSA[r_1+1], invSA[SA[j]+r_1-SA[i]+1]);
                                        r = max ( invSA[r_1+1], invSA[SA[j]+r_1-SA[i]+1]);
                                        INT lce = 0;
                                        INT lSA = SA[l];
                                        INT rSA = SA[r];
                                        
                                        while ( lSA < n && rSA < n && seq[lSA++] == seq[rSA++]  )
                                            lce++;
                                        r_2 = lce + r_1 + 1;
                                    }
                                    else
                                        r_2=r_1+1;
                                    
                                    r_pos[o]=r_1;
                                    r_pos[o+1]=r_2;
                                    
                                    o++;
                                } //end if for the r first mismatch
                                
                                else
                                {
                                    
                                    
                                    r_2= r_pos[o-1]; // to genrlaize it to k
                                    if(r_2<n-1 && (SA[j]+r_2-SA[i])<n-1)
                                    {                                                               //to check not reaching the end of the string
                                        l = min ( invSA[r_2+1], invSA[SA[j]+r_2-SA[i]+1]);
                                        r = max ( invSA[r_2+1], invSA[SA[j]+r_2-SA[i]+1]);
                                        INT lce = 0;
                                        INT lSA = SA[l];
                                        INT rSA = SA[r];
                                        
                                        while ( lSA < n && rSA < n && seq[lSA++] == seq[rSA++]  )
                                            lce++;
                                        
                                        r_3 = lce + r_2 + 1;
                                    }
                                    else
                                        r_3=r_2+1;
                                    r_pos[o]=r_3; // to genrlaize it to k
                                }// if k> 1
                                
                            }//end loop extending to right
                            
                            INT lce = 0; INT lSA = SA[l]; INT rSA = SA[r];
                            
                            for(INT  o=1;o<=k+1;o++) //extending to left with number of mismatches+1
                            {
                                if (o==1){ //if k <= 1
                                    l = min ( i ,  j );
                                    r = max ( i ,  j );
                                    lce = 0;
                                    lSA = SA[l];
                                    rSA = SA[r];
                                    while ( lSA >= 0 && rSA >= 0 && seq[lSA--] == seq[rSA--]  )
                                        lce++;
                                    l_1 =  SA[i] - lce ;   //LCE to the left
                                    if(l_1>0 && (abs(SA[i]-l_1-SA[j]))>0)  // to check there is enough pos berore reaching 0 to serach for L_2
                                    {
                                        
                                        lce = 0;
                                        lSA= l_1-1; // to start from pos l_1-1 before i immeditly
                                        rSA =SA[j]- abs(SA[i]-l_1)-1; // to start from the pos j-l_1-1
                                        
                                        
                                        while ( lSA >= 0 && rSA>=0 && seq[lSA--] == seq[rSA--]  )
                                            lce++;
                                        
                                        l_2 = l_1 - lce -1 ;
                                    }
                                    
                                    else
                                        l_2=l_1-1;
                                    
                                    
                                    l_pos[o]= l_1;
                                    l_pos[o+1]= l_2;
                                    o++;
                                }// end if for l_1
                                
                                
                                else { // if k>1
                                    
                                    l_2= l_pos[o-1]; // genralize it to k
                                    
                                    if(l_2>0 && (abs(SA[i]-l_2-SA[j]))>0)  // to check there is enough pos before reaching 0 to serach for L_2
                                    {
                                        
                                        INT lce = 0;
                                        lSA= l_2-1; // to start from pos l_1-1 before i immeditly
                                        rSA =SA[j]- abs(SA[i]-l_2)-1; // to start from the pos j-l_1-1
                                        while ( lSA >= 0 && seq[lSA--] == seq[rSA--]  )
                                            lce++;
                                        l_3 = l_2 - lce -1 ;
                                    }
                                    else{
                                        l_3=l_2-1;
                                        //  cout << l_3<<endl;
                                    }
                                    
                                    
                                    l_pos[o]=l_3; // genralize it to k
                                    
                                }// end else if k >1
                            }//end loop extending to left
                            
                            ///////////////////////////////////////Prossing k+1 regions with k mismatches ////////////////////////////////////////////////////////////////////////////
                            B_j=0;
                            B_i=0;
                            INT q=0;
                            
                            
                            for (INT o=1;o<=k+1;o++){
                                
                                if (o==1){
                                    
                                    INT mismatch = 1;
                                    for (INT counter=0;counter<k;counter++){
                                        mismatches[counter]=r_pos[mismatch];
                                        mismatch++;}
                                    for (int v=0;v <k;v++){
                                        mismatches_j_pos[v]=mismatches[v]+SA[j]-SA[i];}
                                    q=max(q,l_pos[1]+1);
                                    for (;q<=SA[i];q++){
                                        if ((q+m-1)>=(SA[i]+L-1) && (q+m-1) <r_pos[k+1] && (q+m-1)<n && (q+SA[j]-SA[i]+m-1)<n && (q+SA[j]-SA[i])>=0 ){
                                            B_i=Counting_Blocks_for_K_Mismatch( q, L,  mismatches,  m,k );
                                            B_j=Counting_Blocks_for_K_Mismatch( (q+SA[j]-SA[i]), L, mismatches_j_pos,  m, k);
                                            if ((B_i+B_j)!=0){
                                                C[q]+=1.0/(B_i+B_j);
                                                C[ q+SA[j]-SA[i]]+=1.0/(B_i+B_j);
                                            }
                                        }
                                    }
                                }
                                else{ // other region
                                    B_j=0;B_i=0;
                                    INT p=0,v=0;
                                    int mismatch = o-1;
                                    int r_mismatch=0;
                                    for (INT counter=1;counter<=k;counter++){
                                        if (mismatch>=1){
                                            mismatches[v]=l_pos[mismatch];
                                            v++;
                                            mismatch=mismatch-1;
                                            r_mismatch = mismatch;}
                                        else{
                                            r_mismatch++;
                                            mismatches[v]=r_pos[r_mismatch];
                                            v++;}
                                    }//filling array with mismatches
                                    
                                    for (int v=0;v <k;v++){
                                        mismatches_j_pos[v]=mismatches[v]+SA[j]-SA[i];
                                    }
                                    p=max(p,l_pos[o]+1);
                                    for (;p<=l_pos[o-1];p++){
                                        
                                        if ((p+m-1)>=(SA[i]+L-1) && (p+m-1)< r_pos[k+1-o+1] && (p+m-1)<n && (p+SA[j]-SA[i]+m-1)<n && (p+SA[j]-SA[i])>=0 ){
                                            B_i=Counting_Blocks_for_K_Mismatch( p, L,  mismatches,  m,k );
                                            B_j=Counting_Blocks_for_K_Mismatch( (p+SA[j]-SA[i]), L, mismatches_j_pos,  m, k);
                                            if ((B_i+B_j)!=0 )
                                            {
                                                C[p]+=1.0/(B_i+B_j);
                                            C[ p+SA[j]-SA[i]]+=1.0/(B_i+B_j);
                                              
                                              
                                              
                                            }
                                        }
                                    }
                                    
                                    
                                } //end else
                                
                            }//end o regions loop
                            
                            
                        } //if (i!=j)
                        
                    }// j positions
                    
                    
                    
                } //end if (  SA[i] % L==0)
                
                free(mismatches);
                free(mismatches_j_pos);
                free(r_pos);
                free(l_pos);
                //if ((table.size()> 0))
                
                //   dic->at(tid)= table;
                
            }// loop for i
            i=j;
            
      
            
            
            
     
            
        }//end if
        
        else
            i++;
        
        
    }//end loop
    
    
    
    
    //
    
    
    
    FILE * out_fd;
    if ( ! ( out_fd = fopen ( sw . output_filename, "a") ) )
    {
        fprintf ( stderr, " Error: Cannot open file %s!\n", sw . output_filename );
        return ( 1 );
    }
    fprintf ( out_fd, ">%s\n", ( char * ) seq_id );
    for ( INT i = 0; i < n; i++ )
        fprintf ( out_fd, "%4.00f ", C[i] );
    fprintf( out_fd, "\n" );
    
    if ( fclose ( out_fd ) )
    {
        fprintf( stderr, " Error: file close error!\n");
        return ( 1 );
    }
    
    
    free ( SA );
    free ( LCP );
    free ( invSA );
    
    //free(mismatches);
    //free(mismatches_j_pos);
    //free(r_pos);
    //free(l_pos);
    
    return ( 1 );
}

    




unsigned int Counting_Blocks_for_K_Mismatch ( INT p, INT L,INT * mismatches, INT m, INT k) {

     INT  B=0,c=0,x=0,y=0,d=0;
     INT start_block_of_mismatch=0,end_block_of_mismatch=0;
    

     if (p%L !=0){  //if it is not staring poistion of block  move to the second  block
     x=p%L;
     c=p+L-x;
     }
     else
     c=p;
    
    
     if ((p+m)%L==0)  //if is the end poistion is end of th block otherwise back to the end pos of the block
     d=p+m;
    
     else{
     y=(p+m)%L;
     d=p+m-y;
     }

 if(d>c){
     B=(d-c)/L;
     for (INT i=0;i<k;i++){
         if (i==0){
             if(c<=mismatches[i]&& mismatches[i]<d){
         B=B-1;
                }
         }
         else{
        start_block_of_mismatch= mismatches[i-1]-( mismatches[i-1]%L);
        end_block_of_mismatch= (start_block_of_mismatch+L )-1;
           if(c<=mismatches[i] && mismatches[i]<d ){
               
               if ( mismatches[i] > end_block_of_mismatch)
            B=B-1;
    }
     }
     }
     } // end if d > c
     else {
     B=(c-d)/L;
        // cout << B <<endl;
          for (INT i=0;i<k;i++){
          
              if (i==0){
                  if(d<=mismatches[i] && mismatches[i]<c){
                    B=B-1; ;
                  }
              }
                  else{
                
                      start_block_of_mismatch= mismatches[i-1]-( mismatches[i-1]%L);
                      
                      end_block_of_mismatch= (start_block_of_mismatch+L)-1;
     
     if(d<=mismatches[i] && mismatches[i] < c  )
     if (  mismatches[i] > end_block_of_mismatch)
     
     B=B-1;
              }
              }
     }
    return B;
    
    
}


unsigned int Counting_Blocks_for_Two_Mismatch ( INT p, INT L,INT first_mismatch,INT second_mismatch, INT m) {
   // cout <<"after call"<<endl;
    //cout << first_mismatch << " "<<second_mismatch <<endl;
    
    
    
    INT  B=0,c=0,x=0,y=0,d=0;
    INT start_block_of_mismatch=0,end_block_of_mismatch=0;

    start_block_of_mismatch= first_mismatch-( first_mismatch%L);
    end_block_of_mismatch= (start_block_of_mismatch+L )-1;
    if (p%L !=0){  //if it is not staring poistion of block  move to the second  block
        x=p%L;
        c=p+L-x;
    }
    else
        c=p;
    
    
    if ((p+m)%L==0)  //if is the end poistion is end of th block otherwise back to the end pos of the block
        d=p+m;
    
    else{
        y=(p+m)%L;
        d=p+m-y;
    }
    /* to find the block that the first mistach belog to it*/
    
    
    
    if(d>c){
        
        B=(d-c)/L;
        
        /* for the first mismatch */
        if(c<=first_mismatch && first_mismatch<d){
            B=B-1;
            
            
        }
        
        
        
        if(c<=second_mismatch && second_mismatch<d ){
            
            if (  second_mismatch > end_block_of_mismatch)
                
                B=B-1;
            
            
        }
        
    } // end if d > c
    
    
    
    
    
    else {
        
        
        B=(c-d)/L;
        
        if(d<=first_mismatch && first_mismatch<c)
            B=B-1;
        
        if(d<=second_mismatch && second_mismatch < c  )
            if (  second_mismatch > end_block_of_mismatch)
                
                B=B-1;
        
    }
  //  cout<<"Blocks " << B<<endl;
    return B;
    
    
}

unsigned int Counting_Blocks ( INT p, INT L,INT mismatch, INT m) {
    INT  B=0,c=0,x=0,y=0,d=0;

    if (p%L !=0){  //if it is not staring poistion of block  move to the second  block
        x=p%L;
        c=p+L-x;
    }
    else
        c=p;
    
    
    if ((p+m)%L==0)
        d=p+m;
    
    else{
        y=(p+m)%L;
        d=p+m-y;
    }
    
    
    if(d>c){
        B=(d-c)/L;
        
        if(c<=mismatch && mismatch<d){
            B=B-1;
        }
        
        
    }
    
    else
    {
        B=(c-d)/L;         B=(d-c)/L;

        if(d<=mismatch && mismatch<c){
            B=B-1;
        }
    }
   // cout << "one block "<< B<<endl;
    return B;
    
    
}




