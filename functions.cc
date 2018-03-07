
/**
 K-MAP:   In the k-mappability problem, we are given a string x of length n and integers m and k, and we are asked to count, for each length-m factor y of x, the number of other factors of length m of x that are at Hamming distance at most k from y.
 Copyright (C) 2017 Mai Alzamel .
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
unsigned int compute_At_Most_K_map_simple_mutlthreads ( unsigned char * seq,unsigned char * seq_id, struct TSwitch sw, INT m ,   double  * C ,unsigned int k)
{



      vector<pair<int, int>> chunk;

    INT L= floor ((m/(k+2.0))*1.0);
    int number_of_mismatch=0;
    INT n = strlen ( ( char * ) seq );

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

/**********************************************************************************************************************************************/
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
     if (alpha!=beta)
              chunk.push_back(make_pair(alpha,beta));
     
             

                
             i=j;
        }
            else
                i++;
    }




/****************************************************************************************************************************************/
    


            //************************Processing each set between Alpha-1 and Beta one interval ***********************//
  omp_set_num_threads(4);
  

int interval;

//#pragma omp parallel for schedule(static, chunk.size()/4)




for ( interval =0; interval<chunk.size();interval++){

INT alpha= chunk.at(interval ).first;
 INT beta=chunk.at(interval ).second;
cout << alpha <<" " <<beta<<endl;
#pragma omp parallel for schedule( auto)

           for(INT i=alpha-1;i<=beta;i++)
            {

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
                {                                                               //check if it is a starting position of a block.
                    
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
#pragma omp atomic
                                                C[q]+=1.0/(B_i+B_j);
#pragma omp atomic
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
#pragma omp atomic
                                                C[p]+=1.0/(B_i+B_j);
 #pragma omp atomic
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
           

        
        
   } 
    
    
}
 

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

