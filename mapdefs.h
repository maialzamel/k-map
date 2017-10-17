/**
 K-MAP: In the k-mappability problem, we are given a string x of length n and integers m and k, and we are asked to count, for each length-m factor y of x, the number of other factors of length m of x that are at Hamming distance at most k from y.Copyright (C) 2017 Mai Alzamel.

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

#define ALLOC_SIZE              1048576
#define DEL                     '$'
#define DEL_STR                 "$"

#define DNA                     "ACGTN"                         //DNA alphabet
#define PROT                    "ARNDCQEGHILKMFPSTWYVOUBZJX*"          //Proteins alphabet
#define IUPAC                   "ACGTUWSMKRYBDHVN"          	//IUPAC alphabet
#define max(a,b) ((a) > (b)) ? (a) : (b)
#define min(a,b) ((a) < (b)) ? (a) : (b)


#ifdef _USE_64
typedef int64_t INT;
#endif

#ifdef _USE_32
typedef int32_t INT;
#endif
struct TSwitch
 {
   char               * alphabet;
   char               * input_filename;
   char               * output_filename;
   unsigned int         k,m;
   unsigned int         total_length;
 };
struct Range {
    int i;
    int j;
};


double gettime( void );
int decode_switches ( int argc, char * argv [], struct TSwitch * sw );
void usage ( void );

unsigned int LCParray ( unsigned char *text, INT n, INT * SA, INT * ISA, INT * LCP );
unsigned int LCE(INT * lcp, INT* ISA, INT n,int L, int R);
unsigned int Hamming_Dist ( unsigned char * sub_1, unsigned char * sub_2, INT m );





unsigned int compute_map ( unsigned char * seq, unsigned char * seq_id, struct TSwitch sw, INT m ,   double  * C  );
unsigned int compute_map_simple ( unsigned char * seq, unsigned char * seq_id, struct TSwitch sw, INT m ,   double  * C  );

unsigned int compute_ZeroMapp (INT n, INT m ,INT  * SA ,  INT  * LCP, 	INT * invSA, double * C );


unsigned int compute_At_Most_One_map_simple ( unsigned char * seq,unsigned char * seq_id, struct TSwitch sw, INT m ,   double  * C );

unsigned int compute_At_Most_Two_map_simple ( unsigned char * seq,unsigned char * seq_id, struct TSwitch sw, INT m ,   double  * C );

unsigned int compute_naive_map ( unsigned char * seq, unsigned char * seq_id, struct TSwitch sw,   double * C, INT m, unsigned int k );

unsigned int compute_naive_at_most_k_map ( unsigned char * seq, unsigned char * seq_id, struct TSwitch sw, double * C, INT m , unsigned int k );

unsigned int Counting_Blocks ( INT p, INT L,INT l_1, INT m);
unsigned int Counting_Blocks_for_Two_Mismatch ( INT p, INT L,INT mismatch,INT mismatch_1, INT m);




