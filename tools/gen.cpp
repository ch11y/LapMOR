#include <cstdio>
#include <string>
#include <cstring>
#include <sstream>
#include <vector>
#include <algorithm> 
#include <cstdlib> 
#include <cmath>

using namespace std ;   

#define SZ(x) (int)x.size()

const double eps = 1e-7 ; 
char ss[1000000+10]; 


int main (int argc , char *argv [] ){
    vector< pair< string , string > > vc ;  
    for(int i=1 ; i < argc ; i+=2 )
        vc.push_back( make_pair( string( argv[i]) , string( argv[i+1]) ) ) ; 
     
    sort(vc.begin(), vc.end()) ;  
    
    FILE *fp_in = fopen(vc[0].second.c_str(), "r") ;
    FILE *fp_out = fopen(vc[1].second.c_str() , "r") ; 
    FILE *fp_pa = fopen( vc[2].second.c_str(), "r") ; 

   /* FILE *fp = fopen("../input_s.txt","w") , fp4 = fopen("../output_s.txt", "w")  ; 
    FILE *fp5 = fopen("../input_t.txt","w") , fp6 = fopen("../output_t.txt", "w")  ;    */  

    vector< vector< double > > input , output ;  
    vector< vector< double > > parameter ;  

    while( fgets( ss,  1000000 , fp_in) != NULL ){ 
         string s = string(ss) ;  
         stringstream sin(s) ;  
         double x;
         vector< double > tmp ; 
         while( sin >> x )
             tmp.push_back(x) ; 
         input.push_back( tmp ) ;  
    }
    fclose(fp_in) ;  
    while( fgets( ss ,  1000000, fp_out ) != NULL){ 
        string s = string(ss) ; 
        stringstream sin(s) ; 
        double x ;  
        vector< double > tmp ; 
        while ( sin >> x ) 
            tmp.push_back( x ) ;
        output.push_back(tmp) ; 
    }
    fclose( fp_out ) ; 
    while( fgets(ss , 1000000,  fp_pa ) != NULL ){ 
        string s  = string( ss ) ; 
        stringstream sin(s) ; 
        double x ; 
        vector< double > tmp ;  
        while( sin >> x )
            tmp.push_back(x) ;   
        parameter.push_back( tmp ) ;   
    }
    fclose( fp_pa ) ; 
    for(int p = 0 ;  p < SZ( parameter[0] ) ; ++p){
        stringstream sin ; 
        char buf[100] ; 
        sprintf( buf , "%d", p ); 
        string s_in = vc[0].second +"_";
        string s_out = vc[1].second +"_" ; 
        s_in += string(buf) ; 
        s_out += string(buf) ; 
        fp_in = fopen( s_in.c_str() , "w")  ; 
        fp_out = fopen( s_out.c_str() , "w")  ;   

        for(int i = 0 ; i < SZ(input) ; ++i){
            for(int j=0 ; j < SZ(input[0]) ; ++j){ 
             //   if( fabs( parameter[j][p] ) < eps ) continue ; 
                fprintf( fp_in , "%lf ", input[i][j]) ;  
            }
            fprintf(fp_in , "\n") ;  
            fprintf(fp_out, "%lf\n", output[i][p]) ;
        }
        fclose(fp_in) ; fclose(fp_out) ; 
    }
    return 0; 
}


             
             
