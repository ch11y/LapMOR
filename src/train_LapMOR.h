#ifndef __TRAIN_LAPMOR_H_
#define __TRAIN_LAPMOR_H_ 

#include <cstdio>
#include <cstring>
#include <map>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <set>
#include <bitset>
using namespace std ;  

#define SZ(x) (int)x.size()
const double eps =  1e-5 ;  
typedef vector< double >  VD ;  

bool converge(const  VD &  x , const VD&  y ){
    double res = 0 ; 
    for(int i = 0; i < SZ(x) ; ++i)
        res = res + abs( x[ i ]  - y[ i ] ) ;
    return res < 10*eps ; 
}


void format_input( const vector < VD > & input,  vector < VD > &g_input ){ 

    int n = SZ(input) , p = SZ(input[0]), q = SZ(g_input)/n;  
    VD tmp( p*q ) ;  
    fill(tmp.begin() , tmp.end() , 0 ) ;  

    for(int i = 0 ; i < n ; ++i){         
        for(int j=0 ; j < q;  ++j){  
            for(int k= 0 ;  k < p ; ++k){
                tmp[ j * p + k ] = input[i][k];   
            }
            
            g_input[ i * q + j ] = tmp ;  
            
            for(int k= 0 ; k <p ; ++k){ 
                tmp[ j * p + k ]  = 0 ; 
            }
        }
    }
}

void format_output( const vector < VD > &output, VD & g_output ){ 
    
    int n = SZ(output), q = SZ(output[0]) ; 
    for(int i = 0; i < n ; ++i){
        for(int j=0 ; j < q ; ++j){ 
            g_output[i*q + j ] = output[i][j]; 
        }
    }
}

void format_B( const VD & g_B,  vector< VD > &B ){
    int p = SZ(B) , q = SZ(B[0]) ;  
    for(int i=0 ; i < q ; ++i){
        for(int j= 0; j < p ; ++j){
            B[j][i] = g_B[ i * p + j ]  ;
        }
    }
}

inline double soft_thresholding( double  val, double lambda){
   
    if( fabs( val ) - lambda < 0 ) return 0.; 
    if( val > 0 ) return val - lambda; 
    else return val+ lambda ;    
}

inline double rand_double( ){ 
    return 2.*rand()/RAND_MAX -1 ; 
}
 
static int times = 0 ;
void solve_glasso_sub(const vector < VD > &g_input, const VD &g_output,
        VD & g_B, const VD& A , const double & lambda1,const double & lambda2, const vector< vector<int> >& group){
    if( !times )
    {
        cout <<"init"<<endl;
    for(int i=0; i < SZ(g_B) ; ++i)
       g_B[i] = 0.0;
    }
    ++times;
    int g_n = SZ(g_input), g_p = SZ(g_input[0]), g_num = SZ(group)  ; 
    VD preB  ;    
    set < int > group_set ;  
    VD div = g_input[0] ; 
    for(int i=0; i < SZ(g_input[0]); ++i){
        div[i]=0 ;
        for(int j=0;j <SZ(g_input); ++j){ 
            div[i] = div[i] + g_input[j][i]*g_input[j][i]; 
        }
    }
    int cnt = 0 ; 
    VD tmp_output = g_output ;   
    for(int i=0; i < g_n ; ++i){
        tmp_output[i]=0 ;  
        for(int j= 0 ; j < SZ(group) ; ++j){
            for(int k=0; k < SZ(group[j]);++k){
                tmp_output[i] = tmp_output[i] + g_input[i][ group[j][k]] * g_B[ group[j][k]] ; 
            }
        }
    }
    for(int i=0; i < SZ(group); ++i){
        group_set.insert( i ) ; 
    }

    do{  
        preB = g_B ; 
        vector< int > toRemove ;
        for(set< int > :: iterator it = group_set.begin();  it != group_set.end() ;++it){
            int i = *it ; 
            double norm_tmp = 0 ; 
            for(int j=0; j < SZ(group[i]); ++j){ 
                int id  =  group[i][j];   
                double val = 0 ;
                if(fabs(div[id]) < 1e-5) continue ; 
                for(int k=0; k < g_n ; ++k){ 
                    val = val +  g_input[k][id] *(g_output[k] - (tmp_output[k]- g_input[k][id]*g_B[id]) ) ;  
                } 
                double preB = g_B[id]; 
                double newB = soft_thresholding( val   , A[i]*lambda1*lambda2)/( div[id] + A[i]*(1-lambda2)*lambda1)  ;
                g_B[id] = newB  ;  
                norm_tmp = norm_tmp+ abs(newB); 
                for(int k=0; k < g_n; ++k){ 
                    tmp_output[k] = tmp_output[k] + g_input[k][id] *(newB-preB) ;
                }
            }
            if(fabs(norm_tmp) < eps) {
                toRemove.push_back(i);
            }

        }
        for(int i=0; i < SZ(toRemove); ++i){ 
            group_set.erase( toRemove[i]) ;  
        }
    }while(!converge(preB, g_B)) ;  
}

bool decrease_Laplacian( const VD& normB_1 ,const VD& normB_2, const VD& A, const VD& nextA , const double & lambda1, const double & lambda2, int K){
    double res1 = 0 , res2 = 0 ; 
    for(int i=0; i < SZ(A); ++i){
        assert(A[i]>0 &&A[i]<=1) ;
        assert(nextA[i]>0 && nextA[i]<=1) ; 
        res1 = res1 + A[i]*normB_1[i]*lambda1*lambda2 + 0.5*A[i] * normB_2[i] * lambda1*(1-lambda2)-1.5* log( A[i])*K ; 
        res2 = res2 + nextA[i]*normB_1[i]*lambda1*lambda2+0.5*nextA[i]*normB_2[i]* lambda1*(1-lambda2) - 1.5*log(nextA[i])*K ;  
    }
    return !(res1 < res2 || fabs(res1-res2)< eps) ;   
}

void Generate( VD &A){
    vector< int > x (SZ(A)) ;
    int tot = 0 ; 
    for(int i=0; i < SZ(A); ++i){
        x[i]= rand() % 10000 + 1 ; 
        tot += x[i ] ;  
    }
    for(int i=0; i < SZ(A ) ; ++i){
        A[i] = 1. * x[i]/ tot ; 
    }

}
void update_A_Laplacian( const VD& g_B, VD& A , const double & lambda1, const double& lambda2, int K ){ 
     
     Generate(A) ;
     cout <<"update_A_Laplacian"<<endl; 
     for(int i=0; i < SZ(A); ++i){
         assert(A[i]>0&&A[i]<=1) ; 
     }
     
     double delta[100] ;
     delta[0]=0.0001; 
     for(int i=1; i < 100; ++i){ 
         delta[i]= delta[i-1]*0.90;  
     } 
    
     VD nextA = A ;  
     VD normB_1 = A, normB_2 = A ;  

     for(int i=0; i < SZ(A); ++i){
         normB_1[i]=normB_2[i]=0.; 
         for(int j=0; j < SZ(g_B)/SZ(A); ++j){
            normB_1[i] = normB_1[i] + fabs( g_B[j * SZ(A) + i ] ) ;  
            normB_2[i] = normB_2[i] + g_B[ j * SZ(A)+ i ] * g_B[j*SZ(A)+i] ; 
         }         
     } 
     
     for(int step = 0 ; step < 100 ; ++step ){  
         nextA = A ;
         bool legal = 1 ;
         do{
           A= nextA ;
           for(int i=0; i < SZ(A); ++i){
               nextA[i] = nextA[i] - delta[step] * (  lambda1*lambda2*normB_1[i] + lambda1*(1-lambda2)*normB_2[i] - 1.5*K /(nextA[i]) ) ; 
           }
           double total = 0 ; 
           for(int i=0 ; i <SZ(A) ;++i){
               total = total + nextA[i] ; 
           }
           total -=1 ; 
           for(int i=0; i < SZ(A); ++i){ 
               nextA[i] = nextA[i] - total / SZ(A) ; 
           }
           
           for(int i=0; i <SZ(A); ++i){ 
               if(nextA[i] <=0.0000001 || nextA[i]>=1 ){
                   legal=0 ;
                   break ; 
               }
           }
         } while(legal&& decrease_Laplacian( normB_1,normB_2 , A , nextA, lambda1, lambda2, K) ) ;  
     }
}
void solve_glasso(const vector < VD > &g_input , const VD &g_output , VD& g_B ,  VD & A , const double & lambda1 , const double& lambda2 ){ 
    vector< vector<int> > group(SZ(A)) ; 
    int p = SZ(A) ,q = SZ(g_B)/SZ(A); 
      
    for(int i=0; i < SZ(A) ; ++i){
    
        group[i].clear() ; 
        for(int j=0 ; j < q ; ++j){ 
            group[i].push_back( j * p + i) ; 
        }
    }

    for(int i=0; i < SZ(A); ++i){
        A[i]= 1./ p  ;
    }
    int cnt = 0 ; 
    VD preA  ;                        
    do{
        preA = A ; 
        solve_glasso_sub( g_input , g_output , g_B ,A, lambda1, lambda2, group) ; 
        update_A_Laplacian( g_B , A, lambda1 , lambda2,q );
        ++cnt;
    }while( cnt < 500 && !converge(preA , A ) ) ;
}
void train_LapMOR(const vector< VD > & input,  const vector< VD > & output ,  
        vector< VD > & B , VD& A ,  double & lambda1, double&  lambda2){
    int n = SZ(input), p = SZ(input[0]), q = SZ(output[0]); 
    A = vector< double > (p, 1.0/p);
    B = vector< vector<double> > (p, vector<double>(q, 0.));
    vector < VD > g_input(n*q); 
    VD g_output(n*q);  
    VD g_B( p * q ) ; 
    format_input( input  , g_input ) ; 
    format_output( output , g_output ) ; 
    solve_glasso( g_input , g_output ,g_B, A, lambda1, lambda2) ;  
    format_B( g_B , B) ;   
}
#endif
