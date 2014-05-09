#ifndef __MODEL_LAPMOR_H_
#define __MODEL_LAPMOR_H_

#include <cstdio>
#include <cstring>
#include <map>
#include <iostream>
#include <algorithm>
#include <string>
#include <vector>

#include "train_LapMOR.h"

using namespace std ;   

#define SZ(x) (int)x.size()

typedef  vector< double >  VD ;

double predict_LapMOR(const vector < VD >  & input ,const vector < VD >  & output,
        vector< VD >&  B , VD & A , double lambda1, double lambda2, bool o=false){  
 
    cout <<  "predict_LapMOR" << endl ; 
    double res = 0 ;
    int n  = SZ(input) , p = SZ(input[0]), q = SZ( output[0] )  ;  
    vector< VD > tmp = output ; 
    for(int now =  0 ; now < q ; ++now){

        for(int i=0; i < n ; ++i){
            
            double y_hat = 0 ;
            for(int j= 0 ; j < p ; ++j){ 
                y_hat += input[i][j] * B[j][now]  ; 
            }
            res = res + pow( y_hat  - output[i][now] , 2.) ; 
            tmp[i][now] = y_hat;
        }

    }
    if( o ) { 
        for(int i=0; i < n ; ++i){
            for(int j=0 ; j < q ; ++j){ 
                cout << tmp[i][j] << " " ; 
            }
            cout << endl ; 
        }
    }
    

    return res / n / q  ;
}

double mean( const vector< double > & vc ){ 
    double res = 0  ; 
    for(int i = 0 ; i < SZ(vc) ; ++i){
        res = res + vc [ i  ] ; 
    }
  
    res /= SZ(vc); 
    return res;
}

double variance(const vector< double > &vc ){ 
    double m = mean( vc ) ; 
    double res = 0 ;
    for(int i = 0; i < SZ(vc) ; ++i){ 
        res += pow( vc[ i ] - m , 2 ) ; 
    }
    res/= SZ(vc) ; 
    return res; 
}

void model_LapMOR( const vector< VD > &  input , const vector< VD > & output, vector< VD > & B , VD & A, double & lambda1 , double& lambda2, int div = 5){

    cout <<"model_mtrfs" << endl ; 
   
    pair< double , double >  best_cv  = make_pair(100000000000.,100000000000.) ; 
    pair< double , double >  best_lambda = make_pair( -1 , -1 ) ; 
    int n = SZ(input) ,  p = SZ(input[0]), q = SZ(output[0]) ; 
    int step = n / div ;
   
  
    cout << "model_selection" << endl ;  
    vector< VD > train_input( n - step ), test_input( step ), train_output( n - step ), test_output( step )   ; 
    
    vector< VD > pre_B = B ;  
   // for( lambda1 = 1 ;  lambda1 <= 1 ; lambda1 += 1 ){
      
        bool t = 0;
        B = pre_B ;
       // times= 0 ; 
     //   for( lambda2 = 1 ; lambda2 <= 1 ; lambda2 += 1){
            //lambda1 *= 10 ;lambda2*=10; 
            vector < double > predict_error ; 
            vector < double > train_error ;  
            int cnt = 0 ;
            
            for(int start = 0 ; start < n ; start += step ){
               // times = 0  ; 
              //  if( cnt ++ > 0 ) break ; 
                cout << "CV #: "<< ( ++cnt ) <<endl; 
                for(int i= 0 ; i < n - step ; ++i){ 
                    train_input[i] = input[ (start + i) % n  ]  ;  
                    train_output[i] = output[ (start + i ) % n ] ; 
                }
                for(int i= n - step ; i < n ; ++i){
                    test_input[ i - ( n- step ) ]  = input[ (start+i) % n ] ;   
                    test_output[ i -( n- step ) ]  = output[(start+i) %n  ]; 
                }
                    
                 
                train_LapMOR( train_input, train_output , B , A ,  lambda1 , lambda2 ) ;
                for(int i=0; i < SZ(A) ;++i){
                    cout <<A[i]<<" " ;
                }
                cout << endl; 
                train_error.push_back(  predict_LapMOR( train_input, train_output , B , A , lambda1, lambda2)) ;  
                predict_error.push_back( predict_LapMOR( test_input, test_output , B , A , lambda1, lambda2 )) ;       
                cout <<"lambda1: "<< lambda1 << "  lambda2: " << lambda2<< "  Train Error: "<< train_error[SZ(train_error)-1] << endl ; 
                cout <<"lambda1: "<< lambda1 << "  lambda2: " << lambda2<< "  Test Error:  "<< predict_error[SZ(predict_error)-1]  << endl ;
            }
                     
            cout <<"predict error: " ;
            for(int i=0 ; i < SZ(predict_error); ++i){ 
                cout << predict_error[i] << " " ; 
            }
            cout <<endl; 
            double m = mean ( predict_error  ) , v = variance( predict_error ) ; 
            if( make_pair( m , v ) < best_cv ){

                best_cv = make_pair( m , v ) ;  
                best_lambda = make_pair( lambda1 , lambda2) ;

            }
 
            cout <<"-------------------------------------------------" << endl ; 
            cout <<"Average: "<< lambda1 << " "<< lambda2 << " " << 
                m << " "<< v << endl  ;
      //  }
   // }
    times = 0 ;  
    lambda1 = best_lambda.first , lambda2 = best_lambda.second ; 
    train_LapMOR( input, output,  B , A ,  lambda1, lambda2) ; 
}

#endif


