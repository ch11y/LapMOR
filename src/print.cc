#include <cstdio>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>
#include <sstream>
#include <cassert>

#include "model_amtrfs.h" 

using namespace std ; 

#define SZ(x) (int)x.size()

//"-i input.txt -o outptut.txt -p ans.txt" 



char ss[5000000+ 10] ; 
void init( vector< VD > & data){
    int n = SZ(data) , p = SZ(data[0]) ;
    for(int i=0 ; i< p ; ++i){ 
        double av = 0 , ma= -1 , mi=1000000000.;
        for(int j=0 ; j < n ; ++j){ 
            av += data[j][i];  
            ma = max( data[j][i] , ma) ; 
            mi = min( data[j][i] , mi) ; 
        }
        av /=n ; 
        double var = 0 ; 
        for(int j=0 ; j < n ; ++j){ 
            data[j][i] -= av ;
            var += data[j][i] * data[j][i];  
            //if( fabs( ma - mi ) <  1e-9 ) continue ;  
           // data[j][i] /= (ma -mi) ; 
        }
        var = sqrt( var/n ) ;  
        for(int j=0 ; j < n ; ++j){ 
            //if(fabs(ma - mi) < 1e-7 ) continue ;
            if( fabs(var) < 1e-7 ) continue ; 
            data[j][i] /= var ;  
        }
       
    }
}

int main (int argc , char **argv ){ 
    
    clock_t start = clock() ; 
    vector< pair< string , string > > parameter  ; 
    for(int i=1; i < argc ; i+=2)
        parameter.push_back( make_pair( argv[i] , argv[i+1] ) ) ; 
    sort( parameter.begin() , parameter.end())   ;  
    assert( 5 == SZ(parameter) ) ;      
      
    FILE *fp = fopen( parameter[0].second.c_str() , "r" ) ; 
    assert( NULL != fp ) ;     
    vector< VD > input , output ; 
      
    while ( fgets( ss , 5000000, fp ) != NULL ){ 
        //++lines;  
        string s = string( ss ) ; 
        stringstream sin(s) ;  
        VD tmp ; 
        double x;
        int cnt = 0  ; 
        while( sin >> x ){
    //        if( ++cnt > 300 ) continue ;
            tmp.push_back(x); 
        }
            
        input.push_back( tmp ) ; 
    }
  
    fclose(fp) ;  
    fp = fopen( parameter[1].second.c_str() , "r")  ;  
    assert( NULL != fp ) ; 
    while( fgets( ss , 5000000, fp ) != NULL ){ 
        string s = string( ss )  ; 
        stringstream sin(s) ;  
        VD tmp ;  
        double x ;
        while( sin >> x  ){ 
           tmp.push_back(x) ; 
        }
        output.push_back(tmp) ; 
    }
    
    fclose(fp) ;  
    init( input ) ; 
    init( output ) ; 
    
  
    for(int i = 0 ; i < SZ(input) ; ++i){ 
       int sw=  rand()%(i+1) ;
       
       swap( input[i], input[ sw ] ) ; 
       swap( output[i], output[sw] ) ;
    } 
    for(int i=0 ; i < SZ(input) ; ++i){ 
        for(int j=0 ; j < SZ(input[i]) ;++j){ 
            cout << input[i][j]<< " " ; 
        }
        cout << endl ;
    }

    for(int i=0 ; i < SZ(output); ++i){
        for(int j=0 ;j <SZ(output[i]);++j){
            cout << output[i][j] << " " ; 
        }
        cout << endl ; 
    }
    cout <<"Read Finished" << endl ; 
    for(int i=0 ; i < SZ(output) ; ++i){
        for(int j=0 ; j < SZ(output[i]) ; ++j){ 
           // cout << output[i][j] << " " ; 
        }
        //cout << endl; 
    }
    //exit(0) ; 
    for(int i=0 ; i < SZ(input); ++i){
        input[i].push_back(0.1) ; 
    }
    double lambda1 , lambda2 ; 
    vector < VD > A(SZ(input[0]))  ; 
    VD B(SZ(A)) ;  
    VD tmp( SZ(output[0])) ; 
    for(int i=0 ; i < SZ(A) ;++i){
        A[ i ] = tmp ;  
    }
    lambda1 = atof( parameter[3].second.c_str() ) ; 
    lambda2 = atof( parameter[4].second.c_str() ) ; 
    model_amtrfs( input , output, A , B , lambda1 , lambda2) ; 
    
    fp = fopen( parameter[2].second.c_str() , "w") ;  
    assert( NULL != fp ) ;  
    for(int i = 0 ; i < SZ(A) ; ++i){
        for(int j = 0 ; j < SZ(A[0]) ; ++j){
            fprintf( fp , "%lf " , A[i][j] ) ; 
        }
        fprintf( fp , "\n" ) ; 
    }
    fclose(fp) ; 
    clock_t end = clock() ; 
    cout << ( end - start ) / CLOCKS_PER_SEC << endl ; 
    return 0 ;
}

    


#ifndef __MODEL_AMTRFS_H_
#define __MODEL_AMTRFS_H_

#include <cstdio>
#include <cstring>
#include <map>
#include <iostream>
#include <algorithm>
#include <string>
#include <vector>

#include "train_amtrfs.h"

using namespace std ;   

#define SZ(x) (int)x.size()

typedef  vector< double >  VD ;

double predict_amtrfs(const vector < VD >  & input ,const vector < VD >  & output,
        vector< VD >&  A , VD & B , double lambda1, double lambda2){  
 
    cout <<  "predict_amtrfs" << endl ; 
    double res = 0 ;
    int n  = SZ(input) , p = SZ(input[0]), q = SZ( output[0] )  ;  

    for(int now =  0 ; now < q ; ++now){

        for(int i=0; i < n ; ++i){
            
            double y_hat = 0 ;
            for(int j= 0 ; j < p ; ++j){ 
                y_hat += input[i][j] * A[j][now]  ; 
            }
            res = res + pow( y_hat  - output[i][now] , 2.) ; 
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

void model_amtrfs( const vector< VD > &  input , const vector< VD > & output, vector< VD > & A , VD & B , double & lambda1 , double& lambda2, int div = 5){

    cout <<"model_mtrfs" << endl ; 
   
    pair< double , double >  best_cv  = make_pair(100000000000.,100000000000.) ; 
    pair< double , double >  best_lambda = make_pair( -1 , -1 ) ; 
    int n = SZ(input) ,  p = SZ(input[0]), q = SZ(output[0]) ; 
    int step = n / div ;
   
  
    cout << "model_selection" << endl ;  
    vector< VD > train_input( n - step ), test_input( step ), train_ouput( n - step ), test_output( step )   ; 
    
    vector< VD > pre_A = A ;  
   // for( lambda1 = 1 ;  lambda1 <= 1 ; lambda1 += 1 ){
      
        bool t = 0;
        A = pre_A ;
       // times= 0 ; 
     //   for( lambda2 = 1 ; lambda2 <= 1 ; lambda2 += 1){
            //lambda1 *= 10 ;lambda2*=10; 
            vector < double > predict_error ; 
            vector < double > train_error ;  
            int cnt = 0 ; 
            for(int start = 0 ; start < n ; start += step ){
             
                cout << "CV #: "<< ( ++cnt ) <<endl; 
                for(int i= 0 ; i < n - step ; ++i){ 
                    train_input[i] = input[ (start + i) % n  ]  ;  
                    train_ouput[i] = output[ (start + i ) % n ] ; 
                }
                for(int i= n - step ; i < n ; ++i){
                    test_input[ i - ( n- step ) ]  = input[ (start+i) % n ] ;   
                    test_output[ i -( n- step ) ]  = output[(start+i) %n  ]; 
                }
                    
                 
                train_amtrfs( train_input, train_ouput , A , B , lambda1 , lambda2 ) ; 
                train_error.push_back(  predict_amtrfs( train_input, train_ouput , A , B , lambda1, lambda2)) ;  
                predict_error.push_back( predict_amtrfs( test_input, test_output , A , B ,lambda1, lambda2 )) ;       
                cout <<"lambda1: "<< lambda1 << "  lambda2: " << lambda2<< "  Train Error: "<< train_error[SZ(train_error)-1] << endl ; 
                cout <<"lambda1: "<< lambda1 << "  lambda2: " << lambda2<< "  Test Error:  "<< predict_error[SZ(predict_error)-1]  << endl ;
            }
            if(!t) pre_A = A ; 
            t = 1 ; 
            
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
            cout << lambda1 << " "<< lambda2 << " " << 
                m << " "<< v << endl  ;
      //  }
   // }
    times = 0 ;  
    lambda1 = best_lambda.first , lambda2 = best_lambda.second ; 
    train_amtrfs( input, output,  A , B , lambda1, lambda2) ; 
}

#endif


#ifndef __TRAIN_AMTRFS_H_
#define __TRAIN_AMTRFS_H_ 

#include <cstdio>
#include <cstring>
#include <map>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <set>
using namespace std ;  

#define SZ(x) (int)x.size()
const double eps =  1e-6 ;  
typedef vector< double >  VD ;  

//yep
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

void format_A( const VD & g_A,  vector< VD > &A ){
    int p = SZ(A) , q = SZ(A[0]) ;  
    for(int i=0 ; i < q ; ++i){
        for(int j= 0; j < p ; ++j){
            A[j][i] = g_A[ i * p + j ]  ;
        }
    }
}

inline double rand_double( ){ 
    return 2.*rand()/RAND_MAX -1 ; 
}


inline int cal_t_hat(double x, double y){
    if( fabs( x / y ) <= 1  ) return x/y; 
    else return  ( x / y )  > 0 ? 1:-1  ; 
}

double ternary_search( const vector < VD > g_input, const VD& sub_output , VD& g_A,const VD& B, 
        const double& lambda1, const double& lambda2, const vector< int > group, int id, VD & theta){ 
  
    double l = -3 , r = 3 ;
    int cnt = 0 ;
    double sub_output_sq = 0 , sub_output_z = 0,  z_sq=0;   
    for(int i = 0 ; i < SZ(g_input) ; ++i){
       // sub_output_sq += sub_output[i] * sub_output[i]; 
        sub_output_z += 2*sub_output[i]*g_input[i][ group[id]] ; 
        z_sq += g_input[i][group[id]] * g_input[i][group[id]] ;  
    }
    while( fabs( r -l ) > eps/10 && cnt < 100  ){ 
        double m1=  l + (r-l) /3 ,  m2 = l + 2 * ( r- l ) / 3;  
        double r1 = 0, r2 = 0 ; 
      //  for(int i=0 ; i < SZ(g_input) ; ++i){ 
        //    r1 += pow( sub_output[i] - g_input[i][group[id]] * m1 ,2) ; 
          //  r2 += pow( sub_output[i] - g_input[i][group[id]] * m2, 2 ) ; 
       // }
        r1  =   - sub_output_z * m1 + z_sq * m1*m1 , r2 = -sub_output_z*m2 + z_sq*m2*m2;
        r1/=2 ,r2/=2 ;
        ++cnt; 
        double tmp1 =0 , tmp2 = 0 ;
        /////////////////////////////////////////
        for(int i=0; i < SZ(group) ; ++i){ 
            if( id == i ){ 
                tmp1 += m1 * m1 ; 
                tmp2 += m2 * m2 ; 
                r1 += abs( m1 ) * lambda2 ;
                r2 += abs(m2) * lambda2 ; 
            }else{
                tmp1 += theta[i] * theta[i ] ;
                tmp2 += theta[i] * theta[i]  ;               
            }

        }
        r1 += lambda1 *  sqrt( tmp1) ; r2 += lambda1 * sqrt( tmp2) ; 
        if( r1 < r2 ) r =m2 ; 
        else l = m1 ; 
    }        
    return l ; 
}
    
bool solve_glasso_sub_sub(const vector<VD>&  g_input, const VD &g_output , VD& g_A, const VD& B, 
        const double & lambda1 , const double& lambda2,  const vector< int> group,  set< int >& feature_set ){
    VD sub_output = g_output ;   
    map< int , bool > mp; 
    for(int i=0; i < SZ(group) ; ++i) 
        mp[group[i]] = 1 ; 
    for(int i=0 ; i < SZ(g_input); ++i){
        
        for(set< int > ::iterator it = feature_set.begin() ;  it !=feature_set.end() ; ++it){
            int j = *it; 
            if( mp.find(j) != mp.end()) continue ;
            sub_output[i] -= g_input[i][j] * g_A[j]; 
        }

    }
    
    VD alpha(SZ(group)) ;  
    for(int i=0 ; i < SZ(group) ; ++i){ 
       alpha[ i ] = 0 ; 
       for(int j= 0 ; j < SZ(g_input) ;++j) {
           alpha[ i ] += g_input[j][group[i]] * sub_output[j];  
       }
    }
    double J = 0 ; 
    for(int i=0 ; i < SZ(group) ; ++i)
        J = J +pow( alpha[i]  - lambda2* cal_t_hat( alpha[i], lambda2),2); 
    J /= lambda1*lambda1; 
    if( J < 1 + eps){ 
        for(int i=0; i < SZ(group); ++i){ 
            g_A[ group[i]] = 0.  ; 
        }
        return true ;
    }
    
    VD theta( SZ(group)), pre_theta(SZ(group))  ; 
    for(int i=0 ; i < SZ(group) ; ++i){ 
        theta[i] = g_A[ group[i]] ;
    }
    int cnt = 0 ;  
    while( cnt < 100 ){
        ++cnt; 
        
        if( cnt % 100 == 0 ){
            for(int i=0 ; i < SZ(pre_theta) ; ++i){ 
                cout <<pre_theta[i] << " " ; 
            }
            cout << endl ; 
            for(int i=0 ; i < SZ(theta) ; ++i){ 
                 cout <<theta[i] << " ";  
            }
            cout << endl ; 
            pre_theta = theta ;  
            cout << "Inside Group: "<< cnt << endl ;
        }
        pre_theta = theta ;
        for(int now=0; now < SZ(group) ; ++now){  
            VD w  = sub_output ; 
            for(int i=0 ; i < SZ(g_input) ; ++i){ 

                for(int j=0 ; j < SZ(group) ; ++j){ 
                    if( j == now ) continue ;  
                    w[i] -= g_input[i][ group[j] ] * theta[j] ;  
                }
            }
            double tt = 0  ;
            for(int i=0 ; i < SZ(g_input) ; ++i){
                tt += w[i] * g_input[i][group[now]] ; 
            }
            if( fabs( tt ) < lambda2 ) { 
                theta[now] = 0 ; 
                continue ; 
            }
            else theta[now] = ternary_search( g_input, w , g_A, B , lambda1, lambda2, group,now,theta) ;  
             
        }
        if( converge( pre_theta, theta ) ) break ; 
    }
    double sum = 0 ; 
    for(int i=0; i < SZ(group) ; ++i){ 
        g_A[ group[i]] = theta[i];
        sum += fabs(theta[i]); 
    }
    return fabs(sum) < eps; 
   
}
static int times = 0 ;
void solve_glasso_sub(const vector < VD > &g_input, const VD &g_output, VD & g_A, const VD& B , const double & lambda1,const double & lambda2,
        const vector< vector<int> >& group){ 
    if( !times )
    for(int i=0; i < SZ(g_A) ; ++i)
       g_A[i] = rand_double();
    ++times;
    int g_n = SZ(g_input), g_p = SZ(g_input[0]), g_num = SZ(group)  ; 
    VD preA  ; 
    int cnt =0 ;  
 /*   for(int i=0 ; i < SZ(group) ; ++i){ 
        double sum = 0 ;
        for(int j=0 ; j < SZ(group[i]) ; ++j){ 
            sum += fabs(g_A[ group[i][j] ] ) ; 
        }
        if( sum < eps) state[i] = true ;
    }
*/

    set < int > group_set ;  
    set < int > feature_set; 
    for(int i=0 ; i < SZ(group); ++i){
       group_set.insert(i) ; 
       for(int j=0 ; j <SZ(group[i]); ++j){ 
           feature_set.insert( group[i][j] ) ; 
       }
    }
    
    do{
        preA =g_A ; 
        
        if( cnt % 100 == 0) cout << "Solve Groups: "<< cnt << endl ;
        vector< int > del ; 
        for(set< int > ::iterator it = group_set.begin() ; it!=group_set.end() ; ++it){
            int x = *it ; 
            if(  solve_glasso_sub_sub( g_input , g_output , g_A , B , lambda1, lambda2*B[x], group[x],feature_set) ){
                del.push_back(x) ; 
                for(int i=0 ; i < SZ(group[x]); ++i)
                    feature_set.erase(group[x][i]);
            }
        }      
        for( int i =0 ; i < SZ(del) ; ++i){ 
            group_set.erase( del[i] ) ; 
           /* for(int j=0; j < SZ(group[ del[i]]) ; ++j){ 
                feature_set.erase( group[del[i]][j] ) ; 
            }*/
        }
        ++cnt; 
    }while(cnt < 100 && !converge(preA,g_A)) ; 
}

void update_B( const VD& g_A, VD& B){
    int p = SZ(B) , q = SZ(g_A)/SZ(B) ; 
    double average = 0 ; 
    for(int i= 0 ; i < SZ(g_A); ++i)
        average += fabs(g_A[i]); 

    average /= p ; 
    for(int i=0 ; i < p ; ++i){ 
        B[i]=0 ; 
        for(int j=0 ; j < q ; ++j){ 
            B[i] = B[i] + fabs(g_A[ j * p + i ]); 
        }
        B[i]= ( 1 + average ) / ( 1 + B[i] ) ; 
    }
}

void solve_glasso(const vector < VD > &g_input , const VD &g_output , VD& g_A ,  VD & B , const double & lambda1 , const double& lambda2 ){ 
    vector< vector<int> > group(SZ(B)) ; 
    int p = SZ(B) ; 
       
    for(int i=0; i < SZ(B) ; ++i){
        int q = SZ(g_A)/SZ(B); 
        group[i].clear() ; 
        for(int j=0 ; j < q ; ++j){ 
            group[i].push_back( j * p + i) ; 
        }
        //cout << SZ(group[i]) << endl ;
        
        //exit(0) ;
    }

    for(int i=0; i < SZ(B); ++i){
        B[i]=1.;
    }
    int cnt = 0 ; 
    VD preB  ; 
    do{
        preB = B ; 
        solve_glasso_sub( g_input , g_output , g_A ,B , lambda1, lambda2, group) ; 
        update_B( g_A , B );
        ++cnt; 
        if( cnt % 10 == 0 ) cout <<"Solve Glasso: "<< cnt << endl ;  
    }while( cnt < 100 && !converge(preB , B ) ) ; 
}
 
void train_amtrfs(const vector< VD > & input,  const vector< VD > & output ,  
        vector< VD > & A  , VD B , double & lambda1, double&  lambda2){
    int n = SZ(input), p = SZ(input[0]), q = SZ(output[0]);
   // cout << p << " "<< q << endl; 
   // exit(0) ; 
    vector < VD > g_input(n*q); 
    VD g_output(n*q);  
    VD g_A( p * q ) ; 
    format_input( input  , g_input ) ; 
    format_output( output , g_output ) ; 
    solve_glasso( g_input , g_output ,g_A, B , lambda1, lambda2) ;  
    format_A( g_A , A) ;   
}
#endif
