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
        if( p < 10 ) cout <<av<< " "<< var<<endl ;
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
    assert( 7 == SZ(parameter) ) ;      
      
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
           // if( ++cnt > 500 ) continue ;
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
        int cnt =0 ; 
        while( sin >> x  ){ 
       //    if( ++cnt == 1 )
           tmp.push_back(x) ; 
        }
        output.push_back(tmp) ; 
    }
    
    fclose(fp) ;  
    init( input ) ; 
    init( output ) ; 
    cout <<"fuck"<<endl;  
    int rand_seed = 0;
    srand( rand_seed =  atoi( parameter[6].second.c_str() ) );
  
    for(int i = 0 ; i < SZ(input) ; ++i){ 
       int sw=  rand()%(i+1) ; 
       swap( input[i], input[ sw ] ) ; 
       swap( output[i], output[sw] ) ;
    } 
    for(int i=0 ; i < SZ(input) ; ++i){ 
        for(int j=0 ; j < SZ(input[i]) ;++j){ 
        //    cout << input[i][j]<< " " ; 
        }
       // cout << endl ;
    }

    for(int i=0 ; i < SZ(output); ++i){
        for(int j=0 ;j <SZ(output[i]);++j){
//            cout << output[i][j] << " " ; 
        }
  //      cout << endl ; 
    }
    cout <<"Read Finished" << endl ; 
    //exit(0) ; 
    for(int i=0 ; i < SZ(input); ++i){
        input[i].push_back(0.1) ; 
    }
    double lambda1 , lambda2 ; 
    vector < VD > B(SZ(input[0]))  ; 
    VD A(SZ(B)) ;  
    VD tmp( SZ(output[0])) ; 
    for(int i=0 ; i < SZ(B) ;++i){
        B[ i ] = tmp ;  
    }
    lambda1 = atof( parameter[3].second.c_str() ) ; 
    lambda2 = atof( parameter[4].second.c_str() ) ;
 /*   VD Abun(SZ(input[0])+1)  ;  

    fp = fopen(parameter[5].second.c_str(),"r") ;  
    for(int i=0 ; i < SZ(Abun)-1; ++i)
        fscanf(fp, "%lf",  &Abun[i]);
    fclose(fp); 
    Abun[SZ(Abun)-1]= 1 ;   
*/
    int test = atoi( parameter[5].second.c_str() ) ; 
    vector< VD > train_input(SZ(input) - test ) , train_output( SZ(output) - test) ;  
    vector< VD > test_input(test) , test_output(test) ;
    for(int i = 0; i < SZ(input)-test ; ++i){
        train_input[i] = input[i];
        train_output[i] = output[i]; 
    }
    for(int i= SZ(input) - test  ; i < SZ(input) ; ++i){ 
        test_input[i-(SZ(input) - test)] = input[i]; 
        test_output[i-(SZ(output) - test)] = output[i] ; 
    }
    
    string prefix0 = parameter[0].second + parameter[6].second ;  
    string prefix1 = parameter[1].second + parameter[6].second ; 
    FILE *fp1  = fopen( (prefix0 + ".train.txt").c_str(), "w"); 
    FILE *fp2  = fopen( (prefix1 + ".train.txt").c_str(),"w") ; 
    for(int i=0; i < SZ(train_input) ; ++i){
        for(int j=0 ; j < SZ(train_input[0]) ;++j){ 
            fprintf( fp1, "%lf ", train_input[i][j]) ; 
        }
        for(int j=0 ; j < SZ(train_output[0]) ; ++j){ 
            fprintf( fp2, "%lf ", train_output[i][j]) ;
        }

        fprintf(fp1, "\n") ; 
        fprintf(fp2, "\n") ;
    }
    fclose(fp1) ;
    fclose(fp2) ;
    
    fp1 = fopen( (prefix0 + ".test.txt").c_str() , "w" ) ; 
    fp2 = fopen( (prefix1 + ".test.txt").c_str() , "w" ) ;
    
    for(int i=0 ; i < SZ(test_input) ;++i){
        for(int j=0 ; j < SZ(test_input[0]) ; ++j){
            fprintf( fp1, "%lf " , test_input[i][j]) ;
        }
        for(int j=0 ; j < SZ(test_output[0]) ;++j){ 
            fprintf(fp2, "%lf " ,test_output[i][j] ); 
        } 
        fprintf( fp1, "\n") ; 
        fprintf( fp2, "\n") ;
    }
    fclose(fp1); 
    fclose(fp2); 
    model_amtrfs(train_input , train_output, B , A , lambda1 , lambda2) ;   
    cout <<"I am am am am" <<endl; 
    cout <<predict_amtrfs(test_input , test_output , B , A , lambda1 , lambda2, true) <<endl ;  
    fp = fopen( parameter[2].second.c_str() , "w") ;  
    assert( NULL != fp ) ;  
    for(int i = 0 ; i < SZ(B) ; ++i){
        for(int j = 0 ; j < SZ(B[0]) ; ++j){
            fprintf( fp , "%lf " , B[i][j] ) ; 
        }
        fprintf( fp , "\n" ) ; 
    }
    fclose(fp) ; 
    clock_t end = clock() ; 
    cout << ( end - start ) / CLOCKS_PER_SEC<<" Seconds" << endl ; 
    return 0 ;
}

    


