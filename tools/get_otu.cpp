#include <cstdio>
#include <cstring>
#include <vector>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <cstdlib>

using namespace std ; 

#define SZ(x) (int)x.size()
double data[ 500][ 1000]; 
char ss[1000000+10] ;

double corr(double x[],double y[], int n){
    double meanx =0 , meany =0 ; 
    for(int i=0 ; i < n ; ++i){
        meanx += x[ i ];
        meany += y[ i ]; 
    }
    meanx /= n , meany/=n ; 
    double varx = 0 , vary = 0 ; 
    double covar = 0 ; 
    for(int i=0; i < n ; ++i){
        x[i]-=meanx;
        y[i]-=meany; 
        varx += x[i] * x[i]; 
        vary += y[i] * y[i];
        covar+= x[i] * y[i];
    }
    if(abs(varx) <1e-9 || abs(vary) < 1e-9) return 0; 
    return covar /sqrt(varx*vary); 
}

int main (){
    
    FILE* fp = fopen("id_500.txt","r"); 
    char c; 
    int a; 
    vector< int > vc; 
    while(fscanf (fp,"%c%d\n",&c,&a)==2){
        vc.push_back( a +1 ) ; 
    }
    fclose(fp); 
    
    fp = fopen("parameter_5001000.000000_1.000000_1.txt","r"); 
    double x, y, z; 
    FILE *fp2 = fopen("otu_1000.000000_1.000000_1.txt","w") ; 
    int id = -1 ; 
    vector < int > tmp ; 
    while( fscanf(fp,"%lf%lf%lf", &x,&y,&z) ==3){
       ++id;
       if( abs(x)+abs(y)+abs(z) < 1e-4) continue ;  
       tmp.push_back(id) ;  
       fprintf(fp2,"%d %lf %lf %lf\n", vc[id], x, y ,z ); 
    } 
    fclose(fp) ; 
    fclose(fp2) ; 


    fp = fopen("input_500.txt","r") ;  
    fp2= fopen("otu_matrix.txt","w") ; 
    int samples = 0 ; 
    while( fgets (ss, 1000000, fp ) != NULL ){   
        
        string  s= string( ss ) ;       
        stringstream sin ( s  ) ;  
        double x ; 
        int t = 0 ;
        while( sin >> x ){ 
            data[ samples ][ t ++ ] = x  ;  
        }
        ++samples ;  
    }
    fclose(fp); 
   // fprintf(fp2, "#continuous " ) ;  
    
    for(int i=0 ; i < SZ(tmp ) ; ++i){ 
      //  fprintf(fp2, "OTU%d\t", vc[  tmp[i] ] ) ; 
    }
//    fprintf(fp2, "\n") ; 
   // for(int i=0 ; i < samples;  ++i){ 
     //   fprintf(fp2, "EXP%d ", i) ;  
   // }
   // fprintf(fp2, "\n") ;  

   /* for(int i = 0; i < min(20,SZ(tmp)) ; ++i){
        fprintf(fp2, "OTU%d ", vc[tmp[i]]) ;  
        for(int j=0; j < samples; ++j){  
            fprintf(fp2, "%lf ", 1000*data[ j ] [ tmp[i] ] ) ; 
        }
        fprintf(fp2, "\n") ; 
    }*/
    FILE *otu_name = fopen("otu_name.txt","w");  
    for(int i = 0 ; i < SZ(tmp); ++i){
        fprintf(otu_name, "%d\n", vc[ tmp[i]]);
    }
    fclose(otu_name); 

    for(int i= 0 ; i < samples; ++i){ 
        for(int j=0 ; j <  SZ(tmp) ; ++j){ 
                fprintf(fp2,"%.2f\t", 10000 * data[ i ][ tmp[j]]) ; 
        }
        fprintf(fp2, "\n") ; 
    }
    
    fclose(fp2) ; 

    fp2 = fopen("correlation.txt","w") ;  
    for(int i=0 ; i < SZ(tmp) ; ++i){
        for(int j=1+i ; j < SZ(tmp) ; ++j){ 
            double x[500] ,y[500];
            for(int k =0 ; k < samples; ++k){
                x[k] = data[k][tmp[i]]; 
                y[k] = data[k][tmp[j]];
            }
            double ans = corr(x,y,samples) ; 
           // if(fabs(ans)>0.7){
            fprintf(fp2,"%d %d %lf\n",  vc[tmp[i]] , vc[tmp[j]] ,ans);
          /*  for(int k=0; k < samples ; ++k){
               fprintf(fp2,"%lf ", x[k]);
            }
            fprintf(fp2, "\n"); 
            for(int k=0; k < samples ; ++k){
                fprintf(fp2,"%lf ", y[k]); 
            }*/
            //fprintf(fp2, "\n");
            //}
        }
    }
    fclose(fp2);
    return 0;
}


           
