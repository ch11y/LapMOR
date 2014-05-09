#include <cstdio> 
#include <cstdlib>
#include <map>
#include <iostream>
#include <algorithm>
#include <vector>
using namespace std ;

#define SZ(x) (int)x.size()
int main (){
    char ss[1000]; 
    double a,b,c,d;  
    FILE *fp1 = fopen("log.txt","r");
    FILE *fp2 = fopen("points.txt","w"); 
    vector< pair< pair<double,double> , pair<double,double> > > vc ; 
    while( fscanf (fp1, "%s%lf%lf%lf%lf" , ss, &a, &b, &c,&d)==5){ 
       // fprintf(fp2,"%lf %lf %lf %lf\n", a,b,c,d) ;  
        vc.push_back(make_pair( make_pair(a,b), make_pair(c,d))) ; 
    }
    sort(vc.begin(),vc.end()); 
    for(int i=0 ; i < SZ(vc); ++i){
       fprintf(fp2, "%lf %lf %lf %lf\n", vc[i].first.first, vc[i].first.second, vc[i].second.first, vc[i].second.second) ;
    }

    return 0;
}

