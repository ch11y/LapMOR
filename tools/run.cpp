#include <cstdio>
#include <cstring>
#include <map>
#include <string>
#include <iostream>
#include <cstdlib>

using namespace std ;  

#define SZ(x) (int)x.size()

int main (int argc, char ** argv){ 
     
 
 //   system("mkdir result");  
    
     string seed = argv[1] ;
//    string seed = "1";
     for(double lambda1 = 10000 ; lambda1 <= 10000 ; lambda1  += 100){
        for(double lambda2= 0.1 ; lambda2 <= 1; lambda2 += 0.1){ 

            
            string cmd = "./AMTRFS -i input_1000_genus.txt -z ";
            cmd += seed ; 
            cmd += " -o output_1000_genus.txt  -p parameter_1000_genus" ; 

            char tmp1[20]; 
            sprintf(tmp1, "%lf", lambda1) ; 
            char tmp2[20];
            sprintf(tmp2, "%lf", lambda2) ; 
            string x= string(tmp1) ;                    
            x += "_" ;  
            x += string(tmp2) ; 
            x += "_" ; 
            x += seed ; 
           // x +=".txt " ;
            cmd += x ; 
            cmd += ".txt" ; 

            cmd += " -p1 " ; 
            cmd += string(tmp1) ; 
            cmd += " -p2 " ; 
            cmd += string(tmp2) ; 

            cmd += " -t 71 " ; 
            x = ">log_1000_genus" + x ;   
           // x += seed ; 
           
            cmd += x ; 
            cmd +=".txt &" ;  

           // cmd += "" ; 
            cout << cmd << endl  ;
            system(cmd.c_str()) ; 
        }
    }
    return 0; 
}
