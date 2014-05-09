#include <cstdio>
#include <string>
#include <cstring>
#include <sstream>
#include <algorithm>

using namespace std ;  

int main (){
    for(int i=1; i <= 50; ++i){
        string s = "./run " ; 
        string t = ""; 
        int x = i ; 
        while(x){
            t = t + (char) (x % 10+'0') ;   
            x/=10;
        }
        reverse( t.begin(),t.end()) ;  
        s = s + t ; 
        s = s + " &";
        system( s.c_str() ) ; 
    }
    return 0; 
}
