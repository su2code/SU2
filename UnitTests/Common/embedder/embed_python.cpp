#include <utility>
#include "embed_python.hpp"

using namespace std;
/*Embed_Python(){

};
*/

su2double Embed_Python::sqaure(su2double x, su2double* val){

su2double  sq=x*x;

*val=sq;

exit_code=0;
return exit_code;
    
 }

su2double Embed_Python::addition(su2double a,su2double b, su2double* ans){
 su2double  sum= a+b;

*ans=sum;

exit_code=0;
return exit_code;
    
 }