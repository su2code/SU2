#include "catch.hpp"

#include <sstream>
#include <stdio.h>

#include "../../../Common/include/CConfig.hpp"
#include "../../../Common/include/embed_python.hpp"

TEST_CASE("Embed_python", "[python_embedding]")
{   
    //Embed_Python result();

    su2double x_var=2.2;
    su2double a_var=3.3;
    su2double b_var=4.4;
    su2double val_var;
    su2double ans_var;
    
    Embed_Python sqaure(x_var, &val_var);
    CHECK(val_var==Approx(4.84));

    Embed_Python addition(a_var,b_var,&ans_var);
    CHECK(ans_var==Approx(7.7));


    

}