double fabs(double a){

    double b;
    double epsilon;

    epsilon = 1e-10;

    b = sqrt(a*a + epsilon*epsilon);

    return b;

}
