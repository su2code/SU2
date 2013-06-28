// TAPENADE friendly version of pow()
double pow(double a, double b){

    double c;
    double epsilon;

    epsilon = 1.0e-10;

    if ((fabs(a) < epsilon) && (a >= 0.0)) {
        c = exp( b * log(a + epsilon) );
    } else if ((fabs(a) < epsilon) && (a < 0.0)) {
        c = exp( b * log(a - epsilon) );
    } else {
        c = exp( b * log(a) );
    }

    return c;

}
