// Exact pow_d()
double pow_d(double a, double ad, double b, double *c) {
    double cd;

    c[0] = pow(a, b);
    cd = b*ad*pow(a, (b-1.0));
    return cd;
}
