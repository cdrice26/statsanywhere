#include "t.h"

double tpdf(double t, int df) {
    return gamma((df + 1.0) / 2.0) / (sqrt(M_PI * df) * gamma(df / 2.0)) * pow(1 + pow(t, 2.0) / df, -(df + 1.0) / 2.0);
}

int main() {
    printf("(Should be .0173) tpdf(3, 5) = %f", tpdf(3, 5));
}
