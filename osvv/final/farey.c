#include <math.h>

void farey(long num_in, long den_in, long p, long *num_out, long *den_out) {
    long error;
    long a, b, c, d;
    a = 0;
    b = 1;
    c = ceil(((double) num_in) / den_in);
    d = 1;
    long h, k;

    for (;;) {
        h = a + c;
        k = b + d;

        if (h > p) {
            h = c;
            k = d;
            break;
        }

        error = h * den_in - k * num_in; /* error is positive if h/k > num/den */

        if (error == 0)
            break;

        if (error < 0) {
            a = h;
            b = k;
        }

        if (error > 0) {
            c = h;
            d = k;
        }
    }
    *num_out = h;
    *den_out = k;
}

