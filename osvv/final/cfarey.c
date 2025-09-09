#include <math.h>
#include <stdint.h>

void farey(int64_t num_in, int64_t den_in, int64_t p, int64_t *num_out, int64_t *den_out) {
    int64_t error;
    int64_t a, b, c, d;
    a = 0;
    b = 1;
    c = ceil(((double) num_in) / den_in);
    d = 1;
    int64_t h, k;

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

