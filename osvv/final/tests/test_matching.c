/*
 * Standalone test for the matching (flow decomposition) routines in flow.c.
 * No MATLAB required.
 *
 * Builds a flow problem, runs hipr with route_flag = 1, and checks that the
 * returned matching
 *   - carries exactly the max-flow value,
 *   - only pairs a source-side vertex with a sink-side vertex,
 *   - respects every vertex's source/sink capacity.
 *
 * Usage: test_matching {random|grid} SIZE SEED ALGO
 *   random: SIZE vertices, 5*SIZE random undirected edges; first half of the
 *           vertices attached to the source, second half to the sink.
 *   grid:   SIZE x SIZE grid; left column attached to the source, right column
 *           to the sink, so every flow path has length >= SIZE (long paths).
 *   ALGO:   0 = dinic, 1 = dynamic trees, -1 = max flow only (no matching;
 *           use as a baseline to time the matching step)
 * Prints one line ending in OK or MISMATCH; exit code 1 on MISMATCH.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include "flow.h"

static int64_t *tl, *hd, *w, k, cap_mx;

static void arc(int64_t a, int64_t b, int64_t c) {
    if (k >= cap_mx) {
        cap_mx *= 2;
        tl = realloc(tl, sizeof(int64_t) * cap_mx);
        hd = realloc(hd, sizeof(int64_t) * cap_mx);
        w = realloc(w, sizeof(int64_t) * cap_mx);
    }
    tl[k] = a; hd[k] = b; w[k] = c; k++;
}

int main(int argc, char **argv) {
    if (argc != 5) {
        fprintf(stderr, "usage: %s {random|grid} SIZE SEED ALGO\n", argv[0]);
        return 2;
    }
    int grid = strcmp(argv[1], "grid") == 0;
    int64_t sz = atoll(argv[2]), algo = atoll(argv[4]);
    srand((unsigned) atoi(argv[3]));
    int64_t N = grid ? sz * sz : sz, s = N + 1, t = N + 2;
    char *side = calloc(N + 1, 1);                       /* 1 = source side, 2 = sink side */
    int64_t *vcap = calloc(N + 1, sizeof(int64_t));      /* source/sink capacity of each vertex */
    cap_mx = 1024; k = 0;
    tl = malloc(sizeof(int64_t) * cap_mx);
    hd = malloc(sizeof(int64_t) * cap_mx);
    w = malloc(sizeof(int64_t) * cap_mx);

    if (grid) {
#define ID(r, c) ((r) * sz + (c) + 1)
        for (int64_t r = 0; r < sz; r++) {
            int64_t c0 = 1 + rand() % 3, c1 = 1 + rand() % 3;
            arc(s, ID(r, 0), c0); side[ID(r, 0)] = 1; vcap[ID(r, 0)] = c0;
            arc(ID(r, sz - 1), t, c1); side[ID(r, sz - 1)] = 2; vcap[ID(r, sz - 1)] = c1;
            for (int64_t c = 0; c < sz; c++) {
                if (c + 1 < sz) { int64_t x = 1 + rand() % 4; arc(ID(r, c), ID(r, c + 1), x); arc(ID(r, c + 1), ID(r, c), x); }
                if (r + 1 < sz) { int64_t x = 1 + rand() % 4; arc(ID(r, c), ID(r + 1, c), x); arc(ID(r + 1, c), ID(r, c), x); }
            }
        }
    } else {
        for (int64_t i = 1; i <= N; i++) {
            int64_t c = 1 + rand() % 3;
            if (i <= N / 2) { arc(s, i, c); side[i] = 1; } else { arc(i, t, c); side[i] = 2; }
            vcap[i] = c;
        }
        for (int64_t e = 0; e < 5 * N; e++) {
            int64_t a = 1 + rand() % N, b = 1 + rand() % N, c = 1 + rand() % 4;
            if (a == b) continue;
            arc(a, b, c); arc(b, a, c);
        }
    }

    int64_t *oset = NULL, *mh = NULL, *mt = NULL, *mw = NULL, ne = 0, fl = 0;
    clock_t t0 = clock();
    hipr(N + 2, k, tl, hd, w, s, t, &oset, &mh, &mt, &mw, &ne, &fl, algo >= 0, algo >= 0 ? algo : 0);
    double secs = (double) (clock() - t0) / CLOCKS_PER_SEC;

    int64_t tot = 0, bad = 0, *load = calloc(N + 1, sizeof(int64_t));
    for (int64_t i = 0; i < ne; i++) {
        if (mw[i] <= 0 || mh[i] < 1 || mh[i] > N || mt[i] < 1 || mt[i] > N) { bad++; continue; }
        tot += mw[i];
        load[mh[i]] += mw[i];
        if (side[mh[i]] == 0 || side[mt[i]] == 0 || side[mh[i]] == side[mt[i]]) bad++;
    }
    for (int64_t i = 1; i <= N; i++)
        if (load[i] > vcap[i]) bad++;
    /* the matching is returned symmetrically (u,v) and (v,u) */
    int ok = algo < 0 ? 1 : (tot % 2 == 0) && tot / 2 == fl && bad == 0;
    printf("%-6s size=%-8ld algo=%-7s flow=%-9ld matched=%-9ld violations=%ld hipr+matching=%.3fs %s\n",
           argv[1], (long) sz, algo < 0 ? "none" : algo ? "dynamic" : "dinic", (long) fl, (long) (tot / 2), (long) bad, secs,
           ok ? "OK" : "MISMATCH");
    return ok ? 0 : 1;
}
