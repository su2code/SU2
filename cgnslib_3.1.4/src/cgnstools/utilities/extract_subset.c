/*
 * convert_location.c - convert between vertex and cell-center
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _WIN32
# define unlink _unlink
#else
# include <unistd.h>
#endif

#include "getargs.h"
#include "cgnslib.h"
#include "cgnsutil.h"

#ifndef CG_MODE_MODIFY
# define CG_MODE_MODIFY MODE_MODIFY
#endif

/* command line options */

static char options[] = "fb:z:ijkawB:";

static char *usgmsg[] = {
    "usage  : extract_subset [options] CGNSfile [newCGNSfile]",
    "options:",
    "   -b<base> = use CGNS base number <base> (default 1)",
    "   -z<zone> = read zone number <zone> (default all)",
    "   -i       = subset in I-direction",
    "   -j       = subset in J-direction",
    "   -k       = subset in K-direction",
    "   -a       = subset in all directions (default)",
    "   -w       = use volume weighting for solution",
    "   -B<name> = write subset to base <name> (default same as read)",
    NULL
};

static int weighting = 0;

/*-------------------------------------------------------------------*/

static void field_subset (ZONE *z, SOLUTION *s, int inc[3],
                          double *f, double *w)
{
    cgsize_t i, j, k, ip, jp, kp;
    cgsize_t nn, nf, nw, np;
    ZONE zn;
    double *fn, *wn;

    for (np = 1, nn = 0; nn < 3; nn++) {
        zn.dim[nn] = (z->dim[nn] - 1) / inc[nn] + 1;
        np *= (zn.dim[nn] - 1 + s->rind[nn][0] + s->rind[nn][1]);
    }
    fn = (double *) malloc ((size_t)(2 * np) * sizeof(double));
    if (NULL == fn) FATAL (NULL, "malloc failed for solution field");
    wn = fn + np;
    for (nn = 0; nn < np; nn++) {
        fn[nn] = 0.0;
        wn[nn] = 0.0;
    }

    /* interior */

    for (k = 1; k < z->dim[2]; k += inc[2]) {
        kp = k + 1 < z->dim[2] ? k + 1 : k;
        for (j = 1; j < z->dim[1]; j += inc[1]) {
            jp = j + 1 < z->dim[1] ? j + 1 : j;
            for (i = 1; i < z->dim[0]; i += inc[0]) {
                ip = i + 1 < z->dim[0] ? i + 1 : i;
                nn = solution_index (&zn, s, (i - 1) / inc[0] + 1,
                    (j - 1) / inc[1] + 1, (k - 1) / inc[2] + 1);

                nf = solution_index (z, s, i, j, k);
                nw = cell_index (z, i, j, k);
                fn[nn] += f[nf] * w[nw];
                wn[nn] += w[nw];

                nf = solution_index (z, s, ip, j, k);
                nw = cell_index (z, ip, j, k);
                fn[nn] += f[nf] * w[nw];
                wn[nn] += w[nw];

                nf = solution_index (z, s, i, jp, k);
                nw = cell_index (z, i, jp, k);
                fn[nn] += f[nf] * w[nw];
                wn[nn] += w[nw];

                nf = solution_index (z, s, ip, jp, k);
                nw = cell_index (z, ip, jp, k);
                fn[nn] += f[nf] * w[nw];
                wn[nn] += w[nw];

                nf = solution_index (z, s, i, j, kp);
                nw = cell_index (z, i, j, kp);
                fn[nn] += f[nf] * w[nw];
                wn[nn] += w[nw];

                nf = solution_index (z, s, ip, j, kp);
                nw = cell_index (z, ip, j, kp);
                fn[nn] += f[nf] * w[nw];
                wn[nn] += w[nw];

                nf = solution_index (z, s, i, jp, kp);
                nw = cell_index (z, i, jp, kp);
                fn[nn] += f[nf] * w[nw];
                wn[nn] += w[nw];

                nf = solution_index (z, s, ip, jp, kp);
                nw = cell_index (z, ip, jp, kp);
                fn[nn] += f[nf] * w[nw];
                wn[nn] += w[nw];
            }
        }
    }

    /* rind planes */

    for (i = 1; i <= s->rind[0][0]; i++) {
        ip = 1 - i;
        for (k = 1; k < z->dim[2]; k += inc[2]) {
            kp = k + 1 < z->dim[2] ? k + 1 : k;
            for (j = 1; j < z->dim[1]; j += inc[1]) {
                jp = j + 1 < z->dim[1] ? j + 1 : j;
                nn = solution_index (&zn, s, ip,
                    (j - 1) / inc[1] + 1, (k - 1) / inc[2] + 1);
                fn[nn] += f[solution_index (z, s, ip, j, k)] +
                          f[solution_index (z, s, ip, jp, k)] +
                          f[solution_index (z, s, ip, j, kp)] +
                          f[solution_index (z, s, ip, jp, kp)];
                wn[nn] += 4.0;
            }
        }
    }

    for (i = 1; i <= s->rind[0][1]; i++) {
        ip = z->dim[0] + i - 1;
        for (k = 1; k < z->dim[2]; k += inc[2]) {
            kp = k + 1 < z->dim[2] ? k + 1 : k;
            for (j = 1; j < z->dim[1]; j += inc[1]) {
                jp = j + 1 < z->dim[1] ? j + 1 : j;
                nn = solution_index (&zn, s, zn.dim[0] + i - 1,
                    (j - 1) / inc[1] + 1, (k - 1) / inc[2] + 1);
                fn[nn] += f[solution_index (z, s, ip, j, k)] +
                          f[solution_index (z, s, ip, jp, k)] +
                          f[solution_index (z, s, ip, j, kp)] +
                          f[solution_index (z, s, ip, jp, kp)];
                wn[nn] += 4.0;
            }
        }
    }

    for (j = 1; j <= s->rind[1][0]; j++) {
        jp = 1 - j;
        for (k = 1; k < z->dim[2]; k += inc[2]) {
            kp = k + 1 < z->dim[2] ? k + 1 : k;
            for (i = 1; i < z->dim[0]; i += inc[0]) {
                ip = i + 1 < z->dim[0] ? i + 1 : i;
                nn = solution_index (&zn, s, (i - 1) / inc[0] + 1,
                    jp, (k - 1) / inc[2] + 1);
                fn[nn] += f[solution_index (z, s, i, jp, k)] +
                          f[solution_index (z, s, ip, jp, k)] +
                          f[solution_index (z, s, i, jp, kp)] +
                          f[solution_index (z, s, ip, jp, kp)];
                wn[nn] += 4.0;
            }
        }
    }

    for (j = 1; j <= s->rind[1][1]; j++) {
        jp = z->dim[1] + j - 1;
        for (k = 1; k < z->dim[2]; k += inc[2]) {
            kp = k + 1 < z->dim[2] ? k + 1 : k;
            for (i = 1; i < z->dim[0]; i += inc[0]) {
                ip = i + 1 < z->dim[0] ? i + 1 : i;
                nn = solution_index (&zn, s, (i - 1) / inc[0] + 1,
                    zn.dim[1] + j - 1, (k - 1) / inc[2] + 1);
                fn[nn] += f[solution_index (z, s, i, jp, k)] +
                          f[solution_index (z, s, ip, jp, k)] +
                          f[solution_index (z, s, i, jp, kp)] +
                          f[solution_index (z, s, ip, jp, kp)];
                wn[nn] += 4.0;
            }
        }
    }

    for (k = 1; k <= s->rind[2][0]; k++) {
        kp = 1 - k;
        for (j = 1; j < z->dim[1]; j += inc[1]) {
            jp = j + 1 < z->dim[1] ? j + 1 : j;
            for (i = 1; i < z->dim[0]; i += inc[0]) {
                ip = i + 1 < z->dim[0] ? i + 1 : i;
                nn = solution_index (&zn, s, (i - 1) / inc[0] + 1,
                    (j - 1) / inc[1] + 1, kp);
                fn[nn] += f[solution_index (z, s, i, j, kp)] +
                          f[solution_index (z, s, ip, j, kp)] +
                          f[solution_index (z, s, i, jp, kp)] +
                          f[solution_index (z, s, ip, jp, kp)];
                wn[nn] += 4.0;
            }
        }
    }

    for (k = 1; k <= s->rind[2][1]; k++) {
        kp = z->dim[2] + k - 1;
        for (j = 1; j < z->dim[1]; j += inc[1]) {
            jp = j + 1 < z->dim[1] ? j + 1 : j;
            for (i = 1; i < z->dim[0]; i += inc[0]) {
                ip = i + 1 < z->dim[0] ? i + 1 : i;
                nn = solution_index (&zn, s, (i - 1) / inc[0] + 1,
                    (j - 1) / inc[1] + 1, zn.dim[2] + k - 1);
                fn[nn] += f[solution_index (z, s, i, j, kp)] +
                          f[solution_index (z, s, ip, j, kp)] +
                          f[solution_index (z, s, i, jp, kp)] +
                          f[solution_index (z, s, ip, jp, kp)];
                wn[nn] += 4.0;
            }
        }
    }

    /* rind lines */

    for (i = 1; i <= s->rind[0][0]; i++) {
        ip = 1 - i;
        for (j = 1; j <= s->rind[1][0]; j++) {
            jp = 1 - j;
            for (k = 1; k < zn.dim[2]; k++) {
                nn = solution_index (&zn, s, ip, jp, k);
                nf = solution_index (&zn, s, ip+1, jp, k);
                fn[nn] += fn[nf] * wn[nf];
                wn[nn] += wn[nf];
                nf = solution_index (&zn, s, ip, jp+1, k);
                fn[nn] += fn[nf] * wn[nf];
                wn[nn] += wn[nf];
            }
        }
        for (j = 1; j <= s->rind[1][1]; j++) {
            jp = zn.dim[1] + j - 1;
            for (k = 1; k < zn.dim[2]; k++) {
                nn = solution_index (&zn, s, ip, jp, k);
                nf = solution_index (&zn, s, ip+1, jp, k);
                fn[nn] += fn[nf] * wn[nf];
                wn[nn] += wn[nf];
                nf = solution_index (&zn, s, ip, jp-1, k);
                fn[nn] += fn[nf] * wn[nf];
                wn[nn] += wn[nf];
            }
        }
        for (k = 1; k <= s->rind[2][0]; k++) {
            kp = 1 - k;
            for (j = 1; j < zn.dim[1]; j++) {
                nn = solution_index (&zn, s, ip, j, kp);
                nf = solution_index (&zn, s, ip+1, j, kp);
                fn[nn] += fn[nf] * wn[nf];
                wn[nn] += wn[nf];
                nf = solution_index (&zn, s, ip, j, kp+1);
                fn[nn] += fn[nf] * wn[nf];
                wn[nn] += wn[nf];
            }
        }
        for (k = 1; k <= s->rind[2][1]; k++) {
            kp = zn.dim[2] + k - 1;
            for (j = 1; j < zn.dim[1]; j++) {
                nn = solution_index (&zn, s, ip, j, kp);
                nf = solution_index (&zn, s, ip+1, j, kp);
                fn[nn] += fn[nf] * wn[nf];
                wn[nn] += wn[nf];
                nf = solution_index (&zn, s, ip, j, kp-1);
                fn[nn] += fn[nf] * wn[nf];
                wn[nn] += wn[nf];
            }
        }
    }

    for (i = 1; i <= s->rind[0][1]; i++) {
        ip = zn.dim[0] + i - 1;
        for (j = 1; j <= s->rind[1][0]; j++) {
            jp = 1 - j;
            for (k = 1; k < zn.dim[2]; k++) {
                nn = solution_index (&zn, s, ip, jp, k);
                nf = solution_index (&zn, s, ip-1, jp, k);
                fn[nn] += fn[nf] * wn[nf];
                wn[nn] += wn[nf];
                nf = solution_index (&zn, s, ip, jp+1, k);
                fn[nn] += fn[nf] * wn[nf];
                wn[nn] += wn[nf];
            }
        }
        for (j = 1; j <= s->rind[1][1]; j++) {
            jp = zn.dim[1] + j - 1;
            for (k = 1; k < zn.dim[2]; k++) {
                nn = solution_index (&zn, s, ip, jp, k);
                nf = solution_index (&zn, s, ip-1, jp, k);
                fn[nn] += fn[nf] * wn[nf];
                wn[nn] += wn[nf];
                nf = solution_index (&zn, s, ip, jp-1, k);
                fn[nn] += fn[nf] * wn[nf];
                wn[nn] += wn[nf];
            }
        }
        for (k = 1; k <= s->rind[2][0]; k++) {
            kp = 1 - k;
            for (j = 1; j < zn.dim[1]; j++) {
                nn = solution_index (&zn, s, ip, j, kp);
                nf = solution_index (&zn, s, ip-1, j, kp);
                fn[nn] += fn[nf] * wn[nf];
                wn[nn] += wn[nf];
                nf = solution_index (&zn, s, ip, j, kp+1);
                fn[nn] += fn[nf] * wn[nf];
                wn[nn] += wn[nf];
            }
        }
        for (k = 1; k <= s->rind[2][1]; k++) {
            kp = zn.dim[2] + k - 1;
            for (j = 1; j < zn.dim[1]; j++) {
                nn = solution_index (&zn, s, ip, j, kp);
                nf = solution_index (&zn, s, ip-1, j, kp);
                fn[nn] += fn[nf] * wn[nf];
                wn[nn] += wn[nf];
                nf = solution_index (&zn, s, ip, j, kp-1);
                fn[nn] += fn[nf] * wn[nf];
                wn[nn] += wn[nf];
            }
        }
    }

    for (j = 1; j <= s->rind[1][0]; j++) {
        jp = 1 - j;
        for (k = 1; k <= s->rind[2][0]; k++) {
            kp = 1 - k;
            for (i = 1; i < zn.dim[0]; i++) {
                nn = solution_index (&zn, s, i, jp, kp);
                nf = solution_index (&zn, s, i, jp+1, kp);
                fn[nn] += fn[nf] * wn[nf];
                wn[nn] += wn[nf];
                nf = solution_index (&zn, s, i, jp, kp+1);
                fn[nn] += fn[nf] * wn[nf];
                wn[nn] += wn[nf];
            }
        }
        for (k = 1; k <= s->rind[2][1]; k++) {
            kp = zn.dim[2] + k - 1;
            for (i = 1; i < zn.dim[0]; i++) {
                nn = solution_index (&zn, s, i, jp, kp);
                nf = solution_index (&zn, s, i, jp+1, kp);
                fn[nn] += fn[nf] * wn[nf];
                wn[nn] += wn[nf];
                nf = solution_index (&zn, s, i, jp, kp-1);
                fn[nn] += fn[nf] * wn[nf];
                wn[nn] += wn[nf];
            }
        }
    }

    for (j = 1; j <= s->rind[1][1]; j++) {
        jp = zn.dim[1] + j - 1;
        for (k = 1; k <= s->rind[2][0]; k++) {
            kp = 1 - k;
            for (i = 1; i < zn.dim[0]; i++) {
                nn = solution_index (&zn, s, i, jp, kp);
                nf = solution_index (&zn, s, i, jp-1, kp);
                fn[nn] += fn[nf] * wn[nf];
                wn[nn] += wn[nf];
                nf = solution_index (&zn, s, i, jp, kp+1);
                fn[nn] += fn[nf] * wn[nf];
                wn[nn] += wn[nf];
            }
        }
        for (k = 1; k <= s->rind[2][1]; k++) {
            kp = zn.dim[2] + k - 1;
            for (i = 1; i < zn.dim[0]; i++) {
                nn = solution_index (&zn, s, i, jp, kp);
                nf = solution_index (&zn, s, i, jp-1, kp);
                fn[nn] += fn[nf] * wn[nf];
                wn[nn] += wn[nf];
                nf = solution_index (&zn, s, i, jp, kp-1);
                fn[nn] += fn[nf] * wn[nf];
                wn[nn] += wn[nf];
            }
        }
    }

    /* rind points */

    for (i = 1; i <= s->rind[0][0]; i++) {
        ip = 1 - i;
        for (j = 1; j <= s->rind[1][0]; j++) {
            jp = 1 - j;
            for (k = 1; k <= s->rind[2][0]; k++) {
                kp = 1 - k;
                nn = solution_index (&zn, s, ip, jp, kp);
                nf = solution_index (&zn, s, ip+1, jp, kp);
                fn[nn] += fn[nf] * wn[nf];
                wn[nn] += wn[nf];
                nf = solution_index (&zn, s, ip, jp+1, kp);
                fn[nn] += fn[nf] * wn[nf];
                wn[nn] += wn[nf];
                nf = solution_index (&zn, s, ip, jp, kp+1);
                fn[nn] += fn[nf] * wn[nf];
                wn[nn] += wn[nf];
            }
            for (k = 1; k <= s->rind[2][1]; k++) {
                kp = zn.dim[2] + k - 1;
                nn = solution_index (&zn, s, ip, jp, kp);
                nf = solution_index (&zn, s, ip+1, jp, kp);
                fn[nn] += fn[nf] * wn[nf];
                wn[nn] += wn[nf];
                nf = solution_index (&zn, s, ip, jp+1, kp);
                fn[nn] += fn[nf] * wn[nf];
                wn[nn] += wn[nf];
                nf = solution_index (&zn, s, ip, jp, kp-1);
                fn[nn] += fn[nf] * wn[nf];
                wn[nn] += wn[nf];
            }
        }
        for (j = 1; j <= s->rind[1][1]; j++) {
            jp = zn.dim[1] + j - 1;
            for (k = 1; k <= s->rind[2][0]; k++) {
                kp = 1 - k;
                nn = solution_index (&zn, s, ip, jp, kp);
                nf = solution_index (&zn, s, ip+1, jp, kp);
                fn[nn] += fn[nf] * wn[nf];
                wn[nn] += wn[nf];
                nf = solution_index (&zn, s, ip, jp-1, kp);
                fn[nn] += fn[nf] * wn[nf];
                wn[nn] += wn[nf];
                nf = solution_index (&zn, s, ip, jp, kp+1);
                fn[nn] += fn[nf] * wn[nf];
                wn[nn] += wn[nf];
            }
            for (k = 1; k <= s->rind[2][1]; k++) {
                kp = zn.dim[2] + k - 1;
                nn = solution_index (&zn, s, ip, jp, kp);
                nf = solution_index (&zn, s, ip+1, jp, kp);
                fn[nn] += fn[nf] * wn[nf];
                wn[nn] += wn[nf];
                nf = solution_index (&zn, s, ip, jp-1, kp);
                fn[nn] += fn[nf] * wn[nf];
                wn[nn] += wn[nf];
                nf = solution_index (&zn, s, ip, jp, kp-1);
                fn[nn] += fn[nf] * wn[nf];
                wn[nn] += wn[nf];
            }
        }
    }

    for (i = 1; i <= s->rind[0][1]; i++) {
        ip = zn.dim[0] + i - 1;
        for (j = 1; j <= s->rind[1][0]; j++) {
            jp = 1 - j;
            for (k = 1; k <= s->rind[2][0]; k++) {
                kp = 1 - k;
                nn = solution_index (&zn, s, ip, jp, kp);
                nf = solution_index (&zn, s, ip-1, jp, kp);
                fn[nn] += fn[nf] * wn[nf];
                wn[nn] += wn[nf];
                nf = solution_index (&zn, s, ip, jp+1, kp);
                fn[nn] += fn[nf] * wn[nf];
                wn[nn] += wn[nf];
                nf = solution_index (&zn, s, ip, jp, kp+1);
                fn[nn] += fn[nf] * wn[nf];
                wn[nn] += wn[nf];
            }
            for (k = 1; k <= s->rind[2][1]; k++) {
                kp = zn.dim[2] + k - 1;
                nn = solution_index (&zn, s, ip, jp, kp);
                nf = solution_index (&zn, s, ip-1, jp, kp);
                fn[nn] += fn[nf] * wn[nf];
                wn[nn] += wn[nf];
                nf = solution_index (&zn, s, ip, jp+1, kp);
                fn[nn] += fn[nf] * wn[nf];
                wn[nn] += wn[nf];
                nf = solution_index (&zn, s, ip, jp, kp-1);
                fn[nn] += fn[nf] * wn[nf];
                wn[nn] += wn[nf];
            }
        }
        for (j = 1; j <= s->rind[1][1]; j++) {
            jp = zn.dim[1] + j - 1;
            for (k = 1; k <= s->rind[2][0]; k++) {
                kp = 1 - k;
                nn = solution_index (&zn, s, ip, jp, kp);
                nf = solution_index (&zn, s, ip-1, jp, kp);
                fn[nn] += fn[nf] * wn[nf];
                wn[nn] += wn[nf];
                nf = solution_index (&zn, s, ip, jp-1, kp);
                fn[nn] += fn[nf] * wn[nf];
                wn[nn] += wn[nf];
                nf = solution_index (&zn, s, ip, jp, kp+1);
                fn[nn] += fn[nf] * wn[nf];
                wn[nn] += wn[nf];
            }
            for (k = 1; k <= s->rind[2][1]; k++) {
                kp = zn.dim[2] + k - 1;
                nn = solution_index (&zn, s, ip, jp, kp);
                nf = solution_index (&zn, s, ip-1, jp, kp);
                fn[nn] += fn[nf] * wn[nf];
                wn[nn] += wn[nf];
                nf = solution_index (&zn, s, ip, jp-1, kp);
                fn[nn] += fn[nf] * wn[nf];
                wn[nn] += wn[nf];
                nf = solution_index (&zn, s, ip, jp, kp-1);
                fn[nn] += fn[nf] * wn[nf];
                wn[nn] += wn[nf];
            }
        }
    }

    for (nn = 0; nn < np; nn++) {
        if (wn[nn] == 0.0) wn[nn] = 1.0;
        f[nn] = fn[nn] / wn[nn];
    }
    free (fn);
}

/*-------------------------------------------------------------------*/

static void extract_subset (int nz, int subset[3])
{
    int inc[3];
    cgsize_t n, nn, np, i, j, k;
    ZONE *z = &Zones[nz-1];
    INTERFACE *zi;
    BOCO *zb;
    SOLUTION *zs;
    double *w;

    /* grid coordinates */

    for (n = 0; n < 3; n++)
        inc[n] = subset[n] < z->dim[n] ? subset[n] : 1;

    nn = 0;
    for (k = 1; k <= z->dim[2]; k += inc[2]) {
        for (j = 1; j <= z->dim[1]; j += inc[1]) {
            for (i = 1; i <= z->dim[0]; i += inc[0]) {
                np = vertex_index (z, i, j, k);
                z->verts[nn].x = z->verts[np].x;
                z->verts[nn].y = z->verts[np].y;
                z->verts[nn].z = z->verts[np].z;
                z->verts[nn].id = nn + 1;
                nn++;
            }
        }
    }

    /* interfaces */

    for (zi = z->ints, n = 0; n < z->nints; n++, zi++) {
        for (i = 0; i < 3; i++) {
            nn = np = inc[i];
            if (zi->d_zone)
                np = subset[i] < Zones[zi->d_zone-1].dim[i] ? subset[i] : 1;
            for (j = 0; j < 2; j++) {
                zi->range[i][j] = (zi->range[i][j] - 1) / nn + 1;
                zi->d_range[i][j] = (zi->d_range[i][j] - 1) / np + 1;
            }
        }
    }

    /* connectivities - remove these for now */

    z->nconns = 0;

    /* boundary conditions */

    for (zb = z->bocos, nn = 0; nn < z->nbocos; nn++, zb++) {
        if (zb->ptype == CGNS_ENUMV(PointRange)) {
            for (n = 0, j = 0; j < 2; j++) {
                for (i = 0; i < 3; i++, n++)
                    zb->pnts[n] = (zb->pnts[n] - 1) / inc[i] + 1;
            }
            /* normals ? */
        }
        else
            zb->id = 0;
    }

    /* solutions */

    np = (z->dim[0] - 1) * (z->dim[1] - 1) * (z->dim[2] - 1);
    nn = 0;
    for (zs = z->sols, n = 0; n < z->nsols; n++, zs++)
        nn += zs->nflds;
    if (!np || !nn) {
        z->nsols = 0;
        return;
    }

    w = (double *) malloc ((size_t)np * sizeof(double));
    if (NULL == w)
        FATAL (NULL, "malloc failed for weighting array");
    if (weighting) {
        n = 0;
        for (k = 1; k < z->dim[2]; k++)
            for (j = 1; j < z->dim[1]; j++)
                for (i = 1; i < z->dim[0]; i++)
                    w[n++] = cell_volume (z, i, j, k);
    }
    else {
        for (n = 0; n < np; n++)
            w[n] = 1.0;
    }

    for (zs = z->sols, n = 0; n < z->nsols; n++, zs++) {
        for (n = 0; n < zs->nflds; n++)
            field_subset (z, zs, inc, zs->flds[n].data, w);
    }
    free (w);
}

/*-------------------------------------------------------------------*/

int main (int argc, char *argv[])
{
    int n, nz, celldim, phydim;
    int subset[3], izone = 0;
    char basename[33], *newbase = NULL, *tmpfile;
    ZONE *z;

    if (argc < 2)
        print_usage (usgmsg, NULL);

    for (n = 0; n < 3; n++)
        subset[n] = 1;

    while ((n = getargs (argc, argv, options)) > 0) {
        switch (n) {
            case 'b':
                cgnsbase = atoi (argarg);
                break;
            case 'z':
                izone = atoi (argarg);
                break;
            case 'B':
                newbase = argarg;
                break;
            case 'a':
                for (n = 0; n < 3; n++)
                    subset[n] = 2;
                break;
            case 'i':
            case 'j':
            case 'k':
                subset[n-'i'] = 2;
                break;
            case 'w':
                weighting = 1;
                break;
        }
    }

    if (argind == argc)
        print_usage (usgmsg, "CGNSfile not given");
    if (!file_exists (argv[argind]))
        FATAL (NULL, "CGNSfile does not exist or is not a file");

    if (subset[0] + subset[1] + subset[2] == 3) {
        for (n = 0; n < 3; n++)
            subset[n] = 2;
    }

    /* create a working copy */

    printf ("creating a working copy of %s\n", argv[argind]);
    tmpfile = temporary_file (argv[argind]);
    copy_file (argv[argind], tmpfile);

    /* read CGNS file */

    printf ("reading CGNS file from %s\n", tmpfile);
    if (cg_open (tmpfile, CG_MODE_MODIFY, &cgnsfn) ||
        cg_base_read (cgnsfn, cgnsbase, basename, &celldim, &phydim))
        FATAL (NULL, NULL);
    if (celldim != 3 || phydim != 3)
        FATAL (NULL, "cell and/or physical dimension must be 3");
    printf ("reading zone information for base %d - %s\n",
        cgnsbase, basename);
    read_zones ();

    /* extract the subset */

    for (z = Zones, nz = 1; nz <= nZones; nz++, z++) {
        if (izone && nz != izone)
            z->id = 0;
        else {
            if (z->type != CGNS_ENUMV(Structured)) {
                sprintf (basename, "zone %d is not Structured", nz);
                FATAL (NULL, basename);
            }
            printf ("extracting subset for zone %d - %ldx%ldx%ld ... ", nz,
                (long)z->dim[0], (long)z->dim[1], (long)z->dim[2]);
            fflush (stdout);
            read_zone_data (nz);
            extract_subset (nz, subset);
            puts ("done");
        }
    }

    /* write CGNS file */

    if (newbase != NULL) {
        strncpy (basename, newbase, 32);
        basename[32] = 0;
        if (cg_base_write (cgnsfn, basename, 3, 3, &cgnsbase))
            FATAL (NULL, NULL);
        printf ("output to base %d - %s\n", cgnsbase, basename);
    }

    for (z = Zones, nz = 1; nz <= nZones; nz++, z++) {
        if (z->id) {
            z->nverts = 1;
            for (n = 0; n < 3; n++) {
                if (subset[n] < z->dim[n])
                    z->dim[n] = (z->dim[n] - 1) / subset[n] + 1;
                z->nverts *= z->dim[n];
            }
            printf ("writing zone %d - %ldx%ldx%ld ... ", nz,
                (long)z->dim[0], (long)z->dim[1], (long)z->dim[2]);
            fflush (stdout);
            write_zone_data (nz);
            puts ("done");
        }
    }

    cg_close (cgnsfn);

    if (argind + 1 < argc) argind++;
    printf ("renaming %s to %s\n", tmpfile, argv[argind]);
    unlink (argv[argind]);
    if (rename (tmpfile, argv[argind])) {
        char msg[512];
        sprintf (msg, "rename %s -> %s failed", tmpfile, argv[argind]);
        FATAL (NULL, msg);
        exit (1);
    }
    free (tmpfile);
    return 0;
}
