#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cgnames.h"

static int print_identifier (char *name, int nexps, float *exps, void *user)
{
    int n, i;
    char term[16], num[33], den[33];
    static char label[] = "MLT@a";

    if (nexps <= 0) {
        puts (name);
        return 0;
    }
    printf ("%-32s", name);
    num[0] = den[0] = 0;
    for (n = 0; n < nexps; n++) {
        i = (int)exps[n];
        if (i) {
            if (abs(i) > 1)
                sprintf (term, "%c^%d", label[n], abs(i));
            else
                sprintf (term, "%c", label[n]);
            if (i < 0)
                strcat (den, term);
            else
                strcat (num, term);
        }
    }
    if (!num[0]) strcpy (num, "1");
    printf ("%s", num);
    if (den[0]) printf (" / %s", den);
    putchar ('\n');
    return 0;
}

int main ()
{
    cg_enum_identifier (print_identifier, NULL);
    return 0;
}
