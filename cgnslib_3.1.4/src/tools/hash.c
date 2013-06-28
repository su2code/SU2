/*
 * hash.c - general purpose hash table functions
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hash.h"

/*#define STATISTICS    define this for hit/seeks counts*/
/*#define REORDER       define this for reordering after find*/

typedef struct hash_element_ {
    void *entry;
    struct hash_element_ *next;
    struct hash_element_ **prev;
#ifdef STATISTICS
    size_t hits;
    size_t seeks;
#endif
} BUCKET;

typedef struct hash_tab_ {
    BUCKET **table;     /* pointer to hash table */
    size_t size;        /* Max number of elements in table */
    size_t numsyms;     /* number of elements currently in table */
    BUCKET *lastpos;    /* last bucket accessed by find */
    int (*compare)(     /* entry compare function */
        void *, void *
    );
    size_t (*hash)(     /* pointer to hashing routine */
        void *
    );
} HASH_TAB;

/*---------- HashCreate ----------------------------------------------
 * creates a hash table of indicated size.
 * Should be a prime number for best results.
 * Useful values are:
 *   47 61 89 113 127 157 193 211 257 293 359 401
 *--------------------------------------------------------------------*/

HASH HashCreate (size_t size, int (*compare)(void *entry1, void *entry2),
                 size_t (*hashfunc)(void *entry))
{
    HASH_TAB *tabp;

    if (NULL == compare || NULL == hashfunc)
        return (NULL);
    if (size < 1) size = 11;
    tabp = (HASH_TAB *) malloc (sizeof(HASH_TAB) + size * sizeof(BUCKET *));
    if (NULL != tabp) {
        tabp->table   = (BUCKET **)(tabp + 1);
        tabp->size    = size;
        tabp->numsyms = 0;
        tabp->lastpos = NULL;
        tabp->compare = compare;
        tabp->hash    = hashfunc;
        for (size = 0; size < tabp->size; size++)
            tabp->table[size] = NULL;
    }
    return ((HASH)tabp);
}

/*---------- HashDestroy ---------------------------------------------
 * destroy a hash table
 *--------------------------------------------------------------------*/

void HashDestroy (HASH hash, void (*freeentry)(void *entry))
{
    HASH_TAB *tabp = (HASH_TAB *)hash;
    BUCKET *p, *next;
    size_t i;

    for (i = 0; i < tabp->size; i++) {
        for (p = tabp->table[i]; p != NULL; p = next) {
            next = p->next;
            if (NULL != freeentry)
                (*freeentry) (p->entry);
            free (p);
        }
    }
    free (tabp);
}

/*---------- HashFind -------------------------------------------------
 * finds symbol in hash table or NULL if no match
 * If more than one such entry is in the table, the most
 * recently added one is found (which will be first in the list).
 *---------------------------------------------------------------------*/

void *HashFind (HASH hash, void *entry)
{
    BUCKET *p, **pf;
    HASH_TAB *tabp = (HASH_TAB *)hash;

    pf = &(tabp->table)[(*tabp->hash) (entry) % tabp->size];
    p = *pf;

    while (NULL != p && (*tabp->compare) (entry, p->entry)) {
#ifdef STATISTICS
        p->seeks++;
#endif
        p = p->next;
    }
    if (NULL == p)
        return (NULL);

#ifdef STATISTICS
    p->hits++;
#endif
    tabp->lastpos = p;

    /*
     * this mod places the symbol at the beginning of the list
     * if not already there. This tends to cluster the most often
     * accessed symbols at the beginning of the table, thus
     * reducing the length of the search path
     */

#ifdef REORDER
    if (p != *pf) {
        BUCKET *sym = *pf;
        if ((*(p->prev) = p->next) != NULL)
            p->next->prev = p->prev;
        *pf = p;
        p->prev = pf;
        p->next = sym;
        sym->prev = &p->next;
    }
#endif

    return (p->entry);
}

/*---------- HashAdd -------------------------------------------------
 * add a symbol to the hash table
 * The new symbol is placed first in the bucket.
 *--------------------------------------------------------------------*/

void *HashAdd (HASH hash, void *entry)
{
    BUCKET **p, *tmp;
    BUCKET *sym;
    HASH_TAB *tabp = (HASH_TAB *)hash;

    if ((sym = (BUCKET *) malloc (sizeof(BUCKET))) == NULL)
        return (NULL);

    p = &(tabp->table)[(*tabp->hash) (entry) % tabp->size];

    tmp = *p;
    *p = sym;
    sym->entry = entry;
    sym->prev = p;
    sym->next = tmp;

    if (NULL != tmp)
        tmp->prev = &sym->next;

    tabp->numsyms++;
    tabp->lastpos = NULL;
    return (entry);
}

/*---------- HashDelete ----------------------------------------------
 * removes a symbol from the hash table
 *--------------------------------------------------------------------*/

void *HashDelete (HASH hash, void *entry)
{
    BUCKET *p;
    HASH_TAB *tabp = (HASH_TAB *)hash;

    if (NULL == tabp->lastpos || entry != tabp->lastpos->entry) {
        if (NULL == HashFind (hash, entry))
            return (NULL);
    }
    p = tabp->lastpos;
    tabp->numsyms--;
    if ((*(p->prev) = p->next) != NULL)
        p->next->prev = p->prev;
    free (p);
    return (entry);
}

/*---------- HashList ------------------------------------------------
 * output the hash table
 *--------------------------------------------------------------------*/

size_t HashList (HASH hash, size_t (*listentry)(void *entry, void *userdata),
                 void *userdata)
{
    HASH_TAB *tabp = (HASH_TAB *)hash;
    BUCKET *p;
    size_t i, cnt = 0;

    if (NULL == listentry)
        return (tabp->numsyms);
    for (i = 0; i < tabp->size; i++) {
        for (p = tabp->table[i]; p != NULL; p = p->next)
            cnt += (*listentry) (p->entry, userdata);
    }
    return cnt;
}

/*---------- HashStats ------------------------------------------------
 * print a variety of statistics about the hash table
 *---------------------------------------------------------------------*/

#define MAXINT  ((size_t) ~ 0)
#define MAXLEN  128

void HastStats (HASH hash)
{
    HASH_TAB *tabp = (HASH_TAB *)hash;
    BUCKET *p;
    size_t i;
    size_t chain_len;
    size_t chain_avg;
    size_t maxlen = 0;
    size_t minlen = MAXINT;
    size_t lengths[MAXLEN];
    size_t longer = 0;
#ifdef STATISTICS
    size_t hits = 0;
    size_t seeks = 0;
#endif
    int deviation = 0;

    chain_avg = tabp->numsyms / tabp->size;
    memset (lengths, 0, sizeof(lengths));

    for (i = 0; i < tabp->size; i++) {
        chain_len = 0;
        for (p = tabp->table[i]; p; p = p->next) {
#ifdef STATISTICS
            hits += p->hits;
            seeks += p->seeks;
#endif
            chain_len++;
        }

        if (chain_len >= MAXLEN)
            longer++;
        else
            ++lengths[chain_len];

        if (minlen > chain_len) minlen = chain_len;
        if (maxlen < chain_len) maxlen = chain_len;

        deviation += abs ((int)(chain_len - chain_avg));
    }

    printf ("%ld entries in %ld element hash table, ",
        (long)tabp->numsyms, (long)tabp->size);
    printf ("%ld (%1.0f%%) empty.\n",
        (long)lengths[0], ((double)lengths[0] / tabp->size) * 100.0);
    printf ("Mean chain length = %ld, min = %ld, max = %ld, dev = %ld\n",
        (long)chain_avg, (long)minlen, (long)maxlen,
        (long)(deviation / tabp->size));
#ifdef STATISTICS
    printf ("Table searchs = %ld, found = %ld, percent = %.2f\n",
        (long)seeks, (long)hits,
        (seeks ? (double)hits / seeks * 100.0 : (double)0));
#endif

    for (i = 0; i < MAXLEN; i++)
        if (lengths[i])
            printf ("%ld chains of length %d\n", (long)lengths[i], (int)i);

    if (longer)
        printf ("%ld chains of length %d or longer\n", (long)longer, MAXLEN);
}

