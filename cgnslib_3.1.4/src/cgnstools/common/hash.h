/*
 * hash.h - include file for hash routines
 */

#ifndef _HASH_H_
#define _HASH_H_

#include <stdio.h>  /* for NULL */

typedef void *HASH;

#ifdef __cplusplus
extern "C" {
#endif

HASH   HashCreate  (size_t size, int (*compare) (void *entry1, void *entry2),
                    size_t (*hashfunc) (void *entry));
void   HashDestroy (HASH hash, void (*freeentry) (void *entry));
void  *HashFind    (HASH hash, void *entry);
void  *HashAdd     (HASH hash, void *entry);
void  *HashDelete  (HASH hash, void *entry);
size_t HashList    (HASH hash, size_t (*listentry)(void *entry, void *userdata),
                    void *userdata);
void  HashStats    (HASH hash);

#define HashSize(H) HashList(H,NULL,NULL)

#ifdef __cplusplus
}
#endif

#endif  /* _HASH_H_ */

