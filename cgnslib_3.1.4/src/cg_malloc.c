#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static size_t mem_now = 0;
static size_t mem_max = 0;
static size_t alloc_calls = 0;
static size_t free_calls = 0;

void *cgmalloc(size_t bytes) {
    size_t *data = (size_t *) malloc (bytes + 2 * sizeof(size_t));
    if (data) {
        *data++ = bytes;
        *data++ = (size_t)(data + 1);
        mem_now += bytes;
        if (mem_max < mem_now) mem_max = mem_now;
    }
    alloc_calls++;
    return (void *)data;
}

void *cgrealloc(void *olddata, size_t bytes) {
    size_t *data = (size_t *)olddata;
    if (data && *(data-1) == (size_t)olddata) {
        data -= 2;
        mem_now -= *data;
    }
    data = (size_t *) realloc (data, bytes + 2 * sizeof(size_t));
    if (data) {
        *data++ = bytes;
        *data++ = (size_t)(data + 1);
        mem_now += bytes;
        if (mem_max < mem_now) mem_max = mem_now;
    }
/*    alloc_calls++;*/
    return (void *)data;
}

void *cgcalloc(size_t num, size_t bytes) {
    size_t count = num * bytes;
    size_t *data = (size_t *) malloc (count + 2 * sizeof(size_t));
    if (data) {
        *data++ = count;
        *data++ = (size_t)(data + 1);
        mem_now += count;
        if (mem_max < mem_now) mem_max = mem_now;
        memset(data, 0, count);
    }
    alloc_calls++;
    return (void *)data;
}

void cgfree(void *freedata) {
    size_t *data = (size_t *)freedata;
    if (data && *(data-1) == (size_t)freedata) {
        data -= 2;
        mem_now -= *data;
        free (data);
    }
    free_calls++;
}

size_t cgmemnow() {return mem_now;}
size_t cgmemmax() {return mem_max;}

size_t cgalloccalls() {return alloc_calls;}
size_t cgfreecalls()  {return free_calls;}

