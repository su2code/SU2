


#ifndef __ZZZ_HEADER_H__
#define __ZZZ_HEADER_H__

#include "libm3l.h"
#include "lsipdx.h"


typedef struct d6dof{
	double angles[3];
	double transvec[3], rotcenter[3];
}d6dof_t;

typedef struct conn{
	char name_of_channel_i[MAX_NAME_LENGTH], name_of_channel_o[MAX_NAME_LENGTH];
	char hostname[80];
	int portno;
}conn_t;

#endif
