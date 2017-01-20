/*
 * allocator.cpp
 *
 *  Created on: Oct 6, 2015
 *      Author: cm654063
 */

#include "malloc.h"

void* operator new(size_t size)
{
	return malloc(size);
}
void operator delete(void* p)
{
	free(p);
}
