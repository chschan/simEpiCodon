#define MATHLIB_STANDALONE 1
#include <Rmath.h>

#include <time.h>
#include <limits.h>

void seed_time(void)
{
	time_t timeval1, timeval2;
	unsigned char *ptr;
	unsigned int seed1, seed2;
	size_t i;

	timeval1 = time(NULL);
	ptr = (unsigned char *) &timeval1;
	for (i = 0; i < sizeof(timeval1); i++)
		seed1 = seed1 * (UCHAR_MAX + 2u) + ptr[i];

	timeval2 = time(NULL);
	ptr = (unsigned char *) &timeval2;
	for (i = 0; i < sizeof(timeval2); i++)
		seed2 = seed2 * (UCHAR_MAX + 2u) + ptr[i];

	set_seed(seed1, seed2);
}

