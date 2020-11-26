// (c)2014 Bram Stolk, Vancouver B.C.

// Print out SIMD 8xfloat value.
#define PRT_m256( V ) \
{ \
	float vals[8] __attribute__ ((aligned (32))); \
	 _mm256_store_ps( vals, V ); \
	fprintf \
	( \
		stderr,  \
		"%s: %f/%f/%f/%f/%f/%f/%f/%f\n", \
		#V, \
		vals[0], vals[1], vals[2], vals[3], vals[4], vals[5], vals[6], vals[7] \
	); \
}

// Print out SIMD 8xint value.
#define PRT_m256i( V ) \
{ \
	unsigned int vals[8] __attribute__ ((aligned (32))); \
	 _mm256_store_si256( (__m256i*) vals, V ); \
	fprintf \
	( \
		stderr,  \
		"%s: %08x/%08x/%08x/%08x/%08x/%08x/%08x/%08x\n", \
		#V, \
		vals[0], vals[1], vals[2], vals[3], vals[4], vals[5], vals[6], vals[7] \
	); \
}



// Print out SIMD 16xfloat value.
#define PRT_m512( V ) \
{ \
	float vals[16] __attribute__ ((aligned (64))); \
	 _mm512_store_ps( vals, V ); \
	fprintf \
	( \
		stderr,  \
		"%s: %f/%f/%f/%f/%f/%f/%f/%f %f/%f/%f/%f/%f/%f/%f/%f\n", \
		#V, \
		vals[ 0], vals[ 1], vals[ 2], vals[ 3], vals[ 4], vals[ 5], vals[ 6], vals[ 7], \
		vals[ 8], vals[ 9], vals[10], vals[11], vals[12], vals[13], vals[14], vals[15] \
	); \
}


// Print out SIMD 16xint value.
#define PRT_m512i( V ) \
{ \
	int vals[16] __attribute__ ((aligned (64))); \
	 _mm512_store_epi32( vals, V ); \
	fprintf \
	( \
		stderr,  \
		"%s: %d/%d/%d/%d/%d/%d/%d/%d %d/%d/%d/%d/%d/%d/%d/%d\n", \
		#V, \
		vals[ 0], vals[ 1], vals[ 2], vals[ 3], vals[ 4], vals[ 5], vals[ 6], vals[ 7], \
		vals[ 8], vals[ 9], vals[10], vals[11], vals[12], vals[13], vals[14], vals[15] \
	); \
}

