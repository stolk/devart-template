// Random Space Filling Tiling.
// As described by: http://paulbourke.net/texture_colour/randomtile/
//
// (c)2014 Bram Stolk, Vancouver B.C.

#include <immintrin.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include "hsv.h"	// for colour space conversions.
#include "prt.h"	// for printing __m256 and __m256i values to stderr.


// Scalar left/right test.
// returns >0 for p2 left of line p0-p1
// returns =0 for p2 on the line p0-p1
// returns <0 for p2 right of line p0-p1
float onLeft( float p0x, float p0y, float p1x, float p1y, float p2x, float p2y )
{
	return ( ( p1x - p0x ) * ( p2y - p0y ) - ( p2x - p0x ) * ( p1y - p0y ) );
}


// 8-way SIMD version of testing whether p2 lies to left or right on p0-p1 line.
__m256 onLeft8( __m256 p0x, __m256 p0y, __m256 p1x, __m256 p1y, __m256 p2x, __m256 p2y )
{
	const __m256 t0 = _mm256_sub_ps( p1x, p0x );
	const __m256 t1 = _mm256_sub_ps( p2y, p0y );
	const __m256 t2 = _mm256_sub_ps( p2x, p0x );
	const __m256 t3 = _mm256_sub_ps( p1y, p0y );
	const __m256 f0 = _mm256_mul_ps( t0, t1 );
	const __m256 f1 = _mm256_mul_ps( t2, t3 );
	return _mm256_sub_ps( f0, f1 );
}


// 8-way SIMD calculation of winding number of a point versus a polygon size 8.
int winding_number_8( float px, float py, __m256 vx, __m256 vy, __m256 wx, __m256 wy )
{
	const __m256  zero8  = _mm256_setzero_ps();
	const __m256i zero8i = _mm256_setzero_si256();
	const __m256i one8i = _mm256_set1_epi32( 1 );
	const __m256i minone8i = _mm256_set1_epi32( -1 );

	const __m256 px8 = _mm256_set1_ps( px );
	const __m256 py8 = _mm256_set1_ps( py );

	const __m256 ol = onLeft8( vx, vy, wx, wy, px8, py8 );
	const __m256i olltz = _mm256_cmp_ps( ol, zero8, _CMP_LT_OS );
	const __m256i olgtz = _mm256_cmp_ps( ol, zero8, _CMP_GT_OS );

	const __m256i pred0 = _mm256_cmp_ps( wy, py8, _CMP_GT_OS );
	const __m256i pred1 = _mm256_cmp_ps( wy, py8, _CMP_LE_OS );
	const __m256i useup = _mm256_cmp_ps( vy, py8, _CMP_LE_OS );

	const __m256i up_intersect = _mm256_blendv_ps( zero8i,    one8i, _mm256_and_ps( pred0, olgtz ) );
	const __m256i dn_intersect = _mm256_blendv_ps( zero8i, minone8i, _mm256_and_ps( pred1, olltz ) );

	__m256i addition = _mm256_blendv_ps( dn_intersect, up_intersect, useup );

	// There should be a faster way to sum over the 8 components of a vector value.
	int additionvals[ 8 ] __attribute__ ((aligned (32)));
	 _mm256_store_si256( (__m256i*) additionvals, addition );
	return additionvals[0] + additionvals[1] + additionvals[2] + additionvals[3] +
	       additionvals[4] + additionvals[5] + additionvals[6] + additionvals[7];
}


// Scalar version for calculating winding number of a point versus a polygon size 8.
// after: http://geomalgorithms.com/a03-_inclusion.html
int winding_number( float px, float py, float *vx, float* vy, int n )
{
	int wn=0;
	for ( int i=0; i<n; ++i )
	{
		int j = (i+1)%n;
		const float ol = onLeft( vx[i], vy[i], vx[j], vy[j], px, py );
		int pred0 = (vy[j] >  py );
		int pred1 = (vy[j] <= py );
		printf("vy[j]=%f py=%f\n", vy[j], py );
		printf("pred0:%d pred1:%d\n", pred0, pred1);
		const int up_intersect = ( vy[j] >  py && ol > 0.0f ) ?  1 : 0;
		const int dn_intersect = ( vy[j] <= py && ol < 0.0f ) ? -1 : 0;
		const int addition = ( vy[i] <= py ) ? up_intersect : dn_intersect;
		printf( "up:%d dn:%d addition:%d\n", up_intersect, dn_intersect, addition );
		wn += addition;
	}
	return wn;
}


// Scalar version of edge vs edge intersection.
int edge_vs_edge(float x1,float y1,float x2,float y2,float x3,float y3,float x4,float y4)
{
	float mua,mub;
	float denom,numera,numerb;

	denom  = (y4-y3) * (x2-x1) - (x4-x3) * (y2-y1);
	numera = (x4-x3) * (y1-y3) - (y4-y3) * (x1-x3);
	numerb = (x2-x1) * (y1-y3) - (y2-y1) * (x1-x3);

	// Is the intersection along the the segments
	mua = numera / denom;
	mub = numerb / denom;

	if (mua < 0 || mua > 1 || mub < 0 || mub > 1)
		return 0;
	return 1;
}


// 8-way SIMD version of edge veruses eight edges intersection.
int edge_vs_edge8( float x1, float y1, float x2, float y2, __m256 vx3, __m256 vy3, __m256 vx4, __m256 vy4 )
{
	const __m256 vx1 = _mm256_set1_ps( x1 );
	const __m256 vy1 = _mm256_set1_ps( y1 );
	const __m256 vx2 = _mm256_set1_ps( x2 );
	const __m256 vy2 = _mm256_set1_ps( y2 );

	// calc all inner terms
	const __m256 y4_y3 = _mm256_sub_ps( vy4, vy3 );
	const __m256 x2_x1 = _mm256_sub_ps( vx2, vx1 );
	const __m256 x4_x3 = _mm256_sub_ps( vx4, vx3 );
	const __m256 y2_y1 = _mm256_sub_ps( vy2, vy1 );
	const __m256 y1_y3 = _mm256_sub_ps( vy1, vy3 );
	const __m256 x1_x3 = _mm256_sub_ps( vx1, vx3 );

	// calc all factors
	const __m256 f0 = _mm256_mul_ps( y4_y3, x2_x1 );
	const __m256 f1 = _mm256_mul_ps( x4_x3, y2_y1 );
	const __m256 f2 = _mm256_mul_ps( x4_x3, y1_y3 );
	const __m256 f3 = _mm256_mul_ps( y4_y3, x1_x3 );
	const __m256 f4 = _mm256_mul_ps( x2_x1, y1_y3 );
	const __m256 f5 = _mm256_mul_ps( y2_y1, x1_x3 );

	// calc outer terms
	const __m256 denomi = _mm256_rcp_ps( _mm256_sub_ps( f0, f1 ) );
	const __m256 numera = _mm256_sub_ps( f2, f3 );
	const __m256 numerb = _mm256_sub_ps( f4, f5 );

	const __m256 mua = _mm256_mul_ps( numera, denomi );
	const __m256 mub = _mm256_mul_ps( numerb, denomi );

	const __m256i mua_lo = _mm256_cmp_ps( mua, _mm256_setzero_ps(), _CMP_LT_OS );
	const __m256i mub_lo = _mm256_cmp_ps( mub, _mm256_setzero_ps(), _CMP_LT_OS );
	const __m256i mua_hi = _mm256_cmp_ps( mua, _mm256_set1_ps( 1.0f ), _CMP_GT_OS );
	const __m256i mub_hi = _mm256_cmp_ps( mub, _mm256_set1_ps( 1.0f ), _CMP_GT_OS );

	const __m256i mua_hilo = _mm256_or_ps( mua_lo, mua_hi );
	const __m256i mub_hilo = _mm256_or_ps( mub_lo, mub_hi );

        __m256i hilo = _mm256_or_ps( mua_hilo, mub_hilo );
	hilo = _mm256_xor_ps( hilo, _mm256_set1_epi32(-1) );

	const int t = _mm256_testz_si256( hilo, hilo );
	return !t;
}


// 8-Way SIMD test for intersection between eight edges of one polygon with eight edges of other polygon.
int edges_intersect
(
	__m256 vxa,
	__m256 vya,
	__m256 nxa,
	__m256 nya,

	__m256 vxb,
	__m256 vyb,
	__m256 nxb,
	__m256 nyb
)
{
	float vxvals[8] __attribute__ ((aligned (32)));
	float vyvals[8] __attribute__ ((aligned (32)));
	float nxvals[8] __attribute__ ((aligned (32)));
	float nyvals[8] __attribute__ ((aligned (32)));

 	_mm256_store_si256( (__m256i*) vxvals, vxa );
 	_mm256_store_si256( (__m256i*) vyvals, vya );
 	_mm256_store_si256( (__m256i*) nxvals, nxa );
 	_mm256_store_si256( (__m256i*) nyvals, nya );

	for ( int i=0; i<8; ++i )
	{
		const int rv = edge_vs_edge8( vxvals[i], vyvals[i], nxvals[i], nyvals[i], vxb, vyb, nxb, nyb );
		//printf( "edge nr %d : intersection result %d\n", i, rv );
		if ( rv ) return 1;
	}
	return 0;
}


static float octx[8] __attribute__ ((aligned (32)));
static float octy[8] __attribute__ ((aligned (32)));
static float octX[8] __attribute__ ((aligned (32)));
static float octY[8] __attribute__ ((aligned (32)));
 
static float strx[8] __attribute__ ((aligned (32)));
static float stry[8] __attribute__ ((aligned (32)));
static float strX[8] __attribute__ ((aligned (32)));
static float strY[8] __attribute__ ((aligned (32)));

#define SHPCNT 6
static float shpx[SHPCNT][8] __attribute__ ((aligned(32)));
static float shpy[SHPCNT][8] __attribute__ ((aligned(32)));
static float shpX[SHPCNT][8] __attribute__ ((aligned(32)));
static float shpY[SHPCNT][8] __attribute__ ((aligned(32)));
static __m256 shpx8[SHPCNT];
static __m256 shpy8[SHPCNT];
static __m256 shpX8[SHPCNT];
static __m256 shpY8[SHPCNT];
static float  shparea[SHPCNT];

#define MAXSZ 1600
// db for large polies
static __m256 ldbx[MAXSZ];
static __m256 ldby[MAXSZ];
static __m256 ldbX[MAXSZ];
static __m256 ldbY[MAXSZ];
static int ldbsz=0;
// db for small polies (which we bin in grid buckets for fast spatial lookup.)
static __m256 sdbx[MAXSZ];
static __m256 sdby[MAXSZ];
static __m256 sdbX[MAXSZ];
static __m256 sdbY[MAXSZ];
static int sdbsz=0;

//Spatial indexing
#define GRIDRES 8
#define MAXINCELL MAXSZ
#define OVERLAP 0.4
#define CRITICALSCALE (OVERLAP/GRIDRES)
static int buckets[GRIDRES][GRIDRES][MAXINCELL];
static int bucketsizes[GRIDRES][GRIDRES];


static void add_to_bucket( int x, int y, int idx )
{
	if ( x>=0 && x<GRIDRES && y>=0 && y<GRIDRES )
		buckets[ x ][ y ][ bucketsizes[ x ][ y ]++ ] = idx;
}


// Bin a small polygon wih index 'idx' to one and possibly 1/2/3 neighbours.
static void add_to_buckets( float xo, float yo, int idx )
{
	const int gx = floorf( GRIDRES * xo );
	const int gy = floorf( GRIDRES * yo );
	const float dx = ( GRIDRES * xo ) - gx;
	const float dy = ( GRIDRES * yo ) - gy;
	add_to_bucket( gx, gy, idx );

	const int addl = dx < OVERLAP;
	const int addr = dx > (1.0f - OVERLAP);
	const int addd = dy < OVERLAP;
	const int addu = dy > (1.0f - OVERLAP);

	if ( addl && addd )
	{
		add_to_bucket( gx-1, gy, idx );
		add_to_bucket( gx, gy-1, idx );
		add_to_bucket( gx-1, gy-1, idx );
		return;
	}
	if ( addl && addu )
	{
		add_to_bucket( gx-1, gy, idx );
		add_to_bucket( gx, gy+1, idx );
		add_to_bucket( gx-1, gy+1, idx );
		return;
	}
	if ( addr && addd )
	{
		add_to_bucket( gx+1, gy, idx );
		add_to_bucket( gx, gy-1, idx );
		add_to_bucket( gx+1, gy-1, idx );
		return;
	}
	if ( addr && addu )
	{
		add_to_bucket( gx+1, gy, idx );
		add_to_bucket( gx, gy+1, idx );
		add_to_bucket( gx+1, gy+1, idx );
		return;
	}
	if ( addl )
	{
		add_to_bucket( gx-1, gy, idx );
		return;
	}
	if ( addr )
	{
		add_to_bucket( gx+1, gy, idx );
		return;
	}
	if ( addd )
	{
		add_to_bucket( gx, gy-1, idx );
		return;
	}
	if ( addu )
	{
		add_to_bucket( gx, gy+1, idx );
		return;
	}
}


static int valid_placement( __m256 x8, __m256 y8, __m256 X8, __m256 Y8, float xo, float yo )
{
	static float x[8] __attribute__ ((aligned (32)));
	static float y[8] __attribute__ ((aligned (32)));
	_mm256_store_ps( x, x8 );
	_mm256_store_ps( y, y8 );

	static float X[8] __attribute__ ((aligned (32)));
	static float Y[8] __attribute__ ((aligned (32)));
	_mm256_store_ps( X, X8 );
	_mm256_store_ps( Y, Y8 );

	// consider all currently placed large polygons.
	for ( int i=0; i<ldbsz; ++i )
	{
		// any of the points inside a placed polygon? If so, invalid placement.
		for ( int j=0; j<8; ++j )
		{
			const int pt_inside = winding_number_8( x[j], y[j], ldbx[i], ldby[i], ldbX[i], ldbY[i] );
			if ( pt_inside ) return 0;
		}
		// any of the edges crosses a placed polygon? If so, invalid placement.
		int ei0 = edges_intersect( x8, y8, X8, Y8, ldbx[i], ldby[i], ldbX[i], ldbY[i] );
		if ( ei0 ) return 0;
	}

#if 0
	// consider all currently placed small polygons.
	for ( int i=0; i<sdbsz; ++i )
	{
		// any of the points inside a placed polygon? If so, invalid placement.
		for ( int j=0; j<8; ++j )
		{
			const int pt_inside = winding_number_8( x[j], y[j], sdbx[i], sdby[i], sdbX[i], sdbY[i] );
			if ( pt_inside ) return 0;
		}
		// any of the edges crosses a placed polygon? If so, invalid placement.
		int ei0 = edges_intersect( x8, y8, X8, Y8, sdbx[i], sdby[i], sdbX[i], sdbY[i] );
		if ( ei0 ) return 0;
	}
	return 1;
#endif

	// consider all small polygons that are in the same grid cell
	int gx = floorf( GRIDRES * xo );
	int gy = floorf( GRIDRES * yo );
	gx = gx >= GRIDRES ? GRIDRES-1 : gx;
	gy = gy >= GRIDRES ? GRIDRES-1 : gy;
	assert( gx >= 0 );
	assert( gy >= 0 );
	const int cnt = bucketsizes[gx][gy];
	for ( int k=0; k<cnt; ++k )
	{
		const int i = buckets[gx][gy][k];
		assert( i < MAXSZ );
		// any of the points inside a placed polygon? If so, invalid placement.
		for ( int j=0; j<8; ++j )
		{
			const int pt_inside = winding_number_8( x[j], y[j], sdbx[i], sdby[i], sdbX[i], sdbY[i] );
			if ( pt_inside ) return 0;
		}
		// any of the edges crosses a placed polygon? If so, invalid placement.
		int ei0 = edges_intersect( x8, y8, X8, Y8, sdbx[i], sdby[i], sdbX[i], sdbY[i] );
		if ( ei0 ) return 0;
	}
	return 1;
}


float areasize( __m256 x8, __m256 y8, __m256 X8, __m256 Y8 )
{
	int hit=0;
	int mis=0;
	for ( int y=-32; y<=32; ++y )
	{
		for ( int x=-32; x<=32; ++x )
		{
			const int rv = winding_number_8( x/32.0f, y/32.0f, x8, y8, X8, Y8 );
			hit += rv;
			mis += !rv;
		}
	}
	return  4 * hit / (float)(hit+mis);
}


void rotate_shape( float angle, __m256* x, __m256* y, __m256* X, __m256* Y )
{
	const __m256 co = _mm256_set1_ps( cosf( angle ) );
	const __m256 si = _mm256_set1_ps( sinf( angle ) );
	const __m256 newx = _mm256_sub_ps( _mm256_mul_ps( co, *x ), _mm256_mul_ps( si, *y ) );
	const __m256 newy = _mm256_add_ps( _mm256_mul_ps( si, *x ), _mm256_mul_ps( co, *y ) );
	const __m256 newX = _mm256_sub_ps( _mm256_mul_ps( co, *X ), _mm256_mul_ps( si, *Y ) );
	const __m256 newY = _mm256_add_ps( _mm256_mul_ps( si, *X ), _mm256_mul_ps( co, *Y ) );
	*x = newx;
	*y = newy;
	*X = newX;
	*Y = newY;
}


int main( int argc, char* argv[] )
{
	// init shapes
	for ( int i=0; i<8; ++i )
	{
		int j = (i+1)%8;
		const float a0 = i * M_PI / 4;
		const float a1 = j * M_PI / 4;
		float r = 0.5f;
		shpx[4][ i ] = r * cosf( a0 );
		shpy[4][ i ] = r * sinf( a0 );
		float r0 = (i&1) ? 0.3f : 0.5f;
		float r1 = (j&1) ? 0.3f : 0.5f;
		shpx[5][ i ] = r0 * cosf( a0 );
		shpy[5][ i ] = r0 * sinf( a0 );
	}
	const float sqrx[8] = { -0.5, -0.5, -0.5, 0, 0.5, 0.5, 0.5, 0 };
	const float sqry[8] = { 0.5, 0, -0.5, -0.5, -0.5, 0, 0.5, 0.5 };
	memcpy( shpx[0], sqrx, sizeof(sqrx) );
	memcpy( shpy[0], sqry, sizeof(sqry) );
	const float arrx[8] = { -0, 0.4, 0.2, 0.2,     0, -0.2, -0.2, -0.4 };
	const float arry[8] = { -0.6, -0.2, -0.2, 0.6,    0.4, 0.6, -0.2, -0.2 };
	memcpy( shpx[1], arrx, sizeof(arrx) );
	memcpy( shpy[1], arry, sizeof(arry) );
	const float utnx[8] = { -0.1, -0.2, -0.2, 0.2,    0.2, 0.1, 0.1, -0.1 };
	const float utny[8] = { 0.3, 0.3, -0.3, -0.3,     0.3, 0.3, -0.2, -0.2 };
	memcpy( shpx[2], utnx, sizeof(utnx) );
	memcpy( shpy[2], utny, sizeof(utny) );
	const float jtnx[8] = { -0.1, -0.3, -0.3, -0.1,   0.3, 0.3, 0.1, 0.1 };
	const float jtny[8] = { -0.1, -0.1, -0.3, -0.3,   0.1, 0.3, 0.3, 0.1 };
	memcpy( shpx[3], jtnx, sizeof(jtnx) );
	memcpy( shpy[3], jtny, sizeof(jtny) );

	for ( int s=0; s<SHPCNT; ++s )
	{
		for ( int i=0; i<8; ++i )
		{
			int j = (i+1)%8;
			shpX[s][i] = shpx[s][j];
			shpY[s][i] = shpy[s][j];
		}

		shpx8[s] = _mm256_load_ps( shpx[s] );
		shpy8[s] = _mm256_load_ps( shpy[s] );
		shpX8[s] = _mm256_load_ps( shpX[s] );
		shpY8[s] = _mm256_load_ps( shpY[s] );
		shparea[s] = areasize( shpx8[s], shpy8[s], shpX8[s], shpY8[s] );
		fprintf( stderr, "area of shape %d: %f\n", s, shparea[s] );
	}

	// init grid
	for ( int xx=0; xx<GRIDRES; ++xx )
		for ( int yy=0; yy<GRIDRES; ++yy )
			bucketsizes[xx][yy]=0;

#if 0
	for ( int y=-20; y<=20; ++y )
	{
		for ( int x=-20; x<=20; ++x )
		{
			//const int rv1 = winding_number( x/36.0f, y/36.0f, strx, stry, 8 );
			const int rv = winding_number_8( x/36.0f, y/36.0f, strx8, stry8, strX8, strY8 );
			fprintf( stdout, "%s", (rv==0) ? "--":"[]" );
			//exit(1);
		}
		fprintf( stdout, "\n" );
	}
#endif

#if 0
	int rv = edge_vs_edge8( 0, 0, 0.3, 0.2, octx8, octy8, octX8, octY8 );
	fprintf( stdout, "rv=%d\n", rv );
	rv = edge_vs_edge8( 0, 0.1, 3, 0.1, octx8, octy8, octX8, octY8 );
	fprintf( stdout, "rv=%d\n", rv );
	int ei = edges_intersect( octx8, octy8, octX8, octY8, strx8, stry8, strX8, strY8 );
	printf( "ei = %d\n", ei );
#endif


	fprintf( stdout, "<?xml version=\"1.0\"?>\n" );
	fprintf( stdout, "<svg version=\"1.1\" baseProfile=\"tiny\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" id=\"svg-root\" viewBox=\"0 0 1 1\">\n" );

	const float c = 1.30;
	for ( int i=0; i<MAXSZ; ++i )
	{
		//const int shpnr = 3;
		const int shpnr = i%SHPCNT;
		const float scl = sqrtf( powf( 5+i, -c ) / (2*shparea[shpnr]) );
		fprintf( stderr, "scl %f\n", scl );
		int valid = 0;
		int trials = 0;
		do 
		{
			const float xo = (rand()&65535)/65535.0f;
			const float yo = (rand()&65535)/65535.0f;
			const __m256 xo8 = _mm256_set1_ps( xo );
			const __m256 yo8 = _mm256_set1_ps( yo );
			const __m256 sc8 = _mm256_set1_ps( scl );
			__m256 x8 = _mm256_mul_ps( sc8, shpx8[shpnr] );
			__m256 y8 = _mm256_mul_ps( sc8, shpy8[shpnr] );
			__m256 X8 = _mm256_mul_ps( sc8, shpX8[shpnr] );
			__m256 Y8 = _mm256_mul_ps( sc8, shpY8[shpnr] );
#if 1
			const float angle = (rand()&65535) / 65535.0f * 2.0f * M_PI;
			rotate_shape( angle, &x8, &y8, &X8, &Y8 );
#endif
			x8 = _mm256_add_ps( x8, xo8 );
			y8 = _mm256_add_ps( y8, yo8 );
			X8 = _mm256_add_ps( X8, xo8 );
			Y8 = _mm256_add_ps( Y8, yo8 );
			valid = valid_placement( x8, y8, X8, Y8, xo, yo );
			trials++;
			if ( valid )
			{
				// add to grid buckets
				if ( scl >= CRITICALSCALE )
				{
					// add to db
					ldbx[ldbsz] = x8;
					ldby[ldbsz] = y8;
					ldbX[ldbsz] = X8;
					ldbY[ldbsz] = Y8;
					ldbsz++;
				}
				else
				{
					sdbx[sdbsz] = x8;
					sdby[sdbsz] = y8;
					sdbX[sdbsz] = X8;
					sdbY[sdbsz] = Y8;
					add_to_buckets( xo, yo, sdbsz );
					sdbsz++;
				}
				fprintf( stderr, "Found placement nr %d (shape %d) in %d trials (trials/shapecount=%f).\n", i, shpnr, trials, trials / (float)(sdbsz+ldbsz) );
				const float radius = sqrtf( ( yo-0.5 ) * ( yo-0.5 ) + ( xo-0.5 ) * ( xo-0.5 ) );
				const float h = 150 + shpnr*25.0f;
				const float s = 0.8 - 1.2 * radius;
				const float v = 0.9 - 1.2 * radius;
				float r,g,b;
				hsv2rgb( h, s, v, &r, &g, &b );
				static float x[8] __attribute__ ((aligned (32)));
				static float y[8] __attribute__ ((aligned (32)));
				_mm256_store_ps( x, x8 );
				_mm256_store_ps( y, y8 );
				fprintf( stdout, "<path d=\"M %f %f L %f %f L %f %f L %f %f L %f %f L %f %f L %f %f L %f %f Z\" fill=\"#%02x%02x%02x\" />\n",
					x[0], y[0],
					x[1], y[1],
					x[2], y[2],
					x[3], y[3],
					x[4], y[4],
					x[5], y[5],
					x[6], y[6],
					x[7], y[7],
					(int) (r*255.999f),
					(int) (g*255.999f),
					(int) (b*255.999f)
				);
				fflush(stdout);
			}
		} while ( !valid );
	}

	for ( int y=1; y<GRIDRES; y++ )
	{
		float yc = y / (float) GRIDRES;
		fprintf( stdout, "<path d=\"M %f %f L %f %f\" stroke=\"blue\" stroke-width=\"0.001\" />\n",
			0.0, yc, 1.0, yc );
	}
	for ( int x=1; x<GRIDRES; x++ )
	{
		float xc = x / (float) GRIDRES;
		fprintf( stdout, "<path d=\"M %f %f L %f %f\" stroke=\"blue\" stroke-width=\"0.001\" />\n",
			xc, 0.0, xc, 1.0 );
	}
	fprintf( stdout, "</svg>\n" );
}

