// Random Space Filling Tiling.
// As described by: http://paulbourke.net/texture_colour/randomtile/
//
// (c)2020 Bram Stolk, Vancouver B.C.

#include <immintrin.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include "prt.h"	// for colour space conversions.
#include "hsv.h"	// for colour space conversions.


typedef float fx16 __attribute__((vector_size(4*16)));	// _mm512  : 16 floats of 32b each.
typedef int   ix16 __attribute__((vector_size(4*16)));	// _mm512i : 16 integers of 32b each.


#define CLAMPED(X,LO,HI) \
	( (X) < (LO) ?  \
		(LO) :  \
		( (X) > (HI) ? \
			(HI) : \
			(X) ) \
	)

__inline__ float interpolate
(
 	float y1, float y2,
	float mu
)
{
	const float mu2 = (1-cosf(mu*M_PI))/2;
	return (y1*(1-mu2)+y2*mu2);
}


// returns <0 for p2 right of line p0-p1
static float onLeft( float p0x, float p0y, float p1x, float p1y, float p2x, float p2y )
{
	return ( ( p1x - p0x ) * ( p2y - p0y ) - ( p2x - p0x ) * ( p1y - p0y ) );
}


// 16-way SIMD version of testing whether p2 lies to left or right on p0-p1 line.
static fx16 onLeft16( fx16 p0x, fx16 p0y, fx16 p1x, fx16 p1y, fx16 p2x, fx16 p2y )
{
	return ( ( p1x - p0x ) * ( p2y - p0y ) - ( p2x - p0x ) * ( p1y - p0y ) );
}


// Scalar version for calculating winding number of a point versus a polygon size 16.
// after: http://geomalgorithms.com/a03-_inclusion.html
static int winding_number( float px, float py, float *vx, float* vy, int n )
{
	int wn=0;
	for ( int i=0; i<n; ++i )
	{
		int j = (i+1)%n;
		const float ol = onLeft( vx[i], vy[i], vx[j], vy[j], px, py );
		int pred0 = (vy[j] >  py );
		int pred1 = (vy[j] <= py );
		//printf("vy[j]=%f py=%f\n", vy[j], py );
		//printf("pred0:%d pred1:%d\n", pred0, pred1);
		const int up_intersect = ( vy[j] >  py && ol > 0.0f ) ?  1 : 0;
		const int dn_intersect = ( vy[j] <= py && ol < 0.0f ) ? -1 : 0;
		const int addition = ( vy[i] <= py ) ? up_intersect : dn_intersect;
		//printf( "up:%d dn:%d addition:%d\n", up_intersect, dn_intersect, addition );
		wn += addition;
	}
	return wn;
}


// 16-way SIMD calculation of winding number of a point versus a polygon size 16.
static int winding_number_16( float px, float py, fx16 vx, fx16 vy, fx16 wx, fx16 wy )
{
	const fx16 zero16    = _mm512_set1_ps( 0 );
	const ix16 zero16i   = _mm512_set1_epi32( 0 );
	const ix16 one16i    = _mm512_set1_epi32( 1 );
	const ix16 min16i    = _mm512_set1_epi32(-1 );

	const fx16 px16 = _mm512_set1_ps( px );
	const fx16 py16 = _mm512_set1_ps( py );

	const fx16 ol = onLeft16( vx, vy, wx, wy, px16, py16 );

	const __mmask16 olltz = _mm512_cmp_ps_mask( ol, zero16, _CMP_LT_OS );
	const __mmask16 olgtz = _mm512_cmp_ps_mask( ol, zero16, _CMP_GT_OS );

	const __mmask16 pred0 = _mm512_cmp_ps_mask( wy, py16, _CMP_GT_OS );
	const __mmask16 pred1 = _mm512_cmp_ps_mask( wy, py16, _CMP_LE_OS );
	const __mmask16 useup = _mm512_cmp_ps_mask( vy, py16, _CMP_LE_OS );

	const ix16 up_intersect = _mm512_mask_blend_epi32( pred0 & olgtz, zero16i, one16i );
	const ix16 dn_intersect = _mm512_mask_blend_epi32( pred1 & olltz, zero16i, min16i );

	const ix16 addition     = _mm512_mask_blend_epi32( useup, dn_intersect, up_intersect );

	return _mm512_reduce_add_epi32( addition );
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


// 16-way SIMD version of edge veruses eight edges intersection.
int edge_vs_edge16( float x1, float y1, float x2, float y2, fx16 vx3, fx16 vy3, fx16 vx4, fx16 vy4 )
{
	const fx16 vx1 = _mm512_set1_ps( x1 );
	const fx16 vy1 = _mm512_set1_ps( y1 );
	const fx16 vx2 = _mm512_set1_ps( x2 );
	const fx16 vy2 = _mm512_set1_ps( y2 );

	const fx16 denom16  = (vy4-vy3) * (vx2-vx1) - (vx4-vx3) * (vy2-vy1);
	const fx16 numera16 = (vx4-vx3) * (vy1-vy3) - (vy4-vy3) * (vx1-vx3);
	const fx16 numerb16 = (vx2-vx1) * (vy1-vy3) - (vy2-vy1) * (vx1-vx3);

	const fx16 mua16 = numera16 / denom16;
	const fx16 mub16 = numerb16 / denom16;

	const fx16 zero16 = _mm512_set1_ps( 0 );
	const fx16 one16  = _mm512_set1_ps( 1 );

	const __mmask16 ok0 = _mm512_cmp_ps_mask( mua16, zero16, _CMP_GE_OS );
	const __mmask16 ok1 = _mm512_cmp_ps_mask( mua16,  one16, _CMP_LE_OS );
	const __mmask16 ok2 = _mm512_cmp_ps_mask( mub16, zero16, _CMP_GE_OS );
	const __mmask16 ok3 = _mm512_cmp_ps_mask( mub16,  one16, _CMP_LE_OS );

	return ok0 & ok1 & ok2 & ok3;

//	const ix16 inbounds16 = (mua16 >= zero16) & (mua16 <= one16) & (mub16 >= zero16) & (mub16 <= one16);

//	return _mm512_reduce_and_epi32( inbounds16 );
}


// 16-Way SIMD test for intersection between 16 edges of one polygon with 16 edges of other polygon.
int edges_intersect
(
	fx16 vxa,
	fx16 vya,
	fx16 nxa,
	fx16 nya,

	fx16 vxb,
	fx16 vyb,
	fx16 nxb,
	fx16 nyb
)
{
	float vxvals[16] __attribute__ ((aligned (64)));
	float vyvals[16] __attribute__ ((aligned (64)));
	float nxvals[16] __attribute__ ((aligned (64)));
	float nyvals[16] __attribute__ ((aligned (64)));

 	_mm512_store_si512( (__m512i*) vxvals, vxa );
 	_mm512_store_si512( (__m512i*) vyvals, vya );
 	_mm512_store_si512( (__m512i*) nxvals, nxa );
 	_mm512_store_si512( (__m512i*) nyvals, nya );

	for ( int i=0; i<16; ++i )
	{
		const int rv = edge_vs_edge16( vxvals[i], vyvals[i], nxvals[i], nyvals[i], vxb, vyb, nxb, nyb );
		//printf( "edge nr %d : intersection result %d\n", i, rv );
		if ( rv ) return 1;
	}
	return 0;
}


static float octx[16] __attribute__ ((aligned (64)));
static float octy[16] __attribute__ ((aligned (64)));
static float octX[16] __attribute__ ((aligned (64)));
static float octY[16] __attribute__ ((aligned (64)));
 
static float strx[16] __attribute__ ((aligned (64)));
static float stry[16] __attribute__ ((aligned (64)));
static float strX[16] __attribute__ ((aligned (64)));
static float strY[16] __attribute__ ((aligned (64)));

#define SHPCNT 7
static float shpx[SHPCNT][16] __attribute__ ((aligned(64)));	// verts i0,i1,..,i14,i15
static float shpy[SHPCNT][16] __attribute__ ((aligned(64)));
static float shpX[SHPCNT][16] __attribute__ ((aligned(64)));	// verts i1,i2,.. i15,i0
static float shpY[SHPCNT][16] __attribute__ ((aligned(64)));
static fx16  shpx16[SHPCNT];
static fx16  shpy16[SHPCNT];
static fx16  shpX16[SHPCNT];
static fx16  shpY16[SHPCNT];
static float shparea[SHPCNT];

#define MAXSZ 1600
// db for large polies
static fx16 ldbx[MAXSZ];
static fx16 ldby[MAXSZ];
static fx16 ldbX[MAXSZ];
static fx16 ldbY[MAXSZ];
static int ldbsz=0;
// db for small polies (which we bin in grid buckets for fast spatial lookup.)
static fx16 sdbx[MAXSZ];
static fx16 sdby[MAXSZ];
static fx16 sdbX[MAXSZ];
static fx16 sdbY[MAXSZ];
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


static int numtests=0;
static int valid_placement( fx16 x16, fx16 y16, fx16 X16, fx16 Y16, float xo, float yo )
{
	static float x[16] __attribute__ ((aligned (64)));
	static float y[16] __attribute__ ((aligned (64)));
	_mm512_store_ps( x, x16 );
	_mm512_store_ps( y, y16 );

	static float X[16] __attribute__ ((aligned (64)));
	static float Y[16] __attribute__ ((aligned (64)));
	_mm512_store_ps( X, X16 );
	_mm512_store_ps( Y, Y16 );

	// consider all currently placed large polygons.
	for ( int i=0; i<ldbsz; ++i )
	{
		numtests++;
		// any of the points inside a placed polygon? If so, invalid placement.
		for ( int j=0; j<16; ++j )
		{
			const int pt_inside = winding_number_16( x[j], y[j], ldbx[i], ldby[i], ldbX[i], ldbY[i] );
			if ( pt_inside ) return 0;
		}
		// any of the edges crosses a placed polygon? If so, invalid placement.
		int ei0 = edges_intersect( x16, y16, X16, Y16, ldbx[i], ldby[i], ldbX[i], ldbY[i] );
		if ( ei0 ) return 0;
	}

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
		numtests++;
		const int i = buckets[gx][gy][k];
		assert( i < MAXSZ );
		// any of the points inside a placed polygon? If so, invalid placement.
		for ( int j=0; j<16; ++j )
		{
			const int pt_inside = winding_number_16( x[j], y[j], sdbx[i], sdby[i], sdbX[i], sdbY[i] );
			if ( pt_inside ) return 0;
		}
		// any of the edges crosses a placed polygon? If so, invalid placement.
		int ei0 = edges_intersect( x16, y16, X16, Y16, sdbx[i], sdby[i], sdbX[i], sdbY[i] );
		if ( ei0 ) return 0;
	}
	if ( cnt < sdbsz )
		fprintf( stderr, "Considered %d+%d instead of %d shapes.\n", ldbsz, cnt, ldbsz+sdbsz );
	return 1;
}


static float areasize( fx16 x16, fx16 y16, fx16 X16, fx16 Y16 )
{
	int hit=0;
	int mis=0;
	for ( int y=-32; y<=32; ++y )
	{
		for ( int x=-32; x<=32; ++x )
		{
			const int rv = winding_number_16( x/32.0f, y/32.0f, x16, y16, X16, Y16 );
			hit += rv;
			mis += !rv;
		}
	}
	return  4 * hit / (float)(hit+mis);
}


static void rotate_shape( float angle, fx16* x, fx16* y, fx16* X, fx16* Y )
{
	const fx16 co = _mm512_set1_ps( cosf( angle ) );
	const fx16 si = _mm512_set1_ps( sinf( angle ) );
	const fx16 newx = co * (*x) - si * (*y);
	const fx16 newy = si * (*x) + co * (*y);
	const fx16 newX = co * (*X) - si * (*Y);
	const fx16 newY = si * (*X) + co * (*Y);
	*x = newx;
	*y = newy;
	*X = newX;
	*Y = newY;
}


int main( int argc, char* argv[] )
{
	srand(23);
	// init shapes 1:sphere, 5: sphere, 6: star
	for ( int i=0; i<16; ++i )
	{
		int j = (i+1)%16;
		const float a0 = i * M_PI * 2 / 16;
		const float a1 = j * M_PI * 2 / 16;
		float r = 0.5f;
		shpx[5][ i ] = r * cosf( a0 );
		shpy[5][ i ] = r * sinf( a0 );
		float r0 = (i&1) ? 0.3f : 0.5f;
		float r1 = (j&1) ? 0.3f : 0.5f;
		shpx[6][ i ] = r0 * cosf( a0 );
		shpy[6][ i ] = r0 * sinf( a0 );
		shpx[1][ i ] = r * cosf( a0 );
		shpy[1][ i ] = r * sinf( a0 );
	}

	for ( int s=0; s<=6; ++s )
	{
		for ( int i=0; i<16; ++i )
		{
			int j = (i+1)%16;
			shpX[s][i] = shpx[s][j];
			shpY[s][i] = shpy[s][j];
		}
		shpx16[s] = _mm512_load_ps( shpx[s] );
		shpy16[s] = _mm512_load_ps( shpy[s] );
		shpX16[s] = _mm512_load_ps( shpX[s] );
		shpY16[s] = _mm512_load_ps( shpY[s] );
		shparea[s] = areasize( shpx16[s], shpy16[s], shpX16[s], shpY16[s] );
		fprintf( stderr, "area of shape %d: %f\n", s, shparea[s] );
	}

	// init grid
	for ( int xx=0; xx<GRIDRES; ++xx )
		for ( int yy=0; yy<GRIDRES; ++yy )
			bucketsizes[xx][yy]=0;

	fprintf( stdout, "<?xml version=\"1.0\"?>\n" );
	fprintf( stdout, "<svg version=\"1.1\" baseProfile=\"tiny\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" id=\"svg-root\" viewBox=\"0 0 1 1\">\n" );

//	const float c = 1.31; // for star.
//	const float c = 1.34;
	const float c = 1.29;
	for ( int i=0; i<MAXSZ; ++i )
	{
		const int shpnr = 0;
		float keys[5] = {
			0.15f + 0.09f * ((rand()&1023) / 1023.0f),
			0.15f + 0.09f * ((rand()&1023) / 1023.0f),
			0.15f + 0.09f * ((rand()&1023) / 1023.0f),
			0.15f + 0.09f * ((rand()&1023) / 1023.0f),
			0,
		};
		keys[4] = keys[0];
		for ( int i=0; i<16; ++i )
		{
			const float a0 = i * M_PI * 2 / 16;
			float r;
			if ( i%4 == 0 )
			{
				r = keys[i/4];
			}
			else
			{
				float k0 = keys[i/4+0];
				float k1 = keys[i/4+1];
				float t = (i%4) * 0.25f;
				r = interpolate( k0, k1, t );
			}
			shpx[shpnr][ i ] = r * cosf( a0 );
			shpy[shpnr][ i ] = r * sinf( a0 );
		}
		// Calc area size
		for ( int i=0; i<16; ++i )
		{
			int j = (i+1)%16;
			shpX[shpnr][i] = shpx[shpnr][j];
			shpY[shpnr][i] = shpy[shpnr][j];
		}
		shpx16[shpnr] = _mm512_load_ps( shpx[shpnr] );
		shpy16[shpnr] = _mm512_load_ps( shpy[shpnr] );
		shpX16[shpnr] = _mm512_load_ps( shpX[shpnr] );
		shpY16[shpnr] = _mm512_load_ps( shpY[shpnr] );
		shparea[shpnr] = areasize( shpx16[shpnr], shpy16[shpnr], shpX16[shpnr], shpY16[shpnr] );
		fprintf( stderr, "area of shape %d: %f\n", shpnr, shparea[shpnr] );

		float scl = sqrtf( powf( 5+i, -c ) / (2*shparea[shpnr]) );
		int valid = 0;
		int trials = 0;
		numtests = 0;
		do 
		{
			float xo = (rand()&0xffffff)/(float)0xffffff;
			float yo = (rand()&0xffffff)/(float)0xffffff;

			const fx16 xo16 = _mm512_set1_ps( xo  );
			const fx16 yo16 = _mm512_set1_ps( yo  );
			const fx16 sc16 = _mm512_set1_ps( scl );
			fx16 x16 = sc16 * shpx16[shpnr];
			fx16 y16 = sc16 * shpy16[shpnr];
			fx16 X16 = sc16 * shpX16[shpnr];
			fx16 Y16 = sc16 * shpY16[shpnr];
#if 0
			if ( i>0 )
			{
				const float angle = (rand()&65535) / 65535.0f * 2.0f * M_PI;
				//const float angle = 22.5 * M_PI / 180.0f;
				rotate_shape( angle, &x16, &y16, &X16, &Y16 );
			}
#endif
			x16 = x16 + xo16;
			y16 = y16 + yo16;
			X16 = X16 + xo16;
			Y16 = Y16 + yo16;
			valid = valid_placement( x16, y16, X16, Y16, xo, yo );
			trials++;
			if ( valid )
			{
				// add to grid buckets
				if ( scl >= CRITICALSCALE )
				{
					// add to db
					ldbx[ldbsz] = x16;
					ldby[ldbsz] = y16;
					ldbX[ldbsz] = X16;
					ldbY[ldbsz] = Y16;
					ldbsz++;
				}
				else
				{
					sdbx[sdbsz] = x16;
					sdby[sdbsz] = y16;
					sdbX[sdbsz] = X16;
					sdbY[sdbsz] = Y16;
					add_to_buckets( xo, yo, sdbsz );
					sdbsz++;
				}
				fprintf( stderr, "Found placement nr %d (shape %d scl %f) in %d trials (trials/shapecount=%f) avg num tests/trial = %d.\n", i, shpnr, scl, trials, trials / (float)(sdbsz+ldbsz), numtests / trials );
				const float radius = sqrtf( ( yo-0.5 ) * ( yo-0.5 ) + ( xo-0.5 ) * ( xo-0.5 ) );
				int colnr = (rand()&1023)/1023.0f < yo ? 0:1;
				const float h = colnr ? 260 : 30;
				const float s = 0.9;
				const float v = 0.8;
				float r,g,b;
				hsv2rgb( h, s, v, &r, &g, &b );
				static float x[16] __attribute__ ((aligned (64)));
				static float y[16] __attribute__ ((aligned (64)));
				_mm512_store_ps( x, x16 );
				_mm512_store_ps( y, y16 );
				fprintf( stdout, "<path d=\"M %f %f L %f %f L %f %f L %f %f L %f %f L %f %f L %f %f L %f %f  L %f %f L %f %f L %f %f L %f %f L %f %f L %f %f L %f %f L %f %f  Z\" fill=\"#%02x%02x%02x\" />\n",
					x[0], y[0],
					x[1], y[1],
					x[2], y[2],
					x[3], y[3],
					x[4], y[4],
					x[5], y[5],
					x[6], y[6],
					x[7], y[7],
					x[8], y[8],
					x[9], y[9],
					x[10], y[10],
					x[11], y[11],
					x[12], y[12],
					x[13], y[13],
					x[14], y[14],
					x[15], y[15],
					(int) (r*255.999f),
					(int) (g*255.999f),
					(int) (b*255.999f)
				);
				fflush(stdout);
			}
		} while ( !valid );
	}

	fprintf( stdout, "</svg>\n" );
}

