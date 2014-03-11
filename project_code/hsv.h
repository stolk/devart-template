// From:
// http://stackoverflow.com/a/6930407/301166

inline void rgb2hsv( float r, float g, float b, float* __restrict h, float* __restrict s, float* __restrict v )
{
	float min, max, delta;

	min = r < g ? r : g;
	min = min  < b ? min : b;

	max = r > g ? r : g;
	max = max  > b ? max : b;

	*v = max;
	delta = max - min;
	if( max > 0.0f ) 
	{
		*s = (delta / max);
	} 
	else
	{
		// r = g = b = 0
		*s = 0.0;
		*h = NAN;
		return;
	}
	if ( r >= max )
		*h = ( g - b ) / delta;        // between yellow & magenta
	else
		if ( g >= max )
			*h = 2.0 + ( b - r ) / delta;  // between cyan & yellow
		else
			*h = 4.0 + ( r - g ) / delta;  // between magenta & cyan

	*h *= 60.0; // degrees

	if ( *h < 0.0 )
		*h += 360.0;
}


inline void hsv2rgb( float h, float s, float v, float* __restrict r, float* __restrict g, float* __restrict b )
{
	float hh, p, q, t, ff;
	int i;

	if( s <= 0.0f ) 
	{
		*r = v;
		*g = v;
		*b = v;
		return;
	}
	hh = h;
	if( hh >= 360.0f ) hh = 0.0;
	hh /= 60.0f;
	i = (int)hh;
	ff = hh - i;
	p = v * (1.0f - s);
	q = v * (1.0f - (s * ff));
	t = v * (1.0f - (s * (1.0 - ff)));

	switch( i )
	{
		case 0:
			*r = v;
			*g = t;
			*b = p;
			break;
		case 1:
			*r = q;
			*g = v;
			*b = p;
			break;
		case 2:
			*r = p;
			*g = v;
			*b = t;
			break;

		case 3:
			*r = p;
			*g = q;
			*b = v;
			break;
		case 4:
			*r = t;
			*g = p;
			*b = v;
			break;
		case 5:
		default:
			*r = v;
			*g = p;
			*b = q;
			break;
	}
}


