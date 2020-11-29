#!/usr/bin/python3

import sys


def process_blob( blob, indexoffset ) :
	lines = blob.split( '\n' )
	lines = [ l.strip() for l in lines ]
	name = lines[ 0 ]
	nr = int( name[1:], 16 )
	vlines = [ l for l in lines if "v " in l ]
	verts = []
	faces = []
	trias = []
	width=0
	heigh=0
	lox= 10000
	hix=-10000
	loy= 10000
	hiy=-10000
	for vline in vlines:
		fields = vline.split( " " )
		x = -float( fields[2] )
		y =  float( fields[3] )
		v = ( x, y )
		verts.append( v )
		lox = min( lox, x )
		hix = max( hix, x )
		loy = min( loy, y )
		hiy = max( hiy, y )
	midx = 0.5 * ( lox + hix )
	midy = 0.5 * ( loy + hiy )
	width = hix - lox
	heigh = hiy - loy
	verts = [ (v[0]-midx, -(v[1]-midy)) for v in verts ]
	flines = [ l for l in lines if "f " in l ]
	for fline in flines :
		fields = fline.split( " " )[ 1: ]
		indices = [ int( x.split("/")[0] ) - indexoffset for x in fields ]
		faces.append( indices )
	vstream = []
	assert len(faces) == 1
	for face in faces :
		numv = len( face )
		assert len( face ) == 16
		face.reverse()
		for j in range( len(face) ) :
				i0 = face [ j ]
				v0 = verts[ i0 ]
				vstream.append( v0 )
		return ( nr, width, vstream, indexoffset+len(verts) )


# main entry point:

if len( sys.argv ) != 2 :
	print( "Usage: %s foo.obj" % ( sys.argv[0], )  )
	sys.exit( 1 )

f = open( sys.argv[1], 'r' )
assert f
blob = f.read()
f.close()

blobs = blob.split( "o " )[ 1: ]
#assert len(blobs) == 128-1

widths  = [ 0 for x in range(128) ]
sizes   = [ 0 for x in range(128) ]
streams = [ None for x in range(128 ) ]

indexoffset = 1

for blob in blobs :
	nr, width, vstream, indexoffset = process_blob( blob, indexoffset )
	widths [ nr ] = width
	sizes  [ nr ] = len(vstream)
	streams[ nr ] = vstream


#assert not -1 in widths
#assert not -1 in sizes

totalsize = sum( sizes )

print( "// Machine-generated from %s by tools/fnt.py, do not edit." % ( sys.argv[1], ) )
print( '#pragma clang diagnostic ignored "-Wmissing-braces"' )
print( '#pragma clang diagnostic ignored "-Wconversion"' )
print( "" )
print( "#define VDATASZ %d" % ( totalsize, ) )
print( "#define NUMGLYPHS 128" )

print( "// vertex data for all glyphs combined." )
print( "static float vdatax[] =\n{" )
for i in range(128) :
	sz = sizes[i]
	if sz > 0:
		stream = streams[i]
		for v in stream:
			sys.stdout.write( "%.3f," % ( v[0], ), )
		sys.stdout.write( "\n" )
print( "};" )
print( "static float vdatay[] =\n{" )
for i in range(128) :
	sz = sizes[i]
	if sz > 0:
		stream = streams[i]
		for v in stream:
			sys.stdout.write( "%.3f," % ( v[1], ), )
		sys.stdout.write( "\n" )
print( "};" )


