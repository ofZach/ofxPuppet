// TriangleMesh.h: interface for the TriangleMesh class.
//
//////////////////////////////////////////////////////////////////////

#ifndef _RMS_TRIANGLEMESH_H_
#define _RMS_TRIANGLEMESH_H_



#define WML_ITEM 

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <float.h>
#include <WmlVector2.h>
#include <WmlVector3.h>

#define TM_VERTEX_STRIDE 4    // normal and color same as vertex
#define TM_TEXTURE_STRIDE 2
#define TM_TRIANGLE_STRIDE 3
#define TM_TRITEXCOORD_STRIDE 6		// 3 2D texture coords

#define TM_RESERVE_INCREASE 256

#define TM_VERTEX_INDEX(base,offset)		(((base)*TM_VERTEX_STRIDE)+(offset))
#define TM_TEXTURE_INDEX(base,offset)		(((base)*TM_TEXTURE_STRIDE)+(offset))
#define TM_TRIANGLE_INDEX(base,offset)		(((base)*TM_TRIANGLE_STRIDE)+(offset))
#define TM_TRITEXCOORD_INDEX(base,offset)	(((base)*TM_TRITEXCOORD_STRIDE)+(offset))

#define TM_VERTEX_BIT   1
#define TM_NORMAL_BIT   2
#define TM_COLOR_BIT    4
#define TM_TEXTURE_BIT  8 
#define TM_TRIANGLE_BIT 16
#define TM_TRIANGLE_TEXCOORD_BIT 32



namespace rmsmesh
{
	
	class TriangleMesh
	{
	public:
		enum FileFormat {OBJ_FORMAT, MESHLITE_FORMAT};
		
		typedef unsigned int VertexID;
		typedef unsigned int TriangleID;
		
		// creation, read, write functions
		
		TriangleMesh();
		TriangleMesh(std::string filename, 
					 unsigned int vertex_sizehint = 0, 
					 unsigned int triangle_sizehint = 0);
		~TriangleMesh();
		
		
		/*
		 *  IMesh interface
		 */
		virtual VertexID AppendVertex( const Wml::Vector3f & vVertex, const Wml::Vector3f * pNormal = NULL )
		{ return AppendVertexData( &vVertex, pNormal, false, NULL, &vVertex ); }
		virtual void GetVertex( VertexID vID, Wml::Vector3f & vVertex, Wml::Vector3f * pNormal = NULL ) const;
		virtual void GetNormal( VertexID vID, Wml::Vector3f & vNormal ) const;	
		virtual unsigned int GetVertexCount() const
		{ return (unsigned int)tm_vertices.size() / TM_VERTEX_STRIDE; }
		virtual unsigned int GetMaxVertexID() const
		{ return (unsigned int)tm_vertices.size() / TM_VERTEX_STRIDE; }
		
		virtual TriangleID AppendTriangle( VertexID v1, VertexID v2, VertexID v3 )
		{ unsigned int nTri[3]; nTri[0] = v1; nTri[1] = v2; nTri[2] = v3; return AppendTriangleData(nTri); }
		
		virtual bool SetTriangle( TriangleID tID, VertexID v1, VertexID v2, VertexID v3 )
		{ unsigned int nTri[3]; nTri[0] = v1; nTri[1] = v2; nTri[2] = v3; SetTriangle(tID, nTri); return true; }
		
		virtual void GetTriangle( TriangleID tID, VertexID vTriangle[3]  ) const;
		virtual void GetTriangle( TriangleID tID, Wml::Vector3f vTriangle[3], Wml::Vector3f * pNormals = NULL  ) const;
		virtual unsigned int GetTriangleCount() const
		{ return (unsigned int)tm_triangles.size() / TM_TRIANGLE_STRIDE; }
		virtual unsigned int GetMaxTriangleID() const
		{ return (unsigned int)tm_triangles.size() / TM_TRIANGLE_STRIDE; }
		
		virtual void Clear( bool bFreeMem )
		{ Clear(); }	// never frees mem!
		
		
		bool
		read(const char * pFilename, FileFormat eFormat = OBJ_FORMAT);
		
		bool
		write(const char * pFilename, FileFormat eFormat = OBJ_FORMAT);
		
		// if you are going to use the setX functions, make sure you resize()
		// the mesh so it is large enough to hold all your elements (!)
		
		void
		SetVertexData(unsigned int index, float * vertex = NULL,
					  float * normal = NULL, float * color = NULL,
					  float * texture_coord = NULL);
		
		void
		SetTriangleData(unsigned int index, unsigned int * triangle);
		
		
		// addVertexData/addTriangleData are safe, in that they check to see
		// if the mesh is big enough before they add anything
		
		void
		AddVertexData(unsigned int index, float * vertex = NULL,
					  float * normal = NULL, float * color = NULL,
					  float * texture_coord = NULL);
		
		void
		AddTriangleData(unsigned int index, unsigned int * triangle);
		
		void
		AddTriTexCoordData(unsigned int nIndex, const float * pUV1, const float * pUV2, const float * pUV3 );
		
		
		// appendVertexData/appendTriangleData are safe as well.
		// The return value is the vertex/triangle index. Note that if you
		// are using seperate indices for vertices/normals/colors/textures, 
		// it returns the index for the last non-null data
		
		unsigned int
		AppendVertexData(float * vertex = NULL, float * normal = NULL, 
						 float * color = NULL, float * texture_coord = NULL);
		
		unsigned int
		AppendVertexData( const Wml::Vector3f * pVertex, const Wml::Vector3f * pNormal, bool bFlipNormal, 
						 const Wml::Vector2f * pTextureCoord = NULL, const Wml::Vector3f * pColor = NULL );
		
		unsigned int
		AppendTriangleData(unsigned int * triangle);
		
		
		// resize makes sure that the mesh is big enough to hold new_count
		// reserve is more of a hint that you'll need at least this many elements
		// build the mask from the TM_XX_BIT constants defined above
		
		void 
		Resize(unsigned int mask, unsigned int new_count);
		
		void
		Reserve(unsigned int mask, unsigned int new_count);
		
		void 
		Clear() {
			tm_vertices.resize(0);
			tm_normals.resize(0);
			tm_colors.resize(0);
			tm_texture_coords.resize(0);
			tm_triangles.resize(0);
			tm_triTexCoords.resize(0);
		}
		
		void Clear( unsigned int nBitMask );
		
		
		// other bits
		
		// minx/maxx, miny/maxy, minz/maxz
		void GetBoundingBox( float box[6] );
		
		float GetMaxEdgeLength() const;
		
		// tests
		bool HasVertexTextureCoords() const;
		bool HasTriangleTextureCoords() const;
		
		// data accessor functions
		unsigned int GetNumVertices() const 
		{ return (unsigned int)tm_vertices.size() / TM_VERTEX_STRIDE; }
		unsigned int GetNumTriangles() const 
		{ return (unsigned int)tm_triangles.size() / TM_TRIANGLE_STRIDE; }
		
		//	void GetVertex( unsigned int nVertex, Wml::Vector3f & v ) const;
		//	void GetNormal( unsigned int nVertex, Wml::Vector3f & n ) const;
		void GetTextureCoords( unsigned int nVertex, Wml::Vector2f & vt ) const;
		//void GetTriangle( unsigned int nTriangle, unsigned int triangle[3] ) const;
		//	void GetTriangle( unsigned int nTriangle, Wml::Vector3f * vVertices ) const;
		void GetTextureCoords( unsigned int nTriangle, Wml::Vector2f * vCoords ) const;
		void GetTriangleNormals( unsigned int nTriangle, Wml::Vector3f * vNormals );
		void GetTriTexCoords( unsigned int nTriangle, Wml::Vector2f & uv1, Wml::Vector2f & uv2, Wml::Vector2f & uv3 ) const;
		
		void SetTriangle( unsigned int nTriangle, unsigned int triangle[3] );
		void SetVertex( unsigned int nVertex, const Wml::Vector3f & vVertex );
		void SetNormal( unsigned int nVertex, const Wml::Vector3f & vNormal );
		
		// erase triangles in the bitmap
		void EraseTriangles( const std::vector<bool> & vErase);
		
		
		inline std::vector<float> & 
		GetVertices() { return tm_vertices; } 
		inline const std::vector<float> & 
		GetVertices() const { return tm_vertices; } 
		
		inline std::vector<float> & 
		GetNormals() { return tm_normals; } 
		inline const std::vector<float> & 
		GetNormals() const { return tm_normals; } 
		
		inline std::vector<float> & 
		GetColors() { return tm_colors; } 
		inline const std::vector<float> & 
		GetColors() const { return tm_colors; } 
		
		inline std::vector<float> & 
		GetTextureCoords() { return tm_texture_coords; } 
		inline const std::vector<float> & 
		GetTextureCoords() const { return tm_texture_coords; } 
		
		inline std::vector<unsigned int> &
		GetTriangles() { return tm_triangles; }
		inline const std::vector<unsigned int> &
		GetTriangles() const { return tm_triangles; }
		
		inline std::vector<float> & 
		GetTriTexCoords() { return tm_triTexCoords; } 
		inline const std::vector<float> & 
		GetTriTexCoords() const { return tm_triTexCoords; } 
		
		inline float *
		GetVertexPointer() { return &(*(tm_vertices.begin())); }
		inline const float *
		GetVertexPointer() const { return &(*(tm_vertices.begin())); }
		
		inline float *
		GetNormalPointer() { return &(*(tm_normals.begin())); }
		inline const float *
		GetNormalPointer() const { return &(*(tm_normals.begin())); }
		
		inline float *
		GetColorPointer() { return &(*(tm_colors.begin())); }
		inline const float *
		GetColorPointer() const { return &(*(tm_colors.begin())); }
		
		inline float *
		GetTextureCoordsPointer() { return &(*(tm_texture_coords.begin())); }
		inline const float *
		GetTextureCoordsPointer() const { return &(*(tm_texture_coords.begin())); }
		
		inline unsigned int *
		GetTrianglePointer() { return &(*(tm_triangles.begin())); }
		inline const unsigned int *
		GetTrianglePointer() const { return &(*(tm_triangles.begin())); }
		
		inline float *
		GetTriTexCoordPointer() { return &(*(tm_triTexCoords.begin())); }
		inline const float *
		GetTriTexCoordPointer() const { return &(*(tm_triTexCoords.begin())); }
		
		// status functions
		
		inline std::string const &
		GetFileName() const { return tm_filename; } 
		
		inline std::string const &
		GetError() const { return tm_fileerror; }
		
		inline enum FileFormat 
		GetFileFormat() const { return tm_fileformat; }
		
	protected:
		std::string tm_filename;
		std::string tm_fileerror;
		enum FileFormat tm_fileformat;
		
		std::vector<float> tm_vertices;
		std::vector<float> tm_normals;
		std::vector<float> tm_colors;
		std::vector<float> tm_texture_coords;
		std::vector<unsigned int > tm_triangles;
		
		std::vector<float> tm_triTexCoords;
		
		bool tm_readobj();
		
		bool tm_writeobj();
		bool tm_writeMeshLite();
		
		
		void tm_initialize(unsigned int vertex_sizehint, unsigned int triangle_sizehint);
		
		
		/*
		 * IMesh iterator interface
		 */
		virtual void * ivtx_make_iterator(bool bStart)
		{ return NULL; }
		virtual void * ivtx_make_iterator( void * pFromItr )
		{ return NULL; }
		virtual void ivtx_free_iterator( void * pItr ) 
		{ }
		virtual void ivtx_set( void * pItr, void * pTo )
		{ }
		virtual void ivtx_goto_next( void * pItr )
		{ }
		virtual bool ivtx_equal( void * pItr1, void * pItr2 )
		{ return false; }
		virtual VertexID ivtx_value( void * pItr )
		{ return 0; }
		
		virtual void * itri_make_iterator(bool bStart)
		{ return NULL; }
		virtual void * itri_make_iterator( void * pFromItr )
		{ return NULL; }
		virtual void itri_free_iterator( void * pItr ) 
		{ }
		virtual void itri_set( void * pItr, void * pTo )
		{ }
		virtual void itri_goto_next( void * pItr )
		{ }
		virtual bool itri_equal( void * pItr1, void * pItr2 )
		{ return false; }
		virtual TriangleID itri_value( void * pItr )
		{ return 0; }
	};
	
}  // end namespace rms

#endif // _RMS_TRIANGLEMESH_H_
