// #include "../Src/PreProcessor.h"
#include "PoissonReconLib.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "MyMiscellany.h"
// #include "../Src/CmdLineParser.h"
#include "PPolynomial.h"
#include "FEMTree.h"
// #include "../Src/Ply.h"
#include "VertexFactory.h"
// #include "Image.h"
#include "RegularGrid.h"
using namespace VertexFactory;
template< typename Data >
VectorInputDataStream< Data >::VectorInputDataStream( const std::vector<Data> data ) : _data(data) , _size(data.size()) , _current(0) {}
template< typename Data >
VectorInputDataStream< Data >::~VectorInputDataStream( void ){ ; }
template< typename Data >
void VectorInputDataStream< Data >::reset( void ) { _current=0; }
template< typename Data >
bool VectorInputDataStream< Data >::next( Data &d )
{
	if( _current>=_size ) return false;
	d = _data[_current++];
	return true;
}
struct ParamBool {
    const bool set=false;
    ParamBool(){
    }
};
template<typename Type>
struct Param {
    Type value;
    bool set;
    const char* name="";
    Param(){

    }
    Param( Type v){
        value = v;
    }
};
const Param<float> SamplesPerNode={1.5f};
const Param<int> Iters={8};
const Param<float> DataX={32.f};
const Param<float> CGSolverAccuracy={1e-3f};
const Param<float> Scale={1.1f};
const Param<float> Confidence = {0.f};
const Param<float> ConfidenceBias={0.f};
const Param<float> LowDepthCutOff={0.f};
const Param<float> PointWeight;
enum NormalType
{
	NORMALS_NONE ,
	NORMALS_SAMPLES ,
	NORMALS_GRADIENTS ,
	NORMALS_COUNT
};
const Param<int> BaseVCycles={1};
const Param<int> Normals={NORMALS_NONE};
const Param<int> KernelDepth;
const ParamBool Envelope;
const ParamBool EnvelopeGrid;
const Param<int> EnvelopeDepth;
Param<int> Depth = {8};
Param<float> BaseDepth={0.f};
Param<float> SolveDepth={0.f};
Param<float> FullDepth={0.f};
Param<float> Width={0.f};
const ParamBool Verbose;
const ParamBool Density;
const ParamBool NoDirichletErode;
const ParamBool ShowResidual;
const ParamBool Out;
const ParamBool Grid;
const ParamBool Tree;
const ParamBool PrimalGrid;
const ParamBool ExactInterpolation;
const ParamBool LinearFit;
const ParamBool NonManifold;
const ParamBool PolygonMesh;
const Param<int> Threads={(int)std::thread::hardware_concurrency()};

template< class Real ,typename Index ,typename AuxDataFactory , unsigned int ... FEMSigs >
class PoissonReconLib{
#ifndef NOLIB
	static const int Dim = sizeof ... ( FEMSigs );
#endif
#ifndef NOLIB
	// The input point stream information: First piece of data is the normal; the remainder is the auxiliary data
	typedef InputOrientedPointStreamInfo< Real , Dim , typename AuxDataFactory::VertexType > InputPointStreamInfo;
	// The type of the input sample
	typedef typename InputPointStreamInfo::PointAndDataType InputSampleType;
	// The type of the input sample's data
	typedef typename InputPointStreamInfo::DataType InputSampleDataType;
	typedef InputDataStream< InputSampleType >  InputPointStream;
	typedef TransformedInputDataStream< InputSampleType > XInputPointStream;

    typedef Factory< Real , PositionFactory< Real , Dim > , NormalFactory< Real , Dim > , ValueFactory< Real > , AuxDataFactory > VertexFactory;
#endif
	public:
PoissonReconLib( UIntPack< FEMSigs ... > ){

}
template< typename VertexFactory ,typename OutputIndex , bool UseCharIndex >
	void OutputPolygons(const VertexFactory &vFactory , StreamingMesh< typename VertexFactory::VertexType , Index >* mesh ,std::function< void ( typename VertexFactory::VertexType & ) > xForm,MeshOutputDataStream<typename VertexFactory::VertexType,Index>& output)
	{
        if( mesh->vertexNum()>(size_t)std::numeric_limits< OutputIndex >::max() )
		{
#if NOLIB
			if( std::is_same< Index , OutputIndex >::value ) ERROR_OUT( "more vertices than can be represented using " , Traits< Index >::name );
			WARN( "more vertices than can be represented using " , Traits< OutputIndex >::name , " using " , Traits< Index >::name , " instead" );
			return WritePolygons< VertexFactory , Index , Real , Dim , Index >( fileName , vFactory , mesh , file_type , comments , xForm );
#else
			if( std::is_same< Index , OutputIndex >::value ) throw "more vertices than can be represented using ";
			return WritePolygons< VertexFactory , Index , Real , Dim , Index >( vFactory , mesh, xForm );
#endif
		}
		size_t nr_vertices = mesh->vertexNum();
		size_t nr_faces = mesh->polygonNum();
		float version;
        output.set_vertex_num(mesh->vertexNum);
        output.set_polygon_num(mesh->polygonNum);
		mesh->resetIterator();
		if( vFactory.isStaticallyAllocated() )
		{
			for( size_t i=0; i<mesh->vertexNum() ; i++ )
			{
				typename VertexFactory::VertexType vertex = vFactory();
				mesh->nextVertex( vertex );
				xForm( vertex );
                output.push_vertex(vertex);
			}
		}
		else
		{
			for( size_t i=0; i<mesh->vertexNum() ; i++ )
			{
				typename VertexFactory::VertexType vertex = vFactory();
				mesh->nextVertex( vertex );
				xForm( vertex );
                output.push_vertex(vertex);
			}
		}

	   // write faces
		std::vector< Index > polygon;
		for( size_t i=0 ; i<nr_faces ; i++ )
		{
			//
			// create and fill a struct that the ply code can handle
			//
			mesh->nextPolygon( polygon );
            output.push_polygon_indices(polygon);
		}  // for, write faces
    }
template< typename SetVertexFunction , typename InputSampleDataType , typename VertexFactory>
void ExtractMesh
(
	UIntPack< FEMSigs ... > ,
	FEMTree< sizeof ... ( FEMSigs ) , Real >& tree ,
	const DenseNodeData< Real , UIntPack< FEMSigs ... > >& solution ,
	Real isoValue ,
	const std::vector< typename FEMTree< sizeof ... ( FEMSigs ) , Real >::PointSample > *samples ,
	std::vector< InputSampleDataType > *sampleData ,
	const typename FEMTree< sizeof ... ( FEMSigs ) , Real >::template DensityEstimator< WEIGHT_DEGREE > *density ,
	const VertexFactory &vertexFactory ,
	const InputSampleDataType &zeroInputSampleDataType ,
	SetVertexFunction SetVertex ,
	std::vector< std::string > &comments ,
	XForm< Real , sizeof...(FEMSigs)+1 > unitCubeToModel,
    MeshOutputDataStream<typename VertexFactory::VertexType,Index>& output
)
{
	static const int Dim = sizeof ... ( FEMSigs );
	typedef UIntPack< FEMSigs ... > Sigs;
	typedef typename VertexFactory::VertexType Vertex;

	static const unsigned int DataSig = FEMDegreeAndBType< DATA_DEGREE , BOUNDARY_FREE >::Signature;
	typedef typename FEMTree< Dim , Real >::template DensityEstimator< WEIGHT_DEGREE > DensityEstimator;

	Profiler profiler(20);
#ifdef NOLIB
	char tempHeader[2048];
	{
		char tempPath[1024];
		tempPath[0] = 0;
		if( TempDir.set ) strcpy( tempPath , TempDir.value );
		else SetTempDirectory( tempPath , sizeof(tempPath) );
		if( strlen(tempPath)==0 ) sprintf( tempPath , ".%c" , FileSeparator );
		if( tempPath[ strlen( tempPath )-1 ]==FileSeparator ) sprintf( tempHeader , "%sPR_" , tempPath );
		else                                                  sprintf( tempHeader , "%s%cPR_" , tempPath , FileSeparator );
	}
	StreamingMesh< Vertex , node_index_type > *mesh;
	if( InCore.set ) mesh = new VectorStreamingMesh< Vertex , node_index_type >();
	else             mesh = new FileStreamingMesh< VertexFactory , node_index_type >( vertexFactory , tempHeader );
#else
    StreamingMesh< Vertex , Index > *mesh;
    mesh = new VectorStreamingMesh< Vertex , Index >();
#endif
	profiler.reset();
	typename LevelSetExtractor< Dim , Real , Vertex >::Stats stats;
	if( sampleData )
	{
		SparseNodeData< ProjectiveData< InputSampleDataType , Real > , IsotropicUIntPack< Dim , DataSig > > _sampleData = tree.template setExtrapolatedDataField< DataSig , false >( *samples , *sampleData , (DensityEstimator*)NULL );
		auto nodeFunctor = [&]( const RegularTreeNode< Dim , FEMTreeNodeData , depth_and_offset_type > *n )
		{
			ProjectiveData< InputSampleDataType , Real >* clr = _sampleData( n );
			if( clr ) (*clr) *= (Real)pow( DataX.value , tree.depth( n ) );
		};
		tree.tree().processNodes( nodeFunctor );
		stats = LevelSetExtractor< Dim , Real , Vertex >::template Extract< InputSampleDataType >( Sigs() , UIntPack< WEIGHT_DEGREE >() , UIntPack< DataSig >() , tree , density , &_sampleData , solution , isoValue , *mesh , zeroInputSampleDataType , SetVertex , !LinearFit.set , Normals.value==NORMALS_GRADIENTS , !NonManifold.set , PolygonMesh.set , false );
	}
#if defined( __GNUC__ ) && __GNUC__ < 5
#ifdef SHOW_WARNINGS
#warning "you've got me gcc version<5"
#endif // SHOW_WARNINGS
	else stats = LevelSetExtractor< Dim , Real , Vertex >::template Extract< InputSampleDataType >( Sigs() , UIntPack< WEIGHT_DEGREE >() , UIntPack< DataSig >() , tree , density , (SparseNodeData< ProjectiveData< InputSampleDataType , Real > , IsotropicUIntPack< Dim , DataSig > > *)NULL , solution , isoValue , *mesh , zeroInputSampleDataType , SetVertex , !LinearFit.set , Normals.value==NORMALS_GRADIENTS , !NonManifold.set , PolygonMesh.set , false );
#else // !__GNUC__ || __GNUC__ >=5
	else stats = LevelSetExtractor< Dim , Real , Vertex >::template Extract< InputSampleDataType >( Sigs() , UIntPack< WEIGHT_DEGREE >() , UIntPack< DataSig >() , tree , density , NULL , solution , isoValue , *mesh , zeroInputSampleDataType , SetVertex , !LinearFit.set , Normals.value==NORMALS_GRADIENTS , !NonManifold.set , PolygonMesh.set , false );
#endif // __GNUC__ || __GNUC__ < 4
	if( Verbose.set )
	{
		std::cout << "Vertices / Polygons: " << mesh->vertexNum() << " / " << mesh->polygonNum() << std::endl;
		std::cout << stats.toString() << std::endl;
		if( PolygonMesh.set ) std::cout << "#         Got polygons: " << profiler << std::endl;
		else                  std::cout << "#        Got triangles: " << profiler << std::endl;
	}

	std::vector< std::string > noComments;
	typename VertexFactory::Transform unitCubeToModelTransform( unitCubeToModel );
	auto xForm = [&]( typename VertexFactory::VertexType & v ){ unitCubeToModelTransform.inPlace( v ); };
#ifdef NOLIB
	PLY::WritePolygons< VertexFactory , node_index_type , Real , Dim >( Out.value , vertexFactory , mesh , ASCII.set ? PLY_ASCII : PLY_BINARY_NATIVE , NoComments.set ? noComments : comments , xForm );
#else
	OutputPolygons< VertexFactory , Real , Dim >( vertexFactory , mesh , xForm ,output);
#endif
	delete mesh;
}
void Execute( UIntPack< FEMSigs ... > , const AuxDataFactory &auxDataFactory )
{
#ifdef NOLIB
    static const int Dim = sizeof ... ( FEMSigs );
#endif
	typedef UIntPack< FEMSigs ... > Sigs;
	typedef UIntPack< FEMSignature< FEMSigs >::Degree ... > Degrees;
	typedef UIntPack< FEMDegreeAndBType< NORMAL_DEGREE , DerivativeBoundary< FEMSignature< FEMSigs >::BType , 1 >::BType >::Signature ... > NormalSigs;
	static const unsigned int DataSig = FEMDegreeAndBType< DATA_DEGREE , BOUNDARY_FREE >::Signature;
	typedef typename FEMTree< Dim , Real >::template DensityEstimator< WEIGHT_DEGREE > DensityEstimator;
	typedef typename FEMTree< Dim , Real >::template InterpolationInfo< Real , 0 > InterpolationInfo;
#ifdef NOLIB
	using namespace VertexFactory;
#endif
	// The factory for constructing an input sample
	typedef Factory< Real , PositionFactory< Real , Dim > , Factory< Real , NormalFactory< Real , Dim > , AuxDataFactory > > InputSampleFactory;

	// The factory for constructing an input sample's data
	typedef Factory< Real , NormalFactory< Real , Dim > , AuxDataFactory > InputSampleDataFactory;

#ifdef NOLIB
	// The input point stream information: First piece of data is the normal; the remainder is the auxiliary data
	typedef InputOrientedPointStreamInfo< Real , Dim , typename AuxDataFactory::VertexType > InputPointStreamInfo;

	// The type of the input sample
	typedef typename InputPointStreamInfo::PointAndDataType InputSampleType;

	// The type of the input sample's data
	typedef typename InputPointStreamInfo::DataType InputSampleDataType;

	typedef            InputDataStream< InputSampleType >  InputPointStream;
	typedef TransformedInputDataStream< InputSampleType > XInputPointStream;
#endif
	InputSampleFactory inputSampleFactory( PositionFactory< Real , Dim >() , InputSampleDataFactory( NormalFactory< Real , Dim >() , auxDataFactory ) );
	InputSampleDataFactory inputSampleDataFactory( NormalFactory< Real , Dim >() , auxDataFactory );

	typedef RegularTreeNode< Dim , FEMTreeNodeData , depth_and_offset_type > FEMTreeNode;
	typedef typename FEMTreeInitializer< Dim , Real >::GeometryNodeType GeometryNodeType;
	std::vector< std::string > comments;
	if( Verbose.set )
	{
		std::cout << "*************************************************************" << std::endl;
		std::cout << "*************************************************************" << std::endl;
		std::cout << "** Running Screened Poisson Reconstruction (Version " << VERSION << ") **" << std::endl;
		std::cout << "*************************************************************" << std::endl;
		std::cout << "*************************************************************" << std::endl;
		if( !Threads.set ) std::cout << "Running with " << Threads.value << " threads" << std::endl;
	}

	bool needNormalData = DataX.value>0 && Normals.value;
	bool needAuxData = DataX.value>0 && auxDataFactory.bufferSize();

	XForm< Real , Dim+1 > modelToUnitCube , unitCubeToModel;
#ifdef NOLIB
	if( Transform.set )
	{
		FILE* fp = fopen( Transform.value , "r" );
		if( !fp )
		{
			WARN( "Could not read x-form from: " , Transform.value );
			modelToUnitCube = XForm< Real , Dim+1 >::Identity();
		}
		else
		{
			for( int i=0 ; i<Dim+1 ; i++ ) for( int j=0 ; j<Dim+1 ; j++ )
			{
				float f;
				if( fscanf( fp , " %f " , &f )!=1 ) ERROR_OUT( "Failed to read xform" );
				modelToUnitCube(i,j) = (Real)f;
			}
			fclose( fp );
		}
	}
	else modelToUnitCube = XForm< Real , Dim+1 >::Identity();
#else
    modelToUnitCube = XForm< Real , Dim+1 >::Identity();
#endif

#ifdef NOLIB
	char str[1024];
	for( int i=0 ; params[i] ; i++ )
		if( params[i]->set )
		{
			params[i]->writeValue( str );
			if( Verbose.set )
				if( strlen( str ) ) std::cout << "\t--" << params[i]->name << " " << str << std::endl;
				else                std::cout << "\t--" << params[i]->name << std::endl;
		}
#endif
	double startTime = Time();
	Real isoValue = 0;

	FEMTree< Dim , Real > tree( MEMORY_ALLOCATOR_BLOCK_SIZE );
	Profiler profiler(20);

	if( Depth.set && Width.value>0 )
	{
		WARN( "Both --" , Depth.name  , " and --" , Width.name , " set, ignoring --" , Width.name );
		Width.value = 0;
	}

	size_t pointCount;

	ProjectiveData< Point< Real , 2 > , Real > pointDepthAndWeight;
	std::vector< typename FEMTree< Dim , Real >::PointSample >* samples = new std::vector< typename FEMTree< Dim , Real >::PointSample >();
	DenseNodeData< GeometryNodeType , IsotropicUIntPack< Dim , FEMTrivialSignature > > geometryNodeDesignators;
	std::vector< InputSampleDataType >* sampleData = NULL;
	DensityEstimator* density = NULL;
	SparseNodeData< Point< Real , Dim > , NormalSigs >* normalInfo = NULL;
	Real targetValue = (Real)0.5;

	// Read in the samples (and color data)
	{
		profiler.reset();
		InputPointStream* pointStream;
#ifdef NOLIB
		char* ext = GetFileExtension( In.value );
#endif
		sampleData = new std::vector< InputSampleDataType >();
#ifdef NOLIB
		std::vector< InputSampleType > inCorePoints;

		if( InCore.set )
		{
			InputPointStream *_pointStream;
			if     ( !strcasecmp( ext , "bnpts" ) ) _pointStream = new BinaryInputDataStream< InputSampleFactory >( In.value , inputSampleFactory );
			else if( !strcasecmp( ext , "ply"   ) ) _pointStream = new    PLYInputDataStream< InputSampleFactory >( In.value , inputSampleFactory );
			else                                    _pointStream = new  ASCIIInputDataStream< InputSampleFactory >( In.value , inputSampleFactory );
			InputSampleType p = inputSampleFactory();
			while( _pointStream->next( p ) ) inCorePoints.push_back( p );
			delete _pointStream;

			pointStream = new MemoryInputDataStream< InputSampleType >( inCorePoints.size() , &inCorePoints[0] );
		}
		else
		{
			if     ( !strcasecmp( ext , "bnpts" ) ) pointStream = new BinaryInputDataStream< InputSampleFactory >( In.value , inputSampleFactory );
			else if( !strcasecmp( ext , "ply"   ) ) pointStream = new    PLYInputDataStream< InputSampleFactory >( In.value , inputSampleFactory );
			else                                    pointStream = new  ASCIIInputDataStream< InputSampleFactory >( In.value , inputSampleFactory );
		}
		delete[] ext;
#endif

		typename InputSampleFactory::Transform _modelToUnitCube( modelToUnitCube );
		auto XFormFunctor = [&]( InputSampleType &p ){ _modelToUnitCube.inPlace( p ); };
		XInputPointStream _pointStream( XFormFunctor , *pointStream );

		{
			typename InputSampleDataFactory::VertexType zeroData = inputSampleDataFactory();
			modelToUnitCube = Scale.value>0 ? GetPointXForm< Real , Dim , typename AuxDataFactory::VertexType >( _pointStream , zeroData , (Real)Scale.value ) * modelToUnitCube : modelToUnitCube;
		}
		if( Width.value>0 )
		{
			Real maxScale = 0;
			for( unsigned int i=0 ; i<Dim ; i++ ) maxScale = std::max< Real >( maxScale , (Real)1./modelToUnitCube(i,i) );
			Depth.value = (unsigned int)ceil( std::max< double >( 0. , log( maxScale/Width.value )/log(2.) ) );
		}
		if( SolveDepth.value>Depth.value )
		{
			WARN( "Solution depth cannot exceed system depth: " , SolveDepth.value , " <= " , Depth.value );
			SolveDepth.value = Depth.value;
		}
		if( FullDepth.value>Depth.value )
		{
			WARN( "Full depth cannot exceed system depth: " , FullDepth.value , " <= " , Depth.value );
			FullDepth.value = Depth.value;
		}
		if( BaseDepth.value>FullDepth.value )
		{
			if( BaseDepth.set ) WARN( "Base depth must be smaller than full depth: " , BaseDepth.value , " <= " , FullDepth.value );
			BaseDepth.value = FullDepth.value;
		}

		{
			typename InputSampleFactory::Transform _modelToUnitCube( modelToUnitCube );
			auto XFormFunctor = [&]( InputSampleType &p ){ _modelToUnitCube.inPlace( p ); };
			XInputPointStream _pointStream( XFormFunctor , *pointStream );
			auto ProcessDataWithConfidence = [&]( const Point< Real , Dim > &p , typename InputPointStreamInfo::DataType &d )
			{
				Real l = (Real)Length( d.template get<0>() );
				if( !l || !std::isfinite( l ) ) return (Real)-1.;
				return (Real)pow( l , Confidence.value );
			};
			auto ProcessData = []( const Point< Real , Dim > &p , typename InputPointStreamInfo::DataType &d )
			{
				Real l = (Real)Length( d.template get<0>() );
				if( !l || !std::isfinite( l ) ) return (Real)-1.;
				d.template get<0>() /= l;
				return (Real)1.;
			};

			typename InputSampleDataFactory::VertexType zeroData = inputSampleDataFactory();
			typename FEMTreeInitializer< Dim , Real >::StreamInitializationData sid;
			if( Confidence.value>0 ) pointCount = FEMTreeInitializer< Dim , Real >::template Initialize< InputSampleDataType >( sid , tree.spaceRoot() , _pointStream , zeroData , Depth.value , *samples , *sampleData , true , tree.nodeAllocators[0] , tree.initializer() , ProcessDataWithConfidence );
			else                     pointCount = FEMTreeInitializer< Dim , Real >::template Initialize< InputSampleDataType >( sid , tree.spaceRoot() , _pointStream , zeroData , Depth.value , *samples , *sampleData , true , tree.nodeAllocators[0] , tree.initializer() , ProcessData );
		}

		unitCubeToModel = modelToUnitCube.inverse();
		delete pointStream;

		if( Verbose.set )
		{
			std::cout << "Input Points / Samples: " << pointCount << " / " << samples->size() << std::endl;
			std::cout << "# Read input into tree: " << profiler << std::endl;
		}
	}

	DenseNodeData< Real , Sigs > solution;
	{
		DenseNodeData< Real , Sigs > constraints;
		InterpolationInfo* iInfo = NULL;
		int solveDepth = Depth.value;

		tree.resetNodeIndices( 0 , std::make_tuple() );

		// Get the kernel density estimator
		{
			profiler.reset();
			density = tree.template setDensityEstimator< 1 , WEIGHT_DEGREE >( *samples , KernelDepth.value , SamplesPerNode.value );
			if( Verbose.set ) std::cout << "#   Got kernel density: " << profiler << std::endl;
		}

		// Transform the Hermite samples into a vector field
		{
			profiler.reset();
			normalInfo = new SparseNodeData< Point< Real , Dim > , NormalSigs >();
			std::function< bool ( InputSampleDataType , Point< Real , Dim >& ) > ConversionFunction = []( InputSampleDataType in , Point< Real , Dim > &out )
			{
				Point< Real , Dim > n = in.template get<0>();
				Real l = (Real)Length( n );
				// It is possible that the samples have non-zero normals but there are two co-located samples with negative normals...
				if( !l ) return false;
				out = n / l;
				return true;
			};
			std::function< bool ( InputSampleDataType , Point< Real , Dim >& , Real & ) > ConversionAndBiasFunction = []( InputSampleDataType in , Point< Real , Dim > &out , Real &bias )
			{
				Point< Real , Dim > n = in.template get<0>();
				Real l = (Real)Length( n );
				// It is possible that the samples have non-zero normals but there are two co-located samples with negative normals...
				if( !l ) return false;
				out = n / l;
				bias = (Real)( log( l ) * ConfidenceBias.value / log( 1<<(Dim-1) ) );
				return true;
			};
			if( ConfidenceBias.value>0 ) *normalInfo = tree.setInterpolatedDataField( NormalSigs() , *samples , *sampleData , density , BaseDepth.value , Depth.value , (Real)LowDepthCutOff.value , pointDepthAndWeight , ConversionAndBiasFunction );
			else                         *normalInfo = tree.setInterpolatedDataField( NormalSigs() , *samples , *sampleData , density , BaseDepth.value , Depth.value , (Real)LowDepthCutOff.value , pointDepthAndWeight , ConversionFunction );
			ThreadPool::Parallel_for( 0 , normalInfo->size() , [&]( unsigned int , size_t i ){ (*normalInfo)[i] *= (Real)-1.; } );
			if( Verbose.set )
			{
				std::cout << "#     Got normal field: " << profiler << std::endl;
				std::cout << "Point depth / Point weight / Estimated measure: " << pointDepthAndWeight.value()[0] << " / " << pointDepthAndWeight.value()[1] << " / " << pointCount*pointDepthAndWeight.value()[1] << std::endl;
			}
		}

		// Get the geometry designators indicating if the space node are interior to, exterior to, or contain the envelope boundary
#ifdef NOLIB
		if( Envelope.set )
		{
			profiler.reset();
			{
				// Make the octree complete up to the base depth
				FEMTreeInitializer< Dim , Real >::Initialize( tree.spaceRoot() , BaseDepth.value , []( int , int[] ){ return true; } , tree.nodeAllocators.size() ?  tree.nodeAllocators[0] : NULL , tree.initializer() );

				// Read in the envelope geometry
				std::vector< Point< Real , Dim > > vertices;
				std::vector< SimplexIndex< Dim-1 , node_index_type > > simplices;
				{
					std::vector< typename PositionFactory< Real , Dim >::VertexType > _vertices;
					std::vector< std::vector< int > > polygons;
					std::vector< std::string > comments;
					int file_type;
					PLY::ReadPolygons( Envelope.value , PositionFactory< Real , Dim >() , _vertices , polygons , file_type , comments );
					vertices.resize( _vertices.size() );
					for( int i=0 ; i<vertices.size() ; i++ ) vertices[i] = modelToUnitCube * _vertices[i];
					simplices.resize( polygons.size() );
					for( int i=0 ; i<polygons.size() ; i++ )
						if( polygons[i].size()!=Dim ) ERROR_OUT( "Not a simplex" );
						else for( int j=0 ; j<Dim ; j++ ) simplices[i][j] = polygons[i][j];
				}
				// Get the interior/boundary/exterior designators
				geometryNodeDesignators = FEMTreeInitializer< Dim , Real >::template GetGeometryNodeDesignators( &tree.spaceRoot() , vertices , simplices , BaseDepth.value , EnvelopeDepth.value , tree.nodeAllocators , tree.initializer() );

				// Make nodes in the support of the vector field @{ExactDepth} interior
				if( !NoDirichletErode.set )
				{
					// What to do if we find a node in the support of the vector field
					auto SetScratchFlag = [&]( FEMTreeNode *node )
					{
						if( node )
						{
							while( node->depth()>BaseDepth.value ) node = node->parent;
							node->nodeData.setScratchFlag( true );
						}
					};

					std::function< void ( FEMTreeNode * ) > PropagateToLeaves = [&]( const FEMTreeNode *node )
					{
						geometryNodeDesignators[ node ] = GeometryNodeType::INTERIOR;
						if( node->children ) for( int c=0 ; c<(1<<Dim) ; c++ ) PropagateToLeaves( node->children+c );
					};

					// Flags indicating if a node contains a non-zero vector field coefficient
					std::vector< bool > isVectorFieldElement( tree.nodeCount() , false );

					// Get the set of base nodes
					std::vector< FEMTreeNode * > baseNodes;
					auto nodeFunctor = [&]( FEMTreeNode *node )
					{
						if( node->depth()==BaseDepth.value ) baseNodes.push_back( node );
						return node->depth()<BaseDepth.value;
					};
					tree.spaceRoot().processNodes( nodeFunctor );

					std::vector< node_index_type > vectorFieldElementCounts( baseNodes.size() );
					for( int i=0 ; i<vectorFieldElementCounts.size() ; i++ ) vectorFieldElementCounts[i] = 0;

					// In parallel, iterate over the base nodes and mark the nodes containing non-zero vector field coefficients
					ThreadPool::Parallel_for( 0 , baseNodes.size() , [&]( unsigned int t , size_t  i )
					{
						auto nodeFunctor = [&]( FEMTreeNode *node )
						{
							Point< Real , Dim > *n = (*normalInfo)( node );
							if( n && Point< Real , Dim >::SquareNorm( *n ) ) isVectorFieldElement[ node->nodeData.nodeIndex ] = true , vectorFieldElementCounts[i]++;
						};
						baseNodes[i]->processNodes( nodeFunctor );
					} );
					size_t vectorFieldElementCount = 0;
					for( int i=0 ; i<vectorFieldElementCounts.size() ; i++ ) vectorFieldElementCount += vectorFieldElementCounts[i];

					// Get the subset of nodes containing non-zero vector field coefficients and disable the "scratch" flag
					std::vector< FEMTreeNode * > vectorFieldElements;
					vectorFieldElements.reserve( vectorFieldElementCount );
					{
						std::vector< std::vector< FEMTreeNode * > > _vectorFieldElements( baseNodes.size() );
						for( int i=0 ; i<_vectorFieldElements.size() ; i++ ) _vectorFieldElements[i].reserve( vectorFieldElementCounts[i] );
						ThreadPool::Parallel_for( 0 , baseNodes.size() , [&]( unsigned int t , size_t  i )
						{
							auto nodeFunctor = [&]( FEMTreeNode *node )
							{
								if( isVectorFieldElement[ node->nodeData.nodeIndex ] ) _vectorFieldElements[i].push_back( node );
								node->nodeData.setScratchFlag( false );
							};
							baseNodes[i]->processNodes( nodeFunctor );
						} );
						for( int i=0 ; i<_vectorFieldElements.size() ; i++ ) vectorFieldElements.insert( vectorFieldElements.end() , _vectorFieldElements[i].begin() , _vectorFieldElements[i].end() );
					}

					// Set the scratch flag for the base nodes on which the vector field is supported
#ifdef SHOW_WARNINGS
#pragma message( "[WARNING] In principal, we should unlock finite elements whose support overlaps the vector field" )
#endif // SHOW_WARNINGS
					tree.template processNeighboringLeaves< -BSplineSupportSizes< NORMAL_DEGREE >::SupportStart , BSplineSupportSizes< NORMAL_DEGREE >::SupportEnd >( &vectorFieldElements[0] , vectorFieldElements.size() , SetScratchFlag , false );

					// Set sub-trees rooted at interior nodes @ ExactDepth to interior
					ThreadPool::Parallel_for( 0 , baseNodes.size() , [&]( unsigned int , size_t  i ){ if( baseNodes[i]->nodeData.getScratchFlag() ) PropagateToLeaves( baseNodes[i] ); } );

					// Adjust the coarser node designators in case exterior nodes have become boundary.
					ThreadPool::Parallel_for( 0 , baseNodes.size() , [&]( unsigned int , size_t  i ){ FEMTreeInitializer< Dim , Real >::PullGeometryNodeDesignatorsFromFiner( baseNodes[i] , geometryNodeDesignators ); } );
					FEMTreeInitializer< Dim , Real >::PullGeometryNodeDesignatorsFromFiner( &tree.spaceRoot() , geometryNodeDesignators , BaseDepth.value );
				}
			}
			if( Verbose.set ) std::cout << "#               Initialized envelope constraints: " << profiler << std::endl;
		}
#endif
		if( !Density.set ) delete density , density = NULL;
		if( !needNormalData && !needAuxData ) delete sampleData , sampleData = NULL;

		// Add the interpolation constraints
		if( PointWeight.value>0 )
		{
			profiler.reset();
			if( ExactInterpolation.set ) iInfo = FEMTree< Dim , Real >::template       InitializeExactPointInterpolationInfo< Real , 0 > ( tree , *samples , ConstraintDual< Dim , Real >( targetValue , (Real)PointWeight.value * pointDepthAndWeight.value()[1] ) , SystemDual< Dim , Real >( (Real)PointWeight.value * pointDepthAndWeight.value()[1] ) , true , false );
			else                         iInfo = FEMTree< Dim , Real >::template InitializeApproximatePointInterpolationInfo< Real , 0 > ( tree , *samples , ConstraintDual< Dim , Real >( targetValue , (Real)PointWeight.value * pointDepthAndWeight.value()[1] ) , SystemDual< Dim , Real >( (Real)PointWeight.value * pointDepthAndWeight.value()[1] ) , true , Depth.value , 1 );
			if( Verbose.set ) std::cout <<  "#Initialized point interpolation constraints: " << profiler << std::endl;
		}

		// Trim the tree and prepare for multigrid
		{
			profiler.reset();
			constexpr int MAX_DEGREE = NORMAL_DEGREE > Degrees::Max() ? NORMAL_DEGREE : Degrees::Max();
			typename FEMTree< Dim , Real >::template HasNormalDataFunctor< NormalSigs > hasNormalDataFunctor( *normalInfo );
			auto hasDataFunctor = [&]( const FEMTreeNode *node ){ return hasNormalDataFunctor( node ); };
			auto addNodeFunctor = [&]( int d , const int off[Dim] ){ return d<=FullDepth.value; };
			if( geometryNodeDesignators.size() ) tree.template finalizeForMultigridWithDirichlet< MAX_DEGREE , Degrees::Max() >( BaseDepth.value , addNodeFunctor , hasDataFunctor , [&]( const FEMTreeNode *node ){ return node->nodeData.nodeIndex<(Index)geometryNodeDesignators.size() && geometryNodeDesignators[node]==GeometryNodeType::EXTERIOR; } , std::make_tuple( iInfo ) , std::make_tuple( normalInfo , density , &geometryNodeDesignators ) );
			else                                 tree.template finalizeForMultigrid< MAX_DEGREE , Degrees::Max() >( BaseDepth.value , addNodeFunctor , hasDataFunctor , std::make_tuple( iInfo ) , std::make_tuple( normalInfo , density ) );
#ifdef NOLIB
			if( geometryNodeDesignators.size() && EnvelopeGrid.set )
			{
				FEMTreeInitializer< Dim , Real >::PushGeometryNodeDesignatorsToFiner( &tree.spaceRoot() , geometryNodeDesignators );
				FEMTreeInitializer< Dim , Real >::PullGeometryNodeDesignatorsFromFiner( &tree.spaceRoot() , geometryNodeDesignators );

				auto WriteEnvelopeGrid = [&]( bool showFinest )
				{
					int res = 0;
					DenseNodeData< Real , IsotropicUIntPack< Dim , FEMTrivialSignature > > coefficients = tree.initDenseNodeData( IsotropicUIntPack< Dim , FEMTrivialSignature >() );
					auto nodeFunctor = [&]( const FEMTreeNode *n )
					{
						if( n->nodeData.nodeIndex!=-1 && ( ( showFinest && !n->children ) || ( !showFinest && geometryNodeDesignators[n->parent]==GeometryNodeType::BOUNDARY ) ) )
						{
#if 0 // Randomize the colors
							auto Value = []( double v , double eps ){ return (Real)( v + Random< double >() * 2. * eps - eps ); };
							// Show the octree structure
							if     ( geometryNodeDesignators[n]==GeometryNodeType::INTERIOR ) coefficients[n] = Value(  0.75 , 0.25 );
							else if( geometryNodeDesignators[n]==GeometryNodeType::EXTERIOR ) coefficients[n] = Value( -0.75 , 0.25 );

#else
							if     ( geometryNodeDesignators[n]==GeometryNodeType::INTERIOR ) coefficients[n] = (Real) 1.;
							else if( geometryNodeDesignators[n]==GeometryNodeType::EXTERIOR ) coefficients[n] = (Real)-1.;
#endif
						}
					};
					tree.spaceRoot().processNodes( nodeFunctor );
					Pointer( Real ) values = tree.template regularGridEvaluate< true >( coefficients , res , -1 , false );
					XForm< Real , Dim+1 > voxelToUnitCube = XForm< Real , Dim+1 >::Identity();
					for( int d=0 ; d<Dim ; d++ ) voxelToUnitCube( d , d ) = (Real)( 1. / res ) , voxelToUnitCube( Dim , d ) = (Real)( 0.5 / res );

					WriteGrid< Real , DEFAULT_DIMENSION >( EnvelopeGrid.value , values , res , unitCubeToModel * voxelToUnitCube , Verbose.set );
					DeletePointer( values );
				};

				WriteEnvelopeGrid( true );
			}
#endif
			if( Verbose.set ) std::cout << "#       Finalized tree: " << profiler << std::endl;
		}

		// Add the FEM constraints
		{
			profiler.reset();
			constraints = tree.initDenseNodeData( Sigs() );
			typename FEMIntegrator::template Constraint< Sigs , IsotropicUIntPack< Dim , 1 > , NormalSigs , IsotropicUIntPack< Dim , 0 > , Dim > F;
			unsigned int derivatives2[Dim];
			for( int d=0 ; d<Dim ; d++ ) derivatives2[d] = 0;
			typedef IsotropicUIntPack< Dim , 1 > Derivatives1;
			typedef IsotropicUIntPack< Dim , 0 > Derivatives2;
			for( int d=0 ; d<Dim ; d++ )
			{
				unsigned int derivatives1[Dim];
				for( int dd=0 ; dd<Dim ; dd++ ) derivatives1[dd] = dd==d ?  1 : 0;
				F.weights[d][ TensorDerivatives< Derivatives1 >::Index( derivatives1 ) ][ TensorDerivatives< Derivatives2 >::Index( derivatives2 ) ] = 1;
			}
			tree.addFEMConstraints( F , *normalInfo , constraints , solveDepth );
			if( Verbose.set ) std::cout << "#  Set FEM constraints: " << profiler << std::endl;
		}

		// Free up the normal info
		delete normalInfo , normalInfo = NULL;

		// Add the interpolation constraints
		if( PointWeight.value>0 )
		{
			profiler.reset();
			tree.addInterpolationConstraints( constraints , solveDepth , std::make_tuple( iInfo ) );
			if( Verbose.set ) std::cout << "#Set point constraints: " << profiler << std::endl;
		}

		if( Verbose.set )
		{
			std::cout << "All Nodes / Active Nodes / Ghost Nodes / Dirichlet Supported Nodes: " << tree.allNodes() << " / " << tree.activeNodes() << " / " << tree.ghostNodes() << " / " << tree.dirichletElements() << std::endl;
			std::cout << "Memory Usage: " << float( MemoryInfo::Usage())/(1<<20) << " MB" << std::endl;
		}
		
		// Solve the linear system
		{
			profiler.reset();
			typename FEMTree< Dim , Real >::SolverInfo sInfo;
			sInfo.cgDepth = 0 , sInfo.cascadic = true , sInfo.vCycles = 1 , sInfo.iters = Iters.value , sInfo.cgAccuracy = CGSolverAccuracy.value , sInfo.verbose = Verbose.set , sInfo.showResidual = ShowResidual.set , sInfo.showGlobalResidual = SHOW_GLOBAL_RESIDUAL_NONE , sInfo.sliceBlockSize = 1;
			sInfo.baseVCycles = BaseVCycles.value;
			typename FEMIntegrator::template System< Sigs , IsotropicUIntPack< Dim , 1 > > F( { 0. , 1. } );
			solution = tree.solveSystem( Sigs() , F , constraints , BaseDepth.value , SolveDepth.value , sInfo , std::make_tuple( iInfo ) );
			if( Verbose.set ) std::cout << "# Linear system solved: " << profiler << std::endl;
			if( iInfo ) delete iInfo , iInfo = NULL;
		}
	}

	{
		profiler.reset();
		double valueSum = 0 , weightSum = 0;
		typename FEMTree< Dim , Real >::template MultiThreadedEvaluator< Sigs , 0 > evaluator( &tree , solution );
		std::vector< double > valueSums( ThreadPool::NumThreads() , 0 ) , weightSums( ThreadPool::NumThreads() , 0 );
		ThreadPool::Parallel_for( 0 , samples->size() , [&]( unsigned int thread , size_t j )
		{
			ProjectiveData< Point< Real , Dim > , Real >& sample = (*samples)[j].sample;
			Real w = sample.weight;
			if( w>0 ) weightSums[thread] += w , valueSums[thread] += evaluator.values( sample.data / sample.weight , thread , (*samples)[j].node )[0] * w;
		}
		);
		for( size_t t=0 ; t<valueSums.size() ; t++ ) valueSum += valueSums[t] , weightSum += weightSums[t];
		isoValue = (Real)( valueSum / weightSum );
		if( !needNormalData && !needAuxData ) delete samples , samples = NULL;
		if( Verbose.set )
		{
			std::cout << "Got average: " << profiler << std::endl;
			std::cout << "Iso-Value: " << isoValue << " = " << valueSum << " / " << weightSum << std::endl;
		}
	}
#ifdef NOLIB
	if( Tree.set )
	{
		FILE* fp = fopen( Tree.value , "wb" );
		if( !fp ) ERROR_OUT( "Failed to open file for writing: " , Tree.value );
		FileStream fs(fp);
		FEMTree< Dim , Real >::WriteParameter( fs );
		DenseNodeData< Real , Sigs >::WriteSignatures( fs );
		tree.write( fs , modelToUnitCube , false );
		solution.write( fs );
		if( sampleData )
		{
			SparseNodeData< ProjectiveData< InputSampleDataType , Real > , IsotropicUIntPack< Dim , DataSig > > _sampleData = tree.template setExtrapolatedDataField< DataSig , false >( *samples , *sampleData , (DensityEstimator*)NULL );
			auto nodeFunctor = [&]( const RegularTreeNode< Dim , FEMTreeNodeData , depth_and_offset_type > *n )
			{
				ProjectiveData< InputSampleDataType , Real >* clr = _sampleData( n );
				if( clr ) (*clr) *= (Real)pow( DataX.value , tree.depth( n ) );
			};
			tree.tree().processNodes( nodeFunctor );
			_sampleData.write( fs );
		}
		if( density ) density->write( fs );
		fclose( fp );
	}

	if( Grid.set )
	{
		int res = 0;
		profiler.reset();
		Pointer( Real ) values = tree.template regularGridEvaluate< true >( solution , res , -1 , PrimalGrid.set );
		if( Verbose.set ) std::cout << "Got grid: " << profiler << std::endl;
		XForm< Real , Dim+1 > voxelToUnitCube = XForm< Real , Dim+1 >::Identity();
		if( PrimalGrid.set ) for( int d=0 ; d<Dim ; d++ ) voxelToUnitCube( d , d ) = (Real)( 1. / (res-1) );
		else                 for( int d=0 ; d<Dim ; d++ ) voxelToUnitCube( d , d ) = (Real)( 1. / res ) , voxelToUnitCube( Dim , d ) = (Real)( 0.5 / res );
		WriteGrid< Real , DEFAULT_DIMENSION >( Grid.value , values , res , unitCubeToModel * voxelToUnitCube , Verbose.set );
		DeletePointer( values );
	}
	if( Out.set )
	{
		if( Normals.value )
		{
			if( Density.set )
			{
				typedef Factory< Real , PositionFactory< Real , Dim > , NormalFactory< Real , Dim > , ValueFactory< Real > , AuxDataFactory > VertexFactory;
				VertexFactory vertexFactory( PositionFactory< Real , Dim >() , NormalFactory< Real , Dim >() , ValueFactory< Real >() , auxDataFactory );
				if( Normals.value==NORMALS_SAMPLES )
				{
					auto SetVertex = []( typename VertexFactory::VertexType &v , Point< Real , Dim > p , Point< Real , Dim > g , Real w , InputSampleDataType d ){ v.template get<0>() = p , v.template get<1>() = d.template get<0>() , v.template get<2>() = w , v.template get<3>() = d.template get<1>(); };
					ExtractMesh( UIntPack< FEMSigs ... >() , tree , solution , isoValue , samples , sampleData , density , vertexFactory , inputSampleDataFactory() , SetVertex , comments , unitCubeToModel );
				}
				else if( Normals.value==NORMALS_GRADIENTS )
				{
					auto SetVertex = []( typename VertexFactory::VertexType &v , Point< Real , Dim > p , Point< Real , Dim > g , Real w , InputSampleDataType d ){ v.template get<0>() = p , v.template get<1>() = -g/(1<<Depth.value) , v.template get<2>() = w , v.template get<3>() = d.template get<1>(); };
					ExtractMesh( UIntPack< FEMSigs ... >() , tree , solution , isoValue , samples , sampleData , density , vertexFactory , inputSampleDataFactory() , SetVertex , comments , unitCubeToModel );
				}
			}
			else
			{
				typedef Factory< Real , PositionFactory< Real , Dim > , NormalFactory< Real , Dim > , AuxDataFactory > VertexFactory;
				VertexFactory vertexFactory( PositionFactory< Real , Dim >() , NormalFactory< Real , Dim >() , auxDataFactory );
				if( Normals.value==NORMALS_SAMPLES )
				{
					auto SetVertex = []( typename VertexFactory::VertexType &v , Point< Real , Dim > p , Point< Real , Dim > g , Real w , InputSampleDataType d ){ v.template get<0>() = p                                                 , v.template get<1>() = d.template get<0>() , v.template get<2>() = d.template get<1>(); };
					ExtractMesh( UIntPack< FEMSigs ... >() , tree , solution , isoValue , samples , sampleData , density , vertexFactory , inputSampleDataFactory() , SetVertex , comments , unitCubeToModel );
				}
				else if( Normals.value==NORMALS_GRADIENTS )
				{
					auto SetVertex = []( typename VertexFactory::VertexType &v , Point< Real , Dim > p , Point< Real , Dim > g , Real w , InputSampleDataType d ){ v.template get<0>() = p                                                 , v.template get<1>() = -g/(1<<Depth.value) , v.template get<2>() = d.template get<1>(); };
					ExtractMesh( UIntPack< FEMSigs ... >() , tree , solution , isoValue , samples , sampleData , density , vertexFactory , inputSampleDataFactory() , SetVertex , comments , unitCubeToModel );
				}
			}
		}
		else
		{
			if( Density.set )
			{
				typedef Factory< Real , PositionFactory< Real , Dim > , ValueFactory< Real > , AuxDataFactory > VertexFactory;
				VertexFactory vertexFactory( PositionFactory< Real , Dim >() , ValueFactory< Real >() , auxDataFactory );
				auto SetVertex = []( typename VertexFactory::VertexType &v , Point< Real , Dim > p , Point< Real , Dim > g , Real w , InputSampleDataType d ){ v.template get<0>() = p , v.template get<1>() = w , v.template get<2>() = d.template get<1>(); };
				ExtractMesh( UIntPack< FEMSigs ... >() , tree , solution , isoValue , samples , sampleData , density , vertexFactory , inputSampleDataFactory() , SetVertex , comments , unitCubeToModel );
			}
			else
			{
				typedef Factory< Real , PositionFactory< Real , Dim > , AuxDataFactory > VertexFactory;
				VertexFactory vertexFactory( PositionFactory< Real , Dim >() , auxDataFactory );
				auto SetVertex = []( typename VertexFactory::VertexType &v , Point< Real , Dim > p , Point< Real , Dim > g , Real w , InputSampleDataType d ){ v.template get<0>() = p , v.template get<1>() = d.template get<1>(); };
				ExtractMesh( UIntPack< FEMSigs ... >() , tree , solution , isoValue , samples , sampleData , density , vertexFactory , inputSampleDataFactory() , SetVertex , comments , unitCubeToModel );
			}
		}
		if( sampleData ){ delete sampleData ; sampleData = NULL; }
	}
#else
    // gradient normals and density
    // typedef Factory< Real , PositionFactory< Real , Dim > , NormalFactory< Real , Dim > , ValueFactory< Real > , AuxDataFactory > VertexFactory;
    VertexFactory vertexFactory( PositionFactory< Real , Dim >() , NormalFactory< Real , Dim >() , ValueFactory< Real >() , auxDataFactory );
    auto SetVertex = []( typename VertexFactory::VertexType &v , Point< Real , Dim > p , Point< Real , Dim > g , Real w , InputSampleDataType d ){ v.template get<0>() = p , v.template get<1>() = -g/(1<<Depth.value) , v.template get<2>() = w , v.template get<3>() = d.template get<1>(); };
    ExtractMesh( UIntPack< FEMSigs ... >() , tree , solution , isoValue , samples , sampleData , density , vertexFactory , inputSampleDataFactory() , SetVertex , comments , unitCubeToModel );
#endif
	if( density ) delete density , density = NULL;
	if( Verbose.set ) std::cout << "#          Total Solve: " << Time()-startTime << " (s), " << MemoryInfo::PeakMemoryUsageMB() << " (MB)" << std::endl;
}
};

template< class Real,unsigned int Dim,typename Index>
void Reconstruct(std::vector<VectorTypeUnion<Real,Point<Real,Dim>,VectorTypeUnion<Real,Point<Real,Dim>,Point<Real,Dim>>>> v){
	static const int Degree = DEFAULT_FEM_DEGREE;
	static const BoundaryType BType = DEFAULT_FEM_BOUNDARY;
	typedef IsotropicUIntPack< DEFAULT_DIMENSION , FEMDegreeAndBType< Degree , BType >::Signature > FEMSigs;
	PoissonReconLib<Real,node_index_type,VertexFactory::RGBColorFactory< Real>> lib(FEMSigs());
    MeshOutputDataStream<typename PoissonReconLib<Real,node_index_type,VertexFactory::RGBColorFactory< Real>>::VertexFactory::VertexType,Index> output();
    lib.Execute(FEMSigs() , VertexFactory::RGBColorFactory< Real >(),v);

    MeshOutputDataStream<typename PoissonReconLib<Real,node_index_type,VertexFactory::RGBColorFactory< Real>,5U,5U,5U>::VertexFactory::VertexType,Index> output();
}