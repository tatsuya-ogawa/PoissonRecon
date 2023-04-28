#include "public/PoissonRecon.h"
#include "PoissonReconLib.h"

// #ifdef USE_DOUBLE
// 	typedef double Real;
// #else // !USE_DOUBLE
// 	typedef float  Real;
// #endif // USE_DOUBLE


template<typename Real>
PoissonReconOutput<Real> ExecutePoissonRecon(PoissonReconInput<Real> input){
	static const int Dim = 3;
	typedef InputOrientedPointStreamInfo< Real , Dim , typename VertexFactory::RGBColorFactory< Real>::VertexType > InputPointStreamInfo;
	// The type of the input sample
	typedef typename InputPointStreamInfo::PointAndDataType InputSampleType;
	Point<Real,3U> p(0.f,0.f,0.f);
	VectorTypeUnion<Real,Point<Real,3U>,Point<Real,3U>> nc(p,p);
	VectorTypeUnion<Real,Point<Real,3U>,VectorTypeUnion<Real,Point<Real,3U>,Point<Real,3U>>> all(p,nc);
	InputSampleType a(p,nc);
    std::vector<InputSampleType> v;
    // auto pointStream = new MemoryInputDataStream< InputSampleType >( v.size() , &v[0] );
    auto vectorStream = new VectorInputDataStream< InputSampleType >( v );
    return PoissonReconOutput<Real>();
}
