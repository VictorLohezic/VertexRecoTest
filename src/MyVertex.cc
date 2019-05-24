#include "MyVertex.hh"
using std::vector;
using EVENT::MCParticle;
using IMPL::VertexImpl;
namespace QQbarAnalysis 
{
	MyVertex:: MyVertex()
	{
	}
	
	void MyVertex::__SetMCParticles(const vector< MCParticle * > particles)
	{
		myMCParticles = particles;
	}
	vector< MCParticle * > & MyVertex::__GetMCParticles()
	{
		return myMCParticles;
	}
}
