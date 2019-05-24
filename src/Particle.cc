#include "Particle.hh"
using std::vector;
using std::string;
using EVENT::Vertex;
using EVENT::ReconstructedParticle;
using EVENT::MCParticle;
namespace QQbarAnalysis
{
	Particle::Particle ()
	{
		//myMomentum = NULL;
		myVertex = NULL;
		myMass = 0.0;
		myTruthAngle = -1;
		myVertexAngle = -1;
		myHasDecayed = -1;
		myBtag = -1.0;
		myHits[0] = -1;
		myHits[1] = -1;
		myHits[2] = -1;
		myHits[3] = -1;
		myMCParticle = NULL;
		myRecoParticle = NULL;
	}
	const double * Particle::GetMomentum() const
	{
		return myMomentum;
	}
	void Particle::SetMomentum(const double * momentum)
	{
		//myMomentum = new double[3];
		for (int i = 0; i < 3; i++) 
		{
			myMomentum[i] = momentum[i];
		}
	}
	Vertex * Particle::GetVertex()
	{
		return myVertex;
	}
	float Particle::GetOffset()
	{
		return myOffset;
	}
	void Particle::SetOffset(float offset)
	{
		myOffset = offset;
	}
	void Particle::SetVertex(Vertex * vertex)
	{
		myVertex = vertex;
	}
	float Particle::GetTheta()
	{
		return myTheta;
	}
	float Particle::GetMass()
	{
		return myMass;
	}
	void Particle::SetTheta(float theta)
	{
		myTheta = theta;
	}
	void Particle::SetMass(float mass)
	{
		myMass = mass;
	}
	int Particle::GetGeneration()
	{
		return myGeneration;
	}
	void Particle::SetGeneration(int value)
	{
		myGeneration = value;
	}
	int Particle::GetIsReco()
	{
		return myIsReco;
	}
	void Particle::SetIsReco(int value)
	{
		myIsReco = value;
	}
	const int * Particle::GetHits()
	{
		return myHits;
	}
	void Particle::SetHits(int * value)
	{
		myHits[0] = value[0];
		myHits[1] = value[1];
		myHits[2] = value[2];
		myHits[3] = value[3];
	}
	float Particle::GetChi2()
	{
		return myChi2;
	}
	void Particle::SetChi2(float value)
	{
		myChi2 = value;
	}
	float Particle::GetTruthAngle()
	{
		return myTruthAngle;
	}
	void Particle::SetTruthAngle(float value)
	{
		myTruthAngle = value;
	}
	int Particle::GetInteracted()
	{
		return myHasDecayed;
	}
	void Particle::SetInteracted(int value)
	{
		myHasDecayed = value;
	}
	float Particle::GetJetBtag()
	{
		return myBtag;
	}
	void Particle::SetJetBtag(float value)
	{
		myBtag = value;
	}
	MCParticle * Particle::GetMCParticle()
	{
		return myMCParticle;
	}
	void Particle::SetMCParticle(MCParticle * value)
	{
		myMCParticle = value;
	}
	float Particle::GetVertexAngle()
	{
		return myVertexAngle;
	}
	void Particle::SetVertexAngle(float value)
	{
		myVertexAngle = value;
	}
	ReconstructedParticle * Particle::GetRecoParticle()
	{
		return myRecoParticle;
	}
	void Particle::SetRecoParticle(ReconstructedParticle * particle)
	{
		myRecoParticle = particle;
	}
	
	void Particle::Assign(ReconstructedParticle * initial)
	{
		//myMomentum = new double( *initial->getMomentum());
	}
}
