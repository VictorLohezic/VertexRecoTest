
#include "Jet.hh"
using EVENT::ReconstructedParticle;
using std::vector;
namespace QQbarAnalysis
{
	Jet::Jet (float btag,float ctag, int number,const double * momentum)
	{
		myTrustTag = -1;
		myZeroTag = -1;
		myMinusTag = -1;
		myPlusTag = -1;
		myBTag = btag;
		myCTag = ctag;
		myNumber = number;
		myMomentum = momentum;
		myMCPDG = 0;
		myParticles = NULL;
		myRecoVertices = NULL;
		myTagAngle = -1;
	}
	Jet::Jet ()
	{
		myRecoVertices = NULL;
		myMCPDG = 0;
		myParticles = NULL;
		
	}
	int Jet::GetCharge()
	{
		int result = 0;
		for (unsigned int i = 0; i < myParticles->size(); i++) 
		{
			result+= myParticles->at(i)->getCharge();
		}
		return result;
	}
	void Jet::SetRecoVertices(std::vector<  EVENT::Vertex * > * vertices)
	{
		myRecoVertices = vertices;
	}
	std::vector<  EVENT::Vertex * > * Jet::GetRecoVertices()
	{
		return myRecoVertices;
	}
	void Jet::AddVertexTag(VertexTag * tag)
	{
		myVertexTags.push_back(tag);
	}
	const vector< VertexTag * > & Jet::GetVertexTags() const
	{
		return myVertexTags;
	}
	void Jet::SetChargeTags(float minustag, float zerotag, float plustag)
	{
		myMinusTag = minustag;
		myZeroTag = zerotag;
		myPlusTag = plustag;
	}
	void Jet::SetBTag(float value)
	{
		myBTag = value;
	}
	void Jet::SetCTag(float value)
	{
		myCTag = value;
	}
	float Jet::GetBTag()
	{
		return myBTag;
	}
	float Jet::GetCTag()
	{
		return myCTag;
	}
	float Jet::GetTagAngle()
	{
		return myTagAngle;
	}
	void Jet::SetTagAngle(float angle)
	{
		myTagAngle = angle;
	}
	int Jet::GetNumberOfVertices()
	{
		if (myRecoVertices) 
		{
			return myRecoVertices->size();
		}
		return myNumber;
	}
	int Jet::GetNumberOfVertexParticles()
	{
		int sum = -1;
		if (myRecoVertices) 
		{
			sum = 0;
			for (unsigned int i = 0; i < myRecoVertices->size(); i++) 
			{
				sum += myRecoVertices->at(i)->getAssociatedParticle()->getParticles().size();
			}
		}
		return sum;
	}
	float Jet::GetMinVtxProbability()
	{
		float minprob = 2.;
		if (myRecoVertices) 
		{
			for (unsigned int i = 0; i < myRecoVertices->size(); i++) 
			{
				float probability = myRecoVertices->at(i)->getProbability();
				if (probability < minprob) 
				{
					minprob = probability;
				}
			}
		}
		return minprob;
	}
	float Jet::GetMaxVtxChi2()
	{
		float maxchi = -1.;
		if (myRecoVertices) 
		{
			for (unsigned int i = 0; i < myRecoVertices->size(); i++) 
			{
				float chi = myRecoVertices->at(i)->getChi2();
				if (chi > maxchi) 
				{
					maxchi = chi;
				}
			}
		}
		return maxchi;
	}
	float Jet::GetZeroTag()
	{
		return myZeroTag;
	}
	float Jet::GetPlusTag()
	{
		return myPlusTag;
	}
	float Jet::GetMinusTag()
	{
		return myMinusTag;
	}
	int Jet::__GetGenNumberOfVertexParticles()
	{
		int sum = 0;
		for (int i = 0; i < myVertexTags.size(); i++) 
		{
			sum += myVertexTags[i]->__GetMCTrackNumber();
		}
		return sum;
	}
	int Jet::__GetGenCharge()
	{
		int sum = 0;
		for (int i = 0; i < myVertexTags.size(); i++) 
		{
			if (myVertexTags[i]->__GetMCVertex()->getParameters()[2] == 2) 
			{
				sum = myVertexTags[i]->__GetMCVertex()->getAssociatedParticle()->getCharge();
			}
		}
		return sum;
	}
	float Jet::GetHadronCharge()
	{
		float charge = -5.0;
		if (myRecoVertices) 
		{
			charge = 0.0;
			for (unsigned int i = 0; i < myRecoVertices->size(); i++) 
			{
				charge += myRecoVertices->at(i)->getAssociatedParticle()->getCharge();
			}
		}
		return charge;
	}
	float Jet::GetHadronMomentum()
	{
		float momentum = -1.0;
		if (myRecoVertices) 
		{
			momentum = 0.0;
			for (unsigned int i = 0; i < myRecoVertices->size(); i++)
			{
				momentum += MathOperator::getModule(myRecoVertices->at(i)->getAssociatedParticle()->getMomentum()); // CRUNCH!!!
			}
		}
		return momentum;
	}
	float Jet::GetHadronMass()
	{
		float mass = -1.0;
		if (myRecoVertices) 
		{
			 mass = 0.0;
			 for (int i = 0; i < myRecoVertices->size(); i++) 
			 {
			 	mass += myRecoVertices->at(i)->getAssociatedParticle()->getMass();
			 }
		}
		return mass;
	}
	float Jet::GetHadronCostheta()
	{
		float costheta = -2.0;
		if (myRecoVertices && myRecoVertices->size() > 0 && myRecoVertices->at(0)) 
		{
			//const float * position = myRecoVertices->at(0)->getPosition();
			//std::cout << position[2] << "\n";
			//const double * vertex = MathOperator::toDoubleArray(position,3);
			vector< float > direction = MathOperator::getDirection(myRecoVertices->at(0)->getAssociatedParticle()->getMomentum());
			costheta = std::cos(MathOperator::getAngles(direction)[1]);
		}
		return costheta;
	}
	float Jet::GetHadronDistance()
	{
		float distance = 1000.0;
		if (myRecoVertices && myRecoVertices->size() > 0) 
		{
			for (int i = 0; i < myRecoVertices->size(); i++) 
			{
				float current = MathOperator::getModule(myRecoVertices->at(i)->getPosition());
				if (current < distance) 
				{
					distance = current;
				}
			}
		}
		else 
		{
			distance = -1.0;
		}
		return distance;
	}
	const double * Jet::GetMomentum()
	{
		return myMomentum;
	}
	const vector< ReconstructedParticle * > * Jet::GetParticles() const
	{
		return myParticles;
	}
	void Jet::SetParticles(const vector< ReconstructedParticle * > & particles)
	{
		myParticles = new vector< ReconstructedParticle * >( particles);
	}
	float Jet::GetTrustTag()
	{
		return myTrustTag;
	}
	void Jet::SetTrustTag(float tag)
	{
		myTrustTag = tag;
	}
	int Jet::GetMCPDG()
	{
		return myMCPDG;
	}
	void Jet::SetMCPDG(int pdg)
	{
		myMCPDG = pdg;
	}

}
