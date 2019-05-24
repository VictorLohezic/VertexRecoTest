#include "ParticleOperator.hh"
using std::vector;
using std::string;
using EVENT::LCCollection;
using EVENT::MCParticle;
using EVENT::Vertex;
using EVENT::Track;
using EVENT::ReconstructedParticle;
using IMPL::ReconstructedParticleImpl;
using UTIL::LCRelationNavigator;
using EVENT::LCObject;
namespace QQbarAnalysis 
{
	ParticleOperator:: ParticleOperator()
	{
	}
	bool ParticleOperator::CompareParticles(const ReconstructedParticle * particle1, const ReconstructedParticle * particle2)
	{
		//std::cout << "Comparing particles with modules " << recomodule << " & " << mcmodule << "\n";
		if (particle1 == particle2) 
		{
			//std::cout << "Equal by pointer\n";
			return true;
		}
		if (particle1->getCharge() * particle2->getCharge() < 0.0) 
		{
			//std::cout << "Rejected by charge\n";
			return false;
		}
		float angle = MathOperator::getAngleBtw(particle1->getMomentum(), particle2->getMomentum());
		if (angle > 0.02) 
		{
			//std::cout << "Rejected by angle " << angle << "\n";
			return false;
		}
		float recomodule = MathOperator::getModule(particle2->getMomentum());
		float mcmodule = MathOperator::getModule(particle1->getMomentum());
		float ratio = (1 - recomodule/mcmodule > 0.0) ? 1 - recomodule/mcmodule : recomodule/mcmodule - 1.0;
		if (ratio > 0.04) 
		{
			//std::cout << "Rejected by ratio " << ratio << "\n";
			return false;
		}
		//std::cout << "Particles equal!\n";
		return true;
	}
	bool ParticleOperator::IsDublicate(const ReconstructedParticle * particle, const vector< ReconstructedParticle * > & data)
	{
		bool dublicate = false;
		for (unsigned int j = 0; j < data.size(); j++) 
		{
			if (CompareParticles(particle, data[j])) 
			{
				//std::cout << "Dublicate found!!!!\n";
				dublicate = true;
				break;
			}
		}
		return dublicate;
	}
	bool ParticleOperator::CompareParticles(MCParticle * particle1, ReconstructedParticle * particle2)
	{
		//std::cout << "Comparing particles with modules " << recomodule << " & " << mcmodule << "\n";
		if (particle1->getCharge() * particle2->getCharge() < 0.0) 
		{
			//std::cout << "Rejected by charge\n";
			return false;
		}
		float angle = MathOperator::getAngleBtw(particle1->getMomentum(), particle2->getMomentum());
		if (angle > 0.002) 
		{
			//std::cout << "Rejected by angle " << angle << "\n";
			return false;
		}
		float recomodule = MathOperator::getModule(particle2->getMomentum());
		float mcmodule = MathOperator::getModule(particle1->getMomentum());
		float ratio = (1 - recomodule/mcmodule > 0.0) ? 1 - recomodule/mcmodule : recomodule/mcmodule - 1.0;
		if (ratio > 0.004) 
		{
			//std::cout << "Rejected by ratio " << ratio << "\n";
			return false;
		}
		return true;
	}
	bool ParticleOperator::CheckForVertex(EVENT::ReconstructedParticle * particle1)
	{
		Vertex * vertex = particle1->getStartVertex ();
		if (vertex && !vertex->isPrimary()) 
		{
			return true;
		}
		return false;
	}
	bool ParticleOperator::IsDublicate(MCParticle * particle, vector< MCParticle * > & list)
	{
		for (unsigned int i = 0; i < list.size(); i++) 
		{
			if (particle == list[i]) 
			{
				return true;
			}
		}
		return false;
	}
	MCParticle * ParticleOperator::GetMCParticleTrackRel(const ReconstructedParticle * secondary, LCCollection * trackrel)
	{
		LCRelationNavigator navigator(trackrel);
		MCParticle * mcparticle = NULL;
		if (std::abs(secondary->getCharge()) < 0.1)
		{
			return mcparticle;
		}
		Track * reco = secondary->getTracks()[0];
		vector< LCObject * > obj = navigator.getRelatedToObjects(reco);
		vector< float > weights = navigator.getRelatedToWeights(reco);
		if (obj.size() < 1) 
		{
			std::cout << "ERROR: no Track truthlink for particle either\n";
			return mcparticle;
		}
		float maxweight = 0.5;
		for (int i = 0; i < obj.size(); i++) 
		{
			MCParticle * candidate = dynamic_cast< MCParticle * >(obj[i]);
			if (std::abs(candidate->getCharge()) < 0.09) 
			{
				std::cout << "WARNING: neutral truthlink for particle\n";
				continue;
			}
			if (weights[i] > maxweight) 
			{
				maxweight = weights[i];
				mcparticle = candidate;
			}
		}
		return mcparticle;

	}
	vector< MCParticle * > ParticleOperator::GetMCParticlesRel(const vector< ReconstructedParticle * > & secondaries, LCCollection * rel, LCCollection * trackrel)
	{
		LCRelationNavigator navigator(rel);
		vector< MCParticle * > result;
		for (unsigned int i = 0; i < secondaries.size(); i++) 
		{
			ReconstructedParticle * reco = secondaries[i];
			if (std::abs(reco->getCharge()) < 0.1) 
			{
				continue;
			}
			vector< LCObject * > obj = navigator.getRelatedToObjects(reco);
			vector< float > weights = navigator.getRelatedToWeights(reco);
			if (obj.size() < 1) 
			{
				std::cout << "ERROR: no PFO truthlink for particle here\n";
				MCParticle * candidate = GetMCParticleTrackRel(reco, trackrel);
				if (candidate) 
				{
					result.push_back(candidate);
				}
				continue;
			}

			float maxweight = 0.50;
			MCParticle * mcparticle = NULL;
			for (int j = 0; j < obj.size(); j++) 
			{
				MCParticle * candidate = dynamic_cast< MCParticle * >(obj[j]);
				if (std::abs(candidate->getCharge()) < 0.09) 
				{
					std::cout << "WARNING: neutral truthlink for particle\n";
					continue;
				}
				if (weights[j] > maxweight) 
				{
					maxweight = weights[j];
					mcparticle = candidate;
				}
			}
			if (mcparticle == NULL) 
			{
				std::cout << "ERROR: no charged truthlink for particle\n";
				continue;
			}
			result.push_back(mcparticle);
		}
		return result;
	}
	Track * ParticleOperator::GetTrackRel(MCParticle *  particle, LCCollection * rel)
	{
		if (!particle) 
		{
			return NULL;
		}
		LCRelationNavigator navigator(rel);
		Track * result = NULL;
		vector< LCObject * > obj = navigator.getRelatedFromObjects(particle);
		vector< float > weights = navigator.getRelatedFromWeights(particle);
		float maxweight = 0.0;
		for (int i = 0; i < obj.size(); i++)
		{
			Track * candidate = dynamic_cast< Track *>(obj[i]);
			if (weights[i] > maxweight)
			{
			        maxweight = weights[i];
			        result = candidate;
			}
		}
		return result;
	}
	ReconstructedParticle * ParticleOperator::GetRecoParticle(MCParticle *  particle, LCCollection * rel)
	{
		LCRelationNavigator navigator(rel);
		ReconstructedParticle * result = NULL;
		vector< LCObject * > obj = navigator.getRelatedFromObjects(particle);
		vector< float > weights = navigator.getRelatedFromWeights(particle);
		float maxweight = 0.0;
		for (int i = 0; i < obj.size(); i++) 
		{
			ReconstructedParticle * candidate = dynamic_cast< ReconstructedParticle * >(obj[i]);
			if (std::abs(candidate->getCharge()) < 0.09) 
			{
				std::cout << "WARNING: neutral truthlink for particle\n";
				continue;
			}
			if (weights[i] > maxweight) 
			{
				maxweight = weights[i];
				result = candidate;
			}
		}
		return result;
	}
	ReconstructedParticle * ParticleOperator::ReconstructParticle(EVENT::Track *  track)
	{
		float Bz = 3.5;
		float a = 3.0e-4;
		float omega = track->getOmega();
		float tanl = track->getTanLambda();
		float phi = track->getPhi();
		float pt = a * std::abs(Bz / omega);
		//float p = pt * std::sqrt(1 + tanl * tanl);
		
		double momentum[3];
		momentum[0] = pt * std::cos(phi);
		momentum[1] = pt * std::sin(phi);
		momentum[2] = pt * tanl;

		float charge = Bz / omega / std::abs(Bz / omega);
		ReconstructedParticleImpl * result = new ReconstructedParticleImpl();
		result->setCharge(charge);
		result->setMomentum(momentum);
		result->addTrack(track);
		return result;
	}
	float ParticleOperator::GetError(const ReconstructedParticle * particle)
	{
		float _aParameter = 0.005;
		float _bParameter = 0.01;
		if (!particle) 
		{
			std::cout << "The particle is null!\n";
			return 0.0;
		}
		float p = MathOperator::getModule(particle->getMomentum());
		vector<float> direction = MathOperator::getDirection(particle->getMomentum());
		vector<float> angles = MathOperator::getAngles(direction);
		float accuracy = sqrt(_aParameter*_aParameter + _bParameter*_bParameter /( p * p * pow(sin(angles[1]), 4.0/3.0)) );
		return accuracy;
	}
}
