#include "JetOperator.hh"
#include "TVector3.h"
using namespace lcio ;
using namespace marlin ;
using std::vector;
using std::string;
using EVENT::LCCollection;
using EVENT::MCParticle;
using EVENT::ReconstructedParticle;
namespace QQbarAnalysis 
{
	JetOperator:: JetOperator ( ) 
	{
		myAlgorithmName = "lcfiplus";
		myAngle = 0.6;
	}
	JetOperator:: JetOperator (float angle, string algorithm )
	{
	        myAlgorithmName = algorithm;
	        myAngle = angle;
	}

	
	vector<float> * JetOperator::GetBtags(LCCollection * jets)
	{
		vector<float> * tags = new vector<float>();
		int number = jets->getNumberOfElements();
		PIDHandler pidh(jets);
		int alid = pidh.getAlgorithmID("lcfiplus");
		for (int i = 0; i < number; i++) 
		{
			ReconstructedParticle * jet = dynamic_cast< ReconstructedParticle * >(jets->getElementAt(i));
			const ParticleID& pid = pidh.getParticleID(jet,alid);
			vector<float> params = pid.getParameters();
			float btag = params[pidh.getParameterIndex(alid,"BTag")];
			std::cout << "Jet energy: " << jet->getEnergy() << " b-tag: " << btag <<  '\n';

			tags->push_back(btag);
		}
		return tags;
	}
	vector< Jet * > * JetOperator::GetJets(LCCollection * jets, LCCollection * rel)
	{
		int number = jets->getNumberOfElements();
		vector< Jet * > * result = new vector< Jet * >();
		LCRelationNavigator navigator(rel);
		PIDHandler pidh(jets);
		int alid = pidh.getAlgorithmID(myAlgorithmName);
		for (int i = 0; i < number; i++) 
		{
			ReconstructedParticle * jet = dynamic_cast< ReconstructedParticle * >(jets->getElementAt(i));
			int nvtx = navigator.getRelatedToObjects(jet).size();

			const ParticleID& pid = pidh.getParticleID(jet,alid);
			vector<float> params = pid.getParameters();
			float btag = params[pidh.getParameterIndex(alid,"BTag")];
			float ctag = params[pidh.getParameterIndex(alid,"CTag")];
			std::cout << "Jet energy: " << jet->getEnergy() 
				  << " b-tag: " << btag 
				  << " c-tag: " << ctag 
				  << " # of vtx: " << nvtx 
				  <<  '\n';
			const ReconstructedParticleVec components = jet->getParticles();
			/*for (int j = 0; j < components.size(); j++) 
			{
				std::cout << "\tParticle energy: " << components[j]->getEnergy()
					  << "\tPx: " << components[j]->getMomentum()[0]
					  << "\tPy: " << components[j]->getMomentum()[1]
					  << "\tPz: " << components[j]->getMomentum()[2]
					  <<  '\n';
			}*/
			Jet * tag = new Jet(btag, ctag, nvtx, jet->getMomentum());
			if (nvtx > 0) 
			{
				vector< Vertex * > * vertices = convert(navigator.getRelatedToObjects(jet));
				tag->SetRecoVertices(vertices);
			}
			tag->SetParticles(components);
			result->push_back(tag);
		}
		
		return result;
	}
	void JetOperator::CompareDirection(vector< Jet *> * jets, LCCollection * mc)
	{
		if (!mc || mc->getNumberOfElements() == 0) 
		{
			return;
		}
		int mcnumber = mc->getNumberOfElements();
		/*for (int i = 0; i < jets->size(); i++) 
		{
			Jet * jet = jets->at(i);
			const vector< ReconstructedParticle * > * components = jet->GetParticles();
			for (int k = 0; k < components->size(); k++) 
			{
				ReconstructedParticle * leader = components->at(k);
				bool found = false;
				for (int j = 0; j < mcnumber; j++) 
				{
					MCParticle * mcparticle = dynamic_cast< MCParticle * >( mc->getElementAt(j) ) ;
					float angle = MathOperator::getAngle(mcparticle->getMomentum(), leader->getMomentum());
					std::cout << "Jet angle: " << angle << " particle " << mcparticle->getPDG() << '\n';
					if (angle < myAngle) 
					{
						int pdg = mcparticle->getPDG();
						jet->SetMCPDG(pdg);
						std::cout << "Jet with btag " << jet->GetBTag() << " truth-tagged with " << pdg << '\n';
						found = true;
						break;
					}
				}
				if (found) 
				{
					break;
				}
			}
		}*/
		for (int i = 0; i < mcnumber; i++) 
		{
			int number = -1;
			MCParticle * mcparticle = dynamic_cast< MCParticle * >( mc->getElementAt(i) ) ;
			for (int j = 0; j < jets->size(); j++) 
			{
				if (number == j) 
				{
					continue;
				}
				Jet * jet = jets->at(j);
				const vector< ReconstructedParticle * > * components = jet->GetParticles();
				bool found = false;
				for (int k = 0; k < components->size(); k++) 
				{
					ReconstructedParticle * leader = components->at(k);
					float angle = MathOperator::getAngle(mcparticle->getMomentum(), leader->getMomentum());
					//std::cout << "Jet angle: " << angle << " particle " << mcparticle->getPDG() << '\n';
					if (angle < myAngle) 
					{
						int pdg = mcparticle->getPDG();
						jet->SetMCPDG(pdg);
						jet->SetTagAngle( MathOperator::getAngle(mcparticle->getMomentum(), jet->GetMomentum()));
						std::cout << "Jet with btag " << jet->GetBTag() << " truth-tagged with " << pdg << '\n';
						found = true;
						number = j;
						break;
					}
				}
				if (found) 
				{
					break;
				}
			}
		}
	}
	std::vector< Vertex * > * JetOperator::convert(const std::vector< LCObject * > & objs)
	{
		std::vector< Vertex * > * result = new std::vector< Vertex * >();
		for (int i = 0; i < objs.size(); i++) 
		{
			result->push_back(dynamic_cast< Vertex * >(objs[i]));
		}
		return result;
	}
	
}
