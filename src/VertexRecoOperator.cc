#include "VertexRecoOperator.hh"
using std::vector;
using std::string;
using EVENT::LCCollection;
using EVENT::Vertex;
using EVENT::ReconstructedParticle;
namespace QQbarAnalysis 
{
	VertexRecoOperator:: VertexRecoOperator() 
	{
		myAngleCut = 0.151;
		myPrecisionCut = 1.0;
		ip[0] = 0.0;
		ip[1] = 0.0;
		ip[2] = 0.0;
	}
	VertexRecoOperator:: VertexRecoOperator(float angle, Vertex * primary)
	{
		myAngleCut = angle;
		myPrecisionCut = 1.5;
		myPrimary = primary;
		ip[0] = 0.0;
		ip[1] = 0.0;
		ip[2] = 0.0;
 	}	
	vector< VertexTag * > * VertexRecoOperator::Compare(LCCollection * reconstructed, LCCollection * mc)
	{
		int number = reconstructed->getNumberOfElements();
		vector< VertexTag * > * result = new vector< VertexTag * >();
		for (int i = 0; i < number; i++) 
		{
			Vertex * vertex = dynamic_cast< Vertex * >( reconstructed->getElementAt(i) ) ;
			VertexTag * tagged = Tag(vertex, mc);
			if (tagged) 
			{
				result->push_back(tagged);
			}
			else 
			{
				myUnknownVertexes.push_back(vertex);
			}
		}

		return result;
	}

	vector< VertexTag * > * VertexRecoOperator::CompareDirection(LCCollection * reconstructed, LCCollection * mc)
	{
		int number = reconstructed->getNumberOfElements();
		int mcnumber = mc->getNumberOfElements();
		vector< VertexTag * > * result = new vector< VertexTag * >();
		for (int i = 0; i < number; i++) 
		{
			Vertex * recovertex = dynamic_cast< Vertex * >( reconstructed->getElementAt(i) ) ;
			bool passed = false;
			for (int j = 0; j < mcnumber; j++) 
			{
				Vertex * mcvertex = dynamic_cast< Vertex * >( mc->getElementAt(j) ) ;
				passed = false;
				float angle = MathOperator::getAngle(mcvertex->getPosition(), recovertex->getPosition());
				std::cout << "Angle: " << angle  << '\n';
				if (angle < myAngleCut) 
				{
					std::cout << "Vertex tagged with pdg " << mcvertex->getParameters()[1] << " n-tracks: " << recovertex->getAssociatedParticle()->getParticles().size() << '\n';
					passed = true;
					VertexTag * vtag = new VertexTag(recovertex, mcvertex);
					vtag->SetTruthAngle(angle);
					result->push_back(vtag);
					break;
				}
			}
			if (!passed) 
			{
				myUnknownVertexes.push_back(recovertex);
			}
		}
		/*vector< VertexTag * > * refined = refineTagSelection(result);
		for (int i = 0; i < refined->size(); i++) 
		{
			Vertex * recovertex = refined->at(i)->GetVertex();
			Vertex * mcvertex = refined->at(i)->__GetMCVertex();
			vector< ReconstructedParticle * > * missed = getMissedTracks(recovertex->getAssociatedParticle()->getParticles(), mcvertex->getAssociatedParticle()->getParticles());
		}*/
		return result;
	}

	vector< Particle > * VertexRecoOperator::CompareTracks(LCCollection * reconstructed, LCCollection * mc, vector< ReconstructedParticle * > * recomissed)
	{
		vector< Particle > * missed = new vector< Particle >();
		std::cout << "\n-----------------\nTrack analysis:\n";
		int number = reconstructed->getNumberOfElements();
		int mcnumber = mc->getNumberOfElements();
		vector< ReconstructedParticle * > mctracks;
		vector< ReconstructedParticle * > recotracks;
		for (int i = 0; i < number; i++) 
		{
			Vertex * recovertex = dynamic_cast< Vertex * >( reconstructed->getElementAt(i) ) ;
			recotracks.reserve(recotracks.size()+recovertex->getAssociatedParticle()->getParticles().size());
			recotracks.insert(recotracks.end(),recovertex->getAssociatedParticle()->getParticles().begin(), recovertex->getAssociatedParticle()->getParticles().end());
		}
		std::cout << "All reco tracks number: " << recotracks.size() << '\n';

		for (int i = 0; i < mcnumber; i++) 
		{
			Vertex * mcvertex = dynamic_cast< Vertex * >( mc->getElementAt(i) ) ;
			
			vector< ReconstructedParticle * > * missedvertex = getMissedTracks(recotracks, mcvertex->getAssociatedParticle()->getParticles());
			//double * ip = MathOperator::toDoubleArray(myPrimary->getPosition(), 3);
			for (int j = 0; j < missedvertex->size(); j++) 
			{
				ReconstructedParticle * particle = missedvertex->at(j);
				if (recomissed) 
				{
					recomissed->push_back(particle);
				}
				//ReconstructedParticle * particle = mcvertex->getAssociatedParticle()->getParticles().at(j);
				vector< float > direction = MathOperator::getDirection(particle->getMomentum());
				double * sec = MathOperator::toDoubleArray(mcvertex->getPosition(), 3);
				vector< float > angles = MathOperator::getAngles(direction);
				float offset = MathOperator::getDistanceTo(ip, direction, sec);
				Particle missedParticle;
				missedParticle.SetOffset(offset);
				missedParticle.SetTheta(angles[1]);
				missedParticle.SetMomentum(particle->getMomentum());
				missedParticle.SetVertex(mcvertex);
				missed->push_back(missedParticle);
				std::cout << "INFO: MISSED TRACK OFFSET:  " << offset
					  << " MOMENTUM: " << MathOperator::getModule(missedParticle.GetMomentum())
					  << " THETA: " << missedParticle.GetTheta()
					  << '\n';
			}
		}
		std::cout << "Number of missed tracks: " << missed->size()  << '\n';
		return missed;

	}

	vector< VertexTag * > * VertexRecoOperator::refineTagSelection(vector< VertexTag * > * old)
	{
		vector< VertexTag * > bvtx;
		vector< VertexTag * > bbarvtx;
		vector< VertexTag * > * refined = new vector< VertexTag * >();
		for (int i = 0; i < old->size(); i++) 
		{
			if (old->at(i)->GetInitialPDG() > 0) 
			{
				bvtx.push_back(old->at(i));
			}
			else 
			{
				bbarvtx.push_back(old->at(i));
			}
		}
		refineChain(refined, bvtx); 
		refineChain(refined, bbarvtx); 
		delete old;
		return refined;
	}



	vector< VertexTag * > * VertexRecoOperator::refineChain(vector< VertexTag * > * refined, vector< VertexTag * > & bvtx)
	{
		if (bvtx.size() == 2) 
		{
			std::cout << "Checking:\n";
			Vertex * recovertex1 = bvtx[0]->GetVertex();
			Vertex * mcvertex1 = bvtx[0]->__GetMCVertex();
			float distance1 = MathOperator::getDistance(recovertex1->getPosition(), mcvertex1->getPosition());
			if (distance1 > myPrecisionCut) 
			{
				std::cout << "Distance does not coincide... Refining!\n";
				Vertex * recovertex2 = bvtx[1]->GetVertex();
				Vertex * mcvertex2 = bvtx[1]->__GetMCVertex();
				refined->push_back(new VertexTag(recovertex2, mcvertex1));
				refined->push_back(new VertexTag(recovertex1, mcvertex2));
			}
			else 
			{
				refined->push_back(bvtx[0]);
				refined->push_back(bvtx[1]);
			}
		}
		if (bvtx.size() == 1) 
		{
			refined->push_back(bvtx[0]);
		}
		return refined;
	}

	VertexTag * VertexRecoOperator::Tag(Vertex * vertex, LCCollection * mc)
	{
		int number = mc->getNumberOfElements();
		VertexTag * result = NULL ;
		float min = 1000.0;
		int candidateNumber = -1; 
		for (int i = 0; i < number; i++) 
		{
			Vertex * mcvertex = dynamic_cast< Vertex * >( mc->getElementAt(i) ) ;
			float distance = MathOperator::getDistance(vertex->getPosition(), mcvertex->getPosition());
			myDistances.push_back(distance);
			if (distance < min && !CheckTaken(mcvertex)) 
			{
				min = distance;
				candidateNumber = i;
			}
		}
		if (min < myPrecisionCut) 
		{
			Vertex * mcvertex = dynamic_cast< Vertex * >( mc->getElementAt(candidateNumber) ) ;
			myTakenVertices.push_back(GetID(mcvertex));

			std::cout << "Found "<< mcvertex->getParameters()[2] << " vertex in " << min << "mm radius and mcid " << GetID(mcvertex) << "\n";
			result = new VertexTag(vertex, mcvertex);
			SetMinimalDistance(result, mc);
		}
		return result;
	}
	void VertexRecoOperator::SetMinimalDistance(VertexTag * tag, LCCollection * mc)
	{
		int number = mc->getNumberOfElements();
		float min = 1000.0;
		Vertex * mcvertex = tag->GetVertex();
		for (int j = 0; j < number; j++) 
		{
			Vertex * mcvertex2 = dynamic_cast< Vertex * >( mc->getElementAt(j) ) ;
			if (((int)mcvertex2->getParameters()[1]) != tag->GetInitialPDG() || ((int)mcvertex2->getParameters()[2]) != tag->GetGeneration()) 
			{
				float distance = MathOperator::getDistance(mcvertex->getPosition(), mcvertex2->getPosition());
				if (distance < min) 
				{
					min = distance;
				}
			}
		}
		std::cout << "Minimal distance "<< min << "\n";
		tag->SetMinimalDistance(min);
	}
	vector< ReconstructedParticle * > * VertexRecoOperator::getMissedTracks(const vector< ReconstructedParticle * > & recotracks, const vector< ReconstructedParticle * > & mctracks)
	{
		//recotracks,mctracks
		std::cout << "Comparing the tracks!...\n";
		vector< ReconstructedParticle * > * missed = new vector< ReconstructedParticle * >();
		for (unsigned int i = 0; i < mctracks.size(); i++) 
		{
			bool found = false;
			ReconstructedParticle * mctrack = mctracks[i];
			for (unsigned int j = 0; j < recotracks.size(); j++) 
			{
				ReconstructedParticle * recotrack = recotracks[j];
				float angle = MathOperator::getAngleBtw(mctrack->getMomentum(), recotrack->getMomentum());
				float recomodule = MathOperator::getModule(recotrack->getMomentum());
				float mcmodule = MathOperator::getModule(mctrack->getMomentum());
				float ratio = (1 - recomodule/mcmodule > 0.0) ? 1 - recomodule/mcmodule : recomodule/mcmodule - 1.0;
				//std::cout <<  ratio  << '\n';
				
				if (angle < 0.01 && mctrack->getCharge() == recotrack->getCharge() && ratio < 0.005) 
				{
					//std::cout << "\tFound track with angle: " << angle  << " and ratio " << ratio << '\n';
					found = true;
					break;
				}
				if (recomodule > 20.0) 
				{
					if (angle < 0.1 
					    && mctrack->getCharge() == recotrack->getCharge() 
					    && ratio < 0.007) 
					{
						//std::cout << "\tFound track with angle: " << angle << " and ratio " << ratio << '\n';
						found = true;
						break;
					}
				}
			}
			if (!found) 
			{
				missed->push_back(mctrack);
			}
		}
		return missed;
	}
	bool VertexRecoOperator::CheckTaken(Vertex * mcVertex)
	{
		bool result = false;
		for (int i = 0; i < myTakenVertices.size(); i++) 
		{
			int id = GetID(mcVertex);
			if (myTakenVertices[i] == id) 
			{
				result = true;
				std::cout << "Vertex is taken!\n";
				break;
			}
		}
		return result;
	}
	vector< float >  VertexRecoOperator::GetDistances()
	{
		return myDistances;
	}
	vector< Vertex * > & VertexRecoOperator::GetUnknownVertexes()
	{
		return myUnknownVertexes;
	}
	int VertexRecoOperator::GetID(EVENT::Vertex * mcVertex)
	{
		return ((int)mcVertex->getParameters()[2])*((int)mcVertex->getParameters()[1]);
	}
	vector<float> VertexRecoOperator::getAngles(EVENT::Vertex * vertex)
	{
		ReconstructedParticle * particle = vertex->getAssociatedParticle();
		const float * position = vertex->getPosition();
		double convert[3];
		for (int i = 0; i < 3; i++) 
		{
			convert[i] = position[i];
		}
		vector<float> direction = MathOperator::getDirection(convert);//particle->getMomentum());
		return MathOperator::getAngles(direction);
	}
} /* QQbarAnalysis */
