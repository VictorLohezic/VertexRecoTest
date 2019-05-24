#include "ComparisonOperator.hh"
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
	ComparisonOperator:: ComparisonOperator(LCCollection * pfo, LCCollection * egprongs, LCCollection * trackrel)
	{
		myPFO = pfo;
		myEGProngs = egprongs;
		if (trackrel) 
		{
			myTrackRel = trackrel;
		}
		myAlgorithmName = "lcfiplus";
		myAngleCut = 0.4; // CHANGED
		ip[0] = 0.0;
		ip[1] = 0.0;
		ip[2] = 0.0;
		
	}
	void ComparisonOperator::TagJets(std::vector<RecoJet*> * jets, EVENT::LCCollection * mccol, vector<float> & angleTags)
	{
		int mcnumber = mccol->getNumberOfElements();
		vector<Vertex * > vertices;
		vector<Vertex * > bvertices;
		vector<Vertex * > bbarvertices;
		for (int i = 0; i < mcnumber; i++) 
		{
			Vertex * mcvertex = dynamic_cast< Vertex * >(mccol->getElementAt(i));
			vertices.push_back(mcvertex);
			if (std::abs(mcvertex->getParameters()[1]) >0) 
			{
				bvertices.push_back(mcvertex);
			}
			else 
			{
				bbarvertices.push_back(mcvertex);
			}
		}
		if (bvertices.size() > 0) 
		{
			tagJets(jets, bvertices[0], angleTags);
		}
		if (bbarvertices.size() > 0) 
		{
			tagJets(jets, bbarvertices[0], angleTags);
		}
	}
	void ComparisonOperator::tagJets(vector<RecoJet*> * jets, Vertex * vertex, vector<float> & angleTags)
	{
		int jetnumber = jets->size();
		ReconstructedParticle * mcparticle = vertex->getAssociatedParticle();
		float minangle = myAngleCut;
		int winner = -1;
		for (int i = 0; i < jetnumber; i++)
		{
			RecoJet * jet = jets->at(i);
			float angle = MathOperator::getAngleBtw(mcparticle->getMomentum(), jet->getMomentum());
			if (angle < minangle)
			{
				minangle = angle;
				winner = i;
			}//*/
		}
		if(winner > -1)
		{
			angleTags.push_back(minangle);
			jets->at(winner)->SetMCPDG(vertex->getParameters()[1]);
		}
	}
}
