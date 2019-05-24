#include "VertexTag.hh"
using std::vector;
using std::string;
using EVENT::Vertex;
namespace QQbarAnalysis 
{
	VertexTag:: VertexTag(Vertex * recoVertex, Vertex * mcVertex)
	{
		myRecoVertex = recoVertex;
		myMCVertex = mcVertex;
		myInitialPDG = myMCVertex->getParameters()[1];
		myStatus = UNKNOWN_TAG;
	}
	VertexTag:: VertexTag(Vertex * mcVertex)
	{
		myRecoVertex = NULL;
		myMCVertex = mcVertex;
		myStatus = UNKNOWN_TAG;
	}
	VertexTag:: VertexTag()
	{
		myRecoVertex = NULL;
		myMCVertex = NULL;
		myStatus = UNKNOWN_TAG;
	}
	void VertexTag::SetRecoVertex(EVENT::Vertex * vertex)
	{
		myRecoVertex = vertex;
	}
	void VertexTag::SetStatus(TagStatus status)
	{
		myStatus = status;
	}
	TagStatus VertexTag::GetStatus()
	{
		return myStatus;
	}
	Vertex * VertexTag::GetVertex()
	{
		return myRecoVertex;
	}
	Vertex * VertexTag::__GetMCVertex()
	{
		return myMCVertex;
	}
	int VertexTag::GetGeneration()
	{
		return myMCVertex->getParameters()[2];
	}
	int VertexTag::GetInitialPDG()
	{
		return myMCVertex->getParameters()[1];
	}
	float VertexTag::GetMinimalDistance()
	{
		return myMinimalDistanceMC;
	}
	void VertexTag::SetMinimalDistance(float d)
	{
		myMinimalDistanceMC = d;
	}
	
	float VertexTag::GetTruthAngle()
	{
		return myTruthAngle;
	}
	void VertexTag::SetTruthAngle(float theta)
	{
		myTruthAngle = theta;
	}
	int VertexTag::__GetMCTrackNumber()
	{
		return myMCVertex->getAssociatedParticle()->getParticles().size();
	}
	int VertexTag::GetTrackNumber()
	{
		return myRecoVertex->getAssociatedParticle()->getParticles().size();
	}


} /* QQbarAnalysis */
