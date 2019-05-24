#include "MVAReader.hh"
using std::vector;
using std::string;
namespace QQbarAnalysis 
{
	MVAReader::MVAReader(string weightpath)
	{
		myWeightPath = weightpath;
	}
	float MVAReader::GetTrustTag(MVAInputParameters & parameters)
	{
		string minuspath = myWeightPath+"weights/TMVAClassification_MLP.weights.xml";
		float minustag = GetTag(parameters, minuspath);
		std::cout << "Tag: " << minustag << "\n";
		return minustag;
	}
	vector<float> MVAReader::GetAllTags(MVAInputParameters & parameters)
	{
		vector<float> result;
		string minuspath = myWeightPath+"weights_minus/TMVAClassification_MLP.weights.xml";
		float minustag = GetTag(parameters, minuspath);
		string zeropath = myWeightPath+"weights_zero/TMVAClassification_MLP.weights.xml";
		float zerotag = GetTag(parameters, zeropath);
		string pluspath = myWeightPath+"weights_plus/TMVAClassification_MLP.weights.xml";
		float plustag = GetTag(parameters, pluspath);
		std::cout << "-1tag: " << minustag << " 0tag: " << zerotag << " +1tag: " << plustag << "\n";
		result.push_back(minustag);
		result.push_back(zerotag);
		result.push_back(plustag);
		return result;
	}

	float MVAReader::GetTag(MVAInputParameters & parameters, std::string path)
	{
		float tag = -1;
		TMVA::Reader* reader = new TMVA::Reader("Silent");
		reader->AddVariable( "btag", &parameters.btag );
		reader->AddVariable( "VtxCharge",  &parameters.VtxCharge);
		reader->AddVariable( "nVtxParticles",  &parameters.nVtxParticles);
		reader->AddVariable( "VtxMomentum",  &parameters.VtxMomentum);
		reader->AddVariable( "VtxDistance",  &parameters.VtxDistance);
		reader->AddVariable( "VtxMinHits",  &parameters.VtxMinHits);
		reader->AddVariable( "VtxMinOffset",  &parameters.VtxMinOffset);
		reader->AddVariable( "VtxChargeBalance",  &parameters.VtxChargeBalance);
		
		reader->AddVariable( "JetMaxOffset",  &parameters.JetMaxOffset);
		reader->AddVariable( "JetMaxOffsetMomentum",  &parameters.JetMaxOffsetMomentum);
		reader->AddVariable( "JetMaxOffsetAngle",  &parameters.JetMaxOffsetAngle);
		reader->AddVariable( "JetMaxOffsetCharge",  &parameters.JetMaxOffsetCharge);
		reader->BookMVA( "myMLP", path.c_str() );
		tag = reader->EvaluateMVA( "myMLP" );
		delete reader;
		return tag;
	}
} /* QQbarAnalysis */
