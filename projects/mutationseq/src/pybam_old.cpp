#include <iostream>
#include <iomanip>
#include <vector>
#include <queue>

#include <boost/foreach.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/incremental_components.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/pending/disjoint_sets.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/lexical_cast.hpp>

#include <boost/python/class.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/list.hpp>
#include <boost/python/tuple.hpp>

#include "bamtools/src/api/BamReader.h"
#include "bamtools/src/utils/bamtools_pileup_engine.h"
#include "bamtools/src/utils/bamtools_fasta.h"
#include "bamtools/src/utils/bamtools_utilities.h"

using namespace std;
using namespace boost;
using namespace BamTools;


bool CreatePileupTuple(const PileupPosition& pileupData, python::tuple& tpl, int coverage, bool rmdups, int readQualThreshold)
{
	int ntData[5][6] = {{0}};
	int ambiguous = 0;
	int insertionCount = 0;
	int deletionCount = 0;
	int tempQual = 0;
	int tempMQ = 0;
	float indelReads = 0;
	float allCoverage = 0;

	int aData[2] = {0};
	int cData[2] = {0};
	int gData[2] = {0};
	int tData[2] = {0};
	int allData[2] = {0};
		
	for (vector<PileupAlignment>::const_iterator pileupIter = pileupData.PileupAlignments.begin(); pileupIter != pileupData.PileupAlignments.end(); ++pileupIter)
	{
		const PileupAlignment& pa = (*pileupIter);
		const BamAlignment& ba = pa.Alignment;
		
		// remove duplicates and vendor failed reads
		if ((rmdups && ba.IsDuplicate())|| ba.IsFailedQC())
		{
			continue;
		}
		
		// remove the reads with qual<threshold
		if (ba.MapQuality < readQualThreshold)
		{
			continue;
		}

		//TODO: remove this, this is for mut-174
		//if (ba.Qualities.at(pa.PositionInAlignment) - 33 < 10)
		//	continue;

		// adjacent insertions and deletions
		for (vector<CigarOp>::const_iterator opIter = ba.CigarData.begin(); opIter != ba.CigarData.end(); opIter++)
		{
			if (opIter->Type == 'I')
			{
				insertionCount++;
			}
			else if (opIter->Type == 'D')
			{
				deletionCount++;
			}
		}

		//increment if alignment has an indel 
		allCoverage++;
		for (vector<CigarOp>::const_iterator opIter = ba.CigarData.begin(); opIter != ba.CigarData.end(); opIter++)
		{
			if (opIter->Type == 'I' or opIter->Type == 'D')
			{
				indelReads++;
				break;
			}
		}


		if (pa.IsCurrentDeletion)
		{
			continue;
		}
		
		char base = toupper(ba.QueryBases.at(pa.PositionInAlignment));
		
		if (base == 'N')
		{
			ambiguous++;
			continue;
		}
		
		int baseIdx;
		switch (base)
		{
			case 'A': baseIdx = 0; break;
			case 'C': baseIdx = 1; break;
			case 'G': baseIdx = 2; break;
			case 'T': baseIdx = 3; break;
			default: throw runtime_error("unrecognized base " + string(1, base));
		}
		
		// count
		ntData[baseIdx][0]++;
		ntData[4][0]++;
		
		// base quality, phred quality is scaled by 33
		tempQual = ba.Qualities.at(pa.PositionInAlignment) - 33;
		ntData[baseIdx][1] += tempQual; // * tempQual;
		ntData[4][1] += tempQual; //* tempQual;

		// mapping quality
		tempMQ = ba.MapQuality;
		ntData[baseIdx][2] += tempMQ ; //* tempMQ;
		ntData[4][2] += tempMQ ; //* tempMQ;
	
		// distance
		ntData[baseIdx][3] += (ba.IsReverseStrand()) ? ba.Length - pa.PositionInAlignment - 1 : pa.PositionInAlignment;
		ntData[4][3] += (ba.IsReverseStrand()) ? ba.Length - pa.PositionInAlignment - 1 : pa.PositionInAlignment;
		
		// direction
		ntData[baseIdx][4] += (ba.IsReverseStrand()) ? 1 : 0;
		ntData[4][4] += (ba.IsReverseStrand()) ? 1 : 0;

		//If reverse strand :
		switch(baseIdx){
			case 0 : ba.IsReverseStrand()? aData[1]++:aData[0]++; break;
			case 1 : ba.IsReverseStrand()? cData[1]++:cData[0]++; break;
			case 2 : ba.IsReverseStrand()? gData[1]++:gData[0]++; break;
			case 3 : ba.IsReverseStrand()? tData[1]++:tData[0]++; break;
		}
		ba.IsReverseStrand()?allData[0]++:allData[1]++;
	}
	
	// ignore positions with low or zero coverage 
	if (ntData[4][0] < coverage || ntData[4][0] == 0)
	{
		return false;
	}

	// Identify major base
	int majorBaseIdx = 0;
	for (int baseIdx = 0; baseIdx < 4; baseIdx++)
	{
		if (ntData[baseIdx][0] > ntData[majorBaseIdx][0])
		{
			majorBaseIdx = baseIdx;
		}
	}
	
	// Identify minor base, initialize to base that is not the major base
	int minorBaseIdx = (majorBaseIdx + 1) % 4;
	for (int baseIdx = 0; baseIdx < 4; baseIdx++)
	{
		if (ntData[baseIdx][0] > ntData[minorBaseIdx][0] && baseIdx != majorBaseIdx)
		{
			minorBaseIdx = baseIdx;
		}
	}
	
	// Calculate entropy
	double depth = (double)pileupData.PileupAlignments.size(); // NOTE: COPIED FROM JEFF, THIS MAY BE A WRONG
	double entropy = 0.0;
	for (int baseIdx = 0; baseIdx < 4; baseIdx++)
	{
		double probability = (double)ntData[baseIdx][0] / depth;
		if (probability != 0)
		{
			entropy -= (log(probability) * probability);
		}
	}
	
	// Interface is 1-based, bamtools is 0-based
	int position = pileupData.Position + 1;
	python::tuple atpl = python::make_tuple(aData[0], aData[1]);
	python::tuple ctpl = python::make_tuple(cData[0], cData[1]);
	python::tuple gtpl = python::make_tuple(gData[0], gData[1]);
	python::tuple ttpl = python::make_tuple(tData[0], tData[1]);
	python::tuple alltpl = python::make_tuple(allData[0],allData[1]);

	tpl = python::make_tuple(position,
							  python::make_tuple(ntData[0][0], ntData[0][1], ntData[0][2], ntData[0][3], ntData[0][4], atpl),
							  python::make_tuple(ntData[1][0], ntData[1][1], ntData[1][2], ntData[1][3], ntData[1][4], ctpl),
							  python::make_tuple(ntData[2][0], ntData[2][1], ntData[2][2], ntData[2][3], ntData[2][4], gtpl),
							  python::make_tuple(ntData[3][0], ntData[3][1], ntData[3][2], ntData[3][3], ntData[3][4], ttpl),
							  python::make_tuple(ntData[4][0], ntData[4][1], ntData[4][2], ntData[4][3], ntData[4][4], alltpl),
							  majorBaseIdx,
							  minorBaseIdx,
							  ambiguous,
							  insertionCount,
							  entropy,
							  deletionCount,
							  (indelReads/allCoverage),
							  pileupData.RefId
							  );
	// return success
	return true;
}


struct PileupQueue : PileupVisitor
{
	PileupQueue(int refId=-1, int start=-1, int stop=-1, int c=4, bool r=true, int rq = 10) : RefId(refId), StartPosition(start), StopPosition(stop), Coverage(c), Rmdups(r), ReadQualThreshold(rq)
	{
	}
	
	void Visit(const PileupPosition& pileupData)
	{
		// Reset if we ended up on the wrong chromosome
		if (pileupData.RefId != -1 && pileupData.RefId != RefId)
		{
			return;
		}
		
		// Dont store tuples out of the region of interest
		if (pileupData.Position < StartPosition || (pileupData.Position >= StopPosition && StopPosition != -1))
		{
			return;
		}
		
		python::tuple tpl;
		if(CreatePileupTuple(pileupData, tpl, Coverage, Rmdups, ReadQualThreshold))
		{
			Pileups.push(tpl);
		}
	}
	
	void Clear()
	{
		std::queue<python::tuple> empty;
		swap(Pileups, empty);
	}
	
	std::queue<python::tuple> Pileups;
	int RefId;
	int StartPosition;
	int StopPosition;
	int Coverage;
	bool Rmdups;
	int ReadQualThreshold;
};


class PyPileup
{
public:
	PyPileup(int c=4, bool r=true, int rq = 10) : m_PileupEngine(0), m_PileupQueue(0), RefId(-1), StartPosition(-1), StopPosition(-1), Coverage(c), Rmdups(r), ReadQualThreshold(rq)
	{
	}

	~PyPileup()
	{
		delete m_PileupEngine;
		delete m_PileupQueue;
	}
	
	void Open(const string& bamFilename)
	{
		if (!m_BamReader.Open(bamFilename))
		{
			throw runtime_error("unable to open bam file " + bamFilename);
		}
		
		if (!m_BamReader.LocateIndex())
		{
			throw runtime_error("unable to open index for bam file " + bamFilename);
		}
		
		RefNames = python::list();
		for (RefVector::const_iterator refDataIter = m_BamReader.GetReferenceData().begin(); refDataIter != m_BamReader.GetReferenceData().end(); refDataIter++)
		{
			// refName
			const string refName = refDataIter->RefName;

			// refId
			int refId = m_BamReader.GetReferenceID(refName);

			// keep both refId and its corresponding refName for a chromosome
			python::list temp_list = python::list();
			temp_list.append(refName);
			temp_list.append(refId);

			// creat a list of [refName, refId]
			RefNames.append(temp_list);
		}
		
		SamHeader = m_BamReader.GetHeaderText();
		RestartPileupEngine();
	}
	
	void Rewind()
	{
		m_BamReader.Rewind();
		
		RestartPileupEngine();
	}
	
	void SetChromosome(const string& refName)
	{
		int refId = m_BamReader.GetReferenceID(refName);
		
		if (refId < 0)
		{
			throw runtime_error("invalid ref name " + refName);
		}

		RefId = refId;
		StartPosition = -1;
		StopPosition  = -1;
		
		// set the region to a whole chromosome
	    if(!m_BamReader.SetRegion(BamRegion(refId, 0, refId+1, 1)))
	    	throw runtime_error("failed to set the region to chromosome " + refId);

		RestartPileupEngine();
	}
	
	void SetRegion(const string& refName, int start, int stop)
	{
		int refId = m_BamReader.GetReferenceID(refName);
		
		if (refId < 0)
		{
			throw runtime_error("invalid ref name " + refName);
		}
		
		RefId = refId;

		// Interface is 1-based, bamtools is 0-based
		StartPosition = start - 1;

		// BamRegion is half_open interval, no need for 1 to 0 based conversion
		StopPosition  = stop;
		
		//set the region to the chromosome:start-stop
		if (!m_BamReader.SetRegion(BamRegion(refId, StartPosition, refId, StopPosition)))
			throw runtime_error("failed to set the region to chromosome " + lexical_cast<string>(refId) + ":" + lexical_cast<string>(start) + "-" + lexical_cast<string>(stop));

		RestartPileupEngine();
	}
	
	python::object Next()
	{
		if (m_PileupEngine == 0)
		{
			throw runtime_error("next called before open");
		}
		
		if (!m_PileupQueue->Pileups.empty())
		{
			return PopPileup();
		}


		BamAlignment ba;
		while (m_BamReader.GetNextAlignment(ba)) //&& (ba.RefID == RefId || RefId == -1)
		{
			m_PileupEngine->AddAlignment(ba);
			if (!m_PileupQueue->Pileups.empty())
			{
				return PopPileup();
			}
		}
		m_PileupEngine->Flush();

		if (!m_PileupQueue->Pileups.empty())
		{
			return PopPileup();
		}
		
		return python::object();
	}
	
	python::list RefNames;
	std::string SamHeader;

private:
	void RestartPileupEngine()
	{
		delete m_PileupEngine;
		m_PileupEngine = new PileupEngine();
		
		delete m_PileupQueue;

		m_PileupQueue = new PileupQueue(RefId, StartPosition, StopPosition, Coverage, Rmdups, ReadQualThreshold);
		
		m_PileupEngine->AddVisitor(m_PileupQueue);
	}
	
	python::object PopPileup()
	{
		python::tuple tpl;
		swap(tpl, m_PileupQueue->Pileups.front());
		m_PileupQueue->Pileups.pop();

		return tpl;
	}

	BamReader m_BamReader;
	PileupEngine* m_PileupEngine;
	PileupQueue* m_PileupQueue;

	int RefId;
	int StartPosition;
	int StopPosition;

	// to specify the coverage of a position by which the tuples are filter
	int Coverage;

	// flag to remove/keep duplicates. Used for the deepseq data.
	bool Rmdups;

	//flag to filter out reads with low quality
	int ReadQualThreshold;

};

class PyFasta
{
public:
	PyFasta() : m_IsOpen(false)
	{
	}
	
	~PyFasta()
	{
		if (m_IsOpen)
		{
			m_Fasta.Close();
		}
	}
	
	void Open(const string& fastaFilename)
	{
		if (!Utilities::FileExists(fastaFilename))
		{
			throw runtime_error("invalid fasta file " + fastaFilename);
		}
		
		string indexFilename = fastaFilename + ".fai";
		
		if (!Utilities::FileExists(indexFilename))
		{
			throw runtime_error("index file " + indexFilename + " not found");
		}
		
		if (!m_Fasta.Open(fastaFilename, indexFilename))
		{
			throw runtime_error("unable to open fasta file " + fastaFilename);
		}
		
		// keep both refId and its corresponding refName as well as its length for all chromosomes in the reference 
		vector<string> referenceNames = m_Fasta.GetReferenceNames();
		m_RefLengths = m_Fasta.GetReferenceLengths();
		
		m_RefNames = python::list();
		m_RefLengthsList = python::list();
		for (int refId = 0; refId < referenceNames.size(); refId++)
		{
			python::list temp_list1 = python::list();
			temp_list1.append(referenceNames[refId]);
			temp_list1.append(refId);
			
			python::list temp_list2 = python::list();
			temp_list2.append(referenceNames[refId]);
			temp_list2.append(m_RefLengths[refId]);

			// create a list of [refName, refId]
			m_RefNames.append(temp_list1);
			
			// create a list of [refName, refLength]
			m_RefLengthsList.append(temp_list2);
		}
		
		m_IsOpen = true;
	}

	//python::object CreateFastaTuple(const string& refName, int position)
	python::object CreateFastaTuple(int refId, int position)
	{
		// Interface is 1-based, bamtools is 0-based
		position -= 1;
		
		if (!m_IsOpen)
		{
			throw runtime_error("get called before open");
		}
		
		// reference base
		char referenceBase = 'N';
		if (!m_Fasta.GetBase(refId, position, referenceBase))
		{
			throw runtime_error("unable to get base at " + lexical_cast<string>(refId) + ":" + lexical_cast<string>(position+1));
		}
		referenceBase = toupper(referenceBase);
		
		// get reference sequence for a window of length 500
		string referenceSeq;
		referenceSeq = GetReferenceSequenceByBase(refId, position, 500);

		// gc content
		double gc;
		gc = GcContent(position, referenceSeq);

		// entropy
		double ent;
		ent = Entropy(referenceSeq);

		// forward homopolymer
		int homoForward;
		homoForward = HomoForward(position, refId);

		// backward homopolymer
		int homoBackward;
		homoBackward = HomoBackward(position, refId);

		return python::make_tuple(referenceBase, homoForward, homoBackward, gc, ent);
	}
	
	// this is only used in the API
	string GetReferenceSequence(int refId, int position, int windowLength = 500) 
	{
		int start, stop;
		int refLength = m_RefLengths[refId];
		string referenceSeq;

		if ((position - (windowLength / 2)) < 0)
		{
			start = 0;            //start of the window
			stop = windowLength;  //stop  of the window
		}
		else if ((position + (windowLength / 2)) > refLength)
		{
			start = refLength - windowLength;
			stop = refLength;
		}
		else
		{
			start = (position - (windowLength / 2));
			stop = (position + (windowLength / 2));
		}

		// NOTE: Fasta::GetSequence is super slow, GetReferenceSequenceByBase is an alternative
		if (!m_Fasta.GetSequence(refId, start, stop, referenceSeq))
		{
			throw runtime_error("unable to get reference sequence at " + lexical_cast<string>(position));
		}
		return referenceSeq;
	}

	// an implementation based on the Fasta::GetBase
	string GetReferenceSequenceByBase(int refId, int position, int windowLength = 500)
	{
		string referenceSeq;
		char tempRefSeq[windowLength];
		char currentBase;
		int refLength = m_RefLengths[refId];

		// find the start of the window and set the position to that
		if ((position - (windowLength / 2)) < 0)
			position = 0;
		else if ((position + (windowLength / 2)) > refLength)
			position = refLength - windowLength;
		else
			position = (position - (windowLength / 2));

		// repeat Fasta::GetBase for all the positions in the window
		for(int i=0; i < windowLength; i++)
		{
			if (!m_Fasta.GetBase(refId, position+i, currentBase))
			{
				throw runtime_error("unable to get base at:" + lexical_cast<string>(refId) + ":" + lexical_cast<string>(position+i));
			}
			tempRefSeq[i] = toupper(currentBase);
		}

		tempRefSeq[windowLength] = '\0';
		referenceSeq = tempRefSeq;
		return referenceSeq;
	}

	 // get reference nucleotide at chromosomeId:position
	char GetReferenceBase(int refId, int position)
	{
		// Interface is 1-based, bamtools is 0-based
		position -= 1;
		char referenceBase = 'N';
		if (!m_Fasta.GetBase(refId, position, referenceBase))
		{
			throw runtime_error("unable to get base at " + lexical_cast<string>(refId) + ":" + lexical_cast<string>(position+1));
	    }
	    referenceBase = toupper(referenceBase);

		return referenceBase;
	}

	//unordered_map<string,int> m_RefNameId;
	python::list m_RefNames;
	python::list m_RefLengthsList;
	
private:
	double GcContent(int position, string &referenceSeq)
	{
		int counts[5] = {0};
		// compute the count of each nucleotide
		for (int i = 0; i < referenceSeq.length(); i++)
		{
			switch (referenceSeq[i])
			{
				case 'A': counts[0]++; break;
				case 'C': counts[1]++; break;
				case 'G': counts[2]++; break;
				case 'T': counts[3]++; break;
				default : counts[4]++; break;
			}
		}

		// NOTE: COPIED FROM JEFF, THIS MAY BE WRONG
		double gc, total;
		gc = counts[1] + counts[2];
		total = counts[0] + counts[1] + counts[2] + counts[3];
		if (total == 0)
		{
			return 0.0;
		}
		return gc / total;
	}

	double Entropy(string &referenceSeq)
	{
		int counts[5] = {0};
		// compute the count of each nucleotide
		for (int i = 0; i < referenceSeq.length(); i++)
		{
			switch (referenceSeq[i])
			{
				case 'A': counts[0]++; break;
				case 'C': counts[1]++; break;
				case 'G': counts[2]++; break;
				case 'T': counts[3]++; break;
				default : counts[4]++; break;
			}
		}

		// NOTE: COPIED FROM JEFF, THIS MAY BE WRONG
		double entropy, proba, total;
	    int i;

	    entropy = 0.0;
		total = counts[0] + counts[1] + counts[2] + counts[3];
		if (total == 0)
		{
			return 0;
		}

		for (int i = 0; i < 4; i++)
		{
			proba = (double)counts[i] / total;
			if (proba != 0.0)
			{
				entropy -= ((log(proba)/log(4)) * proba);
			}
		}
		return entropy;
	}

	// NOTE: COPIED FROM JEFF, THIS MAY BE WRONG
	int HomoForward(int position, int refId)
	{
		position += 1; // start after the position of interest
	    char currentBase = 'N';
	    char nextBase    = 'N';
	    int count = 0;

	    if (!m_Fasta.GetBase(refId, position, currentBase))
		{
			throw runtime_error("unable to get base at:" + lexical_cast<string>(position));
		}

	    if (!m_Fasta.GetBase(refId, position+1, nextBase))
		{
			throw runtime_error("unable to get base at:" + lexical_cast<string>(position+1));
		}
		currentBase = toupper(currentBase);
		nextBase = toupper(nextBase);

	    while (currentBase == nextBase)
	    {
	        count++;
	        position++;
	        currentBase = nextBase;
	        if (!m_Fasta.GetBase(refId, position+1, nextBase))
			{
				throw runtime_error("unable to get base at:" + lexical_cast<string>(position+1));
			}
			nextBase = toupper(nextBase);
	    }
	    return count;
	}

	// NOTE: COPIED FROM JEFF, THIS MAY BE WRONG
	int HomoBackward(int position, int refId)
	{
		if (position < 2)
			return 0;

		position -= 1; // start before the position of interest
	    char currentBase = 'N';
	    char previousBase    = 'N';
	    int count = 0;

	    if (!m_Fasta.GetBase(refId, position, currentBase))
		{
			throw runtime_error("unable to get base at:" + lexical_cast<string>(position));
		}

	    if (!m_Fasta.GetBase(refId, position-1, previousBase))
		{
			throw runtime_error("unable to get base at:" + lexical_cast<string>(position-1));
		}
		currentBase = toupper(currentBase);
		previousBase = toupper(previousBase);

	    while (currentBase == previousBase && position > 2)
	    {
	        count++;
	        position--;
	        currentBase = previousBase;
	        if (!m_Fasta.GetBase(refId, position-1, previousBase))
			{
				throw runtime_error("unable to get base at:" + lexical_cast<string>(position-1));
			}
			previousBase = toupper(previousBase);
	    }
	    return count;
	}

	Fasta m_Fasta;
	bool m_IsOpen;
	vector<int> m_RefLengths;
};

BOOST_PYTHON_MODULE(pybam)
{
	using namespace python;
	
	class_<PyPileup>("pileup", init<int, bool, int>())
		.def_readonly("refnames", &PyPileup::RefNames)
		.def_readonly("samheader", &PyPileup::SamHeader)
		.def("open", &PyPileup::Open)
		.def("rewind", &PyPileup::Rewind)
		.def("set_region", &PyPileup::SetChromosome)
		.def("set_region", &PyPileup::SetRegion)
		.def("get_tuple", &PyPileup::Next)
	;
	
	class_<PyFasta>("fasta", init<>())
		.def_readonly("refnames", &PyFasta::m_RefNames)
		.def_readonly("reflenghts", &PyFasta::m_RefLengthsList)
		.def("open", &PyFasta::Open)
		.def("get_tuple", &PyFasta::CreateFastaTuple)
		.def("get_base", &PyFasta::GetReferenceBase)
		.def("get_sequence", &PyFasta::GetReferenceSequence)
		.def("get_sequence_base", &PyFasta::GetReferenceSequenceByBase)
	;
}

