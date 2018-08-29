/*
 *  Common.h
 *  linkexons
 *
 *  Created by Andrew McPherson on 15/09/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef COMMON_H_
#define COMMON_H_

#include <vector>
#include <functional>
#include <numeric>
#include <iostream>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace boost;

enum Strand
{
	PlusStrand = 0,
	MinusStrand = 1,
};

int OtherStrand(int strand);
int InterpretStrand(const string& strand);

struct RefStrand
{
	union
	{
		struct
		{
			unsigned referenceIndex : 31;
			unsigned strand : 1;
		};
		
		uint32_t id;
	};
};

int OtherReadEnd(int readEnd);

struct Region
{
	int start;
	int end;	
};

inline bool operator==(const Region& r1, const Region& r2)
{
	return r1.start == r2.start && r1.end == r2.end;
}

inline int Length(const Region& r)
{
	return r.end - r.start + 1;
}

inline bool Overlap(const Region& r1, const Region& r2)
{
	if (r1.end < r2.start || r1.start > r2.end)
	{
		return false;
	}
	else
	{
		return true;
	}
}

typedef vector<Region> RegionVec;
typedef vector<Region>::iterator RegionVecIter;
typedef vector<Region>::const_iterator RegionVecConstIter;

typedef vector<RegionVec> RegionTable;
typedef vector<RegionVec>::iterator RegionTableIter;
typedef vector<RegionVec>::const_iterator RegionTableConstIter;

typedef vector<string> StringVec;
typedef vector<string>::iterator StringVecIter;
typedef vector<string>::const_iterator StringVecConstIter;

typedef vector<StringVec> StringTable;
typedef vector<StringVec>::iterator StringTableIter;
typedef vector<StringVec>::const_iterator StringTableConstIter;

typedef vector<int> IntegerVec;
typedef vector<int>::iterator IntegerVecIter;
typedef vector<int>::const_iterator IntegerVecConstIter;

typedef vector<IntegerVec> IntegerTable;
typedef vector<IntegerVec>::iterator IntegerTableIter;
typedef vector<IntegerVec>::const_iterator IntegerTableConstIter;

typedef vector<float> FloatVec;
typedef vector<float>::iterator FloatVecIter;
typedef vector<float>::const_iterator FloatVecConstIter;

typedef vector<FloatVec> FloatTable;
typedef vector<FloatVec>::iterator FloatTableIter;
typedef vector<FloatVec>::const_iterator FloatTableConstIter;

typedef vector<double> DoubleVec;
typedef vector<double>::iterator DoubleVecIter;
typedef vector<double>::const_iterator DoubleVecConstIter;

typedef vector<DoubleVec> DoubleTable;
typedef vector<DoubleVec>::iterator DoubleTableIter;
typedef vector<DoubleVec>::const_iterator DoubleTableConstIter;

typedef pair<int,int> IntegerPair;

typedef vector<IntegerPair> IntegerPairVec;
typedef vector<IntegerPair>::iterator IntegerPairVecIter;
typedef vector<IntegerPair>::const_iterator IntegerPairVecConstIter;

typedef vector<IntegerPairVec> IntegerPairTable;
typedef vector<IntegerPairVec>::iterator IntegerPairTableIter;
typedef vector<IntegerPairVec>::const_iterator IntegerPairTableConstIter;

typedef unordered_set<int> IntegerSet;
typedef unordered_set<int>::iterator IntegerSetIter;
typedef unordered_set<int>::const_iterator IntegerSetConstIter;

typedef unordered_map<int,int> IntegerMap;
typedef unordered_map<int,int>::iterator IntegerMapIter;
typedef unordered_map<int,int>::const_iterator IntegerMapConstIter;

typedef unordered_map<int,IntegerVec> IntegerVecMap;
typedef unordered_map<int,IntegerVec>::iterator IntegerVecMapIter;
typedef unordered_map<int,IntegerVec>::const_iterator IntegerVecMapConstIter;

typedef unordered_map<IntegerPair,int> IntegerPairMap;
typedef unordered_map<IntegerPair,int>::iterator IntegerPairMapIter;
typedef unordered_map<IntegerPair,int>::const_iterator IntegerPairMapConstIter;

typedef unordered_map<int,double> DoubleMap;
typedef unordered_map<int,double>::iterator DoubleMapIter;
typedef unordered_map<int,double>::const_iterator DoubleMapConstIter;

typedef pair<IntegerVec,IntegerVec> IntegerVecPair;

typedef vector<IntegerVecPair> IntegerVecPairVec;
typedef vector<IntegerVecPair>::iterator IntegerVecPairVecIter;
typedef vector<IntegerVecPair>::const_iterator IntegerVecPairVecConstIter;

void ReverseComplement(string& sequence);

inline bool operator==(const IntegerVec& vec1, const IntegerVec& vec2)
{
	return vec1.size() == vec2.size() && equal(vec1.begin(), vec1.end(), vec2.begin());
}

inline size_t hash_value(const IntegerVec& vec)
{
	size_t seed = 0;
	for (IntegerVecConstIter elementIter = vec.begin(); elementIter != vec.end(); elementIter++)
	{
		hash_combine(seed, *elementIter);
	}
    return seed;
}

struct ReadID
{
	union
	{
		struct
		{
			unsigned fragmentIndex : 31;
			unsigned readEnd : 1;
		};

		int id;
	};
};

inline bool operator==(const ReadID& readid1, const ReadID& readid2)
{
	return readid1.id == readid2.id;
}

inline size_t hash_value(const ReadID& readid)
{
	size_t seed = 0;
	hash_combine(seed, readid.id);
	return seed;
}

int OtherClusterEnd(int clusterEnd);

struct MatePair
{
	double x;
	double y;
	double u;
	double s;
};

typedef vector<MatePair> MatePairVec;

double normalpdf(double x, double mu, double sigma);

class RandomNumberGenerator
{
public:
	explicit RandomNumberGenerator(int seed = 2014)
	 : mGenerator(seed)
	{}

	template <typename TType>
	TType Next(TType l, TType u)
	{
		return boost::random::uniform_int_distribution<TType>(l, u)(mGenerator);
	}
	
private:
	boost::random::mt19937 mGenerator;
};

class RandomGenomicPositionGenerator
{
public:
	explicit RandomGenomicPositionGenerator(const vector<long>& chromosomeLengths, int seed = 2014)
	 : mChromosomeLengths(chromosomeLengths), mGenerator(seed)
	{
		mGenomeLength = accumulate(chromosomeLengths.begin(), chromosomeLengths.end(), 0L);
		mDistribution = boost::random::uniform_int_distribution<long>(1, mGenomeLength);
	}
	
	void Next(int& chrIdx, long& position)
	{
		position = mDistribution(mGenerator);
		
		chrIdx = 0;
		while (chrIdx < mChromosomeLengths.size())
		{
			if (position <= mChromosomeLengths[chrIdx])
			{
				break;
			}
			
			position -= mChromosomeLengths[chrIdx];
			
			chrIdx++;
		}
	}
	
private:
	const vector<long>& mChromosomeLengths;
	boost::random::mt19937 mGenerator;
	boost::random::uniform_int_distribution<long> mDistribution;
	long mGenomeLength;
};

void CheckFile(const ios& file, const string& filename);

template<typename TType>
TType SafeParseField(const string& field, const char* nameOfType, const char* codeFilename, int codeLine, const string& parseFilename = string(), int parseLine = 0)
{
    using boost::lexical_cast;
    using boost::bad_lexical_cast;
	
	try
	{
		return lexical_cast<TType>(field);
	}
	catch (bad_lexical_cast &)
	{
		if (parseFilename.empty())
		{
			cerr << "error interpreting '" << field << "' as " << nameOfType << endl;
		}
		else
		{
			cerr << "error interpreting '" << field << "' as " << nameOfType << " for " << parseFilename << ":" << parseLine << endl;
		}
		cerr << "parsing failed at " << codeFilename << ":" << codeLine << endl;
		exit(1);
	}
}

#define SAFEPARSE(type, field) SafeParseField<type>(field, #type, __FILE__, __LINE__)

#define SAFEPARSEFIELD(type, field, filename, line) SafeParseField<type>(field, #type, __FILE__, __LINE__, filename, line)

void Print(const IntegerVecMap& clusters);
void Print(const IntegerVec& cluster);
void Print(const IntegerTable& cluster);

bool ReadTSV(istream& file, StringVec& fields);

void ReadFAI(const string& faiFilename, vector<string>& referenceNames, vector<long>& referenceLengths);

#endif
