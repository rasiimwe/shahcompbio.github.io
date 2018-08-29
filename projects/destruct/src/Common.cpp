/*
 *  Common.cpp
 *
 */

#include "Common.h"
#include "DebugCheck.h"

#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>


using namespace std;
using namespace boost;


int OtherStrand(int strand)
{
	return (1 - strand);
}

int OtherReadEnd(int readEnd)
{
	return (1 - readEnd);
}

int OtherClusterEnd(int clusterEnd)
{
	return (1 - clusterEnd);
}

void ReverseComplement(string& sequence)
{
	reverse(sequence.begin(), sequence.end());
	
	for (int seqIndex = 0; seqIndex < sequence.size(); seqIndex++)
	{
		char nucleotide = sequence[seqIndex];
		
		switch (nucleotide)
		{
			case 'A': nucleotide = 'T'; break;
			case 'C': nucleotide = 'G'; break;
			case 'T': nucleotide = 'A'; break;
			case 'G': nucleotide = 'C'; break;
			case 'a': nucleotide = 't'; break;
			case 'c': nucleotide = 'g'; break;
			case 't': nucleotide = 'a'; break;
			case 'g': nucleotide = 'c'; break;
		}
		
		sequence[seqIndex] = nucleotide;
	}
}

double normalpdf(double x, double mu, double sigma)
{
	double coeff = 1.0 / ((double)sigma * sqrt(2 * M_PI));
	
	double dist = (((double)x - (double)mu) / (double)sigma);
	double exponent = -0.5 * dist * dist;
	
	return coeff * exp(exponent);
}

int InterpretStrand(const string& strand)
{
	DebugCheck(strand == "+" || strand == "-");

	if (strand == "+")
	{
		return PlusStrand;
	}
	else
	{
		return MinusStrand;
	}
}

void CheckFile(const ios& file, const string& filename)
{
	if (!file.good())
	{
		cerr << "Error: Unable to open " << filename << endl;
		exit(1);
	}	
}

void Print(const IntegerVecMap& clusters)
{
	for (IntegerVecMapConstIter clusterIter = clusters.begin(); clusterIter != clusters.end(); clusterIter++)
	{
		for (IntegerVecConstIter elementIter = clusterIter->second.begin(); elementIter != clusterIter->second.end(); elementIter++)
		{
			cout << clusterIter->first << "\t" << *elementIter << endl;
		}
	}
}

void Print(const IntegerVec& cluster)
{
	for (IntegerVecConstIter elementIter = cluster.begin(); elementIter != cluster.end(); elementIter++)
	{
		cout << *elementIter << endl;
	}
}

void Print(const IntegerTable& clusters)
{
	for (IntegerTableConstIter clusterIter = clusters.begin(); clusterIter != clusters.end(); clusterIter++)
	{
		for (IntegerVecConstIter elementIter = clusterIter->begin(); elementIter != clusterIter->end(); elementIter++)
		{
			cout << *elementIter << "\t";
		}
		cout << endl;
	}
}

bool ReadTSV(istream& file, StringVec& fields)
{
	string line;
	if (!getline(file, line))
	{
		return false;
	}
	
	split(fields, line, is_any_of("\t"));
	
	return true;
}

void ReadFAI(const string& faiFilename, vector<string>& referenceNames, vector<long>& referenceLengths)
{
	ifstream faiFile(faiFilename.c_str());
	if (!faiFile)
	{
		cerr << "Error: unable to open fai file" << endl;
		exit(1);
	}
	
	string line;
	int lineNumber = 0;
	while (getline(faiFile, line))
	{
		lineNumber++;
		
		if (line.length() == 0)
		{
			cerr << "Error: Empty fai line " << lineNumber << " of " << faiFilename << endl;
			exit(1);
		}
		
		vector<string> faiFields;
		split(faiFields, line, is_any_of("\t"));
		
		if (faiFields.size() < 2)
		{
			cerr << "Error: Format error for fai line " << lineNumber << " of " << faiFilename << endl;
			exit(1);
		}
		
		string referenceName = faiFields[0];
		long referenceLength = SAFEPARSE(long, faiFields[1]);
		
		referenceNames.push_back(referenceName);
		referenceLengths.push_back(referenceLength);
	}
	
	faiFile.close();
}


