/*
 *  AlignmentRecord.cpp
 *
 */

#include "AlignmentRecord.h"
#include "DebugCheck.h"

#include <ostream>

AlignmentKey SpanningAlignmentRecord::GetAlignmentKey() const
{
	AlignmentKey alignKey;
	alignKey.libID = libID;
	alignKey.readID = readID;
	alignKey.readEnd = readEnd;
	alignKey.alignID = alignID;
	return alignKey;
}

std::ostream & operator<<(std::ostream &os, const SpanningAlignmentRecord& record)
{
	os << record.libID << "\t";
	os << record.readID << "\t";
	os << record.readEnd << "\t";
	os << record.alignID << "\t";
	os << record.chromosome << "\t";
	os << record.strand << "\t";
	os << record.position << "\t";
	os << record.alignedLength << "\t";
	os << record.mateLength << "\t";
	os << record.mateScore << std::endl;
	return os;
}

std::istream & operator>>(std::istream &is, SpanningAlignmentRecord& record)
{
	is >> record.libID;
	is >> record.readID;
	is >> record.readEnd;
	is >> record.alignID;
	is >> record.chromosome;
	is >> record.strand;
	is >> record.position;
	is >> record.alignedLength;
	is >> record.mateLength;
	is >> record.mateScore;
	return is;
}

std::ostream & operator<<(std::ostream &os, const SplitAlignmentRecord& record)
{
	os << record.libID << "\t";
	os << record.readID << "\t";
	os << record.readEnd << "\t";
	for (int readEnd = 0; readEnd < 2; readEnd++)
	{
		os << record.alignID[readEnd] << "\t";
		os << record.chromosome[readEnd] << "\t";
		os << record.strand[readEnd] << "\t";
		os << record.position[readEnd] << "\t";
	}
	os << record.homology << "\t";
	if (record.inserted.empty())
	{
		os << ".\t";
	}
	else
	{
		os << record.inserted << "\t";
	}
	os << record.score << std::endl;
	return os;
}

std::istream & operator>>(std::istream &is, SplitAlignmentRecord& record)
{
	is >> record.libID;
	is >> record.readID;
	is >> record.readEnd;
	for (int readEnd = 0; readEnd < 2; readEnd++)
	{
		is >> record.alignID[readEnd];
		is >> record.chromosome[readEnd];
		is >> record.strand[readEnd];
		is >> record.position[readEnd];
	}
	is >> record.homology;
	is >> record.inserted;
	if (record.inserted == ".")
	{
		record.inserted = "";
	}
	is >> record.score;
	return is;
}

AlignmentPairKey SplitAlignmentRecord::GetAlignmentPairKey() const
{
	AlignmentPairKey alignPairKey;
	alignPairKey.libID = libID;
	alignPairKey.readID = readID;
	for (int readEnd = 0; readEnd < 2; readEnd++)
	{
		alignPairKey.alignID[readEnd] = alignID[readEnd];
	}
	return alignPairKey;
}

std::ostream & operator<<(std::ostream &os, const ClusterMemberRecord& record)
{
	os << record.clusterID << "\t";
	os << record.clusterEnd << "\t";
	os << record.libID << "\t";
	os << record.readID << "\t";
	os << record.readEnd << "\t";
	os << record.alignID << std::endl;
	return os;
}

std::istream & operator>>(std::istream &is, ClusterMemberRecord& record)
{
	is >> record.clusterID;
	is >> record.clusterEnd;
	is >> record.libID;
	is >> record.readID;
	is >> record.readEnd;
	is >> record.alignID;
	return is;
}

ReadRecord ClusterMemberRecord::GetReadRecord() const
{
	ReadRecord record;
	record.libID = libID;
	record.readID = readID;
	return record;
}

AlignmentKey ClusterMemberRecord::GetAlignmentKey() const
{
	AlignmentKey alignKey;
	alignKey.libID = libID;
	alignKey.readID = readID;
	alignKey.readEnd = readEnd;
	alignKey.alignID = alignID;
	return alignKey;
}

std::ostream & operator<<(std::ostream &os, const BreakpointRecord& record)
{
	os << record.clusterID << "\t";
	os << record.breakpointID << "\t";
	for (int clusterEnd = 0; clusterEnd < 2; clusterEnd++)
	{
		os << record.chromosome[clusterEnd] << "\t";
		os << record.strand[clusterEnd] << "\t";
		os << record.position[clusterEnd] << "\t";
	}
	os << record.count << "\t";
	os << record.homology << "\t";
	os << record.inserted << "\t";
	os << record.mateScore << std::endl;
	return os;
}

std::istream & operator>>(std::istream &is, BreakpointRecord& record)
{
	is >> record.clusterID;
	is >> record.breakpointID;
	for (int clusterEnd = 0; clusterEnd < 2; clusterEnd++)
	{
		is >> record.chromosome[clusterEnd];
		is >> record.strand[clusterEnd];
		is >> record.position[clusterEnd];
	}
	is >> record.count;
	is >> record.homology;
	is >> record.inserted;
	if (record.inserted == ".")
	{
		record.inserted = "";
	}
	is >> record.mateScore;
	return is;
}

std::ostream & operator<<(std::ostream &os, const BreakAlignScoreRecord& record)
{
	os << record.clusterID << "\t";
	os << record.breakpointID << "\t";
	os << record.clusterEnd << "\t";
	os << record.libID << "\t";
	os << record.readID << "\t";
	os << record.readEnd << "\t";
	os << record.alignID << "\t";
	os << record.alignedLength << "\t";
	os << record.templateLength << "\t";
	os << record.score << std::endl;
	return os;
}

std::istream & operator>>(std::istream &is, BreakAlignScoreRecord& record)
{
	is >> record.clusterID;
	is >> record.breakpointID;
	is >> record.clusterEnd;
	is >> record.libID;
	is >> record.readID;
	is >> record.readEnd;
	is >> record.alignID;
	is >> record.alignedLength;
	is >> record.templateLength;
	is >> record.score;
	return is;
}

