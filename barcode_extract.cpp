#define NOMINMAX

#include "string"
#include "vector"
#include "fstream"
#include "cassert"
#include "map"
#include "set"
#include <algorithm>
#include <iostream>
#include <seqan/seq_io.h>
#include <seqan/find.h>
#include <seqan/index.h>
#include <seqan/arg_parse.h>

using namespace std;

vector<string> mustHave = { "internal_bc_us", "internal_bc_ds" };  //"sample_name", 

//parameters from user
int us_start_search, us_search_range, ds_start_search, ds_search_range;
unsigned int desired_edit_distance;
bool paired;
size_t max_counts_per_file = 100000, batch_read = 10000, stop_after=300000;

std::string sample_sheet; // = "/mnt/c/Users/constantine.chrysost/Documents/fastq_test_files/SampleSheetPipeline.txt"; //"/home/costa/Documents/fastq_test_files/SampleSheetPipeline.txt";
std::string input_r1_read; // = "/mnt/c/Users/constantine.chrysost/Documents/fastq_test_files/Undetermined_S0_L001_R1_001.fastq.gz";  // "/home/costa/Documents/fastq_test_files/Undetermined_S0_L001_R1_001.small.fastq";
std::string input_r2_read; // = "/mnt/c/Users/constantine.chrysost/Documents/fastq_test_files/Undetermined_S0_L001_R2_001.fastq.gz"; // "/home/costa/Documents/fastq_test_files/Undetermined_S0_L001_R2_001.small.fastq";
std::string prefix;

typedef seqan::Dna5String DnaSeq;
typedef seqan::StringSet<DnaSeq> DnaList;
typedef seqan::SeqFileIn MyFileIn;
typedef seqan::SeqFileOut MyFileOut;
typedef seqan::IndexEsa<> MyIndex;

struct Samples {
	/*
		Read tab delimited sample sheet and store each column
	*/
	vector<string> data;
	map<string, int> map_col_to_pos;
	set<string>::iterator dsId, usId;
	DnaSeq myUSSeq, myDSSeq;

	Samples() {}

	Samples(vector<string> _header, vector<string> _data, set<string> & usBc, set<string> & dsBc, set<string> & combinedBc) {
		assert(_header.size() == _data.size());
		for (unsigned int i = 0; i < _header.size(); i++) {			
			map_col_to_pos[_header[i]] = i;
		}
		data = _data;
		myUSSeq = data[map_col_to_pos.at("internal_bc_us")];
		myDSSeq = data[map_col_to_pos.at("internal_bc_ds")];

		usBc.insert(data[map_col_to_pos.at("internal_bc_us")]);
		dsBc.insert(data[map_col_to_pos.at("internal_bc_ds")]);
		combinedBc.insert(data[map_col_to_pos.at("internal_bc_us")] + "+" + data[map_col_to_pos.at("internal_bc_ds")]);

		dsId = dsBc.find(data[map_col_to_pos.at("internal_bc_ds")]);
		usId = usBc.find(data[map_col_to_pos.at("internal_bc_us")]);		
	}
};

struct SampleList {
	vector<Samples> samples;
	vector<string> header;
	set<string> usBc = { "" }, dsBc = { "" }, combinedBc;
	// std::vector<std::vector<std::pair<seqan::CharString, bool>> bc2uniquePair; // make sure a unique set of US and DS sequences exists in sample sheet
	vector < vector<std::pair<seqan::CharString, bool>>> bc2uniquePair;

	SampleList(vector<string> _header) {
		header = _header;
	}

	void InsertSample(const vector<string> & _data) {
		samples.push_back(Samples(header, _data, usBc, dsBc, combinedBc));
	}

	void MapBC2Sample(DnaList bc1, DnaList bc2) {
		seqan::appendValue(bc1, "");
		seqan::appendValue(bc2, "");
		for (int i = 0; i < length(bc1); i++) {
			vector<std::pair<seqan::CharString, bool>> tmp;
			for (int j = 0; j < length(bc2); j++) {
				seqan::CharString storeBC = seqan::CharString(bc1[i]);
				storeBC += "+";
				storeBC += seqan::CharString(bc2[j]);
				tmp.push_back(std::pair<seqan::CharString, bool>(storeBC, false));
				for (auto sample : samples) {
					if ((sample.myUSSeq == bc1[i]) && (sample.myDSSeq == bc2[j])) {
						// this pair exists at least ONCE in sample sheet
						tmp[j].second = true;
					}
				}
			}
			bc2uniquePair.push_back(tmp);
		}
	}
};

struct MultipleHaystack {
	/*
		Create a struct that stores both a finder object for a vetor of barcodes, and also, a mapped position mapping each finder position to a barcode of interest
	*/
	seqan::Finder<seqan::CharString> finder;
	vector<int> pos2index;
	seqan::CharString haystackString, cushion;

	MultipleHaystack() {};
	MultipleHaystack(const DnaList & _barcodes) {
		_InitStruct(_barcodes);
	}
	void _InitStruct(const DnaList & barcodes) {
		for (unsigned int i = 0; i < desired_edit_distance; i++) {
			cushion += ".";
		}

		int currPos = 0, counter = 0;
		for (auto bc : barcodes) {
			haystackString += cushion;
			haystackString += bc;
			haystackString += cushion;
			haystackString += "$";
			currPos = pos2index.size();
			for (int k = currPos; k < length(haystackString); k++) pos2index.push_back(counter);
			counter += 1;
		}

		finder = seqan::Finder<seqan::CharString>(haystackString);
	}
};

struct BatchRecordStore {
	seqan::StringSet<seqan::CharString> header;
	seqan::StringSet<seqan::Dna5QString> seq;

	void append(const seqan::CharString & h, const seqan::Dna5QString & s) {
		appendValue(header, h);
		appendValue(seq, s);		
	}

	void clear() {
		seqan::clear(header);
		seqan::clear(seq);		
	}
};

struct FileWriterHandler {
	std::vector<std::vector<size_t>> currentCounts, currentIter;
	MyFileOut** filePaths;
	MyFileOut noBC, multipleBC;
	size_t maxCountsPerFile;

	std::vector<std::vector<BatchRecordStore>> seqsToWrite;
	std::vector<std::vector<std::pair<string, string>>> filePathNames;

	BatchRecordStore noBCSeqs, multipleBCSeqs;

	FileWriterHandler(
		const SampleList & sampleSheet,
		std::pair<seqan::CharString, seqan::CharString> filePrefixSuffix,
		size_t _maxCountsPerFile
	) {
		maxCountsPerFile = _maxCountsPerFile;
		string filePrefix = toCString(filePrefixSuffix.first), fileSuffix = toCString(filePrefixSuffix.second);
		filePaths = new MyFileOut*[sampleSheet.bc2uniquePair.size()];
		seqan::open(noBC, (filePrefix + ".barcode_notfound" + fileSuffix).c_str());
		seqan::open(multipleBC, (filePrefix + ".barcode_multiplefound" + fileSuffix).c_str());
		for (int xT = 0; xT < sampleSheet.bc2uniquePair.size(); xT++) {
			std::vector<size_t> tmp, tmpIter;
			std::vector<BatchRecordStore> tmpBatch;
			std::vector<std::pair<string, string>> tmpNames;
			filePaths[xT] = new MyFileOut[sampleSheet.bc2uniquePair[xT].size()];
			for (int yT = 0; yT < sampleSheet.bc2uniquePair[xT].size(); yT++) {
				tmp.push_back(0);
				tmpBatch.push_back(BatchRecordStore());
				if (sampleSheet.bc2uniquePair[xT][yT].second) {					
					seqan::open(filePaths[xT][yT], (filePrefix + "_bc_" + toCString(sampleSheet.bc2uniquePair[xT][yT].first) + fileSuffix).c_str());
					tmpIter.push_back(1);
					tmpNames.push_back(std::pair<string, string>((filePrefix + "_bc_" + toCString(sampleSheet.bc2uniquePair[xT][yT].first)).c_str(), fileSuffix.c_str()));
				}
				else {
					tmpIter.push_back(0);
					tmpNames.push_back(std::pair<string, string>("", ""));
				}
			}
			currentIter.push_back(tmpIter);
			currentCounts.push_back(tmp);
			seqsToWrite.push_back(tmpBatch);
			filePathNames.push_back(tmpNames);
		}
	}

	~FileWriterHandler() {		
		for (int xT = 0; xT < currentCounts.size(); xT++) {
			for (int yT = 0; yT < currentCounts[xT].size(); yT++) {				
				seqan::close(filePaths[xT][yT]);
			}			
			delete[] filePaths[xT];
		}		
		delete[] filePaths;
	}

	void addData(size_t i, size_t j, const seqan::CharString & h, const seqan::Dna5QString & s) {
		currentCounts[i][j] += 1;
		seqsToWrite[i][j].append(h, s);
		if (maxCountsPerFile != 0 && currentCounts[i][j] == maxCountsPerFile) {
			// cout << "im writing!" << endl;
			// write all remaining data to file, and then close the file andopen a new one
			seqan::writeRecords(filePaths[i][j], seqsToWrite[i][j].header, seqsToWrite[i][j].seq);
			seqsToWrite[i][j].clear();
			openNewFile(i, j);
		}
	}

	void writeSingleSeq(size_t i, size_t j, const seqan::CharString & h, const seqan::Dna5QString & s) {
		currentCounts[i][j] += 1;
		// seqsToWrite[i][j].append(h, s);
		seqan::writeRecord(filePaths[i][j], h, s);
		if (maxCountsPerFile != 0 && currentCounts[i][j] == maxCountsPerFile) {
			openNewFile(i, j);
		}
	}

	void openNewFile(size_t i, size_t j) {
		currentCounts[i][j] = 0;
		seqan::close(filePaths[i][j]);
		string newFile = filePathNames[i][j].first + "_batch_" + ::to_string(currentIter[i][j]) + filePathNames[i][j].second;
		// cout << "OPENING NEW FILE: " << newFile << endl;
		currentIter[i][j] += 1;
		seqan::open(filePaths[i][j], newFile.c_str());
	}

	void writeData() {
		for (int i = 0; i < currentCounts.size(); i++) {
			for (int j = 0; j < currentCounts[i].size(); j++) {
				if (currentCounts[i][j] > 0) {
					// write data from all arrays in each barcode to file
					seqan::writeRecords(filePaths[i][j], seqsToWrite[i][j].header, seqsToWrite[i][j].seq);
					seqsToWrite[i][j].clear();
				}
			}
		}

		// write not found data to file 
		seqan::writeRecords(noBC, noBCSeqs.header, noBCSeqs.seq);
		seqan::writeRecords(multipleBC, multipleBCSeqs.header, multipleBCSeqs.seq);

		noBCSeqs.clear();
		multipleBCSeqs.clear();
	}
};


std::ostream& operator<<(std::ostream& os, const Samples& m)
{
	// return os << "Data: " << m.data[m.map_col_to_pos.at("sample_name")] << "," << m.data[m.map_col_to_pos.at("internal_bc_us")] << "," << m.data[m.map_col_to_pos.at("internal_bc_ds")];
	return os << "Data: " << m.data[m.map_col_to_pos.at("internal_bc_us")] << "," << m.data[m.map_col_to_pos.at("internal_bc_ds")];
}

vector<string> SplitString(string s, const string & delimiter) {
	size_t pos = 0;
	std::string token;
	vector<string> stringVec;
	while ((pos = s.find(delimiter)) != std::string::npos) {
		token = s.substr(0, pos);
		stringVec.push_back(token);
		s.erase(0, pos + delimiter.length());
	}
	stringVec.push_back(s);
	return stringVec;
}

std::pair<seqan::CharString, seqan::CharString> ExtractReadPrefix(const seqan::CharString  & s) {
	/*
	Extract header string from illumina (ignoring the R1/R2)
	*/
	size_t pos = 0;

	while (pos < length(s)) {
		if (s[pos] == ' ') break;
		pos++;
	}

	return std::pair<seqan::CharString, seqan::CharString>(seqan::prefix(s, pos), seqan::infix(s, pos, length(s)));
}

std::pair<seqan::CharString, seqan::CharString> ExtractFileNamePrefix(const seqan::CharString  & s) {
	/*
	Extract header string from illumina (ignoring the R1/R2)
	*/

	// first try to guess all expected file types....
	vector<seqan::CharString> types = { "fastq", "fq", "fa","fasta" };
	vector<seqan::CharString> ext = { "", ".gz", ".bz2" };
	seqan::CharString suffixString, finalExt, finalPrefix;

	bool found = false;

	for (int i = 0; i < types.size(); i++) {
		for (int j = 0; j < ext.size(); j++) {
			suffixString = ".";
			suffixString += types[i];
			suffixString += ext[j];
			assert(length(seqan::suffix(s, length(s) - length(suffixString))) == length(suffixString));
			if (seqan::suffix(s, length(s) - length(suffixString)) == suffixString) {
				finalExt = suffixString;
				finalPrefix = seqan::prefix(s, length(s) - length(suffixString));
				found = true;
			}
		}
	}

	if (found == false) {
		// coulndt figure it out...default to a fastq file
		std::cerr << "Warning: could not determine the file extension for the input file " << s << endl;
		size_t pos = length(s) - 1;
		while (pos > 0) {
			if (s[pos] == '.') break;
			pos--;
		}

		finalPrefix = seqan::prefix(s, pos);
		finalExt = ".fq";
	}

	return std::pair<seqan::CharString, seqan::CharString>(finalPrefix, finalExt);
}

SampleList ParseSampleSheet(std::string file) {
	ifstream sampleFile;
	sampleFile.open(file);
	assert(sampleFile.is_open());
	string line;

	SampleList *barcodeInfo;

	vector<string> header;
	int line_counter = -1;
	while (!sampleFile.eof()) {
		getline(sampleFile, line);
		// cout << line << endl;
		line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
		line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
		if (line.size() == 0) continue;

		vector<string> data = SplitString(line, "\t");
		line_counter += 1;

		if (line_counter == 0) {
			std::copy(data.begin(), data.end(), back_inserter(header));
			for (unsigned int i = 0; i < mustHave.size(); i++) {
				bool found = false;
				for (unsigned int j = 0; j < header.size(); j++) {
					if (header[j] == mustHave[i]) {
						found = true;
						break;
					}
				}
				if (found == false) {
					cerr << "Error could not find the following column name in sample sheet: " << mustHave[i] << endl;
					exit(1);
				}
			}
			barcodeInfo = new SampleList(header);
			continue;
		}
		barcodeInfo->InsertSample(data);
	}
	return *barcodeInfo;
}

int SearchBarcodesV2(
	seqan::CharString query,
	seqan::StringSet<seqan::Pattern<seqan::CharString, seqan::Myers<>>> barcodePatterns,
	unsigned int editDistance
)
/*
More accurate than searchbarcodesV1 but can be much slower given the length of the barcodepatterns
*/
{
	int minScore = -1 * editDistance;
	int hitPos = -1, numHits = 0;
	seqan::Finder<seqan::CharString> finder(query);
	for (unsigned int i = 0; i < length(barcodePatterns); i++) {
		if (find(finder, barcodePatterns[i], minScore)) {
			numHits += 1;
			hitPos = i;
		}
		if (numHits > 1) {
			hitPos = -2;
			break;
		}

		goBegin(finder);
		clear(finder);
	}

	return hitPos;
}

int SearchBarcodes(seqan::CharString query, MultipleHaystack haystack, unsigned int editDistance) {
	/*
	Search a query string as a pattern against all possible barcodes which are stitched together in a finder/haystack object

	.. note::Accuracy

	this is less accurate than searching each barcode individuall agains tthe query because it will not catch differences we find where the query has a deletion and therefore is an edit distance of 1
	Instead it will think this read has an edit distance of 2 because now the query sequence has an extra base, X, at the end due to the deletion (....-...X). This base should NOT be part of the pattern. So it will miss a few of these solutions
	*/
	int hitPos = -1, numHits = 0, tmpHitPos;
	seqan::Pattern<seqan::CharString, seqan::Myers<>> pattern(query);

	int minScore = -1 * editDistance;

	while (find(haystack.finder, pattern, minScore)) {
		tmpHitPos = haystack.pos2index[position(haystack.finder)];
		if (tmpHitPos != hitPos) {
			numHits += 1;
			hitPos = tmpHitPos;
		}

		if (numHits > 1) {
			hitPos = -2;
			break;
		}
	}
	return hitPos;
}

int SearchExactMatchBarcodes(DnaSeq query, const DnaList & barcodes, int editDistance) {
	int hitPos = -1, numHits = 0;

	seqan::Finder<DnaSeq> finder(query);
	seqan::Pattern<DnaList, seqan::WuManber> pattern(barcodes);

	while (find(finder, pattern)) {
		if (seqan::position(pattern) != hitPos) {
			numHits += 1;
			hitPos = seqan::position(pattern);
		}
		if (numHits > 1) {
			hitPos = -2;
			break;
		}
	}

	return hitPos;
}

int SearchExactMatchBarcodesIndex(DnaSeq query, seqan::Finder<seqan::Index<DnaList, MyIndex>> & barcodeIdx) {
	int hitPos = -1, numHits = 0, tmp;
	seqan::Pair<int, int> found;

	// cout <<"a"<<endl;	
	while (find(barcodeIdx, query)) {
		found = seqan::beginPosition(barcodeIdx);
		if (found.i1 != hitPos) {
			numHits += 1;
			hitPos = found.i1;
		}
		if (numHits > 1) {
			hitPos = -2;
			break;
		}
	}
	seqan::clear(barcodeIdx);
	return hitPos;
}

void FMIndexSearch() {
	// seqan::Dna5String upstreamseqs = "";
	// for (int i = 0; i<length(headerR1List); i++){
	// 	upstreamseqs += infix(seqR1List[i], 0, 10);
	// 	upstreamseqs += "N";
	// }

	// seqan::Dna5String genome = upstreamseqs;
	// seqan::Index<seqan::Dna5String, seqan::BidirectionalIndex<seqan::FMIndex<> > > index(genome);
	// auto delegate = [](auto & iter, seqan::Dna5String const & needle, uint8_t errors)
	// {
	// 	for (auto occ : getOccurrences(iter))
	// 		std::cout << occ << std::endl;
	// };

	// seqan::Dna5String pattern("GGGGTTAT");
	// std::cout << "Hits with up to 2 errors (HammingDistance):" << std::endl;
	// find<0, 2>(delegate, index, pattern, seqan::HammingDistance());
	// find<0, 1>(delegate, index, pattern, seqan::HammingDistance());
	// cout <<"L: " << length(headerR1List)<<endl;
}

void SearchFastqFile(MyFileIn & R1, MyFileIn & R2, DnaList usBarcodes, const DnaList & dsBarcodes, const SampleList & sampleSheet, const std::pair<seqan::CharString, seqan::CharString> & outR1Prefix, const std::pair<seqan::CharString, seqan::CharString> & outR2Prefix) {
	/*
	Read through R1 and R2 files while extracting the position containing internal barcodes
	*/

	std::vector<std::vector<size_t>> occurrenceCounts;
	
	string fileNamePrefixR1 = toCString(outR1Prefix.first), fileNamePrefixR2 = toCString(outR2Prefix.first), fileNameSuffixR1 = toCString(outR1Prefix.second), fileNameSuffixR2 = toCString(outR2Prefix.second);

	FileWriterHandler outDataR1(sampleSheet, outR1Prefix, max_counts_per_file), outDataR2(sampleSheet, outR2Prefix, max_counts_per_file);

	for (int xT = 0; xT < sampleSheet.bc2uniquePair.size(); xT++) {
		std::vector<size_t> tmp;
		for (int yT = 0; yT < sampleSheet.bc2uniquePair[xT].size(); yT++) {
			tmp.push_back(0);
		}
		occurrenceCounts.push_back(tmp);
	}

	seqan::StringSet<seqan::Pattern<seqan::CharString, seqan::Myers<>>> usPattern, dsPattern;

	DnaSeq usSearch, dsSearch;
	size_t minUSSearch = length(usBarcodes[0]), minDSSearch = length(dsBarcodes[0]);

	vector<int> haystackUSMap, haystackDSMap;

	for (auto bc : usBarcodes) {
		// identify the minimum barcode length
		if (length(bc) != minUSSearch) {
			std::cerr << "Warning: 5prime barcode lengths are not the same. We will only search the fastq data for substrings based on the minimum barcode length in the list" << endl;
			minUSSearch = std::min(minUSSearch, length(bc));
		}
		seqan::appendValue(usPattern, seqan::Pattern<seqan::CharString, seqan::Myers<>>(seqan::CharString(bc)));
	}

	for (auto bc : dsBarcodes) {
		// identify the maximum barcode length
		if (length(bc) != minDSSearch) {
			std::cerr << "Warning: 3prime barcode lengths are not the same. We will only search the fastq data for substrings based on the minimum barcode length in the list" << endl;
			minDSSearch = std::min(minDSSearch, length(bc));
		}
		seqan::appendValue(dsPattern, seqan::Pattern<seqan::CharString, seqan::Myers<>>(seqan::CharString(bc)));
	}

	//create a string where all barcodes are stitched together. we can use this to create an online search structure, or a string index
	MultipleHaystack haystackUS(usBarcodes), haystackDS(dsBarcodes);

	//create an index using the barcode sequences. Each barcode will exist in an index position, but can be searched together. 
	seqan::Index<DnaList, MyIndex > usBCIdx(usBarcodes), dsBCIdx(dsBarcodes);
	seqan::Finder<seqan::Index<DnaList, MyIndex > > usBCIdxFinder(usBCIdx), dsBCIdxFinder(dsBCIdx);

	int hitUS, hitDS, breakCounter = 0, totalSeqs = 0;
	seqan::StringSet<seqan::CharString> headerR1List, headerR2List;
	std::pair<seqan::CharString, seqan::CharString> tmpR1Header, tmpR2Header;
	// seqan::CharString headerR1, headerR2;
	// seqan::Dna5QString seqR1, seqR2; // qualR1, qualR2;

	seqan::StringSet<seqan::Dna5QString> seqR1List, seqR2List; //, qualR1List, seqR2List, qualR2List;

	us_search_range = std::max(int(us_search_range), int(minUSSearch));
	ds_search_range = std::max(int(ds_search_range), int(minDSSearch));

	size_t notFound = 0, found = 0, multipleHits = 0;

	while (!seqan::atEnd(R1)) {
		clear(headerR1List);
		clear(headerR2List);
		clear(seqR1List);
		clear(seqR2List);

		seqan::readRecords(headerR1List, seqR1List, R1, batch_read);
		seqan::readRecords(headerR2List, seqR2List, R2, batch_read);

		for (unsigned int i = 0; i < length(headerR1List); i++) {
			totalSeqs += 1;

			// make sure R1/R2 reads match					
			tmpR1Header = ExtractReadPrefix(headerR1List[i]);
			tmpR2Header = ExtractReadPrefix(headerR2List[i]);
			assert(tmpR1Header.first == tmpR2Header.first);
			
			if (desired_edit_distance == 0) {
				// search for exact matches using index
				// search R1 file for upstream sequence
				for (int j = 0; j + minUSSearch <= us_search_range; j++) {
					usSearch = seqan::infix(seqR1List[i], us_start_search + j, us_start_search + j + minUSSearch);
					hitUS = SearchExactMatchBarcodesIndex(usSearch, usBCIdxFinder);
					if (hitUS != -1) {
						// either we found a hit, or too many hits were found...yes we exit if we find a hit even though there might be another hit downstream...we can make this a parameter in the future
						break;
					}
				}

				// search R2 file for downstream sequence
				for (int j = 0; j + minDSSearch <= ds_search_range; j++) {
					dsSearch = seqan::infix(seqR2List[i], ds_start_search + j, ds_start_search + j + minDSSearch);
					hitDS = SearchExactMatchBarcodesIndex(dsSearch, dsBCIdxFinder);
					if (hitDS != -1) {
						break;
					}
				}
			}
			else {
				// do an inexact search where we search all barcodes against full length query 
				usSearch = seqan::infix(seqR1List[i], us_start_search, us_start_search + us_search_range);
				hitUS = SearchBarcodesV2(usSearch, usPattern, desired_edit_distance);
				
				dsSearch = seqan::infix(seqR2List[i], ds_start_search, ds_start_search + ds_search_range);
				hitDS = SearchBarcodesV2(dsSearch, dsPattern, desired_edit_distance);
			}
			
			if (hitUS == -1) hitUS = length(usBarcodes);
			if (hitDS == -1) hitDS = length(dsBarcodes);

			if (hitUS == -2 || hitDS == -2) {
				headerR1List[i] = tmpR1Header.first;
				headerR2List[i] = tmpR1Header.first;
				headerR1List[i] += tmpR1Header.second;
				headerR2List[i] += tmpR2Header.second;
				multipleHits += 1;
				outDataR1.multipleBCSeqs.append(headerR1List[i], seqR1List[i]);
				outDataR2.multipleBCSeqs.append(headerR2List[i], seqR2List[i]);
			}
			else {
				occurrenceCounts[hitUS][hitDS] += 1;
				headerR1List[i] = tmpR1Header.first;
				headerR1List[i] += ":";				
				headerR1List[i] += sampleSheet.bc2uniquePair[hitUS][hitDS].first;
				headerR2List[i] = headerR1List[i];
				headerR1List[i] += tmpR1Header.second;
				headerR2List[i] += tmpR2Header.second;								
				if (sampleSheet.bc2uniquePair[hitUS][hitDS].second) {
					/* this barcode pair exists in the sample sheet*/
					outDataR1.addData(hitUS, hitDS, headerR1List[i], seqR1List[i]);
					outDataR2.addData(hitUS, hitDS, headerR2List[i], seqR2List[i]);
				}
				else {
					/*barcode pair was not found, so print to notfound seqs*/
					outDataR1.noBCSeqs.append(headerR1List[i], seqR1List[i]);
					outDataR2.noBCSeqs.append(headerR2List[i], seqR2List[i]);
				}
			}
			
		}

		outDataR1.writeData();
		outDataR2.writeData();
		
		if (stop_after > 0 && totalSeqs > stop_after) {
			cout << "TOTAL " << totalSeqs << endl;
			break;
		}
	}

	outDataR1.writeData();
	outDataR2.writeData();
	
	std::cout << "Results: " << endl;

	for (int xT = 0; xT < occurrenceCounts.size(); xT++) {
		for (int yT = 0; yT < occurrenceCounts[xT].size(); yT++) {
			if (sampleSheet.bc2uniquePair[xT][yT].second) {
				cout << "Barcode " << sampleSheet.bc2uniquePair[xT][yT].first << ": " << occurrenceCounts[xT][yT] << endl;
				found += occurrenceCounts[xT][yT];
			}
			else {
				notFound += occurrenceCounts[xT][yT];
			}
		}
	}
	cout << endl << endl;
	cout << "Found: " << found << endl;
	cout << "Not Found: " << notFound << endl;
	cout << "Multiple found: " << multipleHits << endl;

	std::cout << endl << endl << "Results from not found: " << endl;
	for (int xT = 0; xT < occurrenceCounts.size(); xT++) {
		for (int yT = 0; yT < occurrenceCounts[xT].size(); yT++) {
			if (!sampleSheet.bc2uniquePair[xT][yT].second) {
				cout << "Barcode " << sampleSheet.bc2uniquePair[xT][yT].first << ": " << occurrenceCounts[xT][yT] << endl;
			}
		}
	}
}

void SearchFastqFile(MyFileIn & R, DnaList usBarcodes, const DnaList & dsBarcodes, const SampleList & sampleSheet, const std::pair<seqan::CharString, seqan::CharString> & outRPrefix) {
	std::vector<std::vector<size_t>> occurrenceCounts;

	string fileNamePrefixR = toCString(outRPrefix.first), fileNameSuffixR = toCString(outRPrefix.second);

	FileWriterHandler outDataR(sampleSheet, outRPrefix, max_counts_per_file);

	for (int xT = 0; xT < sampleSheet.bc2uniquePair.size(); xT++) {
		std::vector<size_t> tmp;
		for (int yT = 0; yT < sampleSheet.bc2uniquePair[xT].size(); yT++) {
			tmp.push_back(0);
		}
		occurrenceCounts.push_back(tmp);
	}

	seqan::StringSet<seqan::Pattern<seqan::CharString, seqan::Myers<>>> usPattern, dsPattern;

	DnaSeq usSearch, dsSearch;
	size_t minUSSearch = length(usBarcodes[0]), minDSSearch = length(dsBarcodes[0]);

	vector<int> haystackUSMap, haystackDSMap;

	for (auto bc : usBarcodes) {
		// identify the minimum barcode length
		if (length(bc) != minUSSearch) {
			std::cerr << "Warning: 5prime barcode lengths are not the same. We will only search the fastq data for substrings based on the minimum barcode length in the list" << endl;
			minUSSearch = std::min(minUSSearch, length(bc));
		}
		seqan::appendValue(usPattern, seqan::Pattern<seqan::CharString, seqan::Myers<>>(seqan::CharString(bc)));
	}

	for (auto bc : dsBarcodes) {
		// identify the maximum barcode length
		if (length(bc) != minDSSearch) {
			std::cerr << "Warning: 3prime barcode lengths are not the same. We will only search the fastq data for substrings based on the minimum barcode length in the list" << endl;
			minDSSearch = std::min(minDSSearch, length(bc));
		}
		seqan::appendValue(dsPattern, seqan::Pattern<seqan::CharString, seqan::Myers<>>(seqan::CharString(bc)));
	}

	//create a string where all barcodes are stitched together. we can use this to create an online search structure, or a string index
	MultipleHaystack haystackUS(usBarcodes), haystackDS(dsBarcodes);

	//create an index using the barcode sequences. Each barcode will exist in an index position, but can be searched together. 
	seqan::Index<DnaList, MyIndex > usBCIdx(usBarcodes), dsBCIdx(dsBarcodes);
	seqan::Finder<seqan::Index<DnaList, MyIndex > > usBCIdxFinder(usBCIdx), dsBCIdxFinder(dsBCIdx);

	int hitUS, hitDS, breakCounter = 0, totalSeqs = 0;
	seqan::StringSet<seqan::CharString> headerR1List;
	std::pair<seqan::CharString, seqan::CharString> tmpRHeader;
	
	seqan::StringSet<seqan::Dna5QString> seqR1List; //, qualR1List, seqR2List, qualR2List;

	us_search_range = std::max(int(us_search_range), int(minUSSearch));
	ds_search_range = std::max(int(ds_search_range), int(minDSSearch));

	size_t notFound = 0, found = 0, multipleHits = 0;

	while (!seqan::atEnd(R)) {
		clear(headerR1List);		
		clear(seqR1List);		

		seqan::readRecords(headerR1List, seqR1List, R, batch_read);		

		for (unsigned int i = 0; i < length(headerR1List); i++) {
			totalSeqs += 1;

			// make sure R1/R2 reads match					
			tmpRHeader = ExtractReadPrefix(headerR1List[i]);


			if (desired_edit_distance == 0) {
				// search for exact matches using index
				// search R1 file for upstream sequence
				for (int j = 0; j + minUSSearch <= us_search_range; j++) {
					usSearch = seqan::infix(seqR1List[i], us_start_search + j, us_start_search + j + minUSSearch);
					hitUS = SearchExactMatchBarcodesIndex(usSearch, usBCIdxFinder);
					if (hitUS != -1) {
						// either we found a hit, or too many hits were found...yes we exit if we find a hit even though there might be another hit downstream...we can make this a parameter in the future
						break;
					}
				}

				// search R1 file for downstream sequence
				for (int j = 0; j + minDSSearch <= ds_search_range; j++) {
					dsSearch = seqan::infix(seqR1List[i], length(seqR1List[i]) - ds_start_search + j, length(seqR1List[i]) - ds_start_search + j + minDSSearch);
					hitDS = SearchExactMatchBarcodesIndex(dsSearch, dsBCIdxFinder);
					if (hitDS != -1) {
						break;
					}
				}
			}
			else {
				// do an inexact search where we search all barcodes against full length query 
				usSearch = seqan::infix(seqR1List[i], us_start_search, us_start_search + us_search_range);
				hitUS = SearchBarcodesV2(usSearch, usPattern, desired_edit_distance);

				dsSearch = seqan::infix(seqR1List[i], length(seqR1List[i]) - ds_start_search, length(seqR1List[i]) - ds_start_search + ds_search_range);
				hitDS = SearchBarcodesV2(dsSearch, dsPattern, desired_edit_distance);
			}


			if (hitUS == -1) hitUS = length(usBarcodes);
			if (hitDS == -1) hitDS = length(dsBarcodes);

			if (hitUS == -2 || hitDS == -2) {
				headerR1List[i] = tmpRHeader.first;				
				headerR1List[i] += tmpRHeader.second;				
				multipleHits += 1;
				outDataR.multipleBCSeqs.append(headerR1List[i], seqR1List[i]);				
			}
			else {
				occurrenceCounts[hitUS][hitDS] += 1;
				headerR1List[i] = tmpRHeader.first;
				headerR1List[i] += ":";
				headerR1List[i] += sampleSheet.bc2uniquePair[hitUS][hitDS].first;				
				headerR1List[i] += tmpRHeader.second;				
				if (sampleSheet.bc2uniquePair[hitUS][hitDS].second) {
					/* this barcode pair exists in the sample sheet*/
					outDataR.addData(hitUS, hitDS, headerR1List[i], seqR1List[i]);					
				}
				else {
					/*barcode pair was not found, so print to notfound seqs*/
					outDataR.noBCSeqs.append(headerR1List[i], seqR1List[i]);					
				}
			}

		}

		outDataR.writeData();	

		if (stop_after > 0 && totalSeqs >= stop_after) {
			cout << "TOTAL " << totalSeqs << endl;
			break;
		}
	}

	outDataR.writeData();	

	std::cout << "Results: " << endl;

	for (int xT = 0; xT < occurrenceCounts.size(); xT++) {
		for (int yT = 0; yT < occurrenceCounts[xT].size(); yT++) {
			if (sampleSheet.bc2uniquePair[xT][yT].second) {
				cout << "Barcode " << sampleSheet.bc2uniquePair[xT][yT].first << ": " << occurrenceCounts[xT][yT] << endl;
				found += occurrenceCounts[xT][yT];
			}
			else {
				notFound += occurrenceCounts[xT][yT];
			}
		}
	}
	cout << endl << endl;
	cout << "Found: " << found << endl;
	cout << "Not Found: " << notFound << endl;
	cout << "Multiple found: " << multipleHits << endl;

	std::cout << endl << endl << "Results from not found: " << endl;
	for (int xT = 0; xT < occurrenceCounts.size(); xT++) {
		for (int yT = 0; yT < occurrenceCounts[xT].size(); yT++) {
			if (!sampleSheet.bc2uniquePair[xT][yT].second) {
				cout << "Barcode " << sampleSheet.bc2uniquePair[xT][yT].first << ": " << occurrenceCounts[xT][yT] << endl;
			}
		}
	}
}

seqan::ArgumentParser::ParseResult ArgumentParse(int argc, char const ** argv) {
	// Setup ArgumentParser.
	seqan::ArgumentParser parser("barcode_extract");

	// set some documentation
	seqan::setShortDescription(parser, "Extract barcodes from NGS");
	seqan::setVersion(parser, "1.0");
	seqan::setDate(parser, "June 2018");

	/*
		Formatting Command Line Documentation
	*/
	seqan::addUsageLine(parser,
		"\"\\fBR1 File\\fP\" \"\\fBR2 File\\fP\" -s \"samplesheet.txt\" [\\fIOPTIONS\\fP]");
	seqan::addDescription(parser,
		"This program splits NGS reads based on barcodes found on the 5' and 3' read ends");

	seqan::addTextSection(parser, "Examples");

	seqan::addListItem(parser,
		"\\fBbarcode_extract\\fP \\fImerged_fastq_file.fq.gz\\fP -s sample_sheet_location.txt",
		"Analyze an unpaired FASTQ file"
	);

	seqan::addListItem(parser,
		"\\fBbarcode_extract\\fP \\fIr1.fq.gz r2.fq.gz\\fP -s sample_sheet_location.txt -p5 5 -p3 3",
		"Analyze a paired R1/R2 read where barcode starts 5bp upstream of read and 3bp downstream of the end of the read/start of R2"
	);

	// set the positional argument to be a list of fastq files (islist = true)
	addArgument(parser, seqan::ArgParseArgument(
		seqan::ArgParseArgument::INPUT_FILE, "FASTQ FILE(S)", true));

	// first argument is location of tab delimited file
	addOption(parser, seqan::ArgParseOption(
		"s", "samplesheet", "Location of TAB delimited file for deconvoluting barcodes",
		seqan::ArgParseArgument::INPUT_FILE, "SAMPLE"));
	// required 
	seqan::setRequired(parser, "s");
	seqan::setValidValues(parser, "s", "txt");

	// define edit distance
	seqan::addOption(parser, seqan::ArgParseOption(
		"e", "editdistance", "Maximum edit distance",
		seqan::ArgParseArgument::INTEGER, "EDIT")
	);
	seqan::setMinValue(parser, "e", "0");
	seqan::setDefaultValue(parser, "e", "0");

	// define where to start searching for the 5' barcode length with respect to the start of an R1 read
	seqan::addOption(parser, seqan::ArgParseOption(
		"p5", "us-start", "Start position to look for 5' barcode",
		seqan::ArgParseArgument::INTEGER, "USSTART")
	);
	seqan::setMinValue(parser, "p5", "0");
	seqan::setDefaultValue(parser, "p5", "0");

	// define where to start searching for the 3' barcode length with respect to the start of an R2 read or END of an R1 read (distance from end)
	seqan::addOption(parser, seqan::ArgParseOption(
		"p3", "ds-start", "Start position (from end) to look for 3' barcode",
		seqan::ArgParseArgument::INTEGER, "DSSTART")
	);
	seqan::setMinValue(parser, "p3", "0");
	seqan::setDefaultValue(parser, "p3", "0");

	// define the range of the search (default of 0 will be determined based on minium barcode length)
	seqan::addOption(parser, seqan::ArgParseOption(
		"l5", "us-start-range", "How many bases from start position can the barcode be found (value of 0 will assume distance must be length of barcode provided)",
		seqan::ArgParseArgument::INTEGER, "USRANGE")
	);
	seqan::setMinValue(parser, "l5", "0");
	seqan::setDefaultValue(parser, "l5", "0");

	// define the range of the search (default of 0 will be determined based on minium barcode length)
	seqan::addOption(parser, seqan::ArgParseOption(
		"l3", "ds-start-range", "How many bases from start position can the 3' barcode be found (value of 0 will assume distance must be length of barcode provided)",
		seqan::ArgParseArgument::INTEGER, "DSRANGE")
	);
	seqan::setMinValue(parser, "l3", "0");
	seqan::setDefaultValue(parser, "l3", "0");

	// define the number of fastq reads we want in each file
	seqan::addOption(parser, seqan::ArgParseOption(
		"r", "num-reads", "Defines the maximum number of reads a fastq file is allowed to have (value of 0 assumes unlimited # of reads allowed)",
		seqan::ArgParseArgument::INTEGER, "MAX_NUM_READS")
	);
	seqan::setMinValue(parser, "r", "0");
	seqan::setDefaultValue(parser, "r", "0");

	// define the number of fastq reads we want in each file
	seqan::addOption(parser, seqan::ArgParseOption(
		"b", "batch-size", "Defines the reads to read at once",
		seqan::ArgParseArgument::INTEGER, "BATCH_SIZE")
	);
	seqan::setMinValue(parser, "r", "1");
	seqan::setDefaultValue(parser, "r", "10000");

	// define the number of fastq reads we want in each file
	seqan::addOption(parser, seqan::ArgParseOption(
		"t", "total-reads", "Force program to stop after reading this many reads (value of 0 will not have a limit on total number of reads)",
		seqan::ArgParseArgument::INTEGER, "TOTAL_READS")
	);
	seqan::setMinValue(parser, "t", "0");
	seqan::setDefaultValue(parser, "t", "0");

	// first argument is location of tab delimited file
	addOption(parser, seqan::ArgParseOption(
		"x", "prefix", "Prefix used for saving output files to a new location (if empty then the prefix will match the filename of the input files)",
		seqan::ArgParseArgument::STRING, "PREFIX"));
	// required 
	seqan::setDefaultValue(parser, "x", "");
	

	// define list of valid values
	// seqan::addOption(parser, seqan::ArgParseOption(
	//	"", "distance-model", "Distance model, either HAMMING or EDIT.",
	//	seqan::ArgParseArgument::STRING, "STR"));
	//seqan::setValidValues(parser, "distance-model", "HAMMING EDIT");

	// defining custom input/output file names
	//seqan::addOption(parser, seqan::ArgParseOption(
	//	"I", "input-file", "Path to the input file",
	//	seqan::ArgParseArgument::INPUT_FILE, "IN"));
	//seqan::addOption(parser, seqan::ArgParseOption(
	//	"O", "output-file", "Path to the output file",
	//	seqan::ArgParseArgument::OUTPUT_FILE, "OUT"));
	// allowed file ext
	//seqan::setValidValues(parser, "input-file", "txt");

	// Parse command line.
	seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

	// If parsing was not successful then exit with code 1 if there were errors.
	// Otherwise, exit with code 0 (e.g. help was printed).
	if (res != seqan::ArgumentParser::PARSE_OK)
		return res;
	
	// make sure the number of read files is <= 2
	unsigned int numFiles = seqan::getArgumentValueCount(parser, 0);
	if (numFiles > 2) {
		std::cerr << "ERROR: You have past too many input files in. Either pass in a single file for unpaired, or two files for paired sequences" << endl;
		return seqan::ArgumentParser::PARSE_ERROR;
	}

	if (numFiles == 1) {
		paired = false;
		seqan::getArgumentValue(input_r1_read, parser, 0, 0);
	}
	else {
		paired = true;
		seqan::getArgumentValue(input_r1_read, parser, 0, 0);
		seqan::getArgumentValue(input_r2_read, parser, 0, 1);
	}

	// extract option values
	seqan::getOptionValue(desired_edit_distance, parser, "editdistance");
	seqan::getOptionValue(sample_sheet, parser, "s");
	seqan::getOptionValue(us_start_search, parser, "p5");
	seqan::getOptionValue(ds_start_search, parser, "p3");
	seqan::getOptionValue(us_search_range, parser, "l5");
	seqan::getOptionValue(ds_search_range, parser, "l3");
	seqan::getOptionValue(max_counts_per_file, parser, "r");
	seqan::getOptionValue(batch_read, parser, "b");
	seqan::getOptionValue(stop_after, parser, "t");
	seqan::getOptionValue(prefix, parser, "x");

	 return seqan::ArgumentParser::PARSE_OK;

}

int main(int argc, char const ** argv) {
	// parse input arguments
	seqan::ArgumentParser::ParseResult res = ArgumentParse(argc, argv);
	// If parsing was not successful then exit with code 1 if there were errors.
	// Otherwise, exit with code 0 (e.g. help was printed).
	if (res != seqan::ArgumentParser::PARSE_OK)
		return res == seqan::ArgumentParser::PARSE_ERROR;


	SampleList barcodeInfo = ParseSampleSheet(sample_sheet);
	DnaList us_barcode, ds_barcode;

	std::pair<seqan::CharString, seqan::CharString> R1FileInfo, R2FileInfo;

	if (paired) {
		R1FileInfo = ExtractFileNamePrefix(input_r1_read);
		R2FileInfo = ExtractFileNamePrefix(input_r2_read);
	}
	else {
		R1FileInfo = ExtractFileNamePrefix(input_r1_read);
	}

	if (prefix != "") {
		//user has predefined the location and prefix for output files
		R1FileInfo.first = prefix;
		R1FileInfo.first += "_R1_001";
		R2FileInfo.first = prefix;
		R2FileInfo.first += "_R2_001";
	}

	for (auto bc : barcodeInfo.usBc) {
		if (bc == "") continue;
		seqan::appendValue(us_barcode, bc);
	}

	for (auto bc : barcodeInfo.dsBc) {
		if (bc == "") continue;
		DnaSeq tmp = DnaSeq(bc);
		if (paired) {
			// R1 and R2 files were provided, so we can search the RC sequence of the ds barcode in the R2 file
			seqan::appendValue(ds_barcode, tmp);
		}
		else {
			// only one file provided, so we assume that the ds barcode is present in its RC direction
			seqan::appendValue(ds_barcode, seqan::Dna5StringReverseComplement(tmp));
		}
	}
	barcodeInfo.MapBC2Sample(us_barcode, ds_barcode);
	seqan::SeqFileIn seqFileR1In, seqFileR2In;

	if (!seqan::open(seqFileR1In, seqan::toCString(input_r1_read))) {
		std::cerr << "ERROR: could not open input R1 file.\n";
		std::cerr << input_r1_read << "\n";
		return 1;
	}

	if (!seqan::open(seqFileR2In, seqan::toCString(input_r2_read))) {
		std::cerr << "ERROR: could not open input R2 file.\n";
		std::cerr << input_r2_read << "\n";
		return 1;
	}

	if (paired) {
		SearchFastqFile(seqFileR1In, seqFileR2In, us_barcode, ds_barcode, barcodeInfo, R1FileInfo, R2FileInfo);
	}
	else {
		SearchFastqFile(seqFileR1In, us_barcode, ds_barcode, barcodeInfo, R1FileInfo);
	}
	
	return 0;
}
