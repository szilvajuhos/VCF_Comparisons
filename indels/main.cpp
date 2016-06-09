//
//  main.cpp
//  VCF_Comparisons
//
//  Created by Pelin Sahlen on 15/06/2015.
//  Copyright (c) 2015 Pelin Sahlen. All rights reserved.
//

#include <stdio.h>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <cstring>
#include <string>
#include <time.h>
#include <vector>
#include <sstream>
#include <boost/unordered/unordered_map.hpp>
#include <boost/unordered/unordered_set.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
using namespace std;

#define BOOST_DISABLE_ASSERTS
#define NDEBUG

vector < string > vcfnames;
string truthsetfile;
string directory;

ofstream ofile("vcf_stats.indels.nofilter.txt");

#include "ProcessVCF.h"
#include "CompareVCFs.h"
#include "HapMap_truth_set.h"


int main(int argc, const char * argv[]) {
   
    vector < string > fname;
    ProcessVCF metavcf;
    
    CompareVCFs comparevcfs;
    
  //  truthsetfile = argv[1];
    truthsetfile = "/pica/v12/b2013064_nobackup/pelin/ANALYSIS/SwedishGenomePilot/Comparisons/genotypes_CEU_r28_nr.b37_fwd.filtered.sortedbyRs";
    vcfnames.push_back(argv[1]);
    vcfnames.push_back(argv[2]);
        
    int vcf_index = 0;
    for (int i; i < 2; ++i) {
        metavcf.readVCF(vcfnames[i], i);
        metavcf.GenerateVCFStats(vcfnames[i], i);
    }
//    TruthSet hapmapgenotypes;
    
//    hapmapgenotypes.CreateTruthSet(metavcf);
//    hapmapgenotypes.PrintConcordance(metavcf);

    comparevcfs.CompareVariants(metavcf,vcf_index);
    
    return 0;
}
