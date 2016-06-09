//
//  CompareVCFs.h
//  VCF_Comparisons
//
//  Created by Pelin Sahlen on 16/06/2015.
//  Copyright (c) 2015 Pelin Sahlen. All rights reserved.
//

#ifndef VCF_Comparisons_CompareVCFs_h
#define VCF_Comparisons_CompareVCFs_h

class CompareVCFs{
friend class ProcessVCF;
public:
    boost::unordered::unordered_map<string, int> sortedvartype;
    void CompareVariants(ProcessVCF,int);
    void CompareTotalReadDepth(ProcessVCF, string, string, int, int, double&);
    void CompareFlagValues(ProcessVCF, string, string, int, int, double&, string);
    void SortVariantTypes(void);
private:
    void GetReadDepths(ProcessVCF, int,int);
    void FillHistogram(vector <long unsigned int>&, vector<double>, int, int);
    void CountDisconcordantVariants(ProcessVCF);
    void DetermineVariantTypes(ProcessVCF, string, string, int, int);
    
};
void CompareVCFs::SortVariantTypes(void){
    
    sortedvartype["AC"] = 0;
    sortedvartype["AG"] = 1;
    sortedvartype["AT"] = 2;
    sortedvartype["CA"] = 3;
    sortedvartype["CG"] = 4;
    sortedvartype["CT"] = 5;
    sortedvartype["GA"] = 6;
    sortedvartype["GC"] = 7;
    sortedvartype["GT"] = 8;
    sortedvartype["TA"] = 9;
    sortedvartype["TC"] = 10;
    sortedvartype["TG"] = 11;
    
}
void CompareVCFs::DetermineVariantTypes(ProcessVCF metavcf, string first, string second, int v1, int v2){
    
    boost::unordered::unordered_map<string, vector < variantfeatures > >::iterator it;
    
    boost::unordered::unordered_map<int, double> vartypes;
    
    string filename;
    filename.append("VariantTypes.");
    filename.append(first);
    filename.append("_");
    filename.append(second);
    filename.append(".txt");
    ofstream outf(filename.c_str());
    
    for (it = metavcf.variantmap.begin(); it != metavcf.variantmap.end(); ++it){
        if (it->second[0].snp_type == v1 && it->second[1].snp_type == v2) {
            string s = metavcf.common_feat_map[it->first].ref;
            s.append(metavcf.common_feat_map[it->first].alt);
            if (vartypes.find(sortedvartype[s]) == vartypes.end())
                    vartypes[sortedvartype[s]] = 0;
                else
                    ++vartypes[sortedvartype[s]];
        }
    }
    
    
    boost::unordered::unordered_map<int, double>::iterator it2;
    for (int i = 0; i < sortedvartype.size(); ++i) {
        for (it2 = vartypes.begin(); it2 != vartypes.end(); ++it2){
            if(it2->first == i){
                outf << it2->first << '\t' << it2->second << endl;
            }
        }
    }
    outf.close();
}
void CompareVCFs::CompareVariants(ProcessVCF metavcf, int nof_vcfs){
    
    double NofCommonKnownSNPs = 0, NofCommonNovelSNPs = 0;
    double NofUnique1_KnownSNPs = 0, NofUnique2_KnownSNPs = 0;
    double NofUnique1_NovelSNPs = 0, NofUnique2_NovelSNPs = 0;
    
    boost::unordered::unordered_map<string, vector < variantfeatures > >::iterator it;
    
    for (it = metavcf.variantmap.begin(); it != metavcf.variantmap.end(); ++it){
        if (it->second[0].snp_type == 0 && it->second[1].snp_type == 1 ) {
            it->second[1].snp_type = 2;
        }
        if (it->second[0].snp_type == 1 && it->second[1].snp_type == 0 ) {
            it->second[0].snp_type = 2;
        }
    }
/*
    SortVariantTypes();
    DetermineVariantTypes(metavcf, "known", "known", 0, 0);
    DetermineVariantTypes(metavcf, "novel", "novel", 1, 1);
    DetermineVariantTypes(metavcf, "known", "notpresent", 0, 2);
    DetermineVariantTypes(metavcf, "notpresent", "known", 2, 0);
    DetermineVariantTypes(metavcf, "novel", "notpresent", 1, 2);
    DetermineVariantTypes(metavcf, "notpresent", "novel", 2, 1);
*/

    CountDisconcordantVariants(metavcf);

    CompareTotalReadDepth(metavcf, "known", "known", 0, 0, NofCommonKnownSNPs);
    CompareTotalReadDepth(metavcf, "novel", "novel", 1, 1, NofCommonNovelSNPs);
    CompareTotalReadDepth(metavcf, "known", "notpresent", 0, 2, NofUnique1_KnownSNPs);
    CompareTotalReadDepth(metavcf, "notpresent", "known", 2, 0, NofUnique2_KnownSNPs);
    CompareTotalReadDepth(metavcf, "novel", "notpresent", 1, 2, NofUnique1_NovelSNPs);
    CompareTotalReadDepth(metavcf, "notpresent", "novel", 2, 1, NofUnique2_NovelSNPs);

/*
    CompareFlagValues(metavcf, "known", "known", 0, 0, NofCommonKnownSNPs, "DP");
    CompareFlagValues(metavcf, "novel", "novel", 1, 1, NofCommonNovelSNPs, "DP");
    CompareFlagValues(metavcf, "known", "notpresent", 0, 2, NofUnique1_KnownSNPs, "DP");
    CompareFlagValues(metavcf, "notpresent", "known", 2, 0, NofUnique2_KnownSNPs, "DP");
    CompareFlagValues(metavcf, "novel", "notpresent", 1, 2, NofUnique1_NovelSNPs, "DP");
    CompareFlagValues(metavcf, "notpresent", "novel", 2, 1, NofUnique2_NovelSNPs, "DP");
    
    CompareFlagValues(metavcf, "known", "known", 0, 0, NofCommonKnownSNPs, "GQ");
    CompareFlagValues(metavcf, "novel", "novel", 1, 1, NofCommonNovelSNPs, "GQ");
    CompareFlagValues(metavcf, "known", "notpresent", 0, 2, NofUnique1_KnownSNPs, "GQ");
    CompareFlagValues(metavcf, "notpresent", "known", 2, 0, NofUnique2_KnownSNPs, "GQ");
    CompareFlagValues(metavcf, "novel", "notpresent", 1, 2, NofUnique1_NovelSNPs, "GQ");
    CompareFlagValues(metavcf, "notpresent", "novel", 2, 1, NofUnique2_NovelSNPs, "GQ");
*/
    cout.setf(ios_base::fixed);
    ofile.setf(ios_base::fixed);
    
    double n1, n2, n3, n4, n5, n6;
    n1 = ((NofCommonKnownSNPs + NofCommonNovelSNPs)/metavcf.vcfstats[0].nof_snps)*100.0;
    n2 = ((NofCommonKnownSNPs + NofCommonNovelSNPs)/metavcf.vcfstats[1].nof_snps)*100.0;
    n3 = (NofCommonKnownSNPs/metavcf.vcfstats[0].nof_dbsnps)*100.0;
    n4 = (NofCommonKnownSNPs/metavcf.vcfstats[1].nof_dbsnps)*100.0;
    n5 = (NofCommonNovelSNPs/metavcf.vcfstats[0].nof_novelsnps)*100.0;
    n6 = (NofCommonNovelSNPs/metavcf.vcfstats[1].nof_novelsnps)*100.0;
    ofile << "Number of Common SNPs      " << '\t' << NofCommonKnownSNPs + NofCommonNovelSNPs << '\t' << n1 << '\t' << n2 << endl
         << "Number of Common Known SNPs" << '\t' << NofCommonKnownSNPs << '\t' << n3 << '\t' << n4 << endl
         << "Number of Common Novel SNPs" << '\t' << NofCommonNovelSNPs << '\t' << n5 << '\t' << n6 << endl
         << "Number of Known SNPs only detected in first"  << '\t' << NofUnique1_KnownSNPs << endl
         << "Number of Known SNPs only detected in second" << '\t' << NofUnique2_KnownSNPs << endl
         << "Number of Novel SNPs only detected in first"  << '\t' << NofUnique1_NovelSNPs << endl
         << "Number of Novel SNPs only detected in second" << '\t' << NofUnique2_NovelSNPs << endl;

   cout << "Number of Common SNPs      " << '\t' << NofCommonKnownSNPs + NofCommonNovelSNPs << '\t' << n1 << '\t' << n2 << endl
        << "Number of Common Known SNPs" << '\t' << NofCommonKnownSNPs << '\t' << n3 << '\t' << n4 << endl
        << "Number of Common Novel SNPs" << '\t' << NofCommonNovelSNPs << '\t' << n5 << '\t' << n6 << endl
        << "Number of Known SNPs only detected in first"  << '\t' << NofUnique1_KnownSNPs << endl
        << "Number of Known SNPs only detected in second" << '\t' << NofUnique2_KnownSNPs << endl
        << "Number of Novel SNPs only detected in first"  << '\t' << NofUnique1_NovelSNPs << endl
        << "Number of Novel SNPs only detected in second" << '\t' << NofUnique2_NovelSNPs << endl;


}
void CompareVCFs::CompareFlagValues(ProcessVCF metavcf, string first, string second, int v1, int v2, double& snpcount, string flag){
    
    boost::unordered::unordered_map<string, vector < variantfeatures > >::iterator it;
    
    vector<double> values1;
    vector<double> values2;
    
    
    for (it = metavcf.variantmap.begin(); it != metavcf.variantmap.end(); ++it){
        if (it->second[0].snp_type == v1 && it->second[1].snp_type == v2) {
            values1.push_back(it->second[0].flag_values[flag]);
            values2.push_back(it->second[1].flag_values[flag]);
            ++snpcount;
        }
    }
    int bin_size = 1;
    int number_of_bins = (int) ceil(1500 / bin_size);
    vector <long unsigned int> histogram1(number_of_bins);
    vector <long unsigned int> histogram2(number_of_bins);
    
    FillHistogram(histogram1, values1, number_of_bins, bin_size);
    FillHistogram(histogram2, values2, number_of_bins, bin_size);
    
    double sum1 = 0, sum2 = 0;
    for (int i = 0; i < histogram1.size(); ++i)
        sum1 += histogram1[i];
    for (int i = 0; i < histogram2.size(); ++i)
        sum2 += histogram2[i];
    
    string filename = "Compare";
    filename.append(flag);
    filename.append("_");
    filename.append(first);
    filename.append("_");
    filename.append(second);
    filename.append(".txt");
    
    ofstream outf(filename.c_str());
    outf << vcfnames[0] << '\t' << "count" << '\t' << "count/sum" << '\t'
    << vcfnames[1] << '\t' << "count" << '\t' << "count/sum" << endl;
    for (int i = 0; i < histogram1.size(); ++i){
        outf << i*bin_size << '\t' << histogram1[i] << '\t' << double (histogram1[i]/sum1) << '\t';
        if (i < histogram2.size()) {
            outf << histogram2[i] << '\t' << double (histogram2[i]/sum1) << endl;
        }
        else{
            outf << 0 << '\t' << 0 << endl;
        }
    }
    if (histogram2.size() > histogram1.size()) {
        for (unsigned long int i = histogram1.size(); i < histogram2.size(); ++i) {
            outf << i*bin_size << '\t' << 0 << '\t' << 0 << '\t';
            outf << histogram2[i] << '\t' << double (histogram2[i]/sum1) << endl;
        }
    }
    cout << filename << "   written" << endl;
}

void CompareVCFs::CompareTotalReadDepth(ProcessVCF metavcf, string first, string second, int v1, int v2, double& snpcount){
    
    boost::unordered::unordered_map<string, vector < variantfeatures > >::iterator it;
    
    vector<double> values1;
    vector<double> values2;

    
    for (it = metavcf.variantmap.begin(); it != metavcf.variantmap.end(); ++it){
        if (it->second[0].snp_type == v1 && it->second[1].snp_type == v2) {
                values1.push_back(it->second[0].ref_read_depth + it->second[0].alt_read_depth);
                values2.push_back(it->second[1].ref_read_depth + it->second[1].alt_read_depth);
                ++snpcount;
        }
    }
    int bin_size = 1;
    int number_of_bins = (int) ceil(1500 / bin_size);
    vector <long unsigned int> histogram1(number_of_bins);
    vector <long unsigned int> histogram2(number_of_bins);
    
    FillHistogram(histogram1, values1, number_of_bins, bin_size);
    FillHistogram(histogram2, values2, number_of_bins, bin_size);
    
    double sum1 = 0, sum2 = 0;
    for (int i = 0; i < histogram1.size(); ++i)
        sum1 += histogram1[i];
    for (int i = 0; i < histogram2.size(); ++i)
        sum2 += histogram2[i];
    
    string filename = "CompareReadDepth_";
    filename.append(first);
    filename.append("_");
    filename.append(second);
    filename.append(".txt");

    ofstream outf(filename.c_str());
    outf << vcfnames[0] << '\t' << "count" << '\t' << "count/sum" << '\t'
         << vcfnames[1] << '\t' << "count" << '\t' << "count/sum" << endl;
    for (int i = 0; i < histogram1.size(); ++i){
        outf << i*bin_size << '\t' << histogram1[i] << '\t' << double (histogram1[i]/sum1) << '\t';
        if (i < histogram2.size()) {
            outf << histogram2[i] << '\t' << double (histogram2[i]/sum1) << endl;
        }
        else{
            outf << 0 << '\t' << 0 << endl;
        }
    }
    if (histogram2.size() > histogram1.size()) {
        for (unsigned long int i = histogram1.size(); i < histogram2.size(); ++i) {
            outf << i*bin_size << '\t' << 0 << '\t' << 0 << '\t';
            outf << histogram2[i] << '\t' << double (histogram2[i]/sum1) << endl;
        }
    }
    cout << filename << "   written" << endl;
}
void CompareVCFs::FillHistogram(vector <long unsigned int>& histogram, vector<double> values, int number_of_bins, int bin_size){

    for (int i = 0; i < number_of_bins; ++i)
        histogram[i] = 0;
    
    vector<double>::iterator iter;
    for (iter = values.begin(); iter != values.end(); ++iter) {
        int bucket = (int)floor((*iter) / bin_size);
        if(bucket > number_of_bins){
            histogram.resize((bucket+1),0);
        }
        histogram[bucket] += 1;
    }
}
void CompareVCFs::CountDisconcordantVariants(ProcessVCF metavcf){
    
    ofstream outf("disconcordant.positions.txt");
    int nkd = 0, nnd = 0, nd = 0;
    
    boost::unordered::unordered_map<string, vector < variantfeatures > >::iterator it;
    for (it = metavcf.variantmap.begin(); it != metavcf.variantmap.end(); ++it){
        if (it->second[0].found && it->second[1].found) {
            if (it->second[0].altgenotype != it->second[1].altgenotype) {
                ++nd;
                if(it->second[0].snp_type == 1)
                    ++nkd;
                else
                    ++nnd;
                outf << it->first << "  " << it->second[0].snp_type << "  " << it->second[1].snp_type << "  " << it->second[0].altgenotype << "  " << it->second[1].altgenotype << endl;
            }
        }
    }
    ofile << "total disconcordant genotypes  " << nd << endl;
    ofile << "dbsnp disconcordant genotypes  " << nkd << endl;
    ofile << "novel disconcordant genotypes  " << nnd << endl;

}
#endif
