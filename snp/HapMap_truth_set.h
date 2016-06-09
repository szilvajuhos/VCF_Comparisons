
class TruthSet{
    friend class ProcessVCF;
public:
    boost::unordered::unordered_map<string, double> vartypes_validated;
    boost::unordered::unordered_map<string, double> vartypes_dbsnp;
    boost::unordered::unordered_map<string, double> vartypes_novel;
    void CreateTruthSet(ProcessVCF&);
    void PrintConcordance(ProcessVCF&);
    void PrintVariantTypes(ProcessVCF&);
private:
    bool processgenotype(string, string);
    
};

int WordOccurrenceCount( std::string const & str, std::string const & word )
{
    int count(0);
    std::string::size_type word_pos( 0 );
    while ( word_pos < word.length() )
    {
        word_pos = str.find(word, word_pos );
        if ( word_pos != std::string::npos )
        {
            ++count;
            
            // start next search after this word
            ++word_pos;
        }
    }
    
    return count;
}
void TruthSet::PrintVariantTypes(ProcessVCF& variant_set){
    
    boost::unordered::unordered_map<string, variant_common_features >::iterator iter;
    
    for(iter = variant_set.common_feat_map.begin(); iter != variant_set.common_feat_map.end();++iter){
        string s = iter->second.ref;
        s.append(iter->second.alt);
        if (iter->second.validated) {
            if (vartypes_validated.find(s) == vartypes_validated.end())
                vartypes_validated[s] = 0;
            else
                ++vartypes_validated[s];
        }
        if (iter->second.dbsnpid.length() > 1) { //dbSNP
            if (vartypes_dbsnp.find(s) == vartypes_dbsnp.end())
                vartypes_dbsnp[s] = 0;
            else
                ++vartypes_dbsnp[s];
        }
        else{
            if (vartypes_novel.find(s) == vartypes_novel.end())
                vartypes_novel[s] = 0;
            else
                ++vartypes_novel[s];
        }
    }
    
    boost::unordered::unordered_map<string, double>::iterator it;
    cout << "validated" << endl;
    for (it = vartypes_validated.begin(); it != vartypes_validated.end(); ++it) {
        cout << it->first << '\t' << it->second << endl;
    }
    cout << "dbsnp" << endl;
    for (it = vartypes_dbsnp.begin(); it != vartypes_dbsnp.end(); ++it) {
        cout << it->first << '\t' << it->second << endl;
    }
    cout << "novel" << endl;
    for (it = vartypes_novel.begin(); it != vartypes_novel.end(); ++it) {
        cout << it->first << '\t' << it->second << endl;
    }
    
}

bool TruthSet::processgenotype(string genotype, string alt){

    if ((genotype.find("N") == std::string::npos)) {
        int c = WordOccurrenceCount(genotype.c_str(), alt.c_str());
        if (c == 2)
            return false;
        else
            return true;
    }
    
    return false;
}

void TruthSet::CreateTruthSet(ProcessVCF& variant_set){
    

    ifstream genotype_file(truthsetfile.c_str());
    
    string temp, rsid, chr, refallele, gen[3], address;
    int pos;
    
    boost::unordered::unordered_map<string, string >::iterator rs_iter;
    boost::unordered::unordered_map<string, variant_common_features >::iterator iter;
//    ofstream outf("truthset.txt");
 
    getline(genotype_file,temp);
    
    genotype_file >> rsid >> temp >> chr >> pos >> temp >> temp >> temp >> temp >> temp >> temp >> gen[0] >> gen[1] >> gen[2];
    do{
        if(chr.find("GL") == string::npos && chr.find("MT") == string::npos){
     //       outf << rsid << "  " << chr << "  " << pos << " " << gen[0] << gen[1] << gen[2] << endl;
            
            rs_iter = variant_set.rsid_to_coord.find(rsid);
            if(rs_iter != variant_set.rsid_to_coord.end()){
                if(variant_set.common_feat_map.find(rs_iter->second) != variant_set.common_feat_map.end()){
                    iter = variant_set.common_feat_map.find(rs_iter->second);
                    variant_set.common_feat_map[rs_iter->second].validated = true;
                    for (int i = 0; i < 3; ++i) {
                        iter->second.truthsnp[i] = processgenotype(gen[i],iter->second.alt);
      //                  outf << gen[i] << '\t' << iter->second.ref << '\t' << iter->second.alt << '\t';
         //               if (iter->second.truthsnp[i])
           //                 outf << 1 << endl;
           //             else
           //                 outf << 0 << endl;
                        
                    }
                }
            }
        }
        genotype_file >> rsid >> temp >> chr >> pos >> temp >> temp >> temp >> temp >> temp >> temp >> gen[0] >> gen[1] >> gen[2];
            
    }while (!genotype_file.eof());
    
}

void TruthSet::PrintConcordance(ProcessVCF& metavcf){
    
    double val1 = 0, val2 = 0, total_val = 0, valboth = 0, valnone = 0;
    boost::unordered::unordered_map<string, variant_common_features >::iterator it;
    boost::unordered::unordered_map<string, vector < variantfeatures > >::iterator it2;

    
    for (it = metavcf.common_feat_map.begin(); it != metavcf.common_feat_map.end(); ++it){
        if (it->second.validated) {
            ++total_val;
            it2 = metavcf.variantmap.find(it->first);
            if (it2->second[0].found == 1){ // dbSNP
                    ++val1;
                if (it2->second[1].found == 1){
                    ++val2;
                    ++valboth;
                }
            }
            else{
                if (it2->second[1].found == 1)
                    ++val2;
                else
                    ++valnone;
            }
        }
    }
    
    cout << "Total number of validated SNPs" << '\t' << total_val << endl
         << "Found in 1" << '\t' << val1 << '\t' << double(val1/total_val) << endl
         << "Found in 2" << '\t' << val2 << '\t' << double(val2/total_val) << endl;
    cout << "Found in both" << '\t' << valboth << '\t' << double(valboth/total_val) << endl
         << "Found in none" << '\t' << valnone << '\t' << double(valnone/total_val) << endl;

    ofile << "Total number of validated SNPs" << '\t' << total_val << endl
          << "Found in 1" << '\t' << val1 << '\t' << double(val1/total_val) << endl
          << "Found in 2" << '\t' << val2 << '\t' << double(val2/total_val) << endl;
    ofile << "Found in both" << '\t' << valboth << '\t' << double(valboth/total_val) << endl
          << "Found in none" << '\t' << valnone << '\t' << double(valnone/total_val) << endl;

}

