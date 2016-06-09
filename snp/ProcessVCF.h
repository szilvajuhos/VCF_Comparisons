//
//  ProcessVCF.h
//  VCF_Comparisons
//
//  Created by Pelin Sahlen on 15/06/2015.
//  Copyright (c) 2015 Pelin Sahlen. All rights reserved.
//

#ifndef VCF_Comparisons_ProcessVCF_h
#define VCF_Comparisons_ProcessVCF_h

struct VCFStruct{
    double nof_snps;
    double nof_dbsnps;
    double nof_novelsnps;
    double transitions;
    double transversions;
};

struct variantfeatures{
    bool found;
    boost::unordered::unordered_map<string, double> flag_values;
    int alt_allele_count; // 0 = hom_alt, 1 = het
    int qual;
    double ref_read_depth;
    double alt_read_depth;
    string refgenotype;
    string altgenotype;
    string filter_value;
    int snp_type;  // 0 = known, 1 = novel, not_present = 2
};

struct variant_common_features{
    bool truthsnp[3]; //if HapMap genotype exist. NA10860 [0], NA11992 [1], NA11993 [3]
    string ref;
    string alt;
    string dbsnpid;
    bool validated;
};
class ProcessVCF{
public:
    //data structures
    boost::unordered::unordered_map<string, string> flags;
    boost::unordered::unordered_map<string, vector < variantfeatures > > variantmap; //index : chr_pos
    boost::unordered::unordered_map<string, variant_common_features > common_feat_map; //index: chr_pos
    boost::unordered::unordered_map<string, string > rsid_to_coord; //index: dbsnpid key: chr_pos
    // vector keeps the different vcf files

    void GetFlags(string);
    void readVCF(string, int); // vcf_file_name, vcf_index
    void GenerateVCFStats(string,int);
    vector< VCFStruct > vcfstats;
private:
    void readheader(ifstream&);
    void GetFlagsfromHeader(ifstream&, boost::iostreams::filtering_stream<boost::iostreams::input>&, bool);
    void fill_snp_map(string, double, int, string, double, bool, stringstream&, stringstream&,stringstream&);
    void AlleleCount(double,double&,double&,double&);
    void GenerateHistogram(int,string,int,int,int, string, string);
    void GenerateHistogram_flags(int, string, int, int, int, string, string);
    void FillHistogram(vector<long unsigned int>&,int, int, double);
    void TiTvRatio(int);
};


void ProcessVCF::GetFlagsfromHeader(ifstream& vcf_file, boost::iostreams::filtering_stream<boost::iostreams::input>& vcf_file_zipped, bool zipped){
    
    string temp;
    if(zipped)
        std::getline(vcf_file_zipped,temp);
    else
        getline(vcf_file,temp);
    while (temp.substr(0,2) == "##") {
        if(zipped)
            std::getline(vcf_file_zipped,temp);
        else
            getline(vcf_file,temp);
        while (temp.substr(0,8) == "##FORMAT") {
            unsigned long int p1 = temp.find("<ID");
            unsigned long int p2 = temp.find("=",p1);
            unsigned long int p3 = temp.find(",", p2);
            string flag = temp.substr(p2+1, p3-p2-1);
            
            p1 = temp.find("Description");
            p2 = temp.find("=",p1);
            p3 = temp.find(">");
            string description = temp.substr(p2+1,p3-p2-1);
            flags[flag] = description;
            if(zipped)
                std::getline(vcf_file_zipped,temp);
            else
                getline(vcf_file,temp);
        }
        while (temp.substr(0,6) == "##INFO") {
            unsigned long int p1 = temp.find("<ID");
            unsigned long int p2 = temp.find("=",p1);
            unsigned long int p3 = temp.find(",", p2);
            string flag = temp.substr(p2+1, p3-p2-1);
            
            p1 = temp.find("Description");
            p2 = temp.find("=",p1);
            p3 = temp.find(">");
            string description = temp.substr(p2+1,p3-p2-1);
            flags[flag] = description;
            if(zipped)
                std::getline(vcf_file_zipped,temp);
            else
                getline(vcf_file,temp);
        }
    }
}
void ProcessVCF::readheader(ifstream& vcf_file ){
    string temp;
    getline(vcf_file,temp);
    while (temp.substr(0,2) == "##") {
        getline(vcf_file,temp);
    }
    getline(vcf_file,temp); //read the column headers, the next line is SNP
}

void ProcessVCF::fill_snp_map(string chr_pos, double pos, int vcf_index, string dbsnp_id, double qual, bool filter, stringstream& info_stream, stringstream& formatdescription, stringstream& format){
    string field, flag,temp;
    double val;
    unsigned long int t;
    
    boost::unordered::unordered_map<string, double> flagmap;
    //Process INFO field
    do{
        getline(info_stream, field, ';');
        t = field.find("=");
        
        flag = (field.substr(0,(t)));
        temp = field.substr(t+1,field.length());
        val = atof(temp.c_str());
        
        flagmap[flag] = val;
    }while (!info_stream.eof());

    if(variantmap.find(chr_pos) == variantmap.end()){
        variantmap[chr_pos].push_back(variantfeatures());
        variantmap[chr_pos].push_back(variantfeatures());
        variantmap[chr_pos][0].found = false;
        variantmap[chr_pos][1].found = false;
        variantmap[chr_pos][0].snp_type = 2; // not_present state is default
        variantmap[chr_pos][1].snp_type = 2; // not_present state is default
        common_feat_map[chr_pos].validated = false;
        common_feat_map[chr_pos].truthsnp[0] = false;
        common_feat_map[chr_pos].truthsnp[1] = false;
        common_feat_map[chr_pos].truthsnp[2] = false;

    }
    //Add the flags from the info field
    //PROCESS FORMAT FIELD
    string value;
    do{
        getline(formatdescription, field, ':');
        t = field.find(":");
        flag = (field.substr(0,(t)));
    
        getline(format,value,':');
        if(value.find("/") == std::string::npos && value.find(",") == std::string::npos){
            val = atof(value.c_str());
            flagmap[flag] = val;
        }
        else{
            if(flag == "GT"){
                if (value == "0/1")
                    variantmap[chr_pos][vcf_index].alt_allele_count = 1;
                if (value == "1/1") {
                    variantmap[chr_pos][vcf_index].alt_allele_count = 2;
                }
            }
            if (flag == "AD") {
                unsigned long int a = value.find(",");
                string refdepth = value.substr(0,a);
                string altdepth = value.substr(a+1,value.size());
                
                variantmap[chr_pos][vcf_index].ref_read_depth = atof(refdepth.c_str());
                variantmap[chr_pos][vcf_index].alt_read_depth = atof(altdepth.c_str());
            }
        }
      //  cout << flag << "  " << value << "  " << variantmap[chr_pos][vcf_index].alt_allele_count << endl;
    }while (!formatdescription.eof());
    
    boost::unordered::unordered_map<string, double>::iterator iter;
    for (iter = flagmap.begin(); iter!= flagmap.end(); ++iter) {
        variantmap[chr_pos][vcf_index].flag_values[iter->first] = iter->second;
    }
    
}
void ProcessVCF::readVCF(string vcf_file_name, int vcf_index){
    
   
    ifstream vcf_file;
    boost::iostreams::filtering_stream<boost::iostreams::input> in;
    bool zipped = 0;
    if (vcf_file_name.find(".gz") == std::string::npos) {  // if not a gzip file
        vcf_file.open(vcf_file_name.c_str());
        cout << vcf_file_name << "  opened" << endl;
        if (!vcf_file.is_open())
            cerr << " File cannot be opened" << endl;
    }
    else{
        vcf_file.open(vcf_file_name.c_str(), std::ios_base::in | std::ios_base::binary);
        zipped = 1;
        try {
            in.push(boost::iostreams::gzip_decompressor());
            in.push(vcf_file);
            cout << vcf_file_name << "  opened" << endl;
        } catch (const boost::iostreams::gzip_error& e) {
            std::cout << e.what() << endl;
        }
    }
    
    string temp1;
    //vcf_file >> temp1;
    getline(vcf_file, temp1);
    GetFlagsfromHeader(vcf_file, in,zipped);
    
    int pos;
    double qual;
    string chr, qualstr, dbsnp_id, ref, alt, formatdescription, filter, info, format;
    bool passed;
    vcfstats.push_back(VCFStruct());
    vcfstats[vcf_index].nof_snps = 0;
    vcfstats[vcf_index].nof_dbsnps = 0;
    vcfstats[vcf_index].nof_novelsnps = 0;
    
//Read SNP lines
    int linecount = 0, n = 0;
    cout.setf(ios_base::fixed);
    do{
   //     if(zipped)
     //       in >> chr >> pos >> dbsnp_id >> ref >> alt >> qualstr >> filter >> info >> formatdescription >> format;
      //  else
            vcf_file >> chr >> pos >> dbsnp_id >> ref >> alt >> qualstr >> filter >> info >> formatdescription >> format;
        ++linecount;
  //      cout << chr << "  " << pos << "  " << filter << endl;
        /*
        if (chr.find("chr") == std::string::npos){
            string chrstr = "chr";
            string newchr = chrstr.append(chr);
            chr.clear();
            chr = newchr;
        }
        */
  //      if (chr != "1" && vcfstats[vcf_index].nof_snps > 0)
    //        break;
        //chr.find("20") != string::npos &&
        string filter_val = "."; //edico
        if (chr.find("chr") == std::string::npos){
          //  filter_val = "PASS"; // piper
            string chrstr = "chr";
            string newchr = chrstr.append(chr);
            chr.clear();
            chr = newchr;
//            n = chr.find("chr");
  //          chr = chr.substr(3, chr.length());
        }
        if (chr.find("GL") == string::npos && chr.find("MT") == string::npos && ref.size() == 1 && alt.size() == 1) { // if it is a variant site and a SNP // NO filter
            string address = "";
            address.append(chr);
            address.append("_");
            string s = boost::lexical_cast< string >(pos); // convert coordinate to string
            address.append(s); //chr_coordinate
  //          cout << chr << "  " << pos << "   " << address << endl;
            stringstream infostream;
            infostream << info;
        
            stringstream formatdesc;
            formatdesc << formatdescription;
        
            stringstream formatstream;
            formatstream << format;

            double postemp = pos;
            qual = atof(qualstr.c_str());
            fill_snp_map(address, postemp, vcf_index, dbsnp_id, qual, passed, infostream, formatdesc, formatstream);
            ++vcfstats[vcf_index].nof_snps;
            variantmap[address][vcf_index].filter_value = filter;
            variantmap[address][vcf_index].refgenotype = ref;
            variantmap[address][vcf_index].altgenotype = alt;
            variantmap[address][vcf_index].qual = qual;
            variantmap[address][vcf_index].found = true;
            
            if (dbsnp_id.length() > 1) { // known
                ++vcfstats[vcf_index].nof_dbsnps;
                variantmap[address][vcf_index].snp_type = 0; //known
                common_feat_map[address].dbsnpid = dbsnp_id;
                common_feat_map[address].ref = ref;
                common_feat_map[address].alt = alt;
                if(rsid_to_coord.find(dbsnp_id) == rsid_to_coord.end())
                    rsid_to_coord[dbsnp_id] = address;
            }
            else{
                ++vcfstats[vcf_index].nof_novelsnps;
                common_feat_map[address].dbsnpid = ".";
                common_feat_map[address].ref = ref;
                common_feat_map[address].alt = alt;
                variantmap[address][vcf_index].snp_type = 1; //novel
            }
            //cout << dbsnp_id << "  " << common_feat_map[address].dbsnpid << "  " << rsid_to_coord[dbsnp_id] << "  " << address << endl;
        }
            if(linecount % 250000 == 0)
                cout << linecount << "    lines read and  " << vcfstats[vcf_index].nof_snps << "  snps filled in the map" << endl;
    }while (!vcf_file.eof());
    
    cout << vcf_file_name << "   Read   " << vcfstats[vcf_index].nof_snps << "  SNPs processed" << endl;
/*
    string ofile;
    ofile.append(vcf_file_name);
    ofile.append(".fields.txt");
    ofstream vcf_out(ofile.c_str());
    
    boost::unordered::unordered_map<string, vector < variantfeatures > >::iterator vit;
    
    for (vit = variantmap.begin(); vit != variantmap.end(); ++vit) {
        vcf_out << vit->first << '\t' << vit->second[vcf_index].filter_value << '\t' << vit->second[vcf_index].altgenotype << '\t' << vit->second[vcf_index].alt_allele_count
        << '\t' << vit->second[vcf_index].alt_read_depth << '\t';
        if(vit->second[vcf_index].found == true)
            vcf_out << 1 << '\t';
        else
            vcf_out << 0 << '\t';
        vcf_out << common_feat_map[vit->first].dbsnpid << '\t' << common_feat_map[vit->first].ref << '\t';
        if(common_feat_map[vit->first].validated == true)
            vcf_out << 1 << '\t';
        else
            vcf_out << 0 << '\t';
        
        for (int z = 0; z<3;++z){
            if(common_feat_map[vit->first].truthsnp[z] == true)
                vcf_out << 1 << '\t';
            else
                vcf_out << 0 << '\t';
        }
        vcf_out << rsid_to_coord[common_feat_map[vit->first].dbsnpid] << endl;
    }
*/
}

void ProcessVCF::GenerateVCFStats(string vcf_name, int vcf_index){

    boost::unordered::unordered_map<string, vector < variantfeatures >>::iterator it;
    
    double nofhets[2], nofhomalts[2], nofmultiallelicsnps[2];
    
    nofhets[0] = 0; nofhomalts[0] = 0; nofmultiallelicsnps[0] = 0; //dbSNP
    nofhets[1] = 0; nofhomalts[1] = 0; nofmultiallelicsnps[1] = 0; //Novel
    
    for (it = variantmap.begin(); it != variantmap.end(); ++it) {
        if(it->second[vcf_index].found){
            if (it->second[vcf_index].snp_type == 0) //if dbSNP
                AlleleCount(it->second[vcf_index].alt_allele_count, nofhets[0], nofhomalts[0], nofmultiallelicsnps[0]);
            if (it->second[vcf_index].snp_type == 1)// if novel
                AlleleCount(it->second[vcf_index].alt_allele_count, nofhets[1], nofhomalts[1], nofmultiallelicsnps[1]);
        }
    }
    cout.setf(ios_base::fixed);
    ofile.setf(ios_base::fixed);
    
    ofile << vcf_name << endl;
    ofile << "Number of SNPs" << '\t' << "Number of dbSNPs" << '\t' << "Number of Novel SNPs" << endl;
    ofile << vcfstats[vcf_index].nof_snps << '\t' << vcfstats[vcf_index].nof_dbsnps << '\t' << vcfstats[vcf_index].nof_novelsnps << endl;
    ofile << "Number of het dbSNPs" << '\t' << "Number of hom alt dbSNPs" << endl;
    ofile << nofhets[0] << '\t' << nofhomalts[0]  << endl;
    ofile << "Number of het Novel SNPs" << '\t' << "Number of hom alt Novel SNPs" << endl;
    ofile << nofhets[1] << '\t' << nofhomalts[1] << endl;
    ofile << "het/hom ratio" << '\t' << (nofhets[0] + nofhets[1])/(nofhomalts[0] + nofhomalts[1]) << endl;

    
    cout << "Number of SNPs" << '\t' << "Number of dbSNPs" << '\t' << "Number of Novel SNPs" << endl;
    cout << vcfstats[vcf_index].nof_snps << '\t' << vcfstats[vcf_index].nof_dbsnps << '\t' << vcfstats[vcf_index].nof_novelsnps << endl;
    cout << "Number of het dbSNPs" << '\t' << "Number of hom alt dbSNPs" << endl;
    cout << nofhets[0] << '\t' << nofhomalts[0]  << endl;
    cout << "Number of het Novel SNPs" << '\t' << "Number of hom alt Novel SNPs" << endl;
    cout << nofhets[1] << '\t' << nofhomalts[1] << endl;
    cout << "het/hom ratio" << '\t' << (nofhets[0] + nofhets[1])/(nofhomalts[0] + nofhomalts[1]) << endl;

    GenerateHistogram_flags(2000,vcf_name,vcf_index,0,1,"GQ","MeanGenotypeQuality_dbSNP_het.txt");
    GenerateHistogram_flags(2000,vcf_name,vcf_index,0,2,"GQ","MeanGenotypeQuality_dbSNP_homalt.txt");
    GenerateHistogram_flags(2000,vcf_name,vcf_index,1,1,"GQ","MeanGenotypeQuality_Novel_het.txt");
    GenerateHistogram_flags(2000,vcf_name,vcf_index,1,2,"GQ","MeanGenotypeQuality_Novel_homalt.txt");
    

    GenerateHistogram(1500,vcf_name,vcf_index,0,1,"AD_ref","ReadDepthRef_dbSNP_het.txt");
    GenerateHistogram(1500,vcf_name,vcf_index,0,2,"AD_ref","ReadDepthRef_dbSNP_homalt.txt");
    GenerateHistogram(1500,vcf_name,vcf_index,1,1,"AD_ref","ReadDepthRef_Novel_het.txt");
    GenerateHistogram(1500,vcf_name,vcf_index,1,2,"AD_ref","ReadDepthRef_Novel_homalt.txt");
    
    GenerateHistogram(1500,vcf_name,vcf_index,0,1,"AD_alt","ReadDepthAlt_dbSNP_het.txt");
    GenerateHistogram(1500,vcf_name,vcf_index,0,2,"AD_alt","ReadDepthAlt_dbSNP_homalt.txt");
    GenerateHistogram(1500,vcf_name,vcf_index,1,1,"AD_alt","ReadDepthAlt_Novel_het.txt");
    GenerateHistogram(1500,vcf_name,vcf_index,1,2,"AD_alt","ReadDepthAlt_Novel_homalt.txt");
    
    GenerateHistogram(1500,vcf_name,vcf_index,0,1,"qual","SNPqual_dbSNP_het.txt");
    GenerateHistogram(1500,vcf_name,vcf_index,0,2,"qual","SNPqual_dbSNP_homalt.txt");
    GenerateHistogram(1500,vcf_name,vcf_index,1,1,"qual","SNPqual_Novel_het.txt");
    GenerateHistogram(1500,vcf_name,vcf_index,1,2,"qual","SNPqual_Novel_homalt.txt");
    
    TiTvRatio(vcf_index);
    cout  << "Ti/Tv" << '\t' << vcfstats[vcf_index].transitions/vcfstats[vcf_index].transversions << endl;
    ofile << "Ti/Tv" << '\t' << vcfstats[vcf_index].transitions/vcfstats[vcf_index].transversions << endl;
}
void ProcessVCF::AlleleCount(double allelecount, double& nofhets, double& nofhomalts, double& nofmultiallelicsnps){
    
    if(allelecount == 1)
        ++nofhets;
    else{
        if(allelecount == 2)
            ++nofhomalts;
        else
            ++nofmultiallelicsnps;
    }
}
void ProcessVCF::GenerateHistogram_flags(int maxval, string vcf_name, int vcf_index, int snp_type, int allelecount, string flag, string outputfilename){
    
    int bin_size =1;
    int number_of_bins = (int) ceil(maxval / bin_size);
    vector <long unsigned int> histogram(number_of_bins);
    for (int i = 0; i < number_of_bins; ++i){
        histogram[i] = 0;
    }
    
    boost::unordered::unordered_map<string, vector < variantfeatures > >::iterator it1;
    boost::unordered::unordered_map<string, double>::iterator iter;
    
    
    for (it1 = variantmap.begin(); it1 != variantmap.end(); ++it1) {
        if(it1->second[vcf_index].snp_type == snp_type && it1->second[vcf_index].alt_allele_count == allelecount)
            FillHistogram(histogram, bin_size, number_of_bins, it1->second[vcf_index].flag_values[flag]);
    }
   
    string temp = vcf_name.substr(0,vcf_name.size()-3);
    temp.append(outputfilename);
    ofstream outf(temp.c_str());
    
    double sum = 0;
    for (int i = 0; i < histogram.size(); ++i)
        sum += histogram[i];
    
    for (int i = 0; i < histogram.size(); ++i)
        outf << i*bin_size << '\t' << histogram[i] << '\t' << double (histogram[i]/sum) << '\t' << endl;
    
    histogram.clear();
    
}

void ProcessVCF::GenerateHistogram(int maxval, string vcf_name, int vcf_index, int snp_type, int allelecount, string flag, string outputfilename){
    
    int bin_size = 1;
    int number_of_bins = (int) ceil(maxval / bin_size);
        vector <long unsigned int> histogram(number_of_bins);
    for (int i = 0; i < number_of_bins; ++i){
        histogram[i] = 0;
    }
    
    boost::unordered::unordered_map<string, vector < variantfeatures > >::iterator it1;
    boost::unordered::unordered_map<string, double>::iterator iter;
    
    
    for (it1 = variantmap.begin(); it1 != variantmap.end(); ++it1) {
        if(it1->second[vcf_index].found){
            if (it1->second[vcf_index].alt_allele_count == allelecount) {
                if(flag == "AD_ref")
                    FillHistogram(histogram, bin_size, number_of_bins, it1->second[vcf_index].ref_read_depth);
                if(flag == "AD_alt")
                    FillHistogram(histogram, bin_size, number_of_bins, it1->second[vcf_index].alt_read_depth);
                if(flag == "qual")
                    FillHistogram(histogram, bin_size, number_of_bins, it1->second[vcf_index].qual);
            }
        }
    }
    string temp = vcf_name.substr(0,vcf_name.size()-3);
    temp.append(outputfilename);
    ofstream outf(temp.c_str());
    
    double sum = 0;
    for (int i = 0; i < histogram.size(); ++i)
        sum += histogram[i];
    
    for (int i = 0; i < histogram.size(); ++i)
        outf << i*bin_size << '\t' << histogram[i] << '\t' << double (histogram[i]/sum) << '\t' << endl;
    
    histogram.clear();
    
}
void ProcessVCF::FillHistogram(vector<long unsigned int> & histogram, int bin_size, int number_of_bins, double value){
    
    int bucket = (int)floor(value / bin_size);
    if(bucket > number_of_bins){
        histogram.resize((bucket+1),0);
    }
    histogram[bucket] += 1;
    
}

void ProcessVCF::TiTvRatio(int vcf_index){
    
    vcfstats[vcf_index].transitions = 0;
    vcfstats[vcf_index].transversions = 0;
    
    boost::unordered::unordered_map<string, vector < variantfeatures > >::iterator it;
    for (it = variantmap.begin(); it != variantmap.end(); ++it) {
        if (it->second[vcf_index].found) {
            if ((it->second[vcf_index].refgenotype == "A" || it->second[vcf_index].refgenotype == "G") && (it->second[vcf_index].altgenotype == "A" || it->second[vcf_index].altgenotype == "G"))
                    ++vcfstats[vcf_index].transitions;
            else{
                if((it->second[vcf_index].refgenotype == "C" || it->second[vcf_index].refgenotype == "T") && (it->second[vcf_index].altgenotype == "C" || it->second[vcf_index].altgenotype == "T"))
                    ++vcfstats[vcf_index].transitions;
                else
                    ++vcfstats[vcf_index].transversions;
            }
        }
    }
}


#endif
