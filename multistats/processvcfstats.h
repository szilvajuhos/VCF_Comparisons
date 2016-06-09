//
//  processvcfstats.h
//  VCF_Comparisons
//
//  Created by Pelin Sahlen on 17/09/2015.
//  Copyright (c) 2015 Pelin Sahlen. All rights reserved.
//
struct stats{
    double x[8][2];
    double z[6];
    
    double discordantgens;
    double discordantdbsnpgens;
    double nofvalsnps;
    double valsnpsfoundinboth;
    double valsnpsfoundinfirst;
    double valsnpsfoundinsecond;
};
vector < stats > statsvector;

void readstatfile(string filename){
    
    statsvector.push_back(stats());
    ifstream infile(filename.c_str());
    string temp;
    getline(infile,temp);
    getline(infile,temp);
    infile >> statsvector.back().x[0][0] >> statsvector.back().x[1][0] >> statsvector.back().x[2][0];
    getline(infile,temp);
    getline(infile,temp);
    getline(infile,temp);
    getline(infile,temp);
    getline(infile,temp);
    infile >> temp >> temp >> statsvector.back().x[3][0] >> temp >> statsvector.back().x[4][0];
    getline(infile,temp);
    getline(infile,temp);
    getline(infile,temp);
    infile >> statsvector.back().x[0][1] >> statsvector.back().x[1][1] >> statsvector.back().x[2][1];
    getline(infile,temp);
    getline(infile,temp);
    getline(infile,temp);
    getline(infile,temp);
    getline(infile,temp);
    infile >> temp >> temp >> statsvector.back().x[3][1] >> temp >> statsvector.back().x[4][1];
/*for indels exclude this   
    infile >> temp >> temp >> temp >> temp >> temp;

    infile >> statsvector.back().z[0] >> temp >> temp >> temp >> temp >> statsvector.back().z[1]
           >> temp >> temp >> temp >> temp >> statsvector.back().z[2]
          >> temp >> temp >> temp >> temp >> statsvector.back().z[3];
    infile >> temp >> temp >> temp >> temp >> temp;
*/ //for indels finished exclusion

    infile >> temp >> temp >> temp;
    infile >> statsvector.back().z[4] >> temp >> temp >> temp >> statsvector.back().z[5];
    infile >> temp >> temp >> temp >> temp;
    infile >> temp >> temp >> temp >> temp >> temp >> statsvector.back().x[5][0] >> statsvector.back().x[5][1];
    infile >> temp >> temp >> temp >> temp >> temp >> temp >> statsvector.back().x[6][0] >> statsvector.back().x[6][1];
    infile >> temp >> temp >> temp >> temp >> temp >> temp >> statsvector.back().x[7][0] >> statsvector.back().x[7][1];
    cout << filename << " stats read" << endl;
    infile.close();
    cout.setf(ios_base::fixed);
    cout.precision(2);
    for (int i = 0; i < 8; ++i) {
        cout << statsvector.back().x[i][0] << "  " << statsvector.back().x[i][1] << endl;
    }
    for (int i = 0; i < 6; ++i) {
        cout << statsvector.back().z[i] << endl;
    }
    
}
void readvcfstats(string infilename, string outfilename){
    
    vector< string > filenames;
    string fname1, fname2,buffer;
    
    ofstream outfile(outfilename.c_str());
    
    outfile.setf(ios_base::fixed);
    outfile.precision(6);
    
    ifstream filelist(infilename.c_str());

    filelist >> fname1 >> fname2;
    int a = fname1.find("/");
    int b = fname2.find("/");
    string comparison1 = fname1.substr(0,a);
    string comparison2 = fname2.substr(0,b);
    filenames.push_back(fname1);
    
    alglib::real_1d_array d1;
    alglib::real_1d_array d2;
    alglib::real_1d_array d3;
    std::vector<double> y1, y2, y3;

    bool samefile = true;

do{
    int i = 0;
    while (samefile){
        ++i;
        filenames.push_back(fname2);
        fname1 = fname2;
        a = fname1.find("/");
        comparison1 = fname1.substr(0,a);
        filelist >> fname2;
        b = fname2.find("/");
        comparison2 = fname2.substr(0,b);
        if(comparison1.compare(comparison2))
            samefile = false;
    }
    
    for (int i = 0; i < filenames.size(); ++i) {
        readstatfile(filenames[i]);
    }
    
    outfile << comparison1 << endl;

    double ms1[8], ms2[8], ms3[6];
    double ss1[8], ss2[8], ss3[6];

    
    
    cout << filenames.size() << endl;
    
    for (int i = 0; i < 8; ++i) {
        y1.resize(0); y2.resize(0);
        for (int j = 0; j < filenames.size(); ++j) {
            y1.push_back(statsvector[j].x[i][0]);
            y2.push_back(statsvector[j].x[i][1]);
        }
        d1.setcontent(y1.size(), &(y1[0]));
        d2.setcontent(y2.size(), &(y2[0]));
        
        ms1[i] = alglib::samplemean(d1);
        alglib::sampleadev(d1, ss1[i]);
        
        ms2[i] = alglib::samplemean(d2);
        alglib::sampleadev(d2, ss2[i]);
        
        if( i < 3)
            outfile.precision(0);
        else
            outfile.precision(6);
        outfile << i << '\t' << ms1[i] << " ± " << ss1[i] << '\t' << ms2[i] << " ± " << ss2[i] << endl;
    }
    for (int i = 0; i < 6; ++i) {
        y3.resize(0);
        for (int j = 0; j < filenames.size(); ++j) {
            y3.push_back(statsvector[j].z[i]);
        }
        d3.setcontent(y3.size(), &(y3[0]));
        
        ms3[i] = alglib::samplemean(d3);
        alglib::sampleadev(d3, ss3[i]);
        
        outfile << i << '\t' << ms3[i] << " ± " << ss3[i] << endl;
    }
    
    statsvector.resize(0);
    filenames.resize(0);
    filenames.push_back(fname2);
    comparison1 = comparison2;
    filelist >> fname2;
    int b = fname2.find("/");
    string comparison2 = fname2.substr(0,b);
    if(comparison1.compare(comparison2))
        samefile = false;
    else
        samefile = true;
    
}while(fname2 != "END");
    outfile << endl;
}
