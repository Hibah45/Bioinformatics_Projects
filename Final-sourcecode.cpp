#include<iostream>
#include<fstream>
#include<tuple>
#include <bits/stdc++.h>
using namespace std;

 int TotalMatchCount = 0;
 int TotalMismatchCount = 0;
 int TotalTransition= 0;
 int TotalTransversion = 0;
 int TotalGaps=0;
 int LenGapCount[1000000];
 int max_gaps = 0;
 unordered_map<string, int> alignment;
int matchSum = 0, mismatchSum=0, gapSum= 0;
int transversionSum = 0, transitionSum=0;


int arr[6][2] = {{265, 13467},
{13467, 21554},{21562,25383},
{26244,26471},
{26522,27190},
{28273,29532}};

int start_index=0, end_index = 0;


void set_gap_length(int currenlen)
{
    LenGapCount[currenlen]++;
    if (currenlen > max_gaps)
        max_gaps = currenlen;
}


void getAlignments(string a, string b)
{
    int matchCount = 0;
    int gaplength = 0;
    int current_gap_len = 0;
    bool continuous_gap = false;
    int mismatchCount = 0;
    int transition =0;
    int transversion =0;
    int len1 = a.length(); 
    int len2 = b.length();
    int maxlen=0;

    if(len2 < len1){
        len1 = len2;
    }
   
    for(int i =0;i< len1;i++){
        char c1 = toupper(a[i]);
        char c2 = toupper(b[i]);
        string str1;
        str1.append(1,c1);
        string str2;
        str2.append(1,c2);

        string key;
        key.append(str1);
        key.append(str2);
        //cout<<key<<endl;
        if((c1 == '-') || (c2 == '-')){
            if (continuous_gap)
            {
                current_gap_len++;
            }
            else
            {
            continuous_gap = true;
            current_gap_len = 1; 
            gaplength ++;
            }
        }
        else
        {   
            if (continuous_gap)
            {
                set_gap_length(current_gap_len);
            }
            if ((c1 != 'A' && c1 != 'G' && c1 != 'T' && c1 != 'C') || (c2 != 'A' && c2 != 'G' && c2 != 'T' && c2 != 'C' && c1!= '-' && c2!='-')){
                continue;
            }
            else if(c1 == c2){
                matchCount++;

            }
            else if((c1== 'A' && c2=='G') || (c1 =='G' && c2 =='A') || (c1 == 'T' && c2 =='C') || (c1=='C' && c2 =='T')){
                transition++;
                mismatchCount++;
            }
            else if((c1== 'A' && c2=='T') || (c1 =='T' && c2 =='A') || (c1 == 'A' && c2 =='C') || (c1=='C' && c2 =='A') || (c1== 'T' && c2=='G') || (c1 =='G' && c2 =='T') || (c1 == 'G' && c2 =='C') || (c1=='C' && c2 =='G')){
                transversion++;
                mismatchCount++;
            }
            else
                mismatchCount++;
                if(alignment.find(key) != alignment.end()){
                    alignment[key] += 1;
                }
                else{
                    alignment[key]=1;
                }
                continuous_gap = false; 
        }
    }
    
    if (continuous_gap)
    {
        set_gap_length(current_gap_len);
    }
    TotalMatchCount += matchCount;
    TotalMismatchCount += mismatchCount;
    TotalTransition += transition;
    TotalTransversion += transversion;
    TotalGaps += gaplength;

    // transition=0;
    // transversion=0;
}

void SetToZero(){
    TotalMatchCount = 0;
    TotalMismatchCount = 0;
    TotalGaps = 0;
    TotalTransition=0;
    TotalTransversion=0;
    for (int j=0; j < 1000000; j++)
    {
    LenGapCount[j] = 0;
    }
    max_gaps = 0;
}

int main(){
    string chr[2];
    int chrCount = 0;
    string line;
    float substitution_rate =0;
    float gap_rate =0;
    int totalGapCount = 0;
    float transRatio = 0;
    string filepaths[1] = {"C:/Users/hibah/projects/termproject/sars2.mers.sing.maf"};
    int flag_sars=0;
    int flag_mers=0;

    for(int i = 0;i<(sizeof(arr)/sizeof(arr[0]));i++){
        // cout<<arr[i][0]<<","<<arr[i][1];
         cout<<"Reading File-"<<i<<"\n\n";
        ifstream files(filepaths[0]);

        while (getline(files, line)){
            char c = line[0];
           if (c == 's')
            {
            std::string delimiter = " ";

            string s = line;
            size_t pos = 0;
            std::string token;
            int token_finder = 0;
            string start_seq, seq_len, genomes;
            while ((pos = s.find(delimiter)) != std::string::npos) {
                token = s.substr(0, pos);
                if(token_finder == 3 && token != " ")
                {
                    seq_len = token;
                    // std::cout << token << std::endl;
                }
                if(token_finder == 2 && token != " ")
                {
                    start_seq = token;
                    // std::cout << token << std::endl;
                }
                if(token_finder == 1 && token != " ")
                {
                    genomes = token;
                    if(genomes == "sars2")
                        flag_sars = 1;
                    if(genomes == "mers")
                            flag_mers = 1;
                    // std::cout << token << std::endl;
                }
                s.erase(0, pos + delimiter.length());
                token_finder++;
            }
            
            token_finder = 0;

            int token_start = stoi(start_seq);
            int token_end = token_start + stoi(seq_len);


            if((token_end>=arr[i][0] && token_start<=arr[i][1]) || genomes=="mers"){
                size_t found = line.find_last_of(" ");
                // chr[chrCount++] = line.substr(found+1);
                if(genomes=="sars2" || (token_end>=arr[i][0] && token_start<=arr[i][1])){
                    start_index=0;
                    end_index=token_end;

                if(arr[i][0]-token_start>0)
                    start_index=arr[i][0]-token_start;
                if(arr[i][1]-token_end>0 && token_end>arr[i][1])
                    end_index = arr[i][1]-token_end;
                if(arr[i][1]-token_end<0)
                    end_index = arr[i][1];
                
                end_index-=token_start;
                }
                cout<< "RANGES::"<<arr[i][0]<<","<<arr[i][1]<<"\n token start:"<<token_start<<"\t"<<"token end:"<<token_end<<"\n";
                cout<<"start index:"<<start_index<<" "<<"end index:"<<end_index<<" "<<"genome:"<<genomes<<"\n";
                if(end_index!=0){
                    if(start_index!=375){
                        chr[chrCount++] = line.substr(found+1).substr(start_index, end_index-start_index);
                    }
                }
                if(genomes=="mers")
                {
                    start_index=0;
                    end_index=0;
                }
                // cout<<"first chromo"<<chr[0]<<endl;
                if (chrCount == 2)
                {
                    getAlignments(chr[0], chr[1]);
                    cout<<""<<endl;
                    cout<<"transition-" << TotalTransition<<endl;
                    cout<<"Transversion-"<< TotalTransversion<<endl;
                    cout<<"Mismatches-"<< TotalMismatchCount<<endl;
                    cout<<"Matches-"<<TotalMatchCount<<endl;
                    cout<<"Gap sum-"<<TotalGaps<<endl;

                    matchSum += TotalMatchCount;
                    mismatchSum += TotalMismatchCount;
                    gapSum += TotalGaps;
                    transitionSum +=TotalTransition;
                    transversionSum +=TotalTransversion;
                    // cout<<"transition sum-"<<transitionSum<<endl;
                    chrCount = 0;
                }
                flag_mers=0;
                flag_sars=0;
                
                cout << "Gap Length \t Gap Count \tGap Frequency" << endl;
                for (int j = 1; j <= max_gaps; j++)
                {
                    totalGapCount += LenGapCount[j];
                    cout << "\t\t" << j << "\t\t" << LenGapCount[j] << "\t\t" << float(LenGapCount[j]) / TotalGaps << endl;
                }
    
                // cout << "\t\tTotal" << "\t\t" << totalGapCount << "\t\t" << float(totalGapCount) / TotalGaps << endl;
                cout << endl;
                
                SetToZero();
            }

            }
        }
       
        cout<<"MAP VALUES OF ALIGNMENT:"<<endl;
        for(auto num:alignment){
            cout<< num.first <<","<<num.second<<endl;
        }
        cout<<"********************FINAL RESULT***********"<<endl;
        cout<<"\n"<<"Mismatches of alignment-"<< mismatchSum<<endl;
        cout<<"Matches of alignment-"<<matchSum<<endl;
        gap_rate = float(gapSum) / (matchSum + mismatchSum + gapSum);
        substitution_rate = float(mismatchSum) /(matchSum+mismatchSum);
        transRatio = float(transitionSum)/ float(transversionSum);

        cout<<"Substitution rate of alignment- "<< substitution_rate<<endl;
        cout<<"gap sum of the alignment-"<<gapSum<<endl;
        cout<<"Gaps rate of alignment- "<< gap_rate<<endl;
        cout<<"Total transition of alignment- " << transitionSum<<endl;
        cout<<"Total Transversion of alignment-"<< transversionSum<<endl;
        cout<<"TransRatio of alignment-"<<transRatio<<endl;
    
        matchSum = 0;
        mismatchSum=0;
        gapSum= 0;
        transitionSum=0;
        transversionSum =0;
        

    }

       
}
