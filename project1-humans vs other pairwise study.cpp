#include<iostream>
#include<fstream>
#include <bits/stdc++.h>
using namespace std;

 int TotalMatchCount = 0;
 int TotalMismatchCount = 0;
 int TotalTransition= 0;
 int TotalTransversion = 0;
 unordered_map<string, int> alignment;

void getAlignments(string a, string b)
{
    int matchCount = 0;
    int mismatchCount = 0;
    int transition =0;
    int transversion =0;
    int len1 = a.length(); 
    int len2 = b.length();
    int maxlen=0;

    if(len2 < len1){
        maxlen = len1;
    }
    else
        maxlen = len2;
    for(int i =0;i< maxlen;i++){
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
        if ((c1 == '-') || (c2 == '-') || (c1 != 'A' && c1 != 'G' && c1 != 'T' && c1 != 'C') || (c2 != 'A' && c2 != 'G' && c2 != 'T' && c2 != 'C')){
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
    }
    
    TotalMatchCount += matchCount;
    TotalMismatchCount += mismatchCount;
    TotalTransition += transition;
    TotalTransversion += transversion;
}

void setZero(){
    TotalMatchCount = 0;
    TotalMismatchCount = 0;
    TotalTransition= 0;
    TotalTransversion = 0; 
}

int main(){
    setZero();
    string chromosomes[2];
    int chromosomesCount = 0;
    string line;
    float substitution_rate =0;
    float transRatio = 0;\
    //uncomment to run for each file."C:\Users\hibah\Downloads\bioinfor\test.txt"
    string filepaths[1] ={"C:/Users/hibah/Downloads/bioinfor/test.txt"};
    // "C:/Users/hibah/Downloads/bioinfor/human.chr22.mouse.maf","C:/Users/hibah/Downloads/bioinfor/human.chr22.dog.maf"};
    for(int i = 0;i<1;i++)
    {
        cout<<"Reading File-"<<i+1<<endl;
        ifstream files(filepaths[i]);
        while (getline(files, line)){
            char c = line[0];
            if (c == 's')
            {
            size_t found = line.find_last_of(" ");
            chromosomes[chromosomesCount++] = line.substr(found+1);
            if (chromosomesCount == 2)
            {
            getAlignments(chromosomes[0], chromosomes[1]);
            chromosomesCount = 0;
            substitution_rate = float(TotalMismatchCount) /(TotalMatchCount+TotalMismatchCount);
            }
            }
        }
        // cout<<"Substitution rate- "<< substitution_rate<<endl;
        //     cout<<"Total transition- " << TotalTransition<<endl;
        //     cout<<"Total Transversion-"<< TotalTransversion<<endl;
        //     cout<<"Mismatches-"<< TotalMismatchCount<<endl;
        //     cout<<"Matches-"<<TotalMatchCount<<endl;
        //     transRatio = float(TotalTransition)/ float(TotalTransversion);
        //     cout<<"TransRatio-"<<transRatio<<endl;
        //setZero();
        files.close();
        }

        cout<<"MAP VALUES:"<<endl;
        for(auto num:alignment){
            cout<< num.first <<","<<num.second<<endl;
        }
        substitution_rate = float(TotalMismatchCount) /(TotalMatchCount+TotalMismatchCount);
        cout<<"Substitution rate- "<< substitution_rate<<endl;
        cout<<"Total transition- " << TotalTransition<<endl;
        cout<<"Total Transversion-"<< TotalTransversion<<endl;
        cout<<"Mismatches-"<< TotalMismatchCount<<endl;
        cout<<"Matches-"<<TotalMatchCount<<endl;
        transRatio = float(TotalTransition)/ float(TotalTransversion);
        cout<<"TransRatio-"<<transRatio<<endl;
       
    return 0;
}




