
#include<iostream>
#include<fstream>
#include <bits/stdc++.h>
using namespace std;
int TotalMatchCount = 0;
int TotalMismatchCount = 0;
int TotalGaps = 0;
int LenGapCount[1000000];
int max_gaps = 0;

void set_gap_length(int currenlen)
{
    LenGapCount[currenlen]++;
    if (currenlen > max_gaps)
        max_gaps = currenlen;
}

void getAlignments(string a, string b)
{
    int matchCount = 0;
    int mismatchCount = 0;
    int gaplength = 0;
    int current_gap_len = 0;
    bool continuous_gap = false;
    int len1 = a.length(); 
    int len2 = b.length();
    
    if (len2 < len1) 
    {
    len1 = len2;
    }

    for (int i =0; i < len1; i++)
    {
        char c1 = toupper(a[i]);
        char c2 = toupper(b[i]);
        if (c1 == '-' || c2 == '-') 
        {
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
            if ((c1 != 'A' && c1 != 'G' && c1 != 'T' && c1 != 'C') || (c2 != 'A' && c2 != 'G' && c2 != 'T' && c2 != 'C' && c1!= '-' && c2!='-'))
            {
                continue;
            }
            else if(c1 == c2)
            {
                matchCount++;
            }
            else
            {
                mismatchCount++;
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
    TotalGaps += gaplength;
}




void SetToZero()
    {
    TotalMatchCount = 0;
    TotalMismatchCount = 0;
    TotalGaps = 0;
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
   //string filepaths[3] ={"C:/Users/hibah/Downloads/bioinfor/human.chr22.chimp.maf","C:/Users/hibah/Downloads/bioinfor/human.chr22.mouse.maf","C:/Users/hibah/Downloads/bioinfor/human.chr22.dog.maf"};
    string filepaths[3] = {"C:/Users/hibah/projects/termproject/sars2.mers.sing.maf", "",""};
    for ( int i=0;i<1;i++)
    {
        cout<<"*********READING FILE:"<<i<<"\n";
        int totalGapCount = 0;
        ifstream files(filepaths[i]);
        while (getline(files, line))
        {
            char c = line[0];
            if (c == 's')
            {
            size_t found = line.find_last_of(" ");
            chr[chrCount++] = line.substr(found+1);
            if (chrCount == 2)
            {
            getAlignments(chr[0], chr[1]);
            chrCount = 0;
            }
        }
        }
    float gapRate = float(TotalGaps) / (TotalMatchCount +TotalMismatchCount + TotalGaps);
    cout << "Gap Rate:" << gapRate << endl;
    cout << "Gap Length \t Gap Count \tGap Frequency" << endl;
    for (int j = 1; j <= max_gaps; j++)
    {
        totalGapCount += LenGapCount[j];
        cout << "\t\t" << j << "\t\t" << LenGapCount[j] << "\t\t" << float(LenGapCount[j]) / TotalGaps << endl;
    }
    cout << "\t\tTotal" << "\t\t" << totalGapCount << "\t\t" << float(totalGapCount) / TotalGaps << endl;
    cout << endl;
    SetToZero();
    files.close();
    }
    return 0;
}