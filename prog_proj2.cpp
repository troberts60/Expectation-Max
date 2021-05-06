#include<iostream>
#include<vector>
#include<sstream>
#include<fstream>
#include<string>
#include<stdlib.h>
#include<math.h>
#include<iomanip>

using namespace std;

int main(){
	
	ifstream seq1File;
	string filename;
	string line1;
	string seq1;
	char decision;
	int seq1_size;
	vector<string> totSeq;
	char currChar;
	int r;
	int randAlignNum;
	int emIteration;
	int moWid;
	
	cout << "Input text file: ";
	cin >> filename;
	
	seq1File.open(filename.c_str());
	if (!seq1File.is_open()){
		cout << "Not correct file name." << endl;
	}
	else{
		while(getline(seq1File,line1)){  
			if (line1.at(0) >= 65 && line1.at(0) <= 90) 
				totSeq.push_back(line1);					 
		}
	}
	seq1File.close();
	
	cout << "\nDefault:\nRandom Starting Alignments: 50\nE-M steps: 500\nMotif length: 6" << endl;
	cout << "Do you wish to change settings[y/n]: ";
	cin >> decision;
	cout << endl;
	if(decision == 'n'){
		randAlignNum = 50;
		emIteration = 500;
		moWid = 6;
	}
	else if(decision == 'y'){
		cout << "How many Random Starting Alignments: ";
		cin >> randAlignNum;
		cout << "How many E-M steps: ";
		cin >> emIteration;
		cout << "Motif length: ";
		cin >> moWid;
	}
	else
		"Incorrect response, run program again";
	
	string motStr;
	int motPos[totSeq.size()];
	int countTab[4][moWid+1] = {0};
	float freqTab[4][moWid+1] = {0};
	float oddTab[4][moWid+1] = {0};
	float logTab[4][moWid+1] = {0};
	float maxScore;
	int maxPos;
	float logScore;
	int newMotifLoc[totSeq.size()];
	float entTab[4][moWid];
	float infoTotal = 0;
	int infoMotif[totSeq.size()];
	float maxInfo = 0;
	int totInfo;
	float maxOddScore[totSeq.size()] = {0};
	int maxMotifLoc[totSeq.size()] = {0};
	for (int l = 0; l < randAlignNum; l++){
		for (int i = 0; i < totSeq.size(); i++){
			seq1_size = totSeq[i].length();
			seq1 = totSeq[i];
			for (int j = 0; j < seq1_size; j++){
				currChar = seq1[j];
				if (currChar == 'A')
					countTab[0][0]++;
				if (currChar == 'C')
					countTab[1][0]++;
				if (currChar == 'G')
					countTab[2][0]++;
				if (currChar == 'T')
					countTab[3][0]++; 
			}
		}

		for (int i = 0; i < totSeq.size(); i++){
			seq1_size = totSeq[i].length();
			r = rand() % (seq1_size - moWid + 1);
			motPos[i] = r;
			motStr =  totSeq[i].substr(r, (moWid));
			for (int j = 0; j < moWid+1; j++){
				currChar = motStr[j];
				if (currChar == 'A'){
					countTab[0][j+1]++;
					countTab[0][0]--;
				}
				if (currChar == 'C'){
					countTab[1][j+1]++;
					countTab[1][0]--;
				}
				if (currChar == 'G'){
					countTab[2][j+1]++;
					countTab[2][0]--;
				}
				if (currChar == 'T'){
					countTab[3][j+1]++;
					countTab[3][0]--;
				}
			}
		}
	
		for (int i = 0; i < 4; i++){
			for (int j = 1; j < moWid+1; j++){
				countTab[i][j]++;
			}
		}

		for (int j= 0; j < moWid+1; j++){
			float totCount = 0;
			for (int i = 0; i < 4; i++){
				totCount = totCount + countTab[i][j];	
			}
			for (int k = 0; k < 4; k++){
				freqTab[k][j] = countTab[k][j]/totCount;
			}
		}
		
		for (int i = 0; i < 4; i++){
			for (int j = 0; j < moWid+1; j++){
				oddTab[i][j] = freqTab[i][j]/freqTab[i][0];
			}
		}
		
		for (int i = 0; i < 4; i++){
			for (int j = 0; j < moWid+1; j++){
				logTab[i][j] = log2(oddTab[i][j]);
			}
		}
		
		for (int t = 0; t < emIteration; t++){
			for (int i = 0; i < totSeq.size(); i++){
				seq1_size = totSeq[i].length();
				maxScore = -1000000000;
				maxPos = -1;
				for (int j = 0; j < seq1_size - moWid + 1; j++){
					motStr = totSeq[i].substr(j, (moWid));
					logScore = 0;
					for (int k = 1; k < moWid + 1; k++){
						if (motStr[k] == 'A')
							logScore = logScore + logTab[0][k];
						if (motStr[k] == 'C')
							logScore = logScore + logTab[1][k];
						if (motStr[k] == 'G')
							logScore = logScore + logTab[2][k];
						if (motStr[k] == 'T')
							logScore = logScore + logTab[3][k];	
					}
					if (logScore > maxScore){
						maxScore = logScore;
						maxPos = j;
					}
					if (maxScore > maxOddScore[i]){
						maxOddScore[i] = maxScore;
						maxMotifLoc[i] = j;
					}
				}
				newMotifLoc[i] = maxPos;
			}
	
			for (int i = 0; i < totSeq.size(); i++){
				motStr =  totSeq[i].substr(motPos[i], (moWid));
				for (int j = 0; j < (moWid + 1); j++){
					currChar = motStr[j];
					if (currChar == 'A'){
						countTab[0][j+1]--;
						countTab[0][0]++;
					}
					if (currChar == 'C'){
						countTab[1][j+1]--;
						countTab[1][0]++;
					}
					if (currChar == 'G'){
						countTab[2][j+1]--;
						countTab[2][0]++;
					}
					if (currChar == 'T'){
						countTab[3][j+1]--;
						countTab[3][0]++;
					}
				}
			}
		
			for (int i = 0; i < totSeq.size(); i++){
				motStr =  totSeq[i].substr(newMotifLoc[i], (moWid));
				for (int j = 0; j < moWid + 1; j++){
					currChar = motStr[j];
					if (currChar == 'A'){
						countTab[0][j+1]++;
						countTab[0][0]--;
					}
					if (currChar == 'C'){
						countTab[1][j+1]++;
						countTab[1][0]--;
					}
					if (currChar == 'G'){
						countTab[2][j+1]++;
						countTab[2][0]--;
					}
					if (currChar == 'T'){
						countTab[3][j+1]++;
						countTab[3][0]--;
					}
				}	
			}
		
			for (int j= 0; j < moWid+1; j++){
				float totCount = 0;
				for (int i = 0; i < 4; i++){
					totCount = totCount + countTab[i][j];	
				}
				for (int k = 0; k < 4; k++){
					freqTab[k][j] = countTab[k][j]/totCount;
				}
			}
			
			for (int i = 0; i < 4; i++){
				for (int j = 0; j < moWid+1; j++){
					oddTab[i][j] = freqTab[i][j]/freqTab[i][0];
				}
			}
			
			for (int i = 0; i < 4; i++){
				for (int j = 0; j < 7; j++){
					logTab[i][j] = log2(oddTab[i][j]);
				}
			}
			
			for (int i = 0; i < 4; i++){
				for (int j = 0; j < moWid; j++){
					entTab[i][j] = freqTab[i][j+1] * (log2(freqTab[i][j+1]));
				}	
			}
			for (int j = 0; j < moWid; j++){
				totInfo = 0;
				for (int i = 0; i < 4; i++){
					totInfo = totInfo + entTab[i][j];
				}
				infoTotal = infoTotal + (2 - (-1 * totInfo));
			}
			if (infoTotal > maxInfo){
				maxInfo = infoTotal;
				for (int i = 0; i < totSeq.size(); i++){
					infoMotif[i] = newMotifLoc[i];
				}
			}
			infoTotal = 0;
			for (int i = 0; i < totSeq.size(); i++){
				motPos[i] = newMotifLoc[i];
			}
		}
		for (int i = 0; i < 4; i++){
			for (int j = 0; j < moWid+1; j++){
				countTab[i][j] = 0;
			}
		}
	}
	cout << "Log Score" << "   Motif Loc" << "   Motif" << endl;
	for (int i = 0; i < totSeq.size(); i++){
		cout << "  ";
		cout << printf("%.1f", maxOddScore[i]) << setw(10);
		cout << maxMotifLoc[i] << setw(12);
		seq1 = totSeq[i].substr(maxMotifLoc[i], moWid);
		cout << seq1 << endl;
	}
	
	return 0;
}
