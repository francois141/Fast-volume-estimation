//
// Created by francois on 29.04.23.
//
#include <iostream>
#include <fstream>

using namespace std;

int main() {

    for(int size = 1; size <= 100;size++) {
        string outputFileName = "simplex_" + to_string(size);

        ofstream outputFile;
        outputFile.open(outputFileName, ios::out | ios::trunc );

        outputFile << size+1 << " " << size << endl;

        for(int i = size; i >= 1;i--) {
            int j = 0;
            for(;j < i;j++) outputFile << "0 ";
            outputFile << "-1 ";
            j++;
            for(; j <= size; j++) outputFile << "0 ";
            outputFile << "\n";
        }

        for(int i = 0; i < size+1;i++) {
            outputFile << "1 ";
        }
        outputFile << "\n";

        outputFile.close();
    }

    return 0;
}