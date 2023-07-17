//
// Created by francois on 29.04.23.
//
#include <iostream>
#include <fstream>

using namespace std;

int main() {

    for(int size = 1; size <= 100;size++) {
        string outputFileName = "cube_" + to_string(size);

        ofstream outputFile;
        outputFile.open(outputFileName, ios::out | ios::trunc );

        outputFile << 2*size << " " << size << endl;

        for(int i = 0; i < size;i++) {
            int j = 0;
            outputFile << "1 ";
            for(;j < i;j++) outputFile << "0 ";
            outputFile << "1 ";
            j++;
            for(; j < size; j++) outputFile << "0 ";
            outputFile << "\n";
        }

        for(int i = size-1; i >= 0;i--) {
            int j = 0;
            outputFile << "1 ";
            for(;j < i;j++) outputFile << "0 ";
            outputFile << "-1 ";
            j++;
            for(; j < size; j++) outputFile << "0 ";
            outputFile << "\n";
        }

        outputFile.close();
    }

    return 0;
}