/**
 * File: PrintDataID.cpp
 * Author: Lili Meng (lilimeng1103@gmail.com)
 * read the line and column data from .txt file
 */

#include "ReadData1.h"


int main() {

    string fileName;
    fileName = "sample_data.txt";
    ReadData l(fileName);


    cout<<"The Data ID from "<<fileName<<" "<<"is the following:"<<endl;
    l.printDataID();


    return 0;
}
