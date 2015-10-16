/**
 * File: ReadData1.cpp
 * Author: Lili Meng (lilimeng1103@gmail.com)
 * read the line and column data from .txt file
 */

#include "ReadData.h"


int main() {

    ReadData l("sample_data.txt");

    cout<<l.get_num_of_dimensions()<<endl;
    cout<<l.get_num_of_elements()<<endl;

    l.printDataID();


    return 0;
}
