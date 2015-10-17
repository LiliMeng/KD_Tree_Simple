/**
 * File:  kNNSearch.cpp
 * Author Lili Meng (lilimeng1103@gmail.com)
 * The required application 2: Read a simple text fle and reports the top 3 nearest neighbors using both exact search and BBF approximate search
 */

#include <iostream>
#include "KD_tree.h"
#include "ReadData1.h"

using namespace std;


int main(int argc, const char * argv[])
{

    int K = 3;

    vector<vector<double> > dataset;
    ReadData rd1("sample_data.txt");
    int N = rd1.get_num_of_elements();
    int dim = rd1.get_num_of_dimensions();
    dataset=rd1.allDataPointsVec;

    //query_point
    vector<double> query_point;
    vector<vector<double> > query_point_dataset;
    ReadData rd2("query_points.txt");
    int N2 = rd2.get_num_of_elements();
    int dim2 = rd2.get_num_of_dimensions();
    query_point_dataset=rd2.allDataPointsVec;

    KD_tree_node* root;
    KD_tree tree;
    int max_leaf_size =1;
    tree.create_tree(dataset,max_leaf_size);


    vector<int> indices;
    vector<double> squared_distances;

    /** Exact kNN Search **/
    cout<<"-------------------------------------------------------------------------------"<<endl;
    cout<<"*********Using Exact k Nearest Neighbor Search, The following are Results******"<<endl;
    for(int i=0; i<query_point_dataset.size(); i++)
    {
        tree.kNN_query(query_point_dataset[i], K, indices, squared_distances);

        for (int j = 0; j<K; j++)
        {
            cout<<"For the number row  "<<i<<"  query point, Using Exact kNN Search 3 Nearest Neigbour : The number "<<j+1<<" nearest neighbor index is  "<<indices[j]<<endl;
        }

    }


    cout<<"--------------------------------------------------------------------------------"<<endl;
    cout<<"--------------------------------------------------------------------------------"<<endl;
    cout<<"--------------------------------------------------------------------------------"<<endl;
    cout<<"***Using BBF Approximate k Nearest Neighbors Search, the following are Results**"<<endl;

    for(int i=0; i<query_point_dataset.size(); i++)
    {
        int max_epoch=1000;
        tree.bbf_kNN_query(query_point_dataset[i], K, indices, squared_distances, max_epoch);

        for (int j = 0; j<K; j++)
        {
            cout<<"For the number row  "<<i<<"  query point, Using Exact kNN Search 3 Nearest Neigbour : The number "<<j+1<<" nearest neighbor index is  "<<indices[j]<<endl;
        }

    }


    return 0;
}
