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
    query_point=query_point_dataset[1];

    KD_tree_node* root;
    KD_tree tree;
    tree.create_tree(dataset, 1, 128);
   

    vector<int> indices;
    vector<double> squared_distances;
    tree.kNN_query(query_point, K, indices, squared_distances);

    for (int i = 0; i<indices.size(); i++) {
        cout<<"index "<<indices[i]<<"\t"<<"distance is "<<squared_distances[i]<<endl;
    }



    // brute force
    vector<double> brute_force_distances;
    for (int i = 0; i< N; i++)
    {
        double dist = KD_tree::distance_sq(query_point, dataset[i]);
        brute_force_distances.push_back(dist);
    }

    std::sort(brute_force_distances.begin(), brute_force_distances.end());
    for(int i = 0; i<K; i++)
    {
        cout<<"brute-force distance is "<<brute_force_distances[i]<<endl;
    }





    return 0;
}
