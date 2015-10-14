#include <iostream>
#include "KD_tree.h"
#include "ReadData.h"

using namespace std;

int main(int argc, const char * argv[])
{
    /*
    // test on random data

    int N = 100;
    int K = 3;
    int dim = 8;

    vector<vector<double> > dataset;
    for (int i = 0; i<N; i++) {
        vector<double> data;
        for (int j = 0; j<dim; j++) {
            double v = rand()%256/256.0;
            data.push_back(v);
        }
        dataset.push_back(data);
    }

    vector<double> query_point;
    for (int j = 0; j<dim; j++) {
        double v = rand()%256/256.0;
        query_point.push_back(v);
    }
    */
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

    KD_tree tree;
    tree.create_tree(dataset, 100, 128);


    vector<int> indices;
    vector<double> squared_distances;
    tree.query(query_point, K, indices, squared_distances);

    for (int i = 0; i<indices.size(); i++) {
        printf("index %d  distance is %f\n", indices[i], squared_distances[i]);
    }



    // brute force
    vector<double> brute_force_distances;
    for (int i = 0; i<N; i++) {
        double dis = KD_tree::distance_sq(query_point, dataset[i]);
        brute_force_distances.push_back(dis);
    }

    std::sort(brute_force_distances.begin(), brute_force_distances.end());
    for (int i = 0; i<K; i++) {
        printf("brute force distance is %f\n", brute_force_distances[i]);
    }




    return 0;
}
