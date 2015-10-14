#include <iostream>
#include "KD_tree.h"

using namespace std;

int main(int argc, const char * argv[])
{
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

    KD_tree tree;
    tree.create_tree(dataset, 4, 4);

    vector<int> indices;
    vector<double> squared_distances;
    tree.query(query_point, K, indices, squared_distances);

    for (int i = 0; i<indices.size(); i++) {
        printf("index, distance are %d %f\n", indices[i], squared_distances[i]);
    }

    // brute force
    vector<double> brute_force_distances;
    for (int i = 0; i<N; i++) {
        double dis = KD_tree::distance_sq(query_point, dataset[i]);
        brute_force_distances.push_back(dis);
    }

    std::sort(brute_force_distances.begin(), brute_force_distances.end());
    for (int i = 0; i<K; i++) {
        printf("brute force distance are %f\n", brute_force_distances[i]);
    }




    return 0;
}

