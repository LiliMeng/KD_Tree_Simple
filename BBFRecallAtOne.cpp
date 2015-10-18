/**
 * File:  BBFRecallAtOne.cpp
 * Author Lili Meng (lilimeng1103@gmail.com)
 * The Recall of BBF Approximate Search at 1 nearest neighbor for 10 query points from a dataset of 1000 points
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
    int max_leaf_size =4;
    tree.create_tree(dataset,max_leaf_size);

    vector<int> indices1;
    vector<double> squared_distances1;
    vector<int> trueSearch(10);
    for(int i=0; i<query_point_dataset.size(); i++)
    {
        tree.kNN_query(query_point_dataset[i], 1, indices1, squared_distances1);
        //cout<<indices1[0]<<endl;
        trueSearch[i]=indices1[0];
       // cout<<trueSearch[i]<<endl;

    }

    vector<int> indices2;
    vector<double> squared_distances2;

    vector<double> correct_number(10);



    for(int max_searched_leaf_number=1; max_searched_leaf_number<=1000; max_searched_leaf_number++)
    //for(int max_searched_leaf_number=10; max_searched_leaf_number<1000; max_searched_leaf_number++)
    /** If the max_epoch is equal or larger than the size of all leaves, the result is the same with the exact kNN search**/
   {
        for(int j=0; j<query_point_dataset.size(); j++)
        {
            tree.bbf_kNN_query(query_point_dataset[j], 1, indices2, squared_distances2, max_searched_leaf_number);
            if(indices2[0]==trueSearch[j])
            {
                correct_number[j]++;
            }
        }
   }

   double sum=0.0;
   for(int i=0; i<query_point_dataset.size(); i++)
   {
        sum+=correct_number[i];
   }
    double RecallAtOne=sum/10;
    cout<<"The Recall of BBF Approximate Search at 1 nearest neighbor for 10 query points from a dataset of 1000 points is: "<<setprecision(3)<<RecallAtOne/1000<<endl;



 return 0;
}
