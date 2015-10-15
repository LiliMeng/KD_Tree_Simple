#ifndef KD_tree_H
#define KD_tree_H

#include <vector>
#include <assert.h>
#include <queue>
#include <stack>
#include <set>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <stdio.h>

using namespace std;

/// key value pair. Not use the std::pair because there is default
/// "<" operator overloading for it and it is too heavy for operation.
template <class T, class V = double>
struct KeyValue
{
    T key;
    V value;
    KeyValue(const T k, const V v) : key(k), value(v){}
};

class Bounding_box
{
private:
    vector<double> min_point_;
    vector<double> max_point_;
public:
    Bounding_box()
    {
    }
    Bounding_box(const vector<double> & min_point, const vector<double> & max_point)
    {
        min_point_ = min_point;
        max_point_ = max_point;
    }

    // squared distance from a point to box
    double distance_sq(const vector<double> & point) const
    {
        double sum_sq = 0;
        unsigned long dim = point.size();
        assert( dim == min_point_.size() && dim == max_point_.size() );

        for ( unsigned int i=0; i<dim; ++i )
        {
            double x_min = min_point_[i];
            double x_max = max_point_[i];
            double x = point[i];
            if ( x < x_min ) {
                sum_sq += ( x_min - x ) * (x_min - x);
            }
            else if ( x > x_max )
            {
                sum_sq += ( x_max - x ) * (x_max - x);
            }
            else
            {
                sum_sq =0;
            }


        }
        return sum_sq;
    }
};

class KD_tree_node
{
public:
    int split_dim_;
    double split_value_;
    bool is_leaf_;
    int level_;
    int index_;
    Bounding_box box_;

    KD_tree_node *left_node_;
    KD_tree_node *right_node_;
    vector<int> data_index_;


    KD_tree_node()
    {
        split_dim_ = 0;
        split_value_ = 0.0;
        is_leaf_ = false;
        left_node_ = NULL;
        right_node_ = NULL;
    }

};

class distance_index
{
public:
    int index_;
    double distance_;

    distance_index(int index, double distance)
    {
        index_ = index;
        distance_ = distance;
    }

    bool operator < (const distance_index & other) const
    {
        return this->distance_ < other.distance_;
    }
};

class KD_tree
{
public:
    KD_tree();
    ~KD_tree();

    bool read_tree_from_file(const char* file_name, KD_tree_node* &root);
    bool create_tree(const vector<vector<double> > & data, unsigned max_leaf_size, unsigned max_level);
    bool save_tree_to_file(const char* file_name);
    bool kNN_query(const vector<double> & query_point, const int K,
               vector<int> & indices,
               vector<double> & squared_distances) const;

    bool bbf_kNN_query(const vector<double> & query_point, const int K,
               vector<int> & indices,
               vector<double> & squared_distances, size_t maxEpoch) const;

private:

    typedef KeyValue<KD_tree_node*> node_bind;
    typedef priority_queue<node_bind, vector<node_bind>, greater<node_bind> > node_minPQ;
    typedef KeyValue<vector<double>> point_bind;
    // build and store data
    vector<vector<double> > data_;
    unsigned long dim_;
    unsigned max_leaf_size_;
    unsigned max_level_;
    KD_tree_node *root_;

    // helper function for reading tree from file
    void read_tree_from_file_helper(FILE *pf, KD_tree_node * & node);

    // helper function for creating the tree using point indices
    bool build_tree(const vector<int> & index, KD_tree_node* & node, unsigned level);

    // helper function for save kdtree
    void save_tree_helper(FILE *pf, KD_tree_node* node);

    const KD_tree_node * get_leaf_node(const vector<double> & query_point, const KD_tree_node * node) const;

    // query update from current node
    bool query(const vector<double> & query_point, const int K,
               stack<KD_tree_node *> & nodes,
               set<KD_tree_node *> & visited_nodes,
               priority_queue<distance_index> & priority_points,
               double & max_distance) const;

    // check points in a left node
    bool explore_leaf_node(const vector<double> & query_point, const int K,
                           priority_queue<distance_index> & priority_points,
                           double & max_distance, const KD_tree_node * cur_node) const;


    Bounding_box build_bounding_box(const vector<int> & indices) const;



public:
    // squared of distance fro brute-force double check
    static double distance_sq(const vector<double> & data0, const vector<double> & data1);
};

KD_tree::KD_tree()
{
    root_ = NULL;
}

KD_tree::~KD_tree()
{

}

void KD_tree::read_tree_from_file_helper(FILE *pf, KD_tree_node * & node)
{
    char lineBuf[1024] = {NULL};
    char *ret = fgets(lineBuf, sizeof(lineBuf), pf);
    if (!ret) {
        node = NULL;
        return;
    }
    if (lineBuf[0] == '#') {
        // empty node
        node = NULL;
        return;
    }
    node = new KD_tree_node();
    assert(node);
    int is_leaf = 0;

    sscanf(lineBuf, "%u\t %d\t\t %d\t\t %lf\n", &node->level_, &node->is_leaf_, &node->split_dim_, &node->split_value_);
    node->is_leaf_ = (is_leaf == 1);

    node->left_node_ = NULL;
    node->right_node_ = NULL;

    read_tree_from_file_helper(pf, node->left_node_);
    read_tree_from_file_helper(pf, node->right_node_);

}

bool KD_tree::read_tree_from_file(const char* file_name, KD_tree_node * & root)
{

    FILE *pf = fopen(file_name, "r");
    if (!pf) {
        cout<<"can not open file "<<file_name<<endl;
        return false;
    }
    //read first line
    char lineBuf[1024] = {NULL};
    fgets(lineBuf, sizeof(lineBuf), pf);
    printf("%s\n", lineBuf);
    read_tree_from_file_helper(pf, root);
    root_ = root;
    fclose(pf);
    return true;
}

bool KD_tree::create_tree(const vector<vector<double> >& data, unsigned max_leaf_size, unsigned max_level)
{
    assert(root_ == NULL);

    data_ = data;
    max_leaf_size_ = max_leaf_size;
    max_level_ = max_level;
    dim_ = data.front().size();

    vector<int> indices;
    for(double i = 0; i< data_.size(); i++)
    {
        indices.push_back(i);
    }

    root_ = new KD_tree_node();
    root_->box_=this->build_bounding_box(indices);
    this->build_tree(indices, root_, 0);

    return true;
}

// helper function for creating the tree using point indices
bool KD_tree::build_tree(const vector<int> & indices, KD_tree_node* & node, unsigned level)
{
    assert(node);
    if(level >= max_level_ || indices.size() <=max_leaf_size_)
    {
        node->is_leaf_ = true;
        node->level_ = level;
        node->data_index_ = indices;
        node->box_ = this->build_bounding_box(indices);
        return true;
    }
    node->level_=level;

    // randomly select split dimensions
    int split_dim = rand() % dim_;
    node->split_dim_ = split_dim;
    vector<double> one_dim_values;
    for(int i=0; i<indices.size(); i++)
    {
        one_dim_values.push_back(data_[indices[i]][split_dim]);
    }

    // median value as split value, every value before the meidan is less than the median, every value after the median is larger than the median;
    size_t n = one_dim_values.size() / 2;
    std::nth_element(one_dim_values.begin(), one_dim_values.begin()+n, one_dim_values.end());

    node->split_value_ = one_dim_values[n];

    node->box_ = this->build_bounding_box(indices);

    vector<int> left_indices;
    vector<int> right_indices;
    for(int i=0; i< indices.size(); i++)
    {
        int index = indices[i];
        double v = data_[index][node->split_dim_];
        if(v < node->split_value_)
        {
            left_indices.push_back(index);
        }
        else
        {
            right_indices.push_back(index);
        }

    }

    assert(left_indices.size() + right_indices.size() == indices.size());

    if(left_indices.size() !=0)
    {
        KD_tree_node* left_node = new KD_tree_node();
        this->build_tree(left_indices, left_node, level+1);
        node->left_node_ = left_node;
    }

    if(right_indices.size() != 0)
    {
        KD_tree_node *right_node = new KD_tree_node();
        this->build_tree(right_indices, right_node, level+1);
        node->right_node_ = right_node;
    }

    return true;
}

void KD_tree::save_tree_helper(FILE *pf, KD_tree_node* node)
{
    if (!node) {
        fprintf(pf, "#\n");
        return;
    }
    // write current node
    fprintf(pf, "%u\t %d\t\t %d\t\t %lf\n", node->level_, (int)node->is_leaf_, (int)node->split_dim_, node->split_value_);

    save_tree_helper(pf, node->left_node_);
    save_tree_helper(pf, node->right_node_);
}

bool KD_tree::save_tree_to_file(const char* file_name)
{
    assert(root_);
    FILE *pf = fopen(file_name, "w");
    if (!pf) {
        cout<<"can not open file "<<file_name<<endl;
        return false;
    }

    fprintf(pf, "level_\t is_leaf_\t split_dim_\t split_value_\n");

    save_tree_helper(pf, root_);
    fclose(pf);
    return true;
}

Bounding_box KD_tree::build_bounding_box(const vector<int> & indices) const
{
    assert(indices.size() >= 1);
    vector<double> min_point = data_[indices.front()];
    vector<double> max_point = min_point;

    for(unsigned int i=1; i< indices.size(); ++i)
    {
        vector<double> point = data_[indices[i]];

        for(unsigned int j=0; j<point.size(); ++j)
        {
            if(point[j] < min_point[j])
            {
                min_point[j] = point[j];
            }

            if(point[j] > max_point[j])
            {
                max_point[j] = point[j];
            }
        }

    }

    return Bounding_box(min_point, max_point);

}

const KD_tree_node* KD_tree::get_leaf_node(const vector<double> & query_point, const KD_tree_node* node) const
{
    assert(node);
    if(node->is_leaf_)
    {
        return node;
    }

    int dim = node->split_dim_;
    double split_value = node->split_value_;
    double value = query_point[dim];

    if(value < split_value && node->left_node_)
    {
        return this->get_leaf_node(query_point, node->left_node_);
    }
    else if(node->right_node_)
    {
        return this->get_leaf_node(query_point, node->right_node_);
    }
    else
    {
        return NULL;
    }
}

bool KD_tree::kNN_query(const vector<double> & query_point, const int K,
                    vector<int> & indices,
                    vector<double> & squared_distances) const
{
    assert(root_);
    assert(K < data_.size());

    double max_sq_distance = numeric_limits<double>::max();
    priority_queue<distance_index> distance_queue;

    stack<KD_tree_node*> candidate_nodes;
    KD_tree_node* cur_node = root_;

    // travel to the leaf node
    while(cur_node!=NULL)
    {
        candidate_nodes.push(cur_node);

        if(cur_node->is_leaf_)
        {
            break;
        }
        int dim = cur_node->split_dim_;
        double split_value = cur_node->split_value_;
        double value = query_point[dim];

        if(value <split_value && cur_node->left_node_!=NULL)
        {
            cur_node = cur_node->left_node_;
        }
        else if(cur_node->right_node_!=NULL)
        {
            cur_node = cur_node->right_node_;
        }
    }

    set<KD_tree_node*> visited_nodes;
    this->query(query_point, K, candidate_nodes, visited_nodes, distance_queue, max_sq_distance);

    indices.resize(K);
    squared_distances.resize(K);
    assert(distance_queue.size()==K);

    int num = K -1;
    while(!distance_queue.empty())
    {
        distance_index top = distance_queue.top();
        distance_queue.pop();
        indices[num] = top.index_;
        squared_distances[num] = top.distance_;
        num--;
    }

    return true;

}

bool KD_tree::query(const vector<double>& query_point, const int K,
                    stack<KD_tree_node*>& nodes,
                    set<KD_tree_node*>& visited_nodes,
                    priority_queue<distance_index>& priority_points,
                    double & max_distance) const
{
    while(!nodes.empty())
    {
        KD_tree_node* cur_node = nodes.top();
        nodes.pop();

        //already visited node
        if(visited_nodes.find(cur_node)!= visited_nodes.end())
        {
            continue;
        }

        visited_nodes.insert(cur_node);
        if(cur_node->is_leaf_)
        {
            this->explore_leaf_node(query_point, K, priority_points, max_distance, cur_node);
        }
        else
        {
            // internal node
            if(cur_node->left_node_)
            {
                double dist_sq = cur_node->left_node_->box_.distance_sq(query_point);
                if(dist_sq < max_distance)
                {
                    nodes.push(cur_node->left_node_);
                }
            }
            if(cur_node->right_node_)
            {
                double dist_sq = cur_node->left_node_->box_.distance_sq(query_point);
                if(dist_sq < max_distance)
                {
                    nodes.push(cur_node->right_node_);
                }
            }
        }

    }
    return true;
}

bool KD_tree::explore_leaf_node(const vector<double>& query_point, const int K,
                                priority_queue<distance_index>& priority_points,
                                double & max_distance, const KD_tree_node* cur_node) const
{
    assert(cur_node->is_leaf_);

    for(int i=0; i<cur_node->data_index_.size(); i++)
    {
        int index = cur_node->data_index_[i];
        const vector<double> & point = data_[index];
        double dist_sq = KD_tree::distance_sq(point, query_point);
        if(priority_points.size() < K || dist_sq < max_distance)
        {
            distance_index di(index, dist_sq);
            priority_points.push(di);
            if(priority_points.size() > K)
            {
                priority_points.pop();
            }
            max_distance = priority_points.top().distance_;
        }
    }
    return true;
}
/*
const KD_tree_node* KD_tree::explore_to_leaf(const vector<double>& query_point, const KD_tree_node* root,  node_minPQ& minPQ)
{
    KD_tree_node* unexplored_node;   // unexplored KD_tree_node*
    KD_tree_node* cur_node = root;

    double split_value;
    double dim;

    while(cur_node!=NULL && !cur_node->is_leaf_)   //currentNode->left!=NULL && currentNode->right!=NULL signifies that currentNode is not a leaf
    {
        //partition dimension and value;
        dim = cur_node->split_dim_;
        split_value = cur_node->split_value_;


        // go to a child and preserve the other
        if(query_point[dim] < split_value)
        {
            unexplored_node = cur_node->right;
            cur_node = cur_node->left;
        }
        else
        {
            unexplored_node = cur_node->left;
            cur_node = cur_node->right;
        }

        if (unexplored_node!=NULL)
        {
            minPQ.push(node_bind(unexplored_node, fabs(query_point[unexplored_node->split_dim_]-unexplored->data_index_[unexplored_node->split_dim_]));
        }

    }

    return cur_node;
}
*/
/*
bool KD_tree::bbf_kNN_query(const vector<double> & query_point, const int K,
               vector<int> & indices,
               vector<double> & squared_distances, size_t max_epoch) const;
{

    size_t epoch = 0;

    KD_tree_node* cur_node = root_;

    priority_queue<distance_index> priority_points;

    double currentBest = numeric_limits<double>::max();

    double dist = 0;

    priority_points.push(root_, 0);

    while(!priority_points.empty() && epoch < maxEpoch)
    {
        cur_node = priority_points.top().key;
        priority_points.pop();

        //find leaf node and push unexplored to minPQ
        currentNode = exploreToLeaf(key, currentNode, priority_points);

        for(int i = 0; i<cur_node->)

    dist = Distance(currentNode->key, key);

    if(dist < currentBest)
    {
        if(kNearestPQ.size()==k)
        {
            //update the currentBest;
            kNearestPQ.pop();
            kNearestPQ.enqueue(currentNode, Distance(currentNode->key, key));
            currentBest = kNearestPQ.best();
        }
        else if(kNearestPQ.size() == k-1)
        {
            kNearestPQ.enqueue(currentNode, Distance(currentNode->key, key));
            currentBest = kNearestPQ.best();
        }
        else
        {
            kNearestPQ.enqueue(currentNode, Distance(currentNode->key, key));
        }
    }
    }
    ++epoch;

    while(!kNearestPQ.empty())
    {
        result.push_back((kNearestPQ.dequeueMin())->value);
    }

    return result;


}*/
//For brute-force comparison
double KD_tree::distance_sq(const vector<double> & data0, const vector<double> & data1)
{
    assert(data0.size() == data1.size());
    double sq = 0.0;
    for (int i=0; i<data0.size(); i++)
    {
        double dif = data0[i]-data1[i];
        sq +=dif*dif;
    }
    return sq;
}


#endif /* KD_TREE_H */
