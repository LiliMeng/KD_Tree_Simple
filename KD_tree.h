#ifndef KD_tree_H
#define KD_tree_H

#include <vector>
#include <assert.h>
#include <queue>
#include <stack>
#include <set>
#include <algorithm>

using namespace std;



class bounding_box
{
private:
    vector<double> min_point_;
    vector<double> max_point_;
public:
    bounding_box()
    {

    }
    bounding_box(const vector<double> & min_point, const vector<double> & max_point)
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

        for ( unsigned int i=0; i<dim; ++i ) {
            double x_min = min_point_[i];
            double x_max = max_point_[i];
            double x = point[i];
            if ( x < x_min ) {
                sum_sq += ( x_min - x ) * (x_min - x);
            }
            else if ( x > x_max ) {
                sum_sq += ( x_max - x ) * (x_max - x);
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
    bounding_box box_;

    KD_tree_node *left_node_;
    KD_tree_node *right_node_;
    vector<int> data_index_;

    KD_tree_node()
    {
        split_dim_ = -1;
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
private:
    // build and store data
    vector<vector<double> > data_;
    unsigned long dim_;
    unsigned bin_num_;
    unsigned max_level_;
    KD_tree_node *root_;


public:
    KD_tree();
    ~KD_tree();

    bool create_tree(const vector<vector<double> > & data, unsigned bin_num, unsigned max_level);
    bool query(const vector<double> & query_point, const int K,
               vector<int> & indices,
               vector<double> & squared_distances) const;

private:

    bool build_tree(const vector<int> & index, KD_tree_node* & node, unsigned level);

    const KD_tree_node * get_leaf_node(const vector<double> & query_point, const KD_tree_node * node) const;

    // query update from current node
    bool query(const vector<double> & query_point, const int K,
               stack<KD_tree_node *> & nodes,
               set<KD_tree_node *> & visited_nodes,
               priority_queue<distance_index> & priority_points,
               double & max_distance) const;

private:

    // check points in a left node
    bool explore_leaf_node(const vector<double> & query_point, const int K,
                           priority_queue<distance_index> & priority_points,
                           double & max_distance, const KD_tree_node * cur_node) const;


    bounding_box build_bounding_box(const vector<int> & indices) const;

public:

    // squared of distance
    static double distance_sq(const vector<double> & data0, const vector<double> & data1);


};

KD_tree::KD_tree()
{
    root_ = NULL;

}
KD_tree::~KD_tree()
{

}

bool KD_tree::create_tree(const vector<vector<double> > & data, unsigned bin_num, unsigned max_level)
{
    assert(root_ == NULL);

    data_ = data;
    bin_num_ = bin_num;
    max_level_ = max_level;
    dim_ = data.front().size();

    vector<int> indices;
    for (int i = 0; i<data_.size(); i++) {
        indices.push_back(i);
    }
    root_ = new KD_tree_node();
    root_->box_ = this->build_bounding_box(indices);
    this->build_tree(indices, root_, 0);

    return true;
}

/*
struct KD_tree_node
{
    int split_dim_;
    double split_value_;
    bool is_leaf_;
    int level_;

    KD_tree_node *left_node_;
    KD_tree_node *right_node_;
    vector<int> data_index_;
};
*/




bool KD_tree::build_tree(const vector<int> & indices, KD_tree_node* & node, unsigned level)
{
    assert(node);
    if (level >= max_level_ || indices.size() <= bin_num_) {
        node->is_leaf_ = true;
        node->level_ = level;
        node->data_index_ = indices;
        node->box_ = this->build_bounding_box(indices);
        return true;
    }

    // randomly select split dimension
    int split_dim = rand()%dim_;
    node->split_dim_ = split_dim;
    vector<double> one_dim_values;
    for (int i = 0; i<indices.size(); i++) {
        one_dim_values.push_back(data_[indices[i]][split_dim]);
    }

    // median value as split value
    size_t n = one_dim_values.size() / 2;
    std::nth_element(one_dim_values.begin(), one_dim_values.begin()+n, one_dim_values.end());
    node->split_value_ = one_dim_values[n];

    node->box_ = this->build_bounding_box(indices);

    vector<int> left_indices;
    vector<int> right_indices;
    for (int i = 0; i<indices.size(); i++) {
        int index = indices[i];
        double v = data_[index][node->split_dim_];
        if (v <= node->split_value_) {
            left_indices.push_back(index);
        }
        else
        {
            right_indices.push_back(index);
        }
    }
    assert(left_indices.size() + right_indices.size() == indices.size());

    if (left_indices.size() != 0) {
        KD_tree_node *left_node = new KD_tree_node();
        this->build_tree(left_indices, left_node, level + 1);
        node->left_node_ = left_node;
    }
    if (right_indices.size() != 0) {
        KD_tree_node * right_node = new KD_tree_node();
        this->build_tree(right_indices, right_node, level + 1);
        node->right_node_ = right_node;
    }

    return true;
}

bounding_box KD_tree::build_bounding_box(const vector<int> & indices) const
{
    assert(indices.size() > 1);
    vector<double> min_point = data_[indices.front()];
    vector<double> max_point = min_point;

    for (unsigned int i=1; i < indices.size(); ++i ) {
        vector<double> point = data_[indices[i]];
        for ( unsigned int j=0; j < point.size(); ++j ) {
            if ( point[j] < min_point[j] ){
                min_point[j] = point[j];
            }
            if ( point[j] > max_point[j] ){
                max_point[j] = point[j];
            }
        }
    }
    return bounding_box(min_point, max_point);
}

const KD_tree_node * KD_tree::get_leaf_node(const vector<double> & query_point, const KD_tree_node * node) const
{
    assert(node);
    if (node->is_leaf_) {
        return node;
    }
    int dim = node->split_dim_;
    double split_value = node->split_value_;
    double value = query_point[dim];
    if (value <= split_value && node->left_node_) {
        return this->get_leaf_node(query_point, node->left_node_);
    }
    else if(node->right_node_)
    {
        return this->get_leaf_node(query_point, node->right_node_);
    }
    return NULL;
}

bool KD_tree::query(const vector<double> & query_point, const int K,
                    vector<int> & indices,
                    vector<double> & squared_distances) const
{
    assert(root_);
    assert(K < data_.size());

    double max_sq_distance = numeric_limits<double>::max();
    priority_queue<distance_index> distance_queue;

    stack<KD_tree_node *> candidate_nodes;
    KD_tree_node * cur_node = root_;
    // travel to leaf node
    while (cur_node) {
        if (cur_node != root_) {
            candidate_nodes.push(cur_node);
        }

        if (cur_node->is_leaf_) {
            break;
        }
        int dim = cur_node->split_dim_;
        double split_value = cur_node->split_value_;
        double value = query_point[dim];
        if (value <= split_value && cur_node->left_node_) {

            cur_node = cur_node->left_node_;
        }
        else if(cur_node->right_node_)
        {
            cur_node = cur_node->right_node_;
        }
    }

    set<KD_tree_node *> visited_nodes;
    this->query(query_point, K, candidate_nodes, visited_nodes, distance_queue, max_sq_distance);

    indices.resize(K);
    squared_distances.resize(K);
    assert(distance_queue.size() == K);

    int num = K-1;
    while (!distance_queue.empty()) {
        distance_index top = distance_queue.top();
        distance_queue.pop();
        indices[num] = top.index_;
        squared_distances[num] = top.distance_;
        num--;
    }

    return true;
}


bool KD_tree::query(const vector<double> & query_point, const int K,
                    stack<KD_tree_node *> & nodes,
                    set<KD_tree_node *> & visited_nodes,
                    priority_queue<distance_index> & priority_points,
                    double & max_distance) const
{
    while (!nodes.empty()) {
        KD_tree_node *cur_node = nodes.top();
        nodes.pop();
        // already visited node
        if (visited_nodes.find(cur_node) != visited_nodes.end()) {
            continue;
        }
        visited_nodes.insert(cur_node);
        if (cur_node->is_leaf_) {
            this->explore_leaf_node(query_point, K, priority_points, max_distance, cur_node);
        }
        else
        {
            // internal node
            if (cur_node->left_node_) {
                double dist_sq = cur_node->left_node_->box_.distance_sq(query_point);
                if (dist_sq < max_distance) {
                    nodes.push(cur_node->left_node_);
                }
            }
            if (cur_node->right_node_) {
                double dist_sq = cur_node->right_node_->box_.distance_sq(query_point);
                if (dist_sq < max_distance) {
                    nodes.push(cur_node->right_node_);
                }
            }
        }
    }
    return true;
}

bool KD_tree::explore_leaf_node(const vector<double> & query_point, const int K,
                                priority_queue<distance_index> & priority_points,
                                double & max_distance, const KD_tree_node * cur_node) const
{
    assert(cur_node->is_leaf_);

    for (int i = 0; i<cur_node->data_index_.size(); i++) {
        int index = cur_node->data_index_[i];
        const vector<double> & point = data_[index];
        double dist_sq = KD_tree::distance_sq(point, query_point);
        if (priority_points.size() < K || dist_sq < max_distance) {
            distance_index di(index, dist_sq);
            priority_points.push(di);
            if (priority_points.size() > K) {
                priority_points.pop();
            }
            max_distance = priority_points.top().distance_;
        }
    }
    return true;
}

double KD_tree::distance_sq(const vector<double> & data0, const vector<double> & data1)
{
    assert(data0.size() == data1.size());
    double sq = 0.0;
    for (int i = 0; i<data0.size(); i++) {
        double dif = data0[i] - data1[i];
        sq += dif * dif;
    }
    return sq;
}

#endif /* KD_TREE_H */
