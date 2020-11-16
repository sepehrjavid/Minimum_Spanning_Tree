#include <iostream>
#include <cmath>
#include <utility>
#include <vector>
#include <fstream>
#include<algorithm>

using namespace std;

int count_edges;

class Node {
public:
    float coordinate[2]{};
    int state;
    Node *left = nullptr;
    Node *right = nullptr;
    bool is_in_fragment = false;
    int name;
    int dimension = 0;


    Node(float x, float y, int s, int n) {
        coordinate[0] = x;
        coordinate[1] = y;
        state = s;
        name = n;

    }

};


class Edge {
public:
    Node *start = nullptr;
    Node *end = nullptr;
    float distance;

    Edge(Node *start, Node *end, float distance) : start(start), end(end), distance(distance) {}

    bool operator()(const Edge a, const Edge b) const {
        return a.distance > b.distance;
    }

};


class KDTree {
public:
    int k;
    Node *root = nullptr;
    Node *candidate_nearest_node = nullptr;
    float candidate_best_distance;

    explicit KDTree(int k) : k(k) {}


    static float find_distance_square(Node *point1, Node *point2) {
        float diff_x = pow(point2->coordinate[0] - point1->coordinate[0], 2);
        float diff_y = pow(point2->coordinate[1] - point1->coordinate[1], 2);
        return sqrt(diff_x + diff_y);
    }

    void insert(Node *current_node, Node *target, int depth) {
        if (root == nullptr) {
            root = target;
            target->dimension = depth % k;
            return;
        }

        int dimension = depth % k;
        if (target->coordinate[dimension] < current_node->coordinate[dimension]) {
            if (current_node->left == nullptr) {
                current_node->left = target;
                target->dimension = (depth + 1) % k;
                return;
            } else {
                return insert(current_node->left, target, depth + 1);
            }
        } else {
            if (current_node->right == nullptr) {
                current_node->right = target;
                target->dimension = (depth + 1) % k;
                return;
            } else {
                return insert(current_node->right, target, depth + 1);
            }
        }

    }


    void find_nearest_candidate(Node *current_node, Node *target) {
        float test_distance = find_distance_square(current_node, target);

        if (candidate_nearest_node == nullptr) {
            candidate_nearest_node = current_node;
            candidate_best_distance = test_distance;
        } else if (test_distance < candidate_best_distance) {
            candidate_best_distance = test_distance;
            candidate_nearest_node = current_node;
        }

        if (target->coordinate[current_node->dimension] < current_node->coordinate[current_node->dimension]) {
            if (current_node->left != nullptr) {
                find_nearest_candidate(current_node->left, target);
            }
        } else {
            if (current_node->right != nullptr) {
                find_nearest_candidate(current_node->right, target);
            }

        }
        float best_possible_distance = abs(
                current_node->coordinate[current_node->dimension] - target->coordinate[current_node->dimension]);
        if (best_possible_distance < candidate_best_distance) {
            if (target->coordinate[current_node->dimension] < current_node->coordinate[current_node->dimension]) {
                if (current_node->right != nullptr) {
                    find_nearest_candidate(current_node->right, target);
                }
            } else {
                if (current_node->left != nullptr) {
                    find_nearest_candidate(current_node->left, target);
                }
            }
        }
    }

    pair<Node *, Node *> find_dimension_min(Node *current_node, Node *parent, int target_dimension) {
        if (current_node->dimension == target_dimension) {
            if (current_node->left == nullptr) {
                return make_pair(parent, current_node);
            }
            return find_dimension_min(current_node->left, current_node, target_dimension);
        } else {
            pair<Node *, Node *> right_min = make_pair(nullptr, nullptr);
            pair<Node *, Node *> left_min = make_pair(nullptr, nullptr);
            if (current_node->right != nullptr) {
                right_min = find_dimension_min(current_node->right, current_node, target_dimension);
            }
            if (current_node->left != nullptr) {
                left_min = find_dimension_min(current_node->left, current_node, target_dimension);
            }
            if (right_min.second == nullptr && left_min.second != nullptr) {
                return left_min.second->coordinate[target_dimension] < current_node->coordinate[target_dimension]
                       ? left_min
                       : make_pair(parent, current_node);
            } else if (left_min.second == nullptr && right_min.second != nullptr) {
                return right_min.second->coordinate[target_dimension] < current_node->coordinate[target_dimension]
                       ? right_min
                       : make_pair(parent, current_node);
            } else if (right_min.second == nullptr && left_min.second == nullptr) {
                return make_pair(parent, current_node);
            } else {
                if (left_min.second->coordinate[target_dimension] <= right_min.second->coordinate[target_dimension]) {
                    return left_min.second->coordinate[target_dimension] < current_node->coordinate[target_dimension]
                           ? left_min
                           : make_pair(parent, current_node);
                } else {
                    return right_min.second->coordinate[target_dimension] < current_node->coordinate[target_dimension]
                           ? right_min
                           : make_pair(parent, current_node);
                }
            }
        }
    }

    void delete_node(Node *current_node, Node *parent, Node *target) {
        if (current_node == target) {
            if (current_node->right != nullptr) {
                pair<Node *, Node *> min = find_dimension_min(current_node->right, current_node,
                                                              current_node->dimension);
                Node *temp = min.second;

                delete_node(min.second, min.first, min.second);

                min.second->dimension = current_node->dimension;
                min.second->right = current_node->right;
                min.second->left = current_node->left;
                if (parent == nullptr) {
                    root = min.second;
                } else {
                    if (parent->left == current_node) {
                        parent->left = min.second;
                    } else if (parent->right == current_node) {
                        parent->right = min.second;
                    }
                }

            } else if (current_node->left != nullptr) {
                pair<Node *, Node *> min = find_dimension_min(current_node->left, current_node,
                                                              current_node->dimension);
                Node *temp = min.second;

                delete_node(min.second, min.first, min.second);

                min.second->dimension = current_node->dimension;
                min.second->right = current_node->left;
                min.second->left = nullptr;
                if (parent == nullptr) {
                    root = min.second;
                } else {
                    if (parent->left == current_node) {
                        parent->left = min.second;
                    } else if (parent->right == current_node) {
                        parent->right = min.second;
                    }
                }

            } else {
                if (parent != nullptr) {
                    if (parent->left == current_node) {
                        parent->left = nullptr;
                    } else if (parent->right == current_node) {
                        parent->right = nullptr;
                    }
                }
            }
        } else {
            if (target->coordinate[current_node->dimension] < current_node->coordinate[current_node->dimension]) {
                delete_node(current_node->left, current_node, target);
            } else {
                delete_node(current_node->right, current_node, target);
            }
        }
    }

};

Node *read_file_and_create_tree(KDTree *tree) {
    float x, y;
    int state, count;
    int i = 0;
    std::ifstream infile("PATH TO THE FILE");
    infile >> count;

    int temp_count = count - 1;
    count_edges = temp_count; // number of edges

    infile >> x;
    infile >> y;
    infile >> state;
    Node *first_node = new Node(x, y, state, i + 1);

    while (temp_count != i) {
        infile >> x;
        infile >> y;
        infile >> state;
        Node *new_node = new Node(x, y, state, i + 2);
        tree->insert(tree->root, new_node, 0);
        i++;
    }

    return first_node;
}


bool compare(const Edge a, const Edge b) {
    return a.distance >= b.distance;
}



int main() {
    float total_distance = 0;
    auto *kdtree = new KDTree(2);
    int i_MST = 0;
    vector<Edge> edges;
    vector<Edge> MST;
    bool (*foo)(const Edge, const Edge);
    foo = &compare;

    clock_t tStart = clock();

    Node *isolated_node = read_file_and_create_tree(kdtree);
    isolated_node->is_in_fragment = true;

    kdtree->find_nearest_candidate(kdtree->root, isolated_node);

    Edge new_edge = Edge(isolated_node, kdtree->candidate_nearest_node, kdtree->candidate_best_distance);
    edges.push_back(new_edge);
    make_heap(edges.begin(), edges.end(), foo);

    while (i_MST != count_edges) {
        Edge candidate_edge = edges.front();

        if (!candidate_edge.end->is_in_fragment) {
            Node *new_node_in_fragment = candidate_edge.end;
            new_node_in_fragment->is_in_fragment = true;
            MST.push_back(candidate_edge);
            i_MST++;
            total_distance = candidate_edge.distance + total_distance;

            kdtree->delete_node(kdtree->root, nullptr, new_node_in_fragment);

            kdtree->candidate_nearest_node = nullptr;
            kdtree->find_nearest_candidate(kdtree->root, new_node_in_fragment);

            new_edge = Edge(new_node_in_fragment, kdtree->candidate_nearest_node, kdtree->candidate_best_distance);
            edges.emplace_back(new_edge);
            push_heap(edges.begin(), edges.end(), foo);

        } else {
            Node *source_node = candidate_edge.start;

            pop_heap(edges.begin(), edges.end(), foo);
            edges.pop_back();

            kdtree->candidate_nearest_node = nullptr;
            kdtree->find_nearest_candidate(kdtree->root, source_node);

            new_edge = Edge(source_node, kdtree->candidate_nearest_node, kdtree->candidate_best_distance);
            edges.emplace_back(new_edge);
            push_heap(edges.begin(), edges.end(), foo);
        }
    }

    printf("Time taken: %.2fs\n", (float) (clock() - tStart) / CLOCKS_PER_SEC);
    cout << total_distance << endl;

    return 0;
}