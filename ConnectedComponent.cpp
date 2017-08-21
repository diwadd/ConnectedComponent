#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <random>
#define FOR(i, N) for(int i = 0; i < N; i++)

#define NUMBER_OF_EDGES_PER_VERTEX 4

const int i_shifts[NUMBER_OF_EDGES_PER_VERTEX] = {1, 0, -1,  0};
const int j_shifts[NUMBER_OF_EDGES_PER_VERTEX] = {0, 1,  0, -1};
#define DOWN_NEIGHBOUR  0  // Corresponds to shift [ 1,  0]
#define RIGHT_NEIGHBOUR 1  // Corresponds to shift [ 0,  1]
#define UP_NEIGHBOUR    2  // Corresponds to shift [-1,  0]
#define LEFT_NEIGHBOUR  3  // Corresponds to shift [ 0, -1]

using namespace std;

class Edge;
class Vertex;

typedef vector<int> vi;
typedef vector<vector<int>> vvi;

typedef vector<Vertex> vx;
typedef vector<Vertex*> vxp;
typedef vector<vector<Vertex>> vvx;
typedef vector<vector<Vertex*>> vvxp;


class Edge {
    public:
        Vertex* m_u;
        Vertex* m_v;

        Edge() { m_u = nullptr; m_v = nullptr; };

        Edge& operator=(Edge e);

};

Edge& Edge::operator=(Edge e) {

    //cerr << "In Edge copy assignment operator, operator=" << endl;
    m_u = e.m_u;
    m_v = e.m_v;

    return *this;
}


class Vertex {
    public:
        int m_key; // The key should be unique for every vertex.
        int m_val; // This is the matrix element != 0;

        // Fields for disjoint-set data structure
        Vertex* m_p; //parent
        int m_r; // rank

        // Each vertex can have at most 4 neighbours (edges).
        Vertex* m_le;
        Vertex* m_re;
        Vertex* m_ue;
        Vertex* m_de;

        Vertex(): m_key(0), m_val(0), m_p(nullptr), m_r(0), m_le(nullptr), m_re(nullptr), m_ue(nullptr), m_de(nullptr) {}

        Vertex(int ik, int iv, Vertex *ip, int ir, Vertex* le = nullptr, Vertex* re = nullptr, Vertex* ue = nullptr, Vertex* de = nullptr): 
                m_key(ik), m_val(iv), m_p(ip), m_r(ir), m_le(le), m_re(re), m_ue(ue), m_de(de) {}

        Vertex& operator=(Vertex v);


};


Vertex& Vertex::operator=(Vertex v) {
    
    //cerr << "In Vertex copy assignment operator, operator=" << endl;
    m_key = v.m_key;
    m_val = v.m_val;
    m_p = v.m_p;
    m_r = v.m_r;

    m_le = v.m_le;
    m_re = v.m_re;
    m_ue = v.m_ue;
    m_de = v.m_de;

    return *this;
}


std::ostream& operator<<(std::ostream& os, const Vertex& v) {

	string s = "key: " + to_string(v.m_key) + ", val: " + to_string(v.m_val);

    os << s;
    return os;
}


void print_vector(vi &v) {

    for(int i = 0; i < v.size(); i++)
        cerr << v[i]  << " ";
    cerr << endl;

}


void print_matrix(vvi &matrix) {
    
    int S = matrix.size();
    for(int i = 0; i < S; i++) {
        for(int j = 0; j < S; j++)
            fprintf(stderr, "%2d ", matrix[i][j]);
        fprintf(stderr, "\n");
    }

}


void print_matrix(vvxp &pointer_vertex_matrix) {
    
    int S = pointer_vertex_matrix.size();
    for(int i = 0; i < S; i++) {
        for(int j = 0; j < S; j++)

            if (pointer_vertex_matrix[i][j]!= nullptr)
                fprintf(stderr, "%2d ", (*pointer_vertex_matrix[i][j]).m_val);
            else
                fprintf(stderr, " 0 ");
        fprintf(stderr, "\n");
    }

}


void array_form_to_matrix(vi &matrix_array_form, 
                          vvi &matrix,
                          int &S) {

    for(int i = 0; i < S; i++)
        for(int j = 0; j < S; j++)
            matrix[i][j] = matrix_array_form[i*S + j];
}


void permute_matrix(vvi &matrix_p, 
                    vvi &matrix, 
                    vi &perm) {

    int S = matrix.size();
    for(int i = 0; i < S; i++)
        for(int j = 0; j < S; j++)
            matrix_p[i][j] = matrix[perm[i]][perm[j]];
}


void permute_pointer_vertex_matrix(vvxp &pointer_vertex_matrix, 
                                   vvxp &pointer_vertex_matrix_permuted, 
                                   vi &perm) {

    int S = pointer_vertex_matrix.size();
    for(int i = 0; i < S; i++)
        for(int j = 0; j < S; j++)
            pointer_vertex_matrix_permuted[i][j] = pointer_vertex_matrix[perm[i]][perm[j]];
}


void permute_pointer_vertex_matrix(vvxp &pointer_vertex_matrix, 
                                   int &pi,
                                   int &pj) {

    // We assume that in the permutation vector only
    // two entries are permuted, e.g., 
    // [p1, p2,... pi,... pj,..., pS] -> [p1, p2,... pj,... pi,..., pS]

    int S = pointer_vertex_matrix.size();
    if (pi < 0 || pi >= S) {
        cerr << "ERROR! In two_permute_pointer_vertex_matrix, pi out of bounds!" << endl;
        exit(0);
    }

    if (pj < 0 || pj >= S) {
        cerr << "ERROR! In two_permute_pointer_vertex_matrix, pj out of bounds!" << endl;
        exit(0);
    }

    //Swap column
    for(int i = 0; i < S; i++) {
        Vertex *t = pointer_vertex_matrix[pj][i];
        pointer_vertex_matrix[pj][i] = pointer_vertex_matrix[pi][i];
        pointer_vertex_matrix[pi][i] = t;
    }

    //Swap row
    for(int i = 0; i < S; i++) {
        Vertex *t = pointer_vertex_matrix[i][pj];
        pointer_vertex_matrix[i][pj] = pointer_vertex_matrix[i][pi];
        pointer_vertex_matrix[i][pi] = t;
    }

}



void initialise_vertex_array(vvi &matrix,
                             vx &vertex_array,
                             vvxp &pointer_vertex_matrix) {

    int n_v = 0;
    int S = matrix.size();

    //Count the number of non zero elements.
    for(int i = 0; i < S; i++)
        for(int j = 0; j < S; j++)
            if (matrix[i][j] != 0)
                n_v++;
                
    vertex_array.resize(n_v, Vertex());
    int vkc = 0; // vertex key counter
    for(int i = 0; i < S; i++) {
        for(int j = 0; j < S; j++) {

            if (matrix[i][j] == 0)
                continue;

            Vertex v(vkc, matrix[i][j], nullptr, 0);
            vertex_array[vkc] = v;

            pointer_vertex_matrix[i][j] = &vertex_array[vkc];

            vkc++;
        }
    }
}


void initialise_vertex_matrix(vvi &matrix,
                              vvx &vertex_matrix) {

    int S = matrix.size();
    int vkc = 0; // vertex key counter
    for(int i = 0; i < S; i++) {
        for(int j = 0; j < S; j++) {

            Vertex v(vkc, matrix[i][j], nullptr, 0);
            vertex_matrix[i][j] = v;

            vkc++;
        }
    }
}

void set_neighbours(int &i, int &j, int &S, vvxp &pointer_vertex_matrix) {

    // Sets the neighbours (edges) of vertex [i, j].

    for(int k = 0; k < NUMBER_OF_EDGES_PER_VERTEX; k++) {

        int i_new = i + i_shifts[k];
        int j_new = j + j_shifts[k];

        if (i_new < 0 || i_new >= S)
            continue;

        if (j_new < 0 || j_new >= S)
            continue;

        //cerr << "S: " << S << " i_new: " << i_new << " j_new: " << j_new << endl;

        // Check if the potential neighbour is a vertex?
        if (pointer_vertex_matrix[i_new][j_new] == nullptr)
            continue;

        if (k == DOWN_NEIGHBOUR)
            pointer_vertex_matrix[i][j]->m_de = pointer_vertex_matrix[i_new][j_new];
        else if (k == RIGHT_NEIGHBOUR) // right neighbours
            pointer_vertex_matrix[i][j]->m_re = pointer_vertex_matrix[i_new][j_new];
        else if (k == UP_NEIGHBOUR) // down neighbours
            pointer_vertex_matrix[i][j]->m_ue = pointer_vertex_matrix[i_new][j_new];
        else if (k == LEFT_NEIGHBOUR) // left neighbours
            pointer_vertex_matrix[i][j]->m_le = pointer_vertex_matrix[i_new][j_new];
        else {}

    }

}


void set_edges(vvxp &pointer_vertex_matrix,
               int &pi,
               int &pj) {

    // Set the edges to vertexes along rows and columns given by pi and pj, O(2xS).

    int S = pointer_vertex_matrix.size();
    vector<int> rc = {pi, pj};

    for(int i = 0; i < rc.size(); i++) {

        int rc_index = rc[i]; // Index of the row/column that we are changing.
        for(int j = 0; j < S; j++) {

            //cerr << "rc_index: " << rc_index << " j: " << j << endl;

            // Scanning horizontally (left -> right).
            if (pointer_vertex_matrix[rc_index][j] != nullptr) {           

                pointer_vertex_matrix[rc_index][j]->m_de = nullptr;
                pointer_vertex_matrix[rc_index][j]->m_re = nullptr;
                pointer_vertex_matrix[rc_index][j]->m_ue = nullptr;
                pointer_vertex_matrix[rc_index][j]->m_le = nullptr;

                set_neighbours(rc_index, j, S, pointer_vertex_matrix);

                // The row is fixed. We are scanning through columns.
                // We have to update the edges of the vertexes 
                // above and below the current vertex.
                int i_up = rc_index - 1;
                int i_down = rc_index + 1;
                if (i_up > 0 && pointer_vertex_matrix[i_up][j] != nullptr)
                    pointer_vertex_matrix[i_up][j]->m_de = pointer_vertex_matrix[rc_index][j];

                if (i_down < S && pointer_vertex_matrix[i_down][j] != nullptr)
                    pointer_vertex_matrix[i_down][j]->m_ue = pointer_vertex_matrix[rc_index][j];

            } else {

                // If there is no vertex at [rc_index][j] we still have to update
                // edges of the neighbours.
                int i_up = rc_index - 1;
                int i_down = rc_index + 1;
                if (i_up > 0 && pointer_vertex_matrix[i_up][j] != nullptr)
                    pointer_vertex_matrix[i_up][j]->m_de = nullptr;

                if (i_down < S && pointer_vertex_matrix[i_down][j] != nullptr)
                    pointer_vertex_matrix[i_down][j]->m_ue = nullptr;

            }

            // Scanning vertically (up -> down).
            if (pointer_vertex_matrix[j][rc_index] != nullptr) {

                pointer_vertex_matrix[j][rc_index]->m_de = nullptr;
                pointer_vertex_matrix[j][rc_index]->m_re = nullptr;
                pointer_vertex_matrix[j][rc_index]->m_ue = nullptr;
                pointer_vertex_matrix[j][rc_index]->m_le = nullptr;

                set_neighbours(j, rc_index, S, pointer_vertex_matrix);

                // The column is fixed. We are scanning through rows.
                // We have to update the edges of the vertexes 
                // on the left and right the current vertex.
                int i_left = rc_index - 1;
                int i_right = rc_index + 1;
                if (i_left > 0 && pointer_vertex_matrix[j][i_left] != nullptr)
                    pointer_vertex_matrix[j][i_left]->m_re = pointer_vertex_matrix[j][rc_index];

                if (i_right < S && pointer_vertex_matrix[j][i_right] != nullptr)
                    pointer_vertex_matrix[j][i_right]->m_le = pointer_vertex_matrix[j][rc_index];

            } else {

                // If there is no vertex at [rc_index][j] we still have to update
                // edges of the neighbours.
                int i_left = rc_index - 1;
                int i_right = rc_index + 1;
                if (i_left > 0 && pointer_vertex_matrix[j][i_left] != nullptr)
                    pointer_vertex_matrix[j][i_left]->m_re = nullptr;

                if (i_right < S && pointer_vertex_matrix[j][i_right] != nullptr)
                    pointer_vertex_matrix[j][i_right]->m_le = nullptr;

            }

        }
    }


}


void set_edges(vvxp &pointer_vertex_matrix) {

    // Sets edges of all the vertexes, O(S^2).

    int S = pointer_vertex_matrix.size();
    for(int i = 0; i < S; i++) {
        for(int j = 0; j < S; j++) {

            //cerr << "i: " << i << " j: " << j << endl;

            // Check if we're pointing to a vertex?
            if (pointer_vertex_matrix[i][j] == nullptr)
                continue;

            pointer_vertex_matrix[i][j]->m_de = nullptr;
            pointer_vertex_matrix[i][j]->m_re = nullptr;
            pointer_vertex_matrix[i][j]->m_ue = nullptr;
            pointer_vertex_matrix[i][j]->m_le = nullptr;

            set_neighbours(i, j, S, pointer_vertex_matrix);
        }
    }

}


void print_vertex_information(Vertex *v) {

    if (v == nullptr)
        cerr << "Pointer is not pointing to any vertex." << endl;
    else {

        cerr << endl << "Basic vertex info: " << endl;
        cerr << v << endl;
        cerr << *v << endl;
        cerr << &v << endl;
        cerr << &(*v) << endl;

        cerr << "Vertex neighbours:" << endl;
        if (v->m_ue != nullptr)        
            cerr << "UP: " << *(v->m_ue) << endl;
        if (v->m_de != nullptr)
            cerr << "DOWN: " << *(v->m_de) << endl;
        if (v->m_re != nullptr) 
            cerr << "RIGHT: " << *(v->m_re) << endl;
        if (v->m_le != nullptr) 
            cerr << "LEFT: " << *(v->m_le) << endl;

    }

}


void make_set(Vertex* v) {
    v->m_p = v;
    v->m_r = 0;
}


void link(Vertex* x, Vertex* y) {

    if (x->m_r > y->m_r)
        y->m_p = x;
    else {
        x->m_p = y;
        if (x->m_r == y->m_r)
            y->m_r = y->m_r + 1;    
    }

}


Vertex* find_set(Vertex* x) {
    if (x != x->m_p)
        x->m_p = find_set(x->m_p);
    return x->m_p;
}


// Union id the disjoint set data structure.
void union_dsds(Vertex *x, Vertex *y) {
    link(find_set(x), find_set(y));
}


void connected_components(vx &vertex_array) {

    int n_v = vertex_array.size();
    for(int i = 0; i < n_v; i++)
        make_set(&vertex_array[i]);

    Vertex* u = nullptr;
    Vertex* v = nullptr;
    for(int i = 0; i < n_v; i++) {

        u = &vertex_array[i];       

        v = vertex_array[i].m_ue;
        if (v != nullptr)   
            if (find_set(u) != find_set(v))
                union_dsds(u, v);

        v = vertex_array[i].m_de;
        if (v != nullptr)
            if (find_set(u) != find_set(v))
                union_dsds(u, v);

        v = vertex_array[i].m_re;
        if (v != nullptr)
            if (find_set(u) != find_set(v))
                union_dsds(u, v);

        v = vertex_array[i].m_le;
        if (v != nullptr)
            if (find_set(u) != find_set(v))
                union_dsds(u, v);
    }
}


bool same_component(Vertex *u, Vertex *v) {

    if (u == nullptr || v == nullptr)
        return false;

    if (find_set(u) == find_set(v))
        return true;
    else
        return false;
}


double count_score(vx &vertex_array) {

    //cerr << "In count_score" << endl;

    int n_v = vertex_array.size();
    vi sum_array(n_v, 0);
    vi count_array(n_v, 0);

    //cerr << "Looping vertexes." << endl;
    for(int i = 0; i < n_v; i++) {
        Vertex* p = find_set(&vertex_array[i]);

        sum_array[p->m_key] = sum_array[p->m_key] + vertex_array[i].m_val;
        count_array[p->m_key] = count_array[p->m_key] + 1;
    }

    //cerr << "Getting vertex scores." << endl;
    vector<double> scores(n_v, 0.0);
    for(int i = 0; i < n_v; i++)
        scores[i] = (double)sum_array[i]*sqrt((double)count_array[i]);

    //cerr << " Finding max score." << endl;
    double s = (*max_element(scores.begin(), scores.end()));

    return s;

}


void apply_n_random_permutations(int &n,
                                 vvxp &pointer_vertex_matrix,
                                 vi &final_perm) {

    int S = pointer_vertex_matrix.size();   
    random_device rd;
    mt19937 gen(rd());

    uniform_int_distribution<> uni_int(0, S - 1);

    for(int i = 0; i < n; i++) {
        int pi = uni_int(gen);
        int pj = uni_int(gen);

        if (pi == pj)
            continue;

        permute_pointer_vertex_matrix(pointer_vertex_matrix, pi, pj);
        set_edges(pointer_vertex_matrix);

        int t = final_perm[pi];
        final_perm[pi] = final_perm[pj];
        final_perm[pj] = t;

    }

}


void metropolis(double &T,
                double &alpha,
                vx &vertex_array,
                vvxp &pointer_vertex_matrix,
                vi &final_perm) {

    cerr << "Running Monte Carlo with Metropolis..." << endl;

    int S = pointer_vertex_matrix.size();   
    random_device rd;
    mt19937 gen(rd());

    double min_cost = numeric_limits<double>::max();

    vi current_perm(S, 0);

    for(int i = 0; i < S; i++)
        current_perm[i] = i;

    uniform_int_distribution<> uni_int(0, S - 1);
    uniform_real_distribution<double> uni_real(0.0, 1.0);

    while(T > 1.0) {

        //cerr << "While beginning." << endl;
        //cerr << "Counting score" << endl;
        connected_components(vertex_array);
        double cost_1 = -count_score(vertex_array);
        //cerr << "cost_1: " << cost_1 << endl;
        //cerr << "Score counted" << endl;

        int pi = uni_int(gen);
        int pj = uni_int(gen);

        if (pi == pj)
            continue;
        
        // Move pointer_vertex_matrix to new state.
        //cerr << "New state generated." << endl;
        permute_pointer_vertex_matrix(pointer_vertex_matrix, pi, pj);
        set_edges(pointer_vertex_matrix);
        //cerr << "New state generated." << endl;

        connected_components(vertex_array);
        double cost_2 = -count_score(vertex_array);
        //cerr << "cost_2: " << cost_2 << endl;

        // Accept new state
        if (cost_2 > cost_1) {
            int t = current_perm[pi];
            current_perm[pi] = current_perm[pj];
            current_perm[pj] = t;

            if (min_cost > cost_2) {
                min_cost = cost_2;
                final_perm = current_perm;
            }

        // Accept new state with probability exp( -(cost_2 - cost_1)/T ).
        } else if ( uni_real(gen) < exp( -(cost_2 - cost_1)/T ) ) {
            int t = current_perm[pi];
            current_perm[pi] = current_perm[pj];
            current_perm[pj] = t;

            if (min_cost > cost_2) {
                min_cost = cost_2;
                final_perm = current_perm;
            }
        // Reject new state.
        } else {
            permute_pointer_vertex_matrix(pointer_vertex_matrix, pj, pi);
            set_edges(pointer_vertex_matrix);
        }

        T = T*alpha;


    }


}


class ConnectedComponent {
public:
    vector<int> permute(vector<int> matrix_array_form) {

        int S = (int)sqrt(matrix_array_form.size());
        vi final_perm(S, 0);

        for(int i = 0; i < S; i++)
            final_perm[i] = i;

        //print_vector(final_perm);

        // S x S matrix
        vvi matrix(S, vi(S, 0));
        array_form_to_matrix(matrix_array_form,
                             matrix,
                             S);

        cerr << "Original matrix: " << endl;
        //print_matrix(matrix);

        vx vertex_array;
        vvxp pointer_vertex_matrix(S, vxp(S, nullptr));
        //vvxp pointer_vertex_matrix_support(S, vxp(S, nullptr));
        initialise_vertex_array(matrix,
                                vertex_array,
                                pointer_vertex_matrix);
        //initialise_vertex_array(matrix,
        //                        vertex_array,
        //                        pointer_vertex_matrix_support);
        set_edges(pointer_vertex_matrix);
        //set_edges(pointer_vertex_matrix_support);

        cerr << "pointer_vertex_matrix" << endl;
        //print_matrix(pointer_vertex_matrix);

        //cerr << "pointer_vertex_matrix_support" << endl;
        //print_matrix(pointer_vertex_matrix_support);

        int iv = 1;
        int jv = 5;
        print_vertex_information(pointer_vertex_matrix[iv][jv]);
        cerr << endl;

        vi perm = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
        //vvxp pointer_vertex_matrix_permuted(S, vxp(S, nullptr));
        //permute_pointer_vertex_matrix(pointer_vertex_matrix, 
        //                              pointer_vertex_matrix_permuted, 
        //                              perm);
        //set_edges(pointer_vertex_matrix_permuted);

        //cerr << "Permutation (pointer matrix): " << endl;        
        //print_matrix(pointer_vertex_matrix_permuted);
        //print_vertex_information(pointer_vertex_matrix_permuted[iv][jv]);
        //cerr << endl;


        cerr << "Small permutation (pointer matrix): " << endl;
        //int pi = 0;
        //int pj = 1;
        //permute_pointer_vertex_matrix(pointer_vertex_matrix, pi, pj);
        //set_edges(pointer_vertex_matrix, pi, pj);
        //pi = 5;
        //pj = 8;
        //permute_pointer_vertex_matrix(pointer_vertex_matrix, pi, pj);
        //set_edges(pointer_vertex_matrix, pi, pj);
        //set_edges(pointer_vertex_matrix);

        //print_matrix(pointer_vertex_matrix);

        //iv = 1; jv = 1;
        //print_vertex_information(pointer_vertex_matrix[iv][jv]);

        //iv = 1; jv = 4;
        //print_vertex_information(pointer_vertex_matrix[iv][jv]);

        //iv = 3; jv = 6;
        //print_vertex_information(pointer_vertex_matrix[iv][jv]);

        //iv = 1; jv = 5;
        //print_vertex_information(pointer_vertex_matrix[iv][jv]);

        //int n = 10;
        //apply_n_random_permutations(n, pointer_vertex_matrix, final_perm);
        //print_vector(final_perm);

        //cerr << "pointer_vertex_matrix after apply_n_random_permutations" << endl;
        //print_matrix(pointer_vertex_matrix);

        //cerr << "pointer_vertex_matrix_support" << endl;
        //print_matrix(pointer_vertex_matrix_support);

        //vvxp pointer_vertex_matrix_permuted(S, vxp(S, nullptr));
        //permute_pointer_vertex_matrix(pointer_vertex_matrix_support, 
        //                              pointer_vertex_matrix_permuted, 
        //                              final_perm);

        //cerr << "pointer_vertex_matrix_permuted" << endl;
        //print_matrix(pointer_vertex_matrix_permuted);
        

        double T = 100.0;
        double alpha = 0.9;
        metropolis(T,
                   alpha,
                   vertex_array,
                   pointer_vertex_matrix,
                   final_perm);


        //cerr << &vertex_array[0] << endl;
        //cerr << vertex_array[0].m_p << endl;
        //connected_components(vertex_array);
        //cerr << vertex_array[0].m_p << endl;

        //Vertex* u = pointer_vertex_matrix[7][0];
        //Vertex* v = pointer_vertex_matrix[9][0];

        vvxp pointer_vertex_matrix_max_score(S, vxp(S, nullptr));
        initialise_vertex_array(matrix,
                                vertex_array,
                                pointer_vertex_matrix_max_score);
        set_edges(pointer_vertex_matrix_max_score);

        vvxp pointer_vertex_matrix_support(S, vxp(S, nullptr));
        initialise_vertex_array(matrix,
                                vertex_array,
                                pointer_vertex_matrix_support);
        set_edges(pointer_vertex_matrix_support);


        permute_pointer_vertex_matrix(pointer_vertex_matrix_support, 
                                      pointer_vertex_matrix_max_score, 
                                      final_perm);
        set_edges(pointer_vertex_matrix_max_score);

        connected_components(vertex_array);

        //cerr << "The vertexes are the same subgraph: " << same_component(u, v) << endl;
        double s = count_score(vertex_array);
        fprintf(stderr, "Internal score: %5.10f\n", s);

        vector<int> ret(S);
        for (int i = 0; i < S; ++i) {
            ret[i] = S - 1 - i;
        }
        return ret;
    }
};
// -------8<------- end of solution submitted to the website -------8<-------

template<class T> void getVector(vector<T>& v) {
    for (int i = 0; i < v.size(); ++i)
        cin >> v[i];
}

int main() {
    ConnectedComponent cc;
    int M;
    cin >> M;
    vector<int> matrix(M);
    getVector(matrix);
    
    vector<int> ret = cc.permute(matrix);
    cout << ret.size() << endl;
    for (int i = 0; i < (int)ret.size(); ++i)
        cout << ret[i] << endl;
    cout.flush();
}
