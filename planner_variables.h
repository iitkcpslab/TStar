#include<bits/stdc++.h>
using namespace std;
namespace planner_info
{
int sz=10,qtot=3;
float prefix_length,pre_len,suf_len=numeric_limits<float>::infinity();
vector<int> dest;

/***djkistra with qtot system states with loops,until,next,few corrections,integer hash improvements,consider node in product-graph for book-keeping***/
//correct the adj updations
struct point
{
    int x,y;
    point(int a,int b)
    {
        x=a;
        y=b;
    }
    point()
    {

    }
};

int nrow,ncol;
vector<struct point> suffix_path,prefix_path,cur_path;
vector<vector<struct point> > stored_suffix_paths;
unordered_map<int,int> fin_aut_state;
int dir=9,total_system_states=0,comp_red=0,tot_comp=0,grid_sz;
int nb[][2]={{0,0},{0,-1},{-1,0},{1,0},{0,1},{-1,-1},{1,-1},{-1,1},{1,1}};


unordered_map<long long int,int> grid;
vector< vector< vector< vector<int> > > > trans(qtot,vector< vector< vector<int> > >(qtot) );
vector< vector< vector< vector<int> > > > negtrans(qtot,vector< vector< vector<int> > >(qtot) );
vector<int> neg_trans_to_neighbour;
vector<vector<int> > pos_system_state;
vector<unordered_map<int,float> > adj;
vector<unordered_map<int,int> > updated;

unordered_map<string,vector<int> > prop;
vector<vector<int> > qsystate;
vector< vector< vector<int> > > qsytrans(qtot,vector< vector<int> >(qtot));
unordered_map<long long int,int> systate_no;
vector<int> path;
vector<int> prevpath;
vector<vector< vector<int> > > prop_sys_states(1000);


long long int key(vector<int> &v)
{
    long long int val;
    if(v.size()==3) 
    val = v[0]*sz*sz+v[1]*sz+v[2];
    else
    val = v[0]*1000+v[1];

    return val;
}
/*unsigned long long int key(vector<int> &vec)
{
    size_t b=0;
    hash_range(b,vec.begin(), vec.end());
    //cout<<"code="<<b<<"\n";
    return b;
}*/


string conv_vec_to_string(vector<int> &vec)
{
    string s="";
    for(int i=0;i<vec.size();i++)
    {
        stringstream ss;
        ss << vec[i];
        s= s+ss.str();
        s+=",";
    }
    return s;
}


struct trans_system_node
{
    float g,h,f;
    point coord;
};


void printvec(vector<int> &vec)
{
    for(int i=0;i<vec.size();i++)
        cout<<vec[i];
}

bool valid(int y,int x)
{
    if(x>=0 && x<ncol && y>=0 && y<nrow && grid[y*sz+x]!=-1)
        return 1;
    return 0;
}

//high level nodes
struct prod_graph_node
{
    prod_graph_node* par;
    prod_graph_node* child;
    point coord;
    float g,h,f;
    bool visit;
    int state,pstate;
};

//actual nodes in path
struct trans_system_path_node
{
    trans_system_path_node* par;
    trans_system_path_node* child;
    point coord;
    float g,h,f;
    bool visit;
    int state;
};


trans_system_path_node* newtrans_system_path_node(int x,int y,int s)
{
    trans_system_path_node* nd = new trans_system_path_node;
    nd->coord.x=x;
    nd->coord.y=y;
    nd->f = 0.0;
    nd->g = 0.0;
    nd->h = 0.0;
    nd->visit = 0;
    nd->state = s;
    nd->par=NULL;
    return nd;
}

prod_graph_node* new_prod_graph_node(int p,int s,int y,int x)
{
    //p-system state s- automata state
    prod_graph_node* nd = new prod_graph_node;
    nd->coord.y=y;
    nd->coord.x=x;
    nd->f = 0.0;
    nd->g = 0.0;
    nd->h = 0.0;
    nd->visit = 0;
    nd->pstate = p;
    nd->state = s;
    nd->par=NULL;
    return nd;
}

struct prod_graph_comparator
{
    bool operator()(prod_graph_node *a , prod_graph_node *b)
    {
        return (a->f > b->f);
    }
};

float cal_heuristic_cost(point a,point b)
{
    float mx = max(abs(a.y-b.y),abs(a.x-b.x));
    float mn = min(abs(a.y-b.y),abs(a.x-b.x));
    return mn*1.5+(mx-mn);

}

//priority_queue for suffixes
priority_queue<prod_graph_node*,vector<prod_graph_node*>,prod_graph_comparator> loopf;

struct trans_system_node_cost_comparator
{
    bool operator()(trans_system_node a , trans_system_node b)
    {
        return (a.f > b.f);
    }
};

struct trans_system_path_node_cost_comparator
{
    bool operator()(trans_system_path_node *a , trans_system_path_node *b)
    {
        return (a->f > b->f);
    }
};

}
