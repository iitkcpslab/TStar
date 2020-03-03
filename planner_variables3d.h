#include<bits/stdc++.h>
#include<boost/functional/hash.hpp>
using namespace boost;
using namespace std;
namespace planner_info
{
int sz=10,qtot=3,lx,ly,lz;
float prefix_length,pre_len,suf_len=numeric_limits<float>::infinity();
vector<int> dest;

/***djkistra with qtot system states with loops,until,next,few corrections,integer hash improvements,consider node in product-graph for book-keeping***/
//correct the adj updations
struct point
{
    int x,y,z;
    point(int a,int b,int c)
    {
        x=a;
        y=b;
        z=c;
    }
    point()
    {

    }
};

vector<struct point> suffix_path,prefix_path,cur_path;
vector<vector<struct point> > stored_suffix_paths;
unordered_map<int,int> fin_aut_state;
int total_system_states=0,comp_red=0,tot_comp=0,grid_sz;


unordered_map<long long int,int> grid;
vector< vector< vector< vector<int> > > > trans(qtot,vector< vector< vector<int> > >(qtot) );
vector< vector< vector< vector<int> > > > negtrans(qtot,vector< vector< vector<int> > >(qtot) );
vector<int> neg_trans_to_neighbour;

vector<vector<int> > pos_system_state;
vector<vector<int> > pos_graph_node;

vector<unordered_map<int,float> > adj;
vector<unordered_map<int,int> > updated;

unordered_map<string,vector<int> > prop;
vector<vector<int> > qsystate;
vector< vector< vector<int> > > qsytrans(qtot,vector< vector<int> >(qtot));
unordered_map<long long int,int> systate_no;
vector<int> path;
vector<int> prevpath;
vector<vector< vector<int> > > prop_sys_states(1000);
int ndir=27;
int nb[][3]={{0,0,0},{0,0,-1},{0,0,1},{0,-1,0},{0,-1,-1},{0,-1,1},{-1,0,0},{-1,0,-1},{-1,0,1},{1,0,0},{1,0,-1},{1,0,1},{0,1,0},{0,1,-1},{0,1,1},{-1,-1,0},{-1,-1,-1},{-1,-1,1},
{1,-1,0},{1,-1,-1},{1,-1,1},{-1,1,0},{-1,1,-1},{-1,1,1},{1,1,0},{1,1,-1},{1,1,1}};


long long int key(vector<int> &v)
{
    long long int val=0;
    if(v.size()==2)
        return v[0]*1000+v[1];

    for(int i=v.size()-1;i>=0;i--)
        val= val+ v[v.size()-i-1]*pow(sz,i);
    return val;
}

long long int key(point pt)
{
    long long int val;
    val = pt.x*sz*sz+pt.y*sz+pt.z;
    return val;
}

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

bool equal_point(struct point point1,struct point point2)
{
    if(point1.x==point2.x && point1.y==point2.y && point1.z==point2.z)
        return 1;
    else
        return 0;
}
void printvec(vector<int> &vec)
{
    for(int i=0;i<vec.size();i++)
        cout<<vec[i]<<", ";
}

bool valid(int x,int y,int z)
{
    if(x>=0 && x<lx && y>=0 && y<ly && z>=0 && z<lz && grid[x*sz*sz+y*sz+z]!=-1)
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


trans_system_path_node* newtrans_system_path_node(int x,int y,int z,int s)
{
    trans_system_path_node* nd = new trans_system_path_node;
    nd->coord.x=x;
    nd->coord.y=y;
    nd->coord.z=z;
    nd->f = 0.0;
    nd->g = 0.0;
    nd->h = 0.0;
    nd->visit = 0;
    nd->state = s;
    nd->par=NULL;
    return nd;
}

prod_graph_node* new_prod_graph_node(int p,int s,int x,int y,int z)
{
    //p-system state s- automata state
    prod_graph_node* nd = new prod_graph_node;
    nd->coord.x=x;
    nd->coord.y=y;
    nd->coord.z=z;
    nd->f = 0.0;
    nd->g = 0.0;
    nd->h = 0.0;
    nd->visit = 0;
    nd->pstate = p;
    nd->state = s;
    nd->par=NULL;
    return nd;
}

void printpt(point a)
{
    cout<<a.x<<","<<a.y<<","<<a.z<<" ";
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
    float disval[] = {abs(a.x-b.x),abs(a.y-b.y),abs(a.z-b.z)};
    for(int i=0;i<3;i++)
    {
      for(int j=0;j<3-i-1;j++)
      {
        if(disval[j]>disval[j+1])
        {
          float temp = disval[j];
          disval[j] = disval[j+1];
          disval[j+1]=temp;
        }
      }
    } 
    float val = disval[0]*1.8+(disval[1]-disval[0])*1.5+(disval[2]-disval[1])*1;  
    return val;
}


bool is_neighbour(point pt,point nbh)
{    
    
    if(abs(nbh.x-pt.x)<=1 && abs(nbh.y-pt.y)<=1 && abs(nbh.z-pt.z)<=1)
        return 1;
    
    return 0;
}

float cal_approx_cost(point a,point b)
{    
    float disval[] = {abs(a.x-b.x),abs(a.y-b.y),abs(a.z-b.z)};
    for(int i=0;i<3;i++)
    {
      for(int j=0;j<3-i-1;j++)
      {
        if(disval[j]>disval[j+1])
        {
          float temp = disval[j];
          disval[j] = disval[j+1];
          disval[j+1]=temp;
        }
      }
    } 
    float val = disval[0]*1.8+(disval[1]-disval[0])*1.5+(disval[2]-disval[1])*1;  
    return val;
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