#include<bits/stdc++.h>
#define lli long long int 
using namespace std;
int sz=10,qtot=3,dir=9,grid_sz;
vector<int> dest;
float pre_len,suf_len=numeric_limits<float>::infinity();
clock_t t;

/**********dijkstra with qtot read from the automaton graph-trans.dat*****************************/
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

int nb[][2]={{0,0},{0,-1},{-1,0},{1,0},{0,1},{-1,-1},{1,-1},{-1,1},{1,1}};

void find_suffix_cycles();
unordered_map<long long int,int> grid;
vector< vector< vector< vector<int> > > > trans(qtot,vector< vector< vector<int> > >(qtot) );
unordered_map<lli,vector<int> > prop;
vector<vector<int> > pstate;
vector<vector<int> > qsystate;
vector<struct point> suffix_path,prefix_path,cur_path;
vector<vector<struct point> > stored_suffix_paths;

long long int key(vector<int> &v)
{
    lli val;
    if(v.size()==3)
    val = v[0]*sz*sz+v[1]*sz+v[2];
    else
    val = v[0]*1000+v[1];
    return val;
}

/*string key(vector<int> &vec)
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
}*/


float calh(point a,point b)
{
    float mx = max(abs(a.y-b.y),abs(a.x-b.x));
    float mn = min(abs(a.y-b.y),abs(a.x-b.x));
    return mn*1.5+(mx-mn);
}

void printvec(vector<int> &vec)
{
for(int i=0;i<vec.size();i++)
cout<<vec[i];
}

bool valid(int y,int x)
{
if(x>=0 && x<grid_sz && y>=0 && y<grid_sz && grid[y*sz+x]!=-1)
return 1;
return 0;
}

struct node
{
node* par;
node* child;
point coord;
float f;
bool visit;
int state,pstate;
};


unordered_map<long long int,vector<node> > node_neighbour;
node* newnode(int x,int y,int s)
{
node* nd = new node;
nd->coord.x=x;
nd->coord.y=y;
nd->f = 0.0;
nd->visit = 0;
nd->state = s;
nd->par=NULL;
return nd;
}

struct comp
{
bool operator()(node *a , node *b)
{
return (a->f > b->f);
}
};

priority_queue<node*,vector<node*>,comp> loopf;

pair<int,int> createtrans(char tpc[])
{

  string s="";
  const int BUFSIZE=128;
  char buf[BUFSIZE];
  char buff[1000];
  strcpy(buff,tpc);
  FILE *fp;
  ofstream outfile;
  outfile.open("trans.dat");
  int i,j,k,cnt=0;
  pair<int,int> rval;

  if ((fp = popen(buff, "r")) == NULL) 
  {
    printf("Error opening pipe!\n");
    return rval;
  }
    /*if ((fp = fopen("result3.txt", "r")) == NULL)
    {
        printf("Error opening pipe!\n");
        return rval;
    }*/

  while (fgets(buf, BUFSIZE, fp) != NULL) 
  {
        // Do whatever you want here...
        //printf("%s", buf);
    s = s+buf;
  }

  if(pclose(fp))  
  {
    printf("Command not found or exited with error status\n");
    return rval;
  }
  cout<<"\n";
  map<string,int> mp;
  map<string,int> state;
  vector<string> st;

    //cout<<s;
  string tmp;
  i=0;
    //cout<<s.size();
  while(s[i]!='T' && s[i]!='a')
    i++;
  int stno=0;
  vector<string> table;
  while(i<s.size())
  {
    tmp="";
    while(i<s.size() && s[i]!=':')
    {
      tmp+=s[i];
      i++;
    }
    if(tmp[0]=='}')
      break;
    if(state.find(tmp)==state.end())
    {
      state[tmp] = stno;
      st.push_back(tmp);
      stno++;
    }

     while(i<s.size() && tmp!="fi;")
    {
    tmp = "";
    while(i<s.size() && s[i]!='\t' && s[i]!=' ' && s[i]!='\n')
    {
     if(s[i]!='(' && s[i]!=')')
       tmp+=s[i];
     i++;
    }
    if(tmp!=" " && tmp!="" && tmp!="if" && tmp!="::" && tmp!="->" && tmp!="goto" &&  tmp!="fi;")
    table.push_back(tmp);
    i++;
    }

  }

    //for(i=0;i<table.size();i++)
    //cout<<table[i]<<"\n";

  for(i=0;i<st.size();i++)
  {
    mp[st[i]]=i;
    if(st[i][0]=='a' && st[i][1]=='c' && st[i][2]=='c')
      dest.push_back(i);
  }

  vector<vector<int> > vec;
  vector<int> tmpvec;
  int x,y;
  k=0;
  i=1;

  while(i<table.size())
  {
  if(table[i][0]==':')
  {
    k++;
  }
  else if(table[i]=="||")
  {
   vec.push_back(tmpvec);
   tmpvec = vector<int>(0);
  }
  else if(table[i]=="&&")
  {
  }
  else if(mp.find(table[i])!=mp.end())
  {
    vec.push_back(tmpvec);
    for(x=0;x<vec.size();x++)
    {
      cnt++;
      outfile<<k<<" "<<mp[table[i]]<<" "<<vec[x].size()<<" ";
      for(y=0;y<vec[x].size();y++)
      {
        outfile<<vec[x][y]<<" ";
      }
      outfile<<"\n";
    }
    vec = vector<vector<int> >(0);
    tmpvec = vector<int>(0);
  }
  else if(table[i]=="1")
  {
    tmpvec.push_back(0);
  }
  else if(table[i][0]=='!')
  {
  int v=0;
  for(x=2;x<table[i].size();x++)
    v = v*10+(table[i][x]-'0');
  tmpvec.push_back(-1*v);
  }
  else if(table[i][0]=='p')
  {
   int v=0;
   for(x=1;x<table[i].size();x++)
     v = v*10+(table[i][x]-'0');
   tmpvec.push_back(v);
  }
  else if(table[i]=="skip")
  {cnt++; outfile<<k<<" "<<k<<" 1 0\n";}

  i++;
}
outfile.close();

rval.first = st.size();
rval.second = cnt;
return rval;
}

void initializegrid(char *filename)
{
    vector<int> v;
    vector<int> p;
    lli s;
    int i,j,k,reg,val,nobs;


    char buff[1000];
    FILE *ft;
    ft = fopen("query.dat", "r");
    fgets(buff, 1000, (FILE*)ft);
    fclose(ft);
    pair<int,int> rval;
    rval = createtrans(buff);

    t = clock();

    ifstream infile;
    infile.open(filename);

    infile>>sz;
    infile>>sz;
    infile>>nobs;
    grid_sz = sz;
    if(sz<100)
      sz=100;
    //grid = vector<vector<int> >(sz,vector<int>(sz));
    //cout<<"obstacles:\n";
    for(i=0;i<nobs;i++)
    {
        infile>>j;
        infile>>k;
        grid[j*sz+k]=-1;
        //cout<<j<<","<<k<<" ";
    }
    cout<<"\n";
    infile>>reg;
    for(i=0;i<reg;i++)
    {
        v = vector<int>(0);
        p = vector<int>(0);
        infile>>val;
        p.push_back(val);
        infile>>val;
        p.push_back(val);
        infile>>val;
        s = key(p);
        prop[s].push_back(val);
        if(prop[s].size()==1)
            pstate.push_back(p);
    }
    infile.close();

    
    //cout<<"query: "<<buff<<"\n";
    
    reg = rval.second;
    qtot = rval.first;
    //cout<<qtot<<","<<reg<<"\n";
    trans = vector< vector< vector< vector<int> > > >(qtot,vector< vector< vector<int> > >(qtot) );
    //cout<<"transitions:\n";

    ifstream ifile;
    ifile.open("trans.dat");
    //cout<<reg<<trans.size()<<trans[0].size();
    //ifile>>reg;
    //cout<<reg;
    //cout<<"fgfg"<<"\n";
    
    for(i=0;i<reg;i++)
    {
        v = vector<int>(0);
        p = vector<int>(0);
        ifile>>val;
        p.push_back(val);
        ifile>>val;
        p.push_back(val);
        int len=0;
        ifile>>len;
        for(j=0;j<len;j++)
        {
            ifile>>val;
            v.push_back(val);
        }
        trans[p[0]][p[1]].push_back(v);
        /*cout<<p[0]<<" "<<p[1]<<" ";
        for(j=0;j<v.size();j++)
            cout<<v[j]<<" ";
        cout<<"\n";*/
    }

}


void calfinstates()
{

  int i,j,k,q,z,cnt=0,nstate;
  vector<int> v;
  lli s;
  qsystate = vector<vector<int> >(dest.size());
  //map<string,int> final_prod_states;
  for(int fs=0;fs<dest.size();fs++)
  {
    nstate = dest[fs];
    for(j=0;j<qtot;j++)
    {
      for(z=0;z<trans[j][nstate].size();z++)
      {
        vector<int> req = trans[j][nstate][z];
        if(req.size()==0)
          continue;

        for(i=0;i<pstate.size();i++)
        {
          point pt;
          //transition
          pt.y =  pstate[i][0];
          pt.x =  pstate[i][1];

          v = vector<int>(0);
          v.push_back(pt.y);
          v.push_back(pt.x);

          vector<int> prop_cur;
          s = key(v);
          if(prop.find(s)==prop.end())
          {
            //nxtstate.push_back(nd->state);
            prop_cur.push_back(0);
            cout<<s<<"...\n";
          }
          else
          {
            prop_cur = prop[s];
            //prop_cur.push_back(0);
          }
          
          cnt=0;
          for(int k=0;k<req.size();k++)
          {
            if(req[k]>=0)
              for(int l=0;l<prop_cur.size();l++)
              {
                if(req[k]==prop_cur[l])
                  {cnt++; break;}
              }
              else
              {
                int tmp=1;
                for(int l=0;l<prop_cur.size();l++)
                {
                  if(-1*req[k]==prop_cur[l])
                    {tmp=0; break;}
                }
                cnt+=tmp;
                if(tmp==0)
                  break;
              }
            }

            //system states from self loop in automata included
            if( cnt==req.size() )
            {
              //cout<<nstate<<" "<<j<<" "<<i<<"\n";
              qsystate[fs].push_back(i);
            }

          }
        }
      }
    }
    
}


void printpath(node* nd)
{
  if(nd!=NULL)
  {
  printpath(nd->par);
  //cout<<"**"<<nd->coord.y<<","<<nd->coord.x<<","<<nd->state<<","<<nd->f<<"\n";
  cur_path.push_back(nd->coord);
  }
}

void expand_node(priority_queue<node*,vector<node*>,comp> &qopen,node* nd,unordered_map<long long int,node*> &init_vertex)
{
    point pt = nd->coord;      
    int nstate = nd->state;  
    point nbh;
    int self_flag=1;    
    vector<int> cur_node{nd->coord.y,nd->coord.x,nd->state};

    if(node_neighbour.find(key(cur_node))==node_neighbour.end())
    {
        for(int j=0;j<qtot;j++)
        {
            for(int z=0;z<trans[nstate][j].size();z++)
            {

                vector<int> req = trans[nstate][j][z];
                vector<int> v;
                for(int i=0;i<dir;i++)
                {
                    if(!valid(pt.y+nb[i][0],pt.x+nb[i][1]))
                        continue;
                    if( abs(nb[i][0])+abs(nb[i][1])==2 && (!valid(pt.y+nb[i][0],pt.x) || !valid(pt.y,pt.x+nb[i][1]) ) )
                      continue;


                    //transition
                    nbh.y =  pt.y+nb[i][0];
                    nbh.x =  pt.x+nb[i][1];
                    v = vector<int>(0);
                    v.push_back(nbh.y);
                    v.push_back(nbh.x);
                    
                    vector<int> prop_cur;
                    unsigned long long st_prop = key(v);
                    if(prop.find(st_prop)==prop.end())
                    {
                        prop_cur.push_back(0);
                    }
                    else
                    {
                        prop_cur = prop[st_prop];
                        prop_cur.push_back(0);
                    }

                    int cnt=0;
                    for(int k=0;k<req.size();k++)
                    {
                        if(req[k]>=0)
                            for(int l=0;l<prop_cur.size();l++)
                            {
                                if(req[k]==prop_cur[l])
                                {
                                    cnt++; break;
                                }
                            }
                            else
                            {
                                int tmp=1;
                                for(int l=0;l<prop_cur.size();l++)
                                {
                                    if(-1*req[k]==prop_cur[l])
                                    {
                                        tmp=0; break;
                                    }
                                }
                                cnt+=tmp;
                            }
                    }
                    if(cnt==req.size())
                    {
                        node neigh_node;
                        neigh_node.coord = nbh;
                        neigh_node.state = j;
                        node_neighbour[key(cur_node)].push_back(neigh_node);
                    }
                }
            }
        }
    }


    vector<node> neigh_set = node_neighbour[key(cur_node)];
    
    for(int i=0;i<neigh_set.size();i++)
    {    
        node* oldtmp;
        nbh = neigh_set[i].coord;
        int automaton_state = neigh_set[i].state;
        float tval = calh(pt,nbh);
        vector<int> tv{nbh.y,nbh.x,automaton_state};
        if(init_vertex.find(key(tv))==init_vertex.end())
        {
            node* tmp = newnode(nbh.x,nbh.y,automaton_state);
            tmp->f = nd->f+tval;
            tmp->par = nd;
            init_vertex[key(tv)]=tmp;
            qopen.push(tmp);
        }
        else
        {
            node* oldtmp = init_vertex[key(tv)];
            if(oldtmp->f > (nd->f+tval) )
            {
                node* tmp = newnode(nbh.x,nbh.y,automaton_state);
                tmp->f = nd->f+tval;
                tmp->par = nd;
                // oldtmp->f = nd->f+tval;
                // oldtmp->par = nd;
                init_vertex[key(tv)]=tmp;
                qopen.push(tmp);
            }

        }
    }         
}

int main(int args,char** argv)
{
  //string s(argv[1]);
  //t = clock();
  initializegrid(argv[1]);
  //return 0;
  calfinstates();
  find_suffix_cycles();

  //cout<<"***inside main*****\n";
  while(!loopf.empty())
  {
    node* finstate = loopf.top();
    loopf.pop();
    int iter,cnt=0;
    //cout<<"..destnode:"<<finstate->coord.y<<" "<<finstate->coord.x<<" "<<finstate->state<<"\n";

    unordered_map<lli,node* > init_vertex;
    unordered_map<lli,int > vis;

    node* st = newnode(0,0,0);

    priority_queue<node*,vector<node*>,comp> qopen;
    vector<node*> closed;
    qopen.push(st);

    vector<int> v;
    v.push_back(0);
    v.push_back(0);
    v.push_back(0);
    init_vertex[key(v)]=st;


    while(!qopen.empty())
    {
      node* nd = qopen.top();
      qopen.pop();
      point pt = nd->coord;

      int nstate = nd->state;
      vector<int> v;
      v.push_back(pt.x);
      v.push_back(pt.y);
      v.push_back(nd->state);

      if(nd->state==finstate->state && nd->coord.x==finstate->coord.x && nd->coord.y==finstate->coord.y)
      {
        cur_path = vector<point>(0);
        printpath(nd);
        prefix_path = cur_path;
        pre_len = nd->f;
        suf_len = finstate->f;
        suffix_path = stored_suffix_paths[finstate->pstate];
        break;
      }

      lli ind;
      ind = key(v);
      if(vis.find(ind)!=vis.end())
        continue;
      else
        vis[ind]=1;

      //cout<<"\n"<<pt.y<<","<<pt.x<<","<<nd->state<<" : ";
      expand_node(qopen,nd,init_vertex);
      
        
        closed.push_back(nd);
      }
      if(!qopen.empty())
        break;
  }

  cout<<"length_prefix="<<pre_len<<"\n";
  cout<<"length_suffix="<<suf_len<<"\n";
  cout<<"\nprefix\n";
  for(int i=0;i<prefix_path.size();i++)
      cout<<"**"<<prefix_path[i].y<<","<<prefix_path[i].x<<"\n";

  cout<<"\nsuffix\n";
  for(int i=0;i<suffix_path.size();i++)
      cout<<"**"<<suffix_path[i].y<<","<<suffix_path[i].x<<"\n";

  t = clock() - t;
  double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds
 
  printf("\nprogram took %f seconds to execute \n", time_taken);
  return 0;
}


void find_suffix_cycles()
{
  int iter,i,suflen=INT_MAX;
  pair<int,int> destnode;
    /*******for every final state of automaton***************/
  for(int fs=0;fs<dest.size();fs++)
  {
    /*******for every sytem state associated with the given automaton final state dest[fs]*********/

    vector<int> vis_sy_state(pstate.size());
    for(int it=0;it<qsystate[fs].size();it++)
    {
      //cout<<"starting for point:\n";
      int systate = qsystate[fs][it];
      if(vis_sy_state[systate]==1)
        continue;
      else
        vis_sy_state[systate]=1;
      unordered_map<long long int,node*> init_vertex;
      unordered_map<long long int,int > vis;
      int destreach=0,nstate = dest[fs];
      node* st = newnode(pstate[systate][1],pstate[systate][0],dest[fs]);
      //cout<<st->coord.y<<","<<st->coord.x<<","<<st->state<<"::\n";
      priority_queue<node*,vector<node*>,comp> qopen;
      vector<node*> closed;
      vector<int> v;
      string s;
      point pt = st->coord;
      point nbh;

      ////////////////////**********************************/////////////////////////

      qopen.push(st);
      int occ=0;
      ///////////////////////////*******************************/////////////////////////////////

      //*********************************************************************************************************\\

      while(!qopen.empty())
      {
        node* nd = qopen.top();
        qopen.pop();
        point pt = nd->coord;

        int nstate = nd->state;
        vector<int> v;
        v.push_back(nd->coord.x);
        v.push_back(nd->coord.y);
        v.push_back(nd->state);

        if(occ && nd->state==st->state && nd->coord.x==st->coord.x && nd->coord.y==st->coord.y)
        {
          cur_path = vector<point>(0);
          printpath(nd);
          //cout<<"\npath_length = "<<nd->f<<"\n";
          stored_suffix_paths.push_back(cur_path);
          nd->pstate = stored_suffix_paths.size()-1;
          loopf.push(nd);
          break;
        }
        if(!occ)
          occ=1;

        lli ind;
        ind = key(v);

        if(vis.find(ind)!=vis.end())
          continue;
        else
          vis[ind]=1;

        
        expand_node(qopen,nd,init_vertex);
             

        }
 
      }
    }
}
