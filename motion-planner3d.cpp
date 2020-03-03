#include<bits/stdc++.h>
#include"planner_variables3d.h"
using namespace std;
using namespace planner_info;

/***djkistra with qtot system states with loops,until,next,few corrections,integer hash improvements,consider node in productgraph for book-keeping***/
//correct the adj updations

unordered_map<unsigned long long int,vector<prod_graph_node> > node_neighbour;


pair<int,int> compute_suffix_cycles();

pair<int,int> create_automata_trans_table(char ltl_query[])
{

    string s="";
    const int BUFSIZE=128;
    char buf[BUFSIZE];
    char buff[1000];
    strcpy(buff,ltl_query);
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
    unordered_map<string,int> mp;
    unordered_map<string,int> state;
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

void initializegrid(char *grid_info_file)
{
    string s;
    int i,j,k,no_of_trans,literal,num_obstacles;
    /////////////////////////////////////
    char buff[1000];
    FILE *ft;
    ft = fopen("query.dat", "r");
    /*reading LTL query*/
    fgets(buff, 1000, (FILE*)ft);
    pair<int,int> automata_info;
    /**reading the automata***/
    automata_info = create_automata_trans_table(buff);
    no_of_trans = automata_info.second;
    qtot = automata_info.first;

    /**vector for storing the automata transition table**/
    trans = vector< vector< vector< vector<int> > > >(qtot,vector< vector< vector<int> > >(qtot) );
    /**storing automata transitions on conjunction of negative literals**/
    negtrans = vector< vector< vector< vector<int> > > >(qtot,vector< vector< vector<int> > >(qtot) );
    neg_trans_to_neighbour = vector<int>(qtot);
    ifstream ifile;
    ifile.open("trans.dat");
    cout<<"Number of automata_states="<<qtot<<"\n";
    cout<<"Number of automata_transitions="<<no_of_trans<<"\n";
    /**reading automata transitions**/
    vector<int> automata_states;
    vector<int> transition_condn;
    
    for(i=0;i<no_of_trans;i++)
    {
        transition_condn = vector<int>(0);
        automata_states = vector<int>(2);
        ifile>>automata_states[0];
        ifile>>automata_states[1];
        int trans_condn_len=0,neg_literals=0,strict_neg_literals=0;
        ifile>>trans_condn_len;
        for(j=0;j<trans_condn_len;j++)
        {
            ifile>>literal;
            transition_condn.push_back(literal);
            if(literal <= 0)
                neg_literals++;
            if(literal < 0)
                strict_neg_literals++;
        }
        /**storing transition condition  in automata**/
        trans[automata_states[0]][automata_states[1]].push_back(transition_condn);
        
        /**if transition condition is conjunction of negative literals**/ 
        if(neg_literals==trans_condn_len)
        {
            if(automata_states[0]!=automata_states[1] && strict_neg_literals==trans_condn_len)
            neg_trans_to_neighbour[automata_states[0]]=1;
            negtrans[automata_states[0]][automata_states[1]].push_back(transition_condn);
        }
    }
    // cout<<"here1";
    
    /////////////////////////////////////
    ifstream grid_file;
    grid_file.open(grid_info_file);
    int x,y,z;
    grid_file>>lx;
    grid_file>>ly;
    grid_file>>lz;
    sz = max(max(lx,ly),lz);
    grid_file>>num_obstacles;
    grid_sz = sz;
    sz = max(sz,qtot);
    /*storing obstacle coordinates*/
    for(i=0;i<num_obstacles;i++)
    {
        grid_file>>x;
        grid_file>>y;
        grid_file>>z;
        grid[x*sz*sz+y*sz+z]=-1;
    }
    
    int num_pos_system_states;
    vector<int> grid_state;
    grid_file>>num_pos_system_states;
    // cout<<"here:"<<num_pos_system_states;
    
    /** reading coordinates of states and proposition true at it **/
    for(i=0;i<num_pos_system_states;i++)
    {
        grid_state = vector<int>(3);
        grid_file>>grid_state[0];
        grid_file>>grid_state[1];
        grid_file>>grid_state[2];
        grid_file>>literal;
        s = conv_vec_to_string(grid_state);
        /**pushing the proposition true at state s**/
        prop[s].push_back(literal);
        /**storing cells associated with a proposition**/
        prop_sys_states[literal].push_back(grid_state);
        
        
        /**storing the list of states with a proposition true at it**/
        if(prop[s].size()==1)
            pos_system_state.push_back(grid_state);
        // cout<<"prop:"<<i<<"\n";
        
    }

    grid_file.close();
    
    adj = vector<unordered_map<int,float> >(100);
    updated = vector<unordered_map<int,int> >(100);
}




void calsystates()
{

    int i,j,k,q,z,cnt=0,tmp_cnt=0;
    vector<int> v;
    string s;
    //TS states satifying an incoming transition to an automaton state
    qsystate = vector< vector<int> >(qtot);
    
    //TS states satifying an automaton transition
    qsytrans = vector< vector< vector<int> > >(qtot,vector< vector<int> >(qtot) );

    for(int nstate=0;nstate<qtot;nstate++)
    {
        for(j=0;j<qtot;j++)
        {
            for(z=0;z<trans[nstate][j].size();z++)
            {
                vector<int> transition_req = trans[nstate][j][z];
                if(transition_req.size()==0)
                    continue;

                for(i=0;i<pos_system_state.size();i++)
                {

                    point pt;
                    //transition
                    pt.x =  pos_system_state[i][0]; 
                    pt.y =  pos_system_state[i][1];
                    pt.z =  pos_system_state[i][2];

                    v = vector<int>(0);
                    v.push_back(pt.x);
                    v.push_back(pt.y);
                    v.push_back(pt.z);

                    vector<int> prop_cur;
                    s = conv_vec_to_string(v);
                    if(prop.find(s)==prop.end())
                    {
                        prop_cur.push_back(0);
                    }
                    else
                    {
                        prop_cur = prop[s];
                    }


                    int pos=0;
                    cnt=0;
                    for(int k=0;k<transition_req.size();k++)
                    {
                        if(transition_req[k]>=0)
                            for(int l=0;l<prop_cur.size();l++)
                            {
                                if(transition_req[k]==prop_cur[l])
                                    {cnt++; pos=1; break;}
                            }
                            else
                            {
                                int tmp=1;
                                for(int l=0;l<prop_cur.size();l++)
                                {
                                    if(-1*transition_req[k]==prop_cur[l])
                                        {tmp=0; break;}
                                }
                                cnt+=tmp;
                                if(tmp==0)
                                    break;
                            }
                    }
                    //if(pos && cnt==transition_req.size() ) added
                    //system states from self loop in automata included
                    if(pos && cnt==transition_req.size())
                    {
                        qsystate[j].push_back(i);
                        qsytrans[nstate][j].push_back(i);
                    }

                }
            }
        }
    }
}


float cal_actual_distance_bw_points(point srcnode,point dest_node)
{
    priority_queue<trans_system_node,vector<trans_system_node>,trans_system_node_cost_comparator> qopen;

    unordered_map<long long int,float> node_dist_from_src;
    unordered_map<long long int,int> init;
    unordered_map<long long int,int> visited;

    trans_system_node start_node;
    start_node.coord  = srcnode;
    start_node.g = 0;
    start_node.h = cal_heuristic_cost(start_node.coord,dest_node);
    start_node.f = start_node.g+start_node.h;

    qopen.push(start_node);
    long long int coord_value = key(srcnode);
    init[coord_value]=1;
    node_dist_from_src[coord_value]=0;

    while(!qopen.empty())
    {
        trans_system_node nd = qopen.top();
        qopen.pop();
        point pt = nd.coord;
        coord_value = key(pt);

        if(equal_point(nd.coord,dest_node))
        return nd.f;

        if(visited[coord_value])
            continue;
        else
            visited[coord_value]=1;

        point nbh;
        
        for(int i=0;i<ndir;i++)
        {
            nbh.x =  pt.x+nb[i][0];
            nbh.y =  pt.y+nb[i][1];
            nbh.z =  pt.z+nb[i][2];
            if(!valid(nbh.x,nbh.y,nbh.z) && (nbh.z!=dest_node.z || nbh.y!=dest_node.y || nbh.x!=dest_node.x))
                continue;

            if(abs(nb[i][0])+abs(nb[i][1])+abs(nb[i][2])>1 && (!valid(pt.x+nb[i][0],pt.y,pt.z) || !valid(pt.x,pt.y+nb[i][1],pt.z) || !valid(pt.x,pt.y,pt.z+nb[i][2])))
                continue;

            
            trans_system_node tmp;
            coord_value = key(nbh);
            tmp.g = nd.g + cal_heuristic_cost(pt,nbh);
            if(!init[coord_value])
            {
                tmp.coord = nbh;
                tmp.h = cal_heuristic_cost(nbh,dest_node);
                tmp.f = tmp.g+tmp.h;
                qopen.push(tmp);
                init[coord_value]=1;
                node_dist_from_src[coord_value]=tmp.g;
            }
            else
            {
                if(tmp.g < node_dist_from_src[coord_value])
                {
                    tmp.coord = nbh;
                    tmp.h = cal_heuristic_cost(nbh,dest_node);
                    tmp.f = tmp.g+tmp.h;
                    qopen.push(tmp);
                    node_dist_from_src[coord_value]=tmp.g;
                }

            }
        }
        

    }
    //pair<float,float> no_path{-1,-1};
    return -1;
}

int print_astar_computed_path(trans_system_path_node* nd)
{
    if(nd!=NULL)
    {
        int path_len = print_astar_computed_path(nd->par);
        cur_path.push_back(nd->coord);
        return path_len+1;
    }
    return 0;
}

/**adding obstacles to the grid**/
void copy_map(unordered_map<long long int,int> &obs,int v)
{
    unordered_map<long long int,int>::iterator it;
    for(it = obs.begin();it!=obs.end();it++)
    {
        if(it->second)
        grid[it->first]=v;
    }
}


int astar_path_src_to_dest(prod_graph_node* src,prod_graph_node* dest)
{
    /////////calculate path from src to dest/////////////////////
    unordered_map<long long int,trans_system_path_node* > init_vertex;
    unordered_map<long long int,int > vis;
    int chng=0,final_motion_state,path_len;

    trans_system_path_node* st = newtrans_system_path_node(src->coord.x,src->coord.y,src->coord.z,src->state);

    priority_queue<trans_system_path_node*,vector<trans_system_path_node*>,trans_system_path_node_cost_comparator> qopen;
    vector<trans_system_path_node*> closed;
    qopen.push(st);

    vector<int> src_vertex{src->coord.x,src->coord.y,src->coord.z,src->state};
    
    init_vertex[key(src_vertex)]=st;
    int v,jt,k,it,automaton_state = src->state,undesired_prop;   
    unordered_map<long long int,int> obs1;
    
    //marking obstacles due to self-trans-on conjunction_of_negative literals on current_automaton_state
    if(negtrans[automaton_state][automaton_state].size() > 0)
    {
        vector<int> transition_req = negtrans[automaton_state][automaton_state][0];
        for(it=0;it<transition_req.size();it++)
        {
            undesired_prop = transition_req[it]*-1;
            for(v=0;v<prop_sys_states[undesired_prop].size();v++)
            {
                obs1[key(prop_sys_states[undesired_prop][v])]=-1;
            }
        }

        for(int ind=1;ind<negtrans[automaton_state][automaton_state].size();ind++)
        {
            unordered_map<long long int,int> obs2;
            vector<int> transition_req = negtrans[automaton_state][automaton_state][ind];
            for(it=0;it<transition_req.size();it++)
            {
                undesired_prop = transition_req[it]*-1;
                /**marking all states where undesired_prop is true as obstacles**/
                for(v=0;v<prop_sys_states[undesired_prop].size();v++)
                {
                    if(obs1[key(prop_sys_states[undesired_prop][v])])
                        obs2[key(prop_sys_states[undesired_prop][v])]=-1;
                }
            }
            obs1 = obs2;
        }
        /**unmarking all states satisfying a conjunction containing positive literal marked as obstacles**/
        for(int z=0;z<qsytrans[automaton_state][automaton_state].size();z++)
        {
            k = qsytrans[automaton_state][automaton_state][z];
            obs1[key(pos_system_state[k])]=0;
        }

        copy_map(obs1,-1);    
    }

    while(!qopen.empty())
    {
        trans_system_path_node* cur_nd = qopen.top();
        qopen.pop();
        point cur_pt = cur_nd->coord;

        int current_automaton_state = cur_nd->state;
        vector<int> cur_node_info{cur_pt.x,cur_pt.y,cur_pt.z,current_automaton_state};

        if(equal_point(cur_nd->coord,dest->coord))
        {
            cur_nd->state = dest->state;
            path_len = print_astar_computed_path(cur_nd);
            break;
        }

        long long int cur_node_k = key(cur_node_info);
        
        if(vis.find(cur_node_k)!=vis.end())
            continue;
        else
            vis[cur_node_k]=1;

        vector<int> neigh_state;
        /** relaxing node cur_node and it's neighbours**/
        point nbh;
        
        
        for(int i=0;i<ndir;i++)
        {
            
            nbh.x =  cur_pt.y+nb[i][0];
            nbh.y =  cur_pt.x+nb[i][1];
            nbh.z =  cur_pt.z+nb[i][2];

            if(!valid(nbh.x,nbh.y,nbh.z) && ((nbh.z)!=dest->coord.z || (nbh.y)!=dest->coord.y || (nbh.x!=dest->coord.x) ))
                continue;
            if(abs(nb[i][0])+abs(nb[i][1])+abs(nb[i][2])>1 && (!valid(cur_pt.x+nb[i][0],cur_pt.y,cur_pt.z) || !valid(cur_pt.x,cur_pt.y+nb[i][1],cur_pt.z) || !valid(cur_pt.x,cur_pt.y,cur_pt.z+nb[i][2])))
                continue;
            
            neigh_state = vector<int>(0);
            neigh_state.push_back(nbh.x);
            neigh_state.push_back(nbh.y);
            neigh_state.push_back(nbh.z);
            neigh_state.push_back(current_automaton_state);


            trans_system_path_node* oldtmp;
            float dis_val = cal_heuristic_cost(cur_pt,nbh);
            if(init_vertex.find(key(neigh_state))==init_vertex.end())
            {   
                trans_system_path_node* neigh_node = newtrans_system_path_node(nbh.x,nbh.y,nbh.z,current_automaton_state);
                neigh_node->g = cur_nd->g+dis_val;
                neigh_node->h = cal_heuristic_cost(nbh,dest->coord);
                neigh_node->f = neigh_node->g+neigh_node->h;
                neigh_node->par = cur_nd;
                init_vertex[key(neigh_state)]=neigh_node;
                qopen.push(neigh_node);
            }
            else
            {
                oldtmp = init_vertex[key(neigh_state)];
                if(oldtmp->g > (cur_nd->g+dis_val) )
                {

                    trans_system_path_node* neigh_node = newtrans_system_path_node(nbh.x,nbh.y,nbh.z,current_automaton_state);
                    neigh_node->g = cur_nd->g+dis_val;
                    neigh_node->h = oldtmp->h;
                    neigh_node->f = neigh_node->g+neigh_node->h;
                    neigh_node->par = cur_nd;
                    init_vertex[key(neigh_state)]=neigh_node;
                    qopen.push(neigh_node);
                }

            }
            
        }
    closed.push_back(cur_nd);
    }

    for(it=0;it<closed.size();it++)
    {
        delete(closed[it]);
    }
 
    copy_map(obs1,0);
    
    return path_len;
}



int printpath(prod_graph_node* nd)
{
    
    if(nd!=NULL)
    {
        int path_len = printpath(nd->par);
        if(nd->par!=NULL)
        {   
            //cout<<"par motion state::="<<nd->par->motion_state<<"\n";
            path_len+=astar_path_src_to_dest(nd->par,nd); 
            //cout<<"node motion state::="<<nd->motion_state<<"\n";
            return path_len;
        }
    }
    return 0;
}


void modifyedges(prod_graph_node* nd,int &total_updates)
{
    
    if(nd!=NULL)
    {
        modifyedges(nd->par,total_updates);

        if(nd->par!=NULL)
        {   
            int par_node=nd->par->pstate;
            int cur_node=nd->pstate;
            
            point par_coord,cur_coord;
            par_coord = nd->par->coord;
            cur_coord = nd->coord;
            

            if(!updated[par_node][cur_node])
            {
                int v,k,it,par_automaton_state = nd->par->state,undesired_prop,grid_state;
                float tmp_dis;

                unordered_map<long long int,int> obs1;
            
                if(negtrans[par_automaton_state][par_automaton_state].size() > 0)
                {
                    //adding obstacles on negative self loops
                    vector<int> transition_req = negtrans[par_automaton_state][par_automaton_state][0];
                    for(it=0;it<transition_req.size();it++)
                    {
                        undesired_prop = transition_req[it]*-1;
                        for(v=0;v<prop_sys_states[undesired_prop].size();v++)
                        {
                            obs1[key(prop_sys_states[undesired_prop][v])]=-1;
                        }
                    }
                    for(int ind=1;ind<negtrans[par_automaton_state][par_automaton_state].size();ind++)
                    {
                        unordered_map<long long int,int> obs2;
                        vector<int> transition_req = negtrans[par_automaton_state][par_automaton_state][ind];
                        for(it=0;it<transition_req.size();it++)
                        {
                            undesired_prop = transition_req[it]*-1;
                            for(v=0;v<prop_sys_states[undesired_prop].size();v++)
                            {
                                if(obs1[key(prop_sys_states[undesired_prop][v])])
                                obs2[key(prop_sys_states[undesired_prop][v])]=-1;
                            }
                        }
                        obs1 = obs2;
                    }
                    for(int z=0;z<qsytrans[par_automaton_state][par_automaton_state].size();z++)
                    {
                        grid_state = qsytrans[par_automaton_state][par_automaton_state][z];
                        obs1[key(pos_system_state[grid_state])]=0;
                    }

                    copy_map(obs1,-1);    
                    
                   
                    tmp_dis = cal_actual_distance_bw_points(par_coord,cur_coord);
                    

                    if(!updated[par_node][cur_node] || tmp_dis < adj[par_node][cur_node])
                    {
                        updated[par_node][cur_node] = 1;
                        adj[par_node][cur_node] = tmp_dis;
                    }
                    
                    copy_map(obs1,0);
                }
                else
                {
                    
                    pair<float,float> info;
                    tmp_dis = cal_actual_distance_bw_points(par_coord,cur_coord);
                    adj[par_node][cur_node] = tmp_dis;
                    updated[par_node][cur_node] = 1;
                }
                total_updates++;
            }
            else if(adj[par_node][cur_node]>=2)
                comp_red++;
            if(adj[par_node][cur_node]>=2)
                tot_comp++;
        }

    }

}

/**assign a id/pstate number to each product graph node**/
int assign_systate_no(long long int coord_value)
{
    systate_no[coord_value] = total_system_states;
    total_system_states++;
    
    if(total_system_states>adj.size())
    {
        int cur_size = adj.size();
        adj.resize(2*cur_size);
        updated.resize(2*cur_size);
    }
    return total_system_states-1;
}


void expand_prod_graph_node(priority_queue<prod_graph_node*,vector<prod_graph_node*>,prod_graph_comparator> &qopen,prod_graph_node* cur_nd,unordered_map<long long int,prod_graph_node*> &init_vertex)
{
    int zeroloop=0,negloop=0;
    int current_automaton_state = cur_nd->state;
    point cur_pt=cur_nd->coord;
    //check whether exists a self transition on 1 of current_automaton_state//
    if(trans[current_automaton_state][current_automaton_state].size()==1 && trans[current_automaton_state][current_automaton_state][0].size()==1 && trans[current_automaton_state][current_automaton_state][0][0]==0)
        zeroloop=1;
    //check whether exists a self transition on conjunction of neg literals of current_automaton_state//
    if(negtrans[current_automaton_state][current_automaton_state].size()>0) 
        negloop=1;

    point nbh;
    int cur_node = cur_nd->pstate,neigh_state;
    vector<int> cur_node_info{cur_pt.x,cur_pt.y,cur_pt.z,cur_nd->state};

    //for every outgoing automaton transition to a neighbour automaton state j//
    if(node_neighbour.find(key(cur_node_info))==node_neighbour.end())
    {
	    for(int j=0;j<qtot;j++)
	    {
	        if(trans[current_automaton_state][j].size()==1 && trans[current_automaton_state][j][0].size()==1 && trans[current_automaton_state][j][0][0]==0)
	        {   
	            //if exists a transition from current_automaton_state to j on 1//
	            prod_graph_node neigh_nd;
	          	neigh_nd.coord = nbh;
	            neigh_nd.state = j;
	            node_neighbour[key(cur_node_info)].push_back(neigh_nd);
	        }
	        else
	        {
	            //for all the grid states satisfying the transition from current_automaton_state to j//

                {
    	            for(int z=0;z<qsytrans[current_automaton_state][j].size();z++)
    	            {
                        neigh_state = qsytrans[current_automaton_state][j][z];
                        vector<int> neigh_node_info{pos_system_state[neigh_state][0],pos_system_state[neigh_state][1],pos_system_state[neigh_state][2],j};
                        nbh.x = neigh_node_info[0];
                        nbh.y = neigh_node_info[1];
                        nbh.z = neigh_node_info[2];
                        
                        
                        
                        if(!valid(nbh.x,nbh.y,nbh.z) || (!negloop && !zeroloop && !is_neighbour(cur_pt,nbh)))
                            continue;
                        prod_graph_node neigh_nd;
    		            neigh_nd.coord = nbh;
    		            neigh_nd.state = j;
    		            node_neighbour[key(cur_node_info)].push_back(neigh_nd);

    	            }
                }
	        }
	    }
	    
        
	    for(int j=0;j<qtot;j++)
	    {
	        if(j==current_automaton_state && !neg_trans_to_neighbour[j])
	            continue;
	        //if(j==nstate)
	           // continue;
	        //for all outgoing automaton transitions on conjunction of negative literals check whether a neighbour of cur_pt in grid satisfies the trans condn//
	        for(int z=0;z<negtrans[current_automaton_state][j].size();z++)
	        {
	            vector<int> transition_req = negtrans[current_automaton_state][j][z];

	            for(int i=0;i<ndir;i++)
	            {
                    
	                nbh.x =  cur_pt.x+nb[i][0];
                    nbh.y =  cur_pt.y+nb[i][1];
                    nbh.z =  cur_pt.z+nb[i][2];
                    if(!valid(nbh.x,nbh.y,nbh.z))
	                    continue;
	                
	                vector<int> neigh_grid_cell{nbh.x,nbh.y,nbh.z};
	                string nbh_s = conv_vec_to_string(neigh_grid_cell);
	                
	                vector<int> prop_nbh;
	                
	                //finding the propositions true at neighbour cell nbh in grid//
	                if(prop.find(nbh_s)==prop.end())
	                {
	                    prop_nbh.push_back(0);
	                }
	                else
	                {
	                    prop_nbh = prop[nbh_s];
	                    prop_nbh.push_back(0);
	                }

	                int trans_literals_satisfied=0;
	                for(int k=0;k<transition_req.size();k++)
	                {
	                    if(transition_req[k]>=0)
	                        for(int l=0;l<prop_nbh.size();l++)
	                        {
	                            if(transition_req[k]==prop_nbh[l])
	                            {
	                            trans_literals_satisfied++; break;
	                            }
	                        }
	                    else
	                    {
	                        int satisfied_neg_literal=1;
	                        for(int l=0;l<prop_nbh.size();l++)
	                        {
	                            if(-1*transition_req[k]==prop_nbh[l])
	                            {
	                                satisfied_neg_literal=0; break;
	                            }
	                        }
	                        trans_literals_satisfied+=satisfied_neg_literal;
	                    }
	                }
	                if(trans_literals_satisfied==transition_req.size())
	                {
	                    prod_graph_node neigh_nd;
			            neigh_nd.coord = nbh;
			            neigh_nd.state = j;
			            node_neighbour[key(cur_node_info)].push_back(neigh_nd);
	                }
	            }
	        }
	    }
	}
	

    vector<prod_graph_node> neigh_set;
    neigh_set = node_neighbour[key(cur_node_info)];
    for(int i=0;i<neigh_set.size();i++)
    {
        vector<int> neigh_nd{neigh_set[i].coord.x,neigh_set[i].coord.y,neigh_set[i].coord.z,neigh_set[i].state};
        unsigned long long int coord_value = key(neigh_nd); 
        if(systate_no.find(coord_value)==systate_no.end())
            assign_systate_no(coord_value);
        point nbh = neigh_set[i].coord;
        long long int neigh_state = systate_no[coord_value];
        if(!updated[cur_node][neigh_state])
            adj[cur_node][neigh_state] = cal_approx_cost(cur_pt,nbh);
        if(adj[cur_node][neigh_state]==-1)
            continue;
        
        float tval = adj[cur_node][neigh_state];

        if(init_vertex.find(key(neigh_nd))==init_vertex.end())
        {
            prod_graph_node* tmp = new_prod_graph_node(neigh_state,neigh_set[i].state,nbh.x,nbh.y,nbh.z);
            tmp->g = cur_nd->g+tval;
            tmp->h = 0;
            tmp->f = tmp->g+tmp->h;
            tmp->par = cur_nd;
            init_vertex[key(neigh_nd)]=tmp;
            qopen.push(tmp);
        }
        else
        {
            prod_graph_node* oldtmp = init_vertex[key(neigh_nd)];
            
            if(oldtmp->g > (cur_nd->g+tval) )
            {
                prod_graph_node* tmp = new_prod_graph_node(neigh_state,neigh_set[i].state,nbh.x,nbh.y,nbh.z);
                tmp->g = cur_nd->g+tval;
                tmp->h = 0;
                tmp->f = tmp->g+tmp->h;
                tmp->par = cur_nd;
                init_vertex[key(neigh_nd)]=tmp;
                qopen.push(tmp);
            }

        }

    }
}


int find_path_in_prod_graph(prod_graph_node* srcnode,prod_graph_node* dest_node,int degen)
{
    int iter;
    //iterate till we get a path where no edge cost has been modified//
    for(iter=1; iter>0 ;iter++)
    {
        unordered_map<long long int ,prod_graph_node*> init_vertex;
        unordered_map<long long int ,int> vis;
        int destreach=0;
        prod_graph_node* st = new_prod_graph_node(srcnode->pstate,srcnode->state,srcnode->coord.x,srcnode->coord.y,srcnode->coord.z);
        priority_queue<prod_graph_node*,vector<prod_graph_node*>,prod_graph_comparator> qopen;
        point cur_pt = st->coord;
        vector<int> cur_node_info{cur_pt.x,cur_pt.y,cur_pt.z,st->state};
        long long int coord_value = key(cur_node_info); 
        
        if(systate_no.find(coord_value)==systate_no.end())
        assign_systate_no(coord_value);
        
        st->pstate = systate_no[coord_value];
        int path_valid=0;
        vector<prod_graph_node*> closed;
        qopen.push(st);
        //if zero length paths are allowed//
        if(degen)
        {
        init_vertex[key(cur_node_info)]=st;
        path_valid=1;
        }

        //path_valid denotes path length is greater than 0 or path length of 0 is allowed in case of finding prefix only//
        while(!qopen.empty())
        {
            prod_graph_node* cur_nd = qopen.top();
            qopen.pop();
            point cur_pt = cur_nd->coord;
            int current_automaton_state = cur_nd->state;
            vector<int> cur_node_info{cur_pt.x,cur_pt.y,cur_pt.z,cur_nd->state};
            // printvec(cur_node_info);
            // cout<<"\n";
            long long int  s;
            s = key(cur_node_info);

            if(path_valid)
            {
            if(vis.find(s)!=vis.end())
                continue;
            else
                vis[s]=1;
            }

            if(path_valid && equal_point(cur_nd->coord,dest_node->coord) && cur_nd->state==dest_node->state )
            {
                int total_updates=0;
                path = vector<int>(0);
                
                modifyedges(cur_nd,total_updates);
                
                destreach=1;
                if(iter>1 && total_updates==0)
                {
                    cur_path = vector<point>(0);
                    iter=-5; 
                    int len = printpath(cur_nd);
                    if(!degen)
                    {
                        stored_suffix_paths.push_back(cur_path);
                        cur_nd->pstate = stored_suffix_paths.size()-1;
                        loopf.push(cur_nd);
                    }
                    else 
                        prefix_length = cur_nd->f;
                }
                break;
            }

            if(!path_valid)
                path_valid=1;

            expand_prod_graph_node(qopen,cur_nd,init_vertex);
            closed.push_back(cur_nd);
        
        }
        for(int i=0;i<closed.size();i++)
            delete(closed[i]);
        if(!destreach)
        break;
    }   
    return iter;
}

int main(int args,char **argv)
{
    
    
    clock_t t;
    t = clock();
    initializegrid(argv[1]);
    
    int cntstate=0;
    calsystates();

    pair<int,int> dest_node;
    cout<<"compute_suffix_cycles\n";

    compute_suffix_cycles();
    cout<<"total suffix_path="<<loopf.size()<<"\n";
    //loopf all suffix cycle starting nodes
    
    while(!loopf.empty())
    {
    prod_graph_node* finstate = loopf.top();
    prod_graph_node* st = new_prod_graph_node(0,0,0,0,0);
    loopf.pop();
    int iter,cnt=0;
    iter = find_path_in_prod_graph(st,finstate,1);
    

    if(iter < 0)
    {
        prefix_path = cur_path;
        suffix_path = stored_suffix_paths[finstate->pstate];
        suf_len = finstate->f;
        break;
    }
    }
    
    cout<<"prefix_length="<<prefix_length<<"\n";
    cout<<"suffix_length="<<suf_len<<"\n";
    cout<<"\nprefix\n";
    vector<int> velocity_primitive;
    
    for(int i=0;i<prefix_path.size();i++)
    {
        printpt(prefix_path[i]);
    }

    cout<<"\nsuffix\n";
    for(int i=0;i<suffix_path.size();i++)
    {
        printpt(suffix_path[i]);
    }

    cout<<"\ntotal_product_graph_nodes="<<total_system_states;
    cout<<"\ncomp_red="<<comp_red;
    cout<<"\ntot_comp="<<tot_comp;
    t = clock() - t;
    double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds
 
    printf("\nprogram took %f seconds to execute \n", time_taken);
    cout<<"\n.......................................................................\n";
}

pair<int,int> compute_suffix_cycles()
{
    int iter,i,suflen=INT_MAX,nsuf=0;
    pair<int,int> dest_node;
    

    /*******for every final state of automaton***************/
    for(int fs=0;fs<dest.size();fs++)
    {
        vector<int> vis_sy_state(pos_system_state.size());
        /*******for every sytem state associated with the given automaton final state dest[fs]*********/
        for(i=0;i<qsystate[dest[fs]].size();i++)
        {
            int systate = qsystate[dest[fs]][i];
            if(vis_sy_state[systate]==1)
                continue;
            else
                vis_sy_state[systate]=1;
            prod_graph_node* st = new_prod_graph_node(systate,dest[fs],pos_system_state[systate][0],pos_system_state[systate][1],pos_system_state[systate][2]);
            nsuf++;
            find_path_in_prod_graph(st,st,0);

        }
    }

    cout<<"attempted="<<nsuf<<"\n";

    return dest_node;
}