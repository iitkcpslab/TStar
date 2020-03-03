# T* ALgorithm


### **Files :**
**1. query.dat**  : LTL query must be specified in this file.

**2. cfile_rec.dat** : 2D workspace descriptor file.
Input format:
nr - number of rows
nc - number of columns
nobs - no of obstacles
nobs coordinates - 'nobs' lines  mentioning the co-odinates of the obstacles
np - no of coordinates where a certain proposition is true
np coordinates - 'np' lines mentioning the coordinates and the proposition true at it.

Example file :
<pre>
10          -----> nr
10          -----> nc
5           -----> nobs
1 2         -----> 'nobs' lines (ex. these is an obstacle at coord (1,2) in the workspace)
3 4
5 6
9 5
3 9
4           -----> np
5 5 1       -----> 'np' lines (ex. at (5,5) , proposition 1 is true)
6 6 2
2 2 3
3 3 4
</pre>

**3. cfile_rec3d.dat** : 3D workspace descriptor file. Same as 2D, we just adding z axis to input and co-ordinates.
Input format:
nx - x axis length 
ny - y axis length 
nz - z axis length 
nobs - no of obstacles followed by coordinates of obstacles 
np - no of coordinates where a certain proposition is true followed by coordinates and the proposition true at it. 
Example file: 
<pre>
10
10
10
5
1 2 1
3 4 1
5 6 2
9 5 3
3 9 5
4
5 5 5 1
6 6 6 2
2 2 2 3
3 3 3 4
</pre>

**4. motion-planner-final.cpp** : T&ast; algorithm for 2D worksapce

**5. motion-planner3d.cpp** : T&ast; algorithm for 3D worksapce

**6. planner_variables.h** : Helper file for 2D workspace

**7. planner_variables3d.h** : Helper file for 3D workspace

**8. ltl2tgba** : LTL to Buchi Automaton converter binary file. 
Install Spot-2.6 tool for LTL2TGBA converter. Copy the ltl2tgba file from the bin folder of spot installation and copy it to the current folder. Check if the tool is working by running following command in the current folder using the command line:
./ltl2tgba \--spin ' [ ] (<>p && <>q) '
This command should give us the Buchi automata transitions for the give quary

**9. djk-optimal-run.cpp**: Optimal_Run algorithm for 2D workspace. This is the Dijkstras  algorithm based solution to the mentioned problem.

----
### Execute Algorithm

**Compile**:  Compile djk-optimal-run.cpp and motion-planner-final.cpp using g++ command with flag -std=c++11 to generate the binaries djk-optimal-run and motion-planner-final(on Ubantu OS) . 
For example -
g++ -std=c++11 djk-optimal-run.cpp -o djk-optimal-run
g++ -std=c++11 motion-planner-final.cpp -o motion-planner-final

**Execute**:  Execute the generated binaries using the following commands
./djk-optimal-run cfile_rec.dat
./motion-planner-final cfile_rec.dat

Similarly, for 3D workspace, execute
./djk-optimal3d cfile_rec3d.dat 
./motion-planner3d cfile_rec3d.dat

------------
