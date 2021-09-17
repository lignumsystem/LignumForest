#include <ForestTrees.h>
#include <cmath>

//Given the  tree_file parse it and  return file names  that match and
//exist in the current directory, associate the given probability with
//the file name
vector<pair<string,double> > FindTreeFiles(const string& tree_file,double tree_age)
{
  QStringList filters;
  vector<double> probabilities;
  ifstream file(tree_file.c_str());
  
  if (!file){
    cout << "FindTreeFiles Cannot open file " << tree_file.c_str() <<endl;
    exit(-1);
  }
  
  string line;
  getline(file,line);
  //Eat the comment lines beginning with '#'
  while (line[0]=='#'){
    getline(file,line);
  }
  //Now this is the first tree
  while (!file.eof()){
    istringstream s(line);
    //Looks like Macs do not like string class with istringstream!!
    int age=0;
    char file_name[255];
    double value = 1.0;
    s >> age >>  file_name >> value;
    if (!s.bad()){
      //Collect the given name matching tree age into a QDir filter
      if (tree_age == age){
	filters << QString(file_name);
	probabilities.push_back(value);
      }
    }
    else{
      cout << "FindTreeFiles Failed reading line " << tree_file.c_str() << " " 
	   << file_name << " " << value <<endl;
      exit(0);
    }
    getline(file,line);    
  }


  QDir dir;
  //Set the directory filter
  QStringList ls1;
  //If not filters do not set such filter (all files will be listed)
  if (filters.size()){
    dir.setNameFilters(filters);
  //This assumes that  the existing files are listed  in the order the
  //appear filter list
    ls1 = dir.entryList(QDir::Files);
  }
  vector<pair<string,double> > ls2;
  //The probability for the target tree
  double self_probability=0.0;
  //Attach probabilities  to files, find  the trees accirding  to tree
  //age
  for (int i = 0; i < ls1.size();i++){
    ls2.push_back(pair<string,double>(ls1[i].toStdString(),probabilities[i]));
    //Update target tree probability, eventually it will hold the highest cumulative probability
    self_probability = probabilities[i];
  }
  //There is a probability for the target tree, which is what is left between the last found tree and 1.0
  if (self_probability && self_probability != 1.0)
    ls2.push_back(pair<string,double>("self",1.0));
  //Check if some of  the the files in the tree file  did not exist in
  //the directory
  if (static_cast<unsigned int>(ls1.size()) !=probabilities.size()){
    cout << "Forest tree files: some tree file(s) cannot be found in the directory!" <<endl;
  }
  return ls2;
}

//Given  a random  number  p from  a  uniform distributed  probability
//[0:1], find the  matching tree to be inserted  into voxel space This
//'<='  actually  returns  "lower  bound",  though  it  is  used  with
//'upper_bound' in InsertFoliageIntoVoxels with light method 4.
bool FindTreeWithProbability::operator()(const double p,const pair<string,double>& e2)const
{
  return p <= e2.second;
}

