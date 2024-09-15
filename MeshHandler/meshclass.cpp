#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <sstream>

int main(int argc, char* argv[]) {
    using namespace std;
    ifstream ifs;
    ifs.open("testmesh",ios::in);

    if(!ifs.is_open()) {
        cout << "Error opening file" << endl;
        return 1;
    }


    string buf;
    map<string,vector<vector<int>>> mesh;
    string key;
    string ind;
    while(getline(ifs,buf))
    {
        istringstream iss(buf);
        if(buf[0] == '$')
        {   
            iss >> ind;
            key = ind.substr(1, ind.size()-1);
            //cout << key<< endl; 
            //cout << ind.size() << endl;
        }
        else
        {
            vector<int> tmp;
            while(iss >> ind)
            {
                tmp.push_back(stoi(ind));
            }
            mesh[key].push_back(tmp);
            //cout << mesh[key][0] << endl; 
        }
    }

    //vector<string>* entity = &mesh["Entities"];
    cout << mesh["Entities"][0][0] << endl;

    return 0;

}