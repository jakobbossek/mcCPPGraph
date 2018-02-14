#include "../graph.h"
#include "../UnionFind.h"
#include "../debug.h"
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <sys/types.h>
#include <dirent.h>
#include <string.h>

using namespace std;

// source: http://www.martinbroadhurst.com/list-the-files-in-a-directory-in-c.html
std::vector<std::string> read_directory(const std::string& name) {
  std::vector<std::string> files;
  DIR* dirp = opendir(name.c_str());
  struct dirent * dp;
  while ((dp = readdir(dirp)) != NULL) {
    files.push_back(dp->d_name);
  }
  closedir(dirp);
  return files;
}

double sum(std::vector<double> x) {
  double sum = 0.0;
  for (auto& element : x) {
    sum += element;
  }
  return sum;
}

double mean(std::vector<double> x) {
  return sum(x) / x.size();
}

int main() {

  // setup
  int nruns = 100;
  vector<string> files = read_directory("../tests/instances/");
  std::string pathToDir = "../tests/instances/";

  fstream logfile("cpp.csv", ios::out | ios::trunc);
  logfile << "prob time.import runtime.mean runtime.total time.total" << endl;

  for (auto pathToFile : files) {
    if (pathToFile == "." || pathToFile == ".." || pathToFile == ".DS_Store")
      continue;

    cout << pathToFile << endl;

    // small benchmark to check how long it takes to import and measure stuff
    std::vector<double> runtimes(nruns);

    std::string pathToFile2 (pathToFile);
    std::string pathToFile3 = pathToDir + pathToFile2;

    cout << pathToFile3 << endl;

    double timeStart = clock();
    Graph g = Graph::importFromGrapheratorFile(pathToFile3);
    double timeImport = ((double)(clock() - timeStart) / (double)CLOCKS_PER_SEC);
    g.print(false);

    for (int i = 0; i < nruns; ++i) {
      cout << ".";
      if (i == nruns - 1)
        cout << endl;
      timeStart = clock();
      Graph mst = g.getMSTKruskal(0.5);
      runtimes[i] = ((double)(clock() - timeStart) / (double)CLOCKS_PER_SEC);
    }
    cout << "===============" << endl << endl;

    double runtimeMean = mean(runtimes);
    double runtimeTotal = sum(runtimes);
    double timeTotal = timeImport + runtimeTotal;

    // print stats
    cout << "Performed " << nruns << " runs on instance " << pathToFile << endl;
    cout << "Time for import     : " << timeImport << endl;
    cout << "Mean running time   : " << runtimeMean << endl;
    cout << "Overall running time: " << runtimeTotal << endl;
    cout << "Overall passed time : " << timeTotal << endl;

    logfile << pathToFile << " " << timeImport << " " << runtimeMean << " " << runtimeTotal << " " << timeTotal << endl;
  } // loop over instances

  // cleanup
  logfile.close();
  exit(0);
}
