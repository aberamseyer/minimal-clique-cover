#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <set>
#include <time.h>

using namespace std;

int main() {
    srand(static_cast<unsigned int>(time(nullptr)));
    ofstream output;
    output.open("25.txt");
    set<pair<int, int>> a;
    int max = 25;
    for(int i = 0; i < max; i++) {
        for(int j = 0; j < max; j++) {
            if(i != j && rand()%100 < 30) {
                pair<int, int> test(i, j);
                pair<int, int> test2(j, i);
                if(a.find(test) != a.end() || a.find(test2) != a.end())
                    continue;
                a.insert(test);
                output << i << " " << j << endl;
            }
        }
    }

    output.close();
    return 0;
}
