#include"aprioriAlgorithm.h"
#include<vector>
int main(){
    std::vector<int> data = {1,2,3,2,1,1,4,3,2,1,5,5,3,2,2};
    Apriori<int> apriori(data);
    std::vector<std::vector<int>> frequentItemsetFirstpass = apriori.frequentItemSetGen(0.1);
    apriori.printIteration(1, frequentItemsetFirstpass);
    return 0;
}
