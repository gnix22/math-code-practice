#include"aprioriAlgorithm.h"
#include<vector>
#include<print>

int main(){
    std::vector<int> data = {1,2,3,2,1,1,4,3,2,1,5,5,3,2,2};
    Apriori<int> apriori(data);
    std::vector<int> frequentItemsetFirstpass = apriori.frequentItemSetGen(0.3);
    std::print(frequentItemsetFirstpass);
    return 0;
}
