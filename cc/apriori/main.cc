#include"aprioriAlgorithm.h"
#include<vector>
int main(){
    std::vector<std::vector<std::string>> transactions = {
        {"bread", "milk"},
        {"bread", "diaper", "beer", "eggs"},
        {"milk", "diaper", "beer", "cola"},
        {"bread", "milk", "diaper", "beer"},
        {"bread", "milk", "diaper", "cola"}
    };
    Apriori<std::string> apriori(transactions);
    auto frequentItemsets = apriori.freqItemsetGeneration(0.3); // 60% min support
    for (const auto& itemset : frequentItemsets){
        for (const auto& item : itemset) std::cout << item << " ";
        std::cout << "\n";
    }
    return 0;
}
