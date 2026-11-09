#include"aprioriAlgorithm.h"
#include<vector>
int main(){
    std::vector<std::vector<std::string>> transactions = {
        {"Milk", "Bread", "Butter"}, {"Milk", "Bread"}, {"Milk", "Butter"}, {"Bread", "Butter"}, {"Milk", "Bread", "Butter"}
    };
    Apriori<std::string> apriori(transactions);
    auto frequentItemsets = apriori.freqItemsetGeneration(0.6); // 60% min support
    for (const auto& itemset : frequentItemsets){
        for (const auto& item : itemset) std::cout << item << " ";
        std::cout << "\n";
    }
    return 0;
}
