#ifndef APRIORI_ALGORITHM_H
#define APRIORI_ALGORITHM_H
#include<iostream>
#include<vector>
#include<algorithm>
#include<iostream>
#include<set>
template<typename T>
class Apriori{
    private:
        std::vector<T> _items;
        bool _isFrequent;
        int _totalSampleSize;
        std::vector<T> _aprioriGen(std::vector<T> prevFrequentItemsets);
        int _kCount;
    public:
        // constructor
        Apriori(std::vector<T> items){
            _items = items;
            _isFrequent = true;
            _totalSampleSize = _items.size();
        }
        std::vector<std::vector<T>> frequentItemSetGen(double minSup){
            std::vector<std::vector<T>> frequentItemsets;
            std::set<T> unique_items(_items.begin(), _items.end());
            for(T item : unique_items){
                if(std::count(_items.begin(), _items.end(), item) >= _totalSampleSize * minSup){
                    frequentItemsets.push_back(std::vector{item});
                }
            }
            return frequentItemsets;
        }
        void printIteration(int k, std::vector<std::vector<T>>currFreqItemset){
            if(currFreqItemset.empty()){
                std::cout << "no itemsets created that matched given support. try lowering support.";
            }else{
                std::cout << "frequent " << k << "-itemsets: ";
                for(int i=0; i<currFreqItemset.size(); i++){
                    std::cout << "{";
                    for(int j=0; j<currFreqItemset[i].size(); j++){
                        std::cout << currFreqItemset[i].at(j);
                        if(j < currFreqItemset[i].size() - 1){
                            std::cout << ", ";
                        }
                    }
                    if(i < currFreqItemset.size()-1){
                        std::cout << "}, ";
                    } else{
                        std::cout << "}";
                    }
                }
            }
        }
};
#endif
