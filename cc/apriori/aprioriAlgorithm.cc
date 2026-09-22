#include"aprioriAlgorithm.h"
#include<vector>
#include<iostream>
#include<algorithm>
template<typename T>
Apriori<T>::Apriori(std::vector<T> items){
    _items = items;
    _isFrequent = true;
    int _totalSampleSize = items.size();
}
template<typename T>
std::vector<T> Apriori<T>::frequentItemSetGen(double minSup){
    std::vector<T> frequentItemsets;
    for(T item : _items){
        if(std::count(_items.begin(), _items.end(), item) >= _totalSampleSize * minSup){
            frequentItemsets.push_back(item);
        }
    }
    return frequentItemsets;
}

