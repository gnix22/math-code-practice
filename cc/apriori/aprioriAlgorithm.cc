#include"aprioriAlgorithm.h"
#include<vector>
#include<iostream>
#include<algorithm>
template<typename T>
//int Apriori::_supCount(){
//    return std::count(_items, _items.size(), item);
//}
template<typename T>
Apriori::Apriori(std::vector<T> items){
    _items = items;
    _isFrequent = true;
    int _totalSampleSize = 0;
}
template<typename T>
std::vector<T> Apriori::frequentItemSetGen(std::vector<T> items, double minSup){
    std::vector<T> frequentItemsets;
    for(T item : _items){
        if(std::count(_items, _items.size(), item) >= _totalSampleSize * minSup){
            
        }
    }

}

