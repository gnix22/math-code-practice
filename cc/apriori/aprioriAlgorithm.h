#ifndef APRIORI_ALGORITHM_H
#define APRIORI_ALGORITHM_H
#include<iostream>
#include<vector>
template<typename T>
class Apriori{
    private:
        std::vector<T> _items;
        bool _isFrequent;
        int _totalSampleSize;
        int _supCount();
        std::vector<T> _aprioriGen(std::vector<T> prevFrequentItemsets);
    public:
        Apriori(std::vector<T> items); // constructor
        std::vector<T> frequentItemSetGen(double minSup);
};
#include"aprioriAlgorithm.cc"
#endif
