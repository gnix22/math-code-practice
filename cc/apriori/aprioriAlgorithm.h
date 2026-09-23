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
        std::vector<std::vector<T>> _data; // the set I containing i itemsets.
        int _totalSampleSize;
        // k = 1
        // F_k = {i|i in I and supp({i})>= N * minsup}
        std::vector<std::vector<T>> firstFreqItemsetsGen(double minSup){
            std::set<T> unique_items;
            for (const auto& transaction : _data)
                unique_items.insert(transaction.begin(), transaction.end());
            std::vector<std::vector<T>> frequentItemsets;
            for (const T& item : unique_items){
                std::vector<T> candidate{item};
                int support = 0;
                for (const auto& t : _data){
                    if (std::find(t.begin(), t.end(), item) != t.end())
                        ++support;
                }
                if (support >= _totalSampleSize * minSup)
                    frequentItemsets.push_back(candidate);
            }
            return frequentItemsets;
        }
        // C_k = apriori_gen(f_{k-1})  -- join + prune only, no support counting
        std::vector<std::vector<T>> aprioriGenCandidatesOnly(std::vector<std::vector<T>> prevFreqItemset){
            for (auto& itemset : prevFreqItemset)
                std::sort(itemset.begin(), itemset.end());
            std::sort(prevFreqItemset.begin(), prevFreqItemset.end());
            std::vector<std::vector<T>> candidates;
            // joining step
            for (size_t i = 0; i < prevFreqItemset.size(); ++i){
                for (size_t j = i + 1; j < prevFreqItemset.size(); ++j){
                    const auto& a = prevFreqItemset[i];
                    const auto& b = prevFreqItemset[j];
                    if (!std::equal(a.begin(), a.end() - 1, b.begin())) continue;
                    if (a.back() >= b.back()) continue;
                    std::vector<T> candidate = a;
                    candidate.push_back(b.back());
                    // candidate pruning step: every (k-1)-subset must already be in prevFreqItemset
                    bool allSubsetsFrequent = true;
                    for (size_t skip = 0; skip < candidate.size(); ++skip){
                        std::vector<T> subset;
                        for (size_t idx = 0; idx < candidate.size(); ++idx)
                            if (idx != skip) subset.push_back(candidate[idx]);
                        if (std::find(prevFreqItemset.begin(), prevFreqItemset.end(), subset)
                            == prevFreqItemset.end()){
                            allSubsetsFrequent = false;
                            break;
                        }
                    }
                    if (allSubsetsFrequent)
                        candidates.push_back(candidate);
                }
            }
            return candidates;
        }
    public:
        // constructor
        Apriori(std::vector<std::vector<T>> itemset){
            _data = itemset;
            _totalSampleSize = _data.size();
        }
        std::vector<std::vector<T>> freqItemsetGeneration(double minSup){
            std::vector<std::vector<T>> resultItemsets;
            int k = 1;
            // F_k
            std::vector<std::vector<T>> firstFreqItemsets = firstFreqItemsetsGen(minSup);
            // insert itemsets at end of result
            resultItemsets.insert(resultItemsets.end(), firstFreqItemsets.begin(), firstFreqItemsets.end());
            do{
                k++; // k = k+1
                // C_k = apriori-gen(F_k-1)
                std::vector<std::vector<T>> candidateItemsets = aprioriGenCandidatesOnly(firstFreqItemsets);
                // for each transaction t in T:
                // for each candidate c in C_k:
                // if c subset of t: support_c += 1
                std::vector<int> supportCounts(candidateItemsets.size(), 0);
                for (const auto& t : _data){
                    std::set<T> tset(t.begin(), t.end());
                    for (size_t idx = 0; idx < candidateItemsets.size(); ++idx){
                        const auto& c = candidateItemsets[idx];
                        bool isSubset = std::all_of(c.begin(), c.end(),
                            [&](const T& item){return tset.count(item) > 0;});
                        if (isSubset) supportCounts[idx]++;
                    }
                }
                // f_k = { c in C_k : support_c >= N*minsup }
                firstFreqItemsets.clear();
                for (size_t idx = 0; idx < candidateItemsets.size(); ++idx){
                    if (supportCounts[idx] >= _totalSampleSize * minSup)
                        candidateItemsets.push_back(candidateItemsets[idx]);
                }
                resultItemsets.insert(resultItemsets.end(), firstFreqItemsets.begin(), firstFreqItemsets.end());
            } while (!firstFreqItemsets.empty());
            return resultItemsets;
        }

};
#endif
