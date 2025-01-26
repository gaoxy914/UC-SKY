#include "BigFloat.h"

template <typename T> class Heap {
// private:
public:
    // unordered_map<int, int> index_map; // id -> index
    vector<int> index_map; // id -> index
    vector<pair<T, int>> heap; // beta, id
    int heap_size = 0;

    void swap(int i, int j) {
        auto item1 = heap[i];
        auto item2 = heap[j];
        heap[i] = item2;
        heap[j] = item1;
        index_map[item2.second] = i;
        index_map[item1.second] = j;
    }

    void adjust_up(int child) {
        int parent = (child - 1)/2;
        while (child > 0) {
            if (heap[parent].first < heap[child].first) {
                swap(parent, child);
                child = parent;
                parent = (child - 1)/2;
            } else {
                break;
            }
        }
    }

    void adjust_down(int parent) {
        int child = parent*2 + 1;
        while (child < heap_size) {
            if (child + 1 < heap_size && heap[child].first < heap[child + 1].first) {
                ++ child;
            }
            if (heap[parent].first < heap[child].first) {
                swap(parent, child);
                parent = child;
                child = parent*2 + 1;
            } else {
                break;
            }
        }
    }

    void resign(const pair<T, int>& item) {
        adjust_up(index_map[item.second]);
        adjust_down(index_map[item.second]);
    }
public:
    Heap(const int& n) {
        index_map.resize(n, -1);
    }
    ~Heap() {
        vector<pair<T, int>> ().swap(heap);
        vector<int> ().swap(index_map);
    }

    void push(const pair<T, int>& item) {
        heap.push_back(item);
        index_map[item.second] = heap_size;
        adjust_up(heap_size ++);
    }

    void pop() {
        if (heap_size > 0) {
            auto& item = heap[0];
            swap(0, heap_size - 1);
            heap_size --;
            // index_map.erase(item.second);
            index_map[item.second] = -1;
            heap.pop_back();
            adjust_down(0);
        }
    }

    void set(const pair<T, int>& item, const pair<T, int>& target) {
        int index = index_map[item.second];
        heap[index] = target;
        adjust_up(index);
        adjust_down(index);
    }

    T get(const int& id) {
        int index = index_map[id];
        return heap[index].first;
    }

    void divide(const int& id, const T& prob) {
        int index = index_map[id];
        heap[index].first = heap[index].first / prob;
        adjust_up(index);
        adjust_down(index);
    }

    void add(const int& id, const T& delta) {
        int index = index_map[id];
        heap[index].first = heap[index].first + delta;
        adjust_up(index);
        adjust_down(index);
    }

    void remove(const int& id) {
        int index = index_map[id];
        swap(index, heap_size - 1);
        heap_size --;
        // index_map.erase(id);
        index_map[id] = -1;
        heap.pop_back();
        if (index == heap_size) return;
        adjust_up(index);
        adjust_down(index);
    }

    pair<T, int> top() {
        if (heap_size > 0) {
            return heap[0];
        }
    }

    int size() {
        return heap_size;
    }

    bool empty() {
        return heap.empty();
    }
};