#ifndef PROJECT_3_VECTOR_H
#define PROJECT_3_VECTOR_H

template <typename T> class vector{
    T* arr;
    unsigned int size;
    unsigned int totalSize;
public:

    vector(){
        size = 0;
        totalSize = 0;
    }

    vector(unsigned int sizeArr){
        arr = new T[sizeArr];
        size = sizeArr;
        totalSize = sizeArr;
    }

    void allocMem(unsigned int wantedSize){
        if (totalSize < wantedSize){
            T* temp_arr = new T[wantedSize];
            totalSize = wantedSize;
            delete[] arr;
            arr = temp_arr;
        }
        size = wantedSize;
    }

    void add(T item){
        T* temp_arr = new T[totalSize + 1];
        for (int i = 0; i < totalSize; ++i) temp_arr[i] = arr[i];
        totalSize++;
        arr = temp_arr;
        arr[totalSize - 1] = item;
    }

    void place(const T& item, int index){
        // if index of place is
        if (index >= totalSize){
            while (index >= totalSize){
                add(item);
            }
        }
        else arr[index] = item;
    }

    T get(unsigned int index) const{
        return arr[index];
    }

    void increment(unsigned int index, int val){
        arr[index] += val;
    }

    unsigned int getSize() const { return size; }

    void changeSize(int val){
        size += val;
    }

    ~vector()= default;
};

#endif //PROJECT_3_VECTOR_H
