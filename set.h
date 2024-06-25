#ifndef PROJECT_3_SET_H
#define PROJECT_3_SET_H


class set{
    unsigned int* arr;
    int size;
    int current;
public:
    set(){
        arr = new unsigned int[100];
        current = 0;
        size = 100;
    }
    set(int newSize){
        arr = new unsigned int[newSize];
        current = 0;
        size = newSize;
    }
    void add(unsigned int newVal){
        // check if value already exists
        for (int i = 0; i < current; ++i) {
            if (arr[i] == newVal) return;
        }
        // if there is not enough place
        if (current == size){
            auto* temp = new unsigned int[size * 2];
            for (int i = 0; i < size; ++i) {
                temp[i] = arr[i];
            }
            size++;
            delete[] arr;
            arr = temp;
        }
        // if not append it
        arr[current] = newVal;
        current++;
    }

    int getSize() { return current; }

    ~set()= default;
};


#endif //PROJECT_3_SET_H
