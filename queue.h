#ifndef PROJECT_3_QUEUE_H
#define PROJECT_3_QUEUE_H

template <typename T> class queue{
    T* arr;
    unsigned int size;
    unsigned int currentGetIndex;
    unsigned int currentPushIndex;
public:
    queue(){
        arr = nullptr;
        size = 0;
        currentGetIndex = 0;
        currentPushIndex = 0;
    }

    queue(unsigned int newSize){
        arr = new T[newSize];
        size = newSize;
        currentGetIndex = 0;
        currentPushIndex = 0;
    }

    void push_back(T item){
        arr[currentPushIndex] = item;
        currentPushIndex++;
    }

    T get_front(){
        currentGetIndex++;
        return arr[currentGetIndex - 1];
    }

    bool is_empty(){
        return currentGetIndex == currentPushIndex;
    }
};

#endif //PROJECT_3_QUEUE_H
