#ifndef PROJECT_3_PRIO_QUEUE_H
#define PROJECT_3_PRIO_QUEUE_H
#include "vector.h"
#include <iostream>

typedef struct Node{
    int rowIndex;
    int saturation;
    int degree;
    Node(int row, int sat, int deg): rowIndex(row), saturation(sat), degree(deg){}
    void changeSat(int val){ saturation = val; }
    void changeDeg(int val){ degree = val; }
} Node;


class PrioQueue{
private:
    vector<Node*>* nodes;
    int* positions;
    int size;

    static bool compare(Node* left, Node* right){
        if (left->saturation == right->saturation) {
            if (left->degree == right->degree)
                return left->rowIndex < right->rowIndex;
            return left->degree > right->degree;
        }
        return left->saturation > right->saturation;
    };

    int parent(int i) { return (i - 1) / 2; }

    int leftChild(int i){ return ((2 * i) + 1); }

    static int rightChild(int i){ return ((2 * i) + 2); }

    void shiftUp(int i){
        while (i > 0 && compare(nodes->get(i), nodes->get(parent(i)))) {
            Node* temp = nodes->get(parent(i));
            nodes->place(nodes->get(i), parent(i));
            nodes->place(temp, i);
            int x = positions[nodes->get(i)->rowIndex];
            positions[nodes->get(i)->rowIndex] = positions[nodes->get(parent(i))->rowIndex];
            positions[nodes->get(parent(i))->rowIndex] = x;
            i = parent(i);
        }
    }

    void shiftDown(int i){
        int maxIndex = i;
        int l = leftChild(i);
        if (l < size && compare(nodes->get(l), nodes->get(maxIndex))) maxIndex = l;
        int r = rightChild(i);
        if (r < size && compare(nodes->get(r), nodes->get(maxIndex))) maxIndex = r;
        if (i != maxIndex) {
            Node* temp = nodes->get(maxIndex);
            nodes->place(nodes->get(i), maxIndex);
            nodes->place(temp, i);
            int x = positions[nodes->get(maxIndex)->rowIndex];
            positions[nodes->get(maxIndex)->rowIndex] = positions[nodes->get(i)->rowIndex];
            positions[nodes->get(i)->rowIndex] = x;
            shiftDown(maxIndex);
        }
    }

public:
    PrioQueue(int newSize){
        this->nodes = new vector<Node*>(newSize);
        this->positions = new int[newSize];
        this->size = 0;
    }

    void add(int row, int sat, int deg){
        Node* newNode = new Node(row, sat, deg);
        if (size == 0){
            nodes->place(newNode, size);
            positions[row] = 0;
        }
        else {
            nodes->place(newNode, size);
            positions[row] = size;
            shiftUp(size);
        }
        size++;
    }

    Node* get(){
        Node* highest = nodes->get(0);
        nodes->get(0)->changeDeg(-1);
        nodes->get(0)->changeSat(-1);
        shiftDown(0);
        return highest;
    }

    void changePriority(int row, int newSat){
        Node oldNode = *nodes->get(positions[row]);
        nodes->get(positions[row])->changeSat(newSat);
        if (compare(nodes->get(positions[row]), &oldNode)) shiftUp(positions[row]);
        else shiftDown(positions[row]);
    }
};

#endif //PROJECT_3_PRIO_QUEUE_H
