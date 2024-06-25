#include <iostream>
#include "vector.h"
#include "set.h"
#include "queue.h"
#include "prio_queue.h"

// 1 (Degree sequence) Begin
void getDegreeSequence(const vector<unsigned int> &graphDegrees, vector<unsigned int> *sortedDegrees, int max,
                       vector<unsigned int> *sortedIndexes) {
    for (int i = 0; i < graphDegrees.getSize(); ++i) {
        sortedDegrees->increment(graphDegrees.get(i), 1);
    }
    for (int i = 1; i <= max; ++i) {
        sortedDegrees->increment(i, sortedDegrees->get(i - 1));
    }
    auto *output = new int[graphDegrees.getSize()];
    for (int i = 0; i < graphDegrees.getSize(); i++) {
        output[sortedDegrees->get(graphDegrees.get(i)) - 1] = graphDegrees.get(i);
        sortedIndexes->place(i, sortedDegrees->get(graphDegrees.get(i)) - 1);
        sortedDegrees->increment(graphDegrees.get(i), -1);
    }
    for (int i = graphDegrees.getSize() - 1; i >= 0; i--) {
        printf("%d ", output[i]);
    }
    printf("\n");
    delete[] output;
}
// End

// 2 (Components) Begin
void DFS(const vector<vector<unsigned int> *> &vertexEdges, int begin, bool *visited) {
    visited[begin] = true;
    for (int i = 0; i < vertexEdges.get(begin)->getSize(); i++) {
        if (!visited[vertexEdges.get(begin)->get(i)]) {
            DFS(vertexEdges, vertexEdges.get(begin)->get(i), visited);
        }
    }
}

void getAmountOfComponents(const vector<vector<unsigned int> *> &vertexEdges, int amountOfVertices) {
    int amount = 0;
    bool *visited = new bool[amountOfVertices];
    for (int i = 0; i < amountOfVertices; ++i) visited[i] = false;

    for (int i = 0; i < amountOfVertices; ++i) {
        if (!visited[i]) {
            DFS(vertexEdges, i, visited);
            amount++;
        }
    }
    printf("%d\n", amount);
    delete[] visited;
}
// End

// 3 (Bipartite) Begin
bool BFS(const vector<vector<unsigned int> *> &graph, int *color, int start, queue<int> *q) {
    int otherVertex;
    q->push_back(start);
    while (!q->is_empty()) {
        int front = q->get_front();
        for (int i = 0; i < graph.get(front)->getSize(); ++i) {
            otherVertex = graph.get(front)->get(i);
            if (color[otherVertex] == -1) {
                color[otherVertex] = 1 - color[front];
                q->push_back(otherVertex);
            } else if (color[otherVertex] == color[front]) return false;
        }
    }
    return true;
}

void isBipartite(const vector<vector<unsigned int> *> &graph) {
    // BFS with some modifications
    int *color = new int[graph.getSize()];
    // set every element to be uncolored
    for (int i = 0; i < graph.getSize(); ++i) {
        color[i] = -1;
    }
    auto *q = new queue<int>(graph.getSize());
    for (int i = 0; i < graph.getSize(); ++i) {
        // do BFS for every vertex which does not have a color yet
        // if it has a color it means that this component has been searched
        // BFS will return false if graph is not bipartite
        if (color[i] == -1) {
            if (!BFS(graph, color, i, q)) {
                printf("F\n");
                return;
            }
        }
    }
    printf("T\n");
    delete q;
    delete[] color;
}
// End

// 6 (Greedy) Begin
void performGreedyColoring(const vector<vector<unsigned int> *> &graph) {
    int *colored = new int[graph.getSize()];
    bool *notColored = new bool[graph.getSize()];
    int color;

    for (int i = 0; i < graph.getSize(); ++i) {
        colored[i] = -1;
        notColored[i] = false;
    }

    colored[0] = 0;

    for (int i = 1; i < graph.getSize(); ++i) {
        for (int j = 0; j < graph.get(i)->getSize(); ++j) {
            if (colored[graph.get(i)->get(j)] != -1)
                notColored[colored[graph.get(i)->get(j)]] = true;
        }
        for (color = 0; color < graph.getSize(); color++) {
            if (!notColored[color]) break;
        }

        colored[i] = color;

        for (int j = 0; j < graph.get(i)->getSize(); ++j) {
            if (colored[graph.get(i)->get(j)] != -1) notColored[colored[graph.get(i)->get(j)]] = false;
        }
    }

    for (int i = 0; i < graph.getSize(); ++i) printf("%d ", colored[i] + 1);
    printf("\n");
    delete[] colored;
    delete[] notColored;
}
// End

// 7 (LF) Begin
void performLFColoring(const vector<vector<unsigned int> *> &graph, const vector<unsigned int> &reversedOrder) {
    int *colored = new int[graph.getSize()];
    bool *notColored = new bool[graph.getSize()];
    int color;
    int node;
    for (int i = 0; i < graph.getSize(); ++i) {
        colored[i] = -1;
        notColored[i] = false;
    }
    colored[reversedOrder.get(reversedOrder.getSize() - 1)] = 0;

    for (int i = reversedOrder.getSize() - 2; i >= 0; i--) {
        node = reversedOrder.get(i);

        for (int j = 0; j < graph.get(node)->getSize(); ++j) {
            if (colored[graph.get(node)->get(j)] != -1)
                notColored[colored[graph.get(node)->get(j)]] = true;
        }

        for (color = 0; color < graph.getSize(); color++) {
            if (!notColored[color]) break;
        }
        colored[node] = color;

        for (int j = 0; j < graph.get(node)->getSize(); ++j) {
            if (colored[graph.get(node)->get(j)] != -1) notColored[colored[graph.get(node)->get(j)]] = false;
        }
    }

    for (int i = 0; i < graph.getSize(); ++i) printf("%d ", colored[i] + 1);
    printf("\n");
    delete[] colored;
    delete[] notColored;
}
// End

// 8 (SLF) Begin
void performSLFColoring(const vector<vector<unsigned int> *> &graph,
                        vector<unsigned int> *degrees) {
    // array with color number of vertex
    int *colored = new int[graph.getSize()];
    // array that says if vetex is 'free'
    bool *notColored = new bool[graph.getSize()];
    auto *sat = new vector<set *>(graph.getSize());
    // color number
    int color;
    int node;
    int x = 0;
    int idx;
    auto *queue = new PrioQueue(graph.getSize());

    for (int i = 0; i < graph.getSize(); ++i) {
        colored[i] = -1;
        notColored[i] = false;
        sat->place(new set(graph.get(i)->getSize()), i);
        queue->add(i, 0, degrees->get(i));
    }
    while (x != graph.getSize()) {
        node = queue->get()->rowIndex;
        x++;
        // go through all adjacent vertices
        for (int j = 0; j < graph.get(node)->getSize(); ++j) {
            idx = graph.get(node)->get(j);
            // set neighbours to true if they are already colored
            if (colored[idx] != -1)
                notColored[colored[idx]] = true;
        }
        // search for the lowest possible color
        for (color = 0; color < graph.get(node)->getSize(); color++) {
            if (!notColored[color]) break;
        }
        // set color
        colored[node] = color;
        // clear neighbouring array
        for (int j = 0; j < graph.get(node)->getSize(); ++j) {
            idx = graph.get(node)->get(j);
            if (colored[idx] != -1) {
                notColored[colored[idx]] = false;
            } else {
                // add saturation to every adjacent vertex
                sat->get(idx)->add(color);
            }
        }
        // update tree
        for (int i = 0; i < graph.get(node)->getSize(); ++i) {
            idx = graph.get(node)->get(i);
            // for every vertex that is not colored
            if (colored[idx] == -1) {
                queue->changePriority(idx, sat->get(idx)->getSize());
            }
        }
    }
    for (int i = 0; i < graph.getSize(); ++i) printf("%d ", colored[i] + 1);
    printf("\n");
    delete queue;
    delete sat;
    delete[] colored;
    delete[] notColored;
}
// End

// 9 cycles Begin
struct C4{
    unsigned int vertex;
    unsigned int degree;

    bool operator<(const C4 &adder) const{
        if (degree == adder.degree) return vertex < adder.vertex;
        return degree < adder.degree;
    }

    bool operator>(const C4 &adder) const{
        if (degree == adder.degree) return vertex > adder.vertex;
        return degree > adder.degree;
    }

    bool compare(const C4 &other){
        if (degree == other.degree) return vertex > other.vertex;
        return degree > other.degree;
    }
};

void quickSort(vector<C4*>* adjNodes, int left, int right){
    if (right <= left) return;
    int i = left - 1;
    int j = right + 1;
    C4* mid = adjNodes->get((left + right) / 2);

    while (true){

        while(i < right && *mid > *adjNodes->get(++i));
//        while(*mid > *adjNodes->get(++i));
        while(j > left && *mid < *adjNodes->get(--j));
//        while(*mid < *adjNodes->get(--j));

        if ( i <= j){
            C4* temp = adjNodes->get(i);
            adjNodes->place(adjNodes->get(j), i);
            adjNodes->place(temp, j);
        }
        else break;
    }

    if (j > left) quickSort(adjNodes, left, j);
    if (i < right) quickSort(adjNodes, i, right);
}

void betterC4Counting(vector<vector<unsigned int> *> &graph) {
    int size = graph.getSize();
    auto* sortedDegrees = new vector<vector<C4*>*>(size);
    for (int i = 0; i < size; ++i) {
        auto* degrees = new vector<C4*>(graph.get(i)->getSize());
        for (int j = 0; j < graph.get(i)->getSize(); ++j) {
            auto* row = new C4{graph.get(i)->get(j), graph.get(graph.get(i)->get(j))->getSize()};
            degrees->place(row, j);
        }
        quickSort(degrees, 0, (int)degrees->getSize() - 1);
        sortedDegrees->place(degrees, i);
    }
    int y;
    int idx;
    int current;
    int degreeV;
    long long numberOfCycles = 0;
    C4 curr;
    int *L = new int[size];
    for (int i = 0; i < size; ++i) L[i] = 0;
    // for every vertex
    for (int v = 0; v < size; v++) {
        curr = C4{(unsigned int)v, sortedDegrees->get(v)->getSize()};
        degreeV = (int)graph.get(v)->getSize();
        // for every neighbouring vertex that is < than i (root of neighbours)
        for (int u = 0; u < degreeV; ++u) {
            if (sortedDegrees->get(v)->get(u)->compare(curr)) break;
            idx = 0;
            // current vertex index
            current = sortedDegrees->get(v)->get(u)->vertex;
            // first neighbour of vertex u
            y = sortedDegrees->get(current)->get(idx)->vertex;
            while (y != v) {
                numberOfCycles += L[y];
                L[y]++;
                idx++;
                y = sortedDegrees->get(current)->get(idx)->vertex;
            }
        }
        for (int u = 0; u < degreeV; ++u) {
//            if (!(*sortedDegrees->get(v)->get(u) < curr)) break;
            if (sortedDegrees->get(v)->get(u)->compare(curr)) break;
            idx = 0;
            current = sortedDegrees->get(v)->get(u)->vertex;
            y = sortedDegrees->get(current)->get(idx)->vertex;
            while (y != v){
                L[y] = 0;
                idx++;
                y = sortedDegrees->get(current)->get(idx)->vertex;
            }
        }
    }
    printf("%lld\n", numberOfCycles);
    for (int i = 0; i < sortedDegrees->getSize(); ++i) {
        delete sortedDegrees->get(i);
    }
    delete sortedDegrees;
    delete[] L;
}
// End


int main() {
    unsigned int numberOfGraphs = 0, numberOfVertices = 0, degreeValue = 0, edge = 0, maxDegree = 0;
    scanf("%d", &numberOfGraphs);
    unsigned long long complementedEdges, presentEdges;
    auto *degreeSequenceVector = new vector<unsigned int>();
    auto *degreeSequenceVectorZeros = new vector<unsigned int>();
    auto *graph = new vector<vector<unsigned int> *>();
    auto *degreeIndexes = new vector<unsigned int>();
    // iterate through all graphs
    for (int i = 0; i < numberOfGraphs; ++i) {
        complementedEdges = 0;
        presentEdges = 0;
        scanf("%d", &numberOfVertices);
        // vector containing graph
        graph->allocMem(numberOfVertices);
        degreeSequenceVector->allocMem(numberOfVertices);
        degreeSequenceVectorZeros->allocMem(numberOfVertices);
        degreeIndexes->allocMem(numberOfVertices);
        // vector containing degrees of vectors
        for (int j = 0; j < numberOfVertices; ++j) {
            // add "row" to graph
            scanf("%d", &degreeValue);
            auto *newV = new vector<unsigned int>();
            newV->allocMem(degreeValue);
            graph->place(newV, j);
            // add degree to vector
            degreeSequenceVector->place(degreeValue, j);
            degreeSequenceVectorZeros->place(0, j);
            // search for max value of degree
            if (maxDegree < degreeValue) maxDegree = degreeValue;
            // add to complemented edges the triangle of ajd matrix
            complementedEdges += numberOfVertices - j - 1;
            presentEdges += degreeValue;
            for (int k = 0; k < degreeValue; ++k) {
                scanf("%d", &edge);
                graph->get(j)->place(edge - 1, k);
            }
        }
        getDegreeSequence(*degreeSequenceVector, degreeSequenceVectorZeros, maxDegree, degreeIndexes); // 1 degree
        getAmountOfComponents(*graph, numberOfVertices); // 2 components
        isBipartite(*graph); // 3 bipartitenes
        performGreedyColoring(*graph); // 6 greedy
        performLFColoring(*graph, *degreeIndexes); // 7 lf
        performSLFColoring(*graph, degreeSequenceVector); // 8 slf
        betterC4Counting(*graph); // 9 count C4
        printf("%llu\n", complementedEdges - (presentEdges / 2)); // 10
    }
    delete degreeSequenceVector;
    delete degreeSequenceVectorZeros;
    delete degreeIndexes;
    for (int i = 0; i < graph->getSize(); ++i) {
        delete graph->get(i);
    }
    delete graph;
    return 0;
}