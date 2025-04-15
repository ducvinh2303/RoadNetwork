#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>
#include <string.h>

#define INF INT_MAX  // Infinite value for comparison
#define MAX_CITIES 26 // Maximum number of cities (A-Z)

// Structure for a road between cities
typedef struct Road {
    int destination;   // Destination city index
    int distance;      // Distance in km
    struct Road* next; // Next road in the list
} Road;

// Structure for a city
typedef struct City {
    char name;                  // City name (A-Z)
    char* description;          // City description
    int population;             // City population
    int landmarks;              // Number of landmarks
    Road* roads;                // Roads from this city (one-to-many relationship)
} City;

// Structure for a region (represents groups of cities)
typedef struct Region {
    char name[50];              // Region name
    City* cities[MAX_CITIES];   // Cities in this region (another one-to-many relationship)
    int cityCount;              // Number of cities in the region
} Region;

// Implementation 1: Adjacency list using linked lists
typedef struct Graph_List {
    int numVertices;
    Road** adjLists;   // Array of linked lists
} Graph_List;

// Implementation 2: Adjacency matrix using arrays 
typedef struct Graph_Matrix {
    int numVertices;
    int** adjMatrix;   // 2D array for edge weights
} Graph_Matrix;

// Structure for MinHeap node used in Dijkstra's algorithm
typedef struct MinHeapNode {
    int vertex;
    int distance;
} MinHeapNode;

// Structure for MinHeap
typedef struct MinHeap {
    int size;
    int capacity;
    int* pos;
    MinHeapNode** array;
} MinHeap;

// ---- City and Road functions ----

// Create a new city
City* createCity(char name, const char* description, int population, int landmarks) {
    City* newCity = (City*)malloc(sizeof(City));
    if (!newCity) return NULL;
    
    newCity->name = name;
    newCity->population = population;
    newCity->landmarks = landmarks;
    newCity->roads = NULL;
    
    // Allocate and copy description
    if (description) {
        newCity->description = (char*)malloc(strlen(description) + 1);
        if (newCity->description) {
            strcpy(newCity->description, description);
        }
    } else {
        newCity->description = NULL;
    }
    
    return newCity;
}

// Free a city and its roads
void freeCity(City* city) {
    if (!city) return;
    
    // Free the description
    if (city->description) {
        free(city->description);
    }
    
    // Free all roads
    Road* currentRoad = city->roads;
    while (currentRoad) {
        Road* nextRoad = currentRoad->next;
        free(currentRoad);
        currentRoad = nextRoad;
    }
    
    free(city);
}

// Create a new road
Road* createRoad(int destination, int distance) {
    Road* newRoad = (Road*)malloc(sizeof(Road));
    if (!newRoad) return NULL;
    
    newRoad->destination = destination;
    newRoad->distance = distance;
    newRoad->next = NULL;
    
    return newRoad;
}

// Add a road to a city
void addRoadToCity(City* city, int destination, int distance) {
    if (!city) return;
    
    Road* newRoad = createRoad(destination, distance);
    if (!newRoad) return;
    
    // Add at the beginning of the list
    newRoad->next = city->roads;
    city->roads = newRoad;
}

// ---- Region functions ----

// Create a new region
Region* createRegion(const char* name) {
    Region* newRegion = (Region*)malloc(sizeof(Region));
    if (!newRegion) return NULL;
    
    strncpy(newRegion->name, name, 49);
    newRegion->name[49] = '\0';
    newRegion->cityCount = 0;
    
    for (int i = 0; i < MAX_CITIES; i++) {
        newRegion->cities[i] = NULL;
    }
    
    return newRegion;
}

// Add a city to a region
void addCityToRegion(Region* region, City* city) {
    if (!region || !city || region->cityCount >= MAX_CITIES) return;
    
    region->cities[region->cityCount++] = city;
}

// Free a region and all its cities
void freeRegion(Region* region) {
    if (!region) return;
    
    for (int i = 0; i < region->cityCount; i++) {
        freeCity(region->cities[i]);
    }
    
    free(region);
}

// ---- Graph implementation using adjacency lists ----

// Create a new graph using adjacency lists
Graph_List* createGraphList(int numVertices) {
    Graph_List* graph = (Graph_List*)malloc(sizeof(Graph_List));
    if (!graph) return NULL;
    
    graph->numVertices = numVertices;
    graph->adjLists = (Road**)malloc(numVertices * sizeof(Road*));
    
    if (!graph->adjLists) {
        free(graph);
        return NULL;
    }
    
    for (int i = 0; i < numVertices; i++) {
        graph->adjLists[i] = NULL;
    }
    
    return graph;
}

// Add an edge to graph using adjacency lists
void addEdgeList(Graph_List* graph, int src, int dest, int weight) {
    if (!graph || src < 0 || src >= graph->numVertices || 
        dest < 0 || dest >= graph->numVertices) return;

    // Add edge from src to dest
    Road* newRoad = createRoad(dest, weight);
    newRoad->next = graph->adjLists[src];
    graph->adjLists[src] = newRoad;
    
    // Add edge from dest to src (undirected graph)
    newRoad = createRoad(src, weight);
    newRoad->next = graph->adjLists[dest];
    graph->adjLists[dest] = newRoad;
}

// Free graph using adjacency lists
void freeGraphList(Graph_List* graph) {
    if (!graph) return;
    
    if (graph->adjLists) {
        for (int i = 0; i < graph->numVertices; i++) {
            Road* current = graph->adjLists[i];
            while (current) {
                Road* next = current->next;
                free(current);
                current = next;
            }
        }
        free(graph->adjLists);
    }
    
    free(graph);
}

// ---- Graph implementation using adjacency matrix ----

// Create a new graph using adjacency matrix
Graph_Matrix* createGraphMatrix(int numVertices) {
    Graph_Matrix* graph = (Graph_Matrix*)malloc(sizeof(Graph_Matrix));
    if (!graph) return NULL;
    
    graph->numVertices = numVertices;
    graph->adjMatrix = (int**)malloc(numVertices * sizeof(int*));
    
    if (!graph->adjMatrix) {
        free(graph);
        return NULL;
    }
    
    for (int i = 0; i < numVertices; i++) {
        graph->adjMatrix[i] = (int*)malloc(numVertices * sizeof(int));
        if (!graph->adjMatrix[i]) {
            // Clean up previously allocated memory
            for (int j = 0; j < i; j++) {
                free(graph->adjMatrix[j]);
            }
            free(graph->adjMatrix);
            free(graph);
            return NULL;
        }
        
        for (int j = 0; j < numVertices; j++) {
            graph->adjMatrix[i][j] = 0; // Initialize with no edges
        }
    }
    
    return graph;
}

// Add an edge to graph using adjacency matrix
void addEdgeMatrix(Graph_Matrix* graph, int src, int dest, int weight) {
    if (!graph || src < 0 || src >= graph->numVertices || 
        dest < 0 || dest >= graph->numVertices) return;

    // Add edge (undirected graph)
    graph->adjMatrix[src][dest] = weight;
    graph->adjMatrix[dest][src] = weight;
}

// Free graph using adjacency matrix
void freeGraphMatrix(Graph_Matrix* graph) {
    if (!graph) return;
    
    if (graph->adjMatrix) {
        for (int i = 0; i < graph->numVertices; i++) {
            free(graph->adjMatrix[i]);
        }
        free(graph->adjMatrix);
    }
    
    free(graph);
}

// ---- MinHeap functions for Dijkstra's algorithm ----

// Create a new MinHeapNode
MinHeapNode* createMinHeapNode(int vertex, int distance) {
    MinHeapNode* node = (MinHeapNode*)malloc(sizeof(MinHeapNode));
    if (!node) return NULL;
    
    node->vertex = vertex;
    node->distance = distance;
    
    return node;
}

// Create a MinHeap
MinHeap* createMinHeap(int capacity) {
    MinHeap* minHeap = (MinHeap*)malloc(sizeof(MinHeap));
    if (!minHeap) return NULL;
    
    minHeap->pos = (int*)malloc(capacity * sizeof(int));
    if (!minHeap->pos) {
        free(minHeap);
        return NULL;
    }
    
    minHeap->size = 0;
    minHeap->capacity = capacity;
    minHeap->array = (MinHeapNode**)malloc(capacity * sizeof(MinHeapNode*));
    
    if (!minHeap->array) {
        free(minHeap->pos);
        free(minHeap);
        return NULL;
    }
    
    return minHeap;
}

// Swap two MinHeap nodes
void swapMinHeapNode(MinHeapNode** a, MinHeapNode** b) {
    MinHeapNode* temp = *a;
    *a = *b;
    *b = temp;
}

// Heapify at a given index
void minHeapify(MinHeap* minHeap, int idx) {
    int smallest, left, right;
    smallest = idx;
    left = 2 * idx + 1;
    right = 2 * idx + 2;

    if (left < minHeap->size &&
        minHeap->array[left]->distance < minHeap->array[smallest]->distance)
        smallest = left;

    if (right < minHeap->size &&
        minHeap->array[right]->distance < minHeap->array[smallest]->distance)
        smallest = right;

    if (smallest != idx) {
        // Update position of nodes in position array
        minHeap->pos[minHeap->array[smallest]->vertex] = idx;
        minHeap->pos[minHeap->array[idx]->vertex] = smallest;

        // Swap nodes
        swapMinHeapNode(&minHeap->array[smallest], &minHeap->array[idx]);

        minHeapify(minHeap, smallest);
    }
}

// Check if MinHeap is empty
int isEmptyMinHeap(MinHeap* minHeap) {
    return minHeap->size == 0;
}

// Extract minimum node from MinHeap
MinHeapNode* extractMin(MinHeap* minHeap) {
    if (isEmptyMinHeap(minHeap))
        return NULL;

    // Store the root node
    MinHeapNode* root = minHeap->array[0];
    MinHeapNode* lastNode = minHeap->array[minHeap->size - 1];
    
    // Replace root with last node
    minHeap->array[0] = lastNode;

    // Update position of last node
    minHeap->pos[root->vertex] = minHeap->size - 1;
    minHeap->pos[lastNode->vertex] = 0;

    // Reduce heap size and heapify root
    --minHeap->size;
    minHeapify(minHeap, 0);

    return root;
}

// Decrease key value of a given vertex
void decreaseKey(MinHeap* minHeap, int vertex, int dist) {
    // Get the index of vertex in heap array
    int i = minHeap->pos[vertex];

    // Update distance value
    minHeap->array[i]->distance = dist;

    // Heapify up until the proper position is found
    while (i && minHeap->array[i]->distance < minHeap->array[(i - 1) / 2]->distance) {
        // Swap with parent
        minHeap->pos[minHeap->array[i]->vertex] = (i - 1) / 2;
        minHeap->pos[minHeap->array[(i - 1) / 2]->vertex] = i;
        swapMinHeapNode(&minHeap->array[i], &minHeap->array[(i - 1) / 2]);

        // Move to parent index
        i = (i - 1) / 2;
    }
}

// Check if a vertex is in MinHeap
int isInMinHeap(MinHeap* minHeap, int vertex) {
    return minHeap->pos[vertex] < minHeap->size;
}

// Free MinHeap
void freeMinHeap(MinHeap* minHeap) {
    if (!minHeap) return;
    
    if (minHeap->array) {
        for (int i = 0; i < minHeap->size; i++) {
            free(minHeap->array[i]);
        }
        free(minHeap->array);
    }
    
    if (minHeap->pos) {
        free(minHeap->pos);
    }
    
    free(minHeap);
}

// ---- Path-finding algorithms ----

// Print the path from source to j using parent array
void printPath(int* parent, int j, char baseChar) {
    if (parent[j] == -1)
        return;

    printPath(parent, parent[j], baseChar);
    printf("%c ", j + baseChar);
}

// Dijkstra's algorithm using adjacency list
void dijkstraList(Graph_List* graph, int src, int dest, char baseChar, double* executionTime) {
    clock_t start, end;
    start = clock();
    
    int V = graph->numVertices;
    int* dist = (int*)malloc(V * sizeof(int));
    int* parent = (int*)malloc(V * sizeof(int));
    
    // Create MinHeap
    MinHeap* minHeap = createMinHeap(V);

    // Initialize distances and MinHeap
    for (int v = 0; v < V; v++) {
        dist[v] = INF;
        parent[v] = -1;
        minHeap->array[v] = createMinHeapNode(v, dist[v]);
        minHeap->pos[v] = v;
    }

    // Set source distance to 0
    dist[src] = 0;
    decreaseKey(minHeap, src, dist[src]);
    minHeap->size = V;

    // Process vertices
    while (!isEmptyMinHeap(minHeap)) {
        MinHeapNode* minNode = extractMin(minHeap);
        int u = minNode->vertex;
        
        // Process adjacent vertices
        Road* road = graph->adjLists[u];
        while (road) {
            int v = road->destination;
            
            // Update distance if shorter path found
            if (isInMinHeap(minHeap, v) && dist[u] != INF && road->distance + dist[u] < dist[v]) {
                dist[v] = dist[u] + road->distance;
                parent[v] = u;
                decreaseKey(minHeap, v, dist[v]);
            }
            road = road->next;
        }
        
        free(minNode);
    }
    
    // Print result
    if (dist[dest] != INF) {
        printf("\nShortest path from %c to %c using Adjacency List: ", src + baseChar, dest + baseChar);
        printf("%c ", src + baseChar);
        printPath(parent, dest, baseChar);
        printf("\nTotal distance: %d km\n", dist[dest]);
    } else {
        printf("\nNo path exists from %c to %c using Adjacency List\n", src + baseChar, dest + baseChar);
    }
    
    end = clock();
    *executionTime = ((double)(end - start)) / CLOCKS_PER_SEC;
    
    // Free memory
    freeMinHeap(minHeap);
    free(dist);
    free(parent);
}

// Dijkstra's algorithm using adjacency matrix
void dijkstraMatrix(Graph_Matrix* graph, int src, int dest, char baseChar, double* executionTime) {
    clock_t start, end;
    start = clock();
    
    int V = graph->numVertices;
    int* dist = (int*)malloc(V * sizeof(int));
    int* parent = (int*)malloc(V * sizeof(int));
    int* visited = (int*)malloc(V * sizeof(int));
    
    // Initialize distances and visited array
    for (int v = 0; v < V; v++) {
        dist[v] = INF;
        parent[v] = -1;
        visited[v] = 0;
    }
    
    dist[src] = 0;
    
    // Find shortest path for all vertices
    for (int count = 0; count < V - 1; count++) {
        // Find minimum distance vertex not yet processed
        int minDist = INF, minIndex = -1;
        for (int v = 0; v < V; v++) {
            if (!visited[v] && dist[v] < minDist) {
                minDist = dist[v];
                minIndex = v;
            }
        }
        
        // If no path exists to remaining vertices
        if (minIndex == -1) break;
        
        // Mark the vertex as processed
        visited[minIndex] = 1;
        
        // Update distances of adjacent vertices
        for (int v = 0; v < V; v++) {
            if (!visited[v] && graph->adjMatrix[minIndex][v] && 
                dist[minIndex] != INF && 
                dist[minIndex] + graph->adjMatrix[minIndex][v] < dist[v]) {
                dist[v] = dist[minIndex] + graph->adjMatrix[minIndex][v];
                parent[v] = minIndex;
            }
        }
    }
    
    // Print result
    if (dist[dest] != INF) {
        printf("\nShortest path from %c to %c using Adjacency Matrix: ", src + baseChar, dest + baseChar);
        printf("%c ", src + baseChar);
        printPath(parent, dest, baseChar);
        printf("\nTotal distance: %d km\n", dist[dest]);
    } else {
        printf("\nNo path exists from %c to %c using Adjacency Matrix\n", src + baseChar, dest + baseChar);
    }
    
    end = clock();
    *executionTime = ((double)(end - start)) / CLOCKS_PER_SEC;
    
    // Free memory
    free(dist);
    free(parent);
    free(visited);
}

// Bellman-Ford algorithm using adjacency list
void bellmanFordList(Graph_List* graph, int src, int dest, char baseChar, double* executionTime) {
    clock_t start, end;
    start = clock();
    
    int V = graph->numVertices;
    int* dist = (int*)malloc(V * sizeof(int));
    int* parent = (int*)malloc(V * sizeof(int));
    
    // Initialize distances
    for (int i = 0; i < V; i++) {
        dist[i] = INF;
        parent[i] = -1;
    }
    
    dist[src] = 0;
    
    // Relax all edges V-1 times
    for (int i = 1; i < V; i++) {
        for (int u = 0; u < V; u++) {
            Road* road = graph->adjLists[u];
            while (road) {
                int v = road->destination;
                int weight = road->distance;
                
                if (dist[u] != INF && dist[u] + weight < dist[v]) {
                    dist[v] = dist[u] + weight;
                    parent[v] = u;
                }
                
                road = road->next;
            }
        }
    }
    
    // Check for negative weight cycles
    for (int u = 0; u < V; u++) {
        Road* road = graph->adjLists[u];
        while (road) {
            int v = road->destination;
            int weight = road->distance;
            
            if (dist[u] != INF && dist[u] + weight < dist[v]) {
                printf("Graph contains negative weight cycle\n");
                free(dist);
                free(parent);
                *executionTime = 0;
                return;
            }
            
            road = road->next;
        }
    }
    
    // Print result
    if (dist[dest] != INF) {
        printf("\nShortest path from %c to %c using Bellman-Ford: ", src + baseChar, dest + baseChar);
        printf("%c ", src + baseChar);
        printPath(parent, dest, baseChar);
        printf("\nTotal distance: %d km\n", dist[dest]);
    } else {
        printf("\nNo path exists from %c to %c using Bellman-Ford\n", src + baseChar, dest + baseChar);
    }
    
    end = clock();
    *executionTime = ((double)(end - start)) / CLOCKS_PER_SEC;
    
    // Free memory
    free(dist);
    free(parent);
}

// ---- File and data loading functions ----

// Load graph from file into both implementations
void loadGraphFromFile(Graph_List* graphList, Graph_Matrix* graphMatrix, const char* filename, char baseChar) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        printf("Could not open file %s for reading\n", filename);
        return;
    }

    char line[256];
    int maxVertex = -1;
    
    // First pass: find the maximum vertex
    while (fgets(line, sizeof(line), file)) {
        char src, dest;
        int weight;
        
        if (sscanf(line, "%c,%c,%dkm", &src, &dest, &weight) == 3) {
            int srcIndex = src - baseChar;
            int destIndex = dest - baseChar;
            
            if (srcIndex > maxVertex) maxVertex = srcIndex;
            if (destIndex > maxVertex) maxVertex = destIndex;
        }
    }
    
    // Ensure both graphs have enough vertices
    graphList->numVertices = maxVertex + 1;
    graphMatrix->numVertices = maxVertex + 1;
    
    // Reset file position to beginning
    rewind(file);
    
    // Second pass: add edges
    while (fgets(line, sizeof(line), file)) {
        char src, dest;
        int weight;
        
        if (sscanf(line, "%c,%c,%dkm", &src, &dest, &weight) == 3) {
            int srcIndex = src - baseChar;
            int destIndex = dest - baseChar;
            
            addEdgeList(graphList, srcIndex, destIndex, weight);
            addEdgeMatrix(graphMatrix, srcIndex, destIndex, weight);
        }
    }
    
    fclose(file);
}

// ---- Testing and performance evaluation functions ----

// Generate a random graph with a specified number of vertices and edges
void generateRandomGraph(Graph_List* graphList, Graph_Matrix* graphMatrix, int V, int E) {
    srand(time(NULL));
    
    for (int i = 0; i < E; i++) {
        int src = rand() % V;
        int dest = rand() % V;
        int weight = rand() % 100 + 1; // Random weight between 1 and 100
        
        if (src != dest) {
            addEdgeList(graphList, src, dest, weight);
            addEdgeMatrix(graphMatrix, src, dest, weight);
        }
    }
}

// Print the adjacency list representation of the graph
void printGraphList(Graph_List* graph, char baseChar) {
    printf("\n--- Graph Adjacency List Representation ---\n");
    for (int v = 0; v < graph->numVertices; v++) {
        Road* road = graph->adjLists[v];
        if (road) {
            printf("Vertex %c: ", v + baseChar);
            while (road) {
                printf("-> %c (%d km) ", road->destination + baseChar, road->distance);
                road = road->next;
            }
            printf("\n");
        }
    }
}

// Print the adjacency matrix representation of the graph
void printGraphMatrix(Graph_Matrix* graph, char baseChar) {
    printf("\n--- Graph Adjacency Matrix Representation ---\n");
    printf("   ");
    for (int v = 0; v < graph->numVertices; v++) {
        printf("%c  ", v + baseChar);
    }
    printf("\n");
    
    for (int i = 0; i < graph->numVertices; i++) {
        printf("%c: ", i + baseChar);
        for (int j = 0; j < graph->numVertices; j++) {
            if (graph->adjMatrix[i][j]) {
                printf("%-2d ", graph->adjMatrix[i][j]);
            } else {
                printf("0  ");
            }
        }
        printf("\n");
    }
}

// Test correctness of the graph implementations
void testGraphImplementations(Graph_List* graphList, Graph_Matrix* graphMatrix, char baseChar) {
    printf("\n--- Testing Graph Implementations ---\n");
    
    // Check if both graphs have the same number of vertices
    if (graphList->numVertices != graphMatrix->numVertices) {
        printf("ERROR: Different number of vertices in the two implementations\n");
        return;
    }
    
    int V = graphList->numVertices;
    int errors = 0;
    
    // Check if both graphs have the same edges with the same weights
    for (int i = 0; i < V; i++) {
        Road* road = graphList->adjLists[i];
        while (road) {
            int j = road->destination;
            int weight = road->distance;
            
            if (graphMatrix->adjMatrix[i][j] != weight) {
                printf("ERROR: Edge (%c,%c) has weight %d in list but %d in matrix\n",
                       i + baseChar, j + baseChar, weight, graphMatrix->adjMatrix[i][j]);
                errors++;
            }
            
            road = road->next;
        }
    }
    
    if (errors == 0) {
        printf("Both graph implementations are consistent\n");
    } else {
        printf("Found %d inconsistencies between the two implementations\n", errors);
    }
}

// Perform performance comparison between the two implementations
void performanceComparison(Graph_List* graphList, Graph_Matrix* graphMatrix, char baseChar) {
    printf("\n--- Performance Comparison ---\n");
    
    int numVertices = graphList->numVertices;
    
    // Test for different source-destination pairs
    int numTests = 5;
    double totalListTime = 0.0;
    double totalMatrixTime = 0.0;
    double totalBellmanFordTime = 0.0;
    
    printf("\nRunning %d tests with random source-destination pairs...\n", numTests);
    printf("--------------------------------------------------------------\n");
    printf("| Test | Source | Dest | Dijkstra (List) | Dijkstra (Matrix) | Bellman-Ford |\n");
    printf("--------------------------------------------------------------\n");
    
    for (int i = 0; i < numTests; i++) {
        int src = rand() % numVertices;
        int dest = rand() % numVertices;
        
        double listTime, matrixTime, bellmanFordTime;
        
        // Run both algorithms
        dijkstraList(graphList, src, dest, baseChar, &listTime);
        dijkstraMatrix(graphMatrix, src, dest, baseChar, &matrixTime);
        bellmanFordList(graphList, src, dest, baseChar, &bellmanFordTime);
        
        // Update totals
        totalListTime += listTime;
        totalMatrixTime += matrixTime;
        totalBellmanFordTime += bellmanFordTime;
        
        // Print individual test results
        printf("| %4d | %6c | %4c | %14.6f | %16.6f | %12.6f |\n",
               i+1, src + baseChar, dest + baseChar, listTime, matrixTime, bellmanFordTime);
    }
    
    printf("--------------------------------------------------------------\n");
    
    // Print averages
    printf("\nAverage execution times over %d tests:\n", numTests);
    printf("Dijkstra (Adjacency List): %.6f seconds\n", totalListTime / numTests);
    printf("Dijkstra (Adjacency Matrix): %.6f seconds\n", totalMatrixTime / numTests);
    printf("Bellman-Ford (Adjacency List): %.6f seconds\n", totalBellmanFordTime / numTests);
    
    // Compare and provide analysis
    printf("\nPerformance Analysis:\n");
    if (totalListTime < totalMatrixTime) {
        printf("- Dijkstra's algorithm using adjacency lists is %.2f%% faster than using adjacency matrix\n",
               (totalMatrixTime - totalListTime) / totalMatrixTime * 100);
    } else {
        printf("- Dijkstra's algorithm using adjacency matrix is %.2f%% faster than using adjacency lists\n",
               (totalListTime - totalMatrixTime) / totalListTime * 100);
    }
    
    if (totalListTime < totalBellmanFordTime) {
        printf("- Dijkstra's algorithm (list) is %.2f%% faster than Bellman-Ford algorithm\n",
               (totalBellmanFordTime - totalListTime) / totalBellmanFordTime * 100);
    } else {
        printf("- Bellman-Ford algorithm is %.2f%% faster than Dijkstra's algorithm (list)\n",
               (totalListTime - totalBellmanFordTime) / totalListTime * 100);
    }
}

int main() {
    Graph_List* graphList = createGraphList(MAX_CITIES);
    Graph_Matrix* graphMatrix = createGraphMatrix(MAX_CITIES);
    char baseChar = 'A';  // Interprets cities as A-Z

    // Load from roads.txt file (ensure it's in the same folder)
    loadGraphFromFile(graphList, graphMatrix, "roads.txt", baseChar);

    // Display graphs
    printGraphList(graphList, baseChar);
    printGraphMatrix(graphMatrix, baseChar);

    // Check consistency
    testGraphImplementations(graphList, graphMatrix, baseChar);

    // Compare performance of all 3 algorithms
    performanceComparison(graphList, graphMatrix, baseChar);

    // Cleanup
    freeGraphList(graphList);
    freeGraphMatrix(graphMatrix);

    return 0;
}    