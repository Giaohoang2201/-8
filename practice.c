#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>

#define V 50
#define SPARSE_E 100
#define DENSE_E 600
#define INF 1e9

// 간선 구조체 (크루스칼 정렬용)
typedef struct {
    int u, v, weight;
} Edge;

// 인접 리스트 노드 구조체 (크루스칼 실험용)
typedef struct AdjNode {
    int dest;
    int weight;
    struct AdjNode* next;
} AdjNode;

// 유니온-파인드 (Union-Find) 구조체
typedef struct {
    int parent;
    int rank;
} UFNode;

int find(UFNode uf[], int i) {
    if (uf[i].parent != i)
        uf[i].parent = find(uf, uf[i].parent);
    return uf[i].parent;
}

void unionSets(UFNode uf[], int x, int y) {
    int rootX = find(uf, x);
    int rootY = find(uf, y);
    if (rootX != rootY) {
        if (uf[rootX].rank < uf[rootY].rank)
            uf[rootX].parent = rootY;
        else if (uf[rootX].rank > uf[rootY].rank)
            uf[rootY].parent = rootX;
        else {
            uf[rootY].parent = rootX;
            uf[rootX].rank++;
        }
    }
}

int compareEdges(const void* a, const void* b) {
    Edge* e1 = (Edge*)a;
    Edge* e2 = (Edge*)b;
    return e1->weight - e2->weight;
}

// DFS를 통한 그래프 연결성 확인 함수 (인접 행렬 기반)
void dfs(int u, int adj_matrix[V][V], bool visited[V]) {
    visited[u] = true;
    for (int v = 0; v < V; v++) {
        if (adj_matrix[u][v] != INF && !visited[v]) {
            dfs(v, adj_matrix, visited);
        }
    }
}

bool isConnected(int adj_matrix[V][V]) {
    bool visited[V];
    for (int i = 0; i < V; i++) visited[i] = false;
    dfs(0, adj_matrix, visited);
    for (int i = 0; i < V; i++) {
        if (!visited[i]) return false;
    }
    return true;
}

// 무작위 연결 그래프 생성 함수
void generateGraph(int E, int adj_matrix[V][V], Edge edges[]) {
    bool connected = false;
    while (!connected) {
        for (int i = 0; i < V; i++) {
            for (int j = 0; j < V; j++) {
                adj_matrix[i][j] = (i == j) ? 0 : INF;
            }
        }

        int edge_count = 0;
        while (edge_count < E) {
            int u = rand() % V;
            int v = rand() % V;
            if (u == v || adj_matrix[u][v] != INF) continue;

            int weight = (rand() % 100) + 1; // 1 ~ 100 가중치
            adj_matrix[u][v] = weight;
            adj_matrix[v][u] = weight;

            edges[edge_count].u = u < v ? u : v;
            edges[edge_count].v = u < v ? v : u;
            edges[edge_count].weight = weight;

            edge_count++;
        }

        if (isConnected(adj_matrix)) {
            connected = true;
        }
    }
}

// 인접 리스트 생성 함수
void buildAdjacencyList(int adj_matrix[V][V], AdjNode* adj_list[V], int E) {
    for (int i = 0; i < V; i++) {
        adj_list[i] = NULL;
    }

    for (int i = 0; i < V; i++) {
        for (int j = i + 1; j < V; j++) {
            if (adj_matrix[i][j] != INF) {
                // i -> j 추가
                AdjNode* newNode1 = (AdjNode*)malloc(sizeof(AdjNode));
                newNode1->dest = j;
                newNode1->weight = adj_matrix[i][j];
                newNode1->next = adj_list[i];
                adj_list[i] = newNode1;

                // j -> i 추가 (무방향 그래프)
                AdjNode* newNode2 = (AdjNode*)malloc(sizeof(AdjNode));
                newNode2->dest = i;
                newNode2->weight = adj_matrix[i][j];
                newNode2->next = adj_list[j];
                adj_list[j] = newNode2;
            }
        }
    }
}

// 인접 리스트 메모리 해제
void freeAdjacencyList(AdjNode* adj_list[V]) {
    for (int i = 0; i < V; i++) {
        AdjNode* curr = adj_list[i];
        while (curr != NULL) {
            AdjNode* temp = curr;
            curr = curr->next;
            free(temp);
        }
    }
}

// [Case 1 & 2] 인접 리스트 + 크루스칼(Kruskal) 알고리즘
long long runKruskal(int E, Edge edges[], int *mst_cost, size_t *memory_used) {
    clock_t start = clock();

    // 간선 배열 복사 및 정렬
    Edge* sorted_edges = (Edge*)malloc(E * sizeof(Edge));
    for (int i = 0; i < E; i++) sorted_edges[i] = edges[i];
    qsort(sorted_edges, E, sizeof(Edge), compareEdges);

    UFNode* uf = (UFNode*)malloc(V * sizeof(UFNode));
    for (int i = 0; i < V; i++) {
        uf[i].parent = i;
        uf[i].rank = 0;
    }

    int mst_edges_count = 0;
    *mst_cost = 0;

    for (int i = 0; i < E && mst_edges_count < V - 1; i++) {
        int u = sorted_edges[i].u;
        int v = sorted_edges[i].v;

        if (find(uf, u) != find(uf, v)) {
            unionSets(uf, u, v);
            *mst_cost += sorted_edges[i].weight;
            mst_edges_count++;
        }
    }

    clock_t end = clock();
    double elapsed_ms = ((double)(end - start) * 1000.0) / CLOCKS_PER_SEC;

    // 메모리 측정: 인접 리스트(헤더 포인터 + 노드) + 간선 배열 + Union-Find
    size_t list_memory = V * sizeof(AdjNode*) + (2 * E) * sizeof(AdjNode);
    size_t aux_memory = E * sizeof(Edge) + V * sizeof(UFNode);
    *memory_used = list_memory + aux_memory;

    free(sorted_edges);
    free(uf);

    return (long long)elapsed_ms;
}

// [Case 3 & 4] 인접 행렬 + 프림(Prim) 알고리즘
long long runPrim(int adj_matrix[V][V], int *mst_cost, size_t *memory_used) {
    clock_t start = clock();

    int parent[V];
    int key[V];
    bool inMST[V];

    for (int i = 0; i < V; i++) {
        key[i] = INF;
        inMST[i] = false;
    }

    key[0] = 0;
    parent[0] = -1;
    *mst_cost = 0;

    for (int count = 0; count < V - 1; count++) {
        int min = INF, u = -1;
        for (int v = 0; v < V; v++) {
            if (!inMST[v] && key[v] < min) {
                min = key[v];
                u = v;
            }
        }

        if (u == -1) break;
        inMST[u] = true;

        for (int v = 0; v < V; v++) {
            if (adj_matrix[u][v] != INF && !inMST[v] && adj_matrix[u][v] < key[v]) {
                parent[v] = u;
                key[v] = adj_matrix[u][v];
            }
        }
    }

    for (int i = 1; i < V; i++) {
        *mst_cost += key[i];
    }

    clock_t end = clock();
    double elapsed_ms = ((double)(end - start) * 1000.0) / CLOCKS_PER_SEC;

    // 메모리 측정: 인접 행렬 + Prim 보조 배열 (parent, key, inMST)
    size_t matrix_memory = V * V * sizeof(int);
    size_t aux_memory = V * sizeof(int) + V * sizeof(int) + V * sizeof(bool);
    *memory_used = matrix_memory + aux_memory;

    return (long long)elapsed_ms;
}

int main() {
    srand((unsigned int)time(NULL));

    // 1. 희소 그래프 생성 (50정점, 100간선)
    int sparse_matrix[V][V];
    Edge sparse_edges[SPARSE_E];
    generateGraph(SPARSE_E, sparse_matrix, sparse_edges);
    
    AdjNode* sparse_list[V];
    buildAdjacencyList(sparse_matrix, sparse_list, SPARSE_E);

    // 2. 밀집 그래프 생성 (50정점, 600간선)
    int dense_matrix[V][V];
    Edge dense_edges[DENSE_E];
    generateGraph(DENSE_E, dense_matrix, dense_edges);

    AdjNode* dense_list[V];
    buildAdjacencyList(dense_matrix, dense_list, DENSE_E);

    int cost1 = 0, cost2 = 0, cost3 = 0, cost4 = 0;
    size_t mem1 = 0, mem2 = 0, mem3 = 0, mem4 = 0;

    // 4가지 실험 케이스 수행
    long long time1 = runKruskal(SPARSE_E, sparse_edges, &cost1, &mem1);
    long long time2 = runKruskal(DENSE_E, dense_edges, &cost2, &mem2);
    long long time3 = runPrim(sparse_matrix, &cost3, &mem3);
    long long time4 = runPrim(dense_matrix, &cost4, &mem4);

    // 지정된 출력 형식에 맞추어 결과 출력
    printf("Case 1 : Sparse Graph - Adjacency List - Kruskal\n");
    printf("Vertices      : %d\n", V);
    printf("Edges         : %d\n", SPARSE_E);
    printf("MST Edges     : %d\n", V - 1);
    printf("MST Cost      : %d\n", cost1);
    printf("Memory        : %zu Bytes\n", mem1);
    printf("Execution Time: %lld ms\n\n", time1);

    printf("Case 2 : Dense Graph - Adjacency List - Kruskal\n");
    printf("Vertices      : %d\n", V);
    printf("Edges         : %d\n", DENSE_E);
    printf("MST Edges     : %d\n", V - 1);
    printf("MST Cost      : %d\n", cost2);
    printf("Memory        : %zu Bytes\n", mem2);
    printf("Execution Time: %lld ms\n\n", time2);

    printf("Case 3 : Sparse Graph - Adjacency Matrix - Prim\n");
    printf("Vertices      : %d\n", V);
    printf("Edges         : %d\n", SPARSE_E);
    printf("MST Edges     : %d\n", V - 1);
    printf("MST Cost      : %d\n", cost3);
    printf("Memory        : %zu Bytes\n", mem3);
    printf("Execution Time: %lld ms\n\n", time3);

    printf("Case 4 : Dense Graph - Adjacency Matrix - Prim\n");
    printf("Vertices      : %d\n", V);
    printf("Edges         : %d\n", DENSE_E);
    printf("MST Edges     : %d\n", V - 1);
    printf("MST Cost      : %d\n", cost4);
    printf("Memory        : %zu Bytes\n", mem4);
    printf("Execution Time: %lld ms\n", time4);

    // 동적 할당 해제
    freeAdjacencyList(sparse_list);
    freeAdjacencyList(dense_list);

    return 0;
}
