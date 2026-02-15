#include <bits/stdc++.h>
using namespace std;

typedef struct{
    int V;
    vector<vector<pair<int,int>>> adj; 
} Graph;

vector<int> dist, parent;
set<int> K;

// min-priority queue: (distance, vertex)
priority_queue<
    pair<int,int>,
    vector<pair<int,int>>,
    greater<pair<int,int>>
> Q;

void initGraph(Graph &G, int V){
    G.V = V;
    G.adj.assign(V, vector<pair<int,int>>());
}

void initGlobalArray(int V){
    dist.assign(V, INT_MAX);
    parent.assign(V, -1);
}

void INITIALIZE_SINGLE_SOURCE(Graph &G, int src){
    initGlobalArray(G.V);
    dist[src] = 0;
}

void RELAX(int u, int v, int w){
    if(dist[v] > dist[u] + w){
        dist[v] = dist[u] + w;
        parent[v] = u;
        Q.push({dist[v], v});
    }
}

void Dijkstra(Graph &G, int src){
    INITIALIZE_SINGLE_SOURCE(G, src);

    // Q ← V
    for(int v = 0; v < G.V; v++){
        Q.push({dist[v], v});
    }

    while(!Q.empty()){
        // u ← EXTRACT-MIN(Q)
        int u = Q.top().second;
        Q.pop();

        if(K.count(u)) continue;

        // K ← K ∪ {u}
        K.insert(u);

        // for each v ∈ Adj[u]
        for(auto edge : G.adj[u]){
            int v = edge.first;
            int w = edge.second;

            if(!K.count(v)){
                RELAX(u, v, w);
            }
        }
    }
}

void printPath(int v){
    if(v == -1) return;      
    printPath(parent[v]); 
    cout << v << " ";
}


int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n, m;
    cin >> n >> m;

    Graph G;
    initGraph(G, n);

    for(int i = 0; i < m; i++){
        int u, v, w;
        cin >> u >> v >> w;
        G.adj[u].push_back({v, w});
    }

    int src;
    cin >> src;

    Dijkstra(G, src);

    // output distances
    for(int i = 0; i < n; i++){
        if(dist[i] == INT_MAX)
            cout << "INF\n";
        else
            cout << dist[i] << "\n";
    }

    return 0;
}
