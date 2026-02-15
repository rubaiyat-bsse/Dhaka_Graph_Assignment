#include<bits/stdc++.h>
using namespace std;

struct Node {
    int id;
    double longitude;
    double latitude;
};

struct Graph {
    int V;
    int E;
    vector<Node> nodes;
    map<pair<double, double>, int> coordToId;
    /*
     while csv parsing we need to know if the parsed node already 
     exists or not. thats why we need a <lon, lat>->id map
     we used red-black tree map because c++ doesn't come with a hash
     function for pair<double, double> -> int
    */
    unordered_map<int, vector<pair<int, double>>> adj; 
    /*
     adjacency list: node -> [(neighbor, distance)]
     for a weighted graph, adjacency list is normally vector<vector<<pair<int,int>>>
     here for efficiency we used hashmap. 
    */
};

//global arrays for dijkstra
vector<double> dist;
vector<int> parent;
set<int> K;

priority_queue<
    pair<double, int>,
    vector<pair<double, int>>,
    greater<pair<double, int>>
> Q;

double haversine(double lon1, double lat1, double lon2, double lat2) {
    const double R = 6371.0;
    
    double dLat = (lat2 - lat1) * M_PI / 180.0;
    double dLon = (lon2 - lon1) * M_PI / 180.0;
    
    double a = sin(dLat / 2) * sin(dLat / 2) +
               cos(lat1 * M_PI / 180.0) * cos(lat2 * M_PI / 180.0) *
               sin(dLon / 2) * sin(dLon / 2);
    
    double c = 2 * atan2(sqrt(a), sqrt(1 - a));
    return R * c;
}

void readCSV(string filename, Graph& G) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Could not open file " << filename << endl;
        return;
    }
    
    string line;
    
    while(getline(file, line)){
        stringstream ss(line);
        string cellValue;
        vector<double> values;
        
        getline(ss, cellValue, ',');
        
        while(getline(ss, cellValue, ',')){
            values.push_back(stod(cellValue));
        }
        
        if(values.size() < 4) continue;
        
        vector<pair<double, double>> roadPoints;
        for(int i = 0; i < values.size() - 2; i += 2){
            double lon = values[i];
            double lat = values[i + 1];
            roadPoints.push_back({lon, lat});
        }
        
        if(roadPoints.size() < 2) continue;
        
        vector<int> nodeIds;
        
        for(auto& coord : roadPoints){
            if (G.coordToId.find(coord) == G.coordToId.end()){
                int newId = G.nodes.size();
                G.coordToId[coord] = newId;
                G.nodes.push_back({newId, coord.first, coord.second});
            }
            nodeIds.push_back(G.coordToId[coord]);
        }
        
        for (int i = 0; i < nodeIds.size() - 1; i++){
            int node1 = nodeIds[i];
            int node2 = nodeIds[i + 1];
            
            double lon1 = G.nodes[node1].longitude;
            double lat1 = G.nodes[node1].latitude;
            double lon2 = G.nodes[node2].longitude;
            double lat2 = G.nodes[node2].latitude;
            
            double distance = haversine(lon1, lat1, lon2, lat2);
            
            G.adj[node1].push_back({node2, distance});
            G.adj[node2].push_back({node1, distance});
        }
    }
    
    file.close();
    
    G.V = G.nodes.size();
    G.E = 0;
    for (auto& entry : G.adj) {
        G.E += entry.second.size();
    }
    G.E /= 2;
}

void initGlobalArray(int V) {
    dist.assign(V, DBL_MAX);
    parent.assign(V, -1);
}

void INITIALIZE_SINGLE_SOURCE(Graph& G, int src) {
    initGlobalArray(G.V);
    dist[src] = 0;
}

void RELAX(int u, int v, double w) {
    if(dist[v] > dist[u] + w) {
        dist[v] = dist[u] + w;
        parent[v] = u;
        Q.push({dist[v], v});
    }
}

void Dijkstra(Graph& G, int src) {
    K.clear();
    while(!Q.empty()) Q.pop();
    
    INITIALIZE_SINGLE_SOURCE(G, src);
    
    for(int v = 0; v < G.V; v++) {
        Q.push({dist[v], v});
    }
    
    while(!Q.empty()) {
        int u = Q.top().second;
        Q.pop();
        
        if(K.count(u)) continue;
        
        K.insert(u);
        
        for(auto edge : G.adj[u]) {
            int v = edge.first;
            double w = edge.second;
            
            if(!K.count(v)) {
                RELAX(u, v, w);
            }
        }
    }
}

void printPath(int v) {
    if(v == -1) return;
    printPath(parent[v]);
    cout << v << " ";
}

vector<int> getPath(int dest) {
    vector<int> path;
    for(int v = dest; v != -1; v = parent[v]) {
        path.push_back(v);
    }
    reverse(path.begin(), path.end());
    return path;
}

void printDetailedPath(Graph& G, vector<int>& path) {
    cout << "Path: " << path.size() << " nodes" << endl;
    for(int i = 0; i < path.size(); i++) {
        int nodeId = path[i];
        cout << "Node " << nodeId << ": (" 
             << fixed << setprecision(6)
             << G.nodes[nodeId].longitude << ", "
             << G.nodes[nodeId].latitude << ")";
        
        if(i < path.size() - 1) {
            int nextNode = path[i + 1];
            for(auto& [neighbor, d] : G.adj[nodeId]) {
                if(neighbor == nextNode) {
                    cout << " -> " << fixed << setprecision(3) << d << " km";
                    break;
                }
            }
        }
        cout << endl;
    }
}

void generateKML(Graph& G, vector<int>& path, string filename) {
    ofstream kml(filename);
    
    kml << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    kml << "<kml xmlns=\"http://www.opengis.net/kml/2.2\">\n";
    kml << "  <Document>\n";
    kml << "    <name>Shortest Path - Dhaka</name>\n";
    kml << "    <Style id=\"routeStyle\">\n";
    kml << "      <LineStyle>\n";
    kml << "        <color>ff0000ff</color>\n";
    kml << "        <width>4</width>\n";
    kml << "      </LineStyle>\n";
    kml << "    </Style>\n";
    
    // source marker
    kml << "    <Placemark>\n";
    kml << "      <name>Source</name>\n";
    kml << "      <Point>\n";
    kml << "        <coordinates>" 
        << fixed << setprecision(6)
        << G.nodes[path[0]].longitude << "," 
        << G.nodes[path[0]].latitude << ",0</coordinates>\n";
    kml << "      </Point>\n";
    kml << "    </Placemark>\n";
    
    // destination marker
    kml << "    <Placemark>\n";
    kml << "      <name>Destination</name>\n";
    kml << "      <Point>\n";
    kml << "        <coordinates>" 
        << fixed << setprecision(6)
        << G.nodes[path.back()].longitude << "," 
        << G.nodes[path.back()].latitude << ",0</coordinates>\n";
    kml << "      </Point>\n";
    kml << "    </Placemark>\n";
    
    // route path
    kml << "    <Placemark>\n";
    kml << "      <name>Route</name>\n";
    kml << "      <styleUrl>#routeStyle</styleUrl>\n";
    kml << "      <LineString>\n";
    kml << "        <coordinates>\n";
    
    for(int nodeId : path) {
        kml << "          " 
            << fixed << setprecision(6)
            << G.nodes[nodeId].longitude << "," 
            << G.nodes[nodeId].latitude << ",0\n";
    }
    
    kml << "        </coordinates>\n";
    kml << "      </LineString>\n";
    kml << "    </Placemark>\n";
    kml << "  </Document>\n";
    kml << "</kml>\n";
    
    kml.close();
}

int main() {
    Graph G;
    
    string filename = "../Dataset/Roadmap-Dhaka.csv";
    readCSV(filename, G);
    
    double src_lon, src_lat, dest_lon, dest_lat;
    cin >> src_lon >> src_lat >> dest_lon >> dest_lat;
    
    // case A: both source and destination are existing nodes
    auto srcCoord = make_pair(src_lon, src_lat);
    auto destCoord = make_pair(dest_lon, dest_lat);
    
    if(G.coordToId.find(srcCoord) == G.coordToId.end()) {
        cerr << "Error: source coordinate not found in graph" << endl;
        return 1;
    }
    
    if(G.coordToId.find(destCoord) == G.coordToId.end()) {
        cerr << "Error: destination coordinate not found in graph" << endl;
        return 1;
    }
    
    int srcNode = G.coordToId[srcCoord];
    int destNode = G.coordToId[destCoord];
    
    Dijkstra(G, srcNode);
    
    if(dist[destNode] == DBL_MAX) {
        cout << "No path exists" << endl;
        return 1;
    }
    
    cout << "Shortest distance: " << fixed << setprecision(3) << dist[destNode] << " km" << endl;
    cout << endl;
    
    vector<int> path = getPath(destNode);
    printDetailedPath(G, path);
    
    generateKML(G, path, "route_caseA.kml");
    cout << "\nKML file generated: route_caseA.kml" << endl;
    
    return 0;
}
