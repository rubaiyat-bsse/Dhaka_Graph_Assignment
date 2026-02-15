#include<bits/stdc++.h>
using namespace std;

struct Node {
    int id;
    double longitude;
    double latitude;
    string type;  // "road" or "metro"
    string name;  // station name for metro, empty for road
};

struct Graph {
    int V;
    int E;
    vector<Node> nodes;
    map<pair<double, double>, int> coordToId;
    unordered_map<int, vector<pair<int, double>>> adj; // weight = COST (not distance)
};

// Global arrays for Dijkstra
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

void readRoadCSV(string filename, Graph& G) {
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
        
        getline(ss, cellValue, ','); // Skip "DhakaStreet"
        
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
                G.nodes.push_back({newId, coord.first, coord.second, "road", ""});
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
            double cost = distance * 20.0; // ৳20/km for car
            
            G.adj[node1].push_back({node2, cost});
            G.adj[node2].push_back({node1, cost});
        }
    }
    
    file.close();
}

void readMetroCSV(string filename, Graph& G) {
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
        
        // Skip first column "DhakaMetroRail"
        getline(ss, cellValue, ',');
        
        // Read all remaining fields
        string remainder;
        getline(ss, remainder);
        
        // Find the last two commas to extract station names
        size_t lastComma = remainder.rfind(',');
        if(lastComma == string::npos) continue;
        
        size_t secondLastComma = remainder.rfind(',', lastComma - 1);
        if(secondLastComma == string::npos) continue;
        
        string endStation = remainder.substr(lastComma + 1);
        string startStation = remainder.substr(secondLastComma + 1, lastComma - secondLastComma - 1);
        string coordsStr = remainder.substr(0, secondLastComma);
        
        // Parse coordinates
        stringstream coordStream(coordsStr);
        while(getline(coordStream, cellValue, ',')){
            values.push_back(stod(cellValue));
        }
        
        if(values.size() < 4) continue;
        
        // Extract coordinate pairs
        vector<pair<double, double>> metroPoints;
        for(int i = 0; i < values.size(); i += 2){
            if(i + 1 < values.size()){
                double lon = values[i];
                double lat = values[i + 1];
                metroPoints.push_back({lon, lat});
            }
        }
        
        if(metroPoints.size() < 2) continue;
        
        vector<int> nodeIds;
        
        // First and last points are stations
        for(int idx = 0; idx < metroPoints.size(); idx++){
            auto& coord = metroPoints[idx];
            
            // Always create new metro node (even if coordinate exists as road node)
            int newId = G.nodes.size();
            
            string stationName = "";
            if(idx == 0) stationName = startStation;
            else if(idx == metroPoints.size() - 1) stationName = endStation;
            
            G.nodes.push_back({newId, coord.first, coord.second, "metro", stationName});
            nodeIds.push_back(newId);
        }
        
        // Create edges between consecutive metro points
        for (int i = 0; i < nodeIds.size() - 1; i++){
            int node1 = nodeIds[i];
            int node2 = nodeIds[i + 1];
            
            double lon1 = G.nodes[node1].longitude;
            double lat1 = G.nodes[node1].latitude;
            double lon2 = G.nodes[node2].longitude;
            double lat2 = G.nodes[node2].latitude;
            
            double distance = haversine(lon1, lat1, lon2, lat2);
            double cost = distance * 5.0; // ৳5/km for metro
            
            G.adj[node1].push_back({node2, cost});
            G.adj[node2].push_back({node1, cost});
        }
    }
    
    file.close();
}

void addTransferEdges(Graph& G, double threshold = 0.5) {
    // For each metro node, find nearest road node within threshold
    for(int i = 0; i < G.nodes.size(); i++){
        if(G.nodes[i].type != "metro") continue;
        
        // Only create transfer edges at STATIONS (nodes with names)
        // Track points (no names) cannot be used for transfers
        if(G.nodes[i].name.empty()) continue;
        
        double metroLon = G.nodes[i].longitude;
        double metroLat = G.nodes[i].latitude;
        
        int nearestRoad = -1;
        double minDist = DBL_MAX;
        
        for(int j = 0; j < G.nodes.size(); j++){
            if(G.nodes[j].type != "road") continue;
            
            double roadLon = G.nodes[j].longitude;
            double roadLat = G.nodes[j].latitude;
            
            double d = haversine(metroLon, metroLat, roadLon, roadLat);
            
            if(d < minDist && d <= threshold){
                minDist = d;
                nearestRoad = j;
            }
        }
        
        if(nearestRoad != -1){
            // Add bidirectional edge with cost = 0 (free walking)
            G.adj[i].push_back({nearestRoad, 0.0});
            G.adj[nearestRoad].push_back({i, 0.0});
        }
    }
}

void initGlobalArray(int V){
    dist.resize(V);
    parent.resize(V);
}

void INITIALIZE_SINGLE_SOURCE(Graph& G, int src){
    K.clear();
    while(!Q.empty()) Q.pop();
    
    for(int i = 0; i < G.V; i++){
        dist[i] = DBL_MAX;
        parent[i] = -1;
    }
    dist[src] = 0;
}

void RELAX(int u, int v, double w){
    if(dist[v] > dist[u] + w){
        dist[v] = dist[u] + w;
        parent[v] = u;
    }
}

void Dijkstra(Graph& G, int src){
    INITIALIZE_SINGLE_SOURCE(G, src);
    Q.push({0, src});
    
    while(!Q.empty()){
        int u = Q.top().second;
        Q.pop();
        
        if(K.find(u) != K.end()) continue;
        K.insert(u);
        
        for(auto& edge : G.adj[u]){
            int v = edge.first;
            double w = edge.second;
            
            double oldDist = dist[v];
            RELAX(u, v, w);
            
            if(dist[v] < oldDist){
                Q.push({dist[v], v});
            }
        }
    }
}

void printPath(Graph& G, int v, vector<int>& path){
    if(v == -1) return;
    printPath(G, parent[v], path);
    path.push_back(v);
}

void generateKML(Graph& G, vector<int>& path, string filename) {
    ofstream kml(filename);
    
    kml << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    kml << "<kml xmlns=\"http://www.opengis.net/kml/2.2\">\n";
    kml << "<Document>\n";
    kml << "  <name>Cheapest Path - Dhaka (Car + Metro)</name>\n";
    
    // Road style (red)
    kml << "  <Style id=\"roadStyle\">\n";
    kml << "    <LineStyle>\n";
    kml << "      <color>ff0000ff</color>\n";
    kml << "      <width>4</width>\n";
    kml << "    </LineStyle>\n";
    kml << "  </Style>\n";
    
    // Metro style (blue)
    kml << "  <Style id=\"metroStyle\">\n";
    kml << "    <LineStyle>\n";
    kml << "      <color>ffff0000</color>\n";
    kml << "      <width>4</width>\n";
    kml << "    </LineStyle>\n";
    kml << "  </Style>\n";
    
    // Source marker
    kml << "  <Placemark>\n";
    kml << "    <name>Source</name>\n";
    kml << "    <Point>\n";
    kml << "      <coordinates>" << G.nodes[path[0]].longitude << "," 
        << G.nodes[path[0]].latitude << ",0</coordinates>\n";
    kml << "    </Point>\n";
    kml << "  </Placemark>\n";
    
    // Destination marker
    kml << "  <Placemark>\n";
    kml << "    <name>Destination</name>\n";
    kml << "    <Point>\n";
    kml << "      <coordinates>" << G.nodes[path.back()].longitude << "," 
        << G.nodes[path.back()].latitude << ",0</coordinates>\n";
    kml << "    </Point>\n";
    kml << "  </Placemark>\n";
    
    // Draw path segments by mode
    for(int i = 0; i < path.size() - 1; i++){
        string style = (G.nodes[path[i]].type == "metro" || G.nodes[path[i+1]].type == "metro") 
                       ? "metroStyle" : "roadStyle";
        
        kml << "  <Placemark>\n";
        kml << "    <styleUrl>#" << style << "</styleUrl>\n";
        kml << "    <LineString>\n";
        kml << "      <coordinates>\n";
        kml << "        " << G.nodes[path[i]].longitude << "," 
            << G.nodes[path[i]].latitude << ",0\n";
        kml << "        " << G.nodes[path[i+1]].longitude << "," 
            << G.nodes[path[i+1]].latitude << ",0\n";
        kml << "      </coordinates>\n";
        kml << "    </LineString>\n";
    kml << "  </Placemark>\n";
    }
    
    kml << "</Document>\n";
    kml << "</kml>\n";
    
    kml.close();
}

int main(){
    Graph G;
    
    cout << "Reading road network..." << endl;
    readRoadCSV("../Dataset/Roadmap-Dhaka.csv", G);
    
    cout << "Reading metro rail network..." << endl;
    readMetroCSV("../Dataset/Routemap-DhakaMetroRail.csv", G);
    
    cout << "Adding transfer edges..." << endl;
    addTransferEdges(G, 0.5); // 500m threshold
    
    G.V = G.nodes.size();
    G.E = 0;
    for(auto& entry : G.adj){
        G.E += entry.second.size();
    }
    G.E /= 2;
    
    cout << "Graph Statistics:" << endl;
    cout << "Total nodes: " << G.V << endl;
    cout << "Total edges: " << G.E << endl;
    
    int roadNodes = 0, metroNodes = 0;
    for(auto& node : G.nodes){
        if(node.type == "road") roadNodes++;
        else metroNodes++;
    }
    cout << "Road nodes: " << roadNodes << endl;
    cout << "Metro nodes: " << metroNodes << endl;
    cout << endl;
    
    // Read source and destination
    double srcLon, srcLat, destLon, destLat;
    cin >> srcLon >> srcLat >> destLon >> destLat;
    
    // Find source and destination nodes
    auto srcCoord = make_pair(srcLon, srcLat);
    auto destCoord = make_pair(destLon, destLat);
    
    if(G.coordToId.find(srcCoord) == G.coordToId.end()){
        cerr << "Error: Source coordinate not found in graph!" << endl;
        return 1;
    }
    
    if(G.coordToId.find(destCoord) == G.coordToId.end()){
        cerr << "Error: Destination coordinate not found in graph!" << endl;
        return 1;
    }
    
    int srcNode = G.coordToId[srcCoord];
    int destNode = G.coordToId[destCoord];
    
    cout << "Source: Node " << srcNode << " (" << srcLon << ", " << srcLat << ")" << endl;
    cout << "Destination: Node " << destNode << " (" << destLon << ", " << destLat << ")" << endl;
    cout << endl;
    
    // Initialize and run Dijkstra
    initGlobalArray(G.V);
    Dijkstra(G, srcNode);
    
    if(dist[destNode] == DBL_MAX){
        cout << "No path found!" << endl;
        return 1;
    }
    
    // Reconstruct path
    vector<int> path;
    printPath(G, destNode, path);
    
    cout << "Minimum cost: ৳" << fixed << setprecision(2) << dist[destNode] << endl;
    cout << endl;
    
    // Analyze path segments
    double totalCarCost = 0, totalMetroCost = 0;
    double carDist = 0, metroDist = 0;
    
    cout << "Path details: " << path.size() << " nodes" << endl;
    for(int i = 0; i < path.size(); i++){
        int nodeId = path[i];
        cout << "Node " << nodeId << " (" << G.nodes[nodeId].type << ")";
        
        if(G.nodes[nodeId].type == "metro" && !G.nodes[nodeId].name.empty()){
            cout << " [" << G.nodes[nodeId].name << "]";
        }
        
        if(i < path.size() - 1){
            int nextNode = path[i + 1];
            double lon1 = G.nodes[nodeId].longitude;
            double lat1 = G.nodes[nodeId].latitude;
            double lon2 = G.nodes[nextNode].longitude;
            double lat2 = G.nodes[nextNode].latitude;
            double segDist = haversine(lon1, lat1, lon2, lat2);
            
            // Determine segment type
            if(G.nodes[nodeId].type == "metro" || G.nodes[nextNode].type == "metro"){
                double segCost = segDist * 5.0;
                totalMetroCost += segCost;
                metroDist += segDist;
                cout << " -> " << fixed << setprecision(3) << segDist << " km (metro, ৳" 
                     << setprecision(2) << segCost << ")";
            } else {
                double segCost = segDist * 20.0;
                totalCarCost += segCost;
                carDist += segDist;
                cout << " -> " << fixed << setprecision(3) << segDist << " km (car, ৳" 
                     << setprecision(2) << segCost << ")";
            }
        }
        cout << endl;
    }
    
    cout << endl;
    cout << "Summary:" << endl;
    cout << "Car distance: " << fixed << setprecision(3) << carDist << " km (Cost: ৳" 
         << setprecision(2) << totalCarCost << ")" << endl;
    cout << "Metro distance: " << fixed << setprecision(3) << metroDist << " km (Cost: ৳" 
         << setprecision(2) << totalMetroCost << ")" << endl;
    cout << "Total cost: ৳" << fixed << setprecision(2) << (totalCarCost + totalMetroCost) << endl;
    
    // Generate KML
    generateKML(G, path, "route_caseA.kml");
    cout << "\nKML file generated: route_caseA.kml" << endl;
    
    return 0;
}
