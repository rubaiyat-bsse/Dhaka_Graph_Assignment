#include<bits/stdc++.h>
using namespace std;

struct Node {
    int id;
    double longitude;
    double latitude;
    string type;  // "road", "metro", or "bus"
    string name;  // station/stop name, empty for track points
};

struct Graph {
    int V;
    int E;
    vector<Node> nodes;
    map<pair<double, double>, int> coordToId;
    unordered_map<int, vector<pair<int, double>>> adj; // weight = TIME in minutes
    unordered_map<int, unordered_map<int, string>> edgeMode;
};

// Global arrays for Dijkstra
vector<double> dist;       // Minimum time to reach each node
vector<int> parent;
vector<double> arrivalTime; // Arrival time at each node (in minutes from midnight)
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
            double timeMinutes = (distance / 10.0) * 60.0; // 10 km/h -> minutes
            
            G.adj[node1].push_back({node2, timeMinutes});
            G.adj[node2].push_back({node1, timeMinutes});
            
            G.edgeMode[node1][node2] = "road";
            G.edgeMode[node2][node1] = "road";
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
        
        getline(ss, cellValue, ',');
        
        string remainder;
        getline(ss, remainder);
        
        size_t lastComma = remainder.rfind(',');
        if(lastComma == string::npos) continue;
        
        size_t secondLastComma = remainder.rfind(',', lastComma - 1);
        if(secondLastComma == string::npos) continue;
        
        string endStation = remainder.substr(lastComma + 1);
        string startStation = remainder.substr(secondLastComma + 1, lastComma - secondLastComma - 1);
        string coordsStr = remainder.substr(0, secondLastComma);
        
        stringstream coordStream(coordsStr);
        while(getline(coordStream, cellValue, ',')){
            values.push_back(stod(cellValue));
        }
        
        if(values.size() < 4) continue;
        
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
        
        for(int idx = 0; idx < metroPoints.size(); idx++){
            auto& coord = metroPoints[idx];
            
            int newId = G.nodes.size();
            
            string stationName = "";
            if(idx == 0) stationName = startStation;
            else if(idx == metroPoints.size() - 1) stationName = endStation;
            
            G.nodes.push_back({newId, coord.first, coord.second, "metro", stationName});
            nodeIds.push_back(newId);
        }
        
        for (int i = 0; i < nodeIds.size() - 1; i++){
            int node1 = nodeIds[i];
            int node2 = nodeIds[i + 1];
            
            double lon1 = G.nodes[node1].longitude;
            double lat1 = G.nodes[node1].latitude;
            double lon2 = G.nodes[node2].longitude;
            double lat2 = G.nodes[node2].latitude;
            
            double distance = haversine(lon1, lat1, lon2, lat2);
            double timeMinutes = (distance / 10.0) * 60.0; // 10 km/h -> minutes
            
            G.adj[node1].push_back({node2, timeMinutes});
            G.adj[node2].push_back({node1, timeMinutes});
            
            G.edgeMode[node1][node2] = "metro";
            G.edgeMode[node2][node1] = "metro";
        }
    }
    
    file.close();
}

void readBusCSV(string filename, Graph& G) {
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
        
        string remainder;
        getline(ss, remainder);
        
        size_t lastComma = remainder.rfind(',');
        if(lastComma == string::npos) continue;
        
        size_t secondLastComma = remainder.rfind(',', lastComma - 1);
        if(secondLastComma == string::npos) continue;
        
        string endStop = remainder.substr(lastComma + 1);
        string startStop = remainder.substr(secondLastComma + 1, lastComma - secondLastComma - 1);
        string coordsStr = remainder.substr(0, secondLastComma);
        
        stringstream coordStream(coordsStr);
        while(getline(coordStream, cellValue, ',')){
            values.push_back(stod(cellValue));
        }
        
        if(values.size() < 4) continue;
        
        vector<pair<double, double>> busPoints;
        for(int i = 0; i < values.size(); i += 2){
            if(i + 1 < values.size()){
                double lon = values[i];
                double lat = values[i + 1];
                busPoints.push_back({lon, lat});
            }
        }
        
        if(busPoints.size() < 2) continue;
        
        vector<int> nodeIds;
        
        for(int idx = 0; idx < busPoints.size(); idx++){
            auto& coord = busPoints[idx];
            
            int newId = G.nodes.size();
            
            string stopName = "";
            if(idx == 0) stopName = startStop;
            else if(idx == busPoints.size() - 1) stopName = endStop;
            
            G.nodes.push_back({newId, coord.first, coord.second, "bus", stopName});
            nodeIds.push_back(newId);
        }
        
        for (int i = 0; i < nodeIds.size() - 1; i++){
            int node1 = nodeIds[i];
            int node2 = nodeIds[i + 1];
            
            double lon1 = G.nodes[node1].longitude;
            double lat1 = G.nodes[node1].latitude;
            double lon2 = G.nodes[node2].longitude;
            double lat2 = G.nodes[node2].latitude;
            
            double distance = haversine(lon1, lat1, lon2, lat2);
            double timeMinutes = (distance / 10.0) * 60.0; // 10 km/h -> minutes
            
            G.adj[node1].push_back({node2, timeMinutes});
            G.adj[node2].push_back({node1, timeMinutes});
            
            G.edgeMode[node1][node2] = "bus";
            G.edgeMode[node2][node1] = "bus";
        }
    }
    
    file.close();
}

void addTransferEdges(Graph& G, double threshold = 0.5) {
    // Add road ↔ metro transfers
    for(int i = 0; i < G.nodes.size(); i++){
        if(G.nodes[i].type != "metro") continue;
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
            double walkTime = (minDist / 2.0) * 60.0; // 2 km/h walking
            G.adj[i].push_back({nearestRoad, walkTime});
            G.adj[nearestRoad].push_back({i, walkTime});
            G.edgeMode[i][nearestRoad] = "walk";
            G.edgeMode[nearestRoad][i] = "walk";
        }
    }
    
    // Add road ↔ bus transfers
    for(int i = 0; i < G.nodes.size(); i++){
        if(G.nodes[i].type != "bus") continue;
        if(G.nodes[i].name.empty()) continue;
        
        double busLon = G.nodes[i].longitude;
        double busLat = G.nodes[i].latitude;
        
        int nearestRoad = -1;
        double minDist = DBL_MAX;
        
        for(int j = 0; j < G.nodes.size(); j++){
            if(G.nodes[j].type != "road") continue;
            
            double roadLon = G.nodes[j].longitude;
            double roadLat = G.nodes[j].latitude;
            
            double d = haversine(busLon, busLat, roadLon, roadLat);
            
            if(d < minDist && d <= threshold){
                minDist = d;
                nearestRoad = j;
            }
        }
        
        if(nearestRoad != -1){
            double walkTime = (minDist / 2.0) * 60.0; // 2 km/h walking
            G.adj[i].push_back({nearestRoad, walkTime});
            G.adj[nearestRoad].push_back({i, walkTime});
            G.edgeMode[i][nearestRoad] = "walk";
            G.edgeMode[nearestRoad][i] = "walk";
        }
    }
    
    // Add metro ↔ bus transfers
    for(int i = 0; i < G.nodes.size(); i++){
        if(G.nodes[i].type != "metro") continue;
        if(G.nodes[i].name.empty()) continue;
        
        double metroLon = G.nodes[i].longitude;
        double metroLat = G.nodes[i].latitude;
        
        int nearestBus = -1;
        double minDist = DBL_MAX;
        
        for(int j = 0; j < G.nodes.size(); j++){
            if(G.nodes[j].type != "bus") continue;
            if(G.nodes[j].name.empty()) continue;
            
            double busLon = G.nodes[j].longitude;
            double busLat = G.nodes[j].latitude;
            
            double d = haversine(metroLon, metroLat, busLon, busLat);
            
            if(d < minDist && d <= threshold){
                minDist = d;
                nearestBus = j;
            }
        }
        
        if(nearestBus != -1){
            double walkTime = (minDist / 2.0) * 60.0; // 2 km/h walking
            G.adj[i].push_back({nearestBus, walkTime});
            G.adj[nearestBus].push_back({i, walkTime});
            G.edgeMode[i][nearestBus] = "walk";
            G.edgeMode[nearestBus][i] = "walk";
        }
    }
}

// Calculate waiting time for metro/bus (in minutes)
double calculateWaitTime(double arrivalTimeMinutes, string mode) {
    if (mode == "road" || mode == "walk") {
        return 0;
    }
    
    int dayMinutes = (int)arrivalTimeMinutes % (24 * 60);
    
    if (dayMinutes < 6 * 60 || dayMinutes >= 23 * 60) {
        return -1; // Outside operating hours
    }
    
    int minutesSince6AM = dayMinutes - 6 * 60;
    int minutesSinceLastDeparture = minutesSince6AM % 15;
    
    if (minutesSinceLastDeparture == 0) {
        return 0;
    }
    
    return 15 - minutesSinceLastDeparture;
}

string formatTime(double minutes) {
    int totalMins = (int)minutes;
    int hours = totalMins / 60;
    int mins = totalMins % 60;
    
    string period = (hours < 12) ? "AM" : "PM";
    int displayHour = hours % 12;
    if (displayHour == 0) displayHour = 12;
    
    ostringstream oss;
    oss << displayHour << ":";
    if (mins < 10) oss << "0";
    oss << mins << " " << period;
    return oss.str();
}

void initGlobalArray(int V){
    dist.resize(V);
    parent.resize(V);
    arrivalTime.resize(V);
}

void INITIALIZE_SINGLE_SOURCE(Graph& G, int src, double startTime){
    K.clear();
    while(!Q.empty()) Q.pop();
    
    for(int i = 0; i < G.V; i++){
        dist[i] = DBL_MAX;
        parent[i] = -1;
        arrivalTime[i] = DBL_MAX;
    }
    dist[src] = 0;
    arrivalTime[src] = startTime;
}

void Dijkstra(Graph& G, int src, double startTime){
    INITIALIZE_SINGLE_SOURCE(G, src, startTime);
    Q.push({0, src});
    
    while(!Q.empty()){
        int u = Q.top().second;
        double currentTime = Q.top().first;
        Q.pop();
        
        if(K.find(u) != K.end()) continue;
        if(currentTime > dist[u]) continue;
        K.insert(u);
        
        double currentArrivalTime = arrivalTime[u];
        
        for(auto& edge : G.adj[u]){
            int v = edge.first;
            double edgeTime = edge.second; // Time to travel this edge
            
            string mode = "road";
            if(G.edgeMode.find(u) != G.edgeMode.end() && 
               G.edgeMode[u].find(v) != G.edgeMode[u].end()){
                mode = G.edgeMode[u][v];
            }
            
            // Calculate waiting time
            double waitTime = 0;
            if(mode == "metro" || mode == "bus"){
                waitTime = calculateWaitTime(currentArrivalTime, mode);
                if(waitTime < 0){
                    continue; // Outside operating hours
                }
            }
            
            // Calculate new arrival time
            double newArrivalTime = currentArrivalTime + waitTime + edgeTime;
            
            // Check if metro/bus edge ends after 11 PM
            if(mode == "metro" || mode == "bus"){
                int dayMins = (int)newArrivalTime % (24 * 60);
                if(dayMins >= 23 * 60){
                    continue;
                }
            }
            
            // Calculate total time (including waiting)
            double totalTime = dist[u] + waitTime + edgeTime;
            
            // Relax edge - MINIMIZE TIME
            if(totalTime < dist[v]){
                dist[v] = totalTime;
                arrivalTime[v] = newArrivalTime;
                parent[v] = u;
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

int findNearestRoadNode(Graph& G, double lon, double lat) {
    int nearestNode = -1;
    double minDist = DBL_MAX;
    
    for(int i = 0; i < G.nodes.size(); i++) {
        if(G.nodes[i].type != "road") continue;
        
        double d = haversine(G.nodes[i].longitude, G.nodes[i].latitude, lon, lat);
        
        if(d < minDist) {
            minDist = d;
            nearestNode = i;
        }
    }
    
    return nearestNode;
}

void generateKML(Graph& G, vector<int>& path, string filename) {
    ofstream kml(filename);
    
    kml << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    kml << "<kml xmlns=\"http://www.opengis.net/kml/2.2\">\n";
    kml << "<Document>\n";
    kml << "  <name>Fastest Path - Dhaka (Problem 5)</name>\n";
    
    kml << "  <Style id=\"roadStyle\">\n";
    kml << "    <LineStyle>\n";
    kml << "      <color>ff0000ff</color>\n";
    kml << "      <width>4</width>\n";
    kml << "    </LineStyle>\n";
    kml << "  </Style>\n";
    
    kml << "  <Style id=\"metroStyle\">\n";
    kml << "    <LineStyle>\n";
    kml << "      <color>ffff0000</color>\n";
    kml << "      <width>4</width>\n";
    kml << "    </LineStyle>\n";
    kml << "  </Style>\n";
    
    kml << "  <Style id=\"busStyle\">\n";
    kml << "    <LineStyle>\n";
    kml << "      <color>ff00ff00</color>\n";
    kml << "      <width>4</width>\n";
    kml << "    </LineStyle>\n";
    kml << "  </Style>\n";
    
    kml << "  <Placemark>\n";
    kml << "    <name>Source</name>\n";
    kml << "    <Point>\n";
    kml << "      <coordinates>" << G.nodes[path[0]].longitude << "," 
        << G.nodes[path[0]].latitude << ",0</coordinates>\n";
    kml << "    </Point>\n";
    kml << "  </Placemark>\n";
    
    kml << "  <Placemark>\n";
    kml << "    <name>Destination</name>\n";
    kml << "    <Point>\n";
    kml << "      <coordinates>" << G.nodes[path.back()].longitude << "," 
        << G.nodes[path.back()].latitude << ",0</coordinates>\n";
    kml << "    </Point>\n";
    kml << "  </Placemark>\n";
    
    // Group path by mode
    for(int i = 0; i < path.size() - 1; i++){
        int u = path[i];
        int v = path[i + 1];
        
        string mode = "road";
        if(G.edgeMode.find(u) != G.edgeMode.end() && 
           G.edgeMode[u].find(v) != G.edgeMode[u].end()){
            mode = G.edgeMode[u][v];
        }
        
        string styleUrl = "#roadStyle";
        if(mode == "metro") styleUrl = "#metroStyle";
        else if(mode == "bus") styleUrl = "#busStyle";
        
        kml << "  <Placemark>\n";
        kml << "    <name>" << mode << "</name>\n";
        kml << "    <styleUrl>" << styleUrl << "</styleUrl>\n";
        kml << "    <LineString>\n";
        kml << "      <coordinates>\n";
        kml << "        " << G.nodes[u].longitude << "," << G.nodes[u].latitude << ",0\n";
        kml << "        " << G.nodes[v].longitude << "," << G.nodes[v].latitude << ",0\n";
        kml << "      </coordinates>\n";
        kml << "    </LineString>\n";
        kml << "  </Placemark>\n";
    }
    
    kml << "</Document>\n";
    kml << "</kml>\n";
    
    kml.close();
}

void printSegmentedPath(Graph& G, vector<int>& path, double srcLon, double srcLat, 
                        double destLon, double destLat, double startTimeMinutes,
                        double walkDistSrc, double walkDistDest) {
    cout << fixed << setprecision(2);
    
    cout << "Problem No: 5\n";
    cout << "Source: (" << srcLon << ", " << srcLat << ")\n";
    cout << "Destination: (" << destLon << ", " << destLat << ")\n";
    cout << "Starting time at source: " << formatTime(startTimeMinutes) << "\n\n";
    
    double currentTime = startTimeMinutes;
    
    // Walking from source to nearest road node
    if(walkDistSrc > 0) {
        double walkTime = (walkDistSrc / 2.0) * 60.0; // 2 km/h
        cout << formatTime(currentTime) << " - " << formatTime(currentTime + walkTime)
             << ": Walk from Source to (" << G.nodes[path[0]].longitude 
             << ", " << G.nodes[path[0]].latitude << ") [" 
             << (int)walkTime << " minutes]\n\n";
        currentTime += walkTime;
    }
    
    for(int i = 0; i < path.size() - 1; i++){
        int u = path[i];
        int v = path[i + 1];
        
        string mode = "road";
        if(G.edgeMode.find(u) != G.edgeMode.end() && 
           G.edgeMode[u].find(v) != G.edgeMode[u].end()){
            mode = G.edgeMode[u][v];
        }
        
        // Calculate waiting time
        double waitTime = 0;
        if(mode == "metro" || mode == "bus"){
            waitTime = calculateWaitTime(currentTime, mode);
        }
        
        if(waitTime > 0){
            cout << formatTime(currentTime) << " - " << formatTime(currentTime + waitTime)
                 << ": Wait for " << mode << " (" << (int)waitTime << " minutes)\n\n";
            currentTime += waitTime;
        }
        
        // Get edge time
        double edgeTime = 0;
        for(auto& edge : G.adj[u]){
            if(edge.first == v){
                edgeTime = edge.second;
                break;
            }
        }
        
        double endTime = currentTime + edgeTime;
        
        string uName = G.nodes[u].name.empty() ? 
                      "(" + to_string(G.nodes[u].longitude) + ", " + to_string(G.nodes[u].latitude) + ")" :
                      G.nodes[u].name + " (" + to_string(G.nodes[u].longitude) + ", " + to_string(G.nodes[u].latitude) + ")";
        
        string vName = G.nodes[v].name.empty() ? 
                      "(" + to_string(G.nodes[v].longitude) + ", " + to_string(G.nodes[v].latitude) + ")" :
                      G.nodes[v].name + " (" + to_string(G.nodes[v].longitude) + ", " + to_string(G.nodes[v].latitude) + ")";
        
        if(mode == "walk"){
            cout << formatTime(currentTime) << " - " << formatTime(endTime)
                 << ": Walk from " << uName << " to " << vName 
                 << " (" << (int)edgeTime << " minutes)\n\n";
        } else if(mode == "road"){
            cout << formatTime(currentTime) << " - " << formatTime(endTime)
                 << ": Ride Car from " << uName << " to " << vName 
                 << " (" << (int)edgeTime << " minutes)\n\n";
        } else if(mode == "metro"){
            cout << formatTime(currentTime) << " - " << formatTime(endTime)
                 << ": Ride Metro from " << uName << " to " << vName 
                 << " (" << (int)edgeTime << " minutes)\n\n";
        } else if(mode == "bus"){
            cout << formatTime(currentTime) << " - " << formatTime(endTime)
                 << ": Ride Bus from " << uName << " to " << vName 
                 << " (" << (int)edgeTime << " minutes)\n\n";
        }
        
        currentTime = endTime;
    }
    
    // Walking from last road node to destination
    if(walkDistDest > 0) {
        double walkTime = (walkDistDest / 2.0) * 60.0; // 2 km/h
        cout << formatTime(currentTime) << " - " << formatTime(currentTime + walkTime)
             << ": Walk from (" << G.nodes[path.back()].longitude 
             << ", " << G.nodes[path.back()].latitude << ") to Destination [" 
             << (int)walkTime << " minutes]\n\n";
        currentTime += walkTime;
    }
    
    cout << "Total travel time: " << (int)(currentTime - startTimeMinutes) << " minutes\n";
    cout << "Arrival time: " << formatTime(currentTime) << "\n";
}

int main(){
    Graph G;
    
    cout << "Reading road network...\n";
    readRoadCSV("../Dataset/Roadmap-Dhaka.csv", G);
    
    cout << "Reading metro network...\n";
    readMetroCSV("../Dataset/Routemap-DhakaMetroRail.csv", G);
    
    cout << "Reading bus network...\n";
    readBusCSV("../Dataset/Routemap-BikolpoBus.csv", G);
    
    G.V = G.nodes.size();
    
    cout << "Adding transfer edges...\n";
    addTransferEdges(G);
    
    cout << "Graph built successfully!\n";
    cout << "Nodes: " << G.V << endl;
    
    initGlobalArray(G.V);
    
    double srcLon, srcLat, destLon, destLat;
    int startHour, startMinute;
    
    cout << "Enter source (lon lat): ";
    cin >> srcLon >> srcLat;
    
    cout << "Enter destination (lon lat): ";
    cin >> destLon >> destLat;
    
    cout << "Enter start time (hour minute in 24h format): ";
    cin >> startHour >> startMinute;
    
    double startTimeMinutes = startHour * 60 + startMinute;
    
    // Find nearest road nodes
    cout << "\nFinding nearest road nodes..." << endl;
    int srcNode = findNearestRoadNode(G, srcLon, srcLat);
    int destNode = findNearestRoadNode(G, destLon, destLat);
    
    if(srcNode == -1){
        cout << "Error: No road node found near source!\n";
        return 1;
    }
    
    if(destNode == -1){
        cout << "Error: No road node found near destination!\n";
        return 1;
    }
    
    // Calculate walking distances
    double walkDistSrc = haversine(srcLon, srcLat, 
                                    G.nodes[srcNode].longitude, 
                                    G.nodes[srcNode].latitude);
    double walkDistDest = haversine(destLon, destLat, 
                                     G.nodes[destNode].longitude, 
                                     G.nodes[destNode].latitude);
    
    // Calculate walking time and adjust start time
    double walkTimeSrc = (walkDistSrc / 2.0) * 60.0; // 2 km/h in minutes
    double adjustedStartTime = startTimeMinutes + walkTimeSrc;
    
    cout << "Walking from source: " << walkDistSrc << " km (" 
         << (int)walkTimeSrc << " minutes)" << endl;
    cout << "Walking to destination: " << walkDistDest << " km (" 
         << (int)((walkDistDest / 2.0) * 60.0) << " minutes)" << endl;
    
    cout << "\nRunning Dijkstra...\n";
    Dijkstra(G, srcNode, adjustedStartTime);
    
    if(dist[destNode] == DBL_MAX){
        cout << "No path found!\n";
        return 1;
    }
    
    vector<int> path;
    printPath(G, destNode, path);
    
    cout << "\n========================================\n\n";
    printSegmentedPath(G, path, srcLon, srcLat, destLon, destLat, startTimeMinutes, 
                       walkDistSrc, walkDistDest);
    cout << "\n========================================\n\n";
    
    generateKML(G, path, "route_caseC.kml");
    cout << "KML file generated: route_caseC.kml\n";
    
    return 0;
}
