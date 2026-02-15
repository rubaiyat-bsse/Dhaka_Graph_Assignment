#include<bits/stdc++.h>
using namespace std;

struct Node {
    int id;
    double longitude;
    double latitude;
    string type;
    string name;
};

struct Graph {
    int V;
    int E;
    vector<Node> nodes;
    map<pair<double, double>, int> coordToId;
    unordered_map<int, vector<pair<int, double>>> adj; // weight = TIME in minutes
    unordered_map<int, unordered_map<int, string>> edgeMode;
};

struct EdgeMatch {
    bool found;
    int node1;
    int node2;
    double distToNode1;
    double distToNode2;
};

vector<double> dist;
vector<int> parent;
vector<double> arrivalTime;
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
            double timeMinutes = (distance / 10.0) * 60.0;
            
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
            double timeMinutes = (distance / 10.0) * 60.0;
            
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
            double timeMinutes = (distance / 10.0) * 60.0;
            
            G.adj[node1].push_back({node2, timeMinutes});
            G.adj[node2].push_back({node1, timeMinutes});
            
            G.edgeMode[node1][node2] = "bus";
            G.edgeMode[node2][node1] = "bus";
        }
    }
    
    file.close();
}

void addTransferEdges(Graph& G, double threshold = 0.5) {
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
            double walkTime = (minDist / 2.0) * 60.0;
            G.adj[i].push_back({nearestRoad, walkTime});
            G.adj[nearestRoad].push_back({i, walkTime});
            G.edgeMode[i][nearestRoad] = "walk";
            G.edgeMode[nearestRoad][i] = "walk";
        }
    }
    
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
            double walkTime = (minDist / 2.0) * 60.0;
            G.adj[i].push_back({nearestRoad, walkTime});
            G.adj[nearestRoad].push_back({i, walkTime});
            G.edgeMode[i][nearestRoad] = "walk";
            G.edgeMode[nearestRoad][i] = "walk";
        }
    }
    
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
            double walkTime = (minDist / 2.0) * 60.0;
            G.adj[i].push_back({nearestBus, walkTime});
            G.adj[nearestBus].push_back({i, walkTime});
            G.edgeMode[i][nearestBus] = "walk";
            G.edgeMode[nearestBus][i] = "walk";
        }
    }
}

double calculateWaitTime(double arrivalTimeMinutes, string mode) {
    if (mode == "road" || mode == "walk") {
        return 0;
    }
    
    int dayMinutes = (int)arrivalTimeMinutes % (24 * 60);
    
    if (dayMinutes < 6 * 60 || dayMinutes >= 23 * 60) {
        return -1;
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

EdgeMatch findEdgeForPoint(Graph& G, double lon, double lat) {
    double tolerance = 0.00001;
    
    for(int u = 0; u < G.nodes.size(); u++) {
        if(G.nodes[u].type != "road") continue;
        
        for(auto& [v, edgeTime] : G.adj[u]) {
            if(G.nodes[v].type != "road") continue;
            if(u < v) {
                double distUV = haversine(G.nodes[u].longitude, G.nodes[u].latitude,
                                         G.nodes[v].longitude, G.nodes[v].latitude);
                
                double distUP = haversine(G.nodes[u].longitude, G.nodes[u].latitude, 
                                         lon, lat);
                double distPV = haversine(lon, lat,
                                         G.nodes[v].longitude, G.nodes[v].latitude);
                
                double diff = abs((distUP + distPV) - distUV);
                
                if(diff < tolerance) {
                    return {true, u, v, distUP, distPV};
                }
            }
        }
    }
    
    return {false, -1, -1, 0, 0};
}

void DijkstraFromEdge(Graph& G, EdgeMatch src, double startTimeMinutes) {
    K.clear();
    while(!Q.empty()) Q.pop();
    
    int V = G.nodes.size();
    dist.assign(V, DBL_MAX);
    parent.assign(V, -1);
    arrivalTime.assign(V, DBL_MAX);
    
    double travelTimeToNode1 = (src.distToNode1 / 10.0) * 60.0;
    double travelTimeToNode2 = (src.distToNode2 / 10.0) * 60.0;
    
    dist[src.node1] = travelTimeToNode1;
    dist[src.node2] = travelTimeToNode2;
    arrivalTime[src.node1] = startTimeMinutes + travelTimeToNode1;
    arrivalTime[src.node2] = startTimeMinutes + travelTimeToNode2;
    
    Q.push({dist[src.node1], src.node1});
    Q.push({dist[src.node2], src.node2});
    
    while(!Q.empty()) {
        int u = Q.top().second;
        double currentTime = Q.top().first;
        Q.pop();
        
        if(K.count(u)) continue;
        if(currentTime > dist[u]) continue;
        K.insert(u);
        
        double currentArrivalTime = arrivalTime[u];
        
        for(auto& edge : G.adj[u]) {
            int v = edge.first;
            double edgeTime = edge.second;
            
            string mode = "road";
            if(G.edgeMode.find(u) != G.edgeMode.end() && 
               G.edgeMode[u].find(v) != G.edgeMode[u].end()){
                mode = G.edgeMode[u][v];
            }
            
            double waitTime = 0;
            if(mode == "metro" || mode == "bus"){
                waitTime = calculateWaitTime(currentArrivalTime, mode);
                if(waitTime < 0) continue;
            }
            
            double newArrivalTime = currentArrivalTime + waitTime + edgeTime;
            
            if(mode == "metro" || mode == "bus"){
                int dayMins = (int)newArrivalTime % (24 * 60);
                if(dayMins >= 23 * 60) continue;
            }
            
            double totalTime = dist[u] + waitTime + edgeTime;
            
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

void generateKML(Graph& G, vector<int>& path, string filename, double srcLon, double srcLat, 
                 double destLon, double destLat) {
    ofstream kml(filename);
    
    kml << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    kml << "<kml xmlns=\"http://www.opengis.net/kml/2.2\">\n";
    kml << "<Document>\n";
    kml << "  <name>Fastest Path - Dhaka (Problem 5 Case B)</name>\n";
    
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
    kml << "      <coordinates>" << srcLon << "," << srcLat << ",0</coordinates>\n";
    kml << "    </Point>\n";
    kml << "  </Placemark>\n";
    
    kml << "  <Placemark>\n";
    kml << "    <name>Destination</name>\n";
    kml << "    <Point>\n";
    kml << "      <coordinates>" << destLon << "," << destLat << ",0</coordinates>\n";
    kml << "    </Point>\n";
    kml << "  </Placemark>\n";
    
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
                        EdgeMatch srcEdge, EdgeMatch destEdge) {
    cout << fixed << setprecision(2);
    
    cout << "Problem No: 5\n";
    cout << "Source: (" << srcLon << ", " << srcLat << ")\n";
    cout << "Destination: (" << destLon << ", " << destLat << ")\n";
    cout << "Starting time at source: " << formatTime(startTimeMinutes) << "\n\n";
    
    double currentTime = startTimeMinutes;
    
    // Initial segment from source to first node
    if(path[0] == srcEdge.node1) {
        double travelTime = (srcEdge.distToNode1 / 10.0) * 60.0;
        cout << formatTime(currentTime) << " - " << formatTime(currentTime + travelTime)
             << ": Ride Car from Source to (" << G.nodes[srcEdge.node1].longitude 
             << ", " << G.nodes[srcEdge.node1].latitude << ") [" 
             << (int)travelTime << " minutes]\n\n";
        currentTime += travelTime;
    } else if(path[0] == srcEdge.node2) {
        double travelTime = (srcEdge.distToNode2 / 10.0) * 60.0;
        cout << formatTime(currentTime) << " - " << formatTime(currentTime + travelTime)
             << ": Ride Car from Source to (" << G.nodes[srcEdge.node2].longitude 
             << ", " << G.nodes[srcEdge.node2].latitude << ") [" 
             << (int)travelTime << " minutes]\n\n";
        currentTime += travelTime;
    }
    
    for(int i = 0; i < path.size() - 1; i++){
        int u = path[i];
        int v = path[i + 1];
        
        string mode = "road";
        if(G.edgeMode.find(u) != G.edgeMode.end() && 
           G.edgeMode[u].find(v) != G.edgeMode[u].end()){
            mode = G.edgeMode[u][v];
        }
        
        double waitTime = 0;
        if(mode == "metro" || mode == "bus"){
            waitTime = calculateWaitTime(currentTime, mode);
        }
        
        if(waitTime > 0){
            cout << formatTime(currentTime) << " - " << formatTime(currentTime + waitTime)
                 << ": Wait for " << mode << " [" << (int)waitTime << " minutes]\n\n";
            currentTime += waitTime;
        }
        
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
                      G.nodes[u].name;
        
        string vName = G.nodes[v].name.empty() ? 
                      "(" + to_string(G.nodes[v].longitude) + ", " + to_string(G.nodes[v].latitude) + ")" :
                      G.nodes[v].name;
        
        if(mode == "walk"){
            cout << formatTime(currentTime) << " - " << formatTime(endTime)
                 << ": Walk from " << uName << " to " << vName 
                 << " [" << (int)edgeTime << " minutes]\n\n";
        } else if(mode == "road"){
            cout << formatTime(currentTime) << " - " << formatTime(endTime)
                 << ": Ride Car from " << uName << " to " << vName 
                 << " [" << (int)edgeTime << " minutes]\n\n";
        } else if(mode == "metro"){
            cout << formatTime(currentTime) << " - " << formatTime(endTime)
                 << ": Ride Metro from " << uName << " to " << vName 
                 << " [" << (int)edgeTime << " minutes]\n\n";
        } else if(mode == "bus"){
            cout << formatTime(currentTime) << " - " << formatTime(endTime)
                 << ": Ride Bus from " << uName << " to " << vName 
                 << " [" << (int)edgeTime << " minutes]\n\n";
        }
        
        currentTime = endTime;
    }
    
    // Final segment to destination
    int lastNode = path.back();
    if(lastNode == destEdge.node1) {
        double travelTime = (destEdge.distToNode1 / 10.0) * 60.0;
        cout << formatTime(currentTime) << " - " << formatTime(currentTime + travelTime)
             << ": Ride Car from (" << G.nodes[destEdge.node1].longitude 
             << ", " << G.nodes[destEdge.node1].latitude << ") to Destination [" 
             << (int)travelTime << " minutes]\n\n";
        currentTime += travelTime;
    } else if(lastNode == destEdge.node2) {
        double travelTime = (destEdge.distToNode2 / 10.0) * 60.0;
        cout << formatTime(currentTime) << " - " << formatTime(currentTime + travelTime)
             << ": Ride Car from (" << G.nodes[destEdge.node2].longitude 
             << ", " << G.nodes[destEdge.node2].latitude << ") to Destination [" 
             << (int)travelTime << " minutes]\n\n";
        currentTime += travelTime;
    }
    
    cout << "Total travel time: " << (int)(currentTime - startTimeMinutes) << " minutes\n";
    cout << "Arrival time: " << formatTime(currentTime) << "\n";
}

int main(){
    Graph G;
    
    cout << "Reading road network..." << endl;
    readRoadCSV("../Dataset/Roadmap-Dhaka.csv", G);
    
    cout << "Reading metro network..." << endl;
    readMetroCSV("../Dataset/Routemap-DhakaMetroRail.csv", G);
    
    cout << "Reading bus network..." << endl;
    readBusCSV("../Dataset/Routemap-BikolpoBus.csv", G);
    
    cout << "Adding transfer edges..." << endl;
    addTransferEdges(G, 0.5);
    
    G.V = G.nodes.size();
    
    cout << "Graph built successfully!\n";
    cout << "Nodes: " << G.V << endl;
    
    double srcLon, srcLat, destLon, destLat;
    int startHour, startMinute;
    
    cout << "Enter source (lon lat): ";
    cin >> srcLon >> srcLat;
    
    cout << "Enter destination (lon lat): ";
    cin >> destLon >> destLat;
    
    cout << "Enter start time (hour minute in 24h format): ";
    cin >> startHour >> startMinute;
    
    double startTimeMinutes = startHour * 60 + startMinute;
    
    cout << "\nFinding edges..." << endl;
    EdgeMatch srcEdge = findEdgeForPoint(G, srcLon, srcLat);
    EdgeMatch destEdge = findEdgeForPoint(G, destLon, destLat);
    
    if(!srcEdge.found){
        cout << "Error: Source point not on any road edge!\n";
        return 1;
    }
    
    if(!destEdge.found){
        cout << "Error: Destination point not on any road edge!\n";
        return 1;
    }
    
    cout << "Running Dijkstra from edge..." << endl;
    DijkstraFromEdge(G, srcEdge, startTimeMinutes);
    
    double finalTime1 = dist[destEdge.node1] + (destEdge.distToNode1 / 10.0) * 60.0;
    double finalTime2 = dist[destEdge.node2] + (destEdge.distToNode2 / 10.0) * 60.0;
    
    int destNode;
    double finalTime;
    if(finalTime1 < finalTime2){
        destNode = destEdge.node1;
        finalTime = finalTime1;
    } else {
        destNode = destEdge.node2;
        finalTime = finalTime2;
    }
    
    if(dist[destNode] == DBL_MAX){
        cout << "No path found!\n";
        return 1;
    }
    
    vector<int> path;
    printPath(G, destNode, path);
    
    cout << "\n========================================\n\n";
    printSegmentedPath(G, path, srcLon, srcLat, destLon, destLat, startTimeMinutes, srcEdge, destEdge);
    cout << "\n========================================\n\n";
    
    generateKML(G, path, "route_caseB.kml", srcLon, srcLat, destLon, destLat);
    cout << "KML file generated: route_caseB.kml\n";
    
    return 0;
}
