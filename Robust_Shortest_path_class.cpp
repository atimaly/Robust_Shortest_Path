#include <bits/stdc++.h>
#include <random>
#include <lemon/list_graph.h>
#include <lemon/time_measure.h>
#include <lemon/dijkstra.h>
#include <lemon/adaptors.h>
#include <lemon/concepts/path.h>
#include <lemon/concepts/graph.h>

using namespace lemon;
using namespace std;

using ll = long long int;

const bool DEBUG = false;
const long long int INF = 10000000000;

#define all(x) begin(x), end(x)
#define FOR(i,n) for(int i = 0; i < (n); ++i)

template <class C>
void Print_vector(const C &Original) {
	for(const auto &v : Original) {
	    cout << v << " ";
	}
	cout << endl;
}

template <class C>
void Print_Matrix(const vector<vector<C>> &M, bool space = true) {
	for(int i = 0; i < (int)M.size(); ++i) {
		for(int j = 0; j < (int)M[i].size(); ++j) {
			cout << M[i][j]; if(space) cout << " ";            
			}
	    cout << endl;
	}
}

template<class T, class C>
void Print_pair(const pair<T,C> &M) {
    cout << "(" << M.first << " , " << M.second << " ) ";
}

template <class C>
void Print_vector_pairs(const C &Original) {
	for(const auto &v : Original) {
	    Print_pair(v);
	}
	cout << endl;
}

template<class T, class C>
void Print_pair_os(const pair<T,C> &M, std::ostream &ot = std::cout) {
    ot << "(" << M.first << " , " << M.second << " ) ";
}

template <class C>
void Print_vector_pairs_os(const C &Original, std::ostream &ot = std::cout) {
	for(const auto &v : Original) {
	    Print_pair_os(v, ot); ot << endl;
	}
	ot << endl;
}


template<class T, class C>
void Print_Matrix_pair(const vector<vector<pair<T,C>>> &M) {
	for(int i = 0; i < (int)M.size(); ++i) {
		cout << i << ": ";
		for(int j = 0; j < (int)M[i].size(); ++j) {
			Print_pair(M[i][j]);          
		}
	    cout << endl;
	}
}

struct less_than_key
{
    inline bool operator() (const pair<double, ListDigraph::Arc>& struct1, const pair<double, ListDigraph::Arc>& struct2)
    {
        return (struct1.first > struct2.first);
    }
};


template <typename T>
pair<T, T> GeneratePoint(std::set<pair<T,T>> &Points) {
	std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0., 1.0);
    pair<T,T> temp = {dis(gen), dis(gen)};
    while(Points.find(temp) != Points.end()) {
		return temp;
	}
    
}

bool Deviation_Chance(double p) {
	//Is the diven edge will deviate
	std::random_device rd;
    std::mt19937 gen(rd());
    std::discrete_distribution<> d({1-p, p});
    int b = d(gen);
    if(b == 1) return true;
    return false;
}

double Deviation_Size(int m) {
	//How many times length[e] is given to length[e]
	std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, (double)m);
    return dis(gen);
}


double calculateDistance(std::pair<double, double> &x, std::pair<double, double> &y)
{
    return sqrt(pow(x.first - y.first, 2) +
                pow(x.second - y.second, 2));
}

class RobustShortestPathProblem{
	
	public:
		ListDigraph graph_;
		ListDigraph::ArcMap<double> length{graph_}; 
		ListDigraph::ArcMap<double> disturbances{graph_};
		int vertex_numb_; 
		double deviation_ch_;
		double erdos_ch_;
		vector<pair<double,double>> points_vec_;
		int edge_counter;
		vector<vector<int>> robust_shortest_paths_;
		vector<double> robust_shortest_paths_length_;
		
		
		RobustShortestPathProblem(int N, double devi, double erdos, std::ostream& os = std::cout) : vertex_numb_{N}, deviation_ch_{devi}, erdos_ch_{erdos} {
			CreateGraph(os);
			int max_gamma = 10;
			int gamma = 4;
			robust_shortest_paths_.resize(max_gamma+1);
			robust_shortest_paths_length_.resize(max_gamma+1, INF);
		}
		
		
		void CreateGraph(std::ostream& os = std::cout) {
			ListDigraph &g = graph_;
			int &N = vertex_numb_;
			
			cout << "N: " << N << endl;
			cout << "ERDOS IN: " << erdos_ch_ << endl;
			
			FOR(i, N) {
				g.addNode();
			}
			FOR(i,N) {
				for(int j = 1; j < N; ++j) {
					if(i != j && Deviation_Chance(erdos_ch_)) {
						//if(i != 0 && j != N-1) 
						ListDigraph::Arc e = g.addArc(g.nodeFromId(i), g.nodeFromId(j));
					}
				}
			}
			
			set<pair<double,double>> Points; Points.insert({0,0}); Points.insert({1,1});
			FOR(i, N-2) {
				Points.insert(GeneratePoint(Points));
			}
			
			if(DEBUG) {cout << "POINTS: \n"; Print_vector_pairs(Points); }
			
			for(auto &v : Points) points_vec_.push_back(v);
			//if(DEBUG) Print_vector_pairs(Points_vec);
			
				
			for(ListDigraph::ArcIt arc(g); arc!=INVALID; ++arc) {
				int u = g.id(g.source(arc)); int v = g.id(g.target(arc));
				length[arc] = calculateDistance(points_vec_[u], points_vec_[v]);
				disturbances[arc] = Deviation_Size(8)*length[arc];
				//if(DEBUG) cout << length[arc] << endl;
			}
		}
		
		void BasicData() {
			ListDigraph &g = graph_;
			int &N = vertex_numb_;
			
			edge_counter = 0;
			for(ListDigraph::ArcIt arc(g); arc!=INVALID; ++arc) 
				++edge_counter;
			cout << "THE GRAPH HAS: " << edge_counter << " edges\n";
			cout << "ERDOS: " << erdos_ch_ << endl;
		}
		
		void EdgesFout(std:: ostream &ot = std::cout) {
			ListDigraph &g = graph_;
			
			for(ListDigraph::ArcIt arc(g); arc!=INVALID; ++arc) {
				int u = g.id(g.source(arc)); int v = g.id(g.target(arc));
				//ot << u << " " << v << endl;
				ot << u << " " << v << " " << length[arc] << " " << disturbances[arc] << endl;
			}
			ot << endl;
		}
		
		void FoutShortestRobustPaths(vector<int> gammas, std:: ostream &ot = std::cout) {
			for(const auto &v : gammas) {
				for(const auto &u : robust_shortest_paths_[v]) {
					ot << u << " ";
				}
				ot << endl;
			}
		}
		
		double ShortestRobustPathSubProb(int l, vector<pair<double, ListDigraph::Arc>> &lengths_p, int gamma) {
			ListDigraph &g = graph_;
			int &N = vertex_numb_;
			
			double import = lengths_p[l].first;
			FOR(i,l+1) {
				length[lengths_p[i].second] += disturbances[lengths_p[i].second];
				length[lengths_p[i].second] -= import; 
			}
			
			if(DEBUG) 
				for(ListDigraph::ArcIt arc(g); arc!=INVALID; ++arc) {
					cout << length[arc] << endl;
				}
			
			ListDigraph::NodeMap<double> dist(g);
			Dijkstra<ListDigraph, ListDigraph::ArcMap<double>>
					  ::Create
					  dijkstra(g, length);
			//Dijkstra<ListGraph> dijkstra(g, length);
			dijkstra.distMap(dist);
			dijkstra.init();
			dijkstra.addSource(g.nodeFromId(0));
			dijkstra.start();
			
			FOR(i,l+1) {
				length[lengths_p[i].second] -= disturbances[lengths_p[i].second];
				length[lengths_p[i].second] += import; 
			}
			
			//Save path if better
				double cost =  gamma*lengths_p[l].first+dijkstra.dist(g.nodeFromId(N-1));
				if(cost < robust_shortest_paths_length_[gamma]) {
					vector<int> indexes;
					lemon::ListDigraph::Node curr = g.nodeFromId(N-1);
					while(curr != INVALID) {
						indexes.push_back(g.id(curr));
						curr = dijkstra.predNode(curr);
					}
					robust_shortest_paths_[gamma] = indexes;
				}
			
			
			return dijkstra.dist(g.nodeFromId(N-1));
		}

		double ShortestRobustPath(int gamma, int &min_ind) {
			ListDigraph &g = graph_;
			int &N = vertex_numb_;
			//d_{i,j} -s in order
				vector<pair<double, ListDigraph::Arc>> lengths_p;
				
				for(ListDigraph::ArcIt arc(g); arc!=INVALID; ++arc) {
					lengths_p.push_back({disturbances[arc], arc});
				}
				
				std::sort(all(lengths_p), less_than_key());
				if(DEBUG) { cout << "LENGTH_P: \n"; FOR(i,(int)lengths_p.size()) cout << lengths_p[i].first << " "; cout << endl; }
			
			
			double minim = 1000000;
			//l = n+1 case
				FOR(i,(int)lengths_p.size()) {
					length[lengths_p[i].second] += disturbances[lengths_p[i].second];
				}
				ListDigraph::NodeMap<double> dist(g);
				Dijkstra<ListDigraph, ListDigraph::ArcMap<double>>
						  ::Create
						  dijkstra(g, length);
				//Dijkstra<ListGraph> dijkstra(g, length);
				dijkstra.distMap(dist);
				dijkstra.init();
				dijkstra.addSource(g.nodeFromId(0));
				dijkstra.start();
				minim =  dijkstra.dist(g.nodeFromId(N-1));
				FOR(i,(int)lengths_p.size()) {
					length[lengths_p[i].second] -= disturbances[lengths_p[i].second];
				}
			
			min_ind = (int)lengths_p.size();
			if(DEBUG) cout << "MINIM: " << minim << endl;
			//Every other case
				FOR(l,(int)lengths_p.size()) {
					double curr_cost = gamma*lengths_p[l].first+ShortestRobustPathSubProb(l, lengths_p, gamma); //ShortestRobustPathSubProb(int l, &lengths_p)
					if(minim > curr_cost) {
						minim = curr_cost;
						min_ind = l;
					}
					if(DEBUG) cout << "MINIM: " << minim << endl;
				}
			
			cout << "MIN INDEX: " << min_ind << endl;
			
			return minim;
			
		}
		
		double ShortestPathRandomInstance(vector<int> &indexes, bool only_for_plot, std::ostream& os = std::cout) {
			ListDigraph &g = graph_;
			int &N = vertex_numb_;
			
			vector<ListDigraph::Arc> deviated;
			for(ListDigraph::ArcIt arc(g); arc!=INVALID; ++arc) {
				if(Deviation_Chance(deviation_ch_)) {
					os << g.id(g.source(arc)) << " " << g.id(g.target(arc)) << ", ";
					length[arc] += disturbances[arc];
					deviated.push_back(arc);
				}
			}
			
			ListDigraph::NodeMap<double> dist(g);
			Dijkstra<ListDigraph, ListDigraph::ArcMap<double>>
					  ::Create
					  dijkstra(g, length);
			//Dijkstra<ListGraph> dijkstra(g, length);
			dijkstra.distMap(dist);
			dijkstra.init();
			dijkstra.addSource(g.nodeFromId(0));
			dijkstra.start();
			
			for(ListDigraph::Arc arc : deviated) {
				length[arc] -= disturbances[arc];
			}
			
			if(only_for_plot) {
				lemon::ListDigraph::Node curr = g.nodeFromId(N-1);
				while(curr != INVALID) {
					indexes.push_back(g.id(curr));
					curr = dijkstra.predNode(curr);
				}
			}
			
			return dijkstra.dist(g.nodeFromId(N-1));
			
		}
		
		
};


int main() {
	ofstream foutp("Points.txt");
	ofstream fouts("Shortest_Paths.txt");
	ofstream foutl("Shortest_Lengths.txt");
	ofstream foutrp("Robust_Paths.txt");
	ofstream foutdev("Deviated_Edges.txt");
	
	int test = 20000; cout << "How many COST tests?\n"; cin >> test;
	int N; cout << "How many VERTICES?\n"; cin >> N; //How many parts should one divide [0,1] intervallum
	double deviation_ch; cout << "What is the chance of the cost deviating?\n"; cin >> deviation_ch;
	double erdos_ch; cout << "What is the chance of an (i,j) directed edge existing?\n"; cin >> erdos_ch;
	bool only_for_plot; cout << "Do you only want to plot?\n"; cin >> only_for_plot;
	
	Timer t;
	RobustShortestPathProblem Test(N, deviation_ch, erdos_ch, foutp);
	Test.BasicData();
	cout << "CREATING GRAPH: " << t << endl; t.restart();
	
	//Points.txt
		Test.BasicData(); cout << "EDGE COUNTER: " << Test.edge_counter << endl;
		foutp << N << endl; foutp << Test.edge_counter << endl;
		Print_vector_pairs_os(Test.points_vec_, foutp);
		Test.EdgesFout(foutp);
		
	
	//Robust Solutions
		int min_ind;
		vector<double> robust_shortest;
		for(int gamma : {0, 3, 6, 10}) {
			robust_shortest.push_back(Test.ShortestRobustPath(gamma, min_ind));
		}
		cout << "ROBUST SOLUTION: " << t << endl; t.restart();
	
		for(const auto v : robust_shortest) foutl << v << " ";
		foutl << endl;
		
		Test.FoutShortestRobustPaths({0, 3, 6, 10}, foutrp);
		
	//Random Instances
		double all = 0;
		int counter = 0;
		FOR(i,test) {
			vector<int> indexes;
			double temp = Test.ShortestPathRandomInstance(indexes, only_for_plot, foutdev);
			all += temp;
			foutl << temp << " ";
			
			if(only_for_plot)
				for(auto &v : indexes) {
					fouts << v << endl;
				}
			if(DEBUG) cout << "TEMP: " << temp << endl;
			
			if(counter >= 1000) {
				cout << "TEMP: " << temp << endl;			
				counter -= 1000;
			}
			++counter;
		}
		cout << "RANDOM INSTANCES: " << t << endl; t.restart();
	//Results
		cout << "ALL: " << all << endl;
		cout << "AVARAGE: " << all/test << endl;
		Print_vector(robust_shortest);
		
		foutl << endl << endl;
	
	
}
