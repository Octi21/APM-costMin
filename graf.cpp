#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
using namespace std;


class Graf
{
	int n, m;
	vector <vector<pair<int, int>>> adiacenta;

	int drumMin(vector <bool> viz, vector <int> cost);

public:
	Graf(int n, int m, vector <vector<pair<int, int>>> adiacenta);

	void algPrim();
	void primPQ();
	void dijkstra();
	void BF();
};

Graf::Graf(int n, int m, vector <vector<pair<int, int>>> adiacenta)
{
	this->n = n;
	this->m = m;
	this->adiacenta = adiacenta;
}


int Graf::drumMin(vector <bool> viz, vector <int> cost) // select drum min
{
	int nod, costMin = 1001;
	for (int i = 1; i <= n; ++i)
	{
		if (viz[i] == false && cost[i] < costMin)
		{
			costMin = cost[i];
			nod = i;
		}
	}
	return nod;
}

void Graf::algPrim()
{
	ofstream out("graf.out");
	vector <bool> viz(n + 1, false);
	vector <int> cost(n + 1, 1001);
	vector <int> par(n + 1, -2);
	cost[1] = 0;
	par[1] = -1;
	int nodp;

	for (int i = 1; i <= n; ++i)
	{
		nodp = drumMin(viz, cost);
		viz[nodp] = true;
		//cout << nodp << " " << cost[nodp] << "\n";
		for (auto nodc : adiacenta[nodp])
		{
			if (viz[nodc.first] == false && cost[nodc.first] > nodc.second)
			{
				par[nodc.first] = nodp;
				cost[nodc.first] = nodc.second;
			}
		}
	}
	int sum = 0;
	for (int i = 1; i <= n; ++i)
	{
		//cout << cost[i] << " ";
		sum += cost[i];
	}
	out << sum << "\n" << n - 1 << "\n";
	for (int i = 2; i <= n; ++i)
		out << i << " " << par[i] << "\n";

}

void p1()
{
	ifstream in("graf.in");

	int n, m, x, y, z;
	in >> n >> m;
	vector <vector<pair<int, int>>> adiacenta(n + 1);

	for (int i = 1; i <= m; ++i)
	{
		in >> x >> y >> z;
		adiacenta[x].push_back(make_pair(y, z));
		adiacenta[y].push_back(make_pair(x, z));
	}
	Graf g(n, m, adiacenta);
	g.algPrim();
}


void Graf::primPQ()
{
	vector <int> par(n + 1, -1);
	vector <bool> viz(n + 1, false);
	vector <int> cost(n + 1, 1001);
	priority_queue <pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;
	int nodp, costf = 0;
	cost[1] = 0;
	pq.push(make_pair(0, 1));

	while (!pq.empty())
	{
		nodp = pq.top().second;
		pq.pop();

		if (viz[nodp] == false)
		{
			viz[nodp] = true;

			for (auto nodc : adiacenta[nodp])
			{
				if (nodc.second < cost[nodc.first] && viz[nodc.first] == false)
				{
					cost[nodc.first] = nodc.second;
					par[nodc.first] = nodp;

					pq.push(make_pair(nodc.second, nodc.first));
				}
			}
		}

	}
	ofstream out("graf.out");
	for (int i = 1; i <= n; ++i)
	{
		costf += cost[i];
	}
	out << costf << "\n";
	out << n - 1 << "\n";
	for (int i = 2; i <= n; ++i)
	{
		out << i << " " << par[i] << "\n";
	}

}

void p2()
{
	ifstream in("graf.in");

	int n, m, x, y, z;
	in >> n >> m;
	vector <vector<pair<int, int>>> adiacenta(n + 1);

	for (int i = 1; i <= m; ++i)
	{
		in >> x >> y >> z;
		adiacenta[x].push_back(make_pair(y, z));
		adiacenta[y].push_back(make_pair(x, z));
	}
	Graf g(n, m, adiacenta);
	
	g.primPQ();
}



void Graf::dijkstra()
{
	vector <int> cost(n + 2, 100001);
	vector <bool> viz(n + 2, false);
	priority_queue <pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;
	cost[1] = 0;
	int nodp;
	pq.push(make_pair(0, 1));
	while (!pq.empty())
	{
		nodp = pq.top().second;
		pq.pop();

		if (viz[nodp] == false)
		{
			viz[nodp] = true;
			for (auto nodc : adiacenta[nodp])
			{
				if (cost[nodp] + nodc.second < cost[nodc.first])
				{
					cost[nodc.first] = cost[nodp] + nodc.second;
					pq.push(make_pair(cost[nodc.first], nodc.first));
				}
			}
		}
		
	}
	ofstream out("dijkstra.out");
	for (int i = 2; i <= n; ++i)
	{
		if (cost[i] != 100001)
			out << cost[i] << " ";
		else
			out << 0 << " ";

	}
		
}

void p3()
{
	ifstream in("dijkstra.in");
	
	int n, m, x, y, z;
	in >> n >> m;
	vector <vector<pair<int, int>>> adiacenta(n + 2);

	for (int i = 1; i <= m; ++i)
	{
		in >> x >> y >> z;
		adiacenta[x].push_back(make_pair(y, z));
	}
	Graf g(n, m, adiacenta);

	g.dijkstra();
	
}


void Graf::BF()
{
	queue <int> que;
	vector <bool> rnd(n + 1, false);
	vector <int> cost(n + 1, 100001);
	vector <int> relax(n + 1, 0);
	int nodp;
	bool cicluN = false;
	cost[1] = 0;
	que.push(1);
	//rnd[1] = true;
	//relax[1] = 1;

	while (!que.empty())
	{
		nodp = que.front();
		que.pop();
		rnd[nodp] = false;

		for (auto nodc : adiacenta[nodp])
		{
			if (cost[nodp] + nodc.second < cost[nodc.first])
			{
				relax[nodc.first] ++;
				cost[nodc.first] = cost[nodp] + nodc.second;

				if (rnd[nodc.first] == false)
				{
					rnd[nodc.first] = true;
					que.push(nodc.first);
				}
				if (relax[nodc.first] == n)
				{
					cicluN = true;
					break;
				}
			}
		}
		if (cicluN == true)
			break;
	}

	ofstream out("bellmanford.out");
	if (cicluN == true)
		out << "Ciclu negativ!";
	else
	{
		for (int i = 2; i <= n; ++i)
		{
			if (cost[i] != 100001)
				out << cost[i] << " ";
			else
				out << 1 << " ";
		}
	}
}

void p4()
{
	ifstream in("bellmanford.in");

	int n, m, x, y, z;
	in >> n >> m;
	vector <vector<pair<int, int>>> adiacenta(n + 1);

	for (int i = 1; i <= m; ++i)
	{
		in >> x >> y >> z;
		adiacenta[x].push_back(make_pair(y, z));
	}
	Graf g(n, m, adiacenta);

	g.BF();
}


class Disjoint
{
private:
	int n, m;
	int rang[100001], par[100001];
	
	int cautR(int x);
	void unite(int x, int y);
public:

	void afis();
};

int Disjoint::cautR(int x)
{
	if (par[x] != x)
		cautR(par[x]);
	else
		return x;
}

void Disjoint::unite(int x, int y)
{
	int px = cautR(x);
	int py = cautR(y);
	if (rang[px] == rang[py])
	{
		par[px] = py;
		rang[py] ++;
	}
	else if (rang[px] > rang[py])
	{
		par[py] = px;
	}
	else
		par[px] = py;
}

void Disjoint:: afis()
{
	int  com,x,y;
	ifstream in("graf.in");
	ofstream out("graf.out");
	in >> n >> m;
	for (int i = 1; i <= n; ++i)
	{
		rang[i] = 1;
		par[i] = i;
	}
	for (int i = 1; i <= m; ++i)
	{
		in >> com >> x >> y;
		if (com == 1)
		{
			cout << "pentru1 : \n";
			unite(x, y);
			cout << "rang = ";
			for (int i = 1; i <= n; ++i)
			{
				cout << rang[i]<<" ";
			}
			cout << "\n";
			cout << "par = ";
			for (int i = 1; i <= n; ++i)
			{
				cout << par[i] << " ";
			}
			cout << "\n";

		}
		else
		{
			cout << "pt2 : \n";
			cout << "rad x = " << cautR(x) << "//// rad y = " << cautR(y);
			cout << "\n";
			if (cautR(x) == cautR(y))
				out << "DA\n";
			else
				out << "NU\n";
		}
	}

}

void p5()
{
	Disjoint d;
	d.afis();
}



int main()
{
	//p1();
	//p2();
	//p3();
	//p4();
	p5();
}










