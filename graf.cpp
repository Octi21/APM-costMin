#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <stack>
#define LIMIT 100020
using namespace std;


class Graf
{
	int n, m;
	vector <vector<pair<int, int>>> adiacenta;

	vector <vector <int>> adiacenta2;
	vector <vector<int>> rev;


	int drumMin(vector <bool> viz, vector <int> cost); // pt alg prim



	void dfs(int nod, bool viz[]); // prob dfs

	int minim(int a, int b);
	void dfsCc(int nodp, vector<int>& prnt, vector<int>& desc, vector<int>& ret, stack<pair<int, int>>& st, int& cc, vector <vector <int>>& compCon); // comp bicon

	void dfsNorm(int& nod, vector <bool>& viz, stack <int>& st);
	void  dfsNorm2(int& nod, vector <bool>& viz, int& nrC, vector <vector<int>>& solutii); // tare conexe

	void dfsMC(int nodp, vector<int>& prnt, vector<int>& desc, vector<int>& ret, stack <pair<int, int>>& st); // muchie critica

public:
	Graf(int n, int m, vector <vector<pair<int, int>>> adiacenta);
	Graf(int n, int m, vector <vector <int>> adiacenta2);
	Graf(int n, int m, vector <vector <int>> adiacenta2, vector <vector <int>> rev);

	void algPrim();
	void primPQ();
	void dijkstra();
	void BF();

	vector<int> bfs(int nod); // prob1 bfs
	int compCon(); // prob2 comn con cu dfs
	void PA(); // prob3 compBicon
	void tareConexC(); // prob4
	void muchieC(); //prob 5 muchie critica
	void topologic(); // prob6 sortTopologica
};

Graf::Graf(int n, int m, vector <vector<pair<int, int>>> adiacenta)
{
	this->n = n;
	this->m = m;
	this->adiacenta = adiacenta;
}

Graf::Graf(int n, int m, vector <vector <int>> adiacenta2)
{
	this->n = n;
	this->m = m;
	this->adiacenta2 = adiacenta2;
}

Graf::Graf(int n, int m, vector <vector <int>> adiacenta2, vector <vector<int>> rev)
{
	this->n = n;
	this->m = m;
	this->adiacenta2 = adiacenta2;
	this->rev = rev;

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


int Graf::minim(int a, int b)
{
	if (a > b)
		return b;
	return a;
}

void Graf::dfsCc(int nodp, vector<int>& prnt, vector<int>& desc, vector<int>& ret, stack<pair<int, int>>& st, int& cc, vector <vector <int>>& compCon)
{
	static int timp = 1;
	desc[nodp] = timp;
	ret[nodp] = timp;
	timp++;
	int nrc = 0;

	for (auto nodc : adiacenta2[nodp])
	{
		if (desc[nodc] == -1)
		{
			nrc++;
			prnt[nodc] = nodp;
			st.push(make_pair(nodp, nodc));
			dfsCc(nodc, prnt, desc, ret, st, cc, compCon);

			ret[nodp] = minim(ret[nodp], ret[nodc]);

			if ((prnt[nodp] == -1 && nrc > 1) || (prnt[nodp] != -1 && desc[nodp] <= ret[nodc]))  // verific noduri sunt ap
			{
				//ap[nodp] = true; => nod ul pe care esti nod de articulatie

				while (st.top().first != nodp && st.top().second != nodc)
				{
					compCon[cc].push_back(st.top().second);
					st.pop();
				}
				compCon[cc].push_back(st.top().second);
				compCon[cc].push_back(st.top().first);

				st.pop();

				cc++;

			}
		}
		else if (nodc != prnt[nodp])
		{
			ret[nodp] = minim(ret[nodp], desc[nodc]);
		}
	}

}

void Graf::PA()
{
	vector<int> prnt(n + 1, -1), desc(n + 1, -1), ret(n + 1, -1);
	vector<bool> ap(n + 1, false);
	stack<pair<int, int>> st;
	vector <vector <int>> compCon(n + 1);
	int cc = 1;
	int chestie;

	for (int i = 1; i <= n; ++i)
	{
		if (desc[i] == -1)
		{
			dfsCc(i, prnt, desc, ret, st, cc, compCon);
		}
		while (!st.empty())
		{
			compCon[cc].push_back(st.top().second);
			chestie = st.top().first;
			st.pop();
		}
	}
	compCon[cc].push_back(chestie);
	ofstream out("graf.out");
	out << cc << "\n";
	for (int i = 1; i <= cc; ++i)
	{
		for (auto j : compCon[i])
		{
			out << j << " ";
		}
		out << "\n";
	}
}

void p6()
{
	ifstream in("graf.in");
	int n, m, x, y;
	in >> n >> m;
	vector <vector <int>> adiacenta2 (n + 1);

	for (int i = 0; i < m; ++i)
	{
		in >> x >> y;
		adiacenta2[x].push_back(y);
		adiacenta2[y].push_back(x);
	}

	Graf g1 = Graf(n, m, adiacenta2);

	g1.PA();
}



vector <int> Graf::bfs(int nod)
{
	vector <int> dist;

	queue <int> coada;

	for (int i = 0; i <= n; ++i)
	{
		dist.push_back(-1);
	}
	coada.push(nod);

	dist[nod] = 0;

	while (!coada.empty())
	{
		int p = coada.front();
		for (int i = 0; i < adiacenta2[p].size(); ++i)
		{
			if (dist[adiacenta2[p][i]] == -1)
			{
				dist[adiacenta2[p][i]] = dist[p] + 1;
				coada.push(adiacenta2[p][i]);
			}
		}
		coada.pop();
	}
	return dist;
}

void p7()
{
	ifstream in("graf.in");
	
	int n, m, nod, x, y;
	in >> n >> m >> nod;
	vector <vector <int>> adiacenta2 (n + 1);



	for (int i = 0; i < m; ++i)
	{
		in >> x >> y;
		adiacenta2 [x].push_back(y);
	}

	Graf g(n, m, adiacenta2);
	vector<int> r = g.bfs(nod);


	ofstream out("graf.out");
	for (int i = 1; i <= n; ++i)
		out << r[i] << " ";
}



void Graf::dfs(int nod, bool viz[])
{
	viz[nod] = true;
	for (int i = 0; i < adiacenta2[nod].size(); ++i)
	{
		if (viz[adiacenta2[nod][i]] == false)
		{
			viz[adiacenta2[nod][i]] = true;
			dfs(adiacenta2[nod][i], viz);
		}
	}
}

int Graf::compCon()
{
	int nrComp = 0;
	bool* viz = new bool[n + 1];
	for (int i = 0; i <= n; ++i)
	{
		viz[i] = false;
	}
	for (int i = 1; i <= n; ++i)
	{
		if (viz[i] == false)
		{
			dfs(i, viz);
			nrComp++;
		}
	}
	return nrComp;

}

void p8()
{

	ifstream in("graf.in");
	int n, m, x, y;
	in >> n >> m;
	vector <vector <int>> adiacenta(n + 1);

	for (int i = 0; i < m; ++i)
	{
		in >> x >> y;
		adiacenta[x].push_back(y);
		adiacenta[y].push_back(x);
	}
	Graf g1 = Graf(n, m, adiacenta);
	int r = g1.compCon();

	ofstream out("graf.out");
	out << r;

}



void Graf::dfsNorm(int& nod, vector <bool>& viz, stack <int>& st)
{
	viz[nod] = true;
	//cout<<nod<<" ";
	for (int i : adiacenta2[nod])
	{
		if (viz[i] == false)
		{
			//viz[nod] = true;
			dfsNorm(i, viz, st);
			//st.push(i);

		}
	}
	st.push(nod);

}

void Graf::dfsNorm2(int& nod, vector <bool>& viz, int& nrC, vector <vector<int>>& solutii)
{
	viz[nod] = true;
	solutii[nrC].push_back(nod);
	//out<<nod<<" ";
	for (int i : rev[nod])
	{
		if (viz[i] == false)
		{
			viz[nod] = true;
			dfsNorm2(i, viz/*,rev*/, nrC, solutii);
			//st.push(i);

		}
	}

}

void Graf::tareConexC()
{
	stack <int> st;
	vector <bool> viz(n + 1, false);
	vector <vector<int>>  solutii(n + 1);
	int nrC = 0;

	for (int i = 1; i <= n; ++i)
	{
		if (viz[i] == false)
		{
			dfsNorm(i, viz, st);
			//st.push(i);
		}
	}


	for (int i = 1; i <= n; ++i)
		viz[i] = false;

	while (!st.empty())
	{
		if (viz[st.top()] == false)
		{
			nrC++;
			dfsNorm2(st.top(), viz, nrC, solutii);
		}
		st.pop();
	}
	ofstream out("graf.out");
	out << nrC << "\n";
	for (int i = 1; i <= nrC; ++i)
	{
		for (auto j : solutii[i])
		{
			out << j << " ";
		}
		out << "\n";
	}

}

void p9()
{
	ifstream in("graf.in");
	
	int n, m, x, y;
	in >> n >> m;
	vector <vector<int>> adiacenta(n + 1), rev(n + 1);
	for (int i = 1; i <= m; ++i)
	{
		in >> x >> y;
		adiacenta[x].push_back(y);
		rev[y].push_back(x);
	}
	Graf g(n, m, adiacenta, rev);
	g.tareConexC();
}



void Graf::dfsMC(int nodp, vector<int>& prnt, vector<int>& desc, vector<int>& ret, stack <pair<int, int>>& st)
{
	static int timp = 1;
	desc[nodp] = timp;
	ret[nodp] = timp;
	timp++;

	for (auto nodc : adiacenta2[nodp])
	{
		if (desc[nodc] == -1)
		{
			prnt[nodc] = nodp;
			dfsMC(nodc, prnt, desc, ret, st);

			ret[nodp] = minim(ret[nodp], ret[nodc]);


			if (desc[nodp] < ret[nodc]) // verific muchie crit
			{
				st.push(make_pair(nodp, nodc));
			}
		}
		else if (nodc != prnt[nodp])
		{
			ret[nodp] = minim(ret[nodp], desc[nodc]);
		}
	}

}

void Graf::muchieC()
{
	vector<int> prnt(n + 1, -1), desc(n + 1, -1), ret(n + 1, -1);
	stack <pair<int, int>> st;

	for (int i = 1; i <= n; ++i)
	{
		if (desc[i] == -1)
		{
			dfsMC(i, prnt, desc, ret, st);
		}
	}
	ofstream out("graf.out");
	out << "[";
	while (!st.empty())
	{
		out << " [" << st.top().first << "," << st.top().second << "] ";
		st.pop();
		if (!st.empty())
			out << " , ";
	}
	out << "]";
}

void p10()
{
	ifstream in("graf.in");
	int n, m, x, y;
	in >> n >> m;
	vector <vector<int>> adiacenta(n + 1);
	for (int i = 1; i <= m; ++i)
	{
		in >> x >> y;
		adiacenta[x].push_back(y);
		adiacenta[y].push_back(x);
	}
	Graf g(n, m, adiacenta);
	g.muchieC();
}



void Graf::topologic()
{
	vector <bool> viz(n + 1, false);
	vector <vector<int>> adiacenta;
	stack <int> st;

	for (int i = 1; i <= n; ++i)
	{
		if (!viz[i])
			dfsNorm(i, viz, st);
	}
	ofstream out("graf.out");
	while (!st.empty())
	{
		out << st.top() << " ";
		st.pop();
	}
}

void p11()
{
	ifstream in("graf.in");
	int n, m, x, y;
	in >> n >> m;
	vector <vector<int>> adiacenta(n + 1);
	for (int i = 1; i <= m; ++i)
	{
		in >> x >> y;
		adiacenta[x].push_back(y);
	}
	Graf g(n, m, adiacenta);
	g.topologic();
}



class Disjoint
{
private:
	int n, m;
	int rang[LIMIT], par[LIMIT];
	
	int cautR(int x);
	void unite(int x, int y);
public:

	void afis();
};

int Disjoint::cautR(int x)
{
	while (x != par[x])
	{
		x = par[x];
	}
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
	ifstream in("disjoint.in");
	ofstream out("disjoint.out");
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
			unite(x, y);
			/*
			cout << "pentru1 : \n";
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
			*/
		}
		else
		{
			//cout << "pt2 : \n";
			//cout << "rad x = " << cautR(x) << "//// rad y = " << cautR(y);
			//cout << "\n";
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


void p12()   // havel hakimi
{
	ifstream in("graf.in");
	ofstream out("graf.out");
	int n, ok = -1, val, x;
	in >> n;
	vector <int> havelH(n);
	for (int i = 0; i < n; ++i)
	{
		in >> x;
		havelH.push_back(x);
	}
	while (ok == -1)
	{
		sort(havelH.begin(), havelH.end(), greater<>());

		if (havelH[0] == 0)
		{
			ok = 1;
			break;
		}
		val = havelH[0];
		havelH.erase(havelH.begin() + 0);

		if (val > havelH.size())
		{
			ok = 0;
			break;
		}

		for (int i = 0; i < val; ++i)
		{
			havelH[i] --;
			if (havelH[i] < 0)
				ok = 0;
		}
	}
	if (ok == 1)
		out << "Da";
	else
		out << "NU";
}  
// havel hakimi


int main()
{
	//p1();
	//p2();
	//p3();
	//p4();
	//p5();
	//p6();
	//p7();
	//p8();
	//p9();
	//p10();
	//p11();
	//p12();
}










