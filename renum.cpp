#define _CRT_SECURE_NO_DEPRECATE
#include "renum.hpp"

renum::renum(){};
renum::~renum(){};

void renum::convert(const int nc, const int nv, const vector<int>& ia, const vector<int>& ja, vector<int>& rcm)
{
	// nc, nv - number of cells, number of vertices
	// ia -array of positions
	// ja -array of nodes

	//working arrays
	int * order = new int[nv]; //number of connections for each node
	for(int k = 0; k < nv; ++k) 
	{
		rcm[k] = -1; //no enumeration provided yet
		order[k] = 0; //no connections yet
	}
	//compute order of each vertex - number of cells that are connected to it
	//std::vector< std::vector<int> > v_ja; // connections from each node to surrounding cells
	for(int k = 0; k < nc; ++k)
	{
		for(int j = ia[k]; j < ia[k+1]; ++j)
		{
			order[ja[j]]++;
			//v_ja[ ja[j] ].push_back(k); // vertex ja[j] is connected to cell k
		}
	}
	int * v_ia = new int[nv+1]; //array of positions for vertices
	int * v_ja = new int[ia[nc]-ia[0]]; //array of connections for vertices
	int * tmp = new int[nv]; //temporal array of positions of current full
	v_ia[0] = 0;
	for(int k = 0; k < nv; ++k)
	{
		tmp[k] = v_ia[k];
		v_ia[k+1] = v_ia[k] + order[k];
	}
	for(int k = 0; k < nc; ++k)
	{
		for(int j = ia[k]; j < ia[k+1]; ++j)
		{
			v_ja[ tmp[ ja[j] ]++ ] = k; 
		}
	}
	delete [] tmp;
	int id = 0; // starting id
	int cur; //selected cell
	do
	{
		cur = -1;
		//find first non-enumerated vertex
		for(int k = 0; k < nv && cur == -1; ++k)
			if( rcm[k] == -1 ) cur = k;
		//find vertex with smallest order
		for(int k = cur; k < nv; ++k)
			if( order[k] < order[cur] ) cur = k;
		assert(cur != -1); //just in case
		//initialize queue
		std::deque<int> q;
		q.push_back(cur);
		//enumerate selected vertex
		rcm[cur] = id++;
		//array to manage connections
		std::vector<int> conns;
		while(!q.empty())
		{
			//retive next node
			cur = q.front();
			q.pop_front();
			//get it's neighbours
			for(int k = v_ia[cur]; k < v_ia[cur+1]; ++k) //go over node's cells
				for(int j = ia[v_ja[k]]; j < ia[v_ja[k]+1]; ++j) //go over k-th cell's nodes
					if( rcm[ja[j]] == -1 ) conns.push_back(ja[j]);
			//connections may be not unique
			std::sort(conns.begin(),conns.end());
			conns.resize(std::unique(conns.begin(),conns.end())-conns.begin());
			//sort according to number of connections
			std::sort(conns.begin(),conns.end(),order_comparator(order));
			//go over connections enumerate them and schedule into queue
			for(int k = 0; k < conns.size(); ++k)
			{
				rcm[conns[k]] = id++; 
				q.push_back(conns[k]);
			}
			//remove any neighbours
			conns.clear();
		}

	} while(id < nv);
	//reverse enumeration
	for(int k = 0; k < nv; ++k)
		rcm[k] = nv - rcm[k] - 1;	
	//write order info
	FILE * fout = fopen("rcm.txt","w");
	for(int k = 0; k < nv; ++k)
		fprintf(fout,"%d %d\n",k, rcm[k]);
	fclose(fout);
	//write matrix in mtx format
	{
		FILE * mtx = fopen("matrix.mtx","w");
		fprintf(mtx,"%%MatrixMarket matrix coordinate real general\n");
		fprintf(mtx,"%d %d %d",nv,nv,ia[nc]-ia[0]);
		std::vector<int> conns;
		for(int k = 0; k < nv; ++k) //go over nodes
		{
			for(int j = v_ia[k]; j < v_ia[k+1]; ++j) //nodes to cells
				for(int l = ia[v_ja[j]]; l < ia[v_ja[j]+1]; ++l) //cells to nodes
					conns.push_back(ja[l]);
			//connections may be not unique
			std::sort(conns.begin(),conns.end());
			conns.resize(std::unique(conns.begin(),conns.end())-conns.begin());
			//write connections
			for(int j = 0; j < (int)conns.size(); ++j)
				fprintf(mtx,"%d %d 1.0\n",rcm[k],rcm[conns[j]]);
			conns.clear();
		}
		fclose(mtx);
	}
	//write matrix in mtx format
	{
		FILE * mtx = fopen("original_matrix.mtx","w");
		fprintf(mtx,"%%MatrixMarket matrix coordinate real general\n");
		fprintf(mtx,"%d %d %d",nv,nv,ia[nc]-ia[0]);
		std::vector<int> conns;
		for(int k = 0; k < nv; ++k) //go over nodes
		{
			for(int j = v_ia[k]; j < v_ia[k+1]; ++j) //nodes to cells
				for(int l = ia[v_ja[j]]; l < ia[v_ja[j]+1]; ++l) //cells to nodes
					conns.push_back(ja[l]);
			//connections may be not unique
			std::sort(conns.begin(),conns.end());
			conns.resize(std::unique(conns.begin(),conns.end())-conns.begin());
			//write connections
			for(int j = 0; j < (int)conns.size(); ++j)
				fprintf(mtx,"%d %d 1.0\n",k,conns[j]);
			conns.clear();
		}
		fclose(mtx);
	}
	delete [] order;
	delete [] v_ia;
	delete [] v_ja;
}