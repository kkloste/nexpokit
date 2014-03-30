#include<iostream>
#include<vector>
#include<unordered_map>
using namespace std;

struct sparse_heap
{
	int type;  // 0 -- min-heap; nonzero -- max-heap
	unordered_map<int, double> values;  // v = values.find(vid)->second
	unordered_map<int, int> lookup;     // tid = lookup.find(vid)->second
	vector<int> heap;                   // vid = heap[tid]
	
	sparse_heap();                      // default min_heap
	sparse_heap(int t);                 // t = 0, min_heap; t != 0, max_heap
	void insert(int vid, double val);
	int pop_top(void);
	void update(int vid, double val);
	int peek_top(void);
};

// default constructor
sparse_heap::sparse_heap(): type(0){}
// set heap type
sparse_heap::sparse_heap(int t)
{
	type = t;
}

void sparse_heap::insert(int vid, double val)
{
	values.insert(std::make_pair(vid, val));
	//n++;
	heap.push_back(vid);
	int tid = heap.size() - 1;
	int ptid = -1;
	//heap.push_back(vid);
	heap[tid] = vid;
	lookup.insert(std::make_pair(vid, tid));
	while(1)
	{
		if (tid == 0)
		{
			return;
		}
		ptid = (tid - 1)/2;
		double tmp1, tmp2;
		if (type == 0){ // min-heap
			tmp1 = values.find(heap[ptid])->second;
			tmp2 = values.find(heap[tid])->second;		
		}
		else{  // max-heap
			tmp1 = -values.find(heap[ptid])->second;
			tmp2 = -values.find(heap[tid])->second;		
		}
		if (tmp1 < tmp2)
		{
			return;
		}
		else{
			int vid_ptid = heap[ptid];
			heap[ptid] = heap[tid];
			lookup.find(heap[tid])->second = ptid;
			heap[tid] = vid_ptid;
			lookup.find(vid_ptid)->second = tid;
			tid = ptid;
		}
	}
}

int sparse_heap::peek_top(void)
{
	if (heap.size() == 0){
		return -1;	
	}
	return heap[0];
}

int sparse_heap::pop_top(void)
{
	if (heap.size() == 1){
		//n --;
		lookup.erase(heap[0]);
		values.erase(heap[0]);
		int rootid = heap[0];
		heap.pop_back();
		return rootid; // return root's vid	
	}
	int vid = heap[0];
	heap[0] = heap[heap.size() - 1]; // put the last element on the root and push down
	heap[heap.size() - 1] = vid;
	lookup.erase(vid);
	values.erase(vid);
	lookup.find(heap[0])->second = 0;
	heap.pop_back();
	// push down the root to a proper position
	int k = 0;
	while(1){
		int i = 2*k+1;
		int kvid = heap[k];
		double tmp1, tmp2;
		if (i > heap.size() - 1){
			return vid;		
		}
		if (i == heap.size() - 1){
			;		
		}
		else{
			int lcvid = heap[i];
			int rcvid = heap[i+1];
			if (type == 0){  // min-heap
				tmp1 = values.find(rcvid)->second;
				tmp2 = values.find(lcvid)->second;
			}
			else{  //max-heap
				tmp1 = -values.find(rcvid)->second;
				tmp2 = -values.find(lcvid)->second;			
			}
			if (tmp1 < tmp2){
				i = i + 1;
			}		
		}
		if (type == 0){ // min-heap
			tmp1 = values.find(kvid)->second;
			tmp2 = values.find(heap[i])->second;
		}
		else{
			tmp1 = -values.find(kvid)->second;
			tmp2 = -values.find(heap[i])->second;
		}
		if (tmp1 < tmp2){
			return vid;		
		}
		else{
			heap[k] = heap[i];
			lookup.find(heap[i])->second = k;
			heap[i] = kvid;
			lookup.find(kvid)->second = i;
			k = i;		
		}	
	}
}

void sparse_heap::update(int vid, double val){
	if (lookup.find(vid) == lookup.end()){
		// cannot find the element
		return;	
	}
	values.find(vid)->second = val;
	int ktid = lookup.find(vid)->second;
	double tmp1, tmp2;
	// push up if the new value is smaller than the old value
	while(1){
		int i = 2*ktid+1;
		int kvid = heap[ktid];
		
		if (i > heap.size() - 1){
			break;		
		}
		if (i == heap.size() - 1){
			;
		}
		else{
			int lcvid = heap[i];
			int rcvid = heap[i+1];
			if (type == 0){ // min-heap
				tmp1 = values.find(rcvid)->second;
				tmp2 = values.find(lcvid)->second;
			}
			else{  // max-heap
				tmp1 = -values.find(rcvid)->second;
				tmp2 = -values.find(lcvid)->second;
			}
			if (tmp1 < tmp2){
				i = i + 1;			
			}		
		}
		if (type == 0){ // min-heap
			tmp1 = values.find(kvid)->second;
			tmp2 = values.find(heap[i])->second;
		}
		else{	// max-heap
			tmp1 = -values.find(kvid)->second;
			tmp2 = -values.find(heap[i])->second;
		}
		if (tmp1 < tmp2){
			break;		
		}
		else{
			heap[ktid] = heap[i];
			lookup.find(heap[i])->second = ktid;
			heap[i] = kvid;
			lookup.find(kvid)->second = i;
			ktid = i;		
		}	
	}
	int jtid = ktid;
	int jvid = heap[jtid];
	int ptid = -1;
	int pvid = -1;
	// push down if the new value is larger than the old value
	while(1){
		if (jtid == 0){
			return;		
		}	
		ptid = (jtid - 1)/2;
		pvid = heap[ptid];
		if (type == 0){ // min-heap
			tmp1 = values.find(pvid)->second;
			tmp2 = values.find(jvid)->second;
		}
		else{  // max-heap
			tmp1 = -values.find(pvid)->second;
			tmp2 = -values.find(jvid)->second;
		}
		if (tmp1 < tmp2){
			return;		
		}
		else{
			heap[ptid] = jvid;
			lookup.find(jvid)->second = ptid;
			heap[jtid] = pvid;
			lookup.find(pvid)->second = jtid;
			jtid = ptid;		
		}
	}
}

void showInfo(sparse_heap &heap)
{
	cout << "Printing heap..." << endl << "[";
	for (int i = 0; i < heap.heap.size(); i++){
		if (i != heap.heap.size()-1){
			cout << heap.heap[i] << ",";		
		}
		else{
			cout << heap.heap[i] << "]" << endl;		
		}	
	}
	cout << "Printing values..." << endl;
	for (int i = 0; i < heap.values.bucket_count(); i ++){
		for(auto it = heap.values.begin(i); it != heap.values.end(i); it ++){
			cout << "[" << it->first << ":" << it->second << "]" << endl;		
		}	
	}
	cout << "Printing lookup..." << endl;
	for (int i = 0; i < heap.lookup.bucket_count(); i ++){
		for (auto it = heap.lookup.begin(i); it != heap.lookup.end(i); it ++){
			cout << "[" << it->first << ":" << it->second << "]" << endl;		
		}	
	}
}
int main()
{
	sparse_heap heap = sparse_heap(1);
	heap.insert(1, 10.5);
	heap.insert(2, 9.8);
	heap.insert(3, 10.0);
	heap.insert(4, 2.5);
	heap.insert(5, 50);
	cout << "Inserting..." << endl;
	showInfo(heap);
	cout << "Popping..." << endl;
	heap.pop_top();
	showInfo(heap);
	cout << "Updating..." << endl;
	heap.update(3, 100.8);
	showInfo(heap);
	return 0;
}
