#include "mvc.h"

std::random_device rand_dev; //these will produce different sequences everytime
std::mt19937 generator(rand_dev());

double elapsed_time_milli(std::chrono::time_point<std::chrono::steady_clock,std::chrono::nanoseconds> start, std::chrono::time_point<std::chrono::steady_clock,std::chrono::nanoseconds> end){
std::chrono::duration<double, std::milli> fp_ms = end-start;
return fp_ms.count();
}

double elapsed_time_seconds(std::chrono::time_point<std::chrono::steady_clock,std::chrono::nanoseconds> start, std::chrono::time_point<std::chrono::steady_clock,std::chrono::nanoseconds> end){
double millis = elapsed_time_milli(start,end);
return millis/1000.0;
}

void print_list(std::list<int> l){
  std::cout << "[";
  for(int j : l){
    std::cout << j << ", ";
  }
  std::cout << "]" << std::endl;
}

void print_vector(std::vector<int> v, std::vector<bool> status){
  std::cout << "[";
  for(int j : v){
    std::cout << j << ": " << status[j] << ", ";
  }
  std::cout << "]" << std::endl;
}

void format_files(std::string sourcefile){
  std::ofstream outfile((sourcefile+".formatted").c_str());
  std::ifstream infile(sourcefile.c_str());
  std::string input;

  std::string c,d;
  int a,b;
  while(std::getline(infile,input)){
    if (input.substr(0,1) == "p"){
      std::stringstream stream(input);
      stream >> c >> d >> a >> b;
      outfile << a << " " << b << std::endl;
    }else if(input.substr(0,1) == "e"){
      std::stringstream stream(input);
      stream >> c >> a >> b;
      outfile << a-1 << " " << b-1 << std::endl;
    }
  }
  outfile.close();
  infile.close();
  }

void swap_vertices(max_heap_node* x, max_heap_node* y, int *indexing){
    max_heap_node temp = *x;
    indexing[x->vertex] = y->current_position;
    indexing[y->vertex] = x->current_position;
    int current_pos_x = x->current_position;
    int current_pos_y = y->current_position;
    *x = *y;
    x->current_position = current_pos_x;
    *y = temp;
    y->current_position = current_pos_y;
}

max_heap::max_heap(){
}


max_heap::max_heap(std::vector<max_heap_node> array, int size, std::vector<int> indexing){
    heap_size = size;
    vertices = array;
    indexes = indexing;

    int i = (heap_size - 1) / 2;
    while (i >= 0){
        max_heapify(i);
        i--;
    }
}

void max_heap::max_heapify(int i){
    int l = left(i);
    int r = right(i);
    int biggest = i;

      if(l < heap_size && vertices[l].degree > vertices[i].degree){
        biggest = l;
    }

      if(r < heap_size && vertices[r].degree > vertices[biggest].degree){
        biggest = r;
    }

    if (biggest != i){
        swap_vertices(&vertices[i], &vertices[biggest], indexes.data());
        max_heapify(biggest);
    }
}

int max_heap::left(int i) {
return (2 * i + 1);
}

int max_heap::right(int i) {
return (2 * i + 2);
}

max_heap_node max_heap::get_max() {
return vertices[0];
}

void max_heap::replace_max(max_heap_node x){
    vertices[0] = x;
    max_heapify(0);
}

void max_heap::change_vertex_degree(int vertex, int degree){
  this->vertices[this->indexes[vertex]].degree = degree;
  max_heapify(this->indexes[vertex]);
}

graph::graph(int vertices){
  this->V = vertices;
  this->E = 0;
  this->mutexes = std::vector<omp_lock_t>(vertices);
  int i = 0;
  while(i < vertices){
    omp_init_lock(&(this->mutexes[i]));
    i++;
  }
  omp_init_lock(&(this->edge_mutex));
  this->adjacency_list = std::vector<std::vector <int> >(vertices);
  this->adjacency_list_sizes = std::vector<int>(vertices,0);
  this->active_vertices = std::vector<bool>(vertices,true);
  this->end = vertices-1;
  this->uniform_distribution = std::uniform_int_distribution<int>(start,end);
  this->heap_vertices = std::vector<max_heap_node>(vertices);
  this->position_in_heap = std::vector<int>(vertices);
}

graph::graph(std::string sourcefile){
  std::ifstream infile(sourcefile.c_str());
  int vertices, edges;
  int v1, v2;
  infile >> vertices >> edges;
  this->V = vertices;

  this->E = 0;
  this->adjacency_list = std::vector< std::vector <int> >(vertices);
  this->heap_vertices = std::vector<max_heap_node>(vertices);

  this->position_in_heap = std::vector<int>(vertices);
  this->adjacency_list_sizes = std::vector<int>(vertices,0);
  this->active_vertices = std::vector<bool>(vertices,true);

  while(infile >> v1 >> v2){
    add_edge(v1,v2);
  }

  build_heap();

  this->end = vertices-1;
  this->uniform_distribution = std::uniform_int_distribution<int>(start,end);
  infile.close();
}

void graph::build_heap(){
  int i = 0;
  while(i < this->V){
    this->heap_vertices[i].vertex = i;
    this->heap_vertices[i].current_position = i;
    this->heap_vertices[i].degree = this->adjacency_list_sizes[i];
    this->position_in_heap[i] = i;
    i++;
  }
  this->heap = max_heap(heap_vertices, this->V, position_in_heap);
}

int graph::v_size(){
  return this->V;
}

long long graph::e_size(){
  return this->E;
}

void graph::add_edge_safely(int v1, int v2){
  omp_set_lock(&(this->mutexes[v1]));
  this->adjacency_list[v1].push_back(v2);
  this->adjacency_list_sizes[v1]++;
  omp_unset_lock(&(this->mutexes[v1]));

  omp_set_lock(&(this->mutexes[v2]));
  this->adjacency_list[v2].push_back(v1);
  this->adjacency_list_sizes[v2]++;
  omp_unset_lock(&(this->mutexes[v2]));

  omp_set_lock(&(this->edge_mutex));
  this->E = this->E + 1;
  omp_unset_lock(&(this->edge_mutex));
}

void graph::add_edge(int v1, int v2){
  if (contains_edge(v1,v2)){
    return;
  }
  this->adjacency_list[v1].push_back(v2); //code is repeated to save a function call, thus decreasing total construction time.
  this->adjacency_list_sizes[v1]++;
  this->adjacency_list[v2].push_back(v1);
  this->adjacency_list_sizes[v2]++;
  this->E = this->E + 1;
}

void graph::destroy_locks(){
  for (int i = 0; i < this->V ; i++){
    omp_destroy_lock(&(this->mutexes[i]));
  }
  omp_destroy_lock(&this->edge_mutex);
}

void graph::copy_state(){
  this->active_vertices_state = this->active_vertices;

  this->position_in_heap_state = this->position_in_heap;

  this->adjacency_list_sizes_state = this->adjacency_list_sizes;

  this->heap_state = this->heap;

  this->heap_vertices_state = this->heap_vertices;

  this->E_state = this->E;
}

void graph::restore_state(){
  this->active_vertices = this->active_vertices_state;

  this->position_in_heap = this->position_in_heap_state;

  this->adjacency_list_sizes = this->adjacency_list_sizes_state;

  this->heap = this->heap_state;

  this->heap_vertices = this->heap_vertices_state;

  this->E = this->E_state;
}

void graph::allocate_adjacency_list_size(int vertex, int neighbours){
  this->adjacency_list[vertex].reserve(neighbours);
}

void graph::shrink_adjacency_list(int vertex){
  omp_set_lock(&(this->mutexes[vertex]));
  this->adjacency_list[vertex].shrink_to_fit();
  omp_unset_lock(&(this->mutexes[vertex]));
}


void graph::remove_neighbour_from_all(int neighbour){
  int i = 0;

  this->active_vertices[neighbour] = false;
  this->E = this->E - this->adjacency_list_sizes[neighbour];

  //rebalance heap
  this->heap.change_vertex_degree(neighbour,0);

  for (int k : this->adjacency_list[neighbour]){ //this cannot be easily parallelized, since we must update the heap sequentially
    if (this->active_vertices[k]){

    this->adjacency_list_sizes[k] = this->adjacency_list_sizes[k]-1;
    if (this->adjacency_list_sizes[k] == 0){
      this->active_vertices[k] = false;
    }

    this->heap.change_vertex_degree(k,this->adjacency_list_sizes[k]);
    }
  }
}

void graph::remove_neighbour_from_all_no_heap(int neighbour){
  this->active_vertices[neighbour] = false;
  this->E = this->E - this->adjacency_list_sizes[neighbour];

  for (int k : this->adjacency_list[neighbour]){ //this is ripe for parallelization, since no neighbour is featured twice in the list
    if (this->active_vertices[k]){

    this->adjacency_list_sizes[k] = this->adjacency_list_sizes[k]-1;
    if (this->adjacency_list_sizes[k] == 0){
      this->active_vertices[k] = false;
    }

    }
  }

}

int graph::get_random_number(){
  return this->uniform_distribution(generator);
}

int graph::pick_random_active_vertex(){
int number = get_random_number();
  while(this->active_vertices[number]==false){
  number = (number+1)%this->V;
  }
return number;
}

bool graph::contains_edge(int v1, int v2){
  std::vector<int>::iterator finder = std::find(this->adjacency_list[v1].begin(), this->adjacency_list[v1].end(), v2);
  if (finder == this->adjacency_list[v1].end()){ // unexpected behaviour when dereferencing end() iterator, so we check
    return false;
  }
  return true;
}

int graph::get_highest_degree_vertex(){
  return this->heap.get_max().vertex;
}

int graph::highest_degree_neighbour(int vertex){
  int max_degree = 0;
  int max_v = 0;
  for (int v : this->adjacency_list[vertex]){ //thread safe since the same vertex won't be twice in the list
    if (this->active_vertices[v]){ //therefore there won't be a race condition to check if it is active
      int size = adjacency_list_sizes[v];
      if (size > max_degree){ //here we will have a race condition, therefore we need a mutex
        max_degree = size;
        max_v = v;
      }
    }
  }
  return max_v;
}

int graph::random_neighbour(int vertex){
  int chosen = this->get_random_number()%adjacency_list[vertex].size();
  while(this->active_vertices[(this->adjacency_list[vertex])[chosen]] == false){
    chosen = (chosen+1)%adjacency_list[vertex].size();
  }
  return (this->adjacency_list[vertex])[chosen];
}

std::list<int> mvc_2_approx(graph *G, double* elapsed_time){
  std::list<int> C;
  auto start = std::chrono::steady_clock::now();
  while(G->e_size() > 0){
    int u = G->pick_random_active_vertex();
    int v = G->random_neighbour(u);
    C.push_back(u);
    C.push_back(v);

    G->remove_neighbour_from_all_no_heap(u);
    G->remove_neighbour_from_all_no_heap(v);
  }
  auto end = std::chrono::steady_clock::now();
  *elapsed_time = elapsed_time_seconds(start,end);
  return C;
}

std::list<int> mvc_highest_degree(graph *G, double* elapsed_time){
  std::list<int> C;
  auto start = std::chrono::steady_clock::now();
  while(G->e_size() > 0){
    int u = G->get_highest_degree_vertex();
    C.push_back(u);
    G->remove_neighbour_from_all(u);
  }
  auto end = std::chrono::steady_clock::now();
  *elapsed_time = elapsed_time_seconds(start,end);
  return C;
}

std::list<int> mvc_2_approx_upgrade(graph *G, double* elapsed_time){
  std::list<int> C;
  auto start = std::chrono::steady_clock::now();
  while(G->e_size() > 0){
    int u = G->get_highest_degree_vertex();
    int v = G->highest_degree_neighbour(u);
    C.push_back(u);
    C.push_back(v);
    G->remove_neighbour_from_all(u);
    G->remove_neighbour_from_all(v);
  }
  auto end = std::chrono::steady_clock::now();
  *elapsed_time = elapsed_time_seconds(start,end);
  return C;
}
