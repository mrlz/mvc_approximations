#include <iostream>
#include <list>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <chrono>
#include <thread> //delete this unnecesary dependency
#include <random>
#include <omp.h>

struct max_heap_node{
    int degree;
    int current_position;
    int vertex;
};

class max_heap{
private:
  std::vector<max_heap_node> vertices;
  int heap_size;
  std::vector<int> indexes;
public:
  max_heap();
  max_heap(std::vector<max_heap_node> array, int size, std::vector<int> indexing);
  void max_heapify(int i);
  int left(int i);
  int right(int i);
  max_heap_node get_max();
  void replace_max(max_heap_node x);
  void change_vertex_degree(int vertex, int degree);
};

class graph{
  int V;
  long long E;
  long long E_state;
  std::vector<bool> active_vertices;
  std::vector<bool> active_vertices_state;

  std::vector<int> position_in_heap;
  std::vector<int> position_in_heap_state;
  std::vector<int> adjacency_list_sizes;
  std::vector<int> adjacency_list_sizes_state;

  max_heap heap;
  max_heap heap_state;

  std::vector<max_heap_node> heap_vertices;
  std::vector<max_heap_node> heap_vertices_state;

  std::vector<omp_lock_t> mutexes;

  omp_lock_t edge_mutex;

  std::vector<std::vector <int> > adjacency_list;

  std::uniform_int_distribution<int> uniform_distribution;
  const int start = 0;
  int end;
public:
  graph(int vertices);
  graph(std::string sourcefile);
  void add_edge(int v1, int v2);
  int vertex_degree(int vertex);
  bool contains_edge(int v1, int v2);
  long long e_size();
  int v_size();
  void remove_neighbour_from_all(int neighbour);
  void remove_neighbour_from_all_no_heap(int neighbour);
  int get_random_number();
  int pick_random_active_vertex();
  int random_neighbour(int vertex);
  int highest_degree_neighbour(int vertex);
  int get_highest_degree_vertex();
  void add_edge_safely(int v1, int v2);
  void build_heap();
  void destroy_locks();
  void copy_state();
  void restore_state();
  void shrink_adjacency_list(int vertex);
  void allocate_adjacency_list_size(int vertex, int neighbours);
};

void print_list(std::list<int> l);
void format_files(std::string sourcefile);
double elapsed_time_milli(std::chrono::time_point<std::chrono::steady_clock,std::chrono::nanoseconds> start, std::chrono::time_point<std::chrono::steady_clock,std::chrono::nanoseconds> end);
double elapsed_time_seconds(std::chrono::time_point<std::chrono::steady_clock,std::chrono::nanoseconds> start, std::chrono::time_point<std::chrono::steady_clock,std::chrono::nanoseconds> end);

std::list<int> mvc_2_approx(graph *G, double* elapsed_time);

std::list<int> mvc_highest_degree(graph *G, double* elapsed_time);

std::list<int> mvc_2_approx_upgrade(graph *G, double* elapsed_time);
