#include "mvc.h"

int approx_average_adj_list_size(int vertices, double p){
  return ((vertices-1)*p + vertices*0.02);
}

void generate_random_graph(int vertices, double p, graph *G){
  int u = 0;
  int v = 0;
  int thread_number = 4;
   //these will produce different sequences everytime, for each separate thread
  std::vector<std::random_device> devices = std::vector<std::random_device>(thread_number);
  std::vector<std::mt19937> generators = std::vector<std::mt19937>(thread_number);
  std::vector<std::bernoulli_distribution> distributions = std::vector<std::bernoulli_distribution>(thread_number);
  for(int b = 0; b < 4; b++){
    distributions[b] = std::bernoulli_distribution(p);
  }

  int avg_neighbours = approx_average_adj_list_size(vertices, p);

  for(int v = 0; v < vertices; v++){
    G->allocate_adjacency_list_size(v, avg_neighbours);
  }
  for (int u = 0; u < vertices; u++){
    #pragma omp parallel for
    for(int v = u+1; v < vertices; v++){
        if ( distributions[omp_get_thread_num()](generators[omp_get_thread_num()]) ){
          G->add_edge_safely(u,v);
        }
    }
  }
  G->build_heap();
}

long long mvc_approximation_iteration_file(std::string filename, int *vertices, double *two_approx_time, double *highest_time, double *upgrade_time, long long *two_approx_sol, long long *highest_approx_sol, long long *upgrade_approx_sol){
  graph G = graph(filename);
  G.copy_state();
  long long edges = G.e_size();
  *vertices = G.v_size();
  double total_time;

  *two_approx_sol = mvc_2_approx(&G,&total_time).size();
  *two_approx_time = total_time;
  G.restore_state();

  *highest_approx_sol = mvc_highest_degree(&G,&total_time).size();
  *highest_time = total_time;
  G.restore_state();

  *upgrade_approx_sol = mvc_2_approx_upgrade(&G,&total_time).size();
  *upgrade_time = total_time;

  G.destroy_locks();
  return edges;
}

long long mvc_approximation_iteration_random(int i, double p, double *two_approx_time, double *highest_time, double *upgrade_time, long long *two_approx_sol, long long *highest_approx_sol, long long *upgrade_approx_sol){
  int vertices = pow(2,i);
  graph G(vertices);
  generate_random_graph(vertices,p,&G);
  G.copy_state();
  long long edges = G.e_size();

  double total_time;

  *two_approx_sol = *two_approx_sol + mvc_2_approx(&G,&total_time).size();
  *two_approx_time = *two_approx_time + total_time;

  G.restore_state();
  *highest_approx_sol = *highest_approx_sol + mvc_highest_degree(&G,&total_time).size();
  *highest_time = *highest_time + total_time;

  G.restore_state();
  *upgrade_approx_sol = *upgrade_approx_sol + mvc_2_approx_upgrade(&G,&total_time).size();
  *upgrade_time = *upgrade_time + total_time;

  G.destroy_locks();
  return edges;
}

void average_values(long long *edges, double *two_approx_time,double *highest_approx_time,double *upgrade_approx_time,long long *two_approx_sol_size,long long *highest_approx_sol_size,long long *upgrade_approx_sol_size, long long iterations){
  *two_approx_time = *two_approx_time/iterations;
  *highest_approx_time = *highest_approx_time/iterations;
  *upgrade_approx_time = *upgrade_approx_time/iterations;

  *two_approx_sol_size = *two_approx_sol_size/iterations;
  *highest_approx_sol_size = *highest_approx_sol_size/iterations;
  *upgrade_approx_sol_size = *upgrade_approx_sol_size/iterations;

  *edges = *edges/iterations;
}

void random_experiments(){
  int iterations = 4;
  double p = 0.18;

  std::string filename = "random_averages.csv";
  omp_lock_t writelock;
  omp_init_lock(&writelock);

  std::ofstream outfile(filename.c_str());
  outfile << "vertices,edges,p,i,2_time,2_C_size,high_time,high_C_size,up_time,up_C_size" << std::endl;
    for(int size = 10; size < 21; size++){
    auto start_size = std::chrono::steady_clock::now();
    std::cout << "i:" << size <<"{" << std::endl;
        for(int p_step = 1; p_step < 6; ++p_step){
          std::cout << "  p:" << p_step*p << std::endl;
        double two_approx_time = 0.0;
        double highest_approx_time = 0.0;
        double upgrade_approx_time = 0.0;

        long long two_approx_sol_size = 0;
        long long highest_approx_sol_size = 0;
        long long upgrade_approx_sol_size = 0;

        long long edges = 0;

        int it = 0;
        auto start_time = std::chrono::steady_clock::now();
        while(it < iterations){
          edges = edges + mvc_approximation_iteration_random(size,p*p_step,&two_approx_time, &highest_approx_time, &upgrade_approx_time, &two_approx_sol_size, &highest_approx_sol_size, &upgrade_approx_sol_size);
          it++;
        }
        auto end_time = std::chrono::steady_clock::now();
        average_values(&edges, &two_approx_time,&highest_approx_time,&upgrade_approx_time,&two_approx_sol_size,&highest_approx_sol_size,&upgrade_approx_sol_size, iterations);

        omp_set_lock(&writelock);
        outfile << pow(2,size) << "," << edges << "," << p*p_step << "," << size << "," << two_approx_time << "," << two_approx_sol_size << "," << highest_approx_time << "," << highest_approx_sol_size << "," << upgrade_approx_time << "," << upgrade_approx_sol_size << std::endl;
        omp_unset_lock(&writelock);
    }
    auto end_size = std::chrono::steady_clock::now();
    omp_set_lock(&writelock);
    std::cout << "} i:" << size << " time: " << elapsed_time_seconds(start_size, end_size) << std::endl;
    omp_unset_lock(&writelock);
  }

  omp_destroy_lock(&writelock);
  outfile.close();
}

void disk_experiments(){
  std::string log_file = "graphs_from_disk.csv";
  std::ofstream outfile(log_file.c_str());
  outfile << "name,vertices,edges,2_time,2_C_size,high_time,high_C_size,up_time,up_C_size" << std::endl;

  std::string folder = "./graphs/";
  std::vector<std::string> filenames = std::vector<std::string>(11);

  filenames[0] = folder+"planar_embedding1000000.pg";
  filenames[1] = folder+"planar_embedding5000000.pg";
  filenames[2] = folder+"planar_embedding10000000.pg";

  filenames[3] = folder+"C1000.9.clq";
  filenames[4] = folder+"C2000.5.col";
  filenames[5] = folder+"C2000.9.clq";
  filenames[6] = folder+"C4000.5.col";

  filenames[7] = folder+"frb35-17-1.mis";
  filenames[8] = folder+"frb50-23-1.mis";
  filenames[9] = folder+"frb59-26-1.mis";
  filenames[10] = folder+"frb100-40.mis";

  double two_approx_time = 0.0;
  double highest_approx_time = 0.0;
  double upgrade_approx_time = 0.0;

  long long two_approx_sol_size = 0;
  long long highest_approx_sol_size = 0;
  long long upgrade_approx_sol_size = 0;

  long long edges;
  int vertices;

  for (int k = 3; k < 11; k++){
    std::cout << "formatting file " << filenames[k] << std::endl;
    format_files(filenames[k]);
    filenames[k] = filenames[k]+".formatted";
  }

  for (int files = 0; files < 11; files++){
    std::cout << "processing file " << filenames[files] << std::endl;
    edges = mvc_approximation_iteration_file(filenames[files], &vertices, &two_approx_time, &highest_approx_time, &upgrade_approx_time, &two_approx_sol_size, &highest_approx_sol_size, &upgrade_approx_sol_size);
    outfile << filenames[files] << ","<< vertices << "," << edges <<  "," << two_approx_time << "," << two_approx_sol_size << "," << highest_approx_time << "," << highest_approx_sol_size << "," << upgrade_approx_time << "," << upgrade_approx_sol_size << std::endl;
  }
}

int main(){

  int experiments_to_run = 1;

  if(experiments_to_run){
    random_experiments();
  }else{
    disk_experiments();
  }

  return 0;
  }
