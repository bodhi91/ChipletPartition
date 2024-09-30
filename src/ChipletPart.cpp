#pragma once
#include "ChipletPart.h"
#include "Hypergraph.h"
#include "Utilities.h"
#include "pugixml.hpp"
#include <Python.h>
#include <omp.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <chrono>
#include <codecvt>
#include <filesystem>
#include <future> // For parallel fitness evaluation
#include <iomanip>
#include <iostream>
#include <locale>
#include <metis.h>
#include <numpy/arrayobject.h>
#include <queue>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// write the main here (accept arguments from command line)
namespace chiplet {

std::vector<int> ChipletPart::CrossBarExpansion(std::vector<int> &crossbars,
                                                int &num_parts) {
  if (crossbars.size() < num_parts) {
    return std::vector<int>();
  }

  // round robin assignment of crossbars to partitions
  std::vector<std::unordered_set<int>> partitions(num_parts);
  std::vector<std::queue<int>> queues(num_parts);
  std::unordered_map<int, int> node_partition_map;
  std::vector<std::unordered_map<int, int>> edge_count(
      hypergraph_->GetNumVertices());

  // Initialize partitions and queues with crossbars
  for (int i = 0; i < num_parts; i++) {
    int partition = i % num_parts;
    int crossbar = crossbars[i];
    partitions[partition].insert(crossbar);
    queues[partition].push(crossbar);
    node_partition_map[crossbar] = partition;
  }

  auto ShouldAddToPartition =
      [](const std::unordered_map<int, int> &edge_counts,
         const std::unordered_map<int, int> &node_partition_map,
         int node) -> bool {
    if (node_partition_map.find(node) != node_partition_map.end()) {
      return false;
    }

    int maxEdges = 0;
    int maxPartition = -1;

    // Find partition with the maximum edge count to this node
    for (const auto &pair : edge_counts) {
      if (pair.second > maxEdges) {
        maxEdges = pair.second;
        maxPartition = pair.first;
      }
    }

    // Only add the node to this partition if it has the majority of its edges
    // here
    return edge_counts.at(maxPartition) > maxEdges / 2;
  };

  auto AllQueuesEmpty = [](const std::vector<std::queue<int>> &queues) -> bool {
    return all_of(queues.begin(), queues.end(),
                  [](const std::queue<int> &q) { return q.empty(); });
  };

  // BFS expansion with weight-based partitioning
  while (!AllQueuesEmpty(queues)) {
    for (int i = 0; i < num_parts; ++i) {
      if (!queues[i].empty()) {
        int current_node = queues[i].front();
        queues[i].pop();

        // get the neighbors of the current node
        std::vector<int> neighbors = hypergraph_->GetNeighbors(current_node);
        for (int neighbor : neighbors) {
          if (partitions[i].count(neighbor) == 0) {
            edge_count[neighbor][i]++;
            if (ShouldAddToPartition(edge_count[neighbor], node_partition_map,
                                     neighbor)) {
              partitions[i].insert(neighbor);
              queues[i].push(neighbor);
              node_partition_map[neighbor] = i;
            }
          }
        }
      }
    }
  }

  std::vector<int> partition(hypergraph_->GetNumVertices());
  for (int i = 0; i < num_parts; i++) {
    for (int node : partitions[i]) {
      partition[node] = i;
    }
  }

  return partition;
}

std::vector<int> ChipletPart::SpectralPartition(const std::string &file_name) {
  // run the python script get_embedding.py
  std::vector<int> partition;
  Py_Initialize();
  // Initialize NumPy array API
  if (_import_array() < 0) {
    PyErr_Print();
    PyErr_SetString(PyExc_ImportError,
                    "numpy.core.multiarray failed to import");
    // Handle error: return an empty vector or throw an exception
    return std::vector<int>(); // Return an empty vector as an error indicator
  }
  if (!Py_IsInitialized()) {
    std::cerr << "Failed to initialize the Python interpreter." << std::endl;
    return partition;
  }

  PyObject *sysPath = PySys_GetObject("path");
  PyList_Append(sysPath,
                PyUnicode_FromString(path_to_embedding_script_.c_str()));
  PyObject *pName = PyUnicode_FromString("get_embedding");
  PyObject *pModule = PyImport_Import(pName);
  Py_DECREF(pName);

  if (!pModule) {
    PyErr_Print();
    std::cerr << "Failed to load the Python module." << std::endl;
    return partition;
  }

  PyObject *pFunc =
      PyObject_GetAttrString(pModule, "get_clusters_from_embedding");
  if (!pFunc || !PyCallable_Check(pFunc)) {
    PyErr_Print();
    std::cerr << "Cannot find function 'get_clusters_from_embedding'."
              << std::endl;
    Py_XDECREF(pFunc);
    Py_DECREF(pModule);
    return partition;
  }

  PyObject *pArgs = PyTuple_New(1);
  PyObject *pValue = PyUnicode_FromString(file_name.c_str());
  PyTuple_SetItem(
      pArgs, 0,
      pValue); // PyTuple_SetItem steals a reference, no need to decref pValue

  PyObject *pResult = PyObject_CallObject(pFunc, pArgs);
  Py_DECREF(pArgs);
  if (!pResult) {
    PyErr_Print();
    std::cerr << "Call to 'get_clusters_from_embedding' failed." << std::endl;
    Py_DECREF(pFunc);
    Py_DECREF(pModule);
    return partition;
  }

  // Ensure pResult is an array and check its type
  if (PyArray_Check(pResult) &&
      PyArray_TYPE((PyArrayObject *)pResult) == NPY_INT) {
    PyArrayObject *arr = reinterpret_cast<PyArrayObject *>(pResult);
    int *data = static_cast<int *>(PyArray_DATA(arr));
    int size = PyArray_SIZE(arr);

    // Copy data to a vector (if necessary)
    partition.assign(data, data + size);
  } else {
    std::cerr << "Returned object is not an int array." << std::endl;
  }

  // Clean up
  Py_DECREF(pResult);
  Py_DECREF(pFunc);
  Py_DECREF(pModule);

  return partition; // Return the vector, even if empty
}

std::vector<int> ChipletPart::KWayCuts(int &num_parts) {
  std::vector<int> partition(hypergraph_->GetNumVertices(), 0);
  std::vector<int> partition_sizes(num_parts, 0);
  int average_size = hypergraph_->GetNumVertices() / num_parts;
  // seed the random number generator
  srand(time(NULL));

  // Random initial partitioning
  for (int i = 0; i < hypergraph_->GetNumVertices(); i++) {
    partition[i] = rand() % num_parts;
    partition_sizes[partition[i]]++;
  }

  // Attempt to balance the partitions
  bool balanced = false;
  while (!balanced) {
    balanced = true;
    for (int i = 0; i < hypergraph_->GetNumVertices(); ++i) {
      int current_partition = partition[i];
      if (partition_sizes[current_partition] > average_size) {
        for (int j = 0; j < num_parts; ++j) {
          if (partition_sizes[j] < average_size) {
            partition[i] = j;
            partition_sizes[current_partition]--;
            partition_sizes[j]++;
            balanced = false;
            break;
          }
        }
      }
    }
  }

  return partition;
}

std::vector<int> ChipletPart::METISPart(int &num_parts) {
  // create the METIS file
  // generate random number
  std::string graph_file = "metis_graph.txt";
  std::ofstream metis_file(graph_file);
  metis_file << hypergraph_->GetNumVertices() << " "
             << hypergraph_->GetNumHyperedges() << std::endl;

  for (int i = 0; i < hypergraph_->GetNumVertices(); i++) {
    for (int j = 0; j < hypergraph_->GetNeighbors(i).size(); j++) {
      metis_file << hypergraph_->GetNeighbors(i)[j] + 1 << " ";
    }
    metis_file << std::endl;
  }

  metis_file.close();

  // run gpmetis to run generate the partition
  std::string command = "./gpmetis " + graph_file + " " +
                        std::to_string(num_parts) + " > /dev/null 2>&1";
  int status = system(command.c_str());

  if (WIFEXITED(status) && WEXITSTATUS(status) != 0) {
    std::cerr << "Error: gpmetis failed with status " << WEXITSTATUS(status)
              << std::endl;
    return std::vector<int>();
  } else {
    // read the partition file
    std::string partition_file =
        "metis_graph.txt.part." + std::to_string(num_parts);
    std::ifstream partition_input(partition_file);
    if (!partition_input.is_open()) {
      std::cerr << "Error: Cannot open partition file " << partition_file
                << std::endl;
      return std::vector<int>();
    }
    std::vector<int> partition(hypergraph_->GetNumVertices(), 0);
    for (int i = 0; i < hypergraph_->GetNumVertices(); i++) {
      partition_input >> partition[i];
    }
    // delete the graph file
    // command = "rm -f " + graph_file;
    // status = system(command.c_str());
    return partition;
  }
}

void ChipletPart::ReadChipletGraph(std::string hypergraph_file,
                                   std::string chiplet_io_file) {
  // read the chiplet_io_file_
  if (!chiplet_io_file.empty()) {
    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load_file(chiplet_io_file.c_str());
    if (!result) {
      std::cerr << "Error: " << result.description() << std::endl;
      return;
    }
    pugi::xml_node ios = doc.child("ios");
    for (pugi::xml_node io = ios.child("io"); io; io = io.next_sibling("io")) {
      // Read the 'type' attribute
      std::string type = io.attribute("type").as_string();

      // Read the 'reach' attribute
      double reach = io.attribute("reach").as_double();

      // Add to hash map
      io_map_[type] = reach;
    }
  }

  // read hypergraph file
  std::ifstream hypergraph_file_input(hypergraph_file);
  if (!hypergraph_file_input.is_open()) {
    std::cerr << "Error: Cannot open hypergraph file " << hypergraph_file
              << std::endl;
    return;
  }
  // Check the number of vertices, number of hyperedges, weight flag
  std::string cur_line;
  std::getline(hypergraph_file_input, cur_line);
  std::istringstream cur_line_buf(cur_line);
  std::vector<int> stats{std::istream_iterator<int>(cur_line_buf),
                         std::istream_iterator<int>()};
  num_hyperedges_ = stats[0];
  num_vertices_ = stats[1];
  bool hyperedge_weight_flag = false;
  bool vertex_weight_flag = false;
  if (stats.size() == 3) {
    if ((stats[2] % 10) == 1) {
      hyperedge_weight_flag = true;
    }
    if (stats[2] >= 10) {
      vertex_weight_flag = true;
    }
  }

  // clear the related vectors
  hyperedges_.clear();
  hyperedge_weights_.clear();
  vertex_weights_.clear();
  hyperedges_.reserve(num_hyperedges_);
  hyperedge_weights_.reserve(num_hyperedges_);
  vertex_weights_.reserve(num_vertices_);
  std::unordered_map<long long, std::vector<int>> distinct_hyperedges;

  // Read hyperedge information
  for (int i = 0; i < num_hyperedges_; i++) {
    std::getline(hypergraph_file_input, cur_line);
    if (hyperedge_weight_flag == true) {
      std::istringstream cur_line_buf(cur_line);
      std::vector<float> hvec{std::istream_iterator<float>(cur_line_buf),
                              std::istream_iterator<float>()};
      // each line is wt reach {vertex list}
      // the first element is the weight, the second element is the reach
      // the remaining elements are the vertex list
      std::vector<float> hwt = {hvec[0]};
      float reach = hvec[1];
      float io_size = hvec[2];
      io_sizes_.push_back(io_size);
      reach_.push_back(reach);
      std::vector<int> hyperedge(hvec.begin() + 3, hvec.end());
      long long hash = 0;
      for (auto &value : hyperedge) {
        value--; // the vertex id starts from 1 in the hypergraph file
        hash += value * value;
      }

      if (distinct_hyperedges.find(hash) != distinct_hyperedges.end()) {
        continue;
      } else {
        distinct_hyperedges[hash] = hyperedge;
      }

      hyperedge_weights_.push_back(hwt);
      hyperedges_.push_back(hyperedge);
    } else {
      std::istringstream cur_line_buf(cur_line);
      std::vector<int> hyperedge{std::istream_iterator<int>(cur_line_buf),
                                 std::istream_iterator<int>()};
      for (auto &value : hyperedge) {
        value--; // the vertex id starts from 1 in the hypergraph file
      }
      std::vector<float> hwts(hyperedge_dimensions_,
                              1.0); // each dimension has the same weight
      hyperedge_weights_.push_back(hwts);
      hyperedges_.push_back(hyperedge);
    }
  }

  // Read weight for vertices
  for (int i = 0; i < num_vertices_; i++) {
    if (vertex_weight_flag == true) {
      std::getline(hypergraph_file_input, cur_line);
      std::istringstream cur_line_buf(cur_line);
      std::vector<float> vwts{std::istream_iterator<float>(cur_line_buf),
                              std::istream_iterator<float>()};
      vertex_weights_.push_back(vwts);
    } else {
      std::vector<float> vwts(vertex_dimensions_, 1.0);
      vertex_weights_.push_back(vwts);
    }
  }

  hypergraph_ = std::make_shared<Hypergraph>(
      vertex_dimensions_, hyperedge_dimensions_, hyperedges_, vertex_weights_,
      hyperedge_weights_, reach_, io_sizes_);
  std::cout << "[INFO] Number of IP blocks in chiplet graph: " << num_vertices_
            << std::endl;
  std::cout << "[INFO] Number of nets in chiplet graph: " << hyperedges_.size()
            << std::endl;
}

void ChipletPart::TechAssignPartition(
    std::string hypergraph_file, std::string chiplet_io_file,
    std::string chiplet_layer_file, std::string chiplet_wafer_process_file,
    std::string chiplet_assembly_process_file, std::string chiplet_test_file,
    std::string chiplet_netlist_file, std::string chiplet_blocks_file,
    float reach, float separation, std::vector<std::string> techs) {
  // seed the rng_ with 0
  rng_.seed(seed_);
  auto start_time = std::chrono::high_resolution_clock::now();
  std::cout << "[INFO] Reading the chiplet graph " << hypergraph_file
            << std::endl;
  std::cout << "[INFO] Reading the chiplet IO file " << chiplet_io_file
            << std::endl;
  std::cout << "[INFO] Reading the chiplet layer file " << chiplet_layer_file
            << std::endl;
  std::cout << "[INFO] Reading the chiplet wafer process file "
            << chiplet_wafer_process_file << std::endl;
  std::cout << "[INFO] Reading the chiplet assembly process file "
            << chiplet_assembly_process_file << std::endl;
  std::cout << "[INFO] Reading the chiplet test file " << chiplet_test_file
            << std::endl;
  std::cout << "[INFO] Reading the chiplet netlist file "
            << chiplet_netlist_file << std::endl;
  std::cout << "[INFO] Reading the chiplet blocks file " << chiplet_blocks_file
            << std::endl;
  std::cout << "[INFO] Reach: " << reach << std::endl;
  std::cout << "[INFO] Separation: " << separation << std::endl;
  std::cout << "[INFO] Techs: ";
  for (auto &tech : techs) {
    std::cout << tech << " ";
  }
  std::cout << std::endl;

  // for each tech node run the partitioner first
  // 45nm, 10nm, 7nm
  std::vector<std::string> optimal_tech_arr;
  std::vector<int> optimal_parts_for_each_tech;
  const std::string separator(60, '-'); // Creates a separator line
  /*for (auto &tech : techs) {
    std::cout << separator << std::endl;
    std::cout << std::setw((separator.size() + 10) / 2) << "TECH: " << tech
              << std::endl;
    std::cout << separator << std::endl;
    Partition(hypergraph_file, chiplet_io_file, chiplet_layer_file,
              chiplet_wafer_process_file, chiplet_assembly_process_file,
              chiplet_test_file, chiplet_netlist_file, chiplet_blocks_file,
              reach, separation, tech);
    optimal_parts_for_each_tech.push_back(tech_parts_);
  }

  // replicate techs based on the number of partitions to estimate number of
  // partitions
  std::vector<std::string> techs_replicated;
  for (int i = 0; i < optimal_parts_for_each_tech.size(); i++) {
    for (int j = 0; j < optimal_parts_for_each_tech[i]; j++) {
      techs_replicated.push_back(techs[i]);
    }
  }

  int estimated_max_parts = techs_replicated.size();*/
  // std::vector<std::string> techs_replicated = techs;
  int estimated_max_parts = techs.size();

  std::cout << "[INFO] Estimated maximum number of partitions: "
            << estimated_max_parts << std::endl;

  // the heavy lifting begins here
  std::cout << separator << std::endl;
  std::cout << std::setw((separator.size() + 10) / 2)
            << "TECH ASSIGNMENT BEGINS: " << std::endl;
  std::cout << separator << std::endl;
  GeneticPart(hypergraph_file, chiplet_io_file, chiplet_layer_file,
              chiplet_wafer_process_file, chiplet_assembly_process_file,
              chiplet_test_file, chiplet_netlist_file, chiplet_blocks_file,
              reach, separation, techs);

  auto end_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_time = end_time - start_time;
  std::cout << "[INFO] Total time taken for tech-aware ChipletPart: "
            << elapsed_time.count() << "s" << std::endl;
}

void ChipletPart::CreateMatingPool(
    const std::vector<std::vector<std::string>> &population,
    const std::vector<float> &fitness,
    std::vector<std::vector<std::string>> &mating_pool) {

  mating_pool.clear(); // Clear existing entries in mating pool
  std::vector<int> indices(population.size());
  for (size_t i = 0; i < indices.size(); ++i) {
    indices[i] = i;
  }

  // Sorting indices based on fitness values (descending order)
  std::sort(indices.begin(), indices.end(), [&](const int &a, const int &b) {
    return fitness[a] < fitness[b];
  });

  // Select the best individuals to fill the mating pool
  mating_pool.reserve(population_size_);

  int idx = 0;
  while (mating_pool.size() < population_size_) {
    mating_pool.push_back(population[indices[idx]]);
    idx = (idx + 1) % indices.size(); // Wrap around if needed
  }
}

void ChipletPart::GeneticPart(
    std::string hypergraph_file, std::string chiplet_io_file,
    std::string chiplet_layer_file, std::string chiplet_wafer_process_file,
    std::string chiplet_assembly_process_file, std::string chiplet_test_file,
    std::string chiplet_netlist_file, std::string chiplet_blocks_file,
    float reach, float separation, std::vector<std::string> &tech_nodes) {

  int no_improvement_counter = 0;
  int max_no_improvement_generations = 2;
  const std::string hseparator(60, '='); // Creates a separator line
  std::cout << hseparator << std::endl;
  std::cout << std::setw((hseparator.size() + 10) / 2)
            << "RUNNING GENETIC ALGORITHM! " << std::endl;
  std::cout << hseparator << std::endl;
  Py_Initialize();
  PyObject *sysPath = PySys_GetObject("path");
  PyList_Append(sysPath,
                PyUnicode_FromString(path_to_embedding_script_.c_str()));

  // the objective of the genetic algorithm is to find
  // a near optimal assignment of tech nodes to partitions
  // such that the cost of the assignment is minimized
  // the cost of running a partition with a given tech
  // node is obtained from the function QuickPart

  // the genetic algorithm is implemented as follows:
  // 1. generate the initial population
  // 2. evaluate the fitness of each individual
  // 3. select the best individuals
  // 4. crossover the best individuals
  // 5. mutate the best individuals
  // 6. repeat steps 2-5 until the stopping criteria is met

  // the stopping criteria is met when the number of generations
  // is reached or the cost of the best individual does not change
  // for a given number of generations or an early threshold is set

  auto CrossOver = [this](std::vector<std::string> &parent1,
                          std::vector<std::string> &parent2) {
    // Determine the minimum size of the two parents
    size_t min_size = std::min(parent1.size(), parent2.size());

    // Create the child with the size of the shorter parent
    std::vector<std::string> child(min_size);

    // Generate a crossover point based on the shorter parent
    std::uniform_int_distribution<int> uni_dist(0, min_size - 1);
    int c = uni_dist(rng_);

    // Perform the crossover, up to the min size
    for (size_t i = 0; i < min_size; ++i) {
      child[i] = (i < c) ? parent1[i] : parent2[i];
    }

    return child;
  };

  auto Mutate = [this](std::vector<std::string> &individual,
                       std::vector<std::string> &tech_nodes) {
    std::uniform_int_distribution<int> part_dist(0, num_parts_ - 1);
    std::uniform_int_distribution<int> tech_dist(0, tech_nodes.size() - 1);
    individual[part_dist(rng_)] = tech_nodes[tech_dist(rng_)];
  };

  // step 1 begins

  std::vector<std::vector<std::string>> population;
  std::vector<float> fitness;
  std::vector<std::string> best_individual;
  std::vector<int> best_partition;
  float best_cost = std::numeric_limits<float>::max();
  int num_generations = 20;
  int num_individuals = 10;
  // seed the random generator
  rng_.seed(std::random_device{}());
  num_parts_ = tech_nodes.size();

  // assign the three technodes to num_parts_
  for (int i = 0; i < num_individuals; i++) {
    std::vector<std::string> individual;
    // set num_part a random number between 1 and 4
    int num_parts = rand() % 8 + 1;
    for (int j = 0; j < num_parts; j++) {
      individual.push_back(tech_nodes[rand() % tech_nodes.size()]);
    }
    population.push_back(individual);
  }

  bool done = false;
  const std::string separator(60, '*'); // Creates a separator line
  for (int i = 0; i < num_generations_; ++i) {
    std::cout << separator << std::endl;
    std::cout << "[INFO] Starting [generation " << i << "] with [population "
              << population.size() << "]" << std::endl;
    std::cout << separator << std::endl;
    // print the technodes for this generation
    for (int j = 0; j < population.size(); ++j) {
      std::cout << "[INFO] Individual " << j << ": ";
      for (auto &tech : population[j]) {
        std::cout << tech << " ";
      }
      std::cout << std::endl;
    }
    if (population.size() < num_individuals) {
      std::cerr << "Error: Population size is less than the number of "
                   "individuals"
                << std::endl;
      return;
    }
    // 1. Evaluate fitness
    for (int j = 0; j < population.size(); ++j) {
      // print the individual
      auto part_tuple =
          QuickPart(hypergraph_file, chiplet_io_file, chiplet_layer_file,
                    chiplet_wafer_process_file, chiplet_assembly_process_file,
                    chiplet_test_file, chiplet_netlist_file,
                    chiplet_blocks_file, reach, separation, population[j]);
      float cost = std::get<0>(part_tuple);
      std::cout << "[INFO] Cost for individual " << j << " is " << cost
                << std::endl;
      std::vector<int> part = std::get<1>(part_tuple);
      int num_parts = std::set<int>(part.begin(), part.end()).size();
      fitness.push_back(cost);
      if (cost < best_cost) {
        best_cost = cost;
        best_partition = std::get<1>(part_tuple);
        best_individual = population[j];
      }
    }

    // print the best cost for this generation and its tech assignment
    std::cout << "[INFO] Best cost for generation " << i << " is " << best_cost
              << " for "
              << *std::max_element(best_partition.begin(),
                                   best_partition.end()) +
                     1
              << " partitions" << std::endl;

    std::cout << "[INFO] Best individual for generation " << i << ": ";
    for (int j = 0; j < best_individual.size(); j++) {
      std::cout << best_individual[j] << " ";
    }

    std::cout << std::endl;

    // 2. Selection
    std::vector<std::vector<std::string>> mating_pool;
    CreateMatingPool(population, fitness, mating_pool);

    // 3. Crossover
    std::vector<std::vector<std::string>> new_population;

    for (size_t j = 0; j < mating_pool.size() / 2; ++j) {
      // Debug prints for mating pool access
      if (2 * j < mating_pool.size() && 2 * j + 1 < mating_pool.size()) {
        std::vector<std::string> parent1 = mating_pool[2 * j];
        std::vector<std::string> parent2 = mating_pool[2 * j + 1];
        std::vector<std::string> child1 = CrossOver(parent1, parent2);
        std::vector<std::string> child2 = CrossOver(parent2, parent1);
        new_population.push_back(child1);
        new_population.push_back(child2);
      } else {
        std::cout << "Error: Out of bounds access to mating pool!" << std::endl;
        break;
      }
    }

    // 4. Mutation
    std::uniform_real_distribution<float> mutation_dist(0.0, 1.0);
    for (auto &individual : new_population) {
      if (mutation_dist(rng_) < mutation_rate_) {
        Mutate(individual, tech_nodes);
      }
    }

    // Optional: Elitism - Retain some of the best individuals from the old
    // population
    int elitism_count = 2; // Keep the top 2 individuals

    // sort the population based on fitness
    std::vector<int> indices(population.size());
    for (size_t j = 0; j < indices.size(); ++j) {
      indices[j] = j;
    }

    // Sorting indices based on fitness values (ascending order)

    std::sort(indices.begin(), indices.end(), [&](const int &a, const int &b) {
      return fitness[a] < fitness[b];
    });

    // now sort the population based on the indices
    std::vector<std::vector<std::string>> sorted_population;
    std::vector<float> sorted_fitness;
    for (size_t j = 0; j < indices.size(); ++j) {
      sorted_population.push_back(population[indices[j]]);
      sorted_fitness.push_back(fitness[indices[j]]);
    }

    population = sorted_population;

    // Keep the elite individuals
    std::vector<std::vector<std::string>> elite_individuals(
        population.begin(), population.begin() + elitism_count);

    // Replace the old population with the new population, but add back the
    // elite
    new_population.insert(new_population.end(), elite_individuals.begin(),
                          elite_individuals.end());
    population = new_population;

    // Sort population by fitness if not already sorted
    std::sort(fitness.begin(), fitness.end(), std::greater<float>());

    if (i > gen_threshold_) {
      if (best_cost <= fitness[0]) {
        no_improvement_counter++;
      } else {
        no_improvement_counter = 0; // Reset if there's improvement
        best_cost = fitness[0];     // Update best cost
      }

      // Stop if no improvement for a certain number of generations
      if (no_improvement_counter >= max_no_improvement_generations) {
        done = true;
        break;
      }
    }

    if (done) {
      break;
    }
  }

  // write the best partition to a file
  std::ofstream partition_output("tech_assignment.chipletpart.parts." +
                                 std::to_string(num_parts_));
  
  for (auto &part : best_partition) {
    partition_output << part << std::endl;
  }

  partition_output.close();

  std::cout << "[INFO] Best cost for the tech assignment is " << best_cost
            << std::endl;
  std::cout << "[INFO] Best individual for the tech assignment: ";
  for (int j = 0; j < best_individual.size(); j++) {
    std::cout << best_individual[j] << " ";
  }
  std::cout << std::endl;

  // write the best individual to a file
  std::ofstream tech_assignment_file("tech_assignment.chipletpart.techs." +
                                     std::to_string(num_parts_));
  for (auto &tech : best_individual) {
    tech_assignment_file << tech << std::endl;
  }

  tech_assignment_file.close();
}

std::tuple<float, std::vector<int>> ChipletPart::InitQuickPart(
    std::string hypergraph_file, std::string chiplet_io_file,
    std::string chiplet_layer_file, std::string chiplet_wafer_process_file,
    std::string chiplet_assembly_process_file, std::string chiplet_test_file,
    std::string chiplet_netlist_file, std::string chiplet_blocks_file,
    float reach, float separation, std::vector<std::string> &tech_nodes) {
  Matrix<float> upper_block_balance;
  Matrix<float> lower_block_balance;
  std::vector<float> zero_vector(vertex_dimensions_, 0.0);
  // removing balance constraints
  for (int i = 0; i < num_parts_; i++) {
    upper_block_balance.emplace_back(hypergraph_->GetTotalVertexWeights());
    lower_block_balance.emplace_back(zero_vector);
  }

  num_parts_ = tech_nodes.size();
  // Do num_parts-way cuts
  std::vector<int> kway_partition = METISPart(num_parts_);
  // remove all part.* files
  std::string command = "rm -f .part.*";
  int status = system(command.c_str());

  std::vector<int> reaches(hypergraph_->GetNumHyperedges(), reach);
  auto chiplet_refiner = std::make_shared<ChipletRefiner>(
      num_parts_, refine_iters_, max_moves_, reaches);
  chiplet_refiner->SetReach(reach);
  chiplet_refiner->SetSeparation(separation);
  // Initialize the cost model here
  chiplet_refiner->InitCostModel(
      chiplet_io_file, chiplet_layer_file, chiplet_wafer_process_file,
      chiplet_assembly_process_file, chiplet_test_file, chiplet_netlist_file,
      chiplet_blocks_file);
  chiplet_refiner->SetTechArray(tech_nodes);
  chiplet_refiner->SetNumParts(num_parts_);
  auto f_tuple = chiplet_refiner->RunFloorplanner(kway_partition, hypergraph_,
                                                  2000, 500, 1.0);
  std::vector<float> aspect_ratios = std::get<0>(f_tuple);
  std::vector<float> x_locations = std::get<1>(f_tuple);
  std::vector<float> y_locations = std::get<2>(f_tuple);
  chiplet_refiner->SetAspectRatios(aspect_ratios);
  chiplet_refiner->SetXLocations(x_locations);
  chiplet_refiner->SetYLocations(y_locations);
  chiplet_refiner->InitSlopes(num_parts_);
  float cost = chiplet_refiner->GetCostFromScratch(kway_partition);
  chiplet_refiner->SetLegacyCost(cost);
  chiplet_refiner->SetRefineIters(1);
  chiplet_refiner->SetMove(10);
  chiplet_refiner->Refine(hypergraph_, upper_block_balance, lower_block_balance,
                          kway_partition);
  int max_part =
      *std::max_element(kway_partition.begin(), kway_partition.end());
  int num_parts = std::set<int>(kway_partition.begin(), kway_partition.end())
                      .size(); // number of unique elements
  cost = chiplet_refiner->GetCostFromScratch(kway_partition);
  std::cout << "[GENETIC-RUNS] [Parts " << num_parts << "] [Chiplet cost "
            << cost << "] ";
  std::cout << "[FP feasible ";
  f_tuple = chiplet_refiner->RunFloorplanner(kway_partition, hypergraph_, 2000,
                                             500, 1.0);
  if (std::get<3>(f_tuple)) {
    std::cout << "YES]" << std::endl;
  } else {
    std::cout << "NO]" << std::endl;
  }

  return std::make_tuple(cost, kway_partition);
}

std::tuple<float, std::vector<int>> ChipletPart::QuickPart(
    std::string hypergraph_file, std::string chiplet_io_file,
    std::string chiplet_layer_file, std::string chiplet_wafer_process_file,
    std::string chiplet_assembly_process_file, std::string chiplet_test_file,
    std::string chiplet_netlist_file, std::string chiplet_blocks_file,
    float reach, float separation, std::vector<std::string> &tech_nodes) {
  rng_.seed(seed_);
  num_parts_ = tech_nodes.size();
  auto start_time_stamp_global = std::chrono::high_resolution_clock::now();

  // Pre-allocate vectors to avoid resizing
  std::vector<std::vector<int>> init_partitions;

  // Prepare upper and lower block balance
  std::vector<float> total_weights = hypergraph_->GetTotalVertexWeights();
  Matrix<float> upper_block_balance(num_parts_, total_weights);
  Matrix<float> lower_block_balance(
      num_parts_, std::vector<float>(vertex_dimensions_, 0.0));

  if (num_parts_ == 1) {
    std::vector<int> partition(hypergraph_->GetNumVertices(), 0);
    init_partitions.push_back(partition);
  } else {
    init_partitions.push_back(METISPart(num_parts_));
  }

  // Remove .part.* files safely
  if (std::filesystem::exists(".part.*")) {
    std::filesystem::remove_all(
        ".part.*"); // C++17 feature, faster than system call
  }

  // Step 4: Run FMRefiner
  auto chiplet_refiner = std::make_shared<ChipletRefiner>(
      num_parts_, refine_iters_, max_moves_,
      std::vector<int>(hypergraph_->GetNumHyperedges(), reach));
  chiplet_refiner->SetReach(reach);
  chiplet_refiner->SetSeparation(separation);
  chiplet_refiner->InitCostModel(
      chiplet_io_file, chiplet_layer_file, chiplet_wafer_process_file,
      chiplet_assembly_process_file, chiplet_test_file, chiplet_netlist_file,
      chiplet_blocks_file);

  // Refactoring loops and minimizing redundant code
  size_t best_partition_idx = 0;
  float best_cost = std::numeric_limits<float>::max();

  for (size_t idx = 0; idx < init_partitions.size(); ++idx) {
    std::vector<int> &partition = init_partitions[idx];
    int max_part = *std::max_element(partition.begin(), partition.end());
    int num_parts = std::set<int>(partition.begin(), partition.end()).size();
    std::vector<std::string> tech_array;
    for (int i = 0; i < num_parts; i++) {
      tech_array.push_back(tech_nodes[i]);
    }

    chiplet_refiner->SetTechArray(tech_array);
    chiplet_refiner->SetNumParts(max_part + 1);

    if (max_part == 0) {
      std::vector<float> aspect_ratios = {1.0};
      std::vector<float> x_locations = {0.0};
      std::vector<float> y_locations = {0.0};
      chiplet_refiner->SetAspectRatios(aspect_ratios);
      chiplet_refiner->SetXLocations(x_locations);
      chiplet_refiner->SetYLocations(y_locations);
      chiplet_refiner->InitSlopes(max_part + 1);
      float cost = chiplet_refiner->GetCostFromScratch(partition);
      int num_parts = std::set<int>(partition.begin(), partition.end()).size();
      // exit(EXIT_SUCCESS);
      best_cost = cost;
      best_partition_idx = idx;
    } else {
      auto f_tuple = chiplet_refiner->RunFloorplanner(partition, hypergraph_,
                                                      2000, 500, 1.0);
      std::vector<float> aspect_ratios = std::get<0>(f_tuple);
      std::vector<float> x_locations = std::get<1>(f_tuple);
      std::vector<float> y_locations = std::get<2>(f_tuple);
      chiplet_refiner->SetAspectRatios(aspect_ratios);
      chiplet_refiner->SetXLocations(x_locations);
      chiplet_refiner->SetYLocations(y_locations);
      chiplet_refiner->InitSlopes(max_part + 1);
      float cost = chiplet_refiner->GetCostFromScratch(partition);
      chiplet_refiner->SetLegacyCost(cost);
      // exit(EXIT_SUCCESS);
      // do the floorplan-aware refinement here
      chiplet_refiner->Refine(hypergraph_, upper_block_balance,
                              lower_block_balance, partition);
      max_part = *std::max_element(partition.begin(), partition.end());
      cost = chiplet_refiner->GetCostFromScratch(partition, true);
      int num_parts = std::set<int>(partition.begin(), partition.end()).size();
      f_tuple = chiplet_refiner->RunFloorplanner(partition, hypergraph_, 2000,
                                                 500, 1.0);
      if (std::get<3>(f_tuple)) {
        best_cost = cost;
        best_partition_idx = idx;
      } else {
        best_cost = std::numeric_limits<float>::max();
        best_partition_idx = idx;
      }
    }
  }

  return std::make_tuple(best_cost, init_partitions[best_partition_idx]);
}

void ChipletPart::Partition(
    std::string hypergraph_file, std::string chiplet_io_file,
    std::string chiplet_layer_file, std::string chiplet_wafer_process_file,
    std::string chiplet_assembly_process_file, std::string chiplet_test_file,
    std::string chiplet_netlist_file, std::string chiplet_blocks_file,
    float reach, float separation, std::string tech) {

  rng_.seed(seed_);

  // Log function (optional: consider reducing verbosity)
  auto log_info = [](const std::string &msg) {
    std::cout << "[INFO] " << msg << std::endl;
  };

  log_info("Reading the chiplet graph " + hypergraph_file);
  log_info("Reading the chiplet IO file " + chiplet_io_file);
  log_info("Reading the chiplet layer file " + chiplet_layer_file);
  log_info("Reading the chiplet wafer process file " +
           chiplet_wafer_process_file);
  log_info("Reading the chiplet assembly process file " +
           chiplet_assembly_process_file);
  log_info("Reading the chiplet test file " + chiplet_test_file);
  log_info("Reading the chiplet netlist file " + chiplet_netlist_file);
  log_info("Reading the chiplet blocks file " + chiplet_blocks_file);
  log_info("Reach: " + std::to_string(reach));
  log_info("Separation: " + std::to_string(separation));
  log_info("Tech: " + tech);

  auto start_time_stamp_global = std::chrono::high_resolution_clock::now();

  // Pre-allocate vectors to avoid resizing
  std::vector<std::vector<int>> init_partitions;
  init_partitions.reserve(10); // assuming we need to store around 10 partitions

  // Prepare upper and lower block balance
  std::vector<float> total_weights = hypergraph_->GetTotalVertexWeights();
  Matrix<float> upper_block_balance(num_parts_, total_weights);
  Matrix<float> lower_block_balance(
      num_parts_, std::vector<float>(vertex_dimensions_, 0.0));

  // Step 1: Run spectral partition
  std::vector<int> init_spec_partition = SpectralPartition(hypergraph_file);
  if (!init_spec_partition.empty() &&
      std::find(init_spec_partition.begin(), init_spec_partition.end(), -1) ==
          init_spec_partition.end()) {
    log_info("Spectral partition obtained!");
    init_partitions.push_back(
        std::move(init_spec_partition)); // Use move to avoid copying
  } else {
    std::cerr << "[ERROR] Spectral partition failed! Ignoring this!"
              << std::endl;
  }

  // Step 2: Run CrossBarExpansion
  float quantile = 0.99;
  std::vector<int> crossbars = FindCrossbars(quantile);
  log_info("High-degree nodes found: " + std::to_string(crossbars.size()));

  for (int i : chiplets_set_) {
    if (i == 1) {
      init_partitions.push_back(std::vector<int>(hypergraph_->GetNumVertices(),
                                                 0)); // create empty partition
    } else {
      std::vector<int> crossbar_partition =
          CrossBarExpansion(crossbars, num_parts_);
      if (!crossbar_partition.empty()) {
        init_partitions.push_back(std::move(crossbar_partition));
      }
    }
  }

  // Step 3: Run KWayCuts
  for (int i : chiplets_set_) {
    if (i != 1) {
      init_partitions.push_back(KWayCuts(i));
      init_partitions.push_back(METISPart(i));
    }
  }

  // Remove .part.* files safely
  if (std::filesystem::exists(".part.*")) {
    std::filesystem::remove_all(
        ".part.*"); // C++17 feature, faster than system call
  }

  log_info("Random and METIS partitions obtained!");
  log_info("Total number of initial partitions: " +
           std::to_string(init_partitions.size()));

  // Step 4: Run FMRefiner
  auto chiplet_refiner = std::make_shared<ChipletRefiner>(
      num_parts_, refine_iters_, max_moves_,
      std::vector<int>(hypergraph_->GetNumHyperedges(), reach));
  chiplet_refiner->SetReach(reach);
  chiplet_refiner->SetSeparation(separation);
  chiplet_refiner->InitCostModel(
      chiplet_io_file, chiplet_layer_file, chiplet_wafer_process_file,
      chiplet_assembly_process_file, chiplet_test_file, chiplet_netlist_file,
      chiplet_blocks_file);

  // Refactoring loops and minimizing redundant code
  size_t best_partition_idx = 0;
  float best_cost = std::numeric_limits<float>::max();

  for (size_t idx = 0; idx < init_partitions.size(); ++idx) {
    std::vector<int> &partition = init_partitions[idx];
    int max_part = *std::max_element(partition.begin(), partition.end());
    std::vector<std::string> tech_array(max_part + 1, tech);

    chiplet_refiner->SetTechArray(tech_array);
    chiplet_refiner->SetNumParts(max_part + 1);

    if (max_part == 0) {
      std::vector<float> aspect_ratios = {1.0};
      std::vector<float> x_locations = {0.0};
      std::vector<float> y_locations = {0.0};
      chiplet_refiner->SetAspectRatios(aspect_ratios);
      chiplet_refiner->SetXLocations(x_locations);
      chiplet_refiner->SetYLocations(y_locations);
      chiplet_refiner->InitSlopes(max_part + 1);
      float cost = chiplet_refiner->GetCostFromScratch(partition);
      int num_parts = std::set<int>(partition.begin(), partition.end()).size();
      std::cout << "[COST-FPL-INFO] [Partition " << idx << "] [Parts "
                << num_parts << "] [Chiplet cost " << cost << "] ";
      std::cout << "[FP feasible YES]" << std::endl;
      // exit(EXIT_SUCCESS);
      best_cost = cost;
      best_partition_idx = idx;
    } else {
      auto f_tuple = chiplet_refiner->RunFloorplanner(partition, hypergraph_,
                                                      2000, 500, 1.0);
      std::vector<float> aspect_ratios = std::get<0>(f_tuple);
      std::vector<float> x_locations = std::get<1>(f_tuple);
      std::vector<float> y_locations = std::get<2>(f_tuple);
      chiplet_refiner->SetAspectRatios(aspect_ratios);
      chiplet_refiner->SetXLocations(x_locations);
      chiplet_refiner->SetYLocations(y_locations);
      chiplet_refiner->InitSlopes(max_part + 1);
      float cost = chiplet_refiner->GetCostFromScratch(partition);
      chiplet_refiner->SetLegacyCost(cost);
      // exit(EXIT_SUCCESS);
      // do the floorplan-aware refinement here
      chiplet_refiner->Refine(hypergraph_, upper_block_balance,
                              lower_block_balance, partition);
      max_part = *std::max_element(partition.begin(), partition.end());
      cost = chiplet_refiner->GetCostFromScratch(partition, true);
      int num_parts = std::set<int>(partition.begin(), partition.end()).size();
      std::cout << "[COST-FPL-INFO] [Partition " << idx << "] [Parts "
                << num_parts << "] [Chiplet cost " << cost << "] ";
      std::cout << "[FP feasible ";
      f_tuple = chiplet_refiner->RunFloorplanner(partition, hypergraph_, 2000,
                                                 500, 1.0);
      if (std::get<3>(f_tuple)) {
        std::cout << "YES]" << std::endl;
        if (cost < best_cost) {
          best_cost = cost;
          best_partition_idx = idx;
        }
      } else {
        std::cout << "NO]" << std::endl;
      }
    }
  }

  log_info("Best partition cost: " + std::to_string(best_cost));
  log_info(
      "Number of chiplets in best partition: " +
      std::to_string(std::set<int>(init_partitions[best_partition_idx].begin(),
                                   init_partitions[best_partition_idx].end())
                         .size()));

  // Write best partition to file
  std::ofstream best_partition_file(hypergraph_file + ".chiplet.part");
  for (const auto &part : init_partitions[best_partition_idx]) {
    best_partition_file << part << std::endl;
  }

  tech_parts_ = *std::max_element(init_partitions[best_partition_idx].begin(),
                                  init_partitions[best_partition_idx].end()) +
                1;

  auto end_time_stamp_global = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_time_global =
      end_time_stamp_global - start_time_stamp_global;

  log_info("Total time taken: " + std::to_string(elapsed_time_global.count()) +
           " seconds");
  log_info("Partitioning completed!");

  float total_fplan_time = chiplet_refiner->GetTotalFloorplanTime();
  log_info("Total floorplanning time: " + std::to_string(total_fplan_time) +
           " seconds");

  float total_cost_model_time = chiplet_refiner->GetCostModelTime();
  log_info("Total cost model time: " + std::to_string(total_cost_model_time) +
           " seconds");
}

void ChipletPart::EvaluatePartition(
    std::string hypergraph_file, std::string hypergraph_part,
    std::string chiplet_io_file, std::string chiplet_layer_file,
    std::string chiplet_wafer_process_file,
    std::string chiplet_assembly_process_file, std::string chiplet_test_file,
    std::string chiplet_netlist_file, std::string chiplet_blocks_file,
    float reach, float separation, std::string tech) {
  Py_Initialize();
  PyObject *sysPath = PySys_GetObject("path");
  PyList_Append(sysPath,
                PyUnicode_FromString(path_to_embedding_script_.c_str()));

  std::cout << "[INFO] Reading the chiplet graph " << hypergraph_file
            << std::endl;
  std::cout << "[INFO] Reading the partition file " << hypergraph_part
            << std::endl;
  std::cout << "[INFO] Reading the chiplet IO file " << chiplet_io_file
            << std::endl;
  std::cout << "[INFO] Reading the chiplet layer file " << chiplet_layer_file
            << std::endl;
  std::cout << "[INFO] Reading the chiplet wafer process file "
            << chiplet_wafer_process_file << std::endl;
  std::cout << "[INFO] Reading the chiplet assembly process file "
            << chiplet_assembly_process_file << std::endl;
  std::cout << "[INFO] Reading the chiplet test file " << chiplet_test_file
            << std::endl;
  std::cout << "[INFO] Reading the chiplet netlist file "
            << chiplet_netlist_file << std::endl;
  std::cout << "[INFO] Reading the chiplet blocks file " << chiplet_blocks_file
            << std::endl;
  std::cout << "[INFO] Reach: " << reach << std::endl;
  std::cout << "[INFO] Separation: " << separation << std::endl;
  std::cout << "[INFO] Tech: " << tech << std::endl;

  auto start_time_stamp_global = std::chrono::high_resolution_clock::now();
  Matrix<float> upper_block_balance;
  Matrix<float> lower_block_balance;
  std::vector<float> zero_vector(vertex_dimensions_, 0.0);
  // removing balance constraints
  for (int i = 0; i < num_parts_; i++) {
    upper_block_balance.emplace_back(hypergraph_->GetTotalVertexWeights());
    lower_block_balance.emplace_back(zero_vector);
  }

  std::vector<int> partition;

  // read hypergraph partitioning file
  std::ifstream hypergraph_part_input(hypergraph_part);
  int num_parts = 0;
  if (!hypergraph_part_input.is_open()) {
    std::cerr << "Error: Cannot open hypergraph partition file "
              << hypergraph_part << std::endl;
    return;
  } else {
    partition = std::vector<int>(hypergraph_->GetNumVertices(), 0);
    for (int i = 0; i < hypergraph_->GetNumVertices(); i++) {
      hypergraph_part_input >> partition[i];
      num_parts = std::max(num_parts, partition[i]);
    }
  }

  num_parts_ = num_parts + 1;
  std::cout << "[INFO] Number of partitions: " << num_parts_ << std::endl;
  std::vector<int> reaches(hypergraph_->GetNumHyperedges(), reach);
  auto chiplet_refiner = std::make_shared<ChipletRefiner>(
      num_parts_, refine_iters_, max_moves_, reaches);
  chiplet_refiner->SetReach(reach);
  chiplet_refiner->SetSeparation(separation);
  // Initialize the cost model here
  chiplet_refiner->InitCostModel(
      chiplet_io_file, chiplet_layer_file, chiplet_wafer_process_file,
      chiplet_assembly_process_file, chiplet_test_file, chiplet_netlist_file,
      chiplet_blocks_file);
  std::vector<std::string> tech_array(num_parts_, tech);
  chiplet_refiner->SetTechArray(tech_array);
  chiplet_refiner->SetNumParts(num_parts_);

  auto f_tuple =
      chiplet_refiner->RunFloorplanner(partition, hypergraph_, 2000, 500, 1.0);
  std::vector<float> aspect_ratios = std::get<0>(f_tuple);
  std::vector<float> x_locations = std::get<1>(f_tuple);
  std::vector<float> y_locations = std::get<2>(f_tuple);
  chiplet_refiner->SetAspectRatios(aspect_ratios);
  chiplet_refiner->SetXLocations(x_locations);
  chiplet_refiner->SetYLocations(y_locations);
  // print the aspect ratios
  std::cout << "[INFO] Aspect ratios: ";
  for (auto &ar : aspect_ratios) {
    std::cout << ar << " ";
  }
  std::cout << std::endl;
  chiplet_refiner->InitSlopes(num_parts_);
  std::cout << "Floorplan is feasible: " << std::get<3>(f_tuple) << std::endl;
  // print the partition
  for (int i = 0; i < partition.size(); i++) {
    std::cout << partition[i] << " ";
  }
  std::cout << std::endl;
  float cost = chiplet_refiner->GetCostFromScratch(partition);
  std::cout << "[COST-FPL-INFO] [Partition has cost " << cost
            << "] [FP feasible " << std::get<3>(f_tuple) << "]" << std::endl;
}

std::vector<int> ChipletPart::FindCrossbars(float &quantile) {
  std::vector<int> degrees(hypergraph_->GetNumVertices(), 0);
  for (int i = 0; i < hypergraph_->GetNumVertices(); i++) {
    degrees[i] = hypergraph_->GetNeighbors(i).size();
  }

  std::vector<int> sorted_degrees = degrees;
  std::sort(sorted_degrees.begin(), sorted_degrees.end());

  int threshold_index = static_cast<int>(quantile * sorted_degrees.size());
  int degree_threshold = sorted_degrees[threshold_index];

  std::vector<int> crossbars;
  for (int i = 0; i < hypergraph_->GetNumVertices(); i++) {
    if (degrees[i] >= degree_threshold) {
      crossbars.push_back(i);
    }
  }

  return crossbars;
}

} // namespace chiplet
