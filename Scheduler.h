#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <climits>
#include <cmath>
#include <numeric> // For std::accumulate
#include <set>     // For DP state


struct Task{
    int id;
    int cost;

    bool operator<(const Task& other) const {
        // Compare based on task number
        return id < other.id;
    }
};

struct Machine{
    int id;
    bool busy;
    int time;
};



class Scheduler{
    public:
        Scheduler(int n_machines, double ptas_K);
        void loadTasks(std::string fname);
        int startLSA();
        int startLPT();
        void display();
        int selectAlgorithm(std::string alg);
        int startPD();
        int startPermute();
        int PTAS(double k);
        int startFPTAS(double K);
        void prune_dp_states(std::set<std::vector<long long>>& dp);
        int startOWN();
    private:
        int tasksCount;
        std::vector<Task> tasks;
        std::vector<Machine> machines;
        int getMinMachineTime();
        void sortTasksDescending();
        double ptas_K = 0;
        void resetMachines();
    };


Scheduler::Scheduler(int n_machines, double ptas_K){
    this->ptas_K = ptas_K;
    for (int i = 0; i < n_machines; i++){
        Machine m;
        m.id = i;
        m.busy = false;
        m.time = 0;
        machines.push_back(m);
    }
}

void Scheduler::loadTasks(std::string filename){
    std::ifstream file;
    file.open(filename);
    if(!file.is_open()){
        std::cerr << "[err] nie znaleziono pliku" << std::endl;
        return;
    }
    file >> tasksCount;
    for(int i = 0; i < tasksCount; i++){
        Task task;
        file >> task.cost;
        task.id = i;
        tasks.push_back(task);
    }
    file.close();
    std::cerr << "[info] wczytano dane z pliku: " << filename << std::endl;
}


void Scheduler::display(){
    for (const auto& task : tasks){
        std::cout << "[" << task.id << "]\t cost: " << task.cost << std::endl;
    }
}

int Scheduler::startLSA(){
    // std::cout << "[info] start LSA" << std::endl;
    int cmax = 0;
    int startTime = 0;
    while (!tasks.empty()){
        for (int i = 0; i < machines.size(); i++){
            Machine &m = machines[i];
            startTime = getMinMachineTime();
            if (!tasks.empty() && m.time <= startTime){
                Task task = tasks.front();
                tasks.erase(tasks.begin());
                // m.busy = true;
                m.time = startTime + task.cost;
                // currentTime = std::min(currentTime, m.time);
                // startTime = getMinMachineTime();
                
                cmax = std::max(cmax, m.time);
                // std::cout << "[info] task: " << task.id << " assigned to machine: " << m.id << std::endl;
                // std::cout << "[info] current time: " << currentTime << std::endl;
                // std::cout << "[info] machine: " << m.id << " time: " << m.time << std::endl;
                // std::cout << "[info] cmax: " << cmax << std::endl;
                // printf("[info] task: %d assigned to machine: %d. current time is %d. cost is %d. will end on %d. \n", task.id, m.id, startTime, task.cost, m.time);
            } else if (m.time >= startTime){
                m.busy = false;
            }
        }
    }
    return cmax;

}

int Scheduler::getMinMachineTime(){
    int minTime = machines[0].time;
    for (int i = 0; i < machines.size(); i++){
        if (machines[i].time < minTime){
            minTime = machines[i].time;
        }
    }
    return minTime;
}


void Scheduler::sortTasksDescending(){
    std::sort(tasks.begin(), tasks.end(), [](const Task &a, const Task &b) {
        return a.cost > b.cost;
    });
}

int Scheduler::startLPT(){
    std::cout << "[info] start LPT" << std::endl;
    sortTasksDescending();
    int cmax = 0;
    cmax = startLSA();
    return cmax;
}


int Scheduler::startPD(){
    std::cout << "[info] start PD dla 2 maszyn" << std::endl;
    if (machines.size() != 2){
        std::cerr << "[err] DP działa tylko dla 2 maszyn" << std::endl;
        return -1;
    }

    int sum = 0;
    for (const auto& t : tasks) sum += t.cost;

    std::vector<bool> dp(sum + 1, false);
    dp[0] = true;

    for (const auto& task : tasks){
        for (int j = sum; j >= task.cost; j--){
            dp[j] = dp[j] || dp[j - task.cost];
        }
    }

    int best = sum;
    for (int i = 0; i <= sum; i++){
        if (dp[i]){
            int other = sum - i;
            int cmax = std::max(i, other);
            best = std::min(best, cmax);
        }
    }

    return best;

}

void Scheduler::resetMachines(){
    for (auto& m : machines){
        m.time = 0;
        m.busy = false;
    }
}

long factorial(const int n)
{
    long f = 1;
    for (int i=1; i<=n; ++i)
        f *= i;
    return f;
}

int Scheduler::startPermute() {
    std::cout << "[info] start permutacji" << std::endl;
    int best_cmax = INT_MAX;
    std::vector<Task> original_tasks = tasks;  // Zachowaj oryginalną listę zadań
    std::sort(tasks.begin(), tasks.end());
    long all_permutations = factorial(tasks.size());
    int permCount = 0;
    do {
        // Przywróć stan początkowy maszyn
        for (auto& m : machines) {
            m.time = 0;
            m.busy = false;
        }
        
        // Wykonaj kopię zadań dla tej permutacji
        std::vector<Task> current_tasks = tasks;
        tasks = current_tasks;  // Upewnij się, że startLSA działa na właściwej liście
        
        int current_cmax = startLSA();
        if (current_cmax < best_cmax) {
            best_cmax = current_cmax;
        }
        
        // Przywróć oryginalną listę zadań dla next_permutation
        tasks = original_tasks;
        std::sort(tasks.begin(), tasks.end());  // Potrzebne dla next_permutation
        // Przesuń do aktualnej permutacji (kosztowne, ale konieczne)
        for (int i = 0; i < original_tasks.size(); ++i) {
            tasks[i] = current_tasks[i];
        }
        if (permCount % 1000 == 0)
            fprintf(stderr, "[%d] progress: %.2f%%\n", permCount, (float)permCount/all_permutations*100);
        permCount++;
    } while (std::next_permutation(tasks.begin(), tasks.end()));
    
    tasks = original_tasks;  // Przywróć oryginalne zadania
    return best_cmax;
}



int Scheduler::selectAlgorithm(std::string s){
    int cmax = 0;
    if (s == "LSA"){
        std::cout << "[info] wybrano algorytm LSA" << std::endl;
        cmax = startLSA();
    } else if (s == "LPT"){
        std::cout << "[info] wybrano algorytm LPT" << std::endl;
        cmax = startLPT();
    } else if (s == "PD"){
        std::cout << "[info] wybrano algorytm PD" << std::endl;
        cmax = startPD();
    } else if (s == "PERMUTE"){
        std::cout << "[info] wybrano przegląd zupełny" << std::endl;
        cmax = startPermute();
    } else if (s == "PTAS"){
        std::cout << "[info] wybrano algorytm PTAS" << std::endl;
        cmax = PTAS(ptas_K);
    } 
    else if (s == "FPTAS"){
        std::cout << "[info] wybrano algorytm FPTAS" << std::endl;
        cmax = startFPTAS(ptas_K);
    } 
    else if (s == "OWN"){
        std::cout << "[info] wybrano algorytm wlasny" << std::endl;
        cmax = startOWN();
    } 

    else {
        std::cerr << "[err] niepoprawny algorytm" << std::endl;
        return -1;
    }
    return cmax;
}

int Scheduler::PTAS(double k) {

    // const double epsilon = k;  // Można dostosować
    // int k_actual = std::min(static_cast<int>(tasks.size()), 
    //                      static_cast<int>(ceil(1.0 / epsilon)));
    int k_actual = static_cast<int>(ceil(k));
    
        // 0. Sprawdź, czy k jest większe od liczby zadań
        if (k_actual > tasks.size()) {
            std::cerr << "[err] k jest większe od liczby zadań" << std::endl;
            return -1;
        }
    
        // 0. Sprawdź, czy są zadania
    if (tasks.empty()) return 0;
    
    
    
        // 1. Posortuj zadania malejąco według czasu wykonania
        sortTasksDescending();
    
        // 2. Podziel zadania na k najdłuższych i resztę
        k_actual = std::min(k_actual, static_cast<int>(tasks.size()));

        std::vector<Task> long_tasks(tasks.begin(), tasks.begin() + k_actual);
        std::vector<Task> short_tasks(tasks.begin() + k_actual, tasks.end());
        std::vector<int> assignment(k_actual, 0); // <- tylko tyle ile trzeba
        
        std::cout << "[info] start PTAS" << std::endl;
        std::cout << "[info] k (eps): " << k << std::endl;
        std::cout << "[info] actual_k: " << k_actual << std::endl;
        std::cout << "[PTAS] Long tasks: " << long_tasks.size() 
                  << ", short tasks: " << short_tasks.size() << std::endl;
    
        // 3. Przegląd wszystkich możliwych przydziałów k najdłuższych zadań
        int best_cmax = INT_MAX;
        
        while (true) {
            // Przygotuj maszyny
            resetMachines();
            
            // Przydziel długie zadania
            for (int i = 0; i < long_tasks.size(); ++i) {
                int machine_idx = assignment[i];
                machines[machine_idx].time += long_tasks[i].cost;
            }
            
            // Przydziel krótkie zadania algorytmem LPT
            for (const auto& task : short_tasks) {
                // Znajdź maszynę z minimalnym obciążeniem
                int min_machine = 0;
                int min_time = machines[0].time;
                for (int i = 1; i < machines.size(); ++i) {
                    if (machines[i].time < min_time) {
                        min_time = machines[i].time;
                        min_machine = i;
                    }
                }
                machines[min_machine].time += task.cost;
            }
            
            // Oblicz Cmax dla tego przydziału
            int current_cmax = 0;
            for (const auto& m : machines) {
                if (m.time > current_cmax) {
                    current_cmax = m.time;
                }
            }
            
            // Aktualizuj najlepsze rozwiązanie
            if (current_cmax < best_cmax) {
                best_cmax = current_cmax;
                std::cout << "[PTAS] New best Cmax: " << best_cmax << std::endl;
            }
            
            // Generuj następny przydział (inkrementacja z przeniesieniem)
            int pos = 0;
            while (pos < k_actual) {
                assignment[pos]++;
                if (assignment[pos] < machines.size()) {
                    break;
                }
                assignment[pos] = 0;
                pos++;
            }
            
            // Jeśli doszliśmy do końca wszystkich przydziałów
            if (pos == k_actual) {
                break;
            }
        }
        
        std::cout << "[PTAS] Final best Cmax: " << best_cmax << std::endl;
        return best_cmax;
}


int Scheduler::startOWN() {
    std::cout << "[info] start algorytmu własnego" << std::endl;
    int total_cost = 0;
    int perfect_cost = 0;
    int cmax = 0;
    std::vector<int> machines_loads(machines.size(), 0);


    for (const auto& task : tasks) {
        total_cost += task.cost;
    }
    perfect_cost = std::ceil(1.0* total_cost / machines.size());

    sortTasksDescending();

    std::vector<bool> task_assigned(tasks.size(), false); // Track assigned tasks
    int tasks_remaining = tasks.size();

    for (int i = 0; i < machines_loads.size() && tasks_remaining > 0; ++i) {
        // Try to fill machine i close to perfect_cost
        bool machine_can_take_more = true;
        while (machine_can_take_more && tasks_remaining > 0) {
            machine_can_take_more = false; // Assume we can't add more in this pass
            int best_task_idx = -1;
            // Find the largest unassigned task that fits
            for (int k = 0; k < tasks.size(); ++k) {
                 // Check if task k exists, is unassigned, and fits
                if (!task_assigned[k] && (tasks[k].cost + machines_loads[i] <= perfect_cost)) {
                    best_task_idx = k; // Found the largest fitting task (since tasks are sorted)
                    break; // Stop searching once the largest fitting task is found
                }
            }

            if (best_task_idx != -1) {
                // Assign the best fitting task found
                machines_loads[i] += tasks[best_task_idx].cost;
                task_assigned[best_task_idx] = true;
                tasks_remaining--;
                machine_can_take_more = true; // We added a task, maybe we can add more
            }
        } // End while trying to fill machine i
    } // End for each machine

    // --- Handle tasks that didn't fit in the first pass ---
    // Use a standard LPT approach for remaining tasks
    std::vector<int> remaining_task_indices;
    for(int k=0; k<tasks.size(); ++k) {
        if (!task_assigned[k]) {
            remaining_task_indices.push_back(k);
        }
    }

    // Assign remaining tasks (largest first) to the least loaded machine
    for (int task_idx : remaining_task_indices) { // Already sorted descending
        int least_loaded_machine = 0;
        for (int m = 1; m < machines_loads.size(); ++m) {
            if (machines_loads[m] < machines_loads[least_loaded_machine]) {
                least_loaded_machine = m;
            }
        }
        machines_loads[least_loaded_machine] += tasks[task_idx].cost;
        task_assigned[task_idx] = true; // Mark as assigned
        tasks_remaining--; // Should reach 0
    }

    // ... (calculate and return cmax - same as before) ...
     int cmax = 0;
     for (int load : machines_loads) {
         if (load > cmax) {
             cmax = load;
         }
     }
     // ... (print machine loads, perfect cost) ...
     return cmax;

}




bool is_component_wise_le(const std::vector<long long>& v1, const std::vector<long long>& v2, int m_size) {
    if (v1.size() != m_size || v2.size() != m_size) return false; // Should not happen
    for (int i = 0; i < m_size; ++i) {
        if (v1[i] > v2[i]) {
            return false;
        }
    }
    return true;
}


// FPTAS (Fully Polynomial-Time Approximation Scheme) for P || Cmax
// This implementation uses DP on scaled task costs.
// The state represents achievable load vectors on machines.
int Scheduler::startFPTAS(double K) {
    // std::cout << "[info] start FPTAS with epsilon = " << epsilon_param << std::endl;

    // if (epsilon_param <= 0) {
    //     std::cerr << "[err] Epsilon must be > 0 for FPTAS." << std::endl;
    //     return -1;
    // }
    if (tasks.empty()) return 0;
    if (machines.empty()) {
        std::cerr << "[err] Brak maszyn do harmonogramowania." << std::endl;
        return -1; // Or some error indicator
    }

    // 1. Calculate sum of task costs (for scaling factor)
    long long sum_costs = 0;
    for (const auto& task : tasks) {
        sum_costs += task.cost;
    }
    if (sum_costs == 0) return 0;

    // 2. Determine scaling factor K
    // K = epsilon * C_est / N, where C_est is an estimate of C*, N is number of tasks.
    // A safe estimate for C* is sum_costs. So K = epsilon * sum_costs / tasks.size()
    // Or K = epsilon * P_sum / (m*n) to bound scaled values or sums better.
    // Let's use K = epsilon * sum_costs / tasks.size() based on some standard presentations.
    // Need to handle tasks.size() == 0, but handled above.
    // K should be a double for the division.
    

    // if (K < 1e-9) { // Avoid division by near zero or very small K which leads to huge scaled costs
    //     // If K is very small, scaled costs are very large. This happens if epsilon*sum_costs is very small.
    //     // This implies tasks are small or epsilon is very small. If epsilon is small, this is expected.
    //     // If tasks are small, running LPT might be fast and give a good enough solution.
    //     // However, to stick to the FPTAS structure, we proceed, but need large types.
    //     // If K is exactly zero (only if epsilon or sum_costs is zero), we handled sum_costs == 0.
    //     // Epsilon is checked > 0. So K > 0 unless sum_costs = 0.
    //      if (K == 0) { // This case should not be reachable if sum_costs > 0 and epsilon > 0
    //           std::cerr << "[err] Scaling factor K calculated as zero unexpectedly." << std::endl;
    //           // Fallback or error
    //           sortTasksDescending();
    //           return startLSA(); // LPT fallback
    //      }
    //     // K is small but positive. Proceed with long long.
    // }


    // 3. Scale task costs
    std::vector<long long> scaled_costs;
    scaled_costs.reserve(tasks.size());
    for (const auto& task : tasks) {
        // scaled_cost = floor(task.cost / K)
        scaled_costs.push_back(static_cast<long long>(std::floor(task.cost / K)));
    }

    // 4. Dynamic Programming
    // dp state: a set of sorted vectors representing achievable scaled loads on machines.
    // std::set guarantees unique vectors and keeps them sorted lexicographically.
    std::set<std::vector<long long>> dp;

    // Initial state: all machines have 0 load
    dp.insert(std::vector<long long>(machines.size(), 0));

    std::cout << "[FPTAS] Starting DP with " << scaled_costs.size() << " scaled tasks." << std::endl;

    // Iterate through scaled tasks
    for (size_t i = 0; i < scaled_costs.size(); ++i) {
        long long p_prime = scaled_costs[i];
        std::set<std::vector<long long>> next_dp;

        // Generate new states by adding p_prime to each machine's load for each existing state
        for (const auto& current_loads : dp) {
            for (size_t m_idx = 0; m_idx < machines.size(); ++m_idx) {
                std::vector<long long> new_loads = current_loads;
                new_loads[m_idx] += p_prime;
                std::sort(new_loads.begin(), new_loads.end()); // Keep vectors sorted
                next_dp.insert(new_loads); // std::set handles uniqueness
            }
        }

        dp = next_dp; // Move to the new set of states

        // Prune dominated states to keep DP set size manageable
        // print size before prune
        // std::cout << "[FPTAS] DP size before prune (task " << i << "): " << dp.size() << std::endl;
        prune_dp_states(dp);
        // print size after prune
        // std::cout << "[FPTAS] DP size after prune (task " << i << "):  " << dp.size() << std::endl;

         if ((i + 1) % 100 == 0) {
             fprintf(stderr, "[FPTAS] Processed %zu/%zu tasks. DP states: %zu\n", i + 1, scaled_costs.size(), dp.size());
         }
    }

    std::cout << "[FPTAS] DP finished. Final DP states: " << dp.size() << std::endl;

    // 5. Find minimum makespan among final states
    long long min_scaled_makespan = -1; // Use -1 as initial sentinel for long long max

    for (const auto& final_loads : dp) {
        // Makespan for a load vector is the maximum load (last element since sorted)
        long long current_scaled_makespan = final_loads.back();
        if (min_scaled_makespan == -1 || current_scaled_makespan < min_scaled_makespan) {
            min_scaled_makespan = current_scaled_makespan;
        }
    }

    // 6. Convert back to original scale
    double final_cmax_double = static_cast<double>(min_scaled_makespan) * K;

    // The actual makespan must be an integer (sum of integer costs).
    // Rounding the double result is necessary.
    int final_cmax = static_cast<int>(std::round(final_cmax_double));

    std::cout << "[FPTAS] Estimated Cmax = " << final_cmax << std::endl;

    // The FPTAS guarantee is C <= (1+epsilon)C*.
    // This DP finds the optimal makespan for the *scaled* problem.
    // The makespan for the original problem using this assignment is then calculated.
    // The error comes from the floor operation during scaling.

    // Note: This specific DP implementation might be slow for large 'm' because
    // the number of vectors in DP can still grow significantly.
    // A more advanced FPTAS DP might use counts of machines at load levels or other techniques.
    // However, this structure is common for demonstrating the DP part of an FPTAS.

    return final_cmax;
}


// Helper for FPTAS pruning - see definition above class
// bool is_component_wise_le(const std::vector<long long>& v1, const std::vector<long long>& v2, int m_size) { ... }
void Scheduler::prune_dp_states(std::set<std::vector<long long>>& dp) {
     if (dp.empty() || dp.size() <= 1) return;

    // Copy states to a vector for easier iteration and pruning logic
    std::vector<std::vector<long long>> states_vec;
    states_vec.reserve(dp.size());
    for(const auto& vec : dp) {
        states_vec.push_back(vec);
    }

    // The states are already sorted lexicographically because they came from a std::set
    // Sorting by sum or makespan might improve pruning efficiency in some cases,
    // but lexicographical is required for the simple prune_set check below.
    // std::sort(states_vec.begin(), states_vec.end()); // Already sorted by set

    std::set<std::vector<long long>> pruned_set;

    for (const auto& current_vec : states_vec) {
        bool is_dominated = false;
        // Check if current_vec is dominated by any vector ALREADY in the pruned_set
        // Because states_vec is processed in lexicographical order, if current_vec
        // is not dominated by any vector currently in pruned_set, it is not dominated
        // by any vector that WILL be added later either (due to non-domination property
        // and sort order). Also, current_vec cannot dominate any vector already in
        // pruned_set, because current_vec comes later in the lexicographical order.
        // This specific pruning logic relies on the lexicographical order of states.
        for (const auto& kept_vec : pruned_set) {
             if (is_component_wise_le(kept_vec, current_vec, machines.size())) {
                 is_dominated = true;
                 break;
             }
        }

        if (!is_dominated) {
            pruned_set.insert(current_vec);
        }
    }

    dp = pruned_set; // Replace the original set with the pruned set
}

