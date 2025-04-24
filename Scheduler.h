#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>

struct Task{
    int id;
    int cost;
};

struct Machine{
    int id;
    bool busy;
    int time;
};



class Scheduler{
    public:
        Scheduler(int n_machines);
        void loadTasks(std::string fname);
        int startLSA();
        int startLPT();
        void display();
        int selectAlgorithm(std::string alg);

    private:
        int tasksCount;
        std::vector<Task> tasks;
        std::vector<Machine> machines;
        int getMinMachineTime();
        void sortTasksDescending();
    };


Scheduler::Scheduler(int n_machines){
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
    std::cout << "[info] start LSA" << std::endl;
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



int Scheduler::selectAlgorithm(std::string s){
    int cmax = 0;
    if (s == "LSA"){
        std::cout << "[info] wybrano algorytm LSA" << std::endl;
        cmax = startLSA();
    } else if (s == "LPT"){
        std::cout << "[info] wybrano algorytm LPT" << std::endl;
        cmax = startLPT();
    } else {
        std::cerr << "[err] niepoprawny algorytm" << std::endl;
        return -1;
    }
    return cmax;
}