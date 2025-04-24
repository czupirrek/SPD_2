#include "Scheduler.h"
#include <iostream>

int main(int argc, char* argv[]){
    // if(argc != 3){
    //     std::cerr << "[err] niepoprawna liczba argumentow" << std::endl;
    //     return 1;
    // }
    int n_machines = std::stoi(argv[1]);
    std::string filename = argv[2];
    std::string algo = argv[3];

    double ptas_K = 0;
    if (argc>4){
        ptas_K = std::stod(argv[4]);
    }

    std::cout << "[main] k: " << ptas_K << std::endl;

    Scheduler scheduler(n_machines, ptas_K);
    scheduler.loadTasks(filename);
    int cmax = scheduler.selectAlgorithm(algo);
    // scheduler.display();
    // int cmax = scheduler.startLSA();
    // int cmax = scheduler.startLPT();

    printf("[result] Cmax: %d\n", cmax);
    // scheduler.display();
    
    return 0;
}