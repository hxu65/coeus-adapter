#include <iostream>
#include <vector>
#include <chrono>
#include <fstream>
#include <string>
#include <sstream>
#include <mpi.h>
#include <adios2.h>
#include <random>
#include <unistd.h>

void mpi_sleep(int time, int rank, std::string reason){
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank == 0) std::cout << reason << std::endl;
  sleep(time);
  if(rank == 0) std::cout << "Done " << reason << std::endl;
  MPI_Barrier(MPI_COMM_WORLD);
}

std::vector<float> generateRandomVector(std::size_t size) {
    // Use the current time as seed for random number generation
    std::mt19937 gen(static_cast<unsigned long>(std::time(nullptr)));
    // Define range for float data type between 1 and 100
    std::uniform_real_distribution<float> dist(0.0f, 100.0f);

    std::vector<float> result(size);

    for(std::size_t i = 0; i < size; ++i) {
        result[i] = dist(gen);
    }

    return result;
}

int main(int argc, char *argv[]) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if(argc < 6) {
        if(rank == 0) {
            std::cout << "Usage: " << argv[0] << " <N> <B> " << std::endl;
        }
        MPI_Finalize();
        return 1;
    }

    const int N = std::stoi(argv[1]);
    const size_t B = std::stoul(argv[2]);
    std::string config_path = argv[3];
    std::string out_file = argv[4];
    int role = std::stoi(argv[5]);

    if(rank==0) {
        std::cout << "Running I/O comparison with " << N << " steps, "
                  << B << " bytes per step, and " << size << " processes."
                  << " with role as " << role << std::endl;
    }
    double localGetTime = 0.0;
    double localPutTime = 0.0;
    std::string engine_name;
    mpi_sleep(10, rank, "INIT");
    if(role == 0 || role == -1){
        adios2::ADIOS adios(config_path, MPI_COMM_WORLD);
        adios2::IO io = adios.DeclareIO("TestIO");

        std::vector<float> data(B);
        auto var_data = io.DefineVariable<float>("U", {size_t(size), B}, {size_t(rank), 0},
                                                 {1, B}, adios2::ConstantDims);
        auto var_data2 = io.DefineVariable<float>("V", {size_t(size), B}, {size_t(rank), 0},
                                                  {1, B}, adios2::ConstantDims);
        auto data_mag = io.DefineDerivedVariable("pdfU",
                                                 "x:U \n"
                                                 "magnitude(x)",
                                                 adios2::DerivedVarType::StoreData);
        auto data_mag2 = io.DefineDerivedVariable("pdfV",
                                                  "x:V \n"
                                                  "magnitude(x)",
                                                  adios2::DerivedVarType::StoreData);

        auto engine = io.Open(out_file, adios2::Mode::Write);
        engine_name = engine.Name();
        mpi_sleep(10, rank, "Write_INIT");
        for (int i = 0; i < N; ++i) {
            data = generateRandomVector(B);
            engine.BeginStep();
            mpi_sleep(5, rank, "BeginStep");
            auto startPut = std::chrono::high_resolution_clock::now();
            engine.Put<float>(var_data, data.data());
            engine.Put<float>(var_data2, data.data());
            mpi_sleep(5, rank, "Put");
            engine.EndStep();
            mpi_sleep(5 , rank, "EndStep");
            auto endPut = std::chrono::high_resolution_clock::now();

            localPutTime += std::chrono::duration<double>(endPut - startPut).count();
        }
        engine.Close();

        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == 0) {
            std::cout << "\tPut done, time: " << localPutTime << std::endl;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) std::cout << "Sleeping" << std::endl;
    sleep(10);
    if (rank == 0) std::cout << "Done Sleeping" << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);

    if(role == 1 || role == -1){
        adios2::ADIOS adios(config_path, MPI_COMM_WORLD);
        adios2::IO io = adios.DeclareIO("TestIO");
        auto readEngine = io.Open(out_file, adios2::Mode::Read);

        std::vector<float> data;
        std::vector<float> derivedData;

        MPI_Barrier(MPI_COMM_WORLD);

        if (rank == 0) std::cout << "BeginStep" << std::endl;
        int i = 0;
        mpi_sleep(5, rank, "Read_INIT");
        while (readEngine.BeginStep() == adios2::StepStatus::OK) {
            mpi_sleep(5, rank, "BeginStep");
//      std::string var_name = "data_" + std::to_string(i) + "_" + std::to_string(rank);
            adios2::Variable<float> readVariable = io.InquireVariable<float>("pdfU");
            adios2::Variable<float> derVariable = io.InquireVariable<float>("pdfV");

            auto startGet = std::chrono::high_resolution_clock::now();
            readEngine.Get<float>(readVariable, data);
            readEngine.Get<float>(derVariable, derivedData);
            readEngine.EndStep();
            mpi_sleep(5, rank, "Gets");
            auto endGet = std::chrono::high_resolution_clock::now();
            localGetTime += std::chrono::duration<double>(endGet - startGet).count();
            i++;
        }
        readEngine.Close();

        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == 0) {
            std::cout << "\tGet done, time: " << localGetTime << std::endl;
        }
    }

    double globalPutTime, globalGetTime;
    MPI_Reduce(&localPutTime, &globalPutTime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&localGetTime, &globalGetTime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if(rank == 0) {
        std::string header = "Size,B,N,GlobalPutTime,GlobalGetTimem rank0Put, rank0Get\n";
        bool needHeader = false;

        auto filename = "io_comp_results.csv";
        std::cout << "Writing results to " << filename << std::endl;
        // Check if the file is empty or doesn't exist
        std::ifstream checkFile(filename);
        if (!checkFile.good() || checkFile.peek() == std::ifstream::traits_type::eof()) {
            needHeader = true;
        }
        checkFile.close();

        // Open the file for appending
        std::ofstream outputFile(filename, std::ios_base::app);

        // Write the header if needed
        if (needHeader) {
            outputFile << header;
        }

        // Append the results
        outputFile << size << "," << B << "," << N << ","
                   << globalPutTime << "," << globalGetTime << ","
                   << localPutTime << "," << localGetTime << std::endl;
        outputFile.close();
    }

    MPI_Finalize();
    return 0;
}
