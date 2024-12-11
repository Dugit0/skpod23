#include "mpi.h"
#include <iostream>
#include <string>
#include <algorithm>
#include <random>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>
#include <assert.h>
#include <math.h>
// #define Max(a, b) (((a) > (b)) ? (a) : (b))
// #define Min(a, b) (((a) < (b)) ? (a) : (b))
#define DEBUG 0

unsigned SECONDS_IN_DAY = 24*60*60;
unsigned SECONDS_IN_HOUR = 60*60;
unsigned SECONDS_IN_MINUTE = 60;
unsigned READ_QUORUM = 3;
unsigned WRITE_QUORUM = 3;
unsigned LOOP_LIMIT = 10;

enum MsgCode {
    REQ_READ,
    RESP_READ,
    REJ_READ,
    REQ_WRITE,
    RESP_WRITE,
    REJ_WRITE,
    SUCC_WRITE,
};


std::vector<unsigned> get_random_threads(unsigned n, unsigned num_threads, unsigned my_num) {
    assert (n < num_threads - 1);
    std::vector<unsigned> container{};
    for (unsigned i = 0; i < num_threads; i++) {
        if (i != my_num) {
            container.push_back(i);
        }
    }
    std::mt19937 g(rand());
    std::shuffle(container.begin(), container.end(), g);
    std::vector<unsigned> res{container.begin(), container.begin() + n};
    return res;
}


void get_timestamp(FILE* stream, int rank, std::string info="") {
    unsigned long cur_time = (unsigned long)time(NULL) % SECONDS_IN_DAY;
    unsigned hours = cur_time / SECONDS_IN_HOUR;
    unsigned minutes = cur_time % SECONDS_IN_HOUR / SECONDS_IN_MINUTE;
    unsigned seconds = cur_time % SECONDS_IN_MINUTE;

    fprintf(stream, "%02u:%02u:%02u||proc: %d||%s\n",
            hours, minutes, seconds, rank, info.c_str());
    // fflush(stream);
    return;
}



// void read() {
//     return;
// }


int main(void) {


    // Init MPI
    int status = MPI_Init(nullptr, nullptr);
    if (status) {
        printf("MPI not supported\nError with code %d\n", status);
        MPI_Abort(MPI_COMM_WORLD, status);
        return status;
    };

    // Get rank and num_threads
    int rank, num_threads;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_threads);

    // Quorum asserts
    assert((int) WRITE_QUORUM > num_threads / 2);
    assert((int) (READ_QUORUM + WRITE_QUORUM) > num_threads);

    // Random lock
    srand(42 ^ rank);

    // Open file
    // size_t buf_size = 50;
    // char filename[buf_size];
    // snprintf(filename, buf_size, "files/file_%d.txt", rank);

    // FILE* out_file = fopen(filename, "w");
    // if (!out_file) {
    //     printf("File %s can not open\n", filename);
    //     MPI_Abort(MPI_COMM_WORLD, 1);
    //     return 1;
    // };

    unsigned file_version = 0;
    int vote = 1;
    std::vector<unsigned> send_message_ids(num_threads, 0);
    std::vector<unsigned> recv_message_ids(num_threads, 0);

    // Sync barrier
    if (DEBUG >= 5)
        printf("Created %d/%d\n", rank, num_threads);
    MPI_Barrier(MPI_COMM_WORLD);

    // Start main loop
    // unsigned i = 0;
    // while (i < LOOP_LIMIT) {
    {
        // unsigned random_action = rand() % 10;
        unsigned random_action = rank;
        if (random_action == 0) {
            get_timestamp(stdout, rank, "READ");
            int buf = MsgCode::REQ_READ;
            // std::vector<MPI_Request> requests(num_threads, {});
            for (int proc = 0; proc < num_threads; proc++) {
                MPI_Send(&buf, 1, MPI_INT, proc, send_message_ids[proc], MPI_COMM_WORLD);
            }
        // } else if (random_action == 1) {
        //     get_timestamp(stdout, rank, "WRITE");
        } else {
            get_timestamp(stdout, rank, "IDLE");
            sleep(1);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        // i++;
    }

    // printf("Proc: %d, Rand: %d\n", rank, rand() % 100);
    // fflush(stdout);

    // sleep(rand() % 5);
    get_timestamp(stdout, rank);


    // Close file
    // fclose(out_file);

    // Finalize MPI
    MPI_Finalize();
    return 0;
}
