#include "mpi.h"
#include <iostream>
#include <string>
#include <algorithm>
#include <random>
#include <vector>
#include <map>
#include <array>

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>
#include <assert.h>
#include <math.h>
// #define Max(a, b) (((a) > (b)) ? (a) : (b))
// #define Min(a, b) (((a) < (b)) ? (a) : (b))
#define DEBUG 5

unsigned SECONDS_IN_DAY = 24*60*60;
unsigned SECONDS_IN_HOUR = 60*60;
unsigned SECONDS_IN_MINUTE = 60;
unsigned READ_QUORUM = 4;
unsigned WRITE_QUORUM = 3;
unsigned LOOP_LIMIT = 10;

enum MsgType {
    REQ,
    RESP,
    RETURN_BACK,
};

enum ReqType {
    READ,
    WRITE,
};

enum RespType {
    REJECT,
    ACCEPT,
};

std::vector<unsigned> get_random_threads(unsigned n, unsigned num_threads, unsigned my_num) {
    assert (n <= num_threads - 1);
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


std::vector<std::array<int, 2>> send_recv_quorum(
        std::string quorum_type,
        int rank,
        int num_threads,
        std::vector<int> REQ_recv_buf,
        std::vector<MPI_Request> REQ_requests,
        std::vector<MPI_Status> REQ_statuses
        ) {

    int buf;
    unsigned quorum;
    if (quorum_type == "READ") {
        buf = ReqType::READ;
        quorum = READ_QUORUM;
    } else {
        buf = ReqType::WRITE;
        quorum = WRITE_QUORUM;
    }

    auto threads = get_random_threads(quorum, num_threads, rank);

    // Send requests to READ/WRITE_QUORUM processes
    for (auto proc: threads) {
        MPI_Send(&buf, 1, MPI_INT,
                 proc, MsgType::REQ, MPI_COMM_WORLD);
    }

    // Open RESP MPI_Irecv
    std::vector<std::array<int, 2>> RESP_recv_buf(threads.size(), {0, 0});
    std::vector<MPI_Request> RESP_requests(threads.size(), {});
    std::vector<MPI_Status> RESP_statuses{};
    for (unsigned i = 0; i < threads.size(); i++) {
        RESP_statuses.push_back({});
    }
    for (unsigned i = 0; i < threads.size(); i++) {
        MPI_Irecv(&(RESP_recv_buf[i].front()), 2, MPI_INT,
                  threads[i], MsgType::RESP, MPI_COMM_WORLD,
                  &RESP_requests[i]);
    }

    // Wait until all responses are received
    // And reject any requests (deadlock protection)
    while (true) {
        int all_resp_received = 0;
        MPI_Testall(RESP_requests.size(), &RESP_requests.front(),
                    &all_resp_received, &RESP_statuses.front());
        if (all_resp_received) {
            break;
        }
        while (true) {
            int any_req_received = 0;
            int index;
            MPI_Testany(REQ_requests.size(), &REQ_requests.front(), &index,
                        &any_req_received, &REQ_statuses.front());
            if (!any_req_received) {
                break;
            } else {
                // TODO, Maybe BSend?
                int send_buf[2] = {RespType::REJECT, 0};
                MPI_Send(&send_buf, 2, MPI_INT, index,
                         MsgType::RESP, MPI_COMM_WORLD);
                MPI_Irecv(&REQ_recv_buf[index], 1, MPI_INT,
                          index, MsgType::REQ, MPI_COMM_WORLD,
                          &REQ_requests[index]);
            }
        }
    }
    if (DEBUG >= 5) {
        sleep(1);
        printf("Proc:\n");
        for(auto proc : threads) {
            printf("%d\n", proc);
        }
        printf("Requests:\n");
        for (auto r : RESP_recv_buf) {
            printf("[%d, %d]\n", r[0], r[1]);
        }
    }
    return RESP_recv_buf;
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

    // Sync barrier
    if (DEBUG >= 5)
        printf("Created %d/%d\n", rank, num_threads);
    MPI_Barrier(MPI_COMM_WORLD);

    // Init process variables
    srand(42 ^ rank); // Random lock
    int file_version = 0;
    int vote = 1;

    // Open REQ MPI_Irecv
    std::vector<int> REQ_recv_buf(num_threads, {});
    std::vector<MPI_Request> REQ_requests(num_threads, {});
    std::vector<MPI_Status> REQ_statuses{};
    for (int i = 0; i < num_threads; i++) {
        REQ_statuses.push_back({});
    }
    for (int i = 0; i < num_threads; i++) {
        MPI_Irecv(&REQ_recv_buf[i], 1, MPI_INT,
                  i, MsgType::REQ, MPI_COMM_WORLD,
                  &REQ_requests[i]);
    }

    // Start main loop
    // unsigned iteration = 0;
    // while (iteration < LOOP_LIMIT) {
    {
        // Choice random action
        // unsigned random_action = rand() % 10;
        unsigned random_action = rank;  // Debug random lock

        if (random_action == 0) {
            // READ
            get_timestamp(stdout, rank, "READ");
            send_recv_quorum(
                    "READ",
                    rank,
                    num_threads,
                    REQ_recv_buf,
                    REQ_requests,
                    REQ_statuses
                    );
        // } else if (random_action == 1) {
        //     // WRITE
        //     get_timestamp(stdout, rank, "WRITE");
        } else {
            // IDLE
            get_timestamp(stdout, rank, "IDLE");
            while (true) {
                int any_req_received = 0;
                int index;
                MPI_Testany(REQ_requests.size(), &REQ_requests.front(), &index,
                            &any_req_received, &REQ_statuses.front());
                if (!any_req_received) {
                    // printf("%d/%d: New cycle\n", rank, num_threads);
                    continue;
                } else {
                    // TODO, Maybe BSend?
                    int send_buf[2] = {RespType::ACCEPT, file_version};
                    printf("%d/%d: REQ from %d is OK\n", rank, num_threads, index);
                    MPI_Send(&send_buf, 2, MPI_INT, index,
                             MsgType::RESP, MPI_COMM_WORLD);
                    MPI_Irecv(&REQ_recv_buf[index], 1, MPI_INT,
                              index, MsgType::REQ, MPI_COMM_WORLD,
                              &REQ_requests[index]);
                }
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        // iteration++;
    }

    // printf("Proc: %d, Rand: %d\n", rank, rand() % 100);
    // fflush(stdout);

    // sleep(rand() % 5);
    MPI_Barrier(MPI_COMM_WORLD);
    get_timestamp(stdout, rank, "END");


    // Close file
    // fclose(out_file);

    // Finalize MPI
    MPI_Finalize();
    return 0;
}
