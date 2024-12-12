#include "mpi.h"
#include <string>
#include <algorithm>
#include <random>
#include <vector>
#include <array>

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <assert.h>
#include <math.h>

#include "pprint.cpp"

#define DEBUG 0


unsigned READ_QUORUM = 3;
unsigned WRITE_QUORUM = 3;


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


void send_recv_quorum(std::string quorum_type,
                      int rank,
                      int num_threads,
                      std::vector<int>& REQ_recv_buf,
                      std::vector<MPI_Request>& REQ_requests,
                      std::vector<MPI_Status>& REQ_statuses) {
    // Choice quorum type
    int buf;
    unsigned quorum;
    if (quorum_type == "READ") {
        buf = ReqType::READ;
        quorum = READ_QUORUM;
    } else {
        buf = ReqType::WRITE;
        quorum = WRITE_QUORUM;
    }

    // Get numbers of threads for quorum
    auto threads = get_random_threads(quorum, num_threads, rank);
    if (DEBUG >= 5)
        printf("%d/%d: threads = %s\n", rank, num_threads, spprint(threads).c_str());

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
        // if (DEBUG >= 5)
        //     printf("%d/%d: Wait until all responses are received\n", rank, num_threads);

        // Test that all responses are received
        int all_resp_received = 0;
        MPI_Testall(RESP_requests.size(), &RESP_requests.front(),
                    &all_resp_received, &RESP_statuses.front());
        if (all_resp_received) {
            break;
        }
        // Reject requests while they exist
        while (true) {
            int any_req_received = 0;
            int index;
            MPI_Testany(REQ_requests.size(), &REQ_requests.front(), &index,
                        &any_req_received, &REQ_statuses.front());
            if (!any_req_received) {
                break;
            } else {
                if (DEBUG >= 5)
                    printf("%d/%d: Reject REQ from %d\n", rank, num_threads, index);
                int send_buf[2] = {RespType::REJECT, 0};
                MPI_Send(&send_buf, 2, MPI_INT, index,
                         MsgType::RESP, MPI_COMM_WORLD);
                MPI_Irecv(&REQ_recv_buf[index], 1, MPI_INT,
                          index, MsgType::REQ, MPI_COMM_WORLD,
                          &REQ_requests[index]);
            }
        }
    }

    //////////// START DEBUG ////////////
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
    ///////////// END DEBUG /////////////

    // Check that quorum was reached
    int quorum_reached = 1;
    for (unsigned i = 0; i < RESP_recv_buf.size(); i++) {
        quorum_reached = quorum_reached && RESP_recv_buf[i][0];
    }

    if (quorum_reached) {
        // Get max (newest) file version
        int newest_file_version = 0;
        for (unsigned i = 0; i < RESP_recv_buf.size(); i++) {
            newest_file_version = std::max(newest_file_version, RESP_recv_buf[i][1]);
        }

        // READ or WRITE
        if (quorum_type == "READ") {
            printf("%d/%d: Succes READ file %d\n",
                   rank, num_threads, newest_file_version);
        } else {
            newest_file_version += 1; // Update file version
            printf("%d/%d: Succes WRITE file %d. New version is %d\n",
                   rank, num_threads, newest_file_version, newest_file_version);
        }

        // RETURN_BACK votes
        for (auto proc: threads) {
            int buf = newest_file_version;
            MPI_Send(&buf, 1, MPI_INT,
                     proc, MsgType::RETURN_BACK, MPI_COMM_WORLD);
        }
        return;
    } else {
        // Reject READ or WRITE
        if (quorum_type == "READ") {
            printf("%d/%d: Reject READ\n",
                   rank, num_threads);
        } else {
            printf("%d/%d: Reject WRITE\n",
                   rank, num_threads);
        }
        // RETURN_BACK votes
        for (unsigned i = 0; i < RESP_recv_buf.size(); i++) {
            if (RESP_recv_buf[i][0]) {
                int buf = 0;
                MPI_Send(&buf, 1, MPI_INT,
                         threads[i], MsgType::RETURN_BACK, MPI_COMM_WORLD);
            }
        }
    }
}


void idle(int rank,
          int num_threads,
          int& vote,
          int& file_version,
          std::vector<int>& REQ_recv_buf,
          std::vector<MPI_Request>& REQ_requests,
          std::vector<MPI_Status>& REQ_statuses) {
    int RETURN_BACK_recv_buf{};
    MPI_Request RETURN_BACK_request{};
    MPI_Status RETURN_BACK_status{};

    int idle_iteration = 0;
    while (true) {
        // If vote == 0, try return back vote
        // And update version if file was written
        if (!vote) {
            int response_received = 0;
            MPI_Test(&RETURN_BACK_request, &response_received,
                     &RETURN_BACK_status);
            if (response_received) {
                if (DEBUG >= 5)
                    printf("%d/%d: Vote was returned back\n", rank, num_threads);
                vote = 1;
                file_version = std::max(file_version, RETURN_BACK_recv_buf);
            }
        }

        // Endless IDLE protection
        if (idle_iteration > 0 && vote) {
            // if (DEBUG >= 5)
            //     printf("%d/%d: Stop IDLE\n", rank, num_threads);
            return;
        }
        idle_iteration++;

        int any_req_received = 0;
        int index;
        MPI_Testany(REQ_requests.size(), &REQ_requests.front(), &index,
                    &any_req_received, &REQ_statuses.front());
        if (!any_req_received) {
            continue;
        } else {
            if (DEBUG >= 5)
                printf("%d/%d: Received REQ\n", rank, num_threads);
            if (vote) {
                vote = 0;
                MPI_Irecv(&RETURN_BACK_recv_buf, 1, MPI_INT,
                          index, MsgType::RETURN_BACK, MPI_COMM_WORLD,
                          &RETURN_BACK_request);
                int send_buf[2] = {RespType::ACCEPT, file_version};
                MPI_Send(&send_buf, 2, MPI_INT,
                         index, MsgType::RESP, MPI_COMM_WORLD);
                printf("%d/%d: REQ from %d. ACCEPT\n",
                       rank, num_threads, index);
            } else {
                int send_buf[2] = {RespType::REJECT, file_version};
                MPI_Send(&send_buf, 2, MPI_INT,
                         index, MsgType::RESP, MPI_COMM_WORLD);
                printf("%d/%d: REQ from %d. REJECT\n",
                       rank, num_threads, index);

            }
            MPI_Irecv(&REQ_recv_buf[index], 1, MPI_INT,
                      index, MsgType::REQ, MPI_COMM_WORLD,
                      &REQ_requests[index]);
        }
    }
}

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

    // Random lock
    int random_seed = 42 ^ rank;
    srand(random_seed);
    if (DEBUG >= 5)
        printf("%d/%d: random_seed = %d\n", rank, num_threads, random_seed);

    // Init file version and vote
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
    while (true) {
        // Sleep for readable output
        usleep(1000);

        // Choice random action
        unsigned random_action = rand() % 2000;
        // unsigned random_action = rank;  // Debug random lock

        if (random_action == 0) {
            // READ
            printf("%d/%d: Want READ\n", rank, num_threads);
            send_recv_quorum("READ",
                             rank,
                             num_threads,
                             REQ_recv_buf,
                             REQ_requests,
                             REQ_statuses);
        } else if (random_action == 1) {
            // WRITE
            printf("%d/%d: Want WRITE\n", rank, num_threads);
            sleep(1);
            send_recv_quorum("WRITE",
                             rank,
                             num_threads,
                             REQ_recv_buf,
                             REQ_requests,
                             REQ_statuses);
        } else {
            // IDLE
            // printf("%d/%d: IDLE\n", rank, num_threads);
            idle(rank,
                 num_threads,
                 vote,
                 file_version,
                 REQ_recv_buf,
                 REQ_requests,
                 REQ_statuses);
        }
    }

    // Finalize MPI
    MPI_Finalize();
    return 0;
}
