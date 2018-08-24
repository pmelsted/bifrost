#ifndef BFG_LOCK_HPP
#define BFG_LOCK_HPP

#include <vector>

#include <thread>
#include <atomic>
#include <mutex>

class LockGraph {

    public:

        LockGraph(const size_t nb_threads);
        LockGraph(LockGraph&& o);

        LockGraph& operator=(LockGraph&& o);

        void clear();

        inline void lock_thread(const size_t thread_id){

            if (!invalid && (thread_id < nb_threads)){

                while (locks_threads[thread_id].test_and_set(std::memory_order_acquire));
            }
        }

        inline void unlock_thread(const size_t thread_id){

            if (!invalid && (thread_id < nb_threads)) locks_threads[thread_id].clear(std::memory_order_release);
        }

        inline void lock_all_threads(){

            if (!invalid){

                for (auto& lck : locks_threads){

                    while (lck.test_and_set(std::memory_order_acquire));
                }
            }
        }

        inline void unlock_all_threads(){

            if (!invalid){

                for (auto& lck : locks_threads) lck.clear(std::memory_order_release);
            }
        }

        inline void lock_unitig(const size_t unitig_id){

            if (!invalid){

                while (locks_unitigs[unitig_id % locks_unitigs.size()].test_and_set(std::memory_order_acquire));
            }
        }

        inline void unlock_unitig(const size_t unitig_id){

            if (!invalid) locks_unitigs[unitig_id % locks_unitigs.size()].clear(std::memory_order_release);
        }

    private :

        bool invalid;

        size_t nb_threads;

        static const size_t nb_locks_per_thread;

        std::vector<std::atomic_flag> locks_threads;
        std::vector<std::atomic_flag> locks_unitigs;
};

#endif
