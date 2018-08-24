#include "Lock.hpp"

LockGraph::LockGraph(const size_t nb_threads_) : nb_threads(nb_threads_) {

    invalid = (nb_threads > std::thread::hardware_concurrency());

    if (!invalid){

        locks_threads = std::vector<std::atomic_flag>(nb_threads);
        locks_unitigs = std::vector<std::atomic_flag>(nb_threads * nb_locks_per_thread);

        for (auto& lck : locks_threads) lck.clear();
        for (auto& lck : locks_unitigs) lck.clear();
    }
}

LockGraph::LockGraph(LockGraph&& o) :   nb_threads(o.nb_threads), invalid(o.invalid),
                                        locks_threads(std::move(o.locks_threads)),
                                        locks_unitigs(std::move(o.locks_unitigs)) {

    o.clear();
}

LockGraph& LockGraph::operator=(LockGraph&& o){

    if (this != &o){

        clear();

        nb_threads = o.nb_threads;
        invalid = o.invalid;

        locks_threads = std::move(o.locks_threads);
        locks_unitigs = std::move(o.locks_unitigs);

        o.clear();
    }

    return *this;
}

void LockGraph::clear() {

    nb_threads = 0;
    invalid = true;

    locks_threads.clear();
    locks_unitigs.clear();
}

const size_t LockGraph::nb_locks_per_thread = 1024;
