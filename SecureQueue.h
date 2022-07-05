#ifndef MULTI_PARALLEL_SECUREQUEUE_H
#define MULTI_PARALLEL_SECUREQUEUE_H

#include <queue>
#include <mutex>
#include <condition_variable>

template<typename Data>
class SecureQueue
{
private:
    std::queue<Data> coda;
    std::mutex mux;
    std::condition_variable cv;
public:
    void push(Data data)
    {
        std::lock_guard<std::mutex> lock(mux);
        coda.push(data);
        cv.notify_one();
    }

    Data pop()
    {
        Data val;
        std::unique_lock<std::mutex> lock(mux);
        while(coda.empty())
        {
            cv.wait(lock);
        }

        val=coda.front();
        coda.pop();

        return val;
    }

};

#endif //MULTI_PARALLEL_SECUREQUEUE_H
