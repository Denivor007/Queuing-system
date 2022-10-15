/*Програмно реалізувати роботу ,наприклад, касового залу банку і показати
графічно стан черги клієнтів і зайнятість касирів(перукарня, гіпермаркет,
автозаправна станція...).
Умови постановки задачі (для банку):
1. Касир - N чол.
2. Час приходу клієнта моделюється за експоненціальним розподілом з
середнім часом 4хв для інтервалів між надходженнями.
3. Час обслуговування клієнта моделюється за експоненціальним
розподілом з середнім часом 0.4хв
4.Обмеження на роботу касового залу: з 04:00 по 17:00.
Необхідно розрахувати:
1. Середнє число заявок, що обслуговуються за одиницю часу кожним
касиром
2. Інтенсивність потоку заявок
3. Середню частку заявок, що обслуговуються системою
4. Інтенсивність потоку обслуговувань
*/

#include <iostream>
//#include <GL/glut.h>
#include <Windows.h>
#include <stdlib.h>
#include <math.h>

#include <iomanip>
#include <math.h>
#include <queue>
#include <vector>     
#include <time.h> 
#include <ctime>
#include <cmath>
#include <array>
#include <list>
#include "profile.h"


using namespace std;


template <typename T>
struct counter
{
    counter()
    {
        objects_created++;
        objects_alive++;
    }

    counter(const counter&)
    {
        objects_created++;
        objects_alive++;
    }

protected:
    virtual ~counter()
    {
        --objects_alive;
    }
    static int objects_created;
    static int objects_alive;
};
template <typename T> int counter<T>::objects_created(0);
template <typename T> int counter<T>::objects_alive(0);


struct Time {
    struct hours
    {
        int value;
        explicit hours(const int& val) : value(val) {}
    };
    struct minutes
    {
        int value;
        explicit minutes(const int& val) : value(val) {}
    };
    int _hours;
    int _minutes;
    int _second;

    explicit operator int() const { return  _hours * 60 + _minutes; }
    explicit operator double() const { return  _hours * 60 + _minutes; }
 

    Time(const double& total) {
        double* tot = new double;
        double sec = modf(total, tot);
        int ttl = (int)*tot;
        delete tot;

        _hours = ttl / 60;
        _minutes = ttl % 60;
        _second = (int)(sec * 60);
    }
    Time(const int& hours_, const int& minutes_, const int& second_ = 0)
        : Time(hours_ * 60 + minutes_ + second_ / 60.0)
    {
    }

    friend ostream& operator << (ostream&, const Time&);
    friend istream& operator >> (istream&, Time&);
};
ostream& operator << (ostream& ost, const Time& t) {
    if (t._hours == 0 && t._minutes == 0) {
        ost << t._second << "sec";
        return ost;
    }
    ost << setfill('0');
    ost << setw(2) << t._hours << ':'
        << setw(2) << t._minutes << ':'
        << setw(2) << t._second;
    return ost;
}
istream& operator >> (istream& ist, Time& t) {
    ist >> t._hours;
    ist.ignore(1);
    ist >> t._minutes;
    ist.ignore(1);
    ist >> t._second;
    return ist;
}





void Start() {
    srand(time(NULL));
    double never_use = rand();
}

double ExpDist(const double& lambda) {
    double Z = (double)rand() / RAND_MAX;
    return -(1 / lambda) * log(1 - Z);
}


class QueueingSystem {
public:
    struct Request;
    class Channel;
protected:
    inline static double maximum_queue_time_ = 30;
    inline static double mu_ = 1; //обработок за минуту 1/service_time
    double lambda_ = 1; //заявок за минуту 1/arrival_time
    vector<Request> requests_; // вектор заявок (время поступления заявки)
    vector<Channel>channels_; // вектор каналов (массивы ссылок, на существующие заявки в векторе заявок)
    size_t needed_ = channels_.size();
    double time_of_start_ = 240;
    double time_of_work_ = 0;

public:
    struct Request : counter<Request> {
        size_t id;
        double time_of_appearance; // время поступления заявки в систему
        double time_of_processing; // время, необходимое на обработку заявки
        double maximum_queue_time; // максимальное время, которое заявка может пробыть в очереди
        double time_in_queue = 0;  // время в очереди
        size_t requests_in_queue=0;
        Channel* channel = nullptr; // канал, в который попала заявка
        size_t index_in_channel = 0;

        Request(const double& toa)
            :counter()
            , id(objects_alive-1)
            , time_of_appearance(toa)
            , time_of_processing(ExpDist(mu_))
            , maximum_queue_time(ExpDist(maximum_queue_time_))
        {
        }
        Request(const double& toa, const double& top)
            :counter()
            , id(objects_alive-1)
            , time_of_appearance(toa)
            , time_of_processing(top)
            , maximum_queue_time(ExpDist(maximum_queue_time_))
        {
        }
        double Done() const {
            return time_of_appearance + time_of_processing + time_in_queue;
        }
        bool Processed() const{
            return !(channel == nullptr);
        }
    }; 

    class Channel {
    protected:
        vector <const Request*> requests_;
        double time_of_working_ = 0;
    public:
        size_t CountOfRequests() {
            return requests_.size();
        }
        double TimeOfWorking() {
            return time_of_working_;
        }
        double RequestsPerTime(const double& unit_of_time) {
            return requests_.size() / unit_of_time;
        }
        void push_back(Request* rq) {
            requests_.push_back(rq);
            time_of_working_ += rq->time_of_processing;
        }
        void pop_back() {
            time_of_working_ -= requests_.back()->time_of_processing;
            requests_.pop_back();
        }
        const Request* back() {
            return requests_.back();
        }
        bool empty() {
            return requests_.empty();
        }
        friend class QueueingSystem;
    };

    QueueingSystem(const size_t& N = 1,
        const size_t& Lambda = 1,
        const double& Mu = 1,
        const double& maximum_queue_time = 30)
        :lambda_(Lambda)
    {
        channels_ = vector<Channel>(N);
        mu_ = Mu;
        maximum_queue_time_ = 1/maximum_queue_time;
    }

    void SimulateClient(const size_t& time = 60) {
        time_of_work_ += time;
        double current_time_ = 0;
        while (current_time_ < time) {
            current_time_ += ExpDist(lambda_);
            requests_.push_back(Request(current_time_));
        }
    }

    void SimulateWork(
        const size_t& max_queue, 
        bool is_static = true,
        const size_t& limit = 2
    ) 
    {
        auto& rq = requests_;
        auto& ch = channels_;

        needed_ = (is_static)? ch.size() : 1;
        ch[0].push_back(&rq[0]);

        for (size_t i = 1; i < rq.size(); i++) {
            double& start = rq[i].time_of_appearance;
            double& queue_time = rq[i].time_in_queue;
            double finish = ch[0].back()->Done();

            size_t c = 0;
            for (size_t j = 1; j < ch.size() && j < needed_; j++) { //выбираем, в каком канале меньше всего очередь
                if (ch[j].empty()) { //если в канале не было заявок то, очевидно, очередь 0
                    finish = 0;
                    c = j;
                    break;
                }
                if (ch[j].back()->Done() < finish) {
                    c = j;
                    finish = ch[j].back()->Done();
                }
            }

            if (finish > start) {
                queue_time = finish - start;
            }

            if (!is_static) {
                if (start > finish && needed_ > 1) needed_--;
                if (queue_time > limit)needed_++;
            }
            if (queue_time > rq[i].maximum_queue_time) continue;

            ch[c].push_back(&rq[i]);
            rq[i].channel = &ch[c];
            rq[i].index_in_channel = ch[c].CountOfRequests()-1 ;

            if (CalculateQueueForRequest(rq[i].id) > max_queue) {
                ch[c].pop_back();
                rq[i].channel = nullptr;
                rq[i].index_in_channel = 0;
            }
        }
    }
    
    void Print(bool get_full_information = false) {
        int Y = CountOfRequests();
        int N = 0;
        for (size_t i = 0; i < CountOfChannels(); i++) {
            N += GetChannel(i).CountOfRequests();
        }

        cout.setf(ios::right);
        cout << setfill(' ');
        cout << setw(40) << "SYSTEM UPTIME - " << Time(time_of_work_) << "\n";
        cout << setfill(' ');
        cout << setw(40) << "TOTAL REQUESTS - " << Y << "\n";
        cout << setfill(' ');
        cout << setw(40) << "REQUESTS PROCESSED - " << N << "\n";
        cout << setfill(' ');
        cout << setw(40) << "RAW REQUESTS - " << Y - N << "\n";
        cout << setfill(' ');
        cout << setw(40) << "PERCENTAGE OF PROCESSED REQUESTS - " << setprecision(4) << 100.0 / Y * N << "%\n";
        cout << setfill(' ');
        cout << setw(40) << "REQUESTS FOR EACH CHANNEL:  \n";

        for (size_t i = 0; i < CountOfChannels(); i++)
        {
            size_t count = GetChannel(i).CountOfRequests();
            cout << setw(35) << "N*" << setw(2) << i << " - " << count << "(" << 100.0 / N * count << "%)\n";
        }

        if (get_full_information) {
            cout << "\n";
            cout << "REQUEST TABLE (number / time of arrival / time in queue / time in processing / done):  \n";
            cout << "------------------------------------------------------------------\n";
            for (auto rq : requests_)
            {
                cout << boolalpha;
                cout << setfill(' ');
                cout << "N*" << rq.id << "(" << setw(5) << rq.Processed() << ")"
                    << "- Start: " << setw(10) << setprecision(5) << Time(rq.time_of_appearance + time_of_start_) << setfill(' ');
                cout.setf(ios::right);
                cout << " | " << setw(5) << Time(rq.time_in_queue)
                    << " | " << setw(5) << Time(rq.time_of_processing)
                    << " | Done: " << setw(5) << Time(rq.Done() + time_of_start_)
                    << " | "<< CalculateQueueForRequest(rq.id) << "\n";
            }
        }
    }

    Channel GetChannel(const size_t& idx) const {
        return channels_[idx];
    }

    Request GetRequest(const size_t& idx) { // вернуть запрос
        if (idx >= requests_.size()) throw - 1;
        return requests_[idx];
    }

    vector<Request> GetRequests() {
        return requests_;
    }


    size_t CalculateQueueForRequest(size_t idx) { //O(N)
        if (idx >= requests_.size()) throw - 1;
        if (idx == 0) return 0;
        if (!requests_[idx].Processed()) return 0;
        double &time = requests_[idx].time_of_appearance; //время поступления текущей записи
        Channel &channel = *requests_[idx].channel; // канал, в который попала предыдущая запись
        size_t& idx_c = requests_[idx].index_in_channel; // индекс этой записи в канале 
        size_t count = 0; // кол-во людей в очереди
        for (size_t i = idx_c;  channel.requests_[i]->Done() > time && i > 0; --i)count++;
        return count-1;
    }

    size_t CountOfRequests() {
        return requests_.size();
    }

    size_t CountOfChannels(bool is_adaptive = false) {
        size_t i = 0;
        if (is_adaptive)return channels_.size();
        for (i; i < channels_.size(); i++)
        {
            if (channels_[i].empty()){
                break;
            }
        }
        return i;
    }
};




 int main() {


    LOG_DURATION("the simulation worked")
    Start();

    size_t count_of_channel = 4;
    double lambda = 1;
    double mu = 1;
    QueueingSystem sys(count_of_channel, lambda, mu);
    
    int wrk = 180 * 1;
    sys.SimulateClient(wrk);
    sys.SimulateWork(3,false);
    
    
    sys.Print(true);
  
}