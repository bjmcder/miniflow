#ifndef __TIMESTEPPER_HPP
#define __TIMESTEPPER_HPP

template<typename T>
class TimeStepper{


    private:
        T _current_time;
        T _dt;
        T _max_time;
        T _step_factor;

    public:

        /**
         * 
        */
       TimeStepper(){}

        /**
         * 
        */
        TimeStepper(const T& t_max, const T& tau, const T& dt_0 = 0){

            _max_time = t_max;
            _step_factor = tau;
            _dt = dt_0;
            _current_time = 0;

        }

        /**
         * 
        */
        const T& dt(){
            return _dt;
        }

        /**
         * 
        */
        const T& current_time(){
            return _current_time;
        }

        /**
         * 
        */
        const T& max_time(){
            return _max_time;
        }

        /**
         * 
        */
        void advance(const std::vector<T>& cell_sizes,
                     const T& reynolds,
                     const std::vector<T>& max_vels){

            auto tau = _step_factor;
            auto h = cell_sizes;
            auto Re = reynolds;

            std::vector<T> h2;
            h2.resize(h.size());

            for(auto& item: h2){
                item *= item;
            }

            auto A1 = 0.5*Re;
            
            auto denom = 0.0;
            for (auto& item: h2){
                denom += (1/item);
            }

            auto A2 = 1/denom;
            auto A = A1*A2;

            std::vector<T> comparands;
            auto max_components = max_vels;
            for (int i=0; i<max_components.size(); i++){
                comparands.push_back(i/max_components[i]);
            }
            
            comparands.push_back(A);

            auto dt_new = tau * (*std::min_element(comparands.begin(), comparands.end()));

            _dt = dt_new; 
        }
};

#endif