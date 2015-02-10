/***********************************************************************
 *                      Kuramoto oscillators                           *
 *             Copyright (c) 2014-2015 Alex Khrabrov                   *
 ***********************************************************************
 * This program is free software. It comes without any warranty, to    *
 * the extent permitted by applicable law. You can redistribute it     *
 * and/or modify it under the terms of the Do What The Fuck You Want   *
 * To Public License, Version 2, as published by Sam Hocevar. See      *
 * license text below.                                                 *
 ***********************************************************************
 *                                                                     *
 *            DO WHAT THE FUCK YOU WANT TO PUBLIC LICENSE              *
 *                    Version 2, December 2004                         *
 *                                                                     *
 * Copyright (C) 2004 Sam Hocevar <sam@hocevar.net>                    *
 *                                                                     *
 * Everyone is permitted to copy and distribute verbatim or modified   *
 * copies of this license document, and changing it is allowed as long *
 * as the name is changed.                                             *
 *                                                                     *
 *            DO WHAT THE FUCK YOU WANT TO PUBLIC LICENSE              *
 *   TERMS AND CONDITIONS FOR COPYING, DISTRIBUTION AND MODIFICATION   *
 *                                                                     *
 *  0. You just DO WHAT THE FUCK YOU WANT TO.                          *
 *                                                                     *
 ***********************************************************************/

typedef double fp_type;

#define _USE_MATH_DEFINES
#include <cmath>
static const fp_type TWO_PI = 2.0 * M_PI;
#undef _USE_MATH_DEFINES

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <map>
#include <random>
#include <algorithm>
#include <numeric>

#include <cassert>
#include <cstdint>
#include <cstdlib>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#if defined(WIN32) || defined(_WIN32)
#define PATH_SEPARATOR "\\"
#else
#define PATH_SEPARATOR "/"
#endif

#ifdef _MSC_VER
#define FORCE_INLINE __forceinline
#else
#define FORCE_INLINE inline
#endif

//--------------------------------------------------------------------------------------------------

template<typename T>
FORCE_INLINE T mean(const std::vector<T> & v)
{
    const auto size = v.size();
    assert(size != 0);
    T sum(0);
    for(typename std::vector<T>::size_type i = 0; i < size; i++)
        sum += v[i];
    return sum / T(size);
}

static std::vector<fp_type> tmp_s;
static std::vector<fp_type> tmp_c;

template<typename T>
T calc_order_parameter(const std::vector<T> & phases)
{
    using namespace std;
    const auto size = phases.size();
    tmp_s.resize(size);
    tmp_c.resize(size);
    for(typename std::vector<T>::size_type i = 0; i < size; i++)
    {
        const T p = phases[i];
        tmp_s[i] = sin(p);
        tmp_c[i] = cos(p);
    }
    const T mean_sin = mean(tmp_s);
    const T mean_cos = mean(tmp_c);
    return sqrt(mean_cos * mean_cos + mean_sin * mean_sin);
}

template<typename T>
FORCE_INLINE T rand_range(T min, T max)
{
    assert(max > min);
    return ((T(rand()) / T(RAND_MAX)) * (max - min)) + min;
}

//--------------------------------------------------------------------------------------------------

template<typename T>
FORCE_INLINE T wrap_phase_1(T p)
{
    if(p > T(TWO_PI))
        p -= T(TWO_PI);
    else if(p < 0.0)
        p += T(TWO_PI);
    return p;
}

template<typename T>
FORCE_INLINE T wrap_phase_2(T p)
{
    while(p < 0.0 || p > T(TWO_PI))
    {
        if(p > T(TWO_PI))
            p -= T(TWO_PI);
        else if(p < 0.0)
            p += T(TWO_PI);
    }
}

//--------------------------------------------------------------------------------------------------

template<typename T>
FORCE_INLINE T kuramoto_classic(const T phi)
{
    return std::sin(phi);
}

template<int n, typename T>
FORCE_INLINE T kuramoto_n(const T phi)
{
    return std::sin(T(n) * phi);
}

template<typename T>
FORCE_INLINE T daido_original(const T phi)
{
    return std::sin(phi) +
        T(0.2) * std::cos(phi) -
        T(0.3) * std::sin(T(2.0) * phi) +
        T(0.6) * std::cos(T(2.0) * phi) +
        T(0.7) * std::sin(T(3.0) * phi) -
        T(0.4) * std::cos(T(3.0) * phi);
}

template<typename T>
FORCE_INLINE T coupling_1(const T phi)
{
    return std::sin(phi) + std::sin(T(3.0) * phi);
}

template<typename T>
FORCE_INLINE T coupling_2(const T phi)
{
    return std::sin(phi) +
           std::sin(T(3.0) * phi) +
           std::sin(T(6.0) * phi);
}

template<typename T>
FORCE_INLINE T coupling_3(const T phi)
{
    return std::sin(phi) + std::cos(T(3.0) * phi);
}

template<typename T>
FORCE_INLINE T coupling_4(const T phi)
{
    return std::cos(phi) + std::cos(T(3.0) * phi);
}

//--------------------------------------------------------------------------------------------------

template<typename T>
void write_vector(std::ostream & os, const std::vector<T> & vec)
{
    for(const auto x : vec)
        os << x << ' ';
}

template<typename T>
void dump_state
(
    const std::string & out_file_name,
    const std::vector<T> & phase,
    const std::vector<T> & vel,
    const unsigned int step,
    const unsigned int N,
    const T r,
    const T mean_phase
)
{
    using namespace std;
    ofstream out(out_file_name, ofstream::trunc);
    out << step << ' ' << N << ' ' << r << ' ' << mean_phase << '\n';
    write_vector(out, phase);
    out << '\n';
    write_vector(out, vel);
}

template<typename T, T(*H)(T)>
void run_simulation
(
    const bool quiet,                    // if true produce no output to console

    const unsigned int N,                // number of oscillators
    const unsigned int N_steps,          // simulations steps

    const T dt,                          // time step

    std::vector<T> phase,                // initial phases
    const std::vector<T> freq,           // oscillator's frequencies
    const std::vector<T> k,              // coupling matrix

    const std::string & name,            // current working set name

    const T noise,                       // noise amplitude

    const T forcing_strength,            // when != 0 all oscillators are forced by external field
    const T forcing_freq,                // external forcing frequency

    const unsigned int dump_interval,    // how often to dump frame snapshots

    const bool r_history_enabled,        // write r values to file
    const bool mean_history_enabled,     // write mean values to file
    const bool mean_vel_history_enabled, // write mean velocity values to file

    const bool freq_modulation_enabled,  // is oscillator frequency modulation enabled
    const std::vector<T> freq_ampl,      // frequency modulation amplitude
    const std::vector<T> freq_freq,      // frequency modulation frequency
    const std::vector<T> freq_offset,    // frequency modulation offset

    const bool k_modulation_enabled,     // is coupling coefficient modulation enabled
    const std::vector<T> k_ampl,         // coupling coefficient modulation amplitude
    const std::vector<T> k_freq,         // coupling coefficient modulation frequency
    const std::vector<T> k_offset,       // coupling coefficient modulation offset

    const std::string & out_dir_name,      // MUST end with directory separator (/ or \)
    const std::string & r_file_path,       // if specified r history will be written to this file
    const std::string & mean_file_path,    // if specified mean phase will be written to this file
    const std::string & mean_vel_file_path // if specified mean velocity will be written to this file
)
{
    using namespace std;

    assert(N > 0 && dt > 0 && phase.size() == N && freq.size() == N && k.size() == N*N);
    assert(freq_modulation_enabled ? (freq_ampl.size() == N && freq_freq.size() == N && freq_offset.size() == N) : true);
    assert(k_modulation_enabled ? (k_ampl.size() == N && k_freq.size() == N*N && k_offset.size() == N*N) : true);

    const T two_dt = T(2.0) * dt;

    const bool add_noise = (noise != 0.0);
    const bool forcing_enabled = (forcing_strength != 0.0);
    const bool dump_enabled = (dump_interval != 0);

    // collect global statistics during simulation
    vector<T> r_hist(N_steps);
    vector<T> mean_hist(N_steps);
    vector<T> mean_vel_hist(N_steps);

    vector<T> phase_old = phase;
    vector<T> phase_old_old = phase;

    vector<T> vel(N);

    unsigned int dump_counter = 0;

    for(unsigned int step = 0; step < N_steps; step++)
    {
        const T t = T(step) * dt;
        for(unsigned int i = 0; i < N; i++)
        {
            T f;
            if(!freq_modulation_enabled)
                f = freq[i];
            else
                f = freq[i] + freq_ampl[i] * sin(T(TWO_PI) * freq_freq[i] * t + freq_offset[i]);

            // NOTE: k[i*N + j] can be interpreted as k[i][j]
            if(!k_modulation_enabled)
            {
                const unsigned int offset = i * N;
                const T phase_old_i = phase_old[i];
                for(unsigned int j = 0; j < N; j++)
                    f += k[offset + j] * H(phase_old[j] - phase_old_i);
            }
            else
            {
                const unsigned int offset = i * N;
                const T phase_old_i = phase_old[i];
                for(unsigned int j = 0; j < N; j++)
                {
                    const unsigned int idx = offset + j;
                    f += (k[idx] + k_ampl[i] * sin(T(TWO_PI) * k_freq[idx] * t + k_offset[idx]))
                        * H(phase_old[j] - phase_old_i);
                }
            }

            if(forcing_enabled)
                f += forcing_strength * sin(forcing_freq * t - phase_old[i]);

            if(add_noise)
                f += rand_range<T>(-noise, noise);

            // calculate new phase
            const T p = phase_old[i] + f * dt;
            phase[i] = wrap_phase_1(p);
        }

        for(unsigned int i = 0; i < N; i++)
            vel[i] = (phase[i] - phase_old_old[i]) / two_dt;

        const T r = calc_order_parameter(phase);
        r_hist[step] = r;

        const T mean_phase = mean(phase);
        mean_hist[step] = mean_phase;

        mean_vel_hist[step] = mean(vel);

        if(dump_enabled)
        {
            if(dump_counter == 0)
            {
                dump_counter = dump_interval;
                const string out_file_name(out_dir_name + to_string(step) + ".txt");
                if(!quiet)
                    cout << "Saving at step: " << step << endl;
                dump_state(out_file_name, phase, vel, step, N, r, mean_phase);
            }
            dump_counter--;
        }
        else if(!(step % 100))
            if(!quiet)
                cout << "Step: " << step << endl;

        phase_old_old = phase_old;
        phase_old = phase;
    }

    if(r_history_enabled)
    {
        ofstream r_out(r_file_path, ofstream::trunc);
        write_vector(r_out, r_hist);
    }

    if(mean_history_enabled)
    {
        ofstream mean_out(mean_file_path, ofstream::trunc);
        write_vector(mean_out, mean_hist);            
    }

    if(mean_vel_history_enabled)
    {
        ofstream mean_vel_out(mean_vel_file_path, ofstream::trunc);
        write_vector(mean_vel_out, mean_vel_hist);            
    }
}

void read_preset(std::ifstream & input,
                 unsigned int & N,
                 std::vector<double> & freq,
                 std::vector<double> & phase,
                 std::vector<double> & k)
{
    using namespace std;

    uint32_t N_;
    input.read(reinterpret_cast<char*>(&N_), sizeof(N_));
    N = N_;

    freq.resize(N);
    phase.resize(N);
    k.resize(N*N);

    input.read(reinterpret_cast<char*>(freq.data()), N*sizeof(double));
    input.read(reinterpret_cast<char*>(phase.data()), N*sizeof(double));
    input.read(reinterpret_cast<char*>(k.data()), N*N*sizeof(double));
}

void read_freq_modulation_data(std::ifstream & input,
                               unsigned int N,
                               std::vector<double> & freq_ampl,
                               std::vector<double> & freq_freq,
                               std::vector<double> & freq_offset)
{
    using namespace std;

    freq_ampl.resize(N);
    freq_freq.resize(N);
    freq_offset.resize(N);

    input.read(reinterpret_cast<char*>(freq_ampl.data()), N*sizeof(double));
    input.read(reinterpret_cast<char*>(freq_freq.data()), N*sizeof(double));
    input.read(reinterpret_cast<char*>(freq_offset.data()), N*sizeof(double));
}

void read_k_modulation_data(std::ifstream & input,
                            unsigned int N,
                            std::vector<double> & k_ampl,
                            std::vector<double> & k_freq,
                            std::vector<double> & k_offset)
{
    using namespace std;

    k_ampl.resize(N);
    k_freq.resize(N*N);
    k_offset.resize(N*N);

    input.read(reinterpret_cast<char*>(k_ampl.data()), N*sizeof(double));
    input.read(reinterpret_cast<char*>(k_freq.data()), N*N*sizeof(double));
    input.read(reinterpret_cast<char*>(k_offset.data()), N*N*sizeof(double));
}

int main(int argc, char * argv[])
{
    using namespace std;
    
    ios_base::sync_with_stdio(false);
    
    try
    {
        namespace po = boost::program_options; 
        po::options_description desc("Options"); 
        desc.add_options() 
            ("help,h", "print help message")
            ("quiet,q", "produce no output to console")
            
            ("preset,p", po::value<std::string>()->required(), "preset name") 
            ("steps,s", po::value<unsigned int>()->default_value(1000), "simulation steps") 
            ("dump-interval,di", po::value<unsigned int>()->default_value(100), "data dump interval") 
            ("dt", po::value<double>()->default_value(0.01), "simulation time step")
            
            ("noise", po::value<double>()->default_value(0.0), "noise level")
            
            ("coupling,c", po::value<std::string>()->default_value("kuramoto"), "coupling type")
            
            ("forcing-strength,fs", po::value<double>()->default_value(0.0), "forcing strength")
            ("forcing-freq,ff", po::value<double>()->default_value(0.0), "forcing frequency")
            
            ("enable-k-modulation", "enable coupling coefficient modulation")
            ("enable-freq-modulation", "enable frequency modulation")
            
            ("r-file", po::value<std::string>()->implicit_value(""), "r file path")
            ("mean-phase-file", po::value<std::string>()->implicit_value(""), "mean phase file path")
            ("mean-vel-file", po::value<std::string>()->implicit_value(""), "mean velocity file path");

        po::variables_map vm; 
        
        try 
        {
            po::store(po::parse_command_line(argc, argv, desc), vm);
            
            if(vm.count("help"))
            {
                std::cout << std::endl 
                          << "kuramoto_simulation - versatile Kuramoto simulation program" << std::endl
                          << std::endl 
                          << desc << std::endl; 
                return EXIT_SUCCESS;
            }
        }
        catch(const po::error & e)
        {
            std::cerr << "ERROR: " << e.what() << std::endl << std::endl; 
            std::cerr << desc << std::endl; 
            return EXIT_FAILURE; 
        } 

        if(!vm.count("preset"))
        {
            std::cerr << "Preset name is required!" << std::endl
                      << "Hint: use -h option to display help." << std::endl; 
            return EXIT_FAILURE;
        }

        const bool quiet = vm.count("quiet");

        const string preset_name = vm["preset"].as<std::string>();
        
        const unsigned int N_steps = vm["steps"].as<unsigned int>();
        unsigned int dump_interval = vm["dump-interval"].as<unsigned int>();
        const fp_type dt = vm["dt"].as<double>();
        const string coupling_type = vm["coupling"].as<std::string>();

        const fp_type noise = vm["noise"].as<double>();
        
        const fp_type forcing_strength = vm["forcing-strength"].as<double>();
        const fp_type forcing_freq = vm["forcing-freq"].as<double>();

        const bool freq_modulation_enabled = vm.count("enable-freq-modulation");
        const bool k_modulation_enabled = vm.count("enable-k-modulation");

        unsigned int N;
        vector<fp_type> freq;
        vector<fp_type> phase;
        vector<fp_type> k;

        // read preset data
        {
            ifstream input(preset_name + ".preset", ifstream::in | ifstream::binary);

            if(!input)
            {
                cerr << "Unable to open input file: " << preset_name << ".preset" << endl;
                return EXIT_FAILURE;
            }

            input.exceptions(ifstream::failbit | ifstream::badbit);

            read_preset(input, N, freq, phase, k);
        }

        vector<fp_type> freq_ampl;
        vector<fp_type> freq_freq;
        vector<fp_type> freq_offset;

        vector<fp_type> k_ampl;
        vector<fp_type> k_freq;
        vector<fp_type> k_offset;

        if(freq_modulation_enabled)
        {
            if(!quiet)
                cout << "freq_modulation_enabled" << endl;
            ifstream input_s(preset_name + ".fm.preset", ifstream::in | ifstream::binary);
            if (!input_s)
            {
                cerr << "Unable to open frequency modulation data file: " << preset_name << ".fm.preset" << endl;
                return EXIT_FAILURE;
            }
            input_s.exceptions(ifstream::failbit | ifstream::badbit);
            read_freq_modulation_data(input_s, N, freq_ampl, freq_freq, freq_offset);
        }

        if(k_modulation_enabled)
        {
            if(!quiet)
                cout << "k_modulation_enabled" << endl;
            ifstream input_s(preset_name + ".km.preset", ifstream::in | ifstream::binary);
            if (!input_s)
            {
                cerr << "Unable to open k modulation data file: " << preset_name << ".km.preset" << endl;
                return EXIT_FAILURE;
            }
            input_s.exceptions(ifstream::failbit | ifstream::badbit);
            read_k_modulation_data(input_s, N, k_ampl, k_freq, k_offset);
        }

        if(!quiet)
            std::cout << "Data loaded." << std::endl;

        const bool r_history_enabled = vm.count("r-file");
        const bool mean_history_enabled = vm.count("mean-phase-file");
        const bool mean_vel_history_enabled = vm.count("mean-vel-file");

        string r_file_path;
        string mean_file_path;
        string mean_vel_file_path;

        string out_dir_name;
        
        // initialize dump directory
        if(dump_interval > 0)
        {
            out_dir_name = "dump_" + preset_name + PATH_SEPARATOR + "steps" + PATH_SEPARATOR;
            // std::cout << "Creating steps dir: " << out_dir_name << std::endl;
            boost::filesystem::create_directories(out_dir_name);
        }

        if(r_history_enabled)
        {
            r_file_path = vm["r-file"].as<std::string>();
            if(r_file_path.empty())
                r_file_path = out_dir_name + "r.txt";
            const string path_str = boost::filesystem::path(r_file_path).parent_path().string();
            if(!path_str.empty())
                boost::filesystem::create_directories(path_str);
        }

        if(mean_history_enabled)
        {
            mean_file_path = vm["mean-phase-file"].as<std::string>();
            if(mean_file_path.empty())
                mean_file_path = out_dir_name + "mean.txt";
            const string path_str = boost::filesystem::path(mean_file_path).parent_path().string();
            if(!path_str.empty())
                boost::filesystem::create_directories(path_str);
        }

        if(mean_vel_history_enabled)
        {
            mean_vel_file_path = vm["mean-vel-file"].as<std::string>();
            if(mean_vel_file_path.empty())
                mean_vel_file_path = out_dir_name + "mean_vel.txt";
            const string path_str = boost::filesystem::path(mean_vel_file_path).parent_path().string();
            if(!path_str.empty())
                boost::filesystem::create_directories(path_str);
        }

        //std::cout << "r_file_path: " << r_file_path << std::endl;
        //std::cout << "mean_file_path: " << mean_file_path << std::endl;
        //std::cout << "mean_vel_file_path: " << mean_vel_file_path << std::endl;

#define RUN_SIMUL(coupling)                  \
        run_simulation<fp_type, coupling>(   \
            quiet,                           \
            N, N_steps, dt,                  \
            phase, freq, k, preset_name,     \
            noise,                           \
            forcing_strength, forcing_freq,  \
            dump_interval,                   \
            r_history_enabled,               \
            mean_history_enabled,            \
            mean_vel_history_enabled,        \
            freq_modulation_enabled,         \
            freq_ampl,                       \
            freq_freq,                       \
            freq_offset,                     \
            k_modulation_enabled,            \
            k_ampl,                          \
            k_freq,                          \
            k_offset,                        \
            out_dir_name,                    \
            r_file_path,                     \
            mean_file_path,                  \
            mean_vel_file_path               \
        )

        if(coupling_type == "kuramoto")
            RUN_SIMUL(kuramoto_classic);
        else if(coupling_type == "kuramoto_2")
            RUN_SIMUL(kuramoto_n<2>);
    	else if(coupling_type == "kuramoto_3")
    		RUN_SIMUL(kuramoto_n<3>);
    	else if(coupling_type == "kuramoto_4")
    		RUN_SIMUL(kuramoto_n<4>);
    	else if(coupling_type == "kuramoto_5")
    		RUN_SIMUL(kuramoto_n<5>);
    	else if(coupling_type == "kuramoto_6")
    		RUN_SIMUL(kuramoto_n<6>);
    	else if(coupling_type == "kuramoto_7")
    		RUN_SIMUL(kuramoto_n<7>);
    	else if(coupling_type == "kuramoto_8")
    		RUN_SIMUL(kuramoto_n<8>);
    	else if(coupling_type == "kuramoto_9")
    		RUN_SIMUL(kuramoto_n<9>);
    	else if(coupling_type == "kuramoto_10")
    		RUN_SIMUL(kuramoto_n<10>);
    	else if(coupling_type == "daido")
            RUN_SIMUL(daido_original);
        else if(coupling_type == "coupling_1")
            RUN_SIMUL(coupling_1);
        else if(coupling_type == "coupling_2")
            RUN_SIMUL(coupling_2);
        else if(coupling_type == "coupling_3")
            RUN_SIMUL(coupling_3);
        else if(coupling_type == "coupling_4")
            RUN_SIMUL(coupling_4);
        else
        {
            std::cerr << "Unknown coupling type: " << coupling_type << std::endl;
            return EXIT_FAILURE;
        }
    }
    catch(std::exception& e) 
    { 
        std::cerr << "Unhandled Exception reached the top of main: " 
                  << e.what() << ", application will now exit." << std::endl; 
        return EXIT_FAILURE;
    } 

    return EXIT_SUCCESS;
}
