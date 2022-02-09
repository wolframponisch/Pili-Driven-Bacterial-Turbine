//
//  main.cpp
//  Turbine simulation
//
//  Created by Wolfram Poenisch on 07/10/2021.
//


// TO DO
// Boost instead of normal random number generator
// mkdir relative path
// start pili only after pili_dynamics_start

#include <stdlib.h>
#include <iomanip>
#include <iostream>
#include <sys/stat.h>
#include <random>
#include <fstream>
#include <sstream>
#include <string>

/// Parameters
    // Parallelization
        //const int no_of_threads                         = 8;
    // Times and Steps
        double time_step                                = 0.00001;
        double total_time                               = 3600;
        const double pili_dynamics_time                 = 0.005;
        const double time_save_results                  = 0.1;
    // Cell properties
        const double cell_radius                        = 0.5;
        //const double cell_radius_std                    = 0.1;
        double mobility_translation                     = 2;
        double mobility_rotation                        = 4;
        double excluded_k_cell_cell                     = 5000;
        double excluded_k_cell_wall                     = 5000;
        double excluded_k_cell_turbine                  = 5000;
        const int N_pili_max                            = 30;
        int N_pili                                      = 15;
    // Pili properties
        const int pili_distribution                      = 0; // 0 homogeneously distributed, 1 randomly distributed
        const double vel_pro                             = 2;
        const double vel_ret                             = 2;
        const double F_stall                             = 180;
        const double pili_k_spring                       = 2000;
        const double pili_av_length                      = 1.5;
    // Walls
        const double wall_x_min                          = -10;
        const double wall_x_max                          =  10;
        const double wall_y_min                          = -10;
        const double wall_y_max                          =  10;
    // Grid properties
        const double grid_length_x                       = 1;
        const double grid_length_y                       = 1;
        const int grid_x_elements                        = 24;
        const int grid_y_elements                        = 24;
        const double grid_start_point_x                  = -12;
        const double grid_start_point_y                  = -12;
    // Rates
        const double rate_pili_ret                       = vel_pro / pili_av_length;
        //double rate_pili_pili_binding                    = 2;
        double rate_sub_attach_first                    = vel_pro / pili_av_length;
        double rate_sub_attach                           = 10;
        double rate_turbine_attach                       = 2;
        double rate_turbine_attach_special               = 2;
    // Detachment rates
        //double time_pili_pili_det                       = 10;
        double time_pili_sub_det1                        = 0.85;
        double time_pili_sub_det2                        = 0.04;
        double time_pili_turbine_de                      = 5;
        double time_pili_turbine_de_special              = 50;
        //double force_pili_pili_det                      = 10;
        double force_pili_sub_det1                       = 1.24;
        double force_pili_sub_det2                       = 33.8;
        double force_pili_turbine_de                     = 180;
        double force_pili_turbine_de_special             = 180;
    // Turbine features
        const int N_turbine_length_max                  = 100;
        const int N_turbine_width_max                   = 50;
        int N_turbine_length                            = 8;
        int N_turbine_width                             = 2;
        double turbine_radius                            = 0.5;  // needs to be smaller than cell radius
        double mobility_turbine                         = 0.005;
/// Settings
    // Initial conditions
        const int max_cell_number = 200;
        int cell_number                     = 40;

/// Variables
    // Times
        int Step = 1;
    // Cell positions
        double cell_pos_x[max_cell_number];
        double cell_pos_y[max_cell_number];
        double cell_disp_x[max_cell_number];
        double cell_disp_y[max_cell_number];
    // Cell forces and torques
        double cell_forces_x[max_cell_number];
        double cell_forces_y[max_cell_number];
        double cell_torques[max_cell_number];
    // Pili positions
        double pili_start_x[N_pili_max * max_cell_number];
        double pili_start_y[N_pili_max * max_cell_number];
        double pili_end_x[N_pili_max * max_cell_number];
        double pili_end_y[N_pili_max * max_cell_number];
    // Pili properties
        int set_pili_ID = 0;
        int pili_state[N_pili_max * max_cell_number];
        int pili_first_attached[N_pili_max * max_cell_number];
        //double pili_contour_length[N_pili * cell_number];
        double pili_free_length[N_pili_max * max_cell_number];
        //double pili_tail_length[N_pili * cell_number];
        int pili_cell_ID[N_pili_max * max_cell_number];
        int pili_attach_turbine_ID[N_pili_max * max_cell_number];
    // Pili forces and torques
        double pili_forces_x[N_pili_max * max_cell_number];
        double pili_forces_y[N_pili_max * max_cell_number];
        double pili_forces[N_pili_max * max_cell_number];
        double pili_torques[N_pili_max * max_cell_number];
    // Cell Grid
        int cell_grid[grid_x_elements][grid_y_elements][max_cell_number];
        int filled_cell_grid[grid_x_elements][grid_y_elements+1];
        int grid_positions_x[max_cell_number+1];
        int grid_positions_y[max_cell_number+1];
    // Turbine Grid
        int turbine_grid_positions_x[N_turbine_length_max*N_turbine_width_max+1];
        int turbine_grid_positions_y[N_turbine_length_max*N_turbine_width_max+1];
    // turbine
        double turbine_circles_x[N_turbine_length_max*N_turbine_width_max];
        double turbine_circles_y[N_turbine_length_max*N_turbine_width_max];
        double turbine_type[N_turbine_length_max*N_turbine_width_max]; // 1 == Normal, 2 == Special Adhesion
        double turbine_torque;
        double turbine_angle;
    // Random number generator
        int seed = 12534;
        std::default_random_engine generator;



void initial_cells()                                                                        // DONE
{
    /// Uniform random number generator
        std::uniform_real_distribution<double> distribution(0.0,1.0);
    /// Initialize cells
        std::cout << "Initialize cells" << std::endl;
    /*
        for (int i = 0; i < cell_number; i++)
        {
                double new_cell_x = (wall_x_min+cell_radius)
                                  + (wall_x_max - wall_x_min - 2 * cell_radius) * distribution(generator);
                double new_cell_y = (wall_y_min+cell_radius)
                                  + (wall_y_max - wall_y_min - 2 * cell_radius) * distribution(generator);
                cell_pos_x[i] = new_cell_x;
                cell_pos_y[i] = new_cell_y;
                cell_disp_x[i] = 0;
                cell_disp_y[i] = 0;
                cell_forces_x[i] = 0;
                cell_forces_y[i] = 0;
                cell_torques[i] = 0;
        };
     */
    
    int i = 0;
    while (i < cell_number)
    {
        double new_cell_x = (wall_x_min+cell_radius)
                          + (wall_x_max - wall_x_min - 2 * cell_radius) * distribution(generator);
        double new_cell_y = (wall_y_min+cell_radius)
                          + (wall_y_max - wall_y_min - 2 * cell_radius) * distribution(generator);
        if ((new_cell_x > N_turbine_width * turbine_radius)
        || (new_cell_x < -N_turbine_width * turbine_radius)
        || (new_cell_y > N_turbine_length * turbine_radius)
        || (new_cell_y < -N_turbine_length * turbine_radius))
        {
            cell_pos_x[i] = new_cell_x;
            cell_pos_y[i] = new_cell_y;
            cell_disp_x[i] = 0;
            cell_disp_y[i] = 0;
            cell_forces_x[i] = 0;
            cell_forces_y[i] = 0;
            cell_torques[i] = 0;
            i++;
        }
    }
    
};

void initial_pili()                                                                        // DONE
{
    /// Uniform random number generator
        std::uniform_real_distribution<double> distribution(0.0,1.0);
    /// Initialize pili
        std::cout << "Initialize pili" << std::endl;
        for (int i = 0; i < cell_number; i++)
        {
            for (int j = 0; j < N_pili; j++)
            {
                double angle;
                if (pili_distribution == 1)
                {
                    angle = 2 * 3.14159265358979323846 * distribution(generator);
                }
                if (pili_distribution == 0)
                {
                    angle = 2 * 3.14159265358979323846 * j/N_pili;
                }
                pili_start_x[i*N_pili+j] = cell_pos_x[i] + cell_radius * std::sin(angle);
                pili_start_y[i*N_pili+j] = cell_pos_y[i] + cell_radius * std::cos(angle);
                pili_end_x[i*N_pili+j] = pili_start_x[i*N_pili+j];
                pili_end_y[i*N_pili+j] = pili_start_y[i*N_pili+j];
                pili_state[i*N_pili+j] = 0;
                pili_free_length[i*N_pili+j] = 0;
                pili_cell_ID[i*N_pili+j] = i;
                pili_forces_x[i*N_pili+j] = 0;
                pili_forces_y[i*N_pili+j] = 0;
                pili_forces[i*N_pili+j] = 0;
                pili_torques[i*N_pili+j] = 0;
                pili_first_attached[i*N_pili+j] = 0;
            }
        };
};

void initial_turbine()
{
    // Intial turbines
    double pos_x;
    double pos_y;
    for (int i = 0; i < N_turbine_width; i++)
    {
        pos_x = -(N_turbine_width-1)*turbine_radius+i*2*turbine_radius;
        for (int j = 0; j < N_turbine_length; j++)
        {
            pos_y = -(N_turbine_length-1)*turbine_radius+j*2*turbine_radius;
            turbine_circles_x[i*N_turbine_length+j] = pos_x;
            turbine_circles_y[i*N_turbine_length+j] = pos_y;
            turbine_type[i*N_turbine_length+j] = 1;
        }
    }
    // Special turbine types
    turbine_type[0*N_turbine_length+0] = 2;
    turbine_type[0*N_turbine_length+1] = 2;
    turbine_type[0*N_turbine_length+2] = 2;
    turbine_type[1*N_turbine_length+N_turbine_length-3] = 2;
    turbine_type[1*N_turbine_length+N_turbine_length-2] = 2;
    turbine_type[1*N_turbine_length+N_turbine_length-1] = 2;
    turbine_angle = 0;
    
}

void set_forces_and_torques_to_zero()
{
    // Cell forces
        for (int i = 0; i < cell_number; i++)
        {
            cell_forces_x[i] = 0;
            cell_forces_y[i] = 0;
            cell_torques[i] = 0;
        }
    // Pili forces
        for (int i = 0; i < cell_number*N_pili; i++)
        {
            pili_forces_x[i] = 0;
            pili_forces_y[i] = 0;
            pili_forces[i] = 0;
            pili_torques[i] = 0;
        }
    // Turbine forces
        turbine_torque = 0;
    
}

void output_data(int step)                                                                            // DONE
{
/*
    /// Cell output
        std::stringstream filename;
        //filename    << "/Users/wolframpoenisch/Dropbox/Work/Projects/Bacteria_Turbine/C++ code/data/cells/Cell_data_" << std::setw(10) << std::setfill('0') << int(step) << ".txt";
        filename    << "data/cells/Cell_data_" << std::setw(10) << std::setfill('0') << int(step) << ".txt";
        std::ofstream output;
        output.open(filename.str().c_str());
        for (int i=0; i<cell_number; i++)
        {
            output  << step *   time_save_results       << " "
                    << cell_pos_x[i]                    << " "
                    << cell_pos_y[i]                    << " "
                    << cell_forces_x[i]                 << " "
                    << cell_forces_y[i]                 << " "
                    << cell_torques[i]                  << " "
                    << "\n";
        };
        output.close();
  
    /// Pili output
        std::stringstream filename2;
        //filename2    << "/Users/wolframpoenisch/Dropbox/Work/Projects/Bacteria_Turbine/C++ code/data/pili/Pili_data_" << std::setw(10) << std::setfill('0') << int(step) << ".txt";
        filename2    << "data/pili/Pili_data_" << std::setw(10) << std::setfill('0') << int(step) << ".txt";
        std::ofstream output2;
        output2.open(filename2.str().c_str());
        for (int i=0; i<cell_number * N_pili; i++)
        {
            output2 << pili_start_x[i]                    << " "
                    << pili_start_y[i]                    << " "
                    << pili_end_x[i]                 << " "
                    << pili_end_y[i]                 << " "
                    << pili_forces[i]               << " "
                    << pili_state[i]                  << " "
                    << pili_cell_ID[i] << " "
                    << "\n";
        };
        output2.close();
    
    /// Turbine output - positions
    
    std::stringstream filename3;
    //filename3    << "/Users/wolframpoenisch/Dropbox/Work/Projects/Bacteria_Turbine/C++ code/data/turbine/Turbine_data_" << std::setw(10) << std::setfill('0') << int(step) << ".txt";
    filename3    << "data/turbine/Turbine_data_" << std::setw(10) << std::setfill('0') << int(step) << ".txt";
    std::ofstream output3;
    output3.open(filename3.str().c_str());
    for (int i=0; i<N_turbine_length * N_turbine_width; i++)
    {
        output3 << turbine_circles_x[i]                    << " "
                << turbine_circles_y[i]                    << " "
                << turbine_type[i]                          << " "
                << "\n";
    };
    output3.close();
*/
    
    /// Turbine output - angles
    std::stringstream filename4;
    //filename4    << "/Users/wolframpoenisch/Dropbox/Work/Projects/Bacteria_Turbine/C++ code/data/Turbine_angle.txt";
    filename4    << "data/Turbine_angle.txt";
    std::ofstream output4;
    if (step == 0)
    {
        output4.open(filename4.str().c_str());
    }
    else
    {
        output4.open(filename4.str().c_str(), std::fstream::app);
    }
    output4 << step * time_save_results         << " "
            << turbine_angle                    << " "
            << "\n";
    output4.close();
    

};

void cell_cell_overlap_forces()
{
    
    double dist_x = 0;
    double dist_y = 0;
    double dist = 0;
    double force_x = 0;
    double force_y = 0;
    double N_compare_cells;
    int ID_x;
    int ID_y;
    int ID_compare_cell;
    
    for (int i = 0; i<cell_number; i++)
    {
        ID_x = grid_positions_x[i];
        ID_y = grid_positions_y[i];
        N_compare_cells = filled_cell_grid[ID_x][ID_y];
        for (int j = 0; j<N_compare_cells; j++)
        {
            ID_compare_cell = cell_grid[ID_x][ID_y][j+1];
            if (ID_compare_cell > i)
            {
                dist_x = cell_pos_x[i] - cell_pos_x[ID_compare_cell];
                dist_y = cell_pos_y[i] - cell_pos_y[ID_compare_cell];
                dist = std::sqrt(dist_x*dist_x + dist_y*dist_y);
                if (dist < 2 * cell_radius)
                {
                    force_x = (2 * cell_radius - dist) * excluded_k_cell_cell * dist_x / dist;
                    force_y = (2 * cell_radius - dist) * excluded_k_cell_cell * dist_y / dist;
                    cell_forces_x[i]               = cell_forces_x[i]               + force_x;
                    cell_forces_x[ID_compare_cell] = cell_forces_x[ID_compare_cell] - force_x;
                    cell_forces_y[i]               = cell_forces_y[i]               + force_y;
                    cell_forces_y[ID_compare_cell] = cell_forces_y[ID_compare_cell] - force_y;
                }
            }
        }
    }
}

void Set_grid_to_zero()
{
    for (int i = 0; i<grid_x_elements; i++)
    {
        for (int j = 0; j<grid_y_elements; j++)
        {
            filled_cell_grid[i][j] = 0;
        }
    }
}

void get_cell_grid_positions()
{
    for (int i = 0; i<cell_number; i++)
    {
        grid_positions_x[i] =  (cell_pos_x[i] - grid_start_point_x)/grid_length_x;
        grid_positions_y[i] =  (cell_pos_y[i] - grid_start_point_y)/grid_length_y;
    }
}

void sort_cells_in_grid()
{
    int ID_x;
    int ID_y;
    for (int i = 0; i<cell_number; i++)
    {
        ID_x = grid_positions_x[i];
        ID_y = grid_positions_y[i];
        filled_cell_grid[ID_x-1][ID_y-1] = filled_cell_grid[ID_x-1][ID_y-1] + 1;
        filled_cell_grid[ID_x-1][ID_y  ] = filled_cell_grid[ID_x-1][ID_y  ] + 1;
        filled_cell_grid[ID_x-1][ID_y+1] = filled_cell_grid[ID_x-1][ID_y+1] + 1;
        filled_cell_grid[ID_x  ][ID_y-1] = filled_cell_grid[ID_x  ][ID_y-1] + 1;
        filled_cell_grid[ID_x  ][ID_y  ] = filled_cell_grid[ID_x  ][ID_y  ] + 1;
        filled_cell_grid[ID_x  ][ID_y+1] = filled_cell_grid[ID_x  ][ID_y+1] + 1;
        filled_cell_grid[ID_x+1][ID_y-1] = filled_cell_grid[ID_x+1][ID_y-1] + 1;
        filled_cell_grid[ID_x+1][ID_y  ] = filled_cell_grid[ID_x+1][ID_y  ] + 1;
        filled_cell_grid[ID_x+1][ID_y+1] = filled_cell_grid[ID_x+1][ID_y+1] + 1;
        cell_grid[ID_x-1][ID_y-1][filled_cell_grid[ID_x-1][ID_y-1]] = i;
        cell_grid[ID_x-1][ID_y  ][filled_cell_grid[ID_x-1][ID_y  ]] = i;
        cell_grid[ID_x-1][ID_y+1][filled_cell_grid[ID_x-1][ID_y+1]] = i;
        cell_grid[ID_x  ][ID_y-1][filled_cell_grid[ID_x  ][ID_y-1]] = i;
        cell_grid[ID_x  ][ID_y  ][filled_cell_grid[ID_x  ][ID_y  ]] = i;
        cell_grid[ID_x  ][ID_y+1][filled_cell_grid[ID_x  ][ID_y+1]] = i;
        cell_grid[ID_x+1][ID_y-1][filled_cell_grid[ID_x+1][ID_y-1]] = i;
        cell_grid[ID_x+1][ID_y  ][filled_cell_grid[ID_x+1][ID_y  ]] = i;
        cell_grid[ID_x+1][ID_y+1][filled_cell_grid[ID_x+1][ID_y+1]] = i;
    }
}

void cell_wall_overlap_forces()
{
    for (int i = 0; i<cell_number; i++)
    {
        if (cell_pos_x[i] > wall_x_max - cell_radius)
        {
            cell_forces_x[i] = cell_forces_x[i] + (wall_x_max - cell_radius - cell_pos_x[i]) * excluded_k_cell_wall;
        }
        if (cell_pos_x[i] < wall_x_min + cell_radius)
        {
            cell_forces_x[i] = cell_forces_x[i] + (wall_x_min + cell_radius - cell_pos_x[i]) * excluded_k_cell_wall;
        }
        if (cell_pos_y[i] > wall_y_max -cell_radius)
        {
            cell_forces_y[i] = cell_forces_y[i] + (wall_y_max - cell_radius - cell_pos_y[i]) * excluded_k_cell_wall;
        }
        if (cell_pos_y[i] < wall_y_min + cell_radius)
        {
            cell_forces_y[i] = cell_forces_y[i] + (wall_y_min + cell_radius - cell_pos_y[i]) * excluded_k_cell_wall;
        }
    }
}

void protrude_pili()
{
    for (int i = 0; i<cell_number*N_pili; i++)
    {
        if (pili_state[i] == 0)
        {
            pili_free_length[i] = pili_free_length[i] + vel_pro * time_step;
        }
    }
}

void get_pili_end_point()
{
    double dist_x;
    double dist_y;
    int cell_ID;
    for (int i = 0; i<cell_number*N_pili; i++)
    {
        if ((pili_state[i] == 0) || (pili_state[i] == 3))
        {
            cell_ID = pili_cell_ID[i];
            dist_x = pili_start_x[i] - cell_pos_x[cell_ID];
            dist_y = pili_start_y[i] - cell_pos_y[cell_ID];
            pili_end_x[i] = pili_start_x[i] + pili_free_length[i] * dist_x / cell_radius;
            pili_end_y[i] = pili_start_y[i] + pili_free_length[i] * dist_y / cell_radius;
        }
    }
}

void attach_pili_to_substrate()
{
    for (int i = 0; i<cell_number*N_pili; i++)
    {
        if (pili_state[i] == 0)
        {
            if (pili_first_attached[i] == 0)
            {
                std::exponential_distribution<> distribution(rate_sub_attach_first);
                if (distribution(generator) < pili_dynamics_time)
                {
                    pili_state[i] = 1;
                    pili_first_attached[i] = 1;
                };
            }
            else
            {
                std::exponential_distribution<> distribution(rate_sub_attach);
                if (distribution(generator) < pili_dynamics_time)
                {
                    pili_state[i] = 1;
                    pili_first_attached[i] = 1;
                };
            }
            
        }
    }
}

void retract_pili()
{
    for (int i = 0; i<cell_number*N_pili; i++)
    {
        if ((pili_state[i] == 1) || (pili_state[i]==2))
        {
            if (pili_forces[i] < F_stall)
            {
                pili_free_length[i] = pili_free_length[i] - vel_ret * time_step * (1 - pili_forces[i] / F_stall);
            }
            
        }
        if (pili_state[i] == 3)
        {
            pili_free_length[i] = pili_free_length[i] - vel_ret * time_step;
            
        }
        if (pili_free_length[i] < 0)
        {
            pili_free_length[i] = 0;
        }
    }
}

void restart_protrusion()
{
    std::uniform_real_distribution<double> distribution(0.0,1.0);
    double angle;
    for (int i = 0; i<cell_number*N_pili; i++)
    {
        if ((pili_state[i] == 3) && (pili_free_length[i] == 0))
        {
            if (pili_distribution == 1)
            {
                angle = 2 * 3.14159265358979323846 * distribution(generator);
                pili_start_x[i] = cell_pos_x[pili_cell_ID[i]] + cell_radius * std::sin(angle);
                pili_start_y[i] = cell_pos_y[pili_cell_ID[i]] + cell_radius * std::cos(angle);
            }
            pili_free_length[i] = 0;
            pili_state[i] = 0;
            pili_end_x[i] = pili_start_x[i];
            pili_end_y[i] = pili_start_y[i];
            pili_forces[i] = 0;
            pili_forces_x[i] = 0;
            pili_forces_y[i] = 0;
            pili_torques[i] = 0;
            pili_first_attached[i] = 0;
        }
    }
}

void get_pili_forces()
{
    double contour_length;
    double dx;
    double dy;
    for (int i = 0; i<cell_number*N_pili; i++)
    {
        if ((pili_state[i] == 1) || (pili_state[i]==2))
        {
            dx = pili_end_x[i] - pili_start_x[i];
            dy = pili_end_y[i] - pili_start_y[i];
            contour_length = std::sqrt(dx * dx + dy * dy);
            if (contour_length > pili_free_length[i])
            {
                pili_forces[i] = (contour_length - pili_free_length[i]) * pili_k_spring;
                pili_forces_x[i] = pili_forces[i] * dx / contour_length;
                pili_forces_y[i] = pili_forces[i] * dy / contour_length;
            }
        }
    }
}

void protrusion_to_retraction()
{
    for (int i = 0; i<cell_number*N_pili; i++)
    {
        if (pili_state[i] == 0)
        {
            std::exponential_distribution<> distribution(rate_pili_ret);
            if (distribution(generator) < pili_dynamics_time)
            {
                pili_state[i] = 3;
            };
        }
    }
}

void detach_pili_from_substrate()
{
    double rate;
    double force;
    for (int i = 0; i<cell_number*N_pili; i++)
    {
        if (pili_state[i] == 1)
        {
            force = pili_forces[i];
            rate = 1 / (time_pili_sub_det1 * std::exp(-force/force_pili_sub_det1)
                      + time_pili_sub_det2 * std::exp(-force/force_pili_sub_det2));
            std::exponential_distribution<> distribution(rate);
            if (distribution(generator) < pili_dynamics_time)
            {
                pili_state[i] = 3;
                pili_forces[i] = 0;
                pili_forces_x[i] = 0;
                pili_forces_y[i] = 0;
                pili_torques[i] = 0;
            };
        }
    }
}

void cell_pili_forces_torques()
{
    double force_x;
    double force_y;
    double dx;
    double dy;
    int ID;
    for (int i = 0; i<cell_number*N_pili; i++)
    {
        if (pili_forces[i] > 0)
        {
            ID = pili_cell_ID[i];
            force_x = pili_forces_x[i];
            force_y = pili_forces_y[i];
            dx = pili_start_x[i] - cell_pos_x[ID];
            dy = pili_start_y[i] - cell_pos_y[ID];
            cell_forces_x[ID] = cell_forces_x[ID]  + force_x;
            cell_forces_y[ID] = cell_forces_y[ID]  + force_y;
            cell_torques[ID] = cell_torques[ID] + dx * force_y - dy * force_x;
        }
    }
}

void translate_cells_and_pili()
{
    // Cells
        for (int i=0; i<cell_number; i++)
        {
            cell_pos_x[i] = cell_pos_x[i] + cell_forces_x[i] * mobility_translation * time_step;
            cell_pos_y[i] = cell_pos_y[i] + cell_forces_y[i] * mobility_translation * time_step;
        }
    // Pili start points
    int ID;
    for (int i=0; i<cell_number*N_pili; i++)
    {
        ID = pili_cell_ID[i];
        pili_start_x[i] = pili_start_x[i] + cell_forces_x[ID] * mobility_translation * time_step;
        pili_start_y[i] = pili_start_y[i] + cell_forces_y[ID] * mobility_translation * time_step;
    }
}

void rotate_cells_and_pili()
{
    double angle;
    double dx;
    double dy;
    int ID;
    double new_dx;
    double new_dy;
    for (int i=0; i<cell_number*N_pili; i++)
    {
        ID = pili_cell_ID[i];
        angle = mobility_rotation * cell_torques[ID] * time_step;
        dx = pili_start_x[i] - cell_pos_x[ID];
        dy = pili_start_y[i] - cell_pos_y[ID];
        new_dx = dx * std::cos(angle) - dy * std::sin(angle);
        new_dy = dx * std::sin(angle) + dy * std::cos(angle);
        pili_start_x[i] = cell_pos_x[ID] + new_dx;
        pili_start_y[i] = cell_pos_y[ID] + new_dy;
    }
}

void get_turbine_grid_positions()
{
    for (int i = 0; i<N_turbine_width*N_turbine_length; i++)
    {
        turbine_grid_positions_x[i] =  (turbine_circles_x[i] - grid_start_point_x)/grid_length_x;
        turbine_grid_positions_y[i] =  (turbine_circles_y[i] - grid_start_point_y)/grid_length_y;
    }
    
}

void cell_turbine_overlap_forces()
{
    int ID_x;
    int ID_y;
    double N_compare_cells;
    double dist_x;
    double dist_y;
    double dist;
    double force_x;
    double force_y;
    int ID_compare_cell;
    for (int i = 0; i<N_turbine_width*N_turbine_length; i++)
    {
        ID_x = turbine_grid_positions_x[i];
        ID_y = turbine_grid_positions_y[i];
        N_compare_cells = filled_cell_grid[ID_x][ID_y];
        for (int j = 0; j<N_compare_cells; j++)
        {
            ID_compare_cell = cell_grid[ID_x][ID_y][j+1];
            dist_x = cell_pos_x[ID_compare_cell] - turbine_circles_x[i];
            dist_y = cell_pos_y[ID_compare_cell] - turbine_circles_y[i];
            dist = std::sqrt(dist_x*dist_x + dist_y*dist_y);
            if (dist < cell_radius + turbine_radius)
            {
                force_x = (cell_radius + turbine_radius - dist) * excluded_k_cell_turbine * dist_x / dist;
                force_y = (cell_radius + turbine_radius - dist) * excluded_k_cell_turbine * dist_y / dist;
                // Cell forces
                    cell_forces_x[ID_compare_cell] = cell_forces_x[ID_compare_cell] + force_x;
                    cell_forces_y[ID_compare_cell] = cell_forces_y[ID_compare_cell] + force_y;
                // Turbine torque
                    turbine_torque = turbine_torque - ( turbine_circles_x[i] * force_y - turbine_circles_y[i] * force_x);
            }
        }
    }
}

void rotate_turbine()
{
    double angle = mobility_turbine * turbine_torque * time_step;
    //double angle = 0.00001;
    double dx;
    double dy;
    double sinangle = std::sin(angle);
    double cosangle = std::cos(angle);
    //std::cout << sinangle << ";" << cosangle << std::endl;
    turbine_angle = turbine_angle + angle;
    for (int i=0; i<N_turbine_width*N_turbine_length; i++)
    {
        dx = turbine_circles_x[i];
        dy = turbine_circles_y[i];
        turbine_circles_x[i] = dx * cosangle - dy * sinangle;
        turbine_circles_y[i] = dx * sinangle + dy * cosangle;
    }
    
    for (int i=0; i<N_pili*cell_number; i++)
    {
        if (pili_state[i] == 2)
        {
            dx = pili_end_x[i];
            dy = pili_end_y[i];
            pili_end_x[i] = dx * cosangle - dy * sinangle;
            pili_end_y[i] = dx * sinangle + dy * cosangle;
        }
    }
}

void attach_pili_to_turbine()
{
    double dist_x;
    double dist_y;
    double dist;
    double turbine_ID;
    double rate;
    for (int i=0; i<cell_number*N_pili; i++)
    {
        for (int j=0; j<N_turbine_width*N_turbine_length; j++)
        {
            dist_x = pili_end_x[i] - turbine_circles_x[j];
            dist_y = pili_end_y[i] - turbine_circles_y[j];
            dist = std::sqrt(dist_x*dist_x+dist_y*dist_y);
            if (dist < turbine_radius)
            {
                if (turbine_type[j] == 1)
                {
                    rate = rate_turbine_attach;
                }
                if (turbine_type[j] == 2)
                {
                    rate = rate_turbine_attach_special;
                }
                std::exponential_distribution<> distribution(rate);
                if (distribution(generator) < pili_dynamics_time)
                {
                    pili_state[i] = 2;
                    pili_first_attached[i] = 1;
                    pili_attach_turbine_ID[i] = j;
                };
            }
        }
    }
}

void detach_pili_from_turbine()
{
    double rate;
    double force;
    double turbine_ID;
    for (int i = 0; i<cell_number*N_pili; i++)
    {
        if (pili_state[i] == 2)
        {
            force = pili_forces[i];
            turbine_ID = pili_attach_turbine_ID[i];
            //std::cout << turbine_ID << " " << std::endl;
            if (turbine_type[pili_attach_turbine_ID[i]] == 1)
            {
                rate = (1 / time_pili_turbine_de) * std::exp(force/force_pili_turbine_de);
            }
            if (turbine_type[pili_attach_turbine_ID[i]] == 2)
            {
                rate = (1 / time_pili_turbine_de_special) * std::exp(force/force_pili_turbine_de_special);
            }
            std::exponential_distribution<> distribution(rate);
            if (distribution(generator) < pili_dynamics_time)
            {
                pili_state[i] = 3;
                pili_forces[i] = 0;
                pili_forces_x[i] = 0;
                pili_forces_y[i] = 0;
                pili_torques[i] = 0;
            };
        }
    }
}

void turbine_pili_torques()
{
    for (int i = 0; i<cell_number*N_pili; i++)
    {
        if (pili_state[i] == 2)
        {
            turbine_torque = turbine_torque + (pili_end_y[i] * pili_forces_x[i] - pili_end_x[i] * pili_forces_y[i]);
        }
    }
}

int main(int argc, const char * argv[]) {
    /// Get parameters from script - only uncomment on the PKS cluster!
        total_time                      = atof(argv[1]);
        time_step                       = atof(argv[2]);
        cell_number                     = atof(argv[3]);
        N_pili                          = atof(argv[4]);
        rate_turbine_attach             = atof(argv[5]);
        rate_turbine_attach_special     = atof(argv[6]);
        time_pili_turbine_de            = atof(argv[7]);
        time_pili_turbine_de_special    = atof(argv[8]);
        force_pili_turbine_de           = atof(argv[9]);
        N_turbine_length                = atof(argv[10]);
        N_turbine_width                 = atof(argv[11]);
        seed                            = atof(argv[12 ]);
    /// Steps
        const int Total_Steps = total_time / time_step;
        const int Pili_dynamics_Steps  = pili_dynamics_time  / time_step + 0.5;
        const int Save_Steps = time_save_results / time_step + 0.5 ;

    /// Make directories
       // mkdir("/Users/wolframpoenisch/Dropbox/Work/Projects/Bacteria_Turbine/C++ code/data", S_IRWXU);
       // mkdir("/Users/wolframpoenisch/Dropbox/Work/Projects/Bacteria_Turbine/C++ code/data/cells", S_IRWXU );
       // mkdir("/Users/wolframpoenisch/Dropbox/Work/Projects/Bacteria_Turbine/C++ code/data/pili", S_IRWXU);
       // mkdir("/Users/wolframpoenisch/Dropbox/Work/Projects/Bacteria_Turbine/C++ code/data/turbine", S_IRWXU);
        mkdir("data", S_IRWXU);
        mkdir("data/cells", S_IRWXU );
        mkdir("data/pili", S_IRWXU);
        mkdir("data/turbine", S_IRWXU);
    
    /// Random number generator
        generator.seed(seed);
    /// Initialize cells
        initial_cells();
        initial_pili();
        initial_turbine();
    /// Loop
    int Save_No = 0;
    for (int it = 0; it<=Total_Steps+1; it++)
    {

        /// Set forces to zero
            set_forces_and_torques_to_zero();
        /// Cell grid
            Set_grid_to_zero();
            get_cell_grid_positions();
            sort_cells_in_grid();
        /// Turbine grid
            get_turbine_grid_positions();
        /// Pili forces
            get_pili_forces();
        /// Pili dynamics
            protrude_pili();
            retract_pili();
            get_pili_end_point();
        /// Cell forces
          cell_pili_forces_torques();
            cell_cell_overlap_forces();
            cell_wall_overlap_forces();
            cell_turbine_overlap_forces();
        /// Turbine torque
            turbine_pili_torques();
        /// Move cells
            translate_cells_and_pili();
            rotate_cells_and_pili();
            get_pili_end_point();
        /// Move turbine
            rotate_turbine();
        /// Pili dynamics step
            if (it % Pili_dynamics_Steps == 0)
            {
                protrusion_to_retraction();
                attach_pili_to_turbine();
                attach_pili_to_substrate();
                detach_pili_from_substrate();
                detach_pili_from_turbine();
                get_pili_end_point();
                restart_protrusion();
            }
         
        /// Save output
            if (it % Save_Steps == 0)
            {
                std::cout << it << " of " << Total_Steps+1 << " Steps" << std::endl;
                output_data(Save_No);
                // output turbine data
                Save_No = Save_No + 1;
            }
    }
    return 0;
}
