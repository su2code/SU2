/*!
 * \file checkpointing.hpp
 * \brief Main file of the SU2 Computational Fluid Dynamics code
 * \author T. Kattmann
 * \version 6.1.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

/* This is a checkpointing implementation which is based upon an existing code
   'revolve' by Walther & Stumm. For more information see http://www.autodiff.org/?module=Tools&tool=Treeverse%20%2F%20Revolve .
   The current implementation is mostly a clone of 'revolve', but with the intent to further adapt/change the code more
   significantly over time.  */

#pragma once

#include <string>
#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>

/*!
 * \enum action
 * \brief Return types of the CP Alg. which defines next be performed action.
 * \author T. Kattmann
 */
namespace ACTION
{
    enum action { primal_step, primal_update, store_full_checkpoint, store_single_state, firsturn,
                  adjoint_step, restore_full_checkpoint, restore_single_state, restore_after_recompute, terminate, error};
}

class CP_state
{
    public:
        ACTION::action whattodo;
        int iCheckpoint;
        int current_timestep;
        bool where_to_put;

        void print_state() {
            std::cout << "Action: " << whattodo << '\t';
            std::cout << "iCheckpoint: " << iCheckpoint << '\t';
            std::cout << "timestep: " << current_timestep << '\t';
            std::cout << "wheretoput: " << where_to_put << std::endl;
        }
};

/*!
 * \class Checkpointing
 * \brief Parent class for checkpointing.
 * \author T. Kattmann
 */
 class Checkpointing_Scheme
 {
     public:
        Checkpointing_Scheme(int input_steps, int input_snaps, int input_snaps_in_RAM);

        ~Checkpointing_Scheme() { };

        virtual ACTION::action revolve() { return ACTION::terminate; };

        /*!
         * \brief Memory (true) or disk (false) checkpoint.
         */
        bool getwhere() { return where_to_put; }

        void setwhere(bool where) {where_to_put = where; }

        /*!
         * \brief Get current (end) timestep.
         */
        int getcapo()  { return capo; }

        /*!
         * \brief Get current (beginning) timestep.
         */
        int getoldcapo() { return oldcapo; }

        /*!
         * \brief Get the current Checkpoint number.
         */
        int getcheck() { return iCheckpoint; }

        void setcheck(int iCP) { iCheckpoint = iCP; }

        int getcounter() { return counter; };

        int getcurrent_timestep() { return current_timestep; };

        void setcurrent_timestep(int timestep) { current_timestep = timestep; }

        unsigned short gettimestepping_order() { return timestepping_order; };

        void increasecounter() { counter++; };

        std::vector<int> CP_actions;  /*!< \brief Vector holding all performed Actions ins chronological order */

        std::vector<CP_state> CP_state_vector; /*!< \brief Vector holding all CP states (ACTION::action, iCheckpoint, current_timestep, where_to_put) in chronological order */

        virtual void Create_vector_with_ACTIONS() { };

     protected:
        int steps;                          /*!< \brief Number of time steps. */
        int counter;                        /*!< \brief Counts the actions taken  */
        int current_timestep;               /*!< \brief stores the physical timestep during forward and reverse sweep */
        int current_timestep_temp;          /*!< \brief For SU2-reimplement counts down timestep and cts is needed for reverse adjoint counting and reloading of timesteps */
        unsigned short timestepping_order;  /*!< \brief Number of timesteps necessary for restart, i.e. DT_2ND = 2. */
        unsigned short depth;               /*!< \brief Number of timesteps necessary for one adjoint iteration, i.e. DT_2ND = 2+1. */
        bool just_advanced;                 /*!< \brief Marks whether a new solution was computed, such that afterwards a store/restore is called */
        bool just_stored_checkpoint;        /*!< \brief Marker is a CP was stored, because a primal_update has to follow. */
        bool primal_sweep;                  /*!< \brief Marker if one is curently in the forward primal sweep */
        bool recompute_primal;              /*!< \brief Marker during adjoint run recompute primal and store CPs */
        int snaps;                          /*!< \brief Number checkpoints. */
        int snaps_in_RAM;                   /*!< \brief Number checkpoints in memory. */
        int capo;                           /*!< \brief Current end timestep. */
        int iCheckpoint;                    /*!< \brief Checkpoint number */
        int oldcapo;                        /*!< \brief Current beginning timestep. */
        bool where_to_put;                  /*!< \brief Checkpoint in RAM (true) or on disk (false) */
 };


/*!
 * \class Equidistant
 * \brief Child class for equidistant checkpointing schemes.
 * \author T. Kattmann
 */
class Equidistant : public Checkpointing_Scheme
{
    protected:
    public:
        /*!
         * \brief Constructor of the class.
         * \param[in] steps - number of time steps.
         * \param[in] snaps_in_RAM - number checkpoints that fit in memory. It is assumed that arbitrary many fit on disk.
         */
        Equidistant(int input_steps, int input_snaps, int input_snaps_in_RAM)
            : Checkpointing_Scheme(input_steps, input_snaps, input_snaps_in_RAM) { };

        /*!
         * \brief Destructor of the class.
         */
        ~Equidistant();

        /*!
         * \brief Decision function for checkpointing.
         */
        ACTION::action revolve();

        /*!
         * \brief Fill vector with all consecutive ACTIONS for the scheme. Later only a loop over
         * the vector returns the appropriate ACTIONS.
         */
        void Create_vector_with_ACTIONS();

};

/*!
 * \class Everything
 * \brief Parent class for checkpointing schemes.
 * \author T. Kattmann
 */
class Everything : public Checkpointing_Scheme
{
    public:
        /*!
         * \brief Constructor of the class.
         * \param[in] steps - number of time steps.
         * \param[in] snaps_in_RAM - number checkpoints that fit in memory. It is assumed that arbitrary many fit on disk.
         */
        Everything(int input_steps, int input_snaps, int input_snaps_in_RAM)
            : Checkpointing_Scheme(input_steps, input_snaps, input_snaps_in_RAM) { };

        /*!
         * \brief Destructor of the class.
         */
        ~Everything();

        /*!
         * \brief Returns the next action to be performed.
         */
        ACTION::action revolve();
};

/*!
 * \class SU2_implementation
 * \brief Parent class for checkpointing schemes.
 * \author T. Kattmann
 */
class SU2_implementation : public Checkpointing_Scheme
{
    public:
        /*!
         * \brief Constructor of the class.
         * \param[in] steps - number of time steps.
         * \param[in] snaps_in_RAM - number checkpoints that fit in memory. It is assumed that arbitrary many fit on disk.
         */
        SU2_implementation(int input_steps, int input_snaps, int input_snaps_in_RAM)
            : Checkpointing_Scheme(input_steps, input_snaps, input_snaps_in_RAM) { };

        /*!
         * \brief Destructor of the class.
         */
        ~SU2_implementation();

        /*!
         * \brief Returns the next action to be performed.
         */
        ACTION::action revolve();
};


/*!
 * \class Checkpointing
 * \brief Parent class for checkpointing. Interface to SU2 solver.
 * \author T. Kattmann
 */
class Checkpointing
{
    public:

        /*!
         * \brief Constructor of the class.
         * \param[in] steps - number of time steps.
         * \param[in] snaps - number checkpoints.
         */
        Checkpointing(int input_steps, int input_snaps, int input_snaps_in_RAM, unsigned short input_CP_type);

        /*!
         * \brief Destructor of the class.
         */
        ~Checkpointing(void);

        /*!
         * \brief Returns the next action to be performed.
         */
        ACTION::action revolve() { return CP_scheme->revolve(); };

        /*!
         * \brief Memory (true) or disk (false) checkpoint.
         */
        bool getwhere() { return CP_scheme->getwhere(); }

        /*!
         * \brief Get current (end) timestep.
         */
        int getcapo()  { return CP_scheme->getcapo(); }

        /*!
         * \brief Get current (beginning) timestep.
         */
        int getoldcapo() { return CP_scheme->getoldcapo(); }

        /*!
         * \brief Get the current Checkpoint number.
         */
        int getcheck() { return CP_scheme->getcheck(); }

        /*!
         * \brief Set amount of output information.
         * \param[in] inf - Int encodes info level {1,2,3}.
         */
        void set_info(int inf) { info=inf; }

        /*!
         * \brief Get level of screen information.
         */
	    int getinfo()  { return info; }

        int getsteps() { return steps; }

        int getcounter() { return CP_scheme->getcounter(); }

        int getcurrent_timestep() { return CP_scheme->getcurrent_timestep(); }

        unsigned short gettimestepping_order() { return CP_scheme->gettimestepping_order(); };

        void setcheck(int iCP) { CP_scheme->setcheck(iCP); }

        Checkpointing_Scheme *CP_scheme; /*!< \brief Checkpointing scheme object, specifies how revolve() behaves */

    protected:
        int steps;                       /*!< \brief Number of time steps. */
        int snaps;                       /*!< \brief Number checkpoints. */
        int snaps_in_RAM;                /*!< \brief Number checkpoints in memory. */
        int info;                        /*!< \brief Amount of screen information, valid {1,2,3} */
        int capo;                        /*!< \brief Current end timestep. */
        int oldcapo;                     /*!< \brief Current beginning timestep. */
        unsigned short CP_type;             /*!< \brief Checkpoint type (EQUIDISTANT, BINOMIAL, QIQIWANG, FULL) */
        ACTION::action schedule[];       /*!< \brief Array holds ordered actions that define the checkpointing */
};
