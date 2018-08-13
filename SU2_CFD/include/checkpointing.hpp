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

#pragma once

#include <string>
#include <iostream>
#include <string>
#include <stdexcept>

/*! 
 * \enum action
 * \brief Return types of the CP Alg. which defines next be performed action.
 * \author T. Kattmann
 */
namespace ACTION
{
    enum action { primal_step, primal_update, store_full_checkpoint, store_single_state, firsturn,
                  adjoint_step, restore_full_checkpoint, restore_single_state, terminate, error};
}

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
        int getcheck() { return check; }
        
        int getcounter() { return counter; };

        void increasecounter() { counter++; };

     protected:
        int steps;                  /*!< \brief Number of time steps. */
        int counter;                /*!< \brief Counts the actions taken  */
        int snaps;                  /*!< \brief Number checkpoints. */
        int snaps_in_RAM;           /*!< \brief Number checkpoints in memory. */
        int capo;                   /*!< \brief Current end timestep. */
        int check;                  /*!< \brief Checkpoint number */
        int oldcapo;                /*!< \brief Current beginning timestep. */
        bool where_to_put;          /*!< \brief Checkpoint in RAM (true) or on disk (false) */
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

        ACTION::action revolve();

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
        ;
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
        
    protected:
        int steps;                       /*!< \brief Number of time steps. */
        int snaps;                       /*!< \brief Number checkpoints. */
        int snaps_in_RAM;                /*!< \brief Number checkpoints in memory. */
        int info;                        /*!< \brief Amount of screen information, valid {1,2,3} */
        int capo;                        /*!< \brief Current end timestep. */
        int check;                       /*!< \brief Checkpoint number */
        int oldcapo;                     /*!< \brief Current beginning timestep. */
        bool where_to_put;               /*!< \brief Checkpoint in RAM (true) or on disk (false) */
        unsigned short CP_type;             /*!< \brief Checkpoint type (EQUIDISTANT, BINOMIAL, QIQIWANG, FULL) */
        ACTION::action schedule[];       /*!< \brief Array holds ordered actions that define the checkpointing */
        Checkpointing_Scheme *CP_scheme; /*!< \brief Checkpointing scheme object, specifies how revolve() behaves */
};