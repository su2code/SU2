/*!
 * \file checkpointing.cpp
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

#include "../include/checkpointing.hpp"

/* Methods of Checkpointing base class */
Checkpointing_Scheme::Checkpointing_Scheme(int input_steps, int input_snaps, int input_snaps_in_RAM)
{
    steps = input_steps;
    snaps = input_snaps;
    snaps_in_RAM = input_snaps_in_RAM;

    // RAM = true, disk = false
    if (snaps_in_RAM == 0) { where_to_put = false; std::cout << "CP's stored on Disk.\n";}
    else if (snaps_in_RAM == snaps) {where_to_put = true; std::cout << "CP's stored in RAM.\n";}
    else throw std::runtime_error("No mixing of snaps in RAM and on Disk currently possible.");

    counter = 0;
    timestepping_order = 2; // here hardcoding needs to be released
    depth = 2+1;
    current_timestep = timestepping_order - 1; // -1 because 0th step is reached after the first timestep
    just_advanced = false;
    just_stored_checkpoint = false;
    primal_sweep = true;
    recompute_primal = false;
    capo = 0;
    oldcapo = 0;
    iCheckpoint = -1;
}

/* Methods of Checkpointing base class */
Checkpointing::Checkpointing(int input_steps, int input_snaps, int input_snaps_in_RAM, unsigned short input_CP_type)
{
    steps = input_steps;
    snaps = input_snaps;
    snaps_in_RAM = input_snaps_in_RAM;
    CP_type = input_CP_type;

    capo = 0;
    oldcapo = 0;
    // must be consistent with Common/include/option_structure.hpp
    if      (CP_type == 1) CP_scheme = new Everything (input_steps, input_snaps, input_snaps_in_RAM);
    else if (CP_type == 2) CP_scheme = new Equidistant(input_steps, input_snaps, input_snaps_in_RAM);
    else if (CP_type == 4) CP_scheme = new SU2_implementation(input_steps, input_snaps, input_snaps_in_RAM);
    else throw std::runtime_error("No valid checkpointing scheme.");

    std::cout <<  "End of Constructor. Steps: " << steps << std::endl;
}

Checkpointing::~Checkpointing(void) { 
    delete CP_scheme; 
    std::cout << "CP_scheme deleted." << std::endl; 
}


/* Methods of CheckpointEverything. Child class of Checkpointing_Scheme */
// DOESNT WORK ANYMORE REPLACED BY SU2 IMPLEMENTATION
ACTION::action Everything::revolve()
{
    ACTION::action whattodo;
    // std::cout << "counter: " << counter << " ";

    if (counter < (steps-1)*2) {
        if (counter%2 == 1) {
            capo++;
            oldcapo = capo - 1;
            iCheckpoint++;
            whattodo = ACTION::primal_step;
        } else {
            whattodo = ACTION::store_full_checkpoint;
        }
    }

    if (counter == (steps-1)*2) whattodo = ACTION::firsturn;

    if (counter > (steps-1)*2 && counter < (steps*2-1)*2-1) {
        if (counter%2 == 1) {
            capo--;
            iCheckpoint--;
            whattodo = ACTION::restore_full_checkpoint;
        } else {
            whattodo = ACTION::adjoint_step;
        }
    }

    if (counter == (steps*2-1)*2-1) whattodo = ACTION::terminate;

    counter++;
    return whattodo;
}

/* Methods of CheckpointEverything child class */
/* General idea is to store full checkpoints on Disk at equidistant Locations 
 * Such that the states in between two disk checkpoints can fit into Ram.
 *  */
ACTION::action Equidistant::revolve()
{
    if (counter > 60) return ACTION::terminate;
    counter++;

    // Initialization: primal_step, store_full_DISK
    //if(counter==1) return ACTION::primal_step;
    if(counter==1) {
        current_timestep++;
        just_stored_checkpoint = true;
        //iCheckpoint += timestepping_order+1;
        where_to_put=false;
        return ACTION::store_full_checkpoint;
    }
    
    if(primal_sweep) {

        if (current_timestep==steps) primal_sweep=false;

        // intermediate CP's
        // do SNAPS_IN_RAM times: primal_update, primal_step, store_single_RAM
        if((current_timestep-timestepping_order)%(snaps_in_RAM+depth) <= snaps_in_RAM && (current_timestep-timestepping_order)%(snaps_in_RAM+depth) != 0 || counter==2) {
            where_to_put=true;
            if (just_stored_checkpoint) {just_stored_checkpoint=false; return ACTION::primal_update;}
            if (!just_advanced && !just_stored_checkpoint) {current_timestep++; just_advanced=true; return ACTION::primal_step;}
            if (just_advanced && !just_stored_checkpoint) {iCheckpoint++; just_stored_checkpoint=true; just_advanced=false; return ACTION::store_single_state;}
            else throw std::runtime_error("Equidistant revolve primal sweep 1."); 
        }
        // full CP's
        // dp DEPTH (or TS_order+1) times: primal_update, primal_step, store_single_DISK
        if((current_timestep-timestepping_order)%(snaps_in_RAM+depth) > snaps_in_RAM || (current_timestep-timestepping_order)%(snaps_in_RAM+depth) == 0) {
            where_to_put=false;
            if (just_stored_checkpoint) {just_stored_checkpoint=false; return ACTION::primal_update;}
            if (!just_advanced && !just_stored_checkpoint) {current_timestep++; just_advanced=true; return ACTION::primal_step;}
            if (just_advanced && !just_stored_checkpoint) {iCheckpoint = -1; just_stored_checkpoint=true; just_advanced=false; return ACTION::store_single_state;}
        }

    } else { // Adjoint computation
        // end condition for whole computation
        if (current_timestep==timestepping_order) return ACTION::terminate;

        if(recompute_primal) { 
            // restore latest full checkpoint
            if((current_timestep-timestepping_order)%(snaps_in_RAM-depth) == 0) {
                where_to_put=false; // from Disk
                current_timestep -= snaps_in_RAM-depth; //set current timestep to last step of full CP
                iCheckpoint = -1;
                return ACTION::restore_full_checkpoint;
            }

            // do SNAPS_IN_RAM times: primal_update, primal_step, store_single_RAM
            if (1) {
                where_to_put=true;
                if (just_stored_checkpoint) {just_stored_checkpoint=false; return ACTION::primal_update;}
                if (!just_advanced && !just_stored_checkpoint) {current_timestep++; just_advanced=true; return ACTION::primal_step;}
                if (just_advanced && !just_stored_checkpoint) {iCheckpoint++; just_stored_checkpoint=true; just_advanced=false; return ACTION::store_single_state;}
                else throw std::runtime_error("Equidistant revolve primal sweep 1.");             
            }
            // Load 2 DEPTH-1 latest entries of latest DISK checkpoint
        } else {
            //std::cout << just_stored_checkpoint << std::endl;
            if (just_stored_checkpoint) {
                just_stored_checkpoint=false; 
                just_advanced=true; 
                if ((current_timestep-timestepping_order)%(snaps_in_RAM+depth) == 0 && current_timestep!=steps) {recompute_primal = true; std::cout << (current_timestep-timestepping_order)%(snaps_in_RAM-depth) << std::endl;}
                if (CP_actions.back()==ACTION::restore_single_state) current_timestep--;
                return ACTION::adjoint_step;
            }
            // do until CP_RAM=0 was loaded: adjoint_step, restore_single
            if(just_advanced && ((current_timestep-timestepping_order)%(snaps_in_RAM+depth) >= depth+1) || ((current_timestep-timestepping_order)%(snaps_in_RAM+depth) == 0))  {
                where_to_put=true;
                just_stored_checkpoint=true;
                just_advanced=false; 
                setcheck((current_timestep-timestepping_order)%(depth+snaps_in_RAM)-depth+1);
                return ACTION::restore_single_state;
            } else if (just_advanced && ((current_timestep-timestepping_order)%(snaps_in_RAM+depth) < depth)) {  // do DEPTH times: adjoint_step, restore_single
                where_to_put=false; // Load CP from disk
                just_stored_checkpoint=true;
                just_advanced=false; 
                return ACTION::restore_single_state;
            } else {
                throw std::runtime_error("Equidistant revolve adjoint reverse."); 
            }
        }
        
        return ACTION::terminate;
    }
}

/* Methods of SU2_implemementation child class */
/* Currently mimics a restart by storing freestream restart_flows in 0 and 1 for DT_2nd 
 * therefore every sim has to be 2 steps longer than the original version which doesn't store freestream  */
ACTION::action SU2_implementation::revolve()
{
    counter++;

    if(primal_sweep) {
        /*--- Always Update the primal time-stepping after a checkpoint was stored. ---*/
        if(just_stored_checkpoint) {
            just_stored_checkpoint = false;
            return ACTION::primal_update;
        }

        /*--- Primal sweep with checkpoints stored ---*/
        if (current_timestep == timestepping_order-1) {
            current_timestep++;
            just_stored_checkpoint = true;
            iCheckpoint += timestepping_order+1; // += depth+1;
            //iCheckpoint += 1; // += depth
            //if(snaps == snaps_in_RAM) where_to_put = false;
            return ACTION::store_full_checkpoint;
        }
        
        // This loop should excecute primal_step, store_single_state, primal_update consecutively
        if (current_timestep < steps-4) { // the -3 is due to dual time stepping
            if (!just_advanced) { // at steps-4 we also need to take a checkpoint
                current_timestep++;
                just_advanced = true;
                return ACTION::primal_step;
            } else {
                just_advanced = false;
                just_stored_checkpoint = true;
                iCheckpoint++;
                return ACTION::store_single_state;
            }
        }

        if(current_timestep == steps-4 && just_advanced == true) { // the last 3 state are stored in the RAM anyway
            just_advanced = false;
            just_stored_checkpoint = true;
            iCheckpoint++;
            return ACTION::store_single_state;
        }
        
        // last (3 for DT_2nd order) are in RAM anyway and therefore don't need to be saved, but should be for comparability
        if (current_timestep < steps-1) {
            current_timestep++;
            if (current_timestep < steps-1) 
                just_stored_checkpoint = true;
            return ACTION::primal_step; std::cout << "hi" << std::endl;
        }
        // first adjoint step here to have correct current_timestep
        if (current_timestep == steps-1) {
            just_advanced = true; 
            primal_sweep = false;
            iCheckpoint++; // artificially add 1 as later iCheckpoint-- is used.
            return ACTION::adjoint_step;
        }
    /*--- Reverse sweep ---*/

    } else {

        if (current_timestep > timestepping_order) {
            if(!just_advanced) {
                current_timestep--;
                just_advanced = true;
                return ACTION::adjoint_step;
            } else {
                just_advanced = false;
                iCheckpoint--;
                return ACTION::restore_single_state;
            }
        }
    }
    
    return ACTION::terminate;
}

/* int main()
{
    std::string CP_type = "TEMP_1";
    Checkpointing *r;
    //r = new Checkpointing(10, 3, 3, "hi");
    r = new Checkpointing(16, 4, 4, 2);
    //r = new Checkpointing(12, 9, 9, 4);
    ACTION::action whattodo;
    do
    {
        whattodo = r->revolve();
        r->CP_scheme->CP_actions.push_back(whattodo);
        //std::cout << r->CP_scheme->CP_actions.back() << std::endl; 
        std::cout << "counter: " << r->getcounter() << '\t' << "current_timestep: " << r->getcurrent_timestep() << '\t' << "iCheckpoint: " << r->getcheck() << '\t';

        if(whattodo == ACTION::primal_step) {
            std::cout << "primal_step " << std::endl; continue; }
        if(whattodo == ACTION::store_full_checkpoint) {
            std::cout << "store_full_checkpoint "; 
            if(r->getwhere()) {
                 std::cout << "in RAM " << r->getcheck();
            }
            std::cout <<  std::endl;
            continue;}
        if(whattodo == ACTION::store_single_state) {
            std::cout << "store_single_state ";
            if(r->getwhere()) {
                 std::cout << "in RAM " << r->getcheck();
            }
            std::cout <<  std::endl;
            continue; }
        if(whattodo == ACTION::adjoint_step) {
            std::cout << "adjoint_step " << std::endl; continue; }
        if(whattodo == ACTION::restore_full_checkpoint) {
            std::cout << "restore_full_checkpoint " << std::endl; continue; }
        if(whattodo == ACTION::restore_single_state) {
            std::cout << "restore_single_state "; 
            if(r->getwhere()) {
                 std::cout << "in RAM " << r->getcheck();
            }
            std::cout <<  std::endl;    
            continue; }
        if(whattodo == ACTION::firsturn){
            std::cout << "firsturn " << std::endl; continue; }
        if(whattodo == ACTION::terminate){
            std::cout << "terminate " << std::endl; continue; }
        if(whattodo == ACTION::primal_update){
            std::cout << "primal_update " << std::endl; continue; }

    } while (whattodo != ACTION::terminate);

    // free allocated memory
    delete r;

    return 0;
} */
