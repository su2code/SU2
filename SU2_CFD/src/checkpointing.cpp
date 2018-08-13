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

#include "../include/checkpointing.hpp"

/* Methods of Checkpointing base class */
Checkpointing_Scheme::Checkpointing_Scheme(int input_steps, int input_snaps, int input_snaps_in_RAM)
{
    steps = input_steps;
    snaps = input_snaps;
    snaps_in_RAM = input_snaps_in_RAM;

    // RAM = true, disk = false
    if (snaps_in_RAM == 0) where_to_put = false;
    else if (snaps_in_RAM == snaps) where_to_put = true;
    else throw std::runtime_error("No mixing of snaps in RAM and on Disk currently possible.");
    counter = 0;
    capo = 0;
    oldcapo = 0;
    check = 0;
}

/* Methods of Checkpointing base class */
Checkpointing::Checkpointing(int input_steps, int input_snaps, int input_snaps_in_RAM, unsigned short input_CP_type)
{
    steps = input_steps;
    snaps = input_snaps;
    snaps_in_RAM = input_snaps_in_RAM;
    CP_type = input_CP_type;

    where_to_put = true;
    capo = 0;
    oldcapo = 0;
    check = 0;
    // must be consistent with Common/include/option_structure.hpp
    if      (CP_type == 1)  CP_scheme = new Everything (input_steps, input_snaps, input_snaps_in_RAM);
    else if (CP_type == 2) CP_scheme = new Equidistant(input_steps, input_snaps, input_snaps_in_RAM);
    else throw std::runtime_error("No valid checkpointing scheme.");

    std::cout <<  "End of Constructor. Steps: " << steps << std::endl;
}

Checkpointing::~Checkpointing(void) { delete CP_scheme; }


/* Methods of CheckpointEverything. Child class of Checkpointing_Scheme */
ACTION::action Everything::revolve()
{
    ACTION::action whattodo;
    // std::cout << "counter: " << counter << " ";

    if (counter < (steps-1)*2) {
        if (counter%2 == 1) {
            capo++;
            oldcapo = capo - 1;
            check++;
            whattodo = ACTION::primal_step;
        } else {
            whattodo = ACTION::store_full_checkpoint;
        }
    }

    if (counter == (steps-1)*2) whattodo = ACTION::firsturn;

    if (counter > (steps-1)*2 && counter < (steps*2-1)*2-1) {
        if (counter%2 == 1) {
            capo--;
            check--;
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
ACTION::action Equidistant::revolve()
{
    ACTION::action whattodo;
    // std::cout << "counter: " << counter << " ";

    if (counter == 0) {whattodo = ACTION::store_full_checkpoint; counter++;}
    else {whattodo = ACTION::terminate;}
    
    return whattodo;
}

/* int main()
{
    std::string CP_type = "EVERYTHING";
    Checkpointing *r;
    //r = new Checkpointing(10, 3, 3, "hi");
    r = new Checkpointing(10, 3, 3, CP_type);
    ACTION::action whattodo;
    do
    {
        std::cout << "getcapo: " << r->getcapo() << std::endl;
        whattodo = r->revolve();
        if(whattodo == ACTION::primal_step) {
            std::cout << "advance " << std::endl; continue; }
        if(whattodo == ACTION::store_full_checkpoint) {
            std::cout << "takeshot " <<  std::endl; continue; }
        if(whattodo == ACTION::adjoint_step) {
            std::cout << "yourturn " << std::endl; continue; }
        if(whattodo == ACTION::restore_full_checkpoint) {
            std::cout << "restore " << std::endl; continue; }
        if(whattodo == ACTION::firsturn){
            std::cout << "firsturn " << std::endl; continue; }
        if(whattodo == ACTION::terminate){
            std::cout << "terminate " << std::endl; continue; }

    } while (whattodo != ACTION::terminate);

    // free allocated memory
    delete r;

    return 0;
} */
