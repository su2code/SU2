/*!
 * \file DBMember.hpp
 * \brief Headers of template class to register classes to a factory
 * \author Mohamed Elwardi Fadeli
 * \version 7.0.6 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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

template<class Type, class FactoryType>
class DBMember
{
protected:
    static bool is_registered;
	// --- This forces is_registered to be initialized without explicitly
	// --- using it elsewhere
	DBMember() { is_registered; }
};

template <class Type, class FactoryType>
bool DBMember<Type, FactoryType>::is_registered = 
FactoryType::Register(Type::GetName(), Type::CreateMethod);
