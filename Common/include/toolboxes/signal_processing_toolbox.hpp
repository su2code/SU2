#pragma once

#include "../datatype_structure.hpp"
#include <vector>


class WindowingTools{
public:
    WindowingTools(){}
    ~WindowingTools(){}

    /*! \brief Returns the value of a windowing function given by fctIdx at CurTimeIdx with given endTimeIdx (i.e. endTimeIdx=nTimeIter, if one starts  windowing at time t =0.) */
    su2double GetWndWeight(int fctIdx, unsigned long CurTimeIdx, unsigned long endTimeIdx);

protected:
    const su2double PI_NUMBER = 4.0 * atan(1.0);	/*!< \brief Pi number. */

// Long time windows
    su2double HannWindow(unsigned long i, unsigned long endTimeIdx);

    su2double HannSquaredWindow(unsigned long i, unsigned long endTimeIdx);

    su2double BumpWindow(unsigned long i, unsigned long endTimeIdx);

// Short time windows
    /*su2double TriangleWindow(unsigned long i, unsigned long endTimeIdx){
        su2double currTime = static_cast<su2double>(i);
        su2double endTime = static_cast<su2double>(endTimeIdx);
        su2double tau = currTime/endTime;
        if(tau < 0.5) return 4*tau;
        else          return 4*(1-tau);
    }

    su2double QuadWindow(unsigned long i, unsigned long endTimeIdx){
        su2double currTime = static_cast<su2double>(i);
        su2double endTime = static_cast<su2double>(endTimeIdx);
        su2double tau = currTime/endTime;
        if(tau < 1.0/3.0) return 13.5*tau*tau;
        if(tau < 2.0/3.0) return 27.0*tau*(1-tau)-4.5;
        else              return 13.5*(1-tau)*(1-tau);
    }
    */
};

class WindowedAverage{
private:
su2double val;
unsigned long count;
//std::vector<su2double> values;


public:
WindowedAverage();

~WindowedAverage(){}

su2double Update(su2double valIn);

su2double Get();

unsigned long Count();

void Reset();

//void addValue(su2double valIn, unsigned long currentIter,unsigned long startIter = 0);

//su2double WindowedUpdate(int fctIdx); //Computes a windowed time average (integral)

private:
//Using Midpoint rule for consistentcy with adjoint solver
//su2double NoWindowing(); // == Square Window

//su2double HannWindowing();

//su2double HannSquaredWindowing();

//su2double BumpWindowing();


};
