#ifndef _PEAKS_H_
#define _PEAKS_H_

/**
 * @file    peaks.h
 * @brief   Defines Position as a way to describe a 2D position;
 *          Defines Peaks as a way to describe a 2D frequency domain entry; 
 *          Also defines Permutation as a way to describe the 2D permutation
 * @author  lixin
 * @date    11/13/2012
 */

#include "common.h"
#include <math.h>
#include <string.h>
#include "utils.h"
#include <list>

/**
 * @brief   A position in the 2D spectrum
 *
 * Note that we describe one position in the 2D frequency domain (NonIntSFFT::n_v x NonIntSFFT::n_h) by its
 * x (vertical) and y (horizontal) position, where 0<=x<NonIntSFFT::n_v, 0<=y<NonIntSFFT::n_h.
 */
template<class Position1D>
class Position {

protected:
    
    /**
     * The x-coordinate of the position,
     * where 0 <= x < NonIntSFFT::n_v
     */
    Position1D x;

    /**
     * The y-coordinate of the position, 
     * where 0 <= y < NonIntSFFT::n_h
     */
    Position1D y;

public:
    
    /*
     * Constructor; specifies the x and y coordinate.
     * @param[in]   _x      the x-coordinate; default value is 0
     * @param[in]   _y      the y-coordinate; default value is 0
     */
    Position(in Position1D _x = 0, in Position1D _y = 0) : 
        x(_x), y(_y) {}

    /**
     * Set the coordinate of the peak. 
     *
     * The function is basically doing:
     *
     * Peaks::x = (_x mod NonIntSFFT::n_v)
     *
     * Peaks::y = (_y mod NonIntSFFT::n_h)
     *
     * @param[in]   _x  the new x-coordinate
     * @param[in]   _y  the new y-coordinate
     */
    inline void setPosition(in Position1D _x, in Position1D _y) forceinline {

        //this might be not a good idea to wrap around
        //because this doesn't work for all cases
        x = fmod(_x + NonIntSFFT::n_v, NonIntSFFT::n_v);
        y = fmod(_y + NonIntSFFT::n_h, NonIntSFFT::n_h);
    }

    /**
     * Set the X and Y coordinate, without wrapping, 
     *
     * Peaks::x = _x
     *
     * Peaks::y = _y
     *
     * @param[in]   _x  the new x-coordinate
     * @param[in]   _y  the new y-coordinate
     */
    inline void setXY(in Position1D _x, in Position1D _y) forceinline {
        x = _x;
        y = _y;
    }

    /**
     * Set the coordinate of the peak. 
     *
     * @param[in]   _pos    the new position of the peak
     */
    inline void setPosition(in Position _pos) forceinline {
        x = _pos.x;
        y = _pos.y;
    }

    /**
     * Get the x coordinate
     *
     * @return  the x coordinate
     */
    inline Position1D getX() forceinline {
        return x;
    }

    /**
     * Get the y coordinate
     *
     * @return  the y coordinate
     */
    inline Position1D getY() forceinline {
        return y;
    }


    /**
     * Map from 2D position to 1D. The formula is
     *
     * index = x * n_h + y
     *
     * where index is often used as index to 1d array. 
     *
     * Note that this doesn't make sense if x and y are not integer. 
     *
     * @return  the 1D index
     */
    inline Position1D map2Index() forceinline {
        return x * NonIntSFFT::n_h + y;
    }

    /**
     * Map from 1D index to 2D. The formular is
     *
     * x = floor(index / n_h); 
     *
     * y = mod(index, n_h);
     *
     * @param[in]   index   the 1d index
     * @return      the corresponding 2d position
     */
    inline static Position<Position1D> index2Position(in int32_t index) forceinline {
        return Position<Position1D>(index / NonIntSFFT::n_h, index % NonIntSFFT::n_h);
/*
        Position1D x = floor(index / NonIntSFFT::n_h);
        Position1D y = index - x * NonIntSFFT::n_h;

        return Position<Position1D>(x, y);
*/
    }

    /**
     * Round the position to integer grid
     *
     * @return the integer grid of that position
     */
    inline Position<int32_t> roundPos() forceinline {

        return Position<int32_t>(Utils::mathModular(round(x), NonIntSFFT::n_v), Utils::mathModular(round(y), NonIntSFFT::n_h));
    }

    /**
     * The circular l1 distance between two peaks
     *
     * @param[in]   target  the other peak
     */
    inline Position1D gridDistance(const Position<Position1D> & target) forceinline {

        Position1D dist_v = fabs(target.x - x);
        Position1D dist_h = fabs(target.y - y);

        if (dist_v > (double)NonIntSFFT::n_v / 2)
            dist_v = NonIntSFFT::n_v - dist_v;

        if (dist_h > (double)NonIntSFFT::n_v / 2)
            dist_h = NonIntSFFT::n_h - dist_h;

        return dist_v + dist_h;
    }

    /**
     * Get the symmetric position
     *
     * @return the symmetric position
     */
    inline Position<Position1D> getSymmetricPosition() {

        return Position<Position1D>(fmod(-x + NonIntSFFT::n_v, NonIntSFFT::n_v), fmod(-y + NonIntSFFT::n_h, NonIntSFFT::n_h));
    }
};

template<>
inline IntPosition IntPosition::getSymmetricPosition() {
    return IntPosition(Utils::mathModular(-x, NonIntSFFT::n_v), Utils::mathModular(-y, NonIntSFFT::n_h));
}

template<> 
inline IntPosition IntPosition::index2Position(int32_t index) {
    
    return IntPosition(index / NonIntSFFT::n_h, index % NonIntSFFT::n_h);
    
}

template<> 
inline void IntPosition::setPosition(in int32_t _x, in int32_t _y) {

    x = Utils::mathModular(_x, NonIntSFFT::n_v);
    y = Utils::mathModular(_y, NonIntSFFT::n_h);
}



/**
 * @brief   A peak in the 2D frequency domain
 *
 * We first specify a peak by its position Position to the 2D spectrum, 
 * Furthermore, we denote its value by v, which is a complex value.
 *
 * @author  lixin
 * @date    11/13/2012
 *
 */
template<class Position1D>
class Peaks {

protected:

    /**
     * The position of this peak
     */
    Position<Position1D> pos;

    /**
     * The complex value of that peak
     */
    Complex v;

public:

    /**
     * Constructor: specifiy the position and value
     *
     * @param[in]   _pos    the position of the peak
     * @param[in]   _v      the value of the peak; default value is 0
     */
    Peaks(in Position<Position1D> _pos = Position<Position1D>(0, 0), in Complex _v = 0) 
        : pos(_pos), v(_v) {};

    /**
     * Constructor: specifiy the coordinates and value
     *
     * @param[in]   _x      the x coordinate of the peak
     * @param[in]   _y      the y coordinate of the peak
     * @param[in]   _v      the value of the peak; default value is 0
     */
    Peaks(in Position1D _x, in Position1D _y, in Complex _v = 0) 
        : pos(_x, _y), v(_v) {};

    /**
     * Set the position of the peak
     *
     * @param[in]   _x  the new x-coordinate
     * @param[in]   _y  the new y-coordinate
     */
    inline void setPosition(in Position1D _x, in Position1D _y) forceinline {
        pos.setPosition(_x, _y);
    }

    /**
     * Set the position of the peak
     *
     * @param[in]   _pos    the new position
     */
    inline void setPosition(in Position<Position1D> _pos) forceinline {
        pos.setPosition(_pos);
    }
    
    /**
     * Set the X and Y coordinate, without wrapping, 
     *
     * Peaks::x = _x
     *
     * Peaks::y = _y
     *
     * @param[in]   _x  the new x-coordinate
     * @param[in]   _y  the new y-coordinate
     */
    inline void setXY(in Position1D _x, in Position1D _y) forceinline {
        pos.setXY(_x, _y);
    }


    /**
     * Set the value of the peak
     *
     * @param[in]   _v the new value of the peak
     */
    inline void setValue(in Complex _v) forceinline {
        v = _v;
    }

    /**
     * Get the x coordinate
     *
     * @return  the x coordinate
     */
    inline Position1D getX() forceinline {
        return pos.getX();
    }

    /**
     * Get the y coordinate
     *
     * @return  the y coordinate
     */
    inline Position1D getY() forceinline {
        return pos.getY();
    }

    /**
     * Get the value of the peak
     *
     * @return  the value of the peak
     */
    inline Complex getV() forceinline {
        return v;
    }

    /**
     * Get the position
     */
    inline Position<Position1D> getPos() {
        return pos;
    }

    /**
     * The circular l1 distance between two peaks
     *
     * @param[in]   target  the other peak
     */
    inline Position1D gridDistance(Peaks<Position1D> & target) forceinline {

        Position<Position1D> target_pos = target.getPos();
        return pos.gridDistance(target_pos);
    }

    /**
     * Given a list of peaks, reconstruct the frequency domain
     *
     * @param[in]   peaks   the list of peaks
     * @param[out]  y       the reconstructed frequency domain
     */
    inline static void reconstructFrequencyDomain(in std::list<Peaks<Position1D> > * peaks, out ComplexPtr y) forceinline {

        memset(y, 0, sizeof(Complex) * NonIntSFFT::n);
        ComplexPtr cur_y = NULL;

        for (PeakIterator it = peaks->begin(); it != peaks->end(); ++ it) {
            Utils::generateSincTemplate(it->getX(), it->getY());
            cur_y = y;
            for (uint32_t f_v = 0; f_v < NonIntSFFT::n_v; ++ f_v)
                for (uint32_t f_h = 0; f_h < NonIntSFFT::n_h; ++ f_h) {
                    *cur_y += it->getV() * Utils::readOffSincTail(f_v, f_h);
                    cur_y ++;
                }
        }
    }

    /**
     * Given a list of peaks, reconstruct the time domain
     *
     * @param[in]   peaks   the list of peaks
     * @param[out]  x       the reconstructed time domain
     */
    inline static void reconstructTimeDomain(in std::list<Peaks<Position1D> > * peaks, out ComplexPtr x) forceinline {

        memset(x, 0, sizeof(Complex) * NonIntSFFT::n);
        ComplexPtr cur_x = NULL;

        for (PeakIterator it = peaks->begin(); it != peaks->end(); ++ it) {
            Utils::generateSinusoidTemplate(it->getX(), it->getY());
            cur_x = x;
            for (uint32_t t_v = 0; t_v < NonIntSFFT::n_v; ++ t_v) 
                for (uint32_t t_h = 0; t_h < NonIntSFFT::n_h; ++ t_h) {
                    *cur_x += it->getV() * Utils::readOffSinusoid(t_v, t_h);
                    cur_x ++;
                }

        }

    }


};



#endif
