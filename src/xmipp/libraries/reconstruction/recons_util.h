/***************************************************************************
 *
 * Authors:     Edgar Garduno Angeles (edgargar@ieee.org)
 * Created on: Sep 6, 2019
 *
 * Department of Computer Science, Institute for Applied Mathematics
 * and Systems Research (IIMAS), National Autonomous University of
 * Mexico (UNAM)
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/
#ifndef REC_UTIL_HH
#define REC_UTIL_HH

#include "superiorization_reconstruct_method.h"

#include <core/xmipp_program.h>

#include <assert.h>

#include "superiorization_proximity_types.h"

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
/********************** Declaration of ART Class ******************************/
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
namespace recu
{
 struct reg_R{
        uint   r,c; // spel Coordinates
        double w;   // Weight --> amount of intersection line and spel
       };
 class pixRay{
 private:
       const double eps_zero = 1e-10;
       double Dh, Dv, Dp, Pc;
       int L,M,N;
       double Lim[4];
       std::vector<double> C,Q;
 public:
       pixRay();
       pixRay(const uint xdim,const uint ydim,double sv,double sh,double sp);
       void init(const uint xdim,const uint ydim,double sv,double sh,double sp);
       std::vector<reg_R> pixray(const int np,const int nr,const std::vector<double>& LA);
       std::vector<reg_R> pixray(const double angle);
      };
}; /* virtual class for ART method */

#endif /* REC_UTIL_HH */