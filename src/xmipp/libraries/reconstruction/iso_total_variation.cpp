/***************************************************************************
 *
 * Authors:     Edgar Garduno Angeles (edgargar@ieee.org)
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

#include "iso_total_variation.h"
#include <core/alglib/ap.h>

#include <functional>
#include <cmath>

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
/****************** Definition of Public Methods ******************************/
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
/**
**
** Computes the Isometric Total Variation
**
*/
double itv::phi(const MultidimArray<double>& v)
{
 double sum = 0.0;
 double dw,dh,dd;
 
// std::cout<<v.xdim; // "physical" horizontal limit (x direction)
// std::cout<<v.ydim; // "physical" horizontal limit (y direction)
// std::cout<<v.zdim; // "physical" horizontal limit (z direction)
 
 for(uint k=0; k < v.zdim;k++){        // Depth
     for(uint j=0;j < v.ydim;j++){     // Height
         for(uint i=0;i < v.xdim;i++){ // Width
             dw = ((i+1) < v.xdim) ? (v.data[k*v.xdim*v.ydim + j*v.xdim + i] - v.data[k*v.xdim*v.ydim + j*v.xdim + i+1]) : 0.0;
             dh = ((j+1) < v.ydim) ? (v.data[k*v.xdim*v.ydim + j*v.xdim + i] - v.data[k*v.xdim*v.ydim + (j+1)*v.xdim + i]) : 0.0;
             dd = ((k+1) < v.zdim) ? (v.data[k*v.xdim*v.ydim + j*v.xdim + i] - v.data[(k+1)*v.xdim*v.ydim + j*v.xdim + i]) : 0.0;
             sum = sum + sqrt(dw*dw + dh*dh + dd*dd);
            }
        }
    }
 
 return sum;
}

/**
**
** Computes the normalized non-ascending vector for the Isometric Total Variation
** TV(x) = SUM of the sqrt( (x(i,j,k) - x(i+1,j,k))^2 + (x(i,j,k) - x(i,j+1,k))^2 +(x(i,j,k) - x(i,j,k+1))^2 ) =
**       = SUM sqrt( (x_i - x_r)^2 + (x_i - x_u)^2 + (x_i - x_b)^2 )
** d/dx(i,j,k) TV / || d/dx(i,j,k) TV ||
**
*/
void itv::nav(const MultidimArray<double>& v, MultidimArray<double>& w)
{
 const double ZERO=pow(10,-15);
 double denom = 0.0;
 double dw,dh,dd;
 
 // std::cout<<v.xdim; // "physical" horizontal limit (x direction)
 // std::cout<<v.ydim; // "physical" horizontal limit (y direction)
 // std::cout<<v.zdim; // "physical" horizontal limit (z direction)
 
 //
 // Computing the gradient of the total variation function
 //
 memset(w.data,0,w.xdim*w.ydim*w.zdim*sizeof(double));
 for(uint k=0; k < v.zdim;k++)         // Depth
     for(uint j=0;j < v.ydim;j++)      // Height
         for(uint i=0;i < v.xdim;i++){ // Width
             //
             // First Case
             // (d/d x_i) of TV
             //
             if(i<(v.xdim-1) && j<(v.ydim-1) && k<(v.zdim-1)){
                dw = v.data[k*v.xdim*v.ydim + j*v.xdim + i] - v.data[k*v.xdim*v.ydim + j*v.xdim + i+1];
                dh = v.data[k*v.xdim*v.ydim + j*v.xdim + i] - v.data[k*v.xdim*v.ydim + (j+1)*v.xdim + i];
                dd = v.data[k*v.xdim*v.ydim + j*v.xdim + i] - v.data[(k+1)*v.xdim*v.ydim + j*v.xdim + i];
                //Computing the denominator
                denom = sqrt(dw*dw + dh*dh + dd*dd);
                if(denom > ZERO)
                   w.data[k*v.xdim*v.ydim + j*v.xdim + i] += (3*v.data[k*v.xdim*v.ydim + j*v.xdim + i] -
                                                                v.data[k*v.xdim*v.ydim + j*v.xdim + i+1] -
                                                                v.data[k*v.xdim*v.ydim + (j+1)*v.xdim + i] -
                                                                v.data[(k+1)*v.xdim*v.ydim + j*v.xdim + i])/denom;
               }
             //
             // Second Case
             // (d/d x_r) of TV (x_r is the base and not x_i)
             //
             if(i>0 && i<v.xdim && j<(v.ydim-1) && k<(v.zdim-1)){
                dw = v.data[k*v.xdim*v.ydim + j*v.xdim + i-1] - v.data[k*v.xdim*v.ydim + j*v.xdim + i];
                dh = v.data[k*v.xdim*v.ydim + j*v.xdim + i-1] - v.data[k*v.xdim*v.ydim + (j+1)*v.xdim + i-1];
                dd = v.data[k*v.xdim*v.ydim + j*v.xdim + i-1] - v.data[(k+1)*v.xdim*v.ydim + j*v.xdim + i-1];
                //Computing the denominator
                denom = sqrt(dw*dw + dh*dh + dd*dd);
                if(denom > ZERO)
                   w.data[k*v.xdim*v.ydim + j*v.xdim + i] += (v.data[k*v.xdim*v.ydim + j*v.xdim + i] - v.data[k*v.xdim*v.ydim + j*v.xdim + i-1])/denom;
               }
             //
             // Third Case
             // (d/d x_u) of TV (x_u is the base and not x_i)
             //
             if(i<(v.xdim-1) && j>0 && j<v.ydim && k<(v.zdim-1)){
                dw = v.data[k*v.xdim*v.ydim + (j-1)*v.xdim + i] - v.data[k*v.xdim*v.ydim + (j-1)*v.xdim + i+1];
                dh = v.data[k*v.xdim*v.ydim + (j-1)*v.xdim + i] - v.data[k*v.xdim*v.ydim + j*v.xdim + i];
                dd = v.data[k*v.xdim*v.ydim + (j-1)*v.xdim + i] - v.data[(k+1)*v.xdim*v.ydim + (j-1)*v.xdim + i];
                //Computing the denominator
                denom = sqrt(dw*dw + dh*dh + dd*dd);
                if(denom > ZERO)
                   w.data[k*v.xdim*v.ydim + j*v.xdim + i] += (v.data[k*v.xdim*v.ydim + j*v.xdim + i] - v.data[k*v.xdim*v.ydim + (j-1)*v.xdim + i])/denom;
               }
             //
             // Fourth Case
             // (d/d x_b) of TV (x_b is the base and not x_i)
             //
             if(i<(v.xdim-1) && j<(v.ydim-1) && k>0 && k<v.zdim){
                dw = v.data[(k-1)*v.xdim*v.ydim + j*v.xdim + i] - v.data[(k-1)*v.xdim*v.ydim + j*v.xdim + i+1];
                dh = v.data[(k-1)*v.xdim*v.ydim + j*v.xdim + i] - v.data[(k-1)*v.xdim*v.ydim + (j+1)*v.xdim + i];
                dd = v.data[(k-1)*v.xdim*v.ydim + j*v.xdim + i] - v.data[k*v.xdim*v.ydim + j*v.xdim + i];
                //Computing the denominator
                denom = sqrt(dw*dw + dh*dh + dd*dd);
                if(denom > ZERO)
                   w.data[k*v.xdim*v.ydim + j*v.xdim + i] += (v.data[k*v.xdim*v.ydim + j*v.xdim + i] - v.data[(k-1)*v.xdim*v.ydim + j*v.xdim + i])/denom;
               }
            }
 
 //
 // Failsafe & Finding the norm of the gradient (vector)
 //
 denom = 0.0;
 for(uint k=0; k < w.zdim;k++)         // Depth
     for(uint j=0;j < w.ydim;j++)      // Height
         for(uint i=0;i < w.xdim;i++){ // Width
             if(std::isnan(w.data[k*w.xdim*w.ydim + j*w.xdim + i]) || fabs(w.data[k*w.xdim*w.ydim + j*w.xdim + i])<=ZERO)
                w.data[k*v.xdim*v.ydim + j*v.xdim + i] = 0.0;
             denom += w.data[k*w.xdim*w.ydim + j*w.xdim + i]*w.data[k*w.xdim*w.ydim + j*w.xdim + i];
            }
 
 //
 // Normalizing the resulting vector
 //
 if(denom <= ZERO)
    memset(w.data,0,w.xdim*w.ydim*w.zdim*sizeof(double));
 else{
    denom = sqrt(denom);
    for(uint k=0; k < w.zdim;k++)         // Depth
        for(uint j=0;j < w.ydim;j++)      // Height
            for(uint i=0;i < w.xdim;i++){ // Width
                w.data[k*w.xdim*w.ydim + j*w.xdim + i] = -1.0 * w.data[k*w.xdim*w.ydim + j*w.xdim + i]/denom;
                if(fabs(w.data[k*w.xdim*w.ydim + j*w.xdim + i]) < ZERO)
                   w.data[k*w.xdim*w.ydim + j*w.xdim + i] = 0.0;
               }
   }
}
#undef DEBUG
