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
#ifndef SUPERIORIZATION_PROXIMITY_HH
#define SUPERIORIZATION_PROXIMITY_HH

#include <core/xmipp_program.h>

// template <typename R, typename ...ARGS> using function = R(*)(ARGS...);
// template< class R, class... Args > class function<R(Args...)>;
class SuperProx
{
 public:
	enum classType{ITV,WTV};

 private:
 //
 // Methods
 //
	const double L2SQ(void);
 public:
	SuperProx();
	SuperProx(String &StrType);
	void set(std::string StrType);
	const double pr(void);
/*    Affinity(Type const t);
    void setAffinity(Type const t);
    string const type()const;
    void clear();
    float const operator()(float const &f_u,float const &f_v,Adjacency::Type const &t)const;
    float const psi(float const &f_u,float const &f_v,Adjacency::Type const &t)const;
    float const g(float const &f_u,float const &f_v,Adjacency::Type const &t)const;
    float const g(float const &f_u,float const &f_v)const;
    void setg(float const val[2],Adjacency::Type const &t);
    float const h(float const &f_u,float const &f_v,Adjacency::Type const &t)const;
    float const h(float const &f_u,float const &f_v)const;
    void seth(float const val[2],Adjacency::Type const &t);
*/
 protected:
	std::function<double(void)> prox;
 private:
};
#endif /* SUPERIORIZATION_PROXIMITY_HH */
