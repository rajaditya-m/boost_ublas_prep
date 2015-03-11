//Written by Rajaditya M
#ifndef _BOOST_LU_GAUSSELIMINATION_
#define _BOOST_LU_GAUSSELIMINATION_

#include <iostream>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/exception/exception.hpp>

namespace boost { namespace numeric { namespace ublas {

	//Returns the row echelon form of the eliminated matrix 
	template<class M> 
		typename M::size_type gauss_elimination(M &m) {
			typename M::size_type singular_ = 0;
			typename M::size_type rows = m.size1();
    	typename M::size_type cols = m.size2();
			typename M::size_type n = (std::min) (rows, cols);
			for (typename M::size_type i = 0; i < n; ++i){
				typename M::size_type maxRow = i;
				typename M::value_type maxElem = std::abs(m(i, i));
				//Search for the maximum in the column 
				for (typename M::size_type k = i + 1; k < n; ++k){
					if (std::abs(m(k, i)) > maxElem){
						maxElem = std::abs(m(k, i));
						maxRow = k;
					}
				}
				//prospective stability issues can occur here ....

				if (m(maxRow, i) != typename M::value_type()){
					//Swap the maximum with the current row 
					for (typename M::size_type k = i; k < n; ++k){
						std::swap(m(i, k), m(maxRow, k));
					}
					//Make all the elements below this zero 
					for (typename M::size_type k = i + 1; k < n; ++k){
						typename M::value_type c = -m(k, i) / m(i, i);
						for (typename M::size_type j = i; j < n; ++j){
							if (i == j){
								m(k, j) = typename M::value_type();
							}
							else {
								m(k, j) += c*m(i, j);
							}
						}
					}
				}
				else {
					singular_ = i + 1;
				}
			}
#ifdef BOOST_UBLAS_TYPE_CHECK
			BOOST_UBLAS_CHECK(singular_ == typename M::size_type(), singular())
#endif
			return singular_;
    }


}}}


#endif