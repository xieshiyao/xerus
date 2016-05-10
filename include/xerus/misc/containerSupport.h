// Xerus - A General Purpose Tensor Library
// Copyright (C) 2014-2016 Benjamin Huber and Sebastian Wolf. 
// 
// Xerus is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
// 
// Xerus is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Affero General Public License for more details.
// 
// You should have received a copy of the GNU Affero General Public License
// along with Xerus. If not, see <http://www.gnu.org/licenses/>.
//
// For further information on Xerus visit https://libXerus.org 
// or contact us at contact@libXerus.org.

/**
* @file
* @brief Header file for the standard contaienr support functions.
*/

#pragma once

#include <algorithm>

#include <vector>
#include <set>
#include <map>
#include <tuple>

#include <iostream>

#include "standard.h"
#include "sfinae.h"
#include "check.h"

namespace xerus {
	namespace misc {
		
		#if __GNUC__ > 4 || defined(__clang__)
			
			GENERATE_HAS_FUNCTION(count)
			GENERATE_HAS_FUNCTION(find)
		
			
			///@brief: Counts how often an element is contained in an arbitary container
			template<template<class, class...> class container_t, class item_t, class... rest_t,
				typename std::enable_if<has_count<container_t<item_t, rest_t...>, item_t>::value, int>::type = 0>
			inline size_t count(const container_t<item_t, rest_t...> &_container, const item_t &_item) { // noexcept(noexcept(_container.count(_itm)))  ?
				return _container.count(_item);
			}
			
			
			///@brief: Counts how often an element is contained in an arbitary container
			template<template<class, class...> class container_t, class item_t, class... rest_t, 
			typename std::enable_if<!has_count<container_t<item_t, rest_t...>, item_t>::value, int>::type = 0>
			inline size_t count(const container_t<item_t, rest_t...> &_container, const item_t &_item) {
				size_t count = 0;
				for(const item_t& otherItem : _container) {
					if(otherItem == _item) { count++; }
				}
				return count;
			}
			
			
			///@brief: Checks whether an arbitary container contains a certain element.
			template<template<class, class...> class container_t, class item_t, class... rest_t,
				typename std::enable_if<has_find<container_t<item_t, rest_t...>, item_t>::value, int>::type = 0>
			inline bool contains(const container_t<item_t, rest_t...> &_container, const item_t &_item) {
				return _container.find(_item) != _container.end();
			}

			
			///@brief: Checks whether an arbitary container contains a certain element.
			template<template<class, class...> class container_t, class item_t, class... rest_t,
				typename std::enable_if<!has_find<container_t<item_t, rest_t...>, item_t>::value, int>::type = 0>
			inline bool contains(const container_t<item_t, rest_t...> &_container, const item_t &_item) {
				return std::find(_container.begin(), _container.end(), _item) != _container.end();
			}

		#else 
			
			GENERATE_HAS_MEMBER(count)
			
			///@brief: Counts how often an element is contained in an arbitary container
			template<template<class, class...> class container_t, class item_t, class... rest_t, 
				typename std::enable_if<!has_member_count<container_t<item_t, rest_t...>>::value, int>::type = 0>
			inline size_t count(const container_t<item_t, rest_t...> &_container, const item_t &_item) {
				size_t count = 0;
				for(const item_t& otherItem : _container) {
					if(otherItem == _item) { count++; }
				}
				return count;
			}
			
			
			///@brief: Counts how often an element is contained in an arbitary container
			template<template<class, class...> class container_t, class item_t, class... rest_t, 
				typename std::enable_if<has_member_count<container_t<item_t, rest_t...>>::value, int>::type = 0>
			inline size_t count(const container_t<item_t, rest_t...> &_container, const item_t &_item) {
				return _container.count(_item);
			}
			
			///@brief: Checks whether an arbitary container contains a certain element.
			template<template<class, class...> class container_t, class item_t, class... rest_t, 
				typename std::enable_if<!has_member_count<container_t<item_t, rest_t...>>::value, int>::type = 0>
			inline bool contains(const container_t<item_t, rest_t...> &_container, const item_t &_item) {
				return std::find(_container.begin(), _container.end(), _item) != _container.end();
			}
			
			
			///@brief: Checks whether an arbitary container contains a certain element.
			template<template<class, class...> class container_t, class item_t, class... rest_t, 
				typename std::enable_if<has_member_count<container_t<item_t, rest_t...>>::value, int>::type = 0>
			inline bool contains(const container_t<item_t, rest_t...> &_container, const item_t &_item) {
				return _container.find(_item) != _container.end();
			}
		
		#endif
		
		///@brief:  Check whether an arbitary container contains all elemets of another arbitary container.
		template<template<class, class...> class containerA_t, template<class, class...> class containerB_t, class item_t, class... restA_t, class... restB_t>
		inline bool contains(const containerA_t<item_t, restA_t...> &_largeContainer, const containerB_t<item_t, restB_t...> &_smallContainer) {
			for(const item_t &item : _smallContainer) {
				if(!contains(_largeContainer, item)) { return false; }
			}
			return true;
		}
		
		
		///@brief: Checks whether two arbitary containers are disjunct, i.e. share no object.
		template<template<class, class...> class container_t, class item_t, class... rest_t>
		inline bool disjunct(const container_t<item_t, rest_t...>& _containerA, const container_t<item_t, rest_t...>& _containerB) {
			for(const item_t& item : _containerA) {
				if(contains(_containerB, item)) { return false; }
			}
			return true;
		}
		
		
		///@brief: Selects the maximal element from an arbitary container.
		template<template<class, class...> class container_t, class item_t, class... rest_t>
		inline item_t max(const container_t<item_t, rest_t...>& _container) {
			REQUIRE(!_container.empty(), "max() must not be invoked with empty container.");
			return *std::max_element(_container.begin(), _container.end());
		}
		

		///@brief: Selects the minimal element from an arbitary container.
		template<template<class, class...> class container_t, class item_t, class... rest_t>
		inline item_t min(const container_t<item_t, rest_t...>& _container) {
			REQUIRE(!_container.empty(), "min() must not be invoked with empty container.");
			return *std::min_element(_container.begin(), _container.end());
		}
		
		
		///@brief: Calculates the sum of all entries of an arbitary container.
		template<template<class, class...> class container_t, class item_t, class... rest_t>
		inline item_t sum(const container_t<item_t, rest_t...>& _container) {
			return std::accumulate(_container.begin(), _container.end(), item_t(0));
		}
		
		///@brief: Calculates the average of all entries of an arbitary container.
		template<template<class, class...> class container_t, class item_t, class... rest_t>
		inline item_t average(const container_t<item_t, rest_t...>& _container) {
			return sum(_container)/_container.size();
		}
		

		///@brief: Calculates the product of all entries of an arbitary container.
		template<template<class, class...> class container_t, class item_t, class... rest_t>
		inline item_t product(const container_t<item_t, rest_t...>& _container) {
			return std::accumulate(_container.begin(), _container.end(), item_t(1), std::multiplies<item_t>());
		}
		
		
		///@brief: Calculates the floating point product of all entries of an arbitary container.
		template<template<class, class...> class container_t, class item_t, class... rest_t>
		inline double fp_product(const container_t<item_t, rest_t...>& _container) {
			return std::accumulate(_container.begin(), _container.end(), double(1.0), std::multiplies<double>());
		}
		
		
		///@brief: Calculates the product the entries in the range [_first, _last).
		template<class item_t, class... rest_t>
		inline item_t product(const std::vector<item_t, rest_t...>& _container, const size_t _first, const size_t _last) {
			REQUIRE(_first <= _last && _last <= _container.size(), "Invalid range " << _first << "-" << _last << " given (Container size " << _container.size() << ")"); 
			return std::accumulate(_container.begin()+_first, _container.begin()+_last, item_t(1), std::multiplies<item_t>());
		}
		
		
		///@brief: Erases all elements specified by @a _rule from the container @a _container.
		template<class rule_t, template<class, class...> class container_t, class item_t, class... rest_t>
		inline void erase(container_t<item_t, rest_t...>& _container, const rule_t& _rule) {
			_container.erase(std::remove_if(_container.begin(), _container.end(), _rule), _container.end());
		}
		
		
		// Internal Helper function to print tuples.
		namespace internal {
			template<size_t n = 0, typename... Tp>
			inline typename std::enable_if<n+1 == sizeof...(Tp), void>::type
			print(std::ostream& _out, const std::tuple<Tp...>& t) {
				_out << std::get<n>(t);
			}

			template<size_t n = 0, typename... Tp>
			inline typename std::enable_if<n+1 < sizeof...(Tp), void>::type
			print(std::ostream& _out, const std::tuple<Tp...>& t) {
				_out << std::get<n>(t) << ", ";
				print<n + 1, Tp...>(_out, t);
			}
		}
	}
}


namespace std {
	
	///@brief Hash function for std::pair
	template<class S, class T>
	struct hash<std::pair<S, T>> {
	public:
		size_t operator()(const std::pair<S, T>& _pair) const {
			// NOTE must not be const to adhere to the standard (eg. n3690) paragraph 8.5(7):
			// "If a program calls for the default initialization of an object of a const-qualified type T, T shall be a class type with a user-provided default constructor."
			// NOTE std::hash has no user provided default constructor
			std::hash<S> hash_fn_s;
			std::hash<T> hash_fn_t;
			return (hash_fn_s(_pair.first) * 0x01000193u) ^ hash_fn_t(_pair.second);
		}
	};
	
	///@brief Add a + operator for iterators and size_t to avoid signed/unsigned errors.
	template<class IteratorType, 
		typename std::enable_if<
			std::is_same<typename std::iterator_traits<IteratorType>::difference_type, long>::value
			&& std::is_class<IteratorType>::value
		, bool>::type = true>
	IteratorType operator+(const IteratorType& _iterator, const size_t _add) {
		return _iterator + typename std::iterator_traits<IteratorType>::difference_type(_add);
	}
	
	
	///@brief: Concatenates two given cointainers.
	template<template<class, class...> class container_t, class item_t, class... rest_t>
	container_t<item_t, rest_t...> operator |(const container_t<item_t, rest_t...> & _left, const container_t<item_t, rest_t...> & _right) {
		container_t<item_t, rest_t...> both(_left);
		both.insert(both.end(), _right.begin(), _right.end());
		return both;
	}
	
	
	///@brief Allow to pipe tuples to ostreams.
	template<class... Tp>
	std::ostream& operator<<(std::ostream& _out, const std::tuple<Tp...>& _tuple) {
		_out << "<";
		xerus::misc::internal::print<0, Tp...>(_out, _tuple);
		_out << ">";
		return _out;
	}
	
	
	///@brief Allow to pipe normal containers to ostreams.
	template<class item_t, class... rest_t>
	std::ostream& operator<<(std::ostream& _out, const std::vector<item_t, rest_t...>& _container) {
		if(_container.size() == 0) { _out << "{ }"; return _out; }
		_out << "{ ";
		for(const item_t& item : _container) { _out << item << ", "; }
		_out << "\b\b }";
		return _out;
	}
	
	template<class item_t, class... rest_t>
	std::ostream& operator<<(std::ostream& _out, const std::set<item_t, rest_t...>& _container) {
		if(_container.size() == 0) { _out << "{ }"; return _out; }
		_out << "{ ";
		for(const item_t& item : _container) { _out << item << ", "; }
		_out << "\b\b }";
		return _out;
	}
	
	
	///@brief Allow to pipe std::maps to ostreams.
	template<class T, class U>
	std::ostream& operator<<(std::ostream& _out, const std::map<T,U>& _set) {
		if(_set.size() == 0) { _out << "{ }"; return _out; }
		_out << "{ ";
		for(const auto& item : _set) {  _out << "(" << item.first << ", " << item.second << "), "; }
		_out << "\b\b }";
		return _out;
	}
}
