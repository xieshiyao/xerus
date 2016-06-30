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
* @brief Header file for the standard container to standard ostream operators.
*/

#pragma once
#include <tuple>
#include <set>
#include <vector>
#include <map>

namespace xerus { namespace misc {
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
	
}} // namespaces xerus::misc