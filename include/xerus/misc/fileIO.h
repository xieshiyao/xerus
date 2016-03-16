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
 * @brief Header file for templates to store and restore objects from / to files / streams.
 */


#pragma once

#include <iostream>
#include <fstream>
#include "stringUtilities.h"
#include "check.h"

namespace xerus { namespace misc {
	
	/**
	 * @brief possible file formats for tensor storage
	 * @details 
	 * BINARY: binary output, exact storage
	 * TSV: tab-separated-value lists
	 */
	enum FileFormat { BINARY, TSV };
	
	// ---------------------------------------- write_to_stream -------------------------------------------------------------
	
	template<class T, typename std::enable_if<std::is_arithmetic<T>::value, bool>::type = true>
	void write_to_stream(std::ostream& _stream, const T &_value, FileFormat _format, const char _space) {
		if(_format == FileFormat::TSV) {
			_stream << _value << _space;
		} else {
			_stream.write(reinterpret_cast<const char*>(&_value), std::streamsize(sizeof(T)));
		}
	}
	
	template<class T,  typename std::enable_if<!std::is_arithmetic<T>::value, bool>::type = true>
	void write_to_stream(std::ostream& _stream, const T &_value, FileFormat _format = BINARY);
	
	template<class T, typename std::enable_if<std::is_arithmetic<T>::value, bool>::type = true>
	void write_to_stream(std::ostream& _stream, const T &_value, FileFormat _format) {
		write_to_stream<T>(_stream, _value, _format, '\t');
	}
	
	
	
	// ---------------------------------------- read_from_stream -------------------------------------------------------------
	
	template<class T, typename std::enable_if<!std::is_arithmetic<T>::value, bool>::type = true>
	void read_from_stream(std::istream& _stream, T &_obj, const FileFormat _format);
	
	template<class T, typename std::enable_if<std::is_arithmetic<T>::value, bool>::type = true>
	void read_from_stream(std::istream& _stream, T &_obj, const FileFormat _format) {
		if(_format == FileFormat::TSV) {
			_stream >> _obj;
		} else {
			_stream.read(reinterpret_cast<char*>(&_obj), sizeof(T));
		}
	}
	
	template<class T, typename std::enable_if<std::is_arithmetic<T>::value, bool>::type = true>
	T read_from_stream(std::istream& _stream, const FileFormat _format) {
		T obj;
		if(_format == FileFormat::TSV) {
			_stream >> obj;
		} else {
			_stream.read(reinterpret_cast<char*>(&obj), sizeof(T));
		}
		return obj;
	}
	
	
	
	// ---------------------------------------- save / load file -------------------------------------------------------------
	
	
	template<class T>
	void save_to_file(const T &_obj, const std::string &_filename, FileFormat _format = BINARY) {
		try {
			if (_format == FileFormat::BINARY) {
				std::ofstream out(_filename, std::ofstream::out | std::ofstream::binary);
				std::string header = std::string("Xerus ") + misc::demangle_cxa(typeid(T).name()) + " datafile.\nFormat: Binary\n";
				out.write(header.c_str(), std::streamsize(header.size()));
				write_to_stream<T>(out, _obj, _format);
			} else {
				std::ofstream out(_filename, std::ofstream::out);
				out << std::string("Xerus ") << misc::demangle_cxa(typeid(T).name()) << " datafile.\nFormat: TSV\n";
				write_to_stream<T>(out, _obj, _format);
			}
		} catch (generic_error &e) {
			throw e << "error occured while writing to file " << _filename << '\n';
		}
	}
	
	
	template<class T>
	void load_from_file(T &_obj, const std::string& _filename) {
		try {
			std::ifstream in(_filename, std::ifstream::in);
			
			std::string firstLine;
			std::getline(in, firstLine);
			
			REQUIRE(in, "Unexpected end of stream in TensorNetwork::load_from_file().");
			
			REQUIRE(firstLine == std::string("Xerus ") + misc::demangle_cxa(typeid(T).name()) + " datafile.", "Invalid binary input file " << _filename << ". DBG: " << firstLine);
			
			std::string formatQual, formatValue;
			in >> formatQual >> formatValue;
			
			REQUIRE(formatQual == std::string("Format:"), "Invalid Sytax detected in file " << _filename << ". DBG: " << formatQual);
			
			FileFormat format;
			if(formatValue == std::string("TSV")) {
				format = FileFormat::TSV;
			} else if(formatValue == std::string("Binary")) {
				format = FileFormat::BINARY;
				
				// Open the stream as binary and skip to the start of the binary part.
				const std::streamoff currPos = in.tellg();
				in.close();
				in.open(_filename, std::ifstream::in | std::ifstream::binary);
				in.seekg(currPos+1); // +1 because of the last \n
			} else {
				LOG(fatal, "Invalid value for format detected. " << formatValue);
			}
			
			read_from_stream<T>(in, _obj, format);
		} catch (generic_error &e) {
			throw e << "error occured in file " << _filename << '\n';
		}
	}
	
	template<class T>
	T load_from_file(const std::string& _filename) {
		T result;
		load_from_file<T>(result, _filename);
		return result;
	}
	
} }


