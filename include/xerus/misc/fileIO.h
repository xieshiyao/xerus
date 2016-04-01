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
#include "exceptions.h"
#include "namedLogger.h"

namespace xerus { namespace misc {
	
	/**
	 * @brief possible file formats for tensor storage
	 * @details 
	 * BINARY: binary output, exact storage
	 * TSV: tab-separated-value lists
	 */
	enum class FileFormat { BINARY, TSV };
	
	// ---------------------------------------- write_to_stream -------------------------------------------------------------
	
	template<class T, class enable = void>
	struct StreamWriter {
		void operator()(std::ostream& _stream, const T& _value, const FileFormat _format) const;
	};
	
	template<class T, typename std::enable_if<std::is_arithmetic<T>::value, bool>::type = true>
	void write_to_stream(std::ostream& _stream, const T &_value, FileFormat _format) {
		if(_format == FileFormat::TSV) {
			_stream << _value << '\t';
		} else {
			_stream.write(reinterpret_cast<const char*>(&_value), std::streamsize(sizeof(T)));
		}
	}
	
	template<class T, typename std::enable_if<!std::is_arithmetic<T>::value, bool>::type = true>
	void write_to_stream(std::ostream& _stream, const T& _value, const FileFormat _format = FileFormat::BINARY) {
		StreamWriter<T> writerObj;
		writerObj(_stream, _value, _format);
	}
	
	
	template<class item_t, class... rest_t>
	struct StreamWriter<std::vector<item_t, rest_t...>> {
		void operator()(std::ostream& _stream, const std::vector<item_t, rest_t...> &_value, const FileFormat _format) const {
			write_to_stream<uint64>(_stream, _value.size(), _format);
			for (size_t i = 0; i < _value.size(); ++i) {
				write_to_stream<item_t>(_stream, _value[i], _format);
			}
		}
	};
	
	
	// ---------------------------------------- read_from_stream -------------------------------------------------------------
	
	template<class T, class enable = void>
	struct StreamReader {
		void operator()(std::istream& _stream, T &_obj, const FileFormat _format) const;
	};
	
	template<class T>
	struct StreamReader<T, typename std::enable_if<std::is_arithmetic<T>::value>::type> {
		void operator()(std::istream& _stream, T &_obj, const FileFormat _format) const {
			if(_format == FileFormat::TSV) {
				_stream >> _obj;
			} else {
				_stream.read(reinterpret_cast<char*>(&_obj), sizeof(T));
			}
		}
	};
	
	template<class T, class enable = void>
	void read_from_stream(std::istream& _stream, T &_obj, const FileFormat _format) {
		static StreamReader<T> callObj;
		callObj(_stream, _obj, _format);
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
	
	
	
	
	template<class item_t, class... rest_t>
	struct StreamReader<std::vector<item_t, rest_t...>> {
		void operator()(std::istream& _stream, std::vector<item_t, rest_t...> &_value, const FileFormat _format) const {
			_value.resize(read_from_stream<uint64>(_stream, _format));
			for (size_t i = 0; i < _value.size(); ++i) {
				read_from_stream<item_t>(_stream, _value[i], _format);
			}
		}
	};
	
	
	// ---------------------------------------- save / load file -------------------------------------------------------------
	
	
	template<class T>
	void save_to_file(const T &_obj, const std::string &_filename, const FileFormat _format = FileFormat::BINARY) {
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
	
}}
