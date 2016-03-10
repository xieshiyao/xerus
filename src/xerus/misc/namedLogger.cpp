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
 * @brief Implementation of the log-buffer functionality and declaration of global logging objects.
 */

#include <time.h>

#ifdef LOG_BUFFER_
    #include <fstream>
#endif

#include <xerus/misc/namedLogger.h>

namespace xerus {
    namespace misc {
        namespace internal {
            //Look if Log shall be streamed into a file
            #ifdef LOGFILE_
                std::ofstream fileStream ( "error.log" , std::ofstream::app | std::ofstream::out);
            #endif
                
                std::mutex namedLoggerMutex;
                std::string logFilePrefix;
                bool silenced = false;
				std::chrono::system_clock::time_point programStartTime;
				static void __attribute__((constructor)) initTime() {
					programStartTime = std::chrono::system_clock::now();
				}
                
            #ifdef LOG_BUFFER_
				namespace buffer {
					std::stringstream current;
					std::stringstream old;
					
					void checkSwitch() {
						if (current.str().size() > 1024*1024) {
							old = std::move(current);
							current = std::stringstream();
						}
					}
					
					void dump_log(std::string _comment)  {
						std::string name = std::string("errors/") + logFilePrefix + std::to_string(std::time(nullptr)) + ".txt";
						std::ofstream out(name, std::ofstream::out);
						out << "Error: " << _comment << std::endl << std::endl;

						// Get callstack
						out << "-------------------------------------------------------------------------------" << std::endl 
							<< "  Callstack : " << std::endl
							<< "-------------------------------------------------------------------------------" << std::endl
							<< get_call_stack() << std::endl;
						
						// Output namedLogger
						current.flush();
						old.flush();
						out << "-------------------------------------------------------------------------------" << std::endl 
							<< "  last " << (current.str().size() + old.str().size()) << " bytes of log:" << std::endl
							<< "-------------------------------------------------------------------------------" << std::endl 
							<< current.str() << old.str() << "horst" << std::endl; 
						out.close();
					}
				}
            #endif
        }
    }
}
