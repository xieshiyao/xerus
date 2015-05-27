// Xerus - A General Purpose Tensor Library
// Copyright (C) 2014-2015 Benjamin Huber and Sebastian Wolf. 
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

#include <time.h>
#include <unistd.h>

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
                
            #ifdef LOG_BUFFER_
                BufferStreams bufferStreams;
                
                void dump_log_buffer(std::string _comment)  {
                    std::string name = std::string("errors/") + logFilePrefix + std::to_string(std::time(0)) + ".txt";
                    std::ofstream out(name, std::ofstream::out);
                    out << "Error: " << _comment << std::endl << std::endl;

                    // get callstack
                    out << "-------------------------------------------------------------------------------" << std::endl 
                        << "  Callstack : " << std::endl
                        << "-------------------------------------------------------------------------------" << std::endl;
                    out.close();
                    // gdb /proc/6345/exe 6345 -batch -ex 'thread apply all backtrace' -q >> errors/138344543.txt
                    std::stringstream tmp;
                    tmp << "gdb /proc/" << getpid() << "/exe " << getpid() << " -batch -ex 'thread apply all backtrace' -q >> " << name << " 2>&1";
                    system(tmp.str().c_str());
                    out.open(name, std::ofstream::out | std::ofstream::app);
                    out << std::endl << std::endl;
                    
                    // output namedLogger
                    err::bufferStreams.curr->flush();
                    err::bufferStreams.old->flush();
                    out << "-------------------------------------------------------------------------------" << std::endl 
                        << "  last " << (err::bufferStreams.curr->str().size() + err::bufferStreams.old->str().size()) << " bytes of log:" << std::endl
                        << "-------------------------------------------------------------------------------" << std::endl 
                        << err::bufferStreams.curr->str() << err::bufferStreams.old->str() << "horst" << std::endl; 
                    out.close();
                }
            #endif
        }
    }
}
