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

#include "../../../include/xerus/misc/callStack.h"

#ifndef NO_FANCY_CALLSTACK

#include <execinfo.h>
#include <signal.h>
#include <bfd.h>
#include <dlfcn.h>
#include <unistd.h>
#include <sstream>
#include <memory>
#include <map>
#include <vector>
#include <xerus/misc/stringUtilities.h>
START_MISC_NAMESPACE

// for the use of dladdr see also https://sourceware.org/git/?p=glibc.git;a=blob;f=debug/backtracesyms.c
struct bfdResolver {
	struct storedBfd {
		typedef bfd_boolean(deleter_t)(bfd*);
		std::unique_ptr<bfd, deleter_t*> abfd;
		std::unique_ptr<asymbol*[]> symbols;
		intptr_t offset;
		storedBfd(bfd *_abfd, deleter_t *_del) : abfd(_abfd, _del) {}
	};
	static std::map<void *, storedBfd> bfds;
	static bool bfd_initialized;
	
	static bool ensure_bfd_loaded(Dl_info &_info) {
		// load the corresponding bfd file (from file or map)
		if (bfds.count(_info.dli_fbase) == 0) {
			std::unique_ptr<storedBfd> newBfd(new storedBfd(bfd_openr(_info.dli_fname, 0), &bfd_close));
			if (!newBfd->abfd) {
				return false;
			}
			bfd_check_format(newBfd->abfd.get(),bfd_object);
			long storageNeeded = bfd_get_symtab_upper_bound(newBfd->abfd.get());
			if (storageNeeded < 0) {
				return false;
			}
			newBfd->symbols.reset(reinterpret_cast<asymbol**>(new char[(size_t)storageNeeded]));
			/*size_t numSymbols = */bfd_canonicalize_symtab(newBfd->abfd.get(), newBfd->symbols.get());
			
			newBfd->offset = (intptr_t)_info.dli_fbase;
			
			bfds.insert(std::pair<void *, storedBfd>(_info.dli_fbase, std::move(*newBfd.get())));
		} 
		return true;
	}
	
	static std::pair<uintptr_t, uintptr_t> get_range_of_section(void * _addr, std::string _name) {
		if (!bfd_initialized) {
			bfd_init();
			bfd_initialized = true;
		}
		
		// get path and offset of shared object that contains this address
		Dl_info info;
		dladdr(_addr, &info);
		if (info.dli_fbase == nullptr) {
			return std::pair<uintptr_t, uintptr_t>(0,0);
		}
		
		if (!ensure_bfd_loaded(info)) {
			return std::pair<uintptr_t, uintptr_t>(0,0);
		}
		storedBfd &currBfd = bfds.at(info.dli_fbase);
		
		asection *section = currBfd.abfd->sections;
		while (section != nullptr) {
			if (_name == section->name) {
				return std::pair<uintptr_t, uintptr_t>(section->vma, section->vma+section->size);
			} else {
				section = section->next;
				continue;
			}
		}
		return std::pair<uintptr_t, uintptr_t>(0,0);
	}
	
	static std::string resolve(void *address) {
		if (!bfd_initialized) {
			bfd_init();
			bfd_initialized = true;
		}
			
		std::stringstream res;
		res << "[0x" << std::setw((int)sizeof(void*)*2) << std::setfill('0') << std::hex << (uintptr_t)address;
		
		// get path and offset of shared object that contains this address
		Dl_info info;
		dladdr(address, &info);
		if (info.dli_fbase == nullptr) {
			return res.str()+" .?] <object to address not found>";
		}
		
		if (!ensure_bfd_loaded(info)) {
			res.str()+" .?] <could not open object file>";
		}
		storedBfd &currBfd = bfds.at(info.dli_fbase);
		
		asection *section = currBfd.abfd->sections;
		bool relative = section->vma < (uintptr_t)currBfd.offset;
// 		std::cout << '\n' << "sections:\n";
		while (section != nullptr) {
			intptr_t offset = ((intptr_t)address) - (relative?currBfd.offset:0) - (intptr_t)section->vma;
// 			std::cout << section->name << " " << section->id << " file: " << section->filepos << " flags: " << section->flags 
// 						<< " vma: " << std::hex << section->vma << " - " << std::hex << (section->vma+section->size) << std::endl;
			
			if (offset < 0 || (size_t)offset > section->size) {
				section = section->next;
				continue;
			}
			res << ' ' << section->name;
			if (!(section->flags | SEC_CODE)) {
				return res.str()+"] <non executable address>";
			}
			// get more info on legal addresses
			const char *file;
			const char *func;
			unsigned line;
			if (bfd_find_nearest_line(currBfd.abfd.get(), section, currBfd.symbols.get(), offset, &file, &func, &line)) {
				if (file) {
					return res.str()+"] "+std::string(file)+":"+to_string(line)+" (inside "+demangle_cxa(func)+")";
				} else {
					if (info.dli_saddr) {
						return res.str()+"] ??:? (inside "+demangle_cxa(func)+" +0x"+std::to_string((intptr_t)address-(intptr_t)info.dli_saddr)+")";
					} else {
						return res.str()+"] ??:? (inside "+demangle_cxa(func)+")";
					}
				}
			} else {
				return res.str()+"] <bfd_error> (inside "+demangle_cxa((info.dli_sname?info.dli_sname:""))+")";
			}
		}
// 		std::cout << " ---- sections end ------ " << std::endl;
		return res.str()+" .none] <not sectioned address>";
	}
};

std::map<void *, bfdResolver::storedBfd> bfdResolver::bfds;
bool bfdResolver::bfd_initialized = false;




std::string get_call_stack() {
	const size_t MAX_FRAMES = 100;
	std::vector<void *> stack(MAX_FRAMES);
	int num = backtrace(&stack[0], MAX_FRAMES);
	if (num <= 0) {
		return "Callstack could not be built.";
	}
	stack.resize((size_t) num);
	std::string res;
	//NOTE i=0 corresponds to get_call_stack and is omitted
	for (size_t i=1; i<(size_t)num; ++i) {
		res += bfdResolver::resolve(stack[i]) + '\n';
	}
	return res;
}

std::pair<uintptr_t, uintptr_t> get_range_of_section(void * _addr, std::string _name) {
	return bfdResolver::get_range_of_section(_addr, _name);
}



#else // no fancy callstack

#include <execinfo.h>
#include <iomanip>
#include <vector>
#include <sstream>

START_MISC_NAMESPACE

std::string get_call_stack() {
	const size_t MAX_FRAMES = 100;
	std::vector<void *> stack(MAX_FRAMES);
	int num = backtrace(&stack[0], MAX_FRAMES);
	if (num <= 0) {
		return "Callstack could not be built.";
	}
	stack.resize((size_t) num);
	std::stringstream res;
	//NOTE i=0 corresponds to get_call_stack and is omitted
	for (size_t i=1; i<(size_t)num; ++i) {
		res << "[0x" << std::setw((int)sizeof(void*)*2) << std::setfill('0') << std::hex << (uintptr_t)stack[i] << " .?] <bfd not loaded, use addr2line to resolve>\n";
	}
	return res.str();
}

std::pair<uintptr_t, uintptr_t> get_range_of_section(void * _addr, std::string _name) {
	return std::pair<uintptr_t, uintptr_t>(0,0);
}




#endif

END_MISC_NAMESPACE
