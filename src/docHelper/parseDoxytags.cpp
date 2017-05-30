#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <map>

using boost::property_tree::ptree;

int main() {
	ptree pt;
	read_xml("xerus.tagfile", pt);
	
	std::map<std::string, std::string> tags;
	for (const auto &compound : pt.get_child("tagfile")) {
		std::string currNamespace = "";
		if (compound.first != "compound") {
			std::cout << "*** skipping " << compound.first << " tag" << std::endl;
			continue;
		}
		std::string kind = compound.second.get<std::string>("<xmlattr>.kind", "N/A");
		if (kind == "page") {
			continue;
		}
		if (kind == "class" || kind == "struct" || kind == "namespace") {
			currNamespace = compound.second.get<std::string>("name");
		}
		
		for (const auto &subc : compound.second) {
			// we are looking for namespaces (in order! -.-) and members
			if (subc.first == "namespace") {
				currNamespace = subc.second.get<std::string>("");
			} else if (subc.first == "member") {
				std::string name = subc.second.get<std::string>("name");
				std::string subkind = subc.second.get<std::string>("<xmlattr>.kind", "N/A");
				std::string anchorfile = subc.second.get<std::string>("anchorfile");
				std::string anchor = subc.second.get<std::string>("anchor");
				std::string args = subc.second.get<std::string>("arglist");
				if (currNamespace!= "" && subkind != "define") {
					name = currNamespace+"::"+name;
				}
				boost::replace_all(name, " ", "");
				tags.insert({name, anchorfile + '#' + anchor});
			} else {
// 				std::cout << "* skipping " << subc.first << std::endl;
			}
		}
	}
	
	std::ofstream outFile("xerus.tags");
	for (const auto &t : tags) {
		outFile << t.first << "\t/doxygen/" << t.second << '\n';
	}
	
	std::cout << "### found " << tags.size() << " unique tags in xerus.tagfile created by doxygen." << std::endl;
}
