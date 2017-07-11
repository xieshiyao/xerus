#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <boost/algorithm/string.hpp>

int main(int argc, char *argv[]) {
	if (argc < 2) {
		std::cerr << "### no argument given to find the tag for..." << std::endl;
		return 1;
	}
	
	std::string memberName(argv[1]);
	std::ifstream inFile("xerus.tags");
	std::string line;
	std::size_t jumpSize = 100;
	// seek rough position
	auto lastPos = inFile.tellg();
	while (!inFile.eof()) {
		lastPos = inFile.tellg();
		for (size_t i=0; i<jumpSize; ++i) {
			std::getline(inFile, line);
		}
		std::vector<std::string> parts;
		boost::split(parts, line, boost::is_any_of("\t"));
		if (parts[0] == memberName) {
			std::cout << parts[1] << std::endl;
			return 0;
		}
		if (parts[0] > memberName) {
			inFile.seekg(lastPos);
			jumpSize /= 10;
			if (jumpSize == 1) {
				break;
			}
		}
	}
	if (inFile.eof()) {
		inFile.clear();
		inFile.seekg(lastPos);
	}
	
	while (std::getline(inFile, line)) {
		std::vector<std::string> parts;
		boost::split(parts, line, boost::is_any_of("\t"));
		if (parts[0] == memberName) {
			std::cout << parts[1] << std::endl;
			return 0;
		}
	}
	
	std::cerr << "### '" << memberName << "' not found in xerus.tags" << std::endl;
	return 2;
}
