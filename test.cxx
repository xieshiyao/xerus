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

    #include <string>
    #include <iostream>

    class A {
    public:
        std::string name;
        A(const A& _other)  : name(_other.name) {
            std::cout << "Copy-Construct with name: " << name << std::endl;
        }
        A(A&& _other)       : name(std::move(_other.name)) {
            std::cout << "Move-Construct with name: " << name << std::endl;
        }
        A(std::string _name): name(_name) { }
    };

    A operator*(const A& _lhs, const A& _rhs) {
        std::cout << "Start Operator Copy with: " << _lhs.name << " " << _rhs.name << std::endl;
        A bla(_lhs.name+" "+_rhs.name);
        return bla;
    }

    A operator*(A&& _lhs, const A& _rhs) {
        std::cout << "Start Operator Move with: " << _lhs.name << " " << _rhs.name << std::endl;
        _lhs.name += " "+_rhs.name;
        return std::move(_lhs);
    }

    void foo( A&& bla) {
        std::cout << "Prelimnary result is: " << bla.name << std::endl;
    }
    int main() {
        A a("a");
        A b("b");
        A c("c");
        A d("d");
//         foo(a*b*A("t1")*c*A("t2")*A("t3")*d);
        
//         
//         A&& x = A("Horst Dieter");
        
        A x = a*b*A("t1")*c*A("t2")*A("t3")*d; 
        
        std::cout << "Final result is: " << x.name << std::endl;
    }
