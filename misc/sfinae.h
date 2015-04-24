#pragma once

#include <type_traits>

// Adds a template argument to a function which can be of only one specified class, but allows & and && types.
#define ALLOW_MOVE(__type, __newTypeName) template<class __newTypeName, typename std::enable_if<std::is_base_of<__type, typename std::decay<__newTypeName>::type>{}, int>::type = 0>

// Adds two template arguments to a function which can each be of only one specified class, but allows & and && types.
#define ALLOW_MOVE_TWO(__newTypeName, __allowedType1, __allowedType2) \
    template<class __newTypeName, \
        typename std::enable_if<\
               std::is_base_of<__allowedType1, typename std::decay<__newTypeName>::type>{} \
            || std::is_base_of<__allowedType2, typename std::decay<__newTypeName>::type>{}, \
            int> \
        ::type = 0 \
    >

// Adds a template arguments whithin an existing template argument list, which can be of only one specified class, but allows & and && types.
#define ADD_MOVE(__type, __newTypeName) class __newTypeName, typename std::enable_if<std::is_base_of<__type, typename std::decay<__newTypeName>::type>{}, int>::type = 0



// Copyright Notice: The following macro is licenced under the CC-BY-SA licence and was created by:
// Wikibooks contributors, "More C++ Idioms/Member Detector," Wikibooks, The Free Textbook Project,
// http://en.wikibooks.org/w/index.php?title=More_C%2B%2B_Idioms/Member_Detector&oldid=2834755 (accessed April 13, 2015).

// Macro to create checks for the existence of member functions
#define GENERATE_HAS_MEMBER(member)                                               \
                                                                                  \
template < class T >                                                              \
class HasMember_##member                                                          \
{                                                                                 \
private:                                                                          \
    using Yes = char[2];                                                          \
    using  No = char[1];                                                          \
                                                                                  \
    struct Fallback { int member; };                                              \
    struct Derived : T, Fallback { };                                             \
                                                                                  \
    template < class U >                                                          \
    static No& test ( decltype(U::member)* );                                     \
    template < typename U >                                                       \
    static Yes& test ( U* );                                                      \
                                                                                  \
public:                                                                           \
    static constexpr bool RESULT = sizeof(test<Derived>(nullptr)) == sizeof(Yes); \
};                                                                                \
                                                                                  \
template < class T >                                                              \
struct has_member_##member                                                        \
: public std::integral_constant<bool, HasMember_##member<T>::RESULT>              \
{ };                                                                              \

