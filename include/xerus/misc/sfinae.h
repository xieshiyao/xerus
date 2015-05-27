#pragma once

#include <type_traits>

// Adds a template argument to a function which can be of only one specified class, but allows & and && types.
#define ALLOW_MOVE(oldTypeName, newTypeName) template<ADD_MOVE(oldTypeName, newTypeName)>

// Adds two template arguments to a function which can each be of only one specified class, but allows & and && types.
#define ALLOW_MOVE_TWO(newTypeName, allowedType1, allowedType2) \
    template<class newTypeName, \
        typename std::enable_if<\
               std::is_base_of<allowedType1, typename std::decay<newTypeName>::type>{} \
            || std::is_base_of<allowedType2, typename std::decay<newTypeName>::type>{}, \
            int> \
        ::type = 0 \
    >

// Adds a template arguments whithin an existing template argument list, which can be of only one specified class, but allows & and && types.
#define ADD_MOVE(oldTypeName, newTypeName) class newTypeName, typename std::enable_if<std::is_base_of<oldTypeName, typename std::decay<newTypeName>::type>{}, int>::type = 0



// Macro to create checks for the existence of member functions
#define GENERATE_HAS_MEMBER(member)                                               \
                                                                                  \
template < class T >                                                              \
class HasMember_##member                                                          \
{                                                                                 \
private:                                                                          \
    typedef char (& yes)[1];							\
    typedef char (& no)[2];							\
										\
    template <typename C> static yes check(decltype(&C::member));		\
    template <typename> static no check(...);					\
										\
public: 									\
    static constexpr bool RESULT = sizeof(check<T>(0)) == sizeof(yes);		\
};                                                                                \
                                                                                  \
template < class T >                                                              \
struct has_member_##member                                                        \
: public std::integral_constant<bool, HasMember_##member<T>::RESULT>              \
{ };                                                                              \

