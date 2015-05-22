#pragma once

#include <type_traits>

// Adds a template argument to a function which can be of only one specified class, but allows & and && types.
#define ALLOW_MOVE(__type, __newTypeName) template<ADD_MOVE(__type, __newTypeName)>

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

