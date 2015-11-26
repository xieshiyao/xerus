# Sets WARNINGS
# Uses USE_CLANG, SUGGEST_ATTRIBUTES

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - Activate/Disable Additinal Warnings - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ##
#General Warning & Error options
WARNINGS += -Wall					# Show "all" warnings
WARNINGS += -Wextra					# Show more warnings
WARNINGS += -pedantic				# Strict ISO C++ checks (warn everywhere the ISO says there should be a warning)


# Additionall Warnings not aktivated by any level
WARNINGS += -Wshadow					    # Warn about variables being shadowed
WARNINGS += -Wconversion				    # warn double -> uint conversions etc.
WARNINGS += -Wswitch-default				# Warn whenever a switch statement does not have a default case. 
#WARNINGS += -Wfloat-equal			    	# Warn if floating-point values are used in equality comparisons. 
WARNINGS += -Wundef					        # Warn if an undefined identifier is evaluated in an ‘#if’ directive. 
WARNINGS += -Wunreachable-code				# Warn about unreachable code
# WARNINGS += -Wframe-larger-than=16384 -Wstack-usage=16384 # Warn if the size of a function frame is larger than 16384 bytes.
# WARNINGS += -Wunsafe-loop-optimizations	    	# Warn if the loop cannot be optimized because the compiler cannot assume anything on the bounds
# WARNINGS += -Wno-aggressive-loop-optimizations 	# Warn if in a loop with constant number of iterations the compiler detects undefined behavior
WARNINGS += -Wdisabled-optimization	                # Warn if a requested optimization pass is disabled
# WARNINGS += -Wpadded -Wpacked			            # misaligned structs (reordering of the fields might help)
# WARNINGS += -Winline			                    # Warn if a function that is declared as inline cannot be inlined
# WARNINGS += -fstack-protector -Wstack-protector	# warns about functions that are not protected against stack smashing
ifdef USE_CLANG
	WARNINGS += -Wlogical-op-parentheses			# Warn about suspicious uses of logical operators in expressions
	WARNINGS += -Wno-error=return-type-c-linkage		# do not warn about c-incompatible return types in "extern-c" blocks
else
	WARNINGS += -Wuseless-cast 				# Warn when an expression is casted to its own type. 
	WARNINGS += -Wno-error=useless-cast			# fails on gcc 5.1.1 otherwise
	WARNINGS += -Wlogical-op				# Warn about suspicious uses of logical operators in expressions
	WARNINGS += -Wtrampolines				# Warn about trampolines
	ifdef SUGGEST_ATTRIBUTES
		WARNINGS += -Wsuggest-attribute=pure 			# Suggest functions that can be pure
		WARNINGS += -Wsuggest-attribute=const 			# Suggest functions that can be const
		WARNINGS += -Wsuggest-attribute=noreturn 		# Suggest functions that do not return
		WARNINGS += -Wno-error=suggest-attribute=pure 		# Do not promote to Error
		WARNINGS += -Wno-error=suggest-attribute=const 		# Do not promote to Error
		WARNINGS += -Wno-error=suggest-attribute=noreturn 	# Do not promote to Error
	endif
endif

# Disabled Warnings
WARNINGS += -Wno-comment				    # Ignore warning about /* inside /*...*/
WARNINGS += -Wno-unknown-pragmas			# Allow unknwon pragmas (e.g. openmp)
WARNINGS += -Wno-unused-parameter			# No warning about unused parameters

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - Error Options - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ##
WARNINGS += -Werror				        	# Promote ALL Warnings to Errors

ifdef USE_CLANG
	WARNINGS += -ferror-limit=3				# Abort at third error
else
	WARNINGS += -fmax-errors=3				# Abort at third error
endif

# Exceptions from Error Promotion
# WARNINGS += -Wno-error=parentheses			        # Do not Promote to Error
WARNINGS += -Wno-error=unused-parameter			    # Do not Promote to Error
WARNINGS += -Wno-error=unused-variable			    # Do not Promote to Error
WARNINGS += -Wno-error=disabled-optimization		# Do not Promote to Error
# WARNINGS += -Wno-error=frame-larger-than		    # Do not Promote to Error
# WARNINGS += -Wno-error=stack-usage			    # Do not Promote to Error

ifdef SUGGEST_ATTRIBUTES
	WARNINGS += -Wno-error=suggest-attribute=pure 		# Do not promote to Error
	WARNINGS += -Wno-error=suggest-attribute=const 		# Do not promote to Error
	WARNINGS += -Wno-error=suggest-attribute=noreturn 	# Do not promote to Error
endif
	
ifndef USE_CLANG
	WARNINGS += -Wno-error=unused-but-set-variable	# Do not Promote to Error
	WARNINGS += -Wno-error=maybe-uninitialized	    # Do not Promote to Error
endif
