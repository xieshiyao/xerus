# ------------------------------------------------------------------------------------------------------
#				Default rule should be the help message
# ------------------------------------------------------------------------------------------------------
help:
	@printf "Possible make targets are:\n \
	\t\thelp \t\t -- Print this help.\n \
	\t\tdoc \t\t -- Build the html documentation for the xerus library.\n \
	\t\tserve \t\t -- Build the html documentation for the xerus library and offer it via 'jekyll serve'.\n \
	\t\tclean \t\t -- Remove all documentation files.\n"

doc:
	-mkdir html
	doxygen doxygen/Doxyfile
	jekyll build --source jekyll/ --destination html/

clean:
	-rm -r html

serve:
	-mkdir html
	doxygen doxygen/Doxyfile
	jekyll serve --source jekyll/ --destination html/
