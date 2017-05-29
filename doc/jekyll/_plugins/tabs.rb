module Jekyll
class TabsConverter < Converter
	safe true
	priority :low

	def matches(ext)
		ext =~ /^\.md$/i
	end

	def output_ext(ext)
		".html"
	end

	def convert(content)
	content.gsub('<p>__tabsInit</p>', "<input id=\"tab1\" type=\"radio\" name=\"tabs\" checked><input id=\"tab2\" type=\"radio\" name=\"tabs\">")
			.gsub('<p>__tabsStart</p>', "<div id=\"tabs\"><label for=\"tab1\">C++</label><label for=\"tab2\">Python</label><div id=\"content\"><section id=\"content1\">")
			.gsub('<p>__tabsMid</p>', "</section><section id=\"content2\">")
			.gsub('<p>__tabsEnd</p>', "</section></div></div>")
			.gsub('<p>__dangerStart</p>', "<div class=\"alert alert-danger\">")
			.gsub('<p>__dangerEnd</p>', "</div>")
			.gsub('<p>__warnStart</p>', "<div class=\"alert alert-warning\">")
			.gsub('<p>__warnEnd</p>', "</div>")
			.gsub('__breakFix1</a></p>', "")
			.gsub('<p>__breakFix2', "</a>")
			.gsub('__version', %x( git describe --tags --always --abbrev=0 ) )
			.gsub(/__doxyref\(([^\)]+)\)/){ |m| %x( ./findDoxytag #{$1} ) }
	end
end
end

