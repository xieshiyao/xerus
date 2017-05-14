module Jekyll
  class TabsConverter < Converter
    safe true
    priority :low
    @@ctr = 0

    def matches(ext)
      ext =~ /^\.md$/i
    end

    def output_ext(ext)
      ".html"
    end

    def convert(content)
      @@ctr += 1
      content.gsub('<p>__tabsInit</p>', "<input id=\"tab1\" type=\"radio\" checked><input id=\"tab2\" type=\"radio\">")
             .gsub('<p>__tabsStart</p>', "<div id=\"tabs\"><label for=\"tab1\">C++</label><label for=\"tab2\">Python</label><div id=\"content\"><section id=\"content1\">")
             .gsub('<p>__tabsMid</p>', "</section><section id=\"content2\">")
             .gsub('<p>__tabsEnd</p>', "</section></div></div>")
    end
  end
end

