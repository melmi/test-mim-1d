all: *.h main.cpp
	g++ main.cpp -o main

doc: readme.md
	pandoc -s --mathjax -c http://kevinburke.bitbucket.org/markdowncss/markdown.css -o readme.html readme.md