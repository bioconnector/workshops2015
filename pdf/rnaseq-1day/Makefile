# make all should build the pdf
all: tmp.pdf

# include the markdown files for all the lessons here
tmp.pdf: ../../lessons/shell/01-intro-unix-shell.md ../../lessons/rnaseq-1day/01-alignment-counting.md ../../lessons/rnaseq-1day/03-differential-expression.md

	# Create front matter in a new file called tmp.md
	echo "---" > tmp.md
	echo "title: RNA-seq Bootcamp Course Materials" >> tmp.md
	echo "author: Stephen D. Turner, Ph.D." >> tmp.md
	echo "date: November 10, 2014" >> tmp.md
	echo "---" >> tmp.md

	# Use a make foreach loop to iterate through $^ (contains all pre-requisites).
	# On each file, run `sed '1,3d' $^ >> tmp.md` to remove the first three lines (front matter)
	$(foreach f,$^,sed '1,3d' $(f) >> tmp.md;)

	# Edit the file in place to remove any image references
	gsed -i 's/.*\/img\/.*/*(image removed from PDF version)*/' tmp.md

	# Create the pdf
	pandoc -s -V fontsize=12pt -V geometry:margin=1in -V documentclass:report --toc tmp.md -o tmp.pdf
	rm tmp.md

clean:
	rm tmp.*

rebuild: clean all
