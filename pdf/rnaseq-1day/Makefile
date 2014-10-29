FILENAME=rnaseq-1day

# make all should build the pdf
all: $(FILENAME).pdf

# include the markdown files for all the lessons here
$(FILENAME).pdf $(FILENAME).epub: ../../lessons/shell/01-intro-unix-shell.md ../../lessons/rnaseq-1day/01-alignment-counting.md ../../lessons/rnaseq-1day/03-differential-expression.md

	# Create front matter in a new file called tmp.md
	echo "---" > $(FILENAME).md
	echo "title: RNA-seq Bootcamp Course Materials" >> $(FILENAME).md
	echo "author: Stephen D. Turner, Ph.D." >> $(FILENAME).md
	echo "date: November 10, 2014" >> $(FILENAME).md
	echo "---" >> $(FILENAME).md

	# Use a make foreach loop to iterate through $^ (contains all pre-requisites).
	# On each file, run `sed '1,3d' $^ >> tmp.md` to remove the first three lines (front matter)
	$(foreach f,$^,sed '1,3d' $(f) >> $(FILENAME).md;)

	# Edit the file in place to remove any image references
	gsed -i 's/.*\/img\/.*/*(image removed from PDF version)*/' $(FILENAME).md

	# Create the pdf and epub
	pandoc -s -V fontsize=12pt -V geometry:margin=1in -V documentclass:report --toc $(FILENAME).md -o $(FILENAME).pdf
	pandoc --toc $(FILENAME).md -o $(FILENAME).epub

clean:
	rm $(FILENAME).*

rebuild: clean all
