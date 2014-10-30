FILENAME=rnaseq-1day

# make all should build the pdf
all: $(FILENAME).pdf

# include the markdown files for all the lessons here.
# note, you cannot include the preface because the foreach loop below will process it, and you don't want this.
# if you change the preface, you need to make rebuild.
$(FILENAME).pdf $(FILENAME).epub: ../../lessons/shell/01-intro-unix-shell.md ../../lessons/rnaseq-1day/01-alignment-counting.md ../../lessons/intro-r-lifesci/01-intro-r.md ../../lessons/rnaseq-1day/03-differential-expression.md

	# Use a make foreach loop to iterate through $^ (contains all pre-requisites).
	# On each file, run `sed '1,3d' $^ >> tmp.md` to remove the first three lines (front matter).
	# Assumes that you have exactly three, and only three YAML lines for each and every md file in the prerequisites.
	$(foreach f,$^,sed '1,3d' $(f) >> $(FILENAME).md;)

	# Edit the file in place to remove any image references
	gsed -i 's/.*\/img\/.*/*(image removed from PDF version)*/' $(FILENAME).md

	# Create the pdf and epub from the preface and the md you created from all the input MD files.
	pandoc -s -V fontsize=12pt -V geometry:margin=1in -V documentclass:report --toc preface.md $(FILENAME).md -o $(FILENAME).pdf
	pandoc --toc $(FILENAME).md -o $(FILENAME).epub

	# Remove the resulting markdown
	rm $(FILENAME).md

clean:
	rm $(FILENAME).*

rebuild: clean all
