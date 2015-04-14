FILENAME=r-workshop-series

# make all should build the pdf
all: $(FILENAME).pdf

# include the markdown files for all the lessons here.
# note, you cannot include the preface because the foreach loop below will process it, and you don't want this.
# if you change the preface, you need to make rebuild.
$(FILENAME).pdf $(FILENAME).epub: ../../lessons/r/r-intro.md ../../lessons/r/r-manipulation.md ../../lessons/r/r-viz.md ../../lessons/r/r-deseq2-demo.md

	# Use a make foreach loop to iterate through $^ (contains all pre-requisites).
	# On each file, run `sed '1,3d' $^ >> tmp.md` to remove the first three lines (front matter).
	# Assumes that you have exactly three, and only three YAML lines for each and every md file in the prerequisites.
	$(foreach f,$^,sed '1,3d' $(f) >> $(FILENAME).md;)

	# Edit the file in place to remove any image references
	gsed -i 's/.*\/img\/.*/*(image removed from PDF version)*/' $(FILENAME).md

	# Remove setup instructions
	gsed -i 's/{% include setup-r.md %}/Please see setup instructions at [bioconnector.org](http:\/\/bioconnector.org)./g' $(FILENAME).md

	# Create the pdf and epub from the preface and the md you created from all the input MD files.
	pandoc -s -V fontsize=12pt -V geometry:margin=1in -V documentclass:report --toc preface.md $(FILENAME).md -o $(FILENAME).pdf
	pandoc --toc $(FILENAME).md -o $(FILENAME).epub

	# Remove the resulting markdown
	rm $(FILENAME).md

clean:
	rm $(FILENAME).*

rebuild: clean all
