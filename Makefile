graph: graph.svg

graph.svg: Snakefile
	snakemake --profile tom \
		--forceall --dag \
		|  grep -v "^[[:space:]+]0" | grep -v "\->[[:space:]]0" \
		| dot -Tsvg \
		> graph.svg
