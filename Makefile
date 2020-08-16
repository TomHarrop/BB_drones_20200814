graph: graph.svg

graph.svg: Snakefile
	snakemake --profile tom \
		--dag \
		|  grep -v "^[[:space:]+]0" | grep -v "\->[[:space:]]0" \
		| dot -Tsvg \
		> graph.svg
