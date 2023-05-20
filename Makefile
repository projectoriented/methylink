deps: ## Install dependencies
	python -m pip install --upgrade pip .

testsuite: ## Run test suite
	pytest

%.link: ## Run test case
	methylink --threads 2 --aln tests/data/$*_aln_test-subsampled.bam --sample $* --methyl_bams tests/data/$*_methylated_test-1.bam --methyl_bams tests/data/$*_methylated_test-2.bam --output $*-linked.bam

.PHONY: clean
clean: ## Clean up
	pip uninstall methylink -y
	rm *-linked.bam*