.PHONY: check demo audit-paths

check:
	bash scripts/smoke_test.sh

demo:
	python3 demo/run_demo.py

audit-paths:
	@echo "Remaining absolute paths in *.R / notebooks (informational):"
	@grep -RInE '/Users/|/mnt/large/' --include='*.R' --include='*.ipynb' . || true
