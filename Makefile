20k-400g-1e/umap.png: cellsim.go main.py
	mkdir -p 20k-400g-1e
	go run cellsim.go -cells 20000 -genes 400 -umis 1200 -num-low-exp 1 -o ./20k-400g-1e/
	uv run main.py ./20k-400g-1e/

20k-400g-2e/umap.png: cellsim.go main.py
	mkdir -p 20k-400g-2e
	go run cellsim.go -cells 20000 -genes 400 -umis 1200 -num-low-exp 2 -o ./20k-400g-2e/
	uv run main.py ./20k-400g-2e/
